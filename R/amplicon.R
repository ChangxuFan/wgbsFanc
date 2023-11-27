amplicon.clean <- function(bams, out.bams = NULL, out.dir = NULL, 
                           flags = BISMARK.TAGS,
                           mapq.th = NULL, target.gr = NULL,
                           # no.indel = F, no.misMatch = F,
                           subset.prob = NULL, # generate a subset for IGV viewing
                           buffers = c(0, 3, 3, 0),
                           threads = 1, samtools = SAMTOOLS) {
  # buffers: left outside, left inside, right inside, right outside.
  # for left mate, its left end must be with in (target.gr$start - left.outside, target.gr$end + left.inside)
  if (is.null(out.bams)) {
    if (is.null(out.dir)) {
      out.dir <- dirname(bams[1])
    }
    rootnames <- basename(bams) %>% sub(".bam$", "", .)
    out.bams <- paste0(out.dir, "/", rootnames, "_cleaned.bam")
  }
  if (length(bams) != length(out.bams)) {
    stop("length(bams) != length(out.bams)")
  }

  filters <- list()
  if (!is.null(mapq.th)) {
    filters$mapq <- function(x) x$mapq > mapq.th
  }
  if (!is.null(target.gr)) {
    if (length(target.gr) != 1) {
      stop("length(target.gr) != 1")
    }
    start.range <- c(start(target.gr) - buffers[1], start(target.gr) + buffers[2])
    end.range <- c(end(target.gr) - buffers[3], end(target.gr) + buffers[4])
    filters$target <- function(x) {
      x$pos.end <- end(unlist(GenomicAlignments::extractAlignmentRangesOnReference(x$cigar, x$pos)))
      res <- x$rname == as.character(seqnames(target.gr)) &
        ((x$mpos >= x$pos & x$pos >= start.range[1] & x$pos <= start.range[2] ) | 
        (x$mpos < x$pos & x$pos.end >= end.range[1] & x$pos.end <= end.range[2] ))
      return(res)
    }
  }
  
  subset.filter <- list(subset = function(x) return(sample(c(T, F), size = nrow(x), replace = T,
                                                           prob = c(subset.prob, 1-subset.prob))))
  if (length(filters) < 1) {
    stop("no filters supplied")
  }
  res <- utilsFanc::safelapply(1:length(bams), function(i) {
    bam <- bams[i]
    out.bam <- out.bams[i]
    Rsamtools::filterBam(file = bam, destination = out.bam, filter = FilterRules(filters)) 
    if (!is.null(subset.prob)) {
      Rsamtools::filterBam(file = out.bam, destination = utilsFanc::insert.name.before.ext(name = out.bam, insert = "subset", delim = "_"),
                           filter = FilterRules(subset.filter)) 
    }
    res <- bamFanc::remove.mate(bam = out.bam)
    return(res)
  }, threads = threads) %>% unlist()
  names(res) <- names(bams)
  return(res)
}


amplicon.pos.read.mat.gen <- function(bams, pos.gr, no.na = T, 
                             threads.bam = 1, threads.pos = 1) {
  if (is.null(names(bams))) {
    names(bams) <- basename(bams)
  }
  mats <- utilsFanc::safelapply(bams, function(bam) {
    df <- utilsFanc::safelapply(1:length(pos.gr), function(i) {
      pos <- pos.gr[i]
      if (width(pos) != 1) {
        stop("width(pos) != 1")
      }
      stack <- stackStringsFromBam(file = bam, 
                                   param = utilsFanc::gr.get.loci(pos),
                                   use.names = T) %>% as.matrix()
      df <- data.frame(qname = rownames(stack), nuc = stack[, 1])
      colnames(df) <- c("qname", paste0(seqnames(pos), ":", start(pos)))
      df <- df[!duplicated(df$qname), ]
      ####
      # this is to deal with scenarios where the 2 mates overlap, but do not agree.
      # from my data, sometimes they do not agree. 
      # but I chose to use the left mate, because the right mate has extensive mismatches
      # around this site. The base quality of the right mate is also low.
      # since the left mate is always in front of the right mate in a sorted bam file,
      # this works...
      ###
      return(df)
    }, threads = threads.pos) %>% Reduce(full_join, .)
    if (no.na) {
      df <- df %>% na.omit()
    }
    rownames(df) <- df$qname
    df$qname <- NULL
    mat <- as.matrix(df)
    return(mat)
  }, threads = threads.bam)
  names(mats) <- names(bams)
  return(mats)
}
# t <- Rsamtools::scanBam(file = "/scratch/fanc/4dn/nk/wgbs/amplicon/v1/bismark_xiaoyu_nd/DHneg_R1_bismark_bt2_pe.sorted.bam")
# 
# filts <- c("peaks", "promoters")
# filters <- FilterRules(filts)
# active(filters) # all TRUE
# 
# df <- DataFrame(peaks = c(TRUE, TRUE, FALSE, FALSE),
#                 promoters = c(TRUE, FALSE, FALSE, TRUE),
#                 introns = c(TRUE, FALSE, FALSE, FALSE))
# eval(filters, df)
# 
# filts <- list(peaks = expression(peaks), promoters = expression(promoters),
#               find_eboxes = function(rd) rep(FALSE, nrow(rd)))
# filters <- FilterRules(filts, active = FALSE)
# eval(filters, df)
# t.bam <- "/scratch/fanc/4dn/nk/wgbs/amplicon/v1/bismark_xiaoyu_nd/DHneg_R1_bismark_bt2_pe.sorted.bam"
# t.out.bam <- "/scratch/fanc/4dn/nk/wgbs/amplicon/v1/bismark_xiaoyu_nd/t_filtered.bam"
# Rsamtools::filterBam(file = t.bam, destination = t.out.bam, 
#                      filter = FilterRules(list(mapq = function(x) x$mapq > 30,
#                                                end = function(x) {
#                                                  end <- end(unlist(GenomicAlignments::extractAlignmentRangesOnReference(x$cigar, x$pos)))
#                                                  return(end < 130129943)
#                                                })))
# 
# t <- GenomicAlignments::extractAlignmentRangesOnReference(c("7M1I6M1I2M3I230M", "7M1I6M1I2M3I230M"), pos = c(1, 100))
# 
amplicon.read.nuc.sum <- function(mats, CT.only = T) {
  df <- lapply(names(mats), function(sample) {
    mat <- mats[[sample]]
    if (CT.only) {
      mat[! mat %in% c("C", "T")] <- NA
      mat <- mat %>% na.omit()
    }
    n <- apply(mat, 1, function(x) sum(x == "C"))
    distro <- table(n) %>% as.data.frame()
    colnames(distro) <- c("n_methylated", sample)
    return(distro)
  }) %>% Reduce(full_join, .)
  df[is.na(df)] <- 0
  return(df)
}
