wgbs.simulate.reads <- function(regions.gr, name.col = NULL, genome, 
                                template.size, read.length, convert.rate, seed = NULL,
                                convert.strand = c("T", "B"), out.root, suffix = ".fq") {
  insert.size <- template.size - 2*read.length
  fa <- paste0("~/genomes/", genome, "/", genome, ".fa")
  if (!file.exists(fa)) {
    stop("!file.exists(fa)")
  }
  dir.create(dirname(out.root), showWarnings = F, recursive = T)
  if (is.null(name.col)) {
    name.col <- "regionID"
    regions.gr$regionID <- paste0("region_", 1:length(regions.gr))
  }
  regions.gr$regionID <- mcols(regions.gr)[, name.col]
  if (any(width(regions.gr) < template.size)) {
    stop("any(width(regions.gr) < template.size)")
  }
  utilsFanc::t.stat("getting sequences")
  seqs.by.region <- lapply(1:length(regions.gr), function(i) {
    region <- regions.gr[i]
    mate1.strand <- "+"
    mate2.strand <- "-"
    df <- data.frame(chr1 = seqnames(region) %>% as.character(),
                     start1 = start(region):(end(region) - template.size + 1))
    df <- df %>% dplyr::mutate(end1 = start1 + read.length - 1,
                               strand1 = mate1.strand,
                               chr2 = chr1, 
                               start2 = start1 + insert.size + read.length,
                               end2 = start1 + template.size - 1,
                               strand2 = mate2.strand)
    df$regionID <- region$regionID
    df <- df %>% dplyr::mutate(rname = paste0(chr1, ":", start1, "-", end2))
    seqs <- lapply(c(1, 2), function(i) {
      df.mate <- df[, c(paste0(c("chr", "start", "end", "strand"), i), "regionID", "rname")]
      colnames(df.mate) <- colnames(df.mate) %>% sub("[12]", "", .)
      gr.mate <- makeGRangesFromDataFrame(df.mate, keep.extra.columns = T)
      gr.mate$regionID <- paste0(gr.mate$regionID, ":", gr.mate$rname)
      seq.vec <- abaFanc2::get.fasta.gr(gr = gr.mate, genome = genome,
                                        as.string = T, id.col = "regionID", 
                                        print.cmd = F) %>% unlist()
      return(seq.vec)
    })
    names(seqs) <- c("R1", "R2")
    return(seqs)
  })
  seqs.by.strand <- lapply(c("R1", "R2"), function(i) {
    lapply(seqs.by.region, function(x) {
      return(x[[i]])
    }) %>% unlist() %>% return()
  })
  names(seqs.by.strand) <- c("R1", "R2")
  rm(seqs.by.region)
  utilsFanc::t.stat("bisulfite convertion...")
  seqs.converted <- lapply(convert.strand, function(c.s) {
    from.list <- list()
    to.list <- list()
    from.list[["R1"]] <- "C"
    from.list[["R2"]] <- "G"
    to.list[["R1"]] <- "T"
    to.list[["R2"]] <- "A"
    
    seq.source <- list()
    if (c.s == "T") {
      seq.source[["R1"]] <- "R1"
      seq.source[["R2"]] <- "R2"
    } else if (c.s == "B") {
      seq.source[["R2"]] <- "R1"
      seq.source[["R1"]] <- "R2"
    } else {
      stop ("convert.strand has to be T or B")
    }
    res <- lapply(c("R1", "R2"), function(r) {
      seqs <- wgbs.convert.core(seqs = seqs.by.strand[[seq.source[[r]]]],
                                from = from.list[[r]], to = to.list[[r]], 
                                conversion.rate = convert.rate, seed = seed)
      return(seqs)
    })
    
    # programming logic see labnotes for 2022-01-11.
    names(res) <- c("R1", "R2")
    return(res)
  })
  names(seqs.converted) <- convert.strand
  
  seqs.out <- seqs.converted
  seqs.out$O <- seqs.by.strand
  utilsFanc::t.stat("writing sequences out ...")
  fastqs.list <- lapply(names(seqs.out), function(type) {
    seqs <- seqs.out[[type]]
    fastqs <- liteRnaSeqFanc::fastq.write(seq.list = seqs, out.root = paste0(out.root, "_", type, "_R"), 
                                out.suffix = ".fq", zip = T)
    return(fastqs)
  })
  names(fastqs.list) <- names(seqs.out)
  return(fastqs.list)
}

wgbs.convert.core <- function(seqs, from, to, conversion.rate,
                              threads = 1, seed = NULL) {
  # reads should be written in the form of strings
  if (length(unique(nchar(seqs))) != 1) {
    stop("unique(nchar(seqs)) != 1")
  }
  if (any(stringr::str_detect(seqs,"[[:lower:]]"))) {
    warning("lower case detected. Convert all to upper cases")
    seqs <- toupper(seqs)
  }
  from <- toupper(from)
  to <- toupper(to)
  # the idea is to put strings into matrices, and use matrix operations
  
  mat <- simplify2array(x = strsplit(seqs, "")) # simplify2array by default add by column
  set.seed(seed = seed)
  bMat.convert <- sample(c(T, F), size = nrow(mat) * ncol(mat), replace = T, 
                         prob = c(conversion.rate, 1-conversion.rate)) %>% 
    matrix(nrow = nrow(mat), ncol = ncol(mat))
  bMat.from <- mat == from
  bMat.convert.from <- bMat.from & bMat.convert
  mat[bMat.convert.from] <- to
  seqs.out <- mat %>% split(., c(col(.)))
  seqs.out <- sapply(seqs.out, function(x) return(paste0(x, collapse = ""))) %>% 
    `names<-`(names(seqs))
  # seqs.out <- utilsFanc::safelapply(seqs, function(seq) {
  #   vec <- seq %>% strsplit("") %>% unlist()
  #   from.pos <- which(vec == from)
  #   convert.pos <- sample(from.pos, floor(conversion.rate * length(from.pos)))
  #   vec[convert.pos] <- to
  #   return(paste0(vec, collapse = ""))
  # }, threads = threads) %>% unlist()
  return(seqs.out)
}

wgbs.parse.qname <- function(qname, return.df = F, field) {
  regionID <- gsub(":.+", "", qname)
  loci <- sub("[^:]+:", "", qname)
  df <- utilsFanc::loci.2.df(loci.vec = loci)
  df$regionID <- regionID
  if (return.df == T)
    return(df)
  res <- df[, field]
  return(res)
}

wgbs.parse.simulated.bam <- function(bam, regions.gr, mode = "start", reverse.coverage = T,
                                     read.length,
                                     write.type = c("on", "off", "miss"),
                                     out.dir) {
  if (mode == "start") {
    bam.gr <- Rsamtools::scanBam(file = bam) %>% as.data.frame() %>% 
      dplyr::filter(mpos > pos) %>% # only take left mate
      dplyr::mutate(chr = rname, start = pos, end = pos, strand = "+") %>% 
      dplyr::select(chr, start, end, strand, qname, flag, mapq) %>% 
      makeGRangesFromDataFrame(keep.extra.columns = T)
    in.range <- subsetByOverlaps(bam.gr, regions.gr, ignore.strand=T)
    bOn.target <- sub("-.+$", "", utilsFanc::gr.get.loci(in.range, keep.strand = F)) == 
      sub(":", "", stringr::str_extract(in.range$qname, ":.+:\\d+"))
  } else if (mode == "coverage") {
    bam.gr <- Rsamtools::scanBam(file = bam) %>% as.data.frame() %>%
      dplyr::mutate(chr = rname, start = pos, end = pos + read.length - 1) %>% 
      dplyr::select(chr, start, end, strand, qname, flag, mapq) %>% 
      makeGRangesFromDataFrame(keep.extra.columns = T)
    in.range <- subsetByOverlaps(bam.gr, regions.gr, ignore.strand=T)
    qname.df <- wgbs.parse.qname(in.range$qname, return.df = T)
    
    bOn.target <- (as.character(seqnames(in.range)) == qname.df$chr) &
                     ( start(in.range) == qname.df$start  | 
                         end(in.range) == qname.df$end)
  }

  bOff.target <- ! bOn.target
  out.grs <- list()
  out.grs$on <- in.range[bOn.target]
  out.grs$off <- in.range[bOff.target]
  out.grs$miss <- setdiff(regions.gr, out.grs$on, ignore.strand=TRUE) 
  dir.create(out.dir, showWarnings = F, recursive = T)
  
  lapply(write.type, function(type) {
    gr <- out.grs[[type]]
    if (type == "on" && length(gr) > 1)
      file.type <- "bedGraph" # wanted to use bw. but not working currently.
    else
      file.type <- "bedGraph"
    
    out.file <- paste0(out.dir, "/", basename(bam), "_", mode, "_", type, ".", file.type)
    
    if (length(gr) > 1) {
      if (mode == "start") {
        mcols(gr) <- NULL
        gr$score <- 1
      } else if (mode == "coverage") {
        gr <- coverage(gr) %>% as("GRanges")
        if (reverse.coverage == T && type == "on")
          gr$score <- 2*read.length - gr$score
      }
      rtracklayer::export(gr, out.file)
    } else {
      df <- data.frame(chr = "chr1", start = 1, end = 2, score = 0) %>% .[F,] %>% 
        write.table(out.file, sep = "\t", quote = F, row.names = F, col.names = F)
    }
    if (file.type == "bedGraph") {
      system(paste0("~/scripts/bed_browser_v2.sh ", out.file))
    }
    return()
  })
  return()
}

gr.bp.resolution <- function(gr) {
  # convert a gr to bp level
  
}
