dss.input <- function(methylC, me.grs = NULL, threads = 4, region = NULL) {
  if (is.null(me.grs)) {
    me.grs <- methylC.import(methylC = methylC, threads = threads, region = region)
  }
  dss <- utilsFanc::safelapply(me.grs, function(gr) {
    df <- gr %>% as.data.frame() %>% 
      dplyr::select(seqnames, start, CG, depth) %>% 
      dplyr::mutate(X = ceiling(CG * depth)) %>% 
      dplyr::mutate(CG = NULL) %>% 
      dplyr::rename(chr = seqnames, pos = start, N = depth)
    return(df)
  }, threads = threads)
  return(dss)
}

dss.quick.check <- function(dml, locus = "chr6:130129523-130130584") {
  locus.gr <- utilsFanc::loci.2.gr(locus)
  if ("pos" %in% colnames(dml)) {
    df <- dml %>% `class<-`("data.frame") %>% 
      dplyr::mutate(start = pos, end = pos) %>% 
      dplyr::mutate(pos = NULL)
  } else
    df <- dml
  df <- df %>% filter(chr == as.character(seqnames(locus.gr)[1]),
                      start >= start(locus.gr),
                      end <= end(locus.gr))
  return(as_tibble(df))
}

dml.filter <- function(dml, dss.in.list = NULL, min.cov = NULL) {
  if (!is.null(dss.in.list)) {
    if (!is.null(min.cov)) {
      print(paste0("filtering for min coverage of ", min.cov, " in all samples"))
      dss.in.list <- lapply(dss.in.list, function(x) {
        return(x[x$N >= min.cov,])
      })
    }
    # ij means inner join
    pos.ij <- dss.in.list %>% lapply(function(x) return(x[, c("chr", "pos")])) %>% 
      Reduce(inner_join, .)
    dml <- dml %>% `class<-`("data.frame") %>% inner_join(pos.ij) %>% 
      `class<-`(c("data.frame", "DMLtest"))
  }
  return(dml)
}

dmr.list.write <- function(dmr.list, items = NULL, 
                           out.dir) {
  if (is.null(items)) {
    items <- names(dmr.list)
  }
  lapply(items, function(item) {
    lapply(names(dmr.list[[item]]), function(comp) {
      dmr <- dmr.list[[item]][[comp]]
      utilsFanc::write.zip.fanc(df = dmr, out.file = paste0(out.dir,"/", item, "_", comp, ".bed"), 
                                bed.shift = T)
      utilsFanc::write.zip.fanc(df = dmr, out.file = paste0(out.dir,"/", item, "_", comp, ".tsv"),
                                zip = F, col.names = T,
                                bed.shift = F)
      return()
    })
  })
}

dmr.pct.in.peaks <- function(dmr, peaks.gr) {
  dmr <- makeGRangesFromDataFrame(dmr)
  o <- findOverlaps(dmr, peaks.gr, ignore.strand = T)
  n <- queryHits(o) %>% unique() %>% length()
  return(n/length(dmr))
}

dmr.bed.import <- function(dmr.bed) {
  dmrs <- read.table(dmr.bed)
  header <- c("chr", "start", "end", "length", "nCG", "meanMethy1", "meanMethy2", "diff.Methy", "areaStat")
  colnames(dmrs) <- header
  return(dmrs)
}

dmr.pileup.assess <- function(dmr.file, methyl.bws, out.dir = NULL, root.name) {
  stop("unfinished")
}