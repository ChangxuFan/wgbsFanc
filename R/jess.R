dmr.mat.gen <- function(bw.vec, dmr.bed, threads = 1) {
  
  dmr.gr <- rtracklayer::import.bed(dmr.bed)
  dmr.gr$dmr.locus <- paste0(seqnames(dmr.gr), ":", start(dmr.gr), "-", end(dmr.gr))
  df <- utilsFanc::safelapply(names(bw.vec), function(sample) {
    bw <- bw.vec[[sample]]
    data.in.dmr <- rtracklayer::import(con = bw, which = dmr.gr)
    j <- plyranges::join_overlap_left(data.in.dmr, dmr.gr)
    
    df <- mcols(j) %>% as.data.frame() %>% 
      dplyr::group_by(dmr.locus) %>% dplyr::summarise(score= mean(score)) %>% 
      dplyr::ungroup() %>% as.data.frame()
    colnames(df) <- c("dmr.locus", sample)
    return(df)
  }, threads = threads) %>% Reduce(left_join, .)
  df <- df %>% na.omit()
  browser()
  rownames(df) <- df$dmr.locus
  df$dmr.locus <- NULL
  mat <- as.matrix(df)
  return(mat)
}