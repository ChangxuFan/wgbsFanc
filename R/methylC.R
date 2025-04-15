methylC.import <- function(methylC, threads = 4, region = NULL) {
  if (is.null(names(methylC))) {
    names(methylC) <- basename(methylC) %>% sub("(\\.CG)*.methylC.*$", "", .)
  }
  grs <- utilsFanc::safelapply(methylC, function(me) {
    gr <- rtracklayer::import.bed(me, which = region)
    mcols(gr) <- mcols(gr) %>% as.data.frame() %>% 
      dplyr::rename(CG = score) %>% 
      dplyr::mutate(depth = thick.start - 1 ) %>% 
      dplyr::select(CG, depth)
    return(gr)
  }, threads = threads)
  return(grs)
}

minmax.dmr <- function(mc.grs, mc.df = NULL, samples,
                       root.name = NULL, work.dir, 
                       depth_min.filter = c(4, 25),
                       depth_cats = c(4, 10, 15, 20, 25),
                       CG_diff_cats = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                       print.tracks = T) {
  if (is.null(root.name)) {
    root.name <- basename(work.dir)
  }
  # if (!is.null())
  mc.dfs <- lapply(samples, function(sample) {
    mc.gr <- mc.grs[[sample]]
    df <- mcols(mc.gr) %>% as.data.frame() %>% 
      dplyr::rename(CG = score) %>% 
      dplyr::mutate(depth = thick.start - 1 ) %>% 
      dplyr::select(CG, depth) %>% `colnames<-`(., paste0(colnames(.), "_", sample))
    df$locus <- utilsFanc::gr.get.loci(mc.gr)
    return(df)
  })
  
  lapply(mc.dfs, head)
  mc.df <- Reduce(inner_join, mc.dfs)
  mc.df$depth_min <- mc.df[, paste0("depth_", samples)] %>% as.matrix() %>% rowMin()
  mc.df$CG_min <- mc.df[, paste0("CG_", samples)] %>% as.matrix() %>% rowMin()
  mc.df$CG_max <- mc.df[, paste0("CG_", samples)] %>% as.matrix() %>% rowMax()
  mc.df$CG_diff <- mc.df$CG_max - mc.df$CG_min 
  mc.df <- mc.df %>% dplyr::filter(depth_min >= depth_min.filter[1], depth_min <= depth_min.filter[2])
  mc.df$depth_cat <- cut(mc.df$depth_min, breaks = depth_cats, 
                         include.lowest = T, ordered_result = T, right = T)
  mc.df$CG_cat <- cut(mc.df$CG_diff, breaks = CG_diff_cats, 
                      include.lowest = T, ordered_result = T, right = T)
  mc.df <- mc.df %>% na.omit()
  pl <- list()
  pl$raw <- ggplot(mc.df, aes(x = depth_cat, fill = CG_cat)) + geom_bar()
  mc.df.pct <- mc.df %>% group_by(CG_cat, depth_cat) %>% summarise(n = n()) %>% 
    ungroup() %>% group_by(depth_cat) %>% dplyr::mutate(pct = n/sum(n)) %>% as.data.frame()
  
  pl$pct <- ggplot(mc.df.pct, aes(x = depth_cat, y = pct, fill = CG_cat)) + geom_bar(stat = "identity")
  
  trash <- scFanc::wrap.plots.fanc(pl, plot.out = paste0(work.dir, "/", root.name, "_minmax_stats.png"))
  
  if (print.tracks) {
    trash <- mc.df %>% split(., f = .$depth_cat) %>% lapply(function(mc.by.depth) {
      mc.by.depth %>% split(., f = .$CG_cat) %>% 
        lapply(function(mc.by.CG) {
          CG.cat <- mc.by.CG$CG_cat[1] %>% gsub("\\(|\\]|\\)|\\[", "", .) %>% gsub(",", "_", .)
          depth.cat <- mc.by.CG$depth_cat[1] %>% gsub("\\(|\\]|\\)|\\[", "", .) %>% gsub(",", "_", .)
          out.file <- paste0(work.dir, "/tracks/", root.name, "_depth", depth.cat, "_CG", CG.cat, ".methylC")
          out.df <- mc.by.CG %>% utilsFanc::loci.2.df(loci.col.name = "locus", remove.loci.col = T) %>% 
            dplyr::mutate(name = "CG", strand = "+") %>% 
            dplyr::select(chr, start, end, name, CG_diff, strand, depth_min) %>% 
            dplyr::arrange(chr, start)
          utilsFanc::write.zip.fanc(out.df, out.file = out.file, bed.shift = T, ez = T)
        })
    })
  }
  # saveRDS(mc.df, paste0(work.dir, "/"))
  mc.o <- list(mc.df = mc.df, mc.df.pct = mc.df.pct, 
               samples = samples,
              work.dir = work.dir, root.name = root.name)
  return(mc.o)
}

mc.plot <- function(mc.o, return.pl = F) {
  pl <- list()
  pl$raw <- ggplot(mc.o$mc.df, aes(x = depth_cat, fill = CG_cat)) + geom_bar()
  mc.o$mc.df.pct <- mc.o$mc.df %>% group_by(CG_cat, depth_cat) %>% summarise(n = n()) %>% 
    ungroup() %>% group_by(depth_cat) %>% dplyr::mutate(pct = n/sum(n)) %>% as.data.frame()
  
  pl$pct <- ggplot(mc.o$mc.df.pct, aes(x = depth_cat, y = pct, fill = CG_cat)) + geom_bar(stat = "identity")
  if (!is.null(mc.o$work.dir)) {
    plot.out <- paste0(mc.o$work.dir, "/", mc.o$root.name, "_minmax_stats.png")
  } else {
    plot.out <- NULL
  }
  if (return.pl) {
    return(pl)
  }
  p <- scFanc::wrap.plots.fanc(pl, plot.out = plot.out)
  invisible(p)
}

mc.smoothing <- function(mc.o, window = 2, new.wd = NULL) {
  mc.o$mc.df$CG_diff <- smooth.by.shift(mc.o$mc.df$CG_diff, win = window)
  mc.o$work.dir <- new.wd
  if (!is.null(new.wd))
    mc.o$root.name <- basename(new.wd)
  return(mc.o)
}

mc.smooth.plot <- function(mc.o, windows, plot.out = NULL) {
  pl <- lapply(windows, function(window) {
    mc.o.smooth <- mc.smoothing(mc.o, window = window)
    # browser()
    pl <- mc.plot(mc.o = mc.o.smooth, return.pl = T) %>% 
      lapply(function(p) {
        return(p + ggtitle(paste0("window: ", window)))
      }) 
    return(pl)
  }) %>% Reduce(c, .)
  p <- scFanc::wrap.plots.fanc(plot.list = pl, plot.out = plot.out, n.split = 2)
  return(p)
}

smooth.by.shift <- function(vec, win) {
  shift.vec <- c(rev(-1 * (1:win)), 1:win)
  n <- length(vec)
  vec.out <- lapply(shift.vec, function(shift) {
    if (shift > 0) {
      vec.shifted <- vec[(shift + 1):n] %>% c(rep(vec[n], shift))
    } else {
      shift <- abs(shift)
      vec.shifted <- vec[1:(n-shift)] %>% c(rep(vec[1], shift), .)
    }
  })  %>% as.data.frame() %>% as.matrix() %>% rowMeans()
  return(vec.out)
}

mc.pca <- function(mc.o, bin.size = NULL, bin.use.genome,
                   ntop = 1000000, depth.cutoff = 5, depth.cap = 30,
                   samples = NULL,
                   rm.sex.chr = T,
                   suffix = NULL) {
  if (is.null(samples))
    samples <- mc.o$samples
  
  df <- mc.o$mc.df 
  
  if (!is.null(bin.size)) {
    browser()
  }
  
  df <- df %>% dplyr::filter(depth_min >= depth.cutoff, depth_max <= depth.cap)
  
  if (rm.sex.chr) df <- df %>% dplyr::filter(!grepl("chrX|chrY", locus))
  
  df <- df %>% .[, c("locus", paste0("CG_", samples))]

  
  if (nrow(df) == 0) {
    stop("nrow(df) == 0")
  }
  
  df$locus <- NULL
  rv <- rowVars(df %>% as.matrix())
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  n.total <- length(select)
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  pc.combinations <- list(
    c(1, 2),
    c(1, 3),
    c(1, 4),
    c(1, 5)
  )
  pl <- lapply(pc.combinations, function(pc.combo) {
    if (max(pc.combo) > ncol(pca$x)) return()
    d <- data.frame(PC1 = pca$x[, pc.combo[1]], PC2 = pca$x[, pc.combo[2]], sample = rownames(pca$x))
    pc.x <- paste0("PC", pc.combo[1])
    pc.y <- paste0("PC", pc.combo[2])
    colnames(d) <- c(pc.x, pc.y, "sample")
    
    p <- ggplot(data = d, aes_string(x = pc.x, y = pc.y, color = "sample", label = "sample")) + 
      geom_point(size = 3) + xlab(paste0(pc.x, ": ", round(percentVar[pc.combo[1]] * 100), "% variance")) + 
      ylab(paste0(pc.y, ": ", round(percentVar[pc.combo[2]] * 100), "% variance")) + coord_fixed() +
      ggrepel::geom_text_repel() +
      ggtitle(paste0("n = ", n.total, "; depth ", depth.cutoff))
  })
  
  trash <- scFanc::wrap.plots.fanc(plot.list = pl, sub.width = 8, sub.height = 8,
                          plot.out = paste0(mc.o$work.dir, "/", mc.o$root.name, "_pca_depth", depth.cutoff,
                                            ifelse(rm.sex.chr, "_rmSex_", "_" ),
                                            "ntop", utilsFanc::so.formatter(ntop), ".png"))
  invisible(pl)
}

mc.cor <- function(mc.o) {
  cor.p <- cor(mc.o$mc.df[, paste0("CG_", mc.o$samples)]) %>% ComplexHeatmap::Heatmap()
  scFanc::save.base.plot(cor.p, file = paste0(mc.o$work.dir, "/", mc.o$root.name, "_cor.png"))
  return()
}

CG.density <- function(gr, genome) {
  stop("never finished coding because nCG column is output by DSS in callDMR()")
  CG.bed <- paste0("~/genomes/", genome, "/cpg/", genome, "_CpG.bed.gz")
  GC.gr <- rtracklayer::import.bed(CG.bed, which = gr)
  t <- plyranges::join_overlap_intersect()
  o <- findOverlaps(query = GC.gr, subject = gr, ignore.strand = T) %>% as.data.frame()
  o <- o %>% group_by(subjectHits) %>% summarise(n = n()) %>% ungroup() %>% 
    as.data.frame()
  
  browser()
  print("miao")
}

