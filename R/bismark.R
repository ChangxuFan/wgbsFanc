bismark.xiaoyu.wrapper <- function(fastqs, suffix = "_R[12].fq.gz", out.dir = "./bismark/",
                                   non.directional = F,
                                   genome, ref.dir = NULL, shift = 0, shift.chr,
                                   min.insert = 0, max.insert = 2000,
                                   parallel.bismark = 2, parallel.bowtie2 = 6, parallel.sample = 1,
                                   add.params = "",
                                   bismark = "/opt/apps/bismark/0.23.1/bismark",
                                   tempt.dir = "/scratch/fanc/tmp/", run = T, sort.only = F,
                                   samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  # written according to Xiaoyu Zhuo's Snakefile
  if (is.null(ref.dir))
    ref.dir <- paste0("~/genomes/", genome, "/bismark")
  if (!dir.exists(ref.dir))
    stop(!dir.exists(ref.dir))
  dir.create(out.dir, showWarnings = F, recursive = T)
  samples <- fastqs %>% basename() %>% sub(suffix, "", .)
  sample.table <- table(samples)
  not.pair <- names(sample.table)[sample.table != 2]
  if (length(not.pair) > 0) {
    stop("some fastqs are not properly paired: \n" %>% 
           paste0(not.pair, collapse = "\n"))
  }
  bams <- fastqs %>% split(f = samples) %>% 
    utilsFanc::safelapply(function(fastq.pair) {
      fastq.pair <- sort(fastq.pair)
      sample <- fastq.pair[1] %>% basename() %>% sub(suffix, "", .)
      if (sort.only == F) {
        log.bismark <- paste0(out.dir, "/", sample, ".log")
        if (non.directional) {
          add.params <- paste0(add.params, " --non_directional")
        }
        cmd <- paste0(bismark, " ", add.params, " -q -I ", min.insert, " -X ", max.insert, 
                      " --parallel ", parallel.bismark, " -p ", parallel.bowtie2, 
                      " --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 ", 
                      " -o ", out.dir, " --temp_dir ", tempt.dir,
                      " --gzip --nucleotide_coverage ", 
                      ref.dir, " -1 ", fastq.pair[1], " -2 ", fastq.pair[2],
                      " 1>", log.bismark, " 2>&1")
        print(cmd)
        if (run == T) {
          system(cmd)
        }
      }
  
      bam <- Sys.glob(paste0(out.dir, "/", sample, "*_bismark_bt2_pe.bam"))
      if (shift !=0 ) {
        bam <- bamFanc::bam.shift.file()
      }
      bam.sorted <- utilsFanc::insert.name.before.ext(name = bam, insert = "sorted", delim = ".")
      cmd <- paste0(samtools, " sort -O bam -m 3G -@ ", parallel.bismark * parallel.bowtie2,
                    " ", bam, " > ", bam.sorted,
                    " && ", samtools, " index ", bam.sorted)
      print(cmd)
      if (run == T) {
        system(cmd)
      }
      bam <- bam.sorted
      if (run == T) {
        liteRnaSeqFanc::bam.browser(bam = bam, thread = parallel.bismark * parallel.bowtie2,
                                    normalization = "None")
      }
      return(bam)
    }, threads = parallel.sample)
  
  return(bams)
}

bismark.filter.noMismatch <- function(bams, regions = NULL, out.dir = NULL, 
                                      out.bams = NULL, threads = 1, index = T,
                                      mode = "noMismatch",
                                      samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  if (is.null(out.bams)) {
    if (is.null(out.dir))
      out.dir <- bams %>% dirname() %>% unique()
    if (length(out.dir) != 1) {
      stop("length(out.dir) != 1")
    }
    basenames <- basename(bams) %>% utilsFanc::insert.name.before.ext(insert = mode, delim = "_") 
    if (! is.null(regions)) {
      basenames <- basenames %>% utilsFanc::insert.name.before.ext(insert = paste0(regions, collapse = "_"), delim = "_") 
    }
    out.bams <- paste0(out.dir, "/", basenames) 
  }
  res <- utilsFanc::safelapply(seq_along(bams), function(i) {
    bam <- bams[i]
    out.bam <- out.bams[i]
    if (is.null(regions)) {
      regions <- ""
    } else {
      regions <- paste0(regions, collapse = " ")
    }
    if (mode == "noMismatch") {
      cmd <- paste0(samtools, " view -h ", bam, " ", regions,  " | ",
                    "awk -F \"\\t\" 'BEGIN{OFS = \"\\t\"} $1 ~ /@/ || ($13 !~ /[AT\\^]/ && $6 !~ /[DI]/ && $5 > 1) {print $0}' | ",
                    samtools, " view -hbo ", out.bam, " -")
      
    } else if (mode == "noIndel") {
      cmd <- paste0(samtools, " view -h ", bam, " ", regions,  " | ",
                    "awk -F \"\\t\" 'BEGIN{OFS = \"\\t\"} $1 ~ /@/ || ($6 !~ /[DI]/ && $5 > 1) {print $0}' | ",
                    samtools, " view -hbo ", out.bam, " -")
    }
    cat(cmd)
    system(cmd)
    if (index == T) {
      cmd <- paste0(samtools, " index ", out.bam)
      cat(cmd)
      system(cmd)
    }
    return(out.bam)
  }, threads = threads)
  return(res)
}

wgbs.bam.cov <- function(bams, smart.base = T, genome, regions, out.dir) {
  if (is.null(names(bams))) {
    names(bams) <- basename(bams)
    if (smart.base) {
      names(bams) <- names(bams) %>% sub("_R1.*$|_bismark.*$", "", .)
    }
  }
  cpgs.gr <- wgbs.get.cpgs(genome = genome, regions = regions)
  lapply(names(bams), function(sample) {
    bam <- bams[sample]
    param <- ScanBamParam(which=cpgs.gr)
    bam.gr <- GenomicAlignments::readGAlignments(file = bam, param = param,
                                                 use.names = T, with.which_label = T) %>% 
      as("GRanges")
    df <- data.frame(read.name = names(bam.gr), cpg = bam.gr$which_label) %>% unique()
    df.cov <- df %>% dplyr::group_by(which_label) %>% dplyr::summarise(cov = n())
    cov.track <- df %>% utilsFanc::loci.2.df(loci.col.name = "cpg", remove.loci.col = T, 
                                             return.gr = T)
    browser()
    print("m")
    print("m")
    print("m")
  })
}

wgbs.get.cpgs <- function(genome, regions) {
  if (is.character(regions)) {
    regions <- utilsFanc::import.bed.fanc(bed = regions, return.gr = T)
  }
  cpgs.bed <- paste0("~/genomes/", genome, "/cpg/", genome, "_CpG.bed.gz")
  cpgs.gr <- lapply(1:length(regions), function(i) {
    region <- regions[i]
    cpgs.gr <- rtracklayer::import.bed(con = cpgs.bed, which = region)
    cpgs.gr$name <- region$forth
    cpgs.gr$name_full <- paste0(cpgs.gr$name, "@", utilsFanc::gr.get.loci(cpgs.gr))
    return(cpgs.gr)
  }) %>% Reduce(c, .) %>% utilsFanc::gr.drop.seqlevel()
  return(cpgs.gr)
}

