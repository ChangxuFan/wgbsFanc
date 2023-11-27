bismark.parse.bam <- function(bams, genome, region.gr) {
  methyl.table <- lapply(bams, function(bam) {
    browser()
    GenomicAlignments::readGAlignmentPairs(bam, use.names = T)
  })
}