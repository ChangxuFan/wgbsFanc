# da <- readRDS("sth/da/da1.2_TpTc41_rep13/bulk.list.Rds")
# da <- deseq2.summary(da)
# names(da$MG$summary)
# DMR.to.da <- function(dmr, 
#                       slot.only = F, slot = "summary", root.name = "DMR",
#                       padj.cutoff = NULL, p.cutoff = NULL, 
#                       log2fc.cutoff = NULL # note that diff.Methy is used as "fake" log2fc.cutoff
#                       ) {
# 
#   up.genes <- dmr[dmr$diff.Methy > 0,]
#   summ <- list(
#     padj.cutoff = padj.cutoff,
#     p.cutoff = p.cutoff, 
#     log2fc.cutoff = log2fc.cutoff,
#     n.de = dmr
#   )
#   da <- list(
#     root.name = root.name,
#     
#   )
# }