# =============================================================================
#  Dot plot, assembly against ARS-UCD1.2
#
#  Same plot as before. Two changes only:
#    1. reference sequences relabelled from RefSeq accessions to 1..29 and X,
#       matching the naming used everywhere else after the rename step in 2.10
#    2. chromosomes ordered by number instead of by size, and the contig
#       labels on the x axis turned off because they overlap
#
#  Neither change touches the filters or drops any alignment.
# =============================================================================

library(pafr)
library(ggplot2)

ali <- read_paf("wsh_vs_ref.paf")
cat("rows read                :", nrow(ali), "\n")

# ---- relabel reference sequences ---------------------------------------
acc <- sprintf("NC_0373%02d.1", 28:57)      # chromosomes 1..29 then X
lab <- c(as.character(1:29), "X")
names(lab) <- acc
hit <- ali$tname %in% acc
ali$tname[hit] <- lab[as.character(ali$tname[hit])]
cat("rows relabelled          :", sum(hit), "\n")

# ---- original filters, unchanged ---------------------------------------
ali <- subset(ali, alen > 500000)
ali <- subset(ali, tlen > 20000000)
ali <- subset(ali, qlen > 5000000)
cat("rows after filters       :", nrow(ali), "\n")

# ---- orderings ----------------------------------------------------------
# contigs keep the size ordering pafr used before
q_ord <- unique(ali$qname[order(-ali$qlen)])
# chromosomes go in number order
t_ord <- intersect(c(as.character(1:29), "X"), unique(ali$tname))

cat("contigs plotted          :", length(q_ord), "\n")
cat("chromosomes plotted      :", length(t_ord), "\n")

p <- dotplot(ali, order_by = "provided", ordering = list(q_ord, t_ord),
             label_seqs = TRUE, dashes = TRUE,
             xlab = "Whitebred Shorthorn assembly contigs",
             ylab = "ARS-UCD1.2 reference chromosomes (Mb)")

# pafr draws the sequence names as text inside the panel, not as axis text,
# so theme() cannot hide them. Drop the text layer that carries the contig
# names and keep the one carrying the chromosome names.
txt <- which(vapply(p$layers, function(L) inherits(L$geom, "GeomText"), logical(1)))
nrows <- vapply(p$layers[txt], function(L) nrow(L$data), numeric(1))
p$layers[[txt[which.max(nrows)]]] <- NULL
cat("dropped contig-label layer of", max(nrows), "labels
")

p <- p + theme(panel.grid.minor = element_blank(),
               panel.grid.major = element_line(color = "grey90", linewidth = 0.2))

ggsave("wsh_vs_ref_named.png", p, width = 14, height = 12, dpi = 150)
cat("written: wsh_vs_ref_named.png\n")
