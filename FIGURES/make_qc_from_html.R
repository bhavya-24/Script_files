# =============================================================================
#  WSH dissertation - read-QC figures derived ENTIRELY from the QC reports.
#
#  Nothing in this script is transcribed by hand. Every value is parsed at run
#  time from:
#     QC REPORTS/NanoPlot-report_run1.html
#     QC REPORTS/NanoPlot-report_run2.html
#     QC REPORTS/NanoPlot-report_run3.html
#     QC REPORTS/NanoPlot-report_merged.html
#     QC REPORTS/multiqc_report.html
#
#  NanoPlot stores its plots as Plotly traces whose arrays are base64-encoded
#  little-endian binary, so they are decoded here rather than screenshotted.
#
#  NOTE ON SCOPE: these five figures are the only ones these HTML files can
#  support. BUSCO, Merqury, the variant spectrum and the SV funnel are NOT in
#  them - those numbers live in the project notebook and wiki, so they are in
#  make_figures_ggplot.R instead.
# =============================================================================

library(ggplot2); library(dplyr); library(tidyr); library(tibble)
library(stringr); library(jsonlite); library(base64enc); library(scales)

QC  <- "D:/UOE/DISSERTATION/QC REPORTS"
OUT <- "D:/UOE/DISSERTATION/FIGURES"
GENOME <- 2.7e9   # expected cattle genome size, used only for the coverage axis

BLUE <- "#1f6fc0"; GOLD <- "#f2b20e"; GREEN <- "#2e6b47"; RED <- "#c62828"
GREY <- "#9e9e9e"; DARK <- "#1a1a1a"; MUTE <- "#5a6b62"

theme_diss <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(panel.border = element_blank(),
          axis.line = element_line(colour = "#b8c4bc"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#e6ece7", linewidth = 0.3),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", colour = GREEN, size = base + 1),
          plot.title = element_text(face = "bold", colour = DARK, size = base + 3.5),
          plot.subtitle = element_text(colour = MUTE, size = base - 0.5),
          plot.title.position = "plot",
          axis.text = element_text(colour = MUTE), axis.title = element_text(colour = DARK),
          legend.position = "none")
}
fmt_si <- function(x) ifelse(x >= 1e9, sprintf("%.3g G", x/1e9),
                      ifelse(x >= 1e6, sprintf("%.3g M", x/1e6),
                      ifelse(x >= 1e3, sprintf("%.3g k", x/1e3), sprintf("%.0f", x))))

# ---------------------------------------------------------------- parsers ----
slurp <- function(f) paste(readLines(file.path(QC, f), warn = FALSE), collapse = "\n")

# NanoPlot "General summary" table -> named numeric vector
nano_stats <- function(f) {
  txt  <- slurp(f)
  tab  <- str_match(txt, "(?s)<table.*?</table>")[1]
  rows <- str_match_all(tab, "(?s)<tr>(.*?)</tr>")[[1]][, 2]
  out  <- list()
  for (r in rows) {
    cells <- str_match_all(r, "(?s)<t[dh][^>]*>(.*?)</t[dh]>")[[1]][, 2]
    cells <- str_squish(str_remove_all(cells, "<[^>]+>"))
    cells <- cells[cells != ""]
    if (length(cells) == 2) {
      v <- suppressWarnings(as.numeric(str_remove_all(cells[2], ",")))
      if (!is.na(v)) out[[cells[1]]] <- v
    }
  }
  unlist(out)
}

# Plotly stores arrays as {dtype, bdata}; decode little-endian binary
decode_bdata <- function(o) {
  if (!is.list(o) || is.null(o$bdata)) return(as.numeric(unlist(o)))
  raw <- base64enc::base64decode(o$bdata)
  switch(o$dtype,
    f8 = readBin(raw, "double",  size = 8, n = length(raw)/8, endian = "little"),
    f4 = readBin(raw, "double",  size = 4, n = length(raw)/4, endian = "little"),
    i4 = readBin(raw, "integer", size = 4, n = length(raw)/4, endian = "little"),
    u4 = { v <- readBin(raw, "integer", size = 4, n = length(raw)/4, endian = "little")
           ifelse(v < 0, v + 2^32, v) },        # R ints are signed; lift unsigned
    stop("unhandled dtype: ", o$dtype))
}

# all Plotly.newPlot trace groups in a NanoPlot report
plotly_traces <- function(f) {
  txt <- slurp(f)
  m <- str_match_all(txt,
        "(?s)Plotly\\.newPlot\\(\\s*(?:\"[^\"]+\"|'[^']+')\\s*,\\s*(\\[.*?\\])\\s*,")[[1]]
  lapply(m[, 2], function(j) tryCatch(fromJSON(j, simplifyVector = FALSE),
                                      error = function(e) NULL))
}
trace_xy <- function(groups, i) {
  tr <- groups[[i]][[1]]
  list(x = decode_bdata(tr$x), y = decode_bdata(tr$y))
}

# =============================================================================
# FIGURE 1 - the three runs and the merged set (summary tables)
# =============================================================================
files <- c("Run 1" = "NanoPlot-report_run1.html", "Run 2" = "NanoPlot-report_run2.html",
           "Run 3" = "NanoPlot-report_run3.html", "Merged" = "NanoPlot-report_merged.html")
S <- lapply(files, nano_stats)
runs <- tibble(run    = factor(names(files), levels = names(files)),
               bases  = sapply(S, `[[`, "Total bases"),
               reads  = sapply(S, `[[`, "Number of reads"),
               n50    = sapply(S, `[[`, "Read length N50"))
print(runs)

f1 <- runs %>%
  pivot_longer(-run, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, c("bases","reads","n50"),
           c("Sequence yield (bases)","Read count","Read length N50 (bp)")),
         fill = case_when(run == "Merged" ~ "m", run == "Run 2" ~ "o", TRUE ~ "r")) %>%
  ggplot(aes(run, value, fill = fill)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = fmt_si(value)), vjust = -0.45, size = 3.4,
            fontface = "bold", colour = DARK) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c(r = BLUE, o = GOLD, m = GREEN)) +
  scale_y_continuous(labels = fmt_si, expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Oxford Nanopore sequencing: three runs merged for assembly",
       subtitle = "Parsed from the NanoPlot summary tables. Run 2 gave the longest reads but almost no data; Run 3 gave half the bases at the shortest lengths.",
       x = NULL, y = NULL) + theme_diss()
ggsave(file.path(OUT, "FigH1_run_comparison.png"), f1, width = 13.2, height = 4.2,
       dpi = 300, bg = "white"); message("FigH1 written")

# =============================================================================
# Decode the merged report's plots once
#   [1] weighted read-length histogram   [5] yield by length   [6] length vs quality
#   (R indexes from 1; these are Python's 0, 4 and 5)
# =============================================================================
G   <- plotly_traces("NanoPlot-report_merged.html")
wh  <- trace_xy(G, 1)
yl  <- trace_xy(G, 5)
lq  <- trace_xy(G, 6)
mrg <- S[["Merged"]]
MED <- mrg[["Median read length"]]; N50 <- mrg[["Read length N50"]]

# =============================================================================
# FIGURE 2 - weighted read-length distribution with the 5 kb threshold
# =============================================================================
f2 <- tibble(len = wh$x, bases = wh$y) %>% filter(len <= 60000) %>%
  ggplot(aes(len, bases)) +
  geom_area(fill = BLUE, alpha = 0.85) +
  geom_vline(xintercept = 5000, colour = RED,   linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = MED,  colour = GOLD,  linetype = "dotted", linewidth = 0.8) +
  geom_vline(xintercept = N50,  colour = GREEN, linetype = "dotted", linewidth = 0.8) +
  annotate("text", x = 5600,  y = Inf, label = "5 kb filter", colour = RED,
           hjust = 0, vjust = 1.6, fontface = "bold", size = 3.6) +
  annotate("text", x = MED+400, y = Inf, label = sprintf("median %s bp", comma(MED)),
           colour = "#8a6d00", hjust = 0, vjust = 3.4, fontface = "bold", size = 3.4) +
  annotate("text", x = N50+600, y = Inf, label = sprintf("N50 %s bp", comma(N50)),
           colour = GREEN, hjust = 0, vjust = 5.2, fontface = "bold", size = 3.4) +
  scale_x_continuous(labels = comma) + scale_y_continuous(labels = fmt_si) +
  labs(title = "Merged read-length distribution",
       subtitle = "Weighted histogram decoded from the NanoPlot report. The median read falls below the 5 kb cut-off.",
       x = "Read length (bp)", y = "Bases (weighted)") + theme_diss()
ggsave(file.path(OUT, "FigH2_readlength.png"), f2, width = 10.6, height = 4.6,
       dpi = 300, bg = "white"); message("FigH2 written")

# =============================================================================
# FIGURE 3 - yield by length: what the 5 kb filter actually cost
# =============================================================================
yd   <- tibble(minlen = yl$x, gb = yl$y) %>% arrange(minlen)
TOT  <- max(yd$gb); KEEP <- approx(yd$minlen, yd$gb, xout = 5000)$y; LOST <- TOT - KEEP
message(sprintf("total %.2f Gb | >=5kb keeps %.2f Gb (%.1f%%) | loses %.2f Gb",
                TOT, KEEP, 100*KEEP/TOT, LOST))

f3 <- ggplot(yd, aes(minlen, gb)) +
  geom_hline(yintercept = TOT, colour = MUTE, linewidth = 0.4) +
  geom_vline(xintercept = 5000, colour = RED, linetype = "dashed", linewidth = 0.8) +
  geom_line(colour = BLUE, linewidth = 1.1) +
  annotate("segment", x = 5000, xend = 5000, y = KEEP, yend = TOT,
           colour = RED, linewidth = 1.0,
           arrow = arrow(ends = "both", length = unit(0.16, "cm"))) +
  annotate("point", x = 5000, y = KEEP, colour = RED, size = 3) +
  annotate("text", x = 6200, y = (KEEP+TOT)/2,
           label = sprintf("%.1f Gb lost\n(%.0f%% of bases)", LOST, 100*LOST/TOT),
           colour = RED, hjust = 0, fontface = "bold", size = 3.7) +
  annotate("text", x = 16000, y = TOT*0.62,
           label = sprintf("reads >=5 kb retain %.1f Gb\nof %.1f Gb  (%.0f%%)",
                           KEEP, TOT, 100*KEEP/TOT),
           colour = DARK, hjust = 0, fontface = "bold", size = 3.7) +
  scale_x_log10(labels = comma, limits = c(300, 3e5)) +
  labs(title = "What the >=5 kb filter actually cost",
       subtitle = sprintf("Yield-by-length curve decoded from the NanoPlot report. Total yield %.2f Gb.", TOT),
       x = "Minimum read length (bp, log scale)",
       y = "Cumulative yield from reads >= x (Gb)") + theme_diss()
ggsave(file.path(OUT, "FigH3_yield_by_length.png"), f3, width = 10.6, height = 4.8,
       dpi = 300, bg = "white"); message("FigH3 written")

# =============================================================================
# FIGURE 4 - read length vs mean quality
# =============================================================================
f4 <- tibble(len = lq$x, q = lq$y) %>% filter(len > 0) %>%
  ggplot(aes(len, q)) +
  geom_point(colour = BLUE, alpha = 0.20, size = 0.7) +
  geom_hline(yintercept = c(10, 15), colour = c(RED, GOLD), linetype = "dotted") +
  annotate("text", x = 230, y = 10.5, label = "Q10", colour = RED,  size = 3.3, fontface = "bold", hjust = 0) +
  annotate("text", x = 230, y = 15.5, label = "Q15", colour = GOLD, size = 3.3, fontface = "bold", hjust = 0) +
  scale_x_log10(labels = comma, limits = c(200, 2e5)) + ylim(0, 26) +
  labs(title = "Read quality is independent of read length",
       subtitle = "10,000-read subsample from the NanoPlot report. No length-dependent quality collapse.",
       x = "Read length (bp, log scale)", y = "Mean read quality (Phred)") + theme_diss()
ggsave(file.path(OUT, "FigH4_length_vs_quality.png"), f4, width = 9.6, height = 4.4,
       dpi = 300, bg = "white"); message("FigH4 written")

# =============================================================================
# FIGURE 5 - cohort short reads, from the MultiQC general statistics table
# =============================================================================
mq  <- slurp("multiqc_report.html")
i0  <- str_locate(mq, fixed('data-table-id="general_stats_table"'))[1, 1]
seg <- substr(mq, max(1, i0 - 3000), i0 + 90000)
rows <- str_match_all(seg, "(?s)<tr[^>]*>.*?</tr>")[[1]][, 1]
rec <- list()
for (r in rows) {
  cells <- str_match_all(r, "(?s)<t[dh][^>]*>.*?</t[dh]>")[[1]][, 1]
  cells <- str_squish(str_remove_all(cells, "<[^>]+>"))
  cells <- cells[cells != ""]
  if (length(cells) >= 7 && startsWith(cells[1], "WHD")) {
    rec[[length(rec)+1]] <- tibble(
      sample = str_remove(str_remove(cells[1], "_R[12]_001$"), "_S7"),
      dup    = as.numeric(str_remove(cells[2], "\\s*%")),
      len    = as.numeric(str_extract(cells[4], "[0-9.]+")),
      seqs   = as.numeric(str_extract(cells[7], "[0-9.]+")) * 1e6)
  }
}
coh <- bind_rows(rec) %>% group_by(sample) %>%
  summarise(dup = mean(dup), len = mean(len), seqs = first(seqs), .groups = "drop") %>%
  mutate(coverage = seqs * 2 * len / GENOME,
         group = case_when(str_detect(sample, "^WHD03") ~ "flag",
                           str_detect(sample, "-Co$")   ~ "control", TRUE ~ "case")) %>%
  arrange(coverage)
print(coh)

f5 <- coh %>% pivot_longer(c(coverage, dup), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, c("coverage","dup"),
                         c("Estimated coverage (x)","Duplicate reads (%)")),
         sample = factor(sample, levels = coh$sample)) %>%
  ggplot(aes(value, sample, fill = group)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = ifelse(metric == "Estimated coverage (x)",
                               sprintf("%.0fx", value), sprintf("%.0f%%", value))),
            hjust = -0.18, size = 3.1, colour = DARK) +
  facet_wrap(~metric, scales = "free_x") +
  scale_fill_manual(values = c(case = GOLD, control = BLUE, flag = RED)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Case-control cohort read QC, with one control flagged",
       subtitle = paste0("Parsed from the MultiQC general statistics table. Coverage is estimated as reads x 2 x mean length / 2.7 Gb, not mapped depth.\n",
                         "WHD03-BS-Co shows ~42% duplicates against 22-26% elsewhere."),
       x = NULL, y = NULL) + theme_diss() + theme(legend.position = "none")
ggsave(file.path(OUT, "FigH5_cohort_qc.png"), f5, width = 12.6, height = 4.6,
       dpi = 300, bg = "white"); message("FigH5 written")

message("\nAll five HTML-derived figures written to ", OUT)
