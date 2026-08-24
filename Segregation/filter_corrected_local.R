# ###########################################################################
#
#   CORRECTED SEGREGATION ANALYSIS  -  runs locally, no Eddie needed
#
#   Input:  graded_candidates.tsv   (the file already on your PC)
#   Run:    Rscript filter_corrected_local.R
#
#   ---------------------------------------------------------------------
#   WHY THE COHORT IS FOUR GROUPS, NOT TWO
#   ---------------------------------------------------------------------
#   Melissa's sample table showed WHD04, WHD05 and WHD08 are MAL,E. The
#   phenotype is female-limited, so the ten-cases-plus-three-controls design
#   does not survive that.
#
#     AFFECTED FEMALES (8)  300042 401868 601898 WHD06 WHD07 WHD09 WHD10 WHD11
#         must be HOMOZYGOUS. They have the disease, so under a recessive
#         model both inherited copies carry the variant.
#
#     SIRES (2)  WHD04, WHD05
#         must carry AT LEAST ONE copy. Both are recorded as sires of affected
#         heifers. A daughter needs two copies, one from each parent, so the
#         sire supplied one. He may carry Tt to appear in. Demanding exactly
#         one is as wrong as demanding two, and the tWO and stay healthy, having no
#         reproductive tract for the defecwo earlier analyses made
#         one of those mistakes each.
#
#     UNAFFECTED FEMALES (2)  WHD01, WHD03
#         must NOT be homozygous. Only these two can exclude anything.
#
#     UNAFFECTED MALE (1)  WHD08
#         ignored. He looks healthy carrying two copies, one, or none.
#
#   ---------------------------------------------------------------------
#   TWO LIMITATIONS OF WORKING FROM THIS FILE, STATED UP FRONT
#   ---------------------------------------------------------------------
#   graded_candidates.tsv holds group COUNTS, not per-animal genotypes, so:
#
#   1. ctl_homalt counts homozygotes across all three unaffected animals
#      together. A value of 1 could be WHD08 (harmless) or WHD01/WHD03
#      (disqualifying), and this file cannot tell them apart. So the script
#      requires ctl_homalt == 0. That is CONSERVATIVE: it discards real
#      candidates where only WHD08 was homozygous, but never admits a false
#      one. Recovering those needs a re-extraction with per-control columns.
#
#   2. The 8/8 rung is exact. Lower rungs are approximate, because when both
#      an affected female and a sire fail at the same site, the counts cannot
#      say which of them was heterozygous and which was reference.
#
#   Neither limitation affects the strictest rung, which is the one that
#   matters most.
#
#   Numbers are never compared as text. In R "10" < "7" is TRUE, so a quoted
#   threshold silently deletes the best-segregating variants. Every count is
#   coerced to integer on read.
#
# ###########################################################################

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

IN_FILE <- "graded_candidates.tsv"
OUT_DIR <- "."

AFFECTED <- c("300042","401868","601898","WHD06","WHD07","WHD09","WHD10","WHD11")
SIRES    <- c("WHD04","WHD05")
N_AFF    <- length(AFFECTED)     # 8

GAP_BP    <- 5e5                 # variants closer than this are one block
SV_MIN_BP <- 50L

say  <- function(...) cat(sprintf(...), "\n", sep = "")
rule <- function(t) cat("\n", strrep("=", 75), "\n ", t, "\n",
                        strrep("=", 75), "\n", sep = "")

# ===========================================================================
#  STEP 1.  READ
# ===========================================================================
rule("STEP 1  read")
stopifnot(file.exists(IN_FILE))

d <- read_tsv(IN_FILE, na = c("", "NA"), show_col_types = FALSE,
              guess_max = 10000)

cnt <- c("case_homalt","case_het","case_homref","case_miss",
         "ctl_homalt","ctl_het","ctl_homref","ctl_miss")
stopifnot(all(cnt %in% names(d)))

d <- d %>% mutate(across(all_of(cnt), as.integer),
                  POS   = as.integer(POS),
                  CHROM = as.character(CHROM),
                  cases_not_homalt = ifelse(is.na(cases_not_homalt), "",
                                            as.character(cases_not_homalt)))
say("rows read : %s", format(nrow(d), big.mark = ","))

# ===========================================================================
#  STEP 2.  REBUILD THE FOUR GROUPS FROM THE COUNTS
#
#  The old "case" group was the eight affected females PLUS the two sires.
#  cases_not_homalt names whichever of those ten were not homozygous, so the
#  groups can be separated again by reading that column.
# ===========================================================================
rule("STEP 2  rebuild groups")

count_in <- function(x, who) {
  vapply(strsplit(x, ",", fixed = TRUE),
         function(v) sum(trimws(v) %in% who), integer(1))
}

d <- d %>%
  mutate(
    n_aff_fail  = count_in(cases_not_homalt, AFFECTED),
    n_sire_fail = count_in(cases_not_homalt, SIRES),
    aff_homalt  = N_AFF - n_aff_fail,

    #  A failing sire is acceptable only if he is HETEROZYGOUS, not reference
    #  and not missing. When nothing in the whole case group is reference or
    #  missing, every failure must be a heterozygote, so both sires carry at
    #  least one copy. This is exact whenever no affected female also fails.
    no_ref_or_missing = (case_homref == 0L & case_miss == 0L),
    sires_ok = (n_sire_fail == 0L) | no_ref_or_missing,

    #  Conservative: WHD08 cannot be separated from WHD01/WHD03 here.
    controls_ok = (ctl_homalt == 0L)
  )

say("affected females homozygous : %d to %d of %d",
    min(d$aff_homalt), max(d$aff_homalt), N_AFF)
say("rows with both sires carrying a copy : %s",
    format(sum(d$sires_ok), big.mark = ","))
say("rows with no unaffected homozygote   : %s",
    format(sum(d$controls_ok), big.mark = ","))

pool <- d %>% filter(sires_ok, controls_ok)
say("\npool after sire and control conditions : %s",
    format(nrow(pool), big.mark = ","))

# ===========================================================================
#  STEP 3.  CLASS, IMPACT, SIFT  -  described and ranked, never filtered
# ===========================================================================
rule("STEP 3  class and impact")

IMPACT_LEVELS <- c("HIGH","MODERATE","LOW","MODIFIER","unknown")

pool <- pool %>%
  mutate(
    symbolic = startsWith(ALT, "<") | grepl("[][]", ALT),
    dlen = abs(nchar(ALT) - nchar(REF)),
    vclass = case_when(symbolic ~ "symbolic_SV",
                       nchar(REF) == 1L & nchar(ALT) == 1L ~ "SNP",
                       dlen >= SV_MIN_BP ~ "SV",
                       dlen >  0L ~ "indel",
                       dlen == 0L & nchar(REF) > 1L ~ "MNP",
                       TRUE ~ "other"),
    impact_stratum = factor(
      if_else(is.na(impact) | !(toupper(impact) %in% IMPACT_LEVELS),
              "unknown", toupper(impact)), levels = IMPACT_LEVELS),
    impact_rank = case_when(impact_stratum == "HIGH" ~ 4L,
                            impact_stratum == "MODERATE" ~ 3L,
                            impact_stratum == "LOW" ~ 2L,
                            impact_stratum == "MODIFIER" ~ 1L, TRUE ~ 0L),
    # unknown ranks ABOVE tolerated: unscored is not the same as cleared
    sift_class = case_when(
      !is.na(sift_term) & grepl("deleterious", sift_term, TRUE) ~ "deleterious",
      !is.na(sift_term) & grepl("tolerated",   sift_term, TRUE) ~ "tolerated",
      TRUE ~ "unknown"))

pool %>% count(vclass, sort = TRUE) %>% as.data.frame() %>% print(row.names = FALSE)

# ===========================================================================
#  STEP 4.  THE LADDER
#  Only the affected-female requirement moves. Sire and control conditions
#  stay fixed, so any change between rungs is attributable to the case side.
# ===========================================================================
rule("STEP 4  the ladder")

TIERS <- 8:6
ladder <- bind_rows(lapply(TIERS, function(k) {
  s <- pool %>% filter(aff_homalt >= k)
  tibble(rung = sprintf("%d/%d", k, N_AFF), k = k,
         n_exactly = sum(pool$aff_homalt == k),
         n_cumulative = nrow(s),
         n_SNP = sum(s$vclass == "SNP"),
         n_SV  = sum(s$vclass %in% c("SV","symbolic_SV")),
         n_HIGH = sum(s$impact_stratum == "HIGH"),
         n_MODERATE = sum(s$impact_stratum == "MODERATE"),
         n_genes = n_distinct(s$symbol[!is.na(s$symbol)]),
         n_chrom = n_distinct(s$CHROM))
}))
print(as.data.frame(ladder), row.names = FALSE)
write_tsv(ladder, file.path(OUT_DIR, "corrected_ladder.tsv"))
say("\n  note: only the 8/8 rung is exact. See the header for why.")

# ===========================================================================
#  STEP 5.  BLOCKS  -  the analysis that actually decides
#
#  A causal variant sits inside ONE inherited stretch of chromosome, so the
#  survivors should gather into a few tight clusters rather than scatter.
#
#  Read n_variants AND per_100kb together. A block of 300 variants spread over
#  40 Mb is background however large it looks. One of 40 across 90 kb is a
#  haplotype.
# ===========================================================================
rule("STEP 5  candidate blocks at 8/8")

blocks_at <- function(k) {
  s <- pool %>% filter(aff_homalt >= k)
  if (!nrow(s)) return(NULL)
  s %>% arrange(CHROM, POS) %>%
    group_by(CHROM) %>%
    mutate(gap = POS - lag(POS),
           new = is.na(gap) | gap > GAP_BP,
           bid = cumsum(new)) %>%
    group_by(CHROM, bid) %>%
    summarise(n_variants = n(),
              start_Mb = round(min(POS)/1e6, 3),
              end_Mb   = round(max(POS)/1e6, 3),
              span_kb  = round((max(POS)-min(POS))/1e3, 1),
              per_100kb = round(n() / pmax((max(POS)-min(POS))/1e5, 0.01), 1),
              n_HIGH = sum(impact_stratum == "HIGH"),
              n_MODERATE = sum(impact_stratum == "MODERATE"),
              n_SV = sum(vclass %in% c("SV","symbolic_SV")),
              genes = paste(unique(na.omit(symbol)), collapse = ","),
              .groups = "drop") %>%
    mutate(rung = k) %>% arrange(desc(n_variants))
}

b8 <- blocks_at(8L)
if (!is.null(b8) && nrow(b8)) {
  say("blocks found: %d\n", nrow(b8))
  b8 %>% head(20) %>% mutate(genes = substr(genes, 1, 40)) %>%
    select(CHROM, n_variants, start_Mb, end_Mb, span_kb, per_100kb,
           n_HIGH, n_MODERATE, n_SV, genes) %>%
    as.data.frame() %>% print(row.names = FALSE)

  tot <- sum(b8$n_variants)
  say("\n  largest block  : %.1f%% of survivors", 100*b8$n_variants[1]/tot)
  say("  top three      : %.1f%%", 100*sum(head(b8$n_variants,3))/tot)
  say("  singletons     : %d of %d blocks",
      sum(b8$n_variants == 1), nrow(b8))
  write_tsv(bind_rows(lapply(TIERS, blocks_at)),
            file.path(OUT_DIR, "corrected_blocks.tsv"))
} else {
  say("no variants survive at 8/8")
}

# ===========================================================================
#  STEP 6.  RUNG AGAINST IMPACT
#  Segregation is evidence measured in these animals. Impact is a prediction
#  made without seeing them. MODIFIER dominating every row is expected, and is
#  exactly why impact must never be used as a filter.
# ===========================================================================
rule("STEP 6  rung against impact")
grid <- bind_rows(lapply(TIERS, function(k) {
  s <- pool %>% filter(aff_homalt >= k)
  tibble(rung = sprintf("%d/%d", k, N_AFF),
         HIGH = sum(s$impact_stratum == "HIGH"),
         MODERATE = sum(s$impact_stratum == "MODERATE"),
         LOW = sum(s$impact_stratum == "LOW"),
         MODIFIER = sum(s$impact_stratum == "MODIFIER"),
         unknown = sum(s$impact_stratum == "unknown"),
         total = nrow(s))
}))
print(as.data.frame(grid), row.names = FALSE)
write_tsv(grid, file.path(OUT_DIR, "corrected_impact_grid.tsv"))

# ===========================================================================
#  STEP 7.  GENES AND SHORTLIST
#  A gene carrying several independently segregating variants is stronger
#  evidence than any single variant's predicted severity.
# ===========================================================================
rule("STEP 7  genes at 8/8")
genes <- pool %>% filter(aff_homalt >= 8L, !is.na(symbol)) %>%
  group_by(symbol) %>%
  summarise(chrom = first(CHROM), pos_Mb = round(min(POS)/1e6, 2),
            n_variants = n(),
            n_HIGH = sum(impact_stratum == "HIGH"),
            n_MODERATE = sum(impact_stratum == "MODERATE"),
            n_delet = sum(sift_class == "deleterious"),
            .groups = "drop") %>%
  arrange(desc(n_HIGH), desc(n_MODERATE), desc(n_variants))

if (nrow(genes)) {
  print(as.data.frame(head(genes, 30)), row.names = FALSE)
  write_tsv(genes, file.path(OUT_DIR, "corrected_genes.tsv"))
} else {
  say("no named genes. Every survivor is intergenic.")
  say("That is a result, not a failure: a long-range regulatory cause looks")
  say("exactly like this. Take the block coordinates from STEP 5 to Ensembl.")
}

short <- pool %>%
  filter(aff_homalt >= 8L, impact_stratum %in% c("HIGH","MODERATE")) %>%
  arrange(desc(impact_rank), desc(sift_class == "deleterious"), CHROM, POS)

rule("STEP 7b  shortlist: HIGH and MODERATE at 8/8")
if (nrow(short)) {
  short %>% select(any_of(c("CHROM","POS","REF","ALT","vclass","consequence",
                            "impact","symbol","sift_term"))) %>%
    head(60) %>% as.data.frame() %>% print(row.names = FALSE)
  write_tsv(short, file.path(OUT_DIR, "corrected_shortlist.tsv"))
} else {
  say("Empty. No protein-altering variant segregates at 8/8.")
  say("Compatible with a regulatory cause, a structural variant VEP cannot")
  say("score, or one poorly genotyped animal. Read the blocks instead.")
}

# ===========================================================================
#  STEP 8.  FILES AND PLOTS
# ===========================================================================
rule("STEP 8  files and plots")
for (k in TIERS) {
  s <- pool %>% filter(aff_homalt >= k) %>% arrange(CHROM, POS)
  write_tsv(s, file.path(OUT_DIR, sprintf("corrected_tier%d.tsv", k)))
  say("corrected_tier%d.tsv : %s rows", k, format(nrow(s), big.mark = ","))
}

theme_set(theme_minimal(base_size = 11))
chr_order <- c(as.character(1:29), "X", "MT")
top <- pool %>% filter(aff_homalt >= 8L) %>%
  mutate(CHROM = factor(CHROM, levels = chr_order)) %>% filter(!is.na(CHROM))

pl <- list()
if (nrow(top)) {
  pl[[1]] <- ggplot(top, aes(POS/1e6, CHROM)) +
    geom_point(alpha = 0.45, size = 1.1, colour = "#2c6e91") +
    scale_y_discrete(limits = rev(levels(top$CHROM))) +
    labs(title = "Surviving variants at 8/8",
         subtitle = "An inherited block appears as a tight vertical stack",
         x = "Position (Mb)", y = "Chromosome")
}
pl[[length(pl)+1]] <- ggplot(ladder, aes(factor(k, levels = TIERS), n_cumulative)) +
  geom_col(fill = "#2c6e91", width = 0.7) +
  geom_text(aes(label = format(n_cumulative, big.mark = ",")),
            vjust = -0.4, size = 3.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Survivors at each rung",
       subtitle = "Sire and control conditions held constant",
       x = "Affected females required homozygous", y = "Variants")

if (!is.null(b8) && nrow(b8)) {
  bb <- b8 %>% head(12) %>%
    mutate(lab = sprintf("chr%s %.2f-%.2f Mb", CHROM, start_Mb, end_Mb))
  pl[[length(pl)+1]] <- ggplot(bb, aes(reorder(lab, n_variants), n_variants)) +
    geom_col(fill = "#7a9e7e", width = 0.7) +
    geom_text(aes(label = sprintf("%.0f/100kb", per_100kb)),
              hjust = -0.15, size = 3) +
    coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(title = "Candidate blocks, ranked",
         subtitle = "Label is density. High count with low density is background.",
         x = NULL, y = "Variants in block")
}

pdf(file.path(OUT_DIR, "corrected_plots.pdf"), width = 9, height = 6)
invisible(lapply(pl, print)); invisible(dev.off())
say("wrote corrected_plots.pdf")

# ===========================================================================
rule("HOW TO READ THIS")
cat("
  1. corrected_blocks.tsv   ranked candidate regions.
     Read n_variants AND per_100kb together. Tight and dense is a haplotype.
     A big count spread over megabases is background, however big it looks.

  2. corrected_genes.tsv, corrected_shortlist.tsv
     What sits in those regions. Where a pattern becomes a hypothesis.

  3. corrected_ladder.tsv
     Whether relaxing from 8/8 recovers anything or only adds noise.

  Several blocks is the EXPECTED result. Eight related heifers share several
  stretches by descent and most cause nothing. The counts cannot rank them.
  The genes inside them can.

  Recessive inheritance is assumed here, not demonstrated. Only WHD01 and
  WHD03 can exclude anything. WHD08 was set aside because a healthy male is
  uninformative for a female-limited trait, not because his data were poor.
")
