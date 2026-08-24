# ###########################################################################
#
#   SCORE EVERY BLOCK, NOT JUST THE BIG ONES
#
#       Rscript score_all_blocks.R
#
#   ---------------------------------------------------------------------
#   WHY THIS EXISTS
#   ---------------------------------------------------------------------
#   Blocks were first ranked by variant count, and the regions inspected by
#   hand were the largest ones. That ranking turned out to be the wrong one.
#
#   The region that survived every test, the SAMD12 promoter, was ninth by
#   count. It won on criteria that were only applied later:
#
#     it barely grew when the threshold was relaxed
#     every one of its variants survived the allele-frequency filter
#     it contained ten variants in a promoter, including a clean 59 bp deletion
#     its variant density was normal rather than elevated
#
#   Selecting regions by count and then discovering the real criteria is
#   backwards, and it means a small block with the same profile could have
#   been missed. This script computes the criteria for EVERY block so the
#   choice of what to inspect is made on the right basis.
#
#   Density needs the full call set and so cannot be computed here. It is
#   checked on Eddie for the top few blocks only.
#
#   Needs corrected_tier8.tsv, corrected_tier7.tsv, corrected_tier6.tsv and
#   tier8_not_nearfixed.tsv in this folder.
#
# ###########################################################################

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

GAP <- 5e5        # same block definition used throughout
MIN_N <- 3        # ignore blocks smaller than this; singletons are not blocks

say  <- function(...) cat(sprintf(...), "\n", sep = "")
rule <- function(t) cat("\n", strrep("=", 100), "\n ", t, "\n", strrep("=", 100), "\n", sep = "")

rd <- function(f) {
  if (!file.exists(f)) stop("missing ", f)
  read_tsv(f, show_col_types = FALSE, na = c("", ".", "NA"), guess_max = 5000) %>%
    mutate(CHROM = as.character(CHROM), POS = as.integer(POS))
}
t8   <- rd("corrected_tier8.tsv")
t7   <- rd("corrected_tier7.tsv")
t6   <- rd("corrected_tier6.tsv")
keep <- rd("tier8_not_nearfixed.tsv")

kept_key <- paste(keep$CHROM, keep$POS)

# ---------------------------------------------------------------------------
#  Define blocks on the 8/8 set, then measure each one
# ---------------------------------------------------------------------------
blocks <- t8 %>% arrange(CHROM, POS) %>% group_by(CHROM) %>%
  mutate(gap = POS - lag(POS), new = is.na(gap) | gap > GAP, bid = cumsum(new)) %>%
  ungroup() %>%
  mutate(
    #  A clean structural variant is one where the length changes by 50 bp or
    #  more without the sequence being a simple repeat. Repeat-length changes
    #  are the class every genotyper handles worst, so they are counted apart.
    dlen    = abs(nchar(ALT) - nchar(REF)),
    isSV    = startsWith(ALT, "<") | grepl("[][]", ALT) | dlen >= 50,
    repeaty = grepl("(AT){4,}|(TA){4,}|(CA){4,}|(AC){4,}|(GT){4,}|(TG){4,}|A{8,}|T{8,}|C{8,}|G{8,}",
                    paste0(REF, ALT)),
    cleanSV = isSV & !repeaty,
    #  Positions where a variant could plausibly act
    functional = impact %in% c("HIGH","MODERATE") |
                 consequence %in% c("upstream_gene_variant","downstream_gene_variant") |
                 grepl("UTR|splice|missense", consequence) | cleanSV,
    promoter   = consequence == "upstream_gene_variant",
    kept       = paste(CHROM, POS) %in% kept_key
  )

count_in <- function(d, chrom, s, e)
  sum(d$CHROM == chrom & d$POS >= s & d$POS <= e)

bs <- blocks %>% group_by(CHROM, bid) %>%
  summarise(n8         = n(),
            start      = min(POS),
            end        = max(POS),
            span_kb    = round((max(POS) - min(POS)) / 1e3, 1),
            n_kept     = sum(kept),
            n_func     = sum(functional),
            n_promoter = sum(promoter),
            n_cleanSV  = sum(cleanSV),
            n_coding   = sum(impact %in% c("HIGH","MODERATE")),
            genes      = paste(unique(na.omit(symbol)), collapse = ","),
            .groups = "drop") %>%
  filter(n8 >= MIN_N)

say("blocks with at least %d variants: %d", MIN_N, nrow(bs))

#  Growth when the threshold is relaxed, measured inside the same interval.
#  A block that barely grows is already saturated: everything there that can
#  segregate does. A block that multiplies only looked good because the strict
#  threshold happened to catch a few of many near-miss variants.
bs <- bs %>% rowwise() %>%
  mutate(n6     = count_in(t6, CHROM, start, end),
         growth = round(n6 / n8, 2),
         pct_kept = round(100 * n_kept / n8, 1)) %>%
  ungroup()

# ---------------------------------------------------------------------------
#  Score
#
#  Deliberately simple and additive, so every point is traceable to a stated
#  reason rather than to a tuned weight.
# ---------------------------------------------------------------------------
bs <- bs %>%
  mutate(
    s_growth   = case_when(growth <= 1.3 ~ 3, growth <= 1.7 ~ 2,
                           growth <= 2.5 ~ 1, TRUE ~ 0),
    s_kept     = case_when(pct_kept == 100 ~ 3, pct_kept >= 70 ~ 2,
                           pct_kept >= 30  ~ 1, TRUE ~ 0),
    s_func     = pmin(n_func, 3),
    s_promoter = ifelse(n_promoter >= 3, 2, ifelse(n_promoter >= 1, 1, 0)),
    s_cleanSV  = ifelse(n_cleanSV >= 1, 2, 0),
    s_coding   = ifelse(n_coding  >= 1, 1, 0),
    score      = s_growth + s_kept + s_func + s_promoter + s_cleanSV + s_coding
  ) %>%
  arrange(desc(score), desc(n_func), growth)

rule("ALL BLOCKS, SCORED")
out <- bs %>%
  transmute(chr = CHROM,
            Mb  = sprintf("%.3f-%.3f", start/1e6, end/1e6),
            kb  = span_kb, n8, n6, growth,
            kept = sprintf("%d (%.0f%%)", n_kept, pct_kept),
            func = n_func, prom = n_promoter, SV = n_cleanSV, cod = n_coding,
            score,
            genes = substr(genes, 1, 32))
print(as.data.frame(head(out, 30)), row.names = FALSE)
write_tsv(bs, "block_scores.tsv")

rule("HOW THE SCORE IS BUILT")
cat("
  growth    <=1.3 gives 3, <=1.7 gives 2, <=2.5 gives 1, else 0
            A block already saturated at 8/8 is more convincing than one that
            multiplies as soon as the threshold moves.

  kept      100% of variants surviving the frequency filter gives 3,
            >=70% gives 2, >=30% gives 1, else 0
            An allele common enough that every unaffected animal carries it
            cannot cause a defect confined to a minority of the breed.

  func      one point per functional variant, capped at 3
  prom      3 or more promoter variants gives 2, at least 1 gives 1
  SV        a clean structural variant, repeats excluded, gives 2
  cod       any HIGH or MODERATE impact variant gives 1

  Maximum is 14. Nothing here weights variant count, deliberately: count
  measures how long the shared haplotype is, not how likely it is to matter.
")

rule("NEXT STEP: DENSITY, WHICH CANNOT BE COMPUTED HERE")
cat("
  Density needs the full call set, which is on Eddie. It was the test that
  separated the SAMD12 promoter from everything else:

    chr13:30.406   4.8x denser than flanking sequence   suspect
    chr12 VWA8     4.6x denser                          suspect
    chr14 SAMD12   0.7x, entirely normal                clean

  Run this for the top blocks printed above:
")
top <- head(bs, 8)
for (i in seq_len(nrow(top))) {
  b <- top[i, ]
  cs <- max(1, b$start - 2000); ce <- b$end + 2000
  w  <- ce - cs
  ctrl_s <- b$end + 1000000
  cat(sprintf(
"echo \"chr%s %.3f-%.3f : $(bcftools view -r %s:%d-%d cohort_vep107_mainchr.vcf.gz | grep -vc '^#')  vs control $(bcftools view -r %s:%d-%d cohort_vep107_mainchr.vcf.gz | grep -vc '^#')\"\n",
    b$CHROM, b$start/1e6, b$end/1e6, b$CHROM, cs, ce, b$CHROM, ctrl_s, ctrl_s + w))
}
cat("
  A block whose count is close to its control window is a normal stretch of
  genome in which the disease haplotype happens to sit. A block several times
  denser is a structurally unusual region, and a signal found inside one is
  much easier to explain away.
")
