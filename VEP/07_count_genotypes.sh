#!/bin/sh
#$ -N wsh_07_counts
#$ -cwd
#$ -m beas
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 2
#$ -o 07_counts.$JOB_ID.out
#$ -e 07_counts.$JOB_ID.err

# ===========================================================================
#  PER-VARIANT GENOTYPE COUNTS FOR THE GRADED FILTER
#
#  Streams the cohort VCF once and writes a compact per-variant summary.
#  Nothing is held in memory, so file size is irrelevant.
#
#  The sample-to-group mapping is built AT RUNTIME from `bcftools query -l`,
#  not hard-coded by column number.  If the sample order in the VCF ever
#  differs from the vault record, this script uses the VCF and says so.
#
#  Outputs
#    graded_full_distribution.txt   every (case_homalt, ctl_homalt) cell -> tiny
#    graded_candidates.tsv          rows with >= MIN_CASE hom-alt cases and
#                                   no hom-alt control -> the R input
#    graded_counts_ALL.tsv.gz       every variant, if KEEP_ALL=yes
#
#  NOTE ON THE MODEL
#  The disease is HYPOTHESISED recessive, not established.  This script only
#  counts genotypes; it makes no claim about mode of inheritance.  Every
#  threshold is a parameter, and the full distribution is written out so the
#  strict tier can be seen in the context of everything below it.
# ===========================================================================

MIN_CASE=6        # write detailed rows from this many hom-alt cases upward
MAX_CTL=2         # allow up to this many hom-alt controls in the detailed rows.
                  # NOT the analysis threshold - it is deliberately loose so R
                  # can apply the sex-limited rule later.  Under sex-limited
                  # inheritance a hom-alt MALE control is a healthy carrier and
                  # must NOT exclude a variant; the sexes of WHD01/WHD03/WHD08
                  # are still outstanding with Melissa, so keep the rows.

. /etc/profile.d/modules.sh
module load igmm/apps/bcftools/1.20
command -v bcftools >/dev/null || { echo "bcftools NOT on PATH" >&2; exit 1; }

# Prefer the annotated main-chromosome file, fall back sensibly.
for f in ./cohort_vep107_mainchr.vcf.gz \
         ./cattle_cohort_genotypes.mainchr.vcf.gz \
         ./cattle_cohort_genotypes.renamed.vcf.gz; do
  [ -f "$f" ] && { IN=$f; break; }
done
[ -n "$IN" ] || { echo "No cohort VCF found in $PWD" >&2; exit 1; }

CASES="300042,401868,601898,WHD04,WHD05,WHD06,WHD07,WHD09,WHD10,WHD11"
CONTROLS="WHD01,WHD03,WHD08"

SAMPLES=$(bcftools query -l "$IN" | tr '\n' ',' | sed 's/,$//')

echo "input    : $IN"
echo "samples  : $SAMPLES"
echo "cases    : $CASES"
echo "controls : $CONTROLS"
echo "min_case : $MIN_CASE"
echo "max_ctl  : $MAX_CTL"
echo

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$IN" | \
awk -F'\t' -v OFS='\t' \
    -v samples="$SAMPLES" -v cases="$CASES" -v controls="$CONTROLS" \
    -v mincase="$MIN_CASE" -v maxctl="$MAX_CTL" -v distfile="graded_full_distribution.txt" '
# --- classify one GT string -----------------------------------------------
#   R = hom ref, H = hom alt (any non-ref allele, doubled), T = het, M = missing
function cls(g,   a, n) {
  if (g == "" || g == "." || g == "./." || g == ".|.") return "M"
  n = split(g, a, /[\/|]/)
  if (n != 2) return "M"
  if (a[1] == "." || a[2] == ".") return "M"
  if (a[1] == a[2]) return (a[1] == "0") ? "R" : "H"
  return "T"
}
BEGIN {
  ns = split(samples,  sn, ",")
  nc = split(cases,    cn, ",")
  nk = split(controls, kn, ",")
  for (i = 1; i <= nc; i++) iscase[cn[i]] = 1
  for (i = 1; i <= nk; i++) isctl[kn[i]]  = 1

  ncase = 0; nctl = 0
  for (i = 1; i <= ns; i++) {
    if (sn[i] in iscase) { ncase++; cidx[ncase] = i; cnm[ncase] = sn[i] }
    else if (sn[i] in isctl) { nctl++; kidx[nctl] = i; knm[nctl] = sn[i] }
    else printf("WARNING: sample %s in VCF is in neither group\n", sn[i]) > "/dev/stderr"
  }
  printf("resolved %d cases, %d controls from the VCF header\n", ncase, nctl) > "/dev/stderr"
  if (ncase != nc || nctl != nk)
    printf("WARNING: expected %d cases and %d controls\n", nc, nk) > "/dev/stderr"

  print "CHROM","POS","REF","ALT","case_homalt","case_het","case_homref","case_miss", \
        "ctl_homalt","ctl_het","ctl_homref","ctl_miss","cases_not_homalt"
}
{
  ch = 0; ct = 0; cr = 0; cm = 0
  kh = 0; kt = 0; kr = 0; km = 0
  nothom = ""

  for (i = 1; i <= ncase; i++) {
    c = cls($(4 + cidx[i]))
    if      (c == "H") ch++
    else if (c == "T") { ct++; nothom = nothom (nothom ? "," : "") cnm[i] }
    else if (c == "R") { cr++; nothom = nothom (nothom ? "," : "") cnm[i] }
    else               { cm++; nothom = nothom (nothom ? "," : "") cnm[i] }
  }
  for (i = 1; i <= nctl; i++) {
    c = cls($(4 + kidx[i]))
    if      (c == "H") kh++
    else if (c == "T") kt++
    else if (c == "R") kr++
    else               km++
  }

  # full joint distribution, held in memory as a tiny table
  dist[ch "\t" kh]++

  # detailed rows only where they matter
  if (ch >= mincase && kh <= maxctl)
    print $1, $2, $3, $4, ch, ct, cr, cm, kh, kt, kr, km, nothom
}
END {
  print "case_homalt\tctl_homalt\tn_variants" > distfile
  for (k in dist) print k "\t" dist[k] > distfile
  close(distfile)
}
' > graded_candidates.tsv

echo
echo "=== full joint distribution (case_homalt x ctl_homalt) ==="
sort -k1,1n -k2,2n graded_full_distribution.txt | head -60
echo

echo "=== candidates written ==="
echo "rows (excl. header): $(( $(wc -l < graded_candidates.tsv) - 1 ))"
ls -lh graded_candidates.tsv graded_full_distribution.txt
echo
echo "Copy both to your PC and run graded_filter.R"
