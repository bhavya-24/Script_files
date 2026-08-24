#!/bin/sh
#$ -N wsh_04_tstv
#$ -cwd
#$ -m beas
#$ -l h_rt=08:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -o 04_tstv.$JOB_ID.out
#$ -e 04_tstv.$JOB_ID.err

# ===========================================================================
#  PER-SAMPLE Ts/Tv ON THE PANGENIE-GENOTYPED COHORT
#
#  Answers James's point from the Stage 2 QC: the graph VCF had no usable
#  per-sample GT, so per-sample Ts/Tv had to wait until after genotyping.
#  It has.
#
#  Uses `bcftools stats -s -`, which emits a PSC (per-sample counts) block
#  in ONE pass over the file.  The earlier failed approach looped
#  `bcftools view -s SAMPLE` per sample, which keeps every site including
#  0/0, so all 13 samples returned the identical cohort-wide 2.16.
#
#  AUTOSOMES ONLY (1-29).  X is excluded because males and females carry
#  different copy numbers, and MT because mitochondrial Ts/Tv is far higher
#  than nuclear and there are only ~62 sites - either would distort a
#  between-sample comparison.
#
#  Output: tstv_per_sample.tsv   (tab-separated, ready for R)
#          cohort_stats_autosomes.txt   (full bcftools stats dump)
# ===========================================================================

. /etc/profile.d/modules.sh
module load igmm/apps/bcftools/1.20
command -v bcftools >/dev/null || { echo "bcftools NOT on PATH" >&2; exit 1; }

# Use the main-chromosome file if step 2 has finished, else the renamed one.
if [ -f ./cattle_cohort_genotypes.mainchr.vcf.gz ]; then
  IN=./cattle_cohort_genotypes.mainchr.vcf.gz
elif [ -f ./cattle_cohort_genotypes.renamed.vcf.gz ]; then
  IN=./cattle_cohort_genotypes.renamed.vcf.gz
else
  echo "No renamed or mainchr VCF found in $PWD - run step 1 first." >&2
  exit 1
fi

AUTOSOMES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29
STATS=cohort_stats_autosomes.txt
TSV=tstv_per_sample.tsv

echo "input     : $IN"
echo "regions   : autosomes 1-29"
echo "samples   : $(bcftools query -l "$IN" | tr '\n' ' ')"
echo

echo "=== running bcftools stats (one pass, all samples) ==="
bcftools stats -s - -r "$AUTOSOMES" "$IN" > "$STATS"
echo "wrote $STATS"
echo

# --- cohort-wide, for reference -------------------------------------------
echo "=== cohort-wide Ts/Tv (all samples pooled) ==="
grep '^TSTV' "$STATS" | head -1 | awk -F'\t' '{printf "  ts=%s  tv=%s  Ts/Tv=%s\n", $3, $4, $5}'
echo

# ---------------------------------------------------------------------------
#  PSC columns (bcftools 1.20):
#   [3]sample [4]nRefHom [5]nNonRefHom [6]nHets [7]nTransitions
#   [8]nTransversions [9]nIndels [10]avg_depth [11]nSingletons
#   [12]nHapRef [13]nHapAlt [14]nMissing
# ---------------------------------------------------------------------------
{
printf "sample\tgroup\tnRefHom\tnNonRefHom\tnHets\tnTs\tnTv\tTsTv\thet_hom\tnIndels\tnMissing\n"
grep '^PSC' "$STATS" | awk -F'\t' '
BEGIN {
  split("300042 401868 601898 WHD04 WHD05 WHD06 WHD07 WHD09 WHD10 WHD11", c, " ")
  for (i in c) grp[c[i]] = "case"
  split("WHD01 WHD03 WHD08", k, " ")
  for (i in k) grp[k[i]] = "control"
}
{
  s = $3
  g = (s in grp) ? grp[s] : "UNKNOWN"
  tstv    = ($8 > 0) ? $7 / $8 : 0
  hethom  = ($5 > 0) ? $6 / $5 : 0
  printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\n", \
         s, g, $4, $5, $6, $7, $8, tstv, hethom, $9, $14
}'
} > "$TSV"

echo "=== per-sample table ==="
column -t -s "$(printf '\t')" "$TSV"
echo
echo "wrote $TSV"
echo

# --- group summary and outlier flag ---------------------------------------
echo "=== group means ==="
awk -F'\t' 'NR>1 {n[$2]++; ts[$2]+=$8; hh[$2]+=$9}
END {for (g in n) printf "  %-8s n=%d  mean Ts/Tv=%.4f  mean het/hom=%.4f\n", g, n[g], ts[g]/n[g], hh[g]/n[g]}' "$TSV"
echo

echo "=== outliers (Ts/Tv more than 2 SD from the cohort mean) ==="
awk -F'\t' 'NR>1 {v[NR]=$8; s+=$8; n++; name[NR]=$1; grp[NR]=$2}
END {
  m = s/n
  for (i in v) ss += (v[i]-m)^2
  sd = (n>1) ? sqrt(ss/(n-1)) : 0
  printf "  cohort mean %.4f, SD %.4f\n", m, sd
  flagged = 0
  for (i in v) if (sd > 0 && (v[i]-m)/sd > 2 || sd > 0 && (m-v[i])/sd > 2) {
    printf "  OUTLIER  %-10s %-8s Ts/Tv=%.4f  (%+.2f SD)\n", name[i], grp[i], v[i], (v[i]-m)/sd
    flagged = 1
  }
  if (!flagged) print "  none"
}' "$TSV"
echo

echo "STEP 4 OK"
ls -lh "$TSV" "$STATS"
