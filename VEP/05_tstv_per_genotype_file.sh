#!/bin/sh
#$ -N wsh_05_tstv_files
#$ -cwd
#$ -m beas
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -o 05_tstv_files.$JOB_ID.out
#$ -e 05_tstv_files.$JOB_ID.err

# ===========================================================================
#  PER-SAMPLE Ts/Tv FROM THE INDIVIDUAL PANGENIE GENOTYPE FILES
#
#  Reads every VCF in pangenie/work/genotypes/ rather than the merged cohort
#  file, so each animal is measured on its own PanGenie output with no merge
#  step in between.
#
#  Two things this handles that a naive loop does not:
#
#  1. PanGenie writes a genotype for EVERY panel site, including 0/0.  Running
#     plain `bcftools stats` would return the panel's Ts/Tv, identical for all
#     animals.  `-s -` emits the PSC block, which counts transitions and
#     transversions from the sample's actual non-reference genotypes.
#
#  2. These files predate the chromosome rename, so contigs are NC_037328.1
#     style.  The autosome list is auto-detected from the first contig, and
#     applied with -t (targets), which streams and needs no index.
#
#  Output: tstv_per_genotype_file.tsv
# ===========================================================================

. /etc/profile.d/modules.sh
module load igmm/apps/bcftools/1.20
command -v bcftools >/dev/null || { echo "bcftools NOT on PATH" >&2; exit 1; }

GTDIR=/exports/eddie/scratch/$USER/Dissertation/pangenie/work/genotypes
TSV=./tstv_per_genotype_file.tsv
RAW=./tstv_raw_psc.txt

[ -d "$GTDIR" ] || { echo "MISSING $GTDIR" >&2; exit 1; }

FILES=$(ls "$GTDIR"/*.vcf.gz "$GTDIR"/*.vcf 2>/dev/null)
[ -n "$FILES" ] || { echo "No VCF files in $GTDIR" >&2; ls -la "$GTDIR" >&2; exit 1; }

echo "directory : $GTDIR"
echo "files     : $(echo "$FILES" | wc -l)"
echo

# --- work out the autosome naming from the first file ----------------------
FIRST=$(echo "$FILES" | head -1)
FIRSTCHR=$(bcftools view -h "$FIRST" | awk -F'[=,]' '/^##contig/ {print $3; exit}')
echo "first contig in $(basename "$FIRST"): $FIRSTCHR"

case "$FIRSTCHR" in
  NC_*)
    REG=$(for i in $(seq 28 56); do printf "NC_0373%02d.1," "$i"; done | sed 's/,$//')
    echo "naming    : RefSeq accessions -> autosomes NC_037328.1 .. NC_037356.1"
    ;;
  *)
    REG=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29
    echo "naming    : plain -> autosomes 1 .. 29"
    ;;
esac
echo

: > "$RAW"
{
printf "sample\tgroup\tfile\tnRefHom\tnNonRefHom\tnHets\tnTs\tnTv\tTsTv\thet_hom\tnIndels\tnMissing\n"

for f in $FILES; do
  SM=$(bcftools query -l "$f" | head -1)
  [ -n "$SM" ] || { echo "  skipped (no sample column): $(basename "$f")" >&2; continue; }
  echo "  processing $SM  <- $(basename "$f")" >&2

  # -t (targets) streams and needs no index, unlike -r (regions).
  bcftools stats -s - -t "$REG" "$f" 2>/dev/null | grep '^PSC' >> "$RAW"

  bcftools stats -s - -t "$REG" "$f" 2>/dev/null | grep '^PSC' | \
  awk -F'\t' -v sm="$SM" -v fn="$(basename "$f")" '
  BEGIN {
    split("300042 401868 601898 WHD04 WHD05 WHD06 WHD07 WHD09 WHD10 WHD11", c, " ")
    for (i in c) grp[c[i]] = "case"
    split("WHD01 WHD03 WHD08", k, " ")
    for (i in k) grp[k[i]] = "control"
  }
  {
    g      = (sm in grp) ? grp[sm] : "UNKNOWN"
    tstv   = ($8 > 0) ? $7 / $8 : 0
    hethom = ($5 > 0) ? $6 / $5 : 0
    printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\n", \
           sm, g, fn, $4, $5, $6, $7, $8, tstv, hethom, $9, $14
  }'
done
} > "$TSV"

echo
echo "=== per-sample table ==="
column -t -s "$(printf '\t')" "$TSV"
echo

echo "=== group means ==="
awk -F'\t' 'NR>1 && $2!="UNKNOWN" {n[$2]++; ts[$2]+=$9; hh[$2]+=$10}
END {for (g in n) printf "  %-8s n=%d  mean Ts/Tv=%.4f  mean het/hom=%.4f\n", g, n[g], ts[g]/n[g], hh[g]/n[g]}' "$TSV"
echo

echo "=== outliers (Ts/Tv more than 2 SD from the cohort mean) ==="
awk -F'\t' 'NR>1 {v[NR]=$9; s+=$9; n++; nm[NR]=$1; gr[NR]=$2}
END {
  if (n < 2) { print "  too few samples"; exit }
  m = s/n
  for (i in v) ss += (v[i]-m)^2
  sd = sqrt(ss/(n-1))
  printf "  cohort mean %.4f, SD %.4f, n=%d\n", m, sd, n
  f = 0
  for (i in v) if (sd > 0 && ((v[i]-m)/sd > 2 || (m-v[i])/sd > 2)) {
    printf "  OUTLIER  %-10s %-8s Ts/Tv=%.4f  (%+.2f SD)\n", nm[i], gr[i], v[i], (v[i]-m)/sd
    f = 1
  }
  if (!f) print "  none"
}' "$TSV"
echo

# Any sample present in the files but not in the case/control key
awk -F'\t' 'NR>1 && $2=="UNKNOWN" {print "  UNRECOGNISED SAMPLE: "$1"  (file "$3")"}' "$TSV"

echo
echo "STEP 5 OK"
ls -lh "$TSV" "$RAW"
