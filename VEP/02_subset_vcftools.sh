#!/bin/sh
#$ -N wsh_02_subset
#$ -cwd
#$ -m beas
#$ -l h_rt=24:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 1
#$ -o 02_subset.$JOB_ID.out
#$ -e 02_subset.$JOB_ID.err

# ===========================================================================
#  STEP 2 of 3 - KEEP ONLY THE MAIN CHROMOSOMES  (vcftools)
#
#  Drops the NW_ unplaced scaffolds, keeping 1-29, X, MT.
#  Cattle are 2n=60: 29 autosomes + X.  ARS-UCD1.2 has no Y (the reference
#  animal, L1 Dominette 01449, is female), so there is no Y to keep.
#
#  vcftools is single-threaded (hence -pe sharedmem 1) and writes UNCOMPRESSED
#  output, so this is the slow step.  bcftools is loaded only to compress and
#  index the result afterwards.
#
#  Input : cattle_cohort_genotypes.renamed.vcf.gz   (from step 1)
#  Output: cattle_cohort_genotypes.mainchr.vcf.gz   (+ .tbi)
#
#  Next: 03_vep.sh
# ===========================================================================

. /etc/profile.d/modules.sh
module load igmm/apps/bcftools/1.20
command -v bcftools >/dev/null || { echo "bcftools NOT on PATH" >&2; exit 1; }

# --- find the vcftools module rather than hard-coding a version -----------
# Eddie sets no default version for most igmm modules, so a bare
# "module load igmm/apps/vcftools" fails.  Ask the module system what exists
# and load the last (highest) match.
VCFT_MOD=$(module avail 2>&1 | tr ' ' '\n' | grep -iE '/vcftools/' | sort -V | tail -1)

if [ -z "$VCFT_MOD" ]; then
  echo "No vcftools module found on this system." >&2
  echo "Candidates seen:" >&2
  module avail 2>&1 | tr ' ' '\n' | grep -i vcftools >&2
  exit 1
fi

echo "loading vcftools module: $VCFT_MOD"
module load "$VCFT_MOD"

command -v vcftools >/dev/null || {
  echo "vcftools still NOT on PATH after loading $VCFT_MOD" >&2; exit 1; }
vcftools --version 2>&1 | head -1
echo

BASE=/exports/eddie/scratch/$USER/Dissertation
export TMPDIR=$BASE/.tmp
mkdir -p "$TMPDIR"

IN=./cattle_cohort_genotypes.renamed.vcf.gz
PREFIX=./cattle_cohort_genotypes.mainchr
OUT=$PREFIX.vcf.gz

[ -f "$IN" ] || { echo "MISSING $IN - run 01_rename_bcftools.sh first" >&2; exit 1; }

N_BEFORE=$(bcftools index --nrecords "$IN")
echo "in       : $IN"
echo "variants : $N_BEFORE"
echo "free     : $(df -h "$BASE" | tail -1 | awk '{print $4}')"
echo

# --- the subset ------------------------------------------------------------
# --recode-INFO-all is NOT optional: without it vcftools strips the INFO
# column entirely.  Harmless here (annotation has not happened yet) but it
# would delete every CSQ record if this were ever run on a VEP output.
#
# vcftools appends ".recode.vcf" to --out, so the file lands as
# cattle_cohort_genotypes.mainchr.recode.vcf, uncompressed.
vcftools --gzvcf "$IN" \
  --chr 1  --chr 2  --chr 3  --chr 4  --chr 5  --chr 6  --chr 7  --chr 8 \
  --chr 9  --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 \
  --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr 23 --chr 24 \
  --chr 25 --chr 26 --chr 27 --chr 28 --chr 29 \
  --chr X  --chr MT \
  --recode --recode-INFO-all \
  --out "$PREFIX"

[ -f "$PREFIX.recode.vcf" ] || { echo "vcftools produced no output" >&2; exit 1; }
echo
echo "uncompressed size: $(du -h "$PREFIX.recode.vcf" | cut -f1)"

# --- compress + index (bcftools, so bgzip/tabix are not needed) ------------
bcftools view "$PREFIX.recode.vcf" -Oz -o "$OUT" --threads 2
bcftools index -t "$OUT"
rm -f "$PREFIX.recode.vcf"

N_AFTER=$(bcftools index --nrecords "$OUT")
echo
echo "--- contigs kept ---"
bcftools index --stats "$OUT" | cut -f1
echo
echo "variants before : $N_BEFORE"
echo "variants after  : $N_AFTER"
echo "dropped         : $((N_BEFORE - N_AFTER))   <-- quote this number in Methods"
echo
echo "STEP 2 OK."
ls -lh "$OUT"*
echo
echo "Next: qsub 03_vep.sh"
