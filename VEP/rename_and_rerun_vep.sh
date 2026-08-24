#!/bin/bash
#$ -N vep_renamed
#$ -cwd
#$ -pe sharedmem 12
#$ -l h_vmem=8G
#$ -l h_rt=06:00:00
#$ -o vep_renamed.$JOB_ID.out
#$ -e vep_renamed.$JOB_ID.err
# Add your project code if your group requires one:
#   #$ -P <project_code>
# =============================================================================
#  Rename chromosomes NC_037328.1 -> 1 ... X, MT, then re-run VEP so the
#  summary HTML is readable (James's request, 2026-07-28).
#
#  Renames the INPUT VCF, not the annotated output, because the summary HTML
#  is written during the VEP run.
# =============================================================================
set -euo pipefail

BASE=/exports/eddie/scratch/$USER/Dissertation
IN=$BASE/pangenie/work/cattle_cohort_genotypes.vcf.gz
MAP=$BASE/vep/chr_map.txt
CACHE=$BASE/vep/vep_cache
RENAMED=$BASE/vep/cattle_cohort_genotypes.renamed.vcf.gz

module load igmm/apps/bcftools
module load igmm/apps/vep/107

echo "=== 1. checking inputs ==="
for f in "$IN" "$MAP"; do
  [ -f "$f" ] || { echo "MISSING: $f" >&2; exit 1; }
done
echo "chromosome names BEFORE rename:"
bcftools index --stats "$IN" | cut -f1 | head -5
echo

echo "=== 2. renaming chromosomes ==="
bcftools annotate --rename-chrs "$MAP" "$IN" -Oz -o "$RENAMED" --threads 4
bcftools index -t "$RENAMED"

echo "chromosome names AFTER rename:"
bcftools index --stats "$RENAMED" | cut -f1 | head -35
echo

# Fail loudly if the rename did nothing (wrong accessions in the map file)
if bcftools index --stats "$RENAMED" | cut -f1 | grep -q '^NC_037'; then
  echo "!!! ERROR: NC_037 names still present. chr_map.txt did not match." >&2
  echo "!!! Compare the map against: bcftools index --stats $IN | cut -f1" >&2
  exit 1
fi
echo "OK: main chromosomes renamed."
echo

echo "=== 3. re-running VEP ==="
cd $BASE/vep
vep -i "$RENAMED" \
    --cache --offline --dir_cache "$CACHE" --cache_version 107 \
    --merged --species bos_taurus \
    --format vcf --vcf \
    --sift b --symbol --regulatory --pick \
    --fork 12 --force_overwrite \
    -o cohort_vep107_renamed.vcf \
    --stats_file cohort_vep107_renamed_summary.html

echo
echo "=== 4. done ==="
ls -lh cohort_vep107_renamed.vcf cohort_vep107_renamed_summary.html
echo "compress + index the output:"
bgzip -f cohort_vep107_renamed.vcf
tabix -p vcf cohort_vep107_renamed.vcf.gz
ls -lh cohort_vep107_renamed.vcf.gz*
