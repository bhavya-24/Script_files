#!/bin/sh
#$ -N wsh_01_rename
#$ -cwd
#$ -m beas
#$ -l h_rt=08:00:00
#$ -l h_vmem=4G
#$ -pe sharedmem 4
#$ -o 01_rename.$JOB_ID.out
#$ -e 01_rename.$JOB_ID.err

# ===========================================================================
#  STEP 1 of 3 - RENAME CHROMOSOMES  (bcftools)
#
#    NC_037328.1 -> 1  ...  NC_037356.1 -> 29
#    NC_037357.1 -> X
#    NC_006853.1 -> MT
#
#  Needs chr_map.txt in the submit directory.
#  Output: cattle_cohort_genotypes.renamed.vcf.gz (+ .tbi)
#
#  Next: 02_subset_vcftools.sh
# ===========================================================================

. /etc/profile.d/modules.sh
module load igmm/apps/bcftools/1.20

command -v bcftools >/dev/null || { echo "bcftools NOT on PATH" >&2; exit 1; }
bcftools --version | head -1
echo

BASE=/exports/eddie/scratch/$USER/Dissertation
VCF=$BASE/pangenie/work/cattle_cohort_genotypes.vcf.gz

MAP=./chr_map.txt
OUT=./cattle_cohort_genotypes.renamed.vcf.gz

[ -f "$VCF" ] || { echo "MISSING input VCF: $VCF" >&2; exit 1; }
[ -f "$MAP" ] || { echo "MISSING chr_map.txt in $PWD" >&2; exit 1; }

echo "in  : $VCF"
echo "map : $MAP"
echo "out : $OUT"
echo

echo "--- chromosome names BEFORE ---"
bcftools index --stats "$VCF" | cut -f1 | head -5
echo "(total contigs: $(bcftools index --stats "$VCF" | wc -l))"
echo

bcftools annotate --rename-chrs "$MAP" "$VCF" -Oz -o "$OUT" --threads 4
bcftools index -t "$OUT"

echo "--- chromosome names AFTER ---"
bcftools index --stats "$OUT" | cut -f1 | head -35
echo "(total contigs: $(bcftools index --stats "$OUT" | wc -l))"
echo

# bcftools does NOT error when the map matches nothing, so check explicitly.
if bcftools index --stats "$OUT" | cut -f1 | grep -q '^NC_037'; then
  echo "ERROR: NC_037 names still present - chr_map.txt did not match the VCF." >&2
  echo "Compare the map against: bcftools index --stats $VCF | cut -f1" >&2
  exit 1
fi

echo "variants: $(bcftools index --nrecords "$OUT")"
echo
echo "STEP 1 OK - chromosomes renamed."
ls -lh "$OUT"*
echo
echo "Next: qsub 02_subset_vcftools.sh"
