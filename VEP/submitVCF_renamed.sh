#!/bin/sh
#$ -N vep_wsh
#$ -cwd
#$ -m beas
#$ -l h_rt=48:00:00
#$ -l h_vmem=32G
#$ -pe sharedmem 16

. /etc/profile.d/modules.sh
module load igmm/apps/vep/107
module load igmm/apps/bcftools

# ---------------------------------------------------------------------------
# Paths.  Only the input VCF and the cache are absolute; chr_map.txt is read
# from, and all output written to, whichever directory you submit from.
# So this script works from any folder, as long as chr_map.txt sits beside it.
# ---------------------------------------------------------------------------
BASE=/exports/eddie/scratch/$USER/Dissertation
VCF=$BASE/pangenie/work/cattle_cohort_genotypes.vcf.gz
CACHE=$BASE/vep/vep_cache

MAP=./chr_map.txt
RENAMED=./cattle_cohort_genotypes.renamed.vcf.gz

echo "working directory : $PWD"
echo "input VCF         : $VCF"
echo "cache             : $CACHE"
echo

# --- check everything exists before burning 45 minutes --------------------
[ -f "$VCF"   ] || { echo "MISSING input VCF: $VCF" >&2; exit 1; }
[ -f "$MAP"   ] || { echo "MISSING chr_map.txt in $PWD - copy it next to this script" >&2; exit 1; }
[ -d "$CACHE" ] || { echo "MISSING cache dir: $CACHE" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 1. Rename chromosomes on the INPUT, so the VEP summary HTML is readable
#    NC_037328.1 -> 1 ... NC_037357.1 -> X, NC_006853.1 -> MT
# ---------------------------------------------------------------------------
echo "chromosomes BEFORE:"
bcftools index --stats $VCF | cut -f1 | head -5

bcftools annotate --rename-chrs $MAP $VCF -Oz -o $RENAMED --threads 4
bcftools index -t $RENAMED

echo "chromosomes AFTER:"
bcftools index --stats $RENAMED | cut -f1 | head -35

# bcftools does NOT error when the map matches nothing, so check explicitly.
if bcftools index --stats $RENAMED | cut -f1 | grep -q '^NC_037'; then
  echo "ERROR: NC_037 names still present - chr_map.txt did not match the VCF." >&2
  exit 1
fi
echo "OK: chromosomes renamed."
echo

# ---------------------------------------------------------------------------
# 2. Annotate the renamed VCF
# ---------------------------------------------------------------------------
vep -i $RENAMED --cache --offline --dir_cache $CACHE --cache_version 107 \
  --merged --species bos_taurus --format vcf --vcf \
  --sift b --symbol --regulatory --pick --fork 16 --force_overwrite \
  -o cohort_vep107_renamed.vcf --stats_file cohort_vep107_renamed_summary.html

bgzip -f cohort_vep107_renamed.vcf
tabix -p vcf cohort_vep107_renamed.vcf.gz

echo
ls -lh cohort_vep107_renamed.vcf.gz* cohort_vep107_renamed_summary.html
