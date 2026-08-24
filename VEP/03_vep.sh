#!/bin/sh
#$ -N wsh_03_vep
#$ -cwd
#$ -m beas
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 16
#$ -o 03_vep.$JOB_ID.out
#$ -e 03_vep.$JOB_ID.err

# ===========================================================================
#  STEP 3 of 3 - ANNOTATE  (Ensembl VEP 107, merged cache, bos_taurus)
#
#  Input : cattle_cohort_genotypes.mainchr.vcf.gz   (from step 2)
#  Output: cohort_vep107_mainchr.vcf.gz             (+ .tbi)
#          cohort_vep107_mainchr_summary.html
#
#  --fork matches -pe sharedmem.  Keep them equal.
# ===========================================================================

. /etc/profile.d/modules.sh
module load igmm/apps/vep/107
module load igmm/apps/bcftools/1.20

command -v vep      >/dev/null || { echo "vep NOT on PATH" >&2; exit 1; }
command -v bcftools >/dev/null || { echo "bcftools NOT on PATH" >&2; exit 1; }

BASE=/exports/eddie/scratch/$USER/Dissertation
CACHE=$BASE/vep/vep_cache
export TMPDIR=$BASE/.tmp
mkdir -p "$TMPDIR"

IN=./cattle_cohort_genotypes.mainchr.vcf.gz
TAG=mainchr

[ -f "$IN"    ] || { echo "MISSING $IN - run 02_subset_vcftools.sh first" >&2; exit 1; }
[ -d "$CACHE" ] || { echo "MISSING cache dir: $CACHE" >&2; exit 1; }

echo "in    : $IN"
echo "cache : $CACHE"
echo "variants going in: $(bcftools index --nrecords "$IN")"
echo

vep -i "$IN" \
    --cache --offline --dir_cache "$CACHE" --cache_version 107 \
    --merged --species bos_taurus \
    --format vcf --vcf \
    --sift b --symbol --regulatory --pick \
    --fork 16 --force_overwrite \
    -o cohort_vep107_$TAG.vcf \
    --stats_file cohort_vep107_${TAG}_summary.html

[ -f cohort_vep107_$TAG.vcf ] || { echo "VEP produced no output" >&2; exit 1; }

# Compress + index with bcftools, so bgzip/tabix are not needed.
bcftools view cohort_vep107_$TAG.vcf -Oz -o cohort_vep107_$TAG.vcf.gz --threads 4
bcftools index -t cohort_vep107_$TAG.vcf.gz
rm -f cohort_vep107_$TAG.vcf

echo
echo "STEP 3 OK."
echo "variants annotated: $(bcftools index --nrecords cohort_vep107_$TAG.vcf.gz)"
ls -lh cohort_vep107_$TAG.vcf.gz* cohort_vep107_${TAG}_summary.html
