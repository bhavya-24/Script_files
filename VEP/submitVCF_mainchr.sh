#!/bin/sh
#$ -N vep_mainchr
#$ -cwd
#$ -m beas
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 16
#$ -o vep_mainchr.$JOB_ID.out
#$ -e vep_mainchr.$JOB_ID.err

# ===========================================================================
#  Two-step preparation James asked for, then VEP:
#
#    1. RENAME    NC_037328.1 -> 1 ... NC_037357.1 -> X, NC_006853.1 -> MT
#    2. SUBSET    keep only 1-29, X, MT   (drops the NW_ unplaced scaffolds)
#    3. ANNOTATE  VEP 107, merged cache, bos_taurus
#
#  bcftools comes from the Eddie module.  The Cactus sandbox was tried first
#  (2026-07-30) but has no bcftools on its PATH and is missing /tmp and
#  /etc/passwd, i.e. it is incomplete - see the log entry.
#
#  Renames the INPUT, not the annotated output, because the summary HTML is
#  written during the VEP run.
#
#  Set SUBSET=no to keep the unplaced scaffolds and annotate everything.
# ===========================================================================

SUBSET=yes

. /etc/profile.d/modules.sh
module load igmm/apps/bcftools/1.20
module load igmm/apps/vep/107

BASE=/exports/eddie/scratch/$USER/Dissertation
VCF=$BASE/pangenie/work/cattle_cohort_genotypes.vcf.gz
CACHE=$BASE/vep/vep_cache

export TMPDIR=$BASE/.tmp
mkdir -p "$TMPDIR"

MAP=./chr_map.txt
RENAMED=./cattle_cohort_genotypes.renamed.vcf.gz
MAINCHR=./cattle_cohort_genotypes.mainchr.vcf.gz

REGIONS=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,X,MT

echo "working directory : $PWD"
echo "input VCF         : $VCF"
echo "cache             : $CACHE"
echo "subset to main    : $SUBSET"
echo

# --- fail in seconds, not after hours in the queue ------------------------
[ -f "$VCF"   ] || { echo "MISSING input VCF: $VCF" >&2; exit 1; }
[ -f "$MAP"   ] || { echo "MISSING chr_map.txt in $PWD" >&2; exit 1; }
[ -d "$CACHE" ] || { echo "MISSING cache dir: $CACHE" >&2; exit 1; }
command -v bcftools >/dev/null || { echo "bcftools NOT on PATH" >&2; exit 1; }
command -v vep      >/dev/null || { echo "vep NOT on PATH"      >&2; exit 1; }
bcftools --version | head -2
echo

# ===========================================================================
# 1. RENAME
# ===========================================================================
echo "=== 1. renaming chromosomes ==="
echo "before:"
bcftools index --stats $VCF | cut -f1 | head -5

bcftools annotate --rename-chrs $MAP $VCF -Oz -o $RENAMED --threads 4
bcftools index -t $RENAMED

echo "after:"
bcftools index --stats $RENAMED | cut -f1 | head -35

# bcftools does NOT error when the map matches nothing, so check explicitly.
if bcftools index --stats $RENAMED | cut -f1 | grep -q '^NC_037'; then
  echo "ERROR: NC_037 names still present - chr_map.txt did not match the VCF." >&2
  exit 1
fi
echo "OK: chromosomes renamed."
echo

# ===========================================================================
# 2. SUBSET to the main chromosomes
# ===========================================================================
if [ "$SUBSET" = "yes" ]; then
  echo "=== 2. restricting to 1-29, X, MT ==="
  N_BEFORE=$(bcftools index --nrecords $RENAMED)

  bcftools view $RENAMED -r $REGIONS -Oz -o $MAINCHR --threads 4
  bcftools index -t $MAINCHR

  N_AFTER=$(bcftools index --nrecords $MAINCHR)
  echo "variants before : $N_BEFORE"
  echo "variants after  : $N_AFTER"
  echo "dropped         : $((N_BEFORE - N_AFTER))   <-- quote this number in Methods"
  echo
  echo "contigs kept:"
  bcftools index --stats $MAINCHR | cut -f1

  INPUT=$MAINCHR
  TAG=mainchr
else
  echo "=== 2. SKIPPED - keeping unplaced scaffolds ==="
  INPUT=$RENAMED
  TAG=renamed
fi
echo

# ===========================================================================
# 3. ANNOTATE
# ===========================================================================
echo "=== 3. running VEP on $INPUT ==="
vep -i $INPUT --cache --offline --dir_cache $CACHE --cache_version 107 \
  --merged --species bos_taurus --format vcf --vcf \
  --sift b --symbol --regulatory --pick --fork 16 --force_overwrite \
  -o cohort_vep107_$TAG.vcf --stats_file cohort_vep107_${TAG}_summary.html

# Compress and index with bcftools, so the job does not
# depend on bgzip/tabix being separately on PATH.
bcftools view cohort_vep107_$TAG.vcf -Oz -o cohort_vep107_$TAG.vcf.gz --threads 4
bcftools index -t cohort_vep107_$TAG.vcf.gz
rm -f cohort_vep107_$TAG.vcf

echo
echo "=== done ==="
ls -lh cohort_vep107_$TAG.vcf.gz* cohort_vep107_${TAG}_summary.html
