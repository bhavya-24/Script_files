#!/bin/sh
# ===========================================================================
#  PREFLIGHT - confirm everything before submitting the real jobs.
#
#  Safe to run on the LOGIN NODE: every check reads metadata or an index,
#  nothing streams the VCF.  Takes seconds.
#
#      sh 00_preflight.sh
#
#  Every line is marked PASS / FAIL / WARN.  Do not qsub anything until
#  there are zero FAILs.
# ===========================================================================

. /etc/profile.d/modules.sh

BASE=/exports/eddie/scratch/$USER/Dissertation
VCF=$BASE/pangenie/work/cattle_cohort_genotypes.vcf.gz
CACHE=$BASE/vep/vep_cache
MAP=./chr_map.txt

FAILS=0
WARNS=0
pass() { echo "  PASS  $*"; }
warn() { echo "  WARN  $*"; WARNS=$((WARNS+1)); }
fail() { echo "  FAIL  $*"; FAILS=$((FAILS+1)); }

echo "==========================================================="
echo " PREFLIGHT   $(date)"
echo " host        $(hostname)"
echo " directory   $PWD"
echo "==========================================================="
echo

# ---------------------------------------------------------------------------
echo "[1] MODULES"
# ---------------------------------------------------------------------------
if module load igmm/apps/bcftools/1.20 2>/dev/null && command -v bcftools >/dev/null; then
  pass "bcftools  -> $(bcftools --version | head -1)"
else
  fail "igmm/apps/bcftools/1.20 did not load"
  echo "        available: $(module avail 2>&1 | tr ' ' '\n' | grep -i '/bcftools/' | tr '\n' ' ')"
fi

if module load igmm/apps/vep/107 2>/dev/null && command -v vep >/dev/null; then
  pass "vep       -> $(command -v vep)"
else
  fail "igmm/apps/vep/107 did not load"
  echo "        available: $(module avail 2>&1 | tr ' ' '\n' | grep -i '/vep/' | tr '\n' ' ')"
fi

VCFT=$(module avail 2>&1 | tr ' ' '\n' | grep -iE '/vcftools/' | sort -V | tail -1)
if [ -n "$VCFT" ]; then
  pass "vcftools module exists -> $VCFT   (optional; bcftools does step 2)"
else
  warn "no vcftools module found - use bcftools view -r for the subset"
fi
echo

# ---------------------------------------------------------------------------
echo "[2] INPUT VCF"
# ---------------------------------------------------------------------------
if [ -f "$VCF" ]; then
  pass "exists  $VCF"
  pass "size    $(du -h "$VCF" | cut -f1)"
else
  fail "MISSING $VCF"
fi

IDX=""
for e in .tbi .csi; do [ -f "$VCF$e" ] && IDX="$VCF$e"; done
if [ -n "$IDX" ]; then
  pass "index   $(basename "$IDX")"
else
  fail "no .tbi or .csi beside the VCF - run: bcftools index -t $VCF"
fi
echo

# ---------------------------------------------------------------------------
echo "[3] chr_map.txt"
# ---------------------------------------------------------------------------
if [ -f "$MAP" ]; then
  NLINE=$(grep -cve '^[[:space:]]*$' "$MAP")
  pass "exists, $NLINE non-blank lines (expect 31)"
  [ "$NLINE" -eq 31 ] || warn "expected 31 lines, found $NLINE"

  if grep -qP '\t' "$MAP" 2>/dev/null || awk '{if (NF!=2) bad=1} END{exit bad+0}' "$MAP"; then
    pass "two whitespace-separated columns throughout"
  else
    fail "some line does not have exactly 2 columns"
  fi

  if awk 'NF==2 && $1 ~ /^NC_[0-9]+\.[0-9]+$/ {n++} END{exit !(n==NR || n==NR-1)}' "$MAP"; then
    pass "all keys look like valid RefSeq accessions"
  else
    warn "a key does not match NC_nnnnnn.n - check for typos"
  fi
else
  fail "MISSING $MAP in $PWD"
fi
echo

# ---------------------------------------------------------------------------
echo "[4] MAP vs VCF  - the silent-failure check"
# ---------------------------------------------------------------------------
# bcftools annotate --rename-chrs does NOT error when a key matches nothing.
# A single wrong accession would pass through unnoticed.  Compare directly.
if [ -f "$VCF" ] && [ -n "$IDX" ] && [ -f "$MAP" ] && command -v bcftools >/dev/null; then
  bcftools index --stats "$VCF" | cut -f1 | sort -u > .pf_vcf_contigs
  awk 'NF==2 {print $1}' "$MAP" | sort -u    > .pf_map_keys

  UNMATCHED=$(comm -23 .pf_map_keys .pf_vcf_contigs)
  if [ -z "$UNMATCHED" ]; then
    pass "all 31 map keys are present in the VCF - rename will fully apply"
  else
    fail "these map keys are NOT in the VCF (they will silently do nothing):"
    echo "$UNMATCHED" | sed 's/^/          /'
  fi

  TOTAL=$(wc -l < .pf_vcf_contigs)
  NW=$(grep -c '^NW_' .pf_vcf_contigs)
  pass "VCF has $TOTAL contigs: 31 main + $NW unplaced (NW_)"

  # Predict the subset result now, from the index alone.
  N_ALL=$(bcftools index --stats "$VCF" | awk '{n+=$3} END{print n+0}')
  N_NW=$(bcftools index --stats "$VCF" | awk '$1 ~ /^NW_/ {n+=$3} END{print n+0}')
  echo
  echo "        PREDICTED OUTCOME (from the index, no run needed):"
  echo "          total variants        $N_ALL"
  echo "          on NW_ scaffolds      $N_NW   <- will be dropped"
  echo "          kept after subset     $((N_ALL - N_NW))"
  echo "        Quote these in Methods."

  rm -f .pf_vcf_contigs .pf_map_keys
else
  warn "skipped - needs the VCF, its index, the map, and bcftools"
fi
echo

# ---------------------------------------------------------------------------
echo "[5] VEP CACHE"
# ---------------------------------------------------------------------------
if [ -d "$CACHE" ]; then
  pass "exists  $CACHE"
  SP=$(ls "$CACHE" 2>/dev/null | grep -i bos_taurus | tr '\n' ' ')
  if [ -n "$SP" ]; then
    pass "species dirs: $SP"
    if ls -d "$CACHE"/bos_taurus_merged/107* >/dev/null 2>&1; then
      pass "merged cache v107 present: $(ls -d "$CACHE"/bos_taurus_merged/107* | xargs -n1 basename | tr '\n' ' ')"
    else
      fail "no 107* under bos_taurus_merged/ - --merged --cache_version 107 will fail"
      echo "        found: $(ls "$CACHE"/bos_taurus_merged 2>/dev/null | tr '\n' ' ')"
    fi
  else
    fail "no bos_taurus directory inside the cache"
  fi
else
  fail "MISSING cache dir $CACHE"
fi
echo

# ---------------------------------------------------------------------------
echo "[6] DISK AND WRITE ACCESS"
# ---------------------------------------------------------------------------
AVAIL=$(df -Ph "$BASE" | tail -1 | awk '{print $4}')
pass "free on scratch: $AVAIL"
df -Ph "$BASE" | tail -1 | awk '{gsub("%","",$5); if ($5+0 > 90) print "  WARN  scratch is "$5"% full"}'

if touch .pf_write_test 2>/dev/null; then
  pass "submit directory is writable"
  rm -f .pf_write_test
else
  fail "cannot write to $PWD"
fi
echo

# ---------------------------------------------------------------------------
echo "==========================================================="
if [ "$FAILS" -eq 0 ]; then
  echo " READY - $WARNS warning(s), 0 failures."
  echo
  echo "   qsub 01_rename_bcftools.sh"
  echo "   (then check 01_rename.*.out says STEP 1 OK)"
  echo "   qsub 02_subset_vcftools.sh"
  echo "   qsub 03_vep.sh"
else
  echo " NOT READY - $FAILS failure(s), $WARNS warning(s). Fix the FAILs above."
fi
echo "==========================================================="
exit $FAILS
