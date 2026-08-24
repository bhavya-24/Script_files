#!/bin/sh
#$ -N wsh_08_filter
#$ -cwd
#$ -m beas
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 2
#$ -o 08_filter.$JOB_ID.out
#$ -e 08_filter.$JOB_ID.err

# ===========================================================================
#  CORRECTED SEGREGATION FILTER  -  four groups, not two
#
#  Melissa's sample table (2026-07-30) showed WHD04, WHD05 and WHD08 are MALE.
#  The phenotype is female-limited, so the cohort is not 10 cases + 3 controls.
#  It is four groups with four different requirements:
#
#    AFFECTED FEMALES (8)   300042 401868 601898 WHD06 WHD07 WHD09 WHD10 WHD11
#        must be HOMOZYGOUS (1/1).  They have the disease, so under a recessive
#        model both inherited copies must carry the variant.
#
#    SIRES (2)              WHD04 WHD05
#        must carry AT LEAST ONE copy (0/1 or 1/1).  Both are recorded as sires
#        of affected heifers.  A daughter needs two copies, one from each
#        parent, so the sire must have supplied one.  He may carry two and
#        still be healthy, because he has no reproductive tract to affect.
#        Demanding exactly one is as wrong as demanding two.
#
#    UNAFFECTED FEMALES (2) WHD01 WHD03
#        must NOT be homozygous.  An unaffected female with two copies
#        contradicts the model outright.  These two are the only animals in
#        the cohort that can exclude anything.  Heterozygous is fine: in a
#        breed this small most healthy animals are expected to be carriers.
#
#    UNAFFECTED MALE (1)    WHD08
#        IGNORED ENTIRELY.  He would look healthy whether he carried two
#        copies, one or none, so his genotype cannot falsify anything.
#        Counting him as an exclusion discards true candidates for free.
#
#  GRADED THRESHOLD
#  The affected-female requirement slides from 8 down to MIN_AFF while the
#  sire and control conditions stay fixed.  A genuine variant can fail the
#  strictest rung on one badly genotyped animal, and only the case side should
#  move so that any change between rungs is attributable to it alone.
#
#  Outputs
#    corrected_candidates.tsv   annotated rows, one per surviving variant
#    corrected_summary.txt      counts per rung, and per 100 kb window
# ===========================================================================

MIN_AFF=6      # lowest rung: this many of the eight affected females homozygous

. /etc/profile.d/modules.sh
module load igmm/apps/bcftools/1.20
command -v bcftools >/dev/null || { echo "bcftools NOT on PATH" >&2; exit 1; }

IN=./cohort_vep107_mainchr.vcf.gz
[ -f "$IN" ] || { echo "MISSING $IN" >&2; exit 1; }

AFFECTED="300042,401868,601898,WHD06,WHD07,WHD09,WHD10,WHD11"
SIRES="WHD04,WHD05"
CTL_FEMALE="WHD01,WHD03"
CTL_MALE="WHD08"

SAMPLES=$(bcftools query -l "$IN" | tr '\n' ',' | sed 's/,$//')
CSQFMT=$(bcftools view -h "$IN" | grep -m1 '^##INFO=<ID=CSQ' \
         | sed 's/.*Format: //; s/">[[:space:]]*$//')
[ -n "$CSQFMT" ] || { echo "No CSQ header - was VEP run with --vcf?" >&2; exit 1; }

echo "input           : $IN"
echo "samples in VCF  : $SAMPLES"
echo "affected females: $AFFECTED"
echo "sires           : $SIRES"
echo "control females : $CTL_FEMALE"
echo "control male    : $CTL_MALE  (ignored)"
echo "lowest rung     : $MIN_AFF of 8"
echo

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ[\t%GT]\n' "$IN" | \
awk -F'\t' -v OFS='\t' \
    -v samples="$SAMPLES" -v affected="$AFFECTED" -v sires="$SIRES" \
    -v ctlf="$CTL_FEMALE" -v ctlm="$CTL_MALE" \
    -v csqfmt="$CSQFMT" -v minaff="$MIN_AFF" '
# R hom-ref | H hom-alt | T het | M missing
function cls(g,   a, n) {
  if (g == "" || g == "." || g == "./." || g == ".|.") return "M"
  n = split(g, a, /[\/|]/)
  if (n != 2)                     return "M"
  if (a[1] == "." || a[2] == ".") return "M"
  if (a[1] == a[2])               return (a[1] == "0") ? "R" : "H"
  return "T"
}
function fld(arr, i) { return (i > 0 && i in arr && arr[i] != "") ? arr[i] : "." }

BEGIN {
  ns = split(samples, sn, ",")
  split(affected, a, ","); for (i in a) isaff[a[i]]  = 1
  split(sires,    s, ","); for (i in s) issire[s[i]] = 1
  split(ctlf,     f, ","); for (i in f) isctlf[f[i]] = 1
  split(ctlm,     m, ","); for (i in m) isctlm[m[i]] = 1

  naff = 0; nsire = 0; nctlf = 0
  for (i = 1; i <= ns; i++) {
    if      (sn[i] in isaff)  { naff++;  aidx[naff]  = i; anm[naff] = sn[i] }
    else if (sn[i] in issire) { nsire++; sidx[nsire] = i; snm[nsire] = sn[i] }
    else if (sn[i] in isctlf) { nctlf++; fidx[nctlf] = i }
    else if (sn[i] in isctlm) { }   # deliberately unused
    else printf("WARNING: %s in neither group\n", sn[i]) > "/dev/stderr"
  }
  printf("resolved: %d affected females, %d sires, %d control females\n",
         naff, nsire, nctlf) > "/dev/stderr"
  if (naff != 8 || nsire != 2 || nctlf != 2)
    printf("WARNING: expected 8 / 2 / 2\n") > "/dev/stderr"

  nf = split(csqfmt, ff, "|")
  for (i = 1; i <= nf; i++) {
    if (ff[i] == "Consequence") i_cons = i
    if (ff[i] == "IMPACT")      i_imp  = i
    if (ff[i] == "SYMBOL")      i_sym  = i
    if (ff[i] == "Gene")        i_gene = i
    if (ff[i] == "BIOTYPE")     i_bio  = i
    if (ff[i] == "SIFT")        i_sift = i
  }
  printf("CSQ: Consequence=%d IMPACT=%d SYMBOL=%d SIFT=%d\n",
         i_cons, i_imp, i_sym, i_sift) > "/dev/stderr"

  print "CHROM","POS","REF","ALT","tier","aff_homalt","aff_het","aff_ref", \
        "aff_miss","aff_not_homalt","sire1_gt","sire2_gt","ctlF_homalt", \
        "consequence","impact","symbol","gene","biotype","sift_term"
}
{
  # --- affected females -----------------------------------------------------
  ah = 0; at = 0; ar = 0; am = 0; notyet = ""
  for (i = 1; i <= naff; i++) {
    c = cls($(5 + aidx[i]))
    if (c == "H") ah++
    else {
      if      (c == "T") at++
      else if (c == "R") ar++
      else               am++
      notyet = notyet (notyet ? "," : "") anm[i]
    }
  }
  if (ah < minaff) next

  # --- sires: at least one copy each ---------------------------------------
  ok = 1
  for (i = 1; i <= nsire; i++) {
    g = $(5 + sidx[i]); c = cls(g)
    sgt[i] = g
    if (c != "H" && c != "T") ok = 0     # R or M means no copy: fails
  }
  if (!ok) next

  # --- unaffected females: must not be homozygous --------------------------
  fh = 0
  for (i = 1; i <= nctlf; i++) if (cls($(5 + fidx[i])) == "H") fh++
  if (fh > 0) next

  # WHD08 never consulted.

  # --- annotation ----------------------------------------------------------
  csq = $5
  if (csq == "." || csq == "") {
    cons = "."; imp = "."; sym = "."; gene = "."; bio = "."; sift = "."
  } else {
    split(csq, blocks, ","); split(blocks[1], cf, "|")
    cons = fld(cf, i_cons); imp  = fld(cf, i_imp)
    sym  = fld(cf, i_sym);  gene = fld(cf, i_gene)
    bio  = fld(cf, i_bio);  sift = fld(cf, i_sift)
  }

  print $1, $2, $3, $4, ah, ah, at, ar, am, notyet, sgt[1], sgt[2], fh, \
        cons, imp, sym, gene, bio, sift
}
' > corrected_candidates.tsv

N=$(( $(wc -l < corrected_candidates.tsv) - 1 ))
echo "surviving variants: $N"
echo

{
echo "=== per rung (affected females homozygous) ==="
printf "%-6s %-12s %-12s\n" "rung" "exactly" "cumulative"
for k in 8 7 6; do
  ex=$(awk -F'\t' -v k=$k 'NR>1 && $5==k' corrected_candidates.tsv | wc -l)
  cu=$(awk -F'\t' -v k=$k 'NR>1 && $5>=k' corrected_candidates.tsv | wc -l)
  printf "%-6s %-12s %-12s\n" "$k/8" "$ex" "$cu"
done

echo
echo "=== strictest rung (8/8): top 100 kb windows ==="
awk -F'\t' 'NR>1 && $5==8 {printf "%s\t%.1f\n", $1, $2/1000000}' corrected_candidates.tsv \
  | sort | uniq -c | sort -rn | head -20

echo
echo "=== strictest rung: impact ==="
awk -F'\t' 'NR>1 && $5==8 {print $15}' corrected_candidates.tsv | sort | uniq -c | sort -rn

echo
echo "=== strictest rung: named genes ==="
awk -F'\t' 'NR>1 && $5==8 && $16!="." {print $16}' corrected_candidates.tsv \
  | sort | uniq -c | sort -rn | head -30

echo
echo "=== strictest rung: HIGH and MODERATE impact variants ==="
awk -F'\t' 'NR==1 || ($5==8 && ($15=="HIGH" || $15=="MODERATE"))' corrected_candidates.tsv \
  | cut -f1,2,3,4,14,15,16,19
} | tee corrected_summary.txt

echo
ls -lh corrected_candidates.tsv corrected_summary.txt
