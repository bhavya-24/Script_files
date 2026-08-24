#!/bin/sh
#$ -N wsh_07_extract
#$ -cwd
#$ -m beas
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 2
#$ -o 07_extract.$JOB_ID.out
#$ -e 07_extract.$JOB_ID.err

# ===========================================================================
#  EXTRACT GENOTYPE COUNTS + VEP ANNOTATION FOR THE GRADED FILTER
#
#  Streams cohort_vep107_mainchr.vcf.gz once.  Nothing is held in memory, so
#  the 13.8 M variants cost nothing but time.
#
#  Two things done at runtime rather than hard-coded, because both have
#  already caused silent errors on this project:
#
#    * The sample-to-group map is built from `bcftools query -l`, not from
#      column numbers.  If the VCF order ever differs from the vault record,
#      this script follows the VCF and prints a warning.
#
#    * The CSQ subfield order is read from the ##INFO=<ID=CSQ ...> header.
#      VEP changes that order depending on which flags were used, so a
#      hard-coded index would quietly pull the wrong field.
#
#  Outputs
#    graded_full_distribution.txt   every (case_homalt, ctl_homalt) cell, tiny
#    graded_candidates.tsv          detailed rows, the R input
# ===========================================================================

MIN_CASE=6        # keep detailed rows from this many hom-alt cases upward
MAX_CTL=2         # keep rows with up to this many hom-alt controls.
                  # Deliberately loose. Under sex-limited inheritance a hom-alt
                  # MALE control is a healthy carrier and must not exclude a
                  # variant; the sexes of WHD01/WHD03/WHD08 are still
                  # outstanding, so the rows are kept and R decides.

. /etc/profile.d/modules.sh
module load igmm/apps/bcftools/1.20
command -v bcftools >/dev/null || { echo "bcftools NOT on PATH" >&2; exit 1; }

for f in ./cohort_vep107_mainchr.vcf.gz ./cohort_vep107_renamed.vcf.gz; do
  [ -f "$f" ] && { IN=$f; break; }
done
[ -n "$IN" ] || { echo "No annotated VCF found in $PWD" >&2; exit 1; }

CASES="300042,401868,601898,WHD04,WHD05,WHD06,WHD07,WHD09,WHD10,WHD11"
CONTROLS="WHD01,WHD03,WHD08"

SAMPLES=$(bcftools query -l "$IN" | tr '\n' ',' | sed 's/,$//')

# --- CSQ subfield order, straight from the header --------------------------
CSQFMT=$(bcftools view -h "$IN" | grep -m1 '^##INFO=<ID=CSQ' \
         | sed 's/.*Format: //; s/">[[:space:]]*$//')
[ -n "$CSQFMT" ] || { echo "No CSQ header found - was VEP run with --vcf?" >&2; exit 1; }

echo "input     : $IN"
echo "samples   : $SAMPLES"
echo "min_case  : $MIN_CASE"
echo "max_ctl   : $MAX_CTL"
echo "CSQ fmt   : $CSQFMT"
echo

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ[\t%GT]\n' "$IN" | \
awk -F'\t' -v OFS='\t' \
    -v samples="$SAMPLES" -v cases="$CASES" -v controls="$CONTROLS" \
    -v csqfmt="$CSQFMT" -v mincase="$MIN_CASE" -v maxctl="$MAX_CTL" \
    -v distfile="graded_full_distribution.txt" '
# --- classify one GT ------------------------------------------------------
#   R hom-ref | H hom-alt (any non-ref allele doubled) | T het | M missing
function cls(g,   a, n) {
  if (g == "" || g == "." || g == "./." || g == ".|.") return "M"
  n = split(g, a, /[\/|]/)
  if (n != 2)                       return "M"
  if (a[1] == "." || a[2] == ".")   return "M"
  if (a[1] == a[2])                 return (a[1] == "0") ? "R" : "H"
  return "T"
}
function fld(arr, i) { return (i > 0 && i in arr && arr[i] != "") ? arr[i] : "." }

BEGIN {
  ns = split(samples,  sn, ",")
  nc = split(cases,    cn, ",")
  nk = split(controls, kn, ",")
  for (i = 1; i <= nc; i++) iscase[cn[i]] = 1
  for (i = 1; i <= nk; i++) isctl[kn[i]]  = 1

  ncase = 0; nctl = 0
  for (i = 1; i <= ns; i++) {
    if (sn[i] in iscase)     { ncase++; cidx[ncase] = i; cnm[ncase] = sn[i] }
    else if (sn[i] in isctl) { nctl++;  kidx[nctl]  = i; knm[nctl]  = sn[i] }
    else printf("WARNING: sample %s is in neither group\n", sn[i]) > "/dev/stderr"
  }
  printf("resolved %d cases, %d controls from the VCF header\n", ncase, nctl) > "/dev/stderr"
  if (ncase != nc || nctl != nk)
    printf("WARNING: expected %d cases / %d controls\n", nc, nk) > "/dev/stderr"

  # locate the CSQ subfields we want, by name
  nf = split(csqfmt, ff, "|")
  for (i = 1; i <= nf; i++) {
    if (ff[i] == "Consequence") i_cons = i
    if (ff[i] == "IMPACT")      i_imp  = i
    if (ff[i] == "SYMBOL")      i_sym  = i
    if (ff[i] == "Gene")        i_gene = i
    if (ff[i] == "BIOTYPE")     i_bio  = i
    if (ff[i] == "SIFT")        i_sift = i
    if (ff[i] == "Feature")     i_feat = i
  }
  printf("CSQ indices: Consequence=%d IMPACT=%d SYMBOL=%d Gene=%d BIOTYPE=%d SIFT=%d\n",
         i_cons, i_imp, i_sym, i_gene, i_bio, i_sift) > "/dev/stderr"
  if (!i_cons || !i_imp)
    printf("WARNING: Consequence or IMPACT not found in the CSQ format\n") > "/dev/stderr"

  print "CHROM","POS","REF","ALT", \
        "case_homalt","case_het","case_homref","case_miss", \
        "ctl_homalt","ctl_het","ctl_homref","ctl_miss","cases_not_homalt", \
        "consequence","impact","symbol","gene","biotype","transcript", \
        "sift_term","sift_score"
}
{
  ch = 0; ct = 0; cr = 0; cm = 0
  kh = 0; kt = 0; kr = 0; km = 0
  nothom = ""

  for (i = 1; i <= ncase; i++) {
    c = cls($(5 + cidx[i]))
    if      (c == "H") ch++
    else {
      if      (c == "T") ct++
      else if (c == "R") cr++
      else               cm++
      nothom = nothom (nothom ? "," : "") cnm[i]
    }
  }
  for (i = 1; i <= nctl; i++) {
    c = cls($(5 + kidx[i]))
    if      (c == "H") kh++
    else if (c == "T") kt++
    else if (c == "R") kr++
    else               km++
  }

  dist[ch "\t" kh]++

  if (ch >= mincase && kh <= maxctl) {
    # --with --pick there is one CSQ block; take the first regardless.
    csq = $5
    if (csq == "." || csq == "") {
      cons = "."; imp = "."; sym = "."; gene = "."; bio = "."; feat = "."
      sterm = "."; sscore = "."
    } else {
      split(csq, blocks, ",")
      split(blocks[1], cf, "|")
      cons = fld(cf, i_cons); imp  = fld(cf, i_imp)
      sym  = fld(cf, i_sym);  gene = fld(cf, i_gene)
      bio  = fld(cf, i_bio);  feat = fld(cf, i_feat)

      # SIFT looks like  deleterious(0.01)  /  tolerated(0.42)  /  empty
      raw = fld(cf, i_sift)
      if (raw == "." ) { sterm = "."; sscore = "." }
      else {
        sterm = raw; sscore = "."
        if (match(raw, /\(/)) {
          sterm  = substr(raw, 1, RSTART - 1)
          sscore = substr(raw, RSTART + 1, length(raw) - RSTART - 1)
        }
      }
    }
    print $1, $2, $3, $4, ch, ct, cr, cm, kh, kt, kr, km, nothom, \
          cons, imp, sym, gene, bio, feat, sterm, sscore
  }
}
END {
  print "case_homalt\tctl_homalt\tn_variants" > distfile
  for (k in dist) print k "\t" dist[k] > distfile
  close(distfile)
}
' > graded_candidates.tsv

echo
echo "=== joint distribution: cases hom-alt x controls hom-alt ==="
printf "case\tctl\tn\n"
sort -k1,1n -k2,2n graded_full_distribution.txt | grep -v '^case_homalt'
echo

echo "=== candidate rows written ==="
echo "rows: $(( $(wc -l < graded_candidates.tsv) - 1 ))"
ls -lh graded_candidates.tsv graded_full_distribution.txt
echo
echo "Copy both to your PC, then: Rscript graded_filter.R"
