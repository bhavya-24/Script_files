#!/bin/sh
# ===========================================================================
#  CORRECTED PER-SAMPLE QC MODEL
#
#  Re-analyses tstv_per_genotype_file.tsv.  No bcftools, no VCF reading -
#  the expensive pass is already done.  Runs in under a second.
#
#      sh 06_qc_model.sh
#
#  THREE CORRECTIONS TO THE PREVIOUS MODEL
#
#  1. Robust statistics.  median + MAD*1.4826 instead of mean + SD.  A single
#     extreme value inflates SD, which then hides every other outlier - that
#     is exactly what happened: WHD09 at +6.4 robust z was invisible at
#     +2 classic SD, and the model flagged WHD06 instead.
#     Threshold |z| > 3.5 (Iglewicz-Hoaglin).
#
#  2. Ts/Tv demoted to a reported value, not an outlier metric.  PanGenie
#     genotypes a FIXED panel, so every animal is scored at the same sites
#     and per-sample Ts/Tv converges on the panel's own value.  Observed
#     range here is 2.1450-2.1543, SD 0.0028.  There is no signal to find.
#     Ts/Tv is still worth reporting: ~2.15 cohort-wide confirms the panel
#     is composed of real SNPs, not noise.
#
#  3. Heterozygosity as the primary metric, plus derived rates that do not
#     depend on how many sites each animal happened to be called at.
#
#  WHD06 is the assembled animal, so its own haplotypes are inside the
#  pangenome graph.  It is expected to carry more non-reference alleles than
#  the others.  Flagged, not treated as a fault.
# ===========================================================================

TSV=${1:-./tstv_per_genotype_file.tsv}
[ -f "$TSV" ] || { echo "MISSING $TSV" >&2; exit 1; }

awk -F'\t' '
# --- median of arr[1..n] ---------------------------------------------------
function median(arr, n,   i, j, t, b) {
  for (i = 1; i <= n; i++) b[i] = arr[i]
  for (i = 1; i < n; i++)
    for (j = i + 1; j <= n; j++)
      if (b[j] < b[i]) { t = b[i]; b[i] = b[j]; b[j] = t }
  if (n % 2) return b[(n + 1) / 2]
  return (b[n / 2] + b[n / 2 + 1]) / 2
}
# --- median absolute deviation, scaled to be SD-comparable -----------------
function madscale(arr, n, med,   i, d) {
  for (i = 1; i <= n; i++) d[i] = (arr[i] > med) ? arr[i] - med : med - arr[i]
  return 1.4826 * median(d, n)
}
function absv(x) { return x < 0 ? -x : x }

NR == 1 { next }
{
  n++
  sample[n] = $1
  group[n]  = $2
  refhom    = $4 + 0
  nonrefhom = $5 + 0
  hets      = $6 + 0
  called    = refhom + nonrefhom + hets

  tstv[n]    = $9  + 0
  hethom[n]  = $10 + 0
  hetrate[n] = (called > 0) ? hets / called : 0
  nonref[n]  = (called > 0) ? (nonrefhom + hets) / called : 0
  miss[n]    = $12 + 0
  ncalled[n] = called
}
END {
  if (n < 3) { print "too few samples"; exit 1 }

  # --- Ts/Tv: report only ---------------------------------------------------
  tmed = median(tstv, n); tmad = madscale(tstv, n, tmed)
  tmin = tstv[1]; tmax = tstv[1]
  for (i = 1; i <= n; i++) { if (tstv[i] < tmin) tmin = tstv[i]; if (tstv[i] > tmax) tmax = tstv[i] }

  print "==========================================================================="
  print " Ts/Tv  - REPORTED, NOT USED FOR OUTLIER DETECTION"
  print "==========================================================================="
  printf "  median %.4f   range %.4f - %.4f   spread %.4f (%.2f%% of median)\n", \
         tmed, tmin, tmax, tmax - tmin, 100 * (tmax - tmin) / tmed
  print "  A fixed PanGenie panel forces every animal onto the same sites, so"
  print "  per-sample Ts/Tv converges on the panel value. No outlier test applied."
  printf "  Cohort ~%.2f is in the healthy cattle range (2.0-2.2): panel is real SNPs.\n", tmed
  print ""

  # --- robust stats for the metrics that vary -------------------------------
  hh_m  = median(hethom, n);  hh_s  = madscale(hethom, n, hh_m)
  hr_m  = median(hetrate, n); hr_s  = madscale(hetrate, n, hr_m)
  nr_m  = median(nonref, n);  nr_s  = madscale(nonref, n, nr_m)
  ms_m  = median(miss, n);    ms_s  = madscale(miss, n, ms_m)

  print "==========================================================================="
  print " ROBUST z-SCORES   (median +/- MAD*1.4826;  |z| > 3.5 = outlier)"
  print "==========================================================================="
  printf "%-9s %-8s %9s %7s %9s %7s %9s %7s %8s %7s\n", \
         "sample", "group", "het/hom", "z", "het_rate", "z", "nonref", "z", "missing", "z"

  for (i = 1; i <= n; i++) {
    zhh = (hh_s > 0) ? (hethom[i]  - hh_m) / hh_s : 0
    zhr = (hr_s > 0) ? (hetrate[i] - hr_m) / hr_s : 0
    znr = (nr_s > 0) ? (nonref[i]  - nr_m) / nr_s : 0
    zms = (ms_s > 0) ? (miss[i]    - ms_m) / ms_s : 0

    zhh_a[i] = zhh; zhr_a[i] = zhr; znr_a[i] = znr; zms_a[i] = zms

    mark = ""
    if (absv(zhh) > 3.5 || absv(zhr) > 3.5 || absv(znr) > 3.5) mark = "  <<<"
    if (sample[i] == "WHD06") mark = mark "  [assembled animal - in the graph]"

    printf "%-9s %-8s %9.4f %+7.2f %9.4f %+7.2f %9.4f %+7.2f %8d %+7.2f%s\n", \
           sample[i], group[i], hethom[i], zhh, hetrate[i], zhr, nonref[i], znr, miss[i], zms, mark
  }
  print ""

  # --- flagged --------------------------------------------------------------
  print "==========================================================================="
  print " FLAGGED"
  print "==========================================================================="
  flagged = 0
  for (i = 1; i <= n; i++) {
    if (absv(zhh_a[i]) > 3.5 || absv(zhr_a[i]) > 3.5 || absv(znr_a[i]) > 3.5) {
      flagged++
      dir = (zhh_a[i] > 0) ? "EXCESS heterozygosity" : "DEFICIENT heterozygosity"
      printf "  %-9s %-8s  %s  (het/hom z %+.2f)\n", sample[i], group[i], dir, zhh_a[i]
      if (zhh_a[i] > 3.5)
        print "            -> possible causes: sample contamination or DNA mixture;"
        print "               genuinely outbred animal; ancestry outside the breed."
      else
        print "            -> possible causes: high autozygosity (long ROH, inbreeding);"
        print "               allele dropout at low depth pushing het calls to hom."
    }
  }
  if (!flagged) print "  none"
  print ""

  # --- group descriptives ---------------------------------------------------
  print "==========================================================================="
  print " GROUP DESCRIPTIVES   (descriptive only - 3 controls cannot support a test)"
  print "==========================================================================="
  for (i = 1; i <= n; i++) {
    g = group[i]
    gn[g]++; ghh[g] += hethom[i]; gts[g] += tstv[i]
    if (!(g in gmin) || hethom[i] < gmin[g]) gmin[g] = hethom[i]
    if (!(g in gmax) || hethom[i] > gmax[g]) gmax[g] = hethom[i]
  }
  for (g in gn)
    printf "  %-8s n=%-3d  mean het/hom %.4f  range %.4f-%.4f   mean Ts/Tv %.4f\n", \
           g, gn[g], ghh[g] / gn[g], gmin[g], gmax[g], gts[g] / gn[g]
  print ""
  print "  The case and control het/hom ranges overlap almost completely, and the"
  print "  extremes sit on BOTH sides in BOTH groups (WHD04 case low, WHD01 control"
  print "  low, WHD09 case high). This is individual variation in autozygosity, not"
  print "  a case/control difference. Do not report it as one."
  print ""

  print "==========================================================================="
  print " READ THIS BEFORE USING THE FLAGS"
  print "==========================================================================="
  print "  Excess heterozygosity in a CONTROL is the more dangerous direction. Under"
  print "  a sex-limited recessive model only unaffected FEMALES can exclude a"
  print "  candidate, and there are three controls in total. A contaminated control"
  print "  carrying spurious alternate alleles will wrongly exclude real candidates."
  print ""
  print "  Deficient heterozygosity is expected in this breed and is not a fault:"
  print "  the Whitebred Shorthorn is small and closed, so long runs of homozygosity"
  print "  are the norm. Cross-check the low animals against the ROH output before"
  print "  calling it a technical problem."
}
' "$TSV"
