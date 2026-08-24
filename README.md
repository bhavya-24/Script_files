# Whitebred Shorthorn pangenome analysis

Scripts for an MSc dissertation project searching for the variant behind an inherited infertility that affects only female Whitebred Shorthorn cattle.

The approach was to assemble one affected animal from nanopore long reads, build a pangenome graph containing that assembly alongside seven other cattle genomes, genotype all thirteen study animals against the graph, and then filter the resulting call set under a recessive model.

**This repository holds code only.** No sequence data, variant calls or genotype tables are included. Data for the thirteen animals is not ours to release. `.gitignore` blocks those file types so they cannot be committed by accident.

Every script here produced something that appears in the dissertation. Earlier and superseded versions have been removed so that what remains matches the report.

Analyses ran on the University of Edinburgh Eddie cluster under Sun Grid Engine. Paths in the shell scripts point at scratch and group storage on that system and will need changing to run anywhere else.

## The chain that produced the reported numbers

```
stagein.sh, stage_longreads.sh    raw data from group storage onto scratch
        |
01_rename_bcftools.sh             RefSeq contig names to the Ensembl convention
02_subset_vcftools.sh             keep the 29 autosomes, the X and the mitochondrion
03_vep.sh                         VEP 107 offline against the merged ARS-UCD1.2 cache
        |
07_extract_candidates.sh          genotypes and the VEP CSQ field  ->  graded_candidates.tsv
08_filter_corrected.sh            the segregation filter           ->  corrected_candidates.tsv
        |
filter_corrected_local.R          candidates at each threshold     ->  corrected_tier8/7/6.tsv
score_all_blocks.R                blocks on the 500 kb gap rule    ->  block_scores.tsv
```

`corrected_tier8/7/6.tsv` are behind Table 8 of the dissertation. `block_scores.tsv` is behind Table 9.

## Layout

```
VEP/          cluster jobs, annotation and the segregation filter
Segregation/  block grouping and candidate prioritisation, run locally in R
DOTPLOT/      assembly against reference dot plot
FIGURES/      read quality figure
*.sh (root)   staging raw data from group storage onto scratch
```

## VEP

Numbered shell scripts run in order as Grid Engine jobs.

| Script | What it does | Dissertation |
|---|---|---|
| `00_preflight.sh` | checks inputs and tool versions before anything else runs | 2.3 |
| `01_rename_bcftools.sh` | renames RefSeq contigs to the Ensembl convention | 2.10 |
| `02_subset_vcftools.sh` | keeps the twenty-nine autosomes, the X and the mitochondrion | 2.10, 4.4 |
| `03_vep.sh` | runs VEP 107 offline against the merged *Bos taurus* ARS-UCD1.2 cache | 2.10, 4.4 |
| `04_tstv_per_sample.sh` | per-animal transition-to-transversion ratios from the cohort VCF | 2.9, Table 7 |
| `05_tstv_per_genotype_file.sh` | the same figures taken from the individual per-animal files | 2.9 |
| `06_qc_model.sh` | genotype quality checks across the cohort | 4.3 |
| `07_count_genotypes.sh` | genotype class counts per animal | 4.3 |
| `07_extract_candidates.sh` | pulls genotypes and the VEP CSQ field out of the annotated call set | 2.11 |
| `08_filter_corrected.sh` | applies the segregation filter | 2.11, 4.5 |
| `submitVCF_mainchr.sh`, `submitVCF_renamed.sh` | job wrappers for the VCF stages | 2.10 |
| `rename_and_rerun_vep.sh` | reruns annotation after the contig rename | 2.10 |

## Segregation

Run locally in R on the output of `08_filter_corrected.sh`.

| Script | What it does | Dissertation |
|---|---|---|
| `filter_corrected_local.R` | applies the segregation conditions at eight, seven and six of eight affected females | 2.11, Table 8 |
| `score_all_blocks.R` | groups candidates into blocks on a 500 kb gap rule and measures each block | 2.12, Table 9 |

Note on `score_all_blocks.R`: an earlier version combined growth and frequency-filter survival with counts of functional, promoter, structural and coding variants into a single score. That score was dropped because those annotation counts overlap, so the same variants were being counted more than once. The two measures reported in the dissertation are the ones that are independent of each other.

## DOTPLOT

| Script | Produces |
|---|---|
| `make_dotplot_named.R` | Figure 4, the polished assembly against ARS-UCD1.2, built with pafr from `wsh_vs_ref.paf` |

Reference sequences are relabelled from RefSeq accessions to chromosome numbers so the figure matches the naming used after the contig rename, chromosomes are ordered by number rather than by size, and the contig-name layer is dropped because 161 labels overlap. Filters are unchanged: alignments over 500 kb, target sequences over 20 Mb, contigs over 5 Mb.

## FIGURES

| Script | Produces |
|---|---|
| `make_qc_from_html.R` | Figure 3, read length and read quality, parsed from the NanoPlot report |

Every value in that figure is parsed from the report at run time. Nothing is transcribed by hand.

## Software

Versions as used. Full details are in the dissertation methods.

hifiasm 0.25.0, NextPolish 1.4.1, Trim Galore, Merqury and meryl, BUSCO, minimap2, pafr 0.0.2, minigraph-Cactus 2.9.9, odgi, PanGenie 4.2.1, bcftools 1.20, VEP 107, Singularity 4.1.3.
