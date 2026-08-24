#!/bin/bash
#$ -N stage_multi
#$ -cwd
#$ -q staging
#$ -l h_rt=24:00:00

DEST=/exports/eddie/scratch/$USER/Dissertation
mkdir -p "$DEST"

# Only WH-prefixed files from 00_fastq
rsync -av --include='WH*' --exclude='*' \
  /exports/cmvm/datastore/eb/groups/prendergast_grp/jamesp/Data/Annogen_WHD_WGS/Fastq/40-670701359/00_fastq/ \
  "$DEST/"

# Whole second_prep folder
rsync -av \
  /exports/cmvm/datastore/eb/groups/prendergast_grp/jamesp/Data/Annogen_WHD_WGS/Fastq/40-670701359/second_prep \
  "$DEST/"