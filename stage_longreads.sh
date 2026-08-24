#!/bin/bash
#$ -N stage_longreads
#$ -cwd
#$ -q staging
#$ -l h_rt=12:00:00
#$ -o stage_longreads.$JOB_ID.out
#$ -e stage_longreads.$JOB_ID.err
# NOTE: no -pe and no big -l h_vmem. Transfers are I/O-bound, not CPU or memory
# bound, so extra cores/RAM only make the job queue longer for no benefit.
# NOTE: if your group requires a project code, add it and CONFIRM it first:
#   #$ -P <your_project_code>
# =============================================================================
#  Stage ONT long reads between Eddie scratch and DataStore.
#
#  DIRECTION is explicit and DRYRUN is on by default, because rsync in the
#  wrong direction can overwrite good data with old data.
#
#  Usage:
#     qsub stage_longreads.sh                 # dry run, backup direction
#     DRYRUN=no qsub stage_longreads.sh       # actually copy
#     DIRECTION=restore DRYRUN=no qsub ...    # pull from DataStore to scratch
# =============================================================================
set -euo pipefail

# ----------------------------------------------------------------- settings --
SCRATCH="/exports/eddie/scratch/$USER/Dissertation/long_reads"
DATASTORE="/exports/cmvm/datastore/eb/groups/prendergast_grp/$USER/long_reads"

DIRECTION="${DIRECTION:-backup}"   # backup  = scratch    -> DataStore
                                   # restore = DataStore  -> scratch
DRYRUN="${DRYRUN:-yes}"            # yes = show what would happen, copy nothing

# ------------------------------------------------------------------ resolve --
case "$DIRECTION" in
  backup)  SRC="$SCRATCH/";   DST="$DATASTORE/" ;;
  restore) SRC="$DATASTORE/"; DST="$SCRATCH/"   ;;
  *) echo "DIRECTION must be 'backup' or 'restore', got '$DIRECTION'" >&2; exit 2 ;;
esac

echo "=================================================================="
echo " direction : $DIRECTION"
echo " source    : $SRC"
echo " target    : $DST"
echo " dry run   : $DRYRUN"
echo " started   : $(date)"
echo "=================================================================="

[ -d "$SRC" ] || { echo "SOURCE does not exist: $SRC" >&2; exit 1; }

echo "--- source contents ---"
du -sh "$SRC" 2>/dev/null || true
find "$SRC" -maxdepth 2 -name '*.fastq.gz' -printf '%10s  %p\n' 2>/dev/null | sort -k2 || true
echo

mkdir -p "$DST"

# ---------------------------------------------------------------- transfer --
#  -a  preserve times/perms      --partial  resume interrupted large files
#  -h  human sizes               --info=progress2  one overall progress line
#  --exclude tmp/work dirs so we copy data, not scratch clutter
RSYNC_OPTS=(-a -h --partial --info=progress2
            --exclude='tmp/' --exclude='work/' --exclude='*.tmp' --exclude='js/')

if [ "$DRYRUN" = "yes" ]; then
  echo ">>> DRY RUN. Nothing will be written. Re-run with DRYRUN=no to transfer."
  rsync --dry-run --itemize-changes "${RSYNC_OPTS[@]}" "$SRC" "$DST"
  echo
  echo ">>> dry run complete. Review the list above, then:"
  echo "    DRYRUN=no qsub stage_longreads.sh"
  exit 0
fi

echo ">>> transferring..."
rsync "${RSYNC_OPTS[@]}" "$SRC" "$DST"

# ------------------------------------------------------------- verification --
echo
echo "--- verifying (checksum comparison, no data copied) ---"
# rsync -c compares checksums; any output here means a file did NOT match.
if rsync -rc --dry-run --itemize-changes "${RSYNC_OPTS[@]}" "$SRC" "$DST" \
     | grep -v '^$' | grep . ; then
  echo "!!! WARNING: the files above differ after transfer. Investigate before"
  echo "!!! deleting anything from the source."
  exit 1
else
  echo "OK: all files match by checksum."
fi

echo
echo "--- final sizes ---"
du -sh "$SRC" "$DST"
echo "finished : $(date)"
