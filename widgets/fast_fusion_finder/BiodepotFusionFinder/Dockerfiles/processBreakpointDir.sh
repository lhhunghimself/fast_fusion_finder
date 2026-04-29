#!/usr/bin/env bash
# processBreakpointDir.sh
#
# Build fusion consensus sequences for every breakpoint cluster in a
# directory whose layout matches the BFF "breakpoint_files/" output —
# pairs of {name}.nearby.sam (alignment dump, no @HD/@SQ) plus
# {name}.nearby.readnamescoords (TSV: read_name<TAB>chr:s:e;chr:s:e).
#
# For each cluster we:
#   1. Prepend a header derived from the reference's .fai so samtools
#      can sort the SAM into a coordinate-sorted, indexed BAM.
#   2. Derive the reads list from column 1 of the .readnamescoords.
#   3. Pick the modal encoding from column 2 as the breakpoint file row
#      (matches the canonical breakpoint of the cluster).
#   4. Run breakpointPileup.pl, writing {name}.fa (consensus) and
#      {name}.tsv (per-position counts) into the output directory.
#
# If the SAM uses "chr"-prefixed contig names but the reference does not
# (e.g. cellranger's Ensembl-style FASTA), we strip the "chr" prefix on
# the fly from both the SAM RNAME/RNEXT and the breakpoint encoding.
#
# Usage:
#   processBreakpointDir.sh INPUT_DIR REFERENCE_FASTA OUTPUT_DIR [WINDOW]
#
# Requires samtools, perl, and breakpointPileup.pl on PATH.

set -euo pipefail

inDir=${1:-}
ref=${2:-}
outDir=${3:-}
window=${4:-100}

if [[ -z "$inDir" || -z "$ref" || -z "$outDir" ]]; then
    echo "usage: $0 INPUT_DIR REFERENCE_FASTA OUTPUT_DIR [WINDOW]" >&2
    exit 1
fi

[[ -d "$inDir" ]] || { echo "input dir not found: $inDir" >&2; exit 1; }
[[ -f "$ref"  ]]  || { echo "reference not found: $ref"   >&2; exit 1; }

if [[ ! -f "$ref.fai" ]]; then
    echo "indexing $ref (one-time)..." >&2
    samtools faidx "$ref"
fi

mkdir -p "$outDir"

# --- decide whether SAM's "chr"-style names need stripping --------------------

faiHasChr=0
firstFaiName=$(awk 'NR==1{print $1; exit}' "$ref.fai")
[[ "$firstFaiName" == chr* ]] && faiHasChr=1

# --- header (one @SQ per FAI contig, no @PG) ---------------------------------

header=$(mktemp)
trap 'rm -f "$header"' EXIT
{
    printf '@HD\tVN:1.6\tSO:coordinate\n'
    awk -v OFS='\t' '{print "@SQ", "SN:" $1, "LN:" $2}' "$ref.fai"
} >"$header"

# --- per-cluster work -------------------------------------------------------

shopt -s nullglob
processed=0
skipped=0
for samFile in "$inDir"/*.nearby.sam; do
    base=$(basename "$samFile" .nearby.sam)
    coordsFile="$inDir/${base}.nearby.readnamescoords"
    if [[ ! -f "$coordsFile" ]]; then
        echo "warn: missing $coordsFile, skipping $base" >&2
        skipped=$((skipped+1))
        continue
    fi

    # Decide whether this cluster needs chr-strip (per-cluster in case the
    # input dir mixes naming, though normally it won't).
    samHasChr=0
    samFirstRname=$(awk -F'\t' '!/^@/{print $3; exit}' "$samFile")
    [[ "$samFirstRname" == chr* ]] && samHasChr=1
    stripChr=0
    if [[ $samHasChr -eq 1 && $faiHasChr -eq 0 ]]; then stripChr=1; fi

    work=$(mktemp -d)

    # 1. header + (optionally renamed) records -> sorted+indexed BAM
    if [[ $stripChr -eq 1 ]]; then
        {
            cat "$header"
            awk -v FS='\t' -v OFS='\t' \
                '{ sub(/^chr/, "", $3); if ($7 != "*" && $7 != "=") sub(/^chr/, "", $7); print }' \
                "$samFile"
        }
    else
        cat "$header" "$samFile"
    fi | samtools sort -O bam -o "$work/cluster.bam" -
    samtools index "$work/cluster.bam"

    # 2. reads list
    cut -f1 "$coordsFile" | sort -u >"$work/reads.txt"

    # 3. breakpoint file: prefer the row whose bp matches the filename's
    #    canonical (chrA-posA-chrB-posB) bp pair; fall back to the modal
    #    encoding across rows when no row matches (e.g. clusters named
    #    differently or filenames that don't follow the convention).
    chosen=""
    if [[ "$base" =~ ^([^-]+)-([0-9]+)-([^-]+)-([0-9]+)$ ]]; then
        fnPosA="${BASH_REMATCH[2]}"
        fnPosB="${BASH_REMATCH[4]}"
        chosen=$(awk -v FS='\t' -v pA="$fnPosA" -v pB="$fnPosB" '
            {
                n = split($2, sides, ";")
                if (n != 2) next
                split(sides[1], p1, ":")
                split(sides[2], p2, ":")
                if (p1[3] == pA && p2[2] == pB) { print $2; exit }
            }' "$coordsFile")
    fi
    if [[ -z "$chosen" ]]; then
        chosen=$(cut -f2 "$coordsFile" | sort | uniq -c | sort -rn \
                 | head -1 | awk '{print $2}')
        echo "warn: no row in $coordsFile matched filename's bp; using modal" >&2
    fi
    if [[ $stripChr -eq 1 ]]; then
        printf '%s\n' "$chosen" | sed 's/\(^\|;\)chr/\1/g' >"$work/bp.txt"
    else
        printf '%s\n' "$chosen" >"$work/bp.txt"
    fi

    if [[ ! -s "$work/bp.txt" ]]; then
        echo "warn: no breakpoint encoding in $coordsFile, skipping $base" >&2
        rm -rf "$work"
        skipped=$((skipped+1))
        continue
    fi

    # 4. consensus + variant table
    breakpointPileup.pl \
        --reads     "$work/reads.txt" \
        --bam       "$work/cluster.bam" \
        --ref       "$ref" \
        -b          "$work/bp.txt" \
        --window    "$window" \
        --pileupOut "$outDir/${base}.tsv" \
        >"$outDir/${base}.fa"

    rm -rf "$work"
    echo "ok: $base -> $outDir/${base}.fa, $outDir/${base}.tsv"
    processed=$((processed+1))
done

echo "processed=$processed skipped=$skipped output=$outDir"
