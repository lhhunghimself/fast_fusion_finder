#!/usr/bin/env bash
# polishFusionConsensus.sh
#
# For each cluster in a BFF-style breakpoint_files/ directory, build a
# chimeric reference (5' partner window + 3' partner window in transcript
# order, RC where the breakpoint encoding says -), realign the supporting
# reads to it with minimap2 (map-ont preset), then polish with racon to
# produce a primer-design-grade consensus FASTA. Output is one record
# per cluster, named after the cluster.
#
# Why this exists: pileup-against-genome (processBreakpointDir.sh) is fast
# and gives a sensible per-side sequence at the breakpoint, but at low
# cluster depth on noisy long reads it doesn't denoise indels well, and
# bases right at the junction tend to be soft-clipped in the genome
# alignment. Realigning to a chimeric ref + racon polishing gives one
# clean record across the junction with indels handled coherently.
#
# Usage:
#   polishFusionConsensus.sh INPUT_DIR REFERENCE_FASTA OUTPUT_DIR [WINDOW]
#
# Defaults: WINDOW=300 (each side of the junction). Bigger than the
# pileup script's 100 because minimap2 needs more anchor to map noisy
# long reads to a small chimeric contig.
#
# Requires samtools, minimap2, racon, bash, awk on PATH. The
# biodepot/bff:test image bundles all of these.

set -euo pipefail

inDir=${1:-}
ref=${2:-}
outDir=${3:-}
window=${4:-300}

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

# Detect contig naming: do we need to strip "chr" off the cluster SAMs to
# match the reference?
faiHasChr=0
firstFaiName=$(awk 'NR==1{print $1; exit}' "$ref.fai")
[[ "$firstFaiName" == chr* ]] && faiHasChr=1

# SAM header for sorting cluster alignments.
header=$(mktemp)
trap 'rm -f "$header"' EXIT
{
    printf '@HD\tVN:1.6\tSO:coordinate\n'
    awk -v OFS='\t' '{print "@SQ", "SN:" $1, "LN:" $2}' "$ref.fai"
} >"$header"

revcomp() {
    # stdin -> reverse complement on stdout. Preserves N case.
    tr 'ACGTacgtNn' 'TGCAtgcaNn' | rev
}

faExtract() {
    # samtools faidx, strip the header, join into one line.
    samtools faidx "$ref" "$1" | awk 'NR>1' | tr -d '\n'
}

shopt -s nullglob
processed=0
skipped=0
fellback=0
for samFile in "$inDir"/*.nearby.sam; do
    base=$(basename "$samFile" .nearby.sam)
    coordsFile="$inDir/${base}.nearby.readnamescoords"
    if [[ ! -f "$coordsFile" ]]; then
        echo "warn: missing $coordsFile, skipping $base" >&2
        skipped=$((skipped+1))
        continue
    fi

    # Parse canonical bp from filename: chrA-posA-chrB-posB
    if [[ ! "$base" =~ ^([^-]+)-([0-9]+)-([^-]+)-([0-9]+)$ ]]; then
        echo "warn: $base does not match chrA-posA-chrB-posB, skipping" >&2
        skipped=$((skipped+1))
        continue
    fi
    chrA="${BASH_REMATCH[1]}"; posA="${BASH_REMATCH[2]}"
    chrB="${BASH_REMATCH[3]}"; posB="${BASH_REMATCH[4]}"

    # Pull a representative encoding to recover strand. Prefer the row that
    # matches the filename's bp, fall back to the first row.
    encoding=$(awk -v FS='\t' -v pA="$posA" -v pB="$posB" '
        {
            n = split($2, sides, ";")
            if (n != 2) next
            split(sides[1], p1, ":")
            split(sides[2], p2, ":")
            if (p1[3] == pA && p2[2] == pB) { print $2; exit }
        }' "$coordsFile")
    [[ -z "$encoding" ]] && encoding=$(awk 'NR==1{print $2; exit}' "$coordsFile")

    side1Encoded="${encoding%%;*}"
    side2Encoded="${encoding#*;}"
    s1Start=$(echo "$side1Encoded" | cut -d: -f2)
    s1End=$(echo   "$side1Encoded" | cut -d: -f3)
    s2Start=$(echo "$side2Encoded" | cut -d: -f2)
    s2End=$(echo   "$side2Encoded" | cut -d: -f3)
    [[ "$s1Start" -gt "$s1End" ]] && dirA='-' || dirA='+'
    [[ "$s2Start" -gt "$s2End" ]] && dirB='-' || dirB='+'

    # Resolve contig names against the reference (strip chr if FAI is bare).
    refChrA="$chrA"; refChrB="$chrB"
    if [[ $faiHasChr -eq 0 ]]; then
        [[ "$chrA" == chr* ]] && refChrA="${chrA#chr}"
        [[ "$chrB" == chr* ]] && refChrB="${chrB#chr}"
    fi

    # Window placement (same logic as breakpointPileup.pl):
    #   5'+ or 3'-: window ends at bp going up genomically  -> [bp-W+1, bp]
    #   5'- or 3'+: window starts at bp going up genomically -> [bp, bp+W-1]
    if [[ $dirA == '+' ]]; then loA=$((posA - window + 1)); hiA=$posA
    else                        loA=$posA; hiA=$((posA + window - 1)); fi
    if [[ $dirB == '-' ]]; then loB=$((posB - window + 1)); hiB=$posB
    else                        loB=$posB; hiB=$((posB + window - 1)); fi
    [[ $loA -lt 1 ]] && loA=1
    [[ $loB -lt 1 ]] && loB=1

    work=$(mktemp -d)

    # 1. Build chimeric reference (5' + 3', in transcript order).
    sideA=$(faExtract "${refChrA}:${loA}-${hiA}")
    sideB=$(faExtract "${refChrB}:${loB}-${hiB}")
    [[ $dirA == '-' ]] && sideA=$(printf '%s' "$sideA" | revcomp)
    [[ $dirB == '-' ]] && sideB=$(printf '%s' "$sideB" | revcomp)
    chimericName="${chrA}_${posA}${dirA}_${chrB}_${posB}${dirB}"
    {
        printf '>%s\n' "$chimericName"
        printf '%s%s\n' "$sideA" "$sideB"
    } >"$work/chimeric.fa"
    samtools faidx "$work/chimeric.fa"

    # 2. Pull supporting reads' sequences out of the cluster SAM. The SAM
    # is headerless; prepend the FAI-derived header, sort, then samtools
    # fasta gives reads back in their original orientation.
    if [[ $faiHasChr -eq 0 ]]; then
        {
            cat "$header"
            awk -v FS='\t' -v OFS='\t' \
                '{ sub(/^chr/, "", $3); if ($7 != "*" && $7 != "=") sub(/^chr/, "", $7); print }' \
                "$samFile"
        }
    else
        cat "$header" "$samFile"
    fi | samtools sort -O bam -o "$work/cluster.bam" -
    samtools fasta "$work/cluster.bam" 2>/dev/null >"$work/reads.fa" || true
    if [[ ! -s "$work/reads.fa" ]]; then
        echo "warn: no reads extracted for $base, skipping" >&2
        rm -rf "$work"
        skipped=$((skipped+1))
        continue
    fi

    # 3. Realign reads to the chimeric ref. SAM (not BAM) for racon.
    minimap2 -ax map-ont --secondary=no \
        "$work/chimeric.fa" "$work/reads.fa" 2>"$work/minimap2.log" \
        >"$work/realign.sam"

    # 4. Polish. Racon can fail with very low coverage / no overlaps; if
    # it does, we fall back to the unpolished chimeric draft and flag it.
    if racon -t 1 "$work/reads.fa" "$work/realign.sam" "$work/chimeric.fa" \
            2>"$work/racon.log" >"$outDir/${base}.fa"; then
        echo "ok: $base -> $outDir/${base}.fa (polished)"
    else
        cp "$work/chimeric.fa" "$outDir/${base}.fa"
        echo "warn: racon failed for $base (likely low overlap); kept draft chimeric ref" >&2
        echo "      racon log: $(cat "$work/racon.log" | tr '\n' ' ' | head -c 200)" >&2
        fellback=$((fellback+1))
    fi

    rm -rf "$work"
    processed=$((processed+1))
done

echo "processed=$processed skipped=$skipped fallback_to_draft=$fellback output=$outDir"
