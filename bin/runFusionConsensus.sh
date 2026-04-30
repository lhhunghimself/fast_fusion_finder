#!/usr/bin/env bash
# runFusionConsensus.sh
#
# Host-side driver. Builds biodepot/bff:test if it isn't already present,
# then runs either the polish (default) or the per-side pileup pipeline
# inside the container, and emits a TSV table of consensus sequences to
# stdout (or a file via -o).
#
# Two input modes:
#   -d INPUT_DIR        process every *.nearby.sam in INPUT_DIR
#   SAM [SAM ...]       process the given *.nearby.sam files explicitly
#                       (each must have its matching *.nearby.readnamescoords
#                        file alongside it)
#
# Output table (tab-separated):
#   cluster    side     length    sequence
# - polish mode: side = "junction", one row per cluster (chimeric +
#   racon, primer-design grade).
# - pileup mode: side = "5p" or "3p", two rows per cluster
#   (per-side mpileup consensus, lowercase variants preserved).
#
# Usage:
#   runFusionConsensus.sh [opts] -d INPUT_DIR REFERENCE.fa
#   runFusionConsensus.sh [opts]    REFERENCE.fa SAM [SAM ...]
#
# Options:
#   -m polish|pileup   pipeline (default: polish)
#   -w N               window size (default 300 polish / 100 pileup)
#   -o PATH            write TSV to PATH instead of stdout
#   -D DIR             keep the per-cluster FASTAs (and pileup TSVs) in
#                      DIR after the run; DIR is created if missing.
#                      Default: per-cluster artefacts are dropped.
#   -i IMAGE           override the docker image tag (default biodepot/bff:test)
#   -h                 this help

set -euo pipefail

MODE=polish
WINDOW=
OUT=-
SEARCH_DIR=
KEEP_DIR=
IMAGE=biodepot/bff:test

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DOCKERFILES_DIR="$SCRIPT_DIR/../widgets/fast_fusion_finder/BiodepotFusionFinder/Dockerfiles"

usage() {
    cat <<'EOF' >&2
Usage:
  runFusionConsensus.sh [opts] -d INPUT_DIR REFERENCE.fa
  runFusionConsensus.sh [opts]    REFERENCE.fa SAM [SAM ...]

Options:
  -m polish|pileup   pipeline (default: polish)
  -w N               window size (default 300 polish / 100 pileup)
  -o PATH            write TSV to PATH instead of stdout
  -D DIR             keep per-cluster FASTAs (and pileup TSVs) in DIR
  -d DIR             use every *.nearby.sam in DIR
  -i IMAGE           docker image tag (default biodepot/bff:test)
  -h                 show this help
EOF
    exit "${1:-1}"
}

while getopts "m:w:o:d:D:i:h" opt; do
    case "$opt" in
        m) MODE=$OPTARG;;
        w) WINDOW=$OPTARG;;
        o) OUT=$OPTARG;;
        d) SEARCH_DIR=$OPTARG;;
        D) KEEP_DIR=$OPTARG;;
        i) IMAGE=$OPTARG;;
        h) usage 0;;
        *) usage;;
    esac
done
shift $((OPTIND - 1))

case "$MODE" in
    polish) WRAPPER=polishFusionConsensus.sh; WINDOW=${WINDOW:-300};;
    pileup) WRAPPER=processBreakpointDir.sh;  WINDOW=${WINDOW:-100};;
    *) echo "unknown mode: $MODE (use polish or pileup)" >&2; exit 1;;
esac

command -v docker >/dev/null 2>&1 \
    || { echo "docker not found on PATH" >&2; exit 1; }

# Build the image if missing (and possible).
if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    if [[ -d "$DOCKERFILES_DIR" ]]; then
        echo "image $IMAGE not found, building from $DOCKERFILES_DIR ..." >&2
        docker build -t "$IMAGE" "$DOCKERFILES_DIR" >&2
    else
        echo "image $IMAGE not found and Dockerfile dir not at $DOCKERFILES_DIR" >&2
        echo "build the image first:" >&2
        echo "  docker build -t $IMAGE <path-to-Dockerfiles>" >&2
        exit 1
    fi
fi

# Cleanup tracking
cleanup_dirs=()
cleanup() {
    local d
    for d in "${cleanup_dirs[@]+"${cleanup_dirs[@]}"}"; do
        [[ -n "$d" && -d "$d" ]] && rm -rf "$d"
    done
}
trap cleanup EXIT

# Resolve input source.
if [[ -n "$SEARCH_DIR" ]]; then
    [[ $# -lt 1 ]] && usage
    REF=$1; shift
    [[ $# -gt 0 ]] && { echo "extra args after reference: $*" >&2; usage; }
    [[ -d "$SEARCH_DIR" ]] || { echo "directory not found: $SEARCH_DIR" >&2; exit 1; }
    INPUT_DIR=$(realpath "$SEARCH_DIR")
else
    [[ $# -lt 2 ]] && usage
    REF=$1; shift
    INPUT_DIR=$(mktemp -d)
    cleanup_dirs+=("$INPUT_DIR")
    for sam in "$@"; do
        [[ -f "$sam" ]] || { echo "sam not found: $sam" >&2; exit 1; }
        sam=$(realpath "$sam")
        base=$(basename "$sam" .nearby.sam)
        if [[ "$(basename "$sam")" == "$base" ]]; then
            echo "expected *.nearby.sam, got: $sam" >&2; exit 1
        fi
        coords="$(dirname "$sam")/${base}.nearby.readnamescoords"
        [[ -f "$coords" ]] \
            || { echo "missing $coords (expected next to $sam)" >&2; exit 1; }
        # Hardlink first (fast, same FS); fall back to copy across FSes.
        ln "$sam"    "$INPUT_DIR/${base}.nearby.sam"               2>/dev/null \
            || cp "$sam"    "$INPUT_DIR/${base}.nearby.sam"
        ln "$coords" "$INPUT_DIR/${base}.nearby.readnamescoords"   2>/dev/null \
            || cp "$coords" "$INPUT_DIR/${base}.nearby.readnamescoords"
    done
fi
[[ -f "$REF" ]] || { echo "reference not found: $REF" >&2; exit 1; }
REF=$(realpath "$REF")
REF_DIR=$(dirname "$REF")
REF_NAME=$(basename "$REF")

# Per-cluster output dir (kept on -D, otherwise scratch).
if [[ -n "$KEEP_DIR" ]]; then
    mkdir -p "$KEEP_DIR"
    OUT_DIR=$(realpath "$KEEP_DIR")
else
    OUT_DIR=$(mktemp -d)
    cleanup_dirs+=("$OUT_DIR")
fi

# Run inside the container. All wrapper progress lines go to stderr so
# stdout stays clean for the table.
docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$INPUT_DIR:/in:ro" \
    -v "$REF_DIR:/ref:ro" \
    -v "$OUT_DIR:/out" \
    "$IMAGE" \
    "$WRAPPER" /in "/ref/$REF_NAME" /out "$WINDOW" >&2

# Emit the TSV table.
[[ "$OUT" == "-" ]] && OUT=/dev/stdout
{
    printf 'cluster\tside\tlength\tsequence\n'
    shopt -s nullglob
    for fa in "$OUT_DIR"/*.fa; do
        base=$(basename "$fa" .fa)
        awk -v cluster="$base" -v mode="$MODE" '
            BEGIN { rec = 0; seq = ""; have = 0 }
            function emit() {
                side = (mode == "polish") ? "junction" : (rec == 1 ? "5p" : "3p")
                printf "%s\t%s\t%d\t%s\n", cluster, side, length(seq), seq
            }
            /^>/ {
                if (have) emit()
                have = 1; rec++; seq = ""
                next
            }
            { seq = seq $0 }
            END { if (have) emit() }
        ' "$fa"
    done
} > "$OUT"
