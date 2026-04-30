#!/usr/bin/env bash
# runFusionConsensus.sh
#
# Top-level host-side driver for the BFF fusion-consensus pipeline.
#
# What it does:
#   1. Makes sure the biodepot/bff:test Docker image exists. If it
#      doesn't, builds it from the in-repo Dockerfiles directory.
#   2. Stages the input *.nearby.sam files (and matching
#      *.nearby.readnamescoords) into a directory the container can see.
#   3. Runs either the polish pipeline (chimeric-ref + minimap2 + racon
#      → primer-grade junction sequence) or the per-side mpileup
#      pipeline (variant-aware, two halves per cluster) inside the
#      container.
#   4. Walks the per-cluster FASTAs the container produced and prints a
#      tab-separated table to stdout (or to -o PATH).
#
# Two input modes:
#   -d INPUT_DIR        process every *.nearby.sam in INPUT_DIR
#   SAM [SAM ...]       process the given *.nearby.sam files explicitly
#                       (each must have its matching *.nearby.readnamescoords
#                        file alongside it)
#
# Output table (tab-separated):
#   cluster    side     length    sequence
# - polish mode: side = "junction", one row per cluster.
# - pileup mode: side = "5p" or "3p", two rows per cluster
#   (lowercase letters mark non-reference consensus calls).
#
# Run with --help (or -h) for the usage block plus an example.

set -euo pipefail

# ---- defaults --------------------------------------------------------------

MODE=polish                      # which container script to run
WINDOW=                          # half-window in bp; mode-specific default below
OUT=-                            # "-" means stdout; otherwise a file path
SEARCH_DIR=                      # set by -d; overrides positional SAM list
KEEP_DIR=                        # set by -D; preserves per-cluster FASTAs
IMAGE=biodepot/bff:test          # container tag; -i overrides

# Locate the in-repo Dockerfiles dir relative to this script. Used only
# if the image is missing and we need to build it. If the script is
# moved out of the repo (eg /usr/local/bin), this path won't resolve and
# the script will refuse to auto-build (telling the user how instead).
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DOCKERFILES_DIR="$SCRIPT_DIR/../widgets/fast_fusion_finder/BiodepotFusionFinder/Dockerfiles"

# ---- help ------------------------------------------------------------------

usage() {
    # Always written to stderr so the TSV stdout stays clean if someone
    # accidentally pipes runFusionConsensus.sh -h | column -t.
    cat <<'EOF' >&2
Usage:
  runFusionConsensus.sh [opts] -d INPUT_DIR REFERENCE.fa
  runFusionConsensus.sh [opts]    REFERENCE.fa SAM [SAM ...]

Options:
  -m polish|pileup   pipeline (default: polish)
                       polish: chimeric-ref + minimap2 + racon
                               -> one primer-grade record per cluster
                       pileup: per-side mpileup vs the genome
                               -> two records (5'/3') per cluster,
                                  lowercase letters = non-reference calls
  -w N               half-window size in bp
                       defaults: 300 (polish) / 100 (pileup)
                       polish needs more anchor for minimap2 to map
                       noisy long reads onto the small chimeric ref.
  -o PATH            write the TSV table to PATH instead of stdout
  -D DIR             keep the per-cluster FASTAs (and pileup .tsv files)
                     in DIR after the run; DIR is created if missing.
                     Default: per-cluster artefacts are dropped, only
                     the table survives.
  -d DIR             process every *.nearby.sam under DIR.
                     Mutually exclusive with passing SAM files
                     positionally.
  -i IMAGE           override the docker image tag
                     (default: biodepot/bff:test)
  -h, --help         show this help and exit

Output (TSV, on stdout unless -o):
  cluster<TAB>side<TAB>length<TAB>sequence
    polish: side = "junction", one row per cluster
    pileup: side = "5p" / "3p", two rows per cluster

Example:
  # Polish every cluster in a breakpoint_files/ dir using the
  # cellranger Ensembl reference; redirect the table to a file.
  bin/runFusionConsensus.sh \
      -d /mnt/pikachu/bff-test/breakpoint_files \
      -o consensus.tsv \
      /mnt/pikachu/autoindex_110_44/cellranger_ref_cache/Homo_sapiens.GRCh38.dna.primary_assembly.fa

Notes:
  * Only docker is required on the host; the image carries samtools,
    perl, minimap2, racon, and the wrapper scripts.
  * The reference FASTA must have a .fai (the in-container wrappers
    will create one lazily if its directory is writable).
  * The SAM files passed positionally must end in .nearby.sam and have
    their .nearby.readnamescoords companion next to them.
EOF
    exit "${1:-1}"
}

# Translate --help into -h before getopts sees it. getopts only handles
# single-char short options.
for arg in "$@"; do
    [[ "$arg" == "--help" ]] && usage 0
done

# ---- option parsing --------------------------------------------------------

while getopts "m:w:o:d:D:i:h" opt; do
    case "$opt" in
        m) MODE=$OPTARG;;        # polish | pileup
        w) WINDOW=$OPTARG;;      # half-window override
        o) OUT=$OPTARG;;         # TSV destination
        d) SEARCH_DIR=$OPTARG;;  # directory-mode input
        D) KEEP_DIR=$OPTARG;;    # preserve per-cluster artefacts
        i) IMAGE=$OPTARG;;       # custom image tag
        h) usage 0;;
        *) usage;;
    esac
done
shift $((OPTIND - 1))

# Pick the in-container wrapper script and apply mode-specific window.
case "$MODE" in
    polish) WRAPPER=polishFusionConsensus.sh; WINDOW=${WINDOW:-300};;
    pileup) WRAPPER=processBreakpointDir.sh;  WINDOW=${WINDOW:-100};;
    *) echo "unknown mode: $MODE (use polish or pileup)" >&2; exit 1;;
esac

# ---- preflight: docker + image --------------------------------------------

command -v docker >/dev/null 2>&1 \
    || { echo "docker not found on PATH" >&2; exit 1; }

# Build the image if missing AND we know where the Dockerfiles live.
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

# ---- cleanup tracking ------------------------------------------------------

# We may create up to two tempdirs (an input staging dir for explicit
# SAM-file mode, and an output dir if -D wasn't given). Track them so
# the EXIT trap can clean them up after the table has been emitted.
cleanup_dirs=()
cleanup() {
    local d
    for d in "${cleanup_dirs[@]+"${cleanup_dirs[@]}"}"; do
        [[ -n "$d" && -d "$d" ]] && rm -rf "$d"
    done
}
trap cleanup EXIT

# ---- resolve inputs --------------------------------------------------------

if [[ -n "$SEARCH_DIR" ]]; then
    # Directory mode: -d INPUT_DIR REFERENCE.fa
    [[ $# -lt 1 ]] && usage
    REF=$1; shift
    [[ $# -gt 0 ]] && { echo "extra args after reference: $*" >&2; usage; }
    [[ -d "$SEARCH_DIR" ]] || { echo "directory not found: $SEARCH_DIR" >&2; exit 1; }
    INPUT_DIR=$(realpath "$SEARCH_DIR")
else
    # Explicit SAM mode: REFERENCE.fa SAM [SAM ...]
    [[ $# -lt 2 ]] && usage
    REF=$1; shift

    # Stage the SAM/coords pair for each cluster into a tempdir so the
    # in-container wrapper sees a well-shaped directory. Hardlink first
    # (no copy when on the same FS), fall back to cp across filesystems.
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
        ln "$sam"    "$INPUT_DIR/${base}.nearby.sam"             2>/dev/null \
            || cp "$sam"    "$INPUT_DIR/${base}.nearby.sam"
        ln "$coords" "$INPUT_DIR/${base}.nearby.readnamescoords" 2>/dev/null \
            || cp "$coords" "$INPUT_DIR/${base}.nearby.readnamescoords"
    done
fi

# Reference: must exist; resolve to absolute path so docker -v sees it.
[[ -f "$REF" ]] || { echo "reference not found: $REF" >&2; exit 1; }
REF=$(realpath "$REF")
REF_DIR=$(dirname "$REF")
REF_NAME=$(basename "$REF")

# ---- per-cluster artefact dir ---------------------------------------------

# OUT_DIR holds the .fa (and, in pileup mode, .tsv) files the wrapper
# emits. With -D we keep them; otherwise we make a tempdir that the
# EXIT trap removes after we've built the table from them.
if [[ -n "$KEEP_DIR" ]]; then
    mkdir -p "$KEEP_DIR"
    OUT_DIR=$(realpath "$KEEP_DIR")
else
    OUT_DIR=$(mktemp -d)
    cleanup_dirs+=("$OUT_DIR")
fi

# ---- run the in-container pipeline ----------------------------------------

# All wrapper progress goes to stderr so stdout is reserved for the TSV.
# -u keeps the host user as the owner of files written into OUT_DIR.
docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$INPUT_DIR:/in:ro" \
    -v "$REF_DIR:/ref:ro" \
    -v "$OUT_DIR:/out" \
    "$IMAGE" \
    "$WRAPPER" /in "/ref/$REF_NAME" /out "$WINDOW" >&2

# ---- build and emit the TSV table -----------------------------------------

# "-" -> /dev/stdout so the same { ... } > "$OUT" works either way.
[[ "$OUT" == "-" ]] && OUT=/dev/stdout

{
    printf 'cluster\tside\tlength\tsequence\n'

    # Iterate over each cluster's FASTA. Each file has 1 record (polish)
    # or 2 records (pileup, 5'/3'). awk concatenates the wrapped sequence
    # lines into a single string and emits one TSV row per record.
    shopt -s nullglob
    for fa in "$OUT_DIR"/*.fa; do
        base=$(basename "$fa" .fa)
        awk -v cluster="$base" -v mode="$MODE" '
            BEGIN { rec = 0; seq = ""; have = 0 }

            # Emit the most-recently-collected record (called when we
            # see the next ">" or hit EOF).
            function emit() {
                side = (mode == "polish") ? "junction" : (rec == 1 ? "5p" : "3p")
                printf "%s\t%s\t%d\t%s\n", cluster, side, length(seq), seq
            }

            /^>/ {
                if (have) emit()
                have = 1; rec++; seq = ""
                next
            }

            # All non-header lines are sequence; concatenate.
            { seq = seq $0 }

            END { if (have) emit() }
        ' "$fa"
    done
} > "$OUT"
