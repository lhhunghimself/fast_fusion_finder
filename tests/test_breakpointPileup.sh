#!/usr/bin/env bash
# Self-contained smoke test for bin/breakpointPileup.pl.
#
# Builds a tiny synthetic reference + BAM and exercises:
#   1. forward+forward pair with no variants
#   2. forward+forward pair with a majority-rule variant on side 1
#   3. forward+reverse pair, with a majority-rule variant on the reverse side
#      (so we can verify RC preserves variant lowercase and ends up in the
#       right transcript-order position)
#   4. tightening --window further than the breakpoint range and checking the
#      output anchors on the breakpoint position, not the range bounds
#   5. that --pileupOut emits a per-position TSV with correct counts
#
# Run with `bash tests/test_breakpointPileup.sh` from the repo root.
# Requires samtools 1.12+ on PATH.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PILEUP="$REPO_ROOT/bin/breakpointPileup.pl"

if [[ ! -x "$PILEUP" ]]; then
    echo "FAIL: $PILEUP is not executable" >&2
    exit 1
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# ---- helpers -----------------------------------------------------------------

repeat() { # repeat CHAR N
    local ch="$1" n="$2"
    head -c "$n" </dev/zero | tr '\0' "$ch"
}

assert_eq() { # assert_eq LABEL EXPECTED ACTUAL
    local label="$1" expected="$2" actual="$3"
    if [[ "$expected" == "$actual" ]]; then
        printf 'PASS  %s\n' "$label"
    else
        printf 'FAIL  %s\n  expected: %s\n  got:      %s\n' "$label" "$expected" "$actual"
        exit 1
    fi
}

fasta_seq() { # fasta_seq PATH HEADER_NUM(1-based)  -> seq with no newlines
    awk -v want="$2" 'BEGIN{n=0} /^>/{n++; next} n==want{printf "%s", $0}' "$1"
}

# ---- fixtures ----------------------------------------------------------------

# chr1 = 200 A, chr2 = 200 G. Picking different bases per chromosome makes RC
# unambiguous (RC of all-G is all-C, distinguishable from chr1's A).
{
    printf '>chr1\n%s\n' "$(repeat A 200)"
    printf '>chr2\n%s\n' "$(repeat G 200)"
} >ref.fa
samtools faidx ref.fa

QUAL50=$(repeat I 50)
A50=$(repeat A 50)
A49=$(repeat A 49)
G50=$(repeat G 50)
G49=$(repeat G 49)

emit_sam() {
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    printf '@SQ\tSN:chr1\tLN:200\n'
    printf '@SQ\tSN:chr2\tLN:200\n'

    # 10 reads on chr1 covering 101-150. 8 carry a T at position 150 (last
    # base), 2 carry the reference A.
    for i in 1 2 3 4 5 6 7 8; do
        printf 'chr1var%d\t0\tchr1\t101\t60\t50M\t*\t0\t0\t%sT\t%s\n' "$i" "$A49" "$QUAL50"
    done
    for i in 9 10; do
        printf 'chr1ref%d\t0\tchr1\t101\t60\t50M\t*\t0\t0\t%s\t%s\n' "$i" "$A50" "$QUAL50"
    done

    # 10 reads on chr2 covering 50-99 with no variants.
    for i in 1 2 3 4 5 6 7 8 9 10; do
        printf 'chr2plus%d\t0\tchr2\t50\t60\t50M\t*\t0\t0\t%s\t%s\n' "$i" "$G50" "$QUAL50"
    done

    # 10 reads on chr2 covering 100-149. 8 carry a T at position 100 (first
    # base), 2 carry the reference G. Used by the forward+reverse test.
    for i in 1 2 3 4 5 6 7 8; do
        printf 'chr2var%d\t0\tchr2\t100\t60\t50M\t*\t0\t0\tT%s\t%s\n' "$i" "$G49" "$QUAL50"
    done
    for i in 9 10; do
        printf 'chr2ref%d\t0\tchr2\t100\t60\t50M\t*\t0\t0\t%s\t%s\n' "$i" "$G50" "$QUAL50"
    done
}

emit_sam >reads.sam
samtools sort -O bam reads.sam -o reads.bam
samtools index reads.bam

# Reads list — every synthetic read participates.
{
    for i in 1 2 3 4 5 6 7 8;  do echo "chr1var$i";   done
    for i in 9 10;             do echo "chr1ref$i";   done
    for i in 1 2 3 4 5 6 7 8 9 10; do echo "chr2plus$i"; done
    for i in 1 2 3 4 5 6 7 8;  do echo "chr2var$i";   done
    for i in 9 10;             do echo "chr2ref$i";   done
} >reads.txt

# ---- test 1: forward+forward, variant on side 1 ------------------------------

echo "chr1:51:150;chr2:50:149" >bp1.txt
"$PILEUP" --reads reads.txt --bam reads.bam --ref ref.fa -b bp1.txt \
          --window 50 --pileupOut t1.tsv >t1.fa

EXP1A="$(repeat A 49)t"       # side 1 (5'): 49 As + lowercase t variant
EXP1B="$G50"                  # side 2 (3'): 50 Gs
GOT1A=$(fasta_seq t1.fa 1)
GOT1B=$(fasta_seq t1.fa 2)
assert_eq "test1 fwd side-1 seq with variant" "$EXP1A" "$GOT1A"
assert_eq "test1 fwd side-2 seq" "$EXP1B" "$GOT1B"

# TSV at chr1:150 should show consensus 't' with A=2, T=8.
ROW=$(awk '$1=="chr1" && $2==150' t1.tsv)
assert_eq "test1 TSV chr1:150" "chr1	150	A	t	10	2	0	0	8" "$ROW"

# Sanity: a non-variant column.
ROW=$(awk '$1=="chr2" && $2==50' t1.tsv)
assert_eq "test1 TSV chr2:50 (no variant)" "chr2	50	G	G	10	0	0	10	0" "$ROW"

# ---- test 2: tighter --window anchors on the breakpoint, not the range -------

# Same breakpoint file, --window 5 -> sides should be chr1:146-150 + chr2:50-54.
"$PILEUP" --reads reads.txt --bam reads.bam --ref ref.fa -b bp1.txt \
          --window 5 >t2.fa

EXP2A="AAAAt"
EXP2B="$(repeat G 5)"
assert_eq "test2 narrow window side 1" "$EXP2A" "$(fasta_seq t2.fa 1)"
assert_eq "test2 narrow window side 2" "$EXP2B" "$(fasta_seq t2.fa 2)"

# ---- test 3: forward+reverse, variant on the reverse side --------------------

# Reverse side 2 by encoding it start>end. Side 2 breakpoint = start = 149,
# pileup window [100, 149], reads have 8/10 T at chr2:100. After RC the
# variant should appear at the *last* position of the fusion (RC moves
# position 100 to the tail and complements T -> a).
echo "chr1:51:150;chr2:149:100" >bp3.txt
"$PILEUP" --reads reads.txt --bam reads.bam --ref ref.fa -b bp3.txt \
          --window 50 --pileupOut t3.tsv >t3.fa

EXP3A="${A49}t"                 # side 1 (5' fwd): 49 As + t variant
EXP3B="$(repeat C 49)a"         # side 2 (3' rev): 49 Cs + a (variant RC'd to lowercase a)
assert_eq "test3 fwd side-1 seq" "$EXP3A" "$(fasta_seq t3.fa 1)"
assert_eq "test3 rev side-2 seq RC preserves variant case" "$EXP3B" "$(fasta_seq t3.fa 2)"

# Pileup TSV stays in genomic order, so chr2:100 row should record 't' alt.
ROW=$(awk '$1=="chr2" && $2==100' t3.tsv)
assert_eq "test3 TSV chr2:100 (variant on rev side, genomic order)" \
          "chr2	100	G	t	10	0	0	2	8" "$ROW"

# Each side's FASTA header should reflect its resolved direction.
H1=$(awk '/^>/{print; if (++n==1) exit}' t3.fa)
H2=$(awk '/^>/{n++; if (n==2) {print; exit}}' t3.fa)
[[ "$H1" == *"chr1:101-150+"* ]] || { echo "FAIL  test3 header 1: $H1"; exit 1; }
[[ "$H2" == *"chr2:100-149-"* ]] || { echo "FAIL  test3 header 2: $H2"; exit 1; }
echo "PASS  test3 per-side FASTA headers record strands"

# ---- test 4: zero-coverage falls back to reference ---------------------------

# Empty reads list -> every column has depth 0 -> output should be the plain
# reference (no lowercase, no RC effect on case).
:>empty_reads.txt
"$PILEUP" --reads empty_reads.txt --bam reads.bam --ref ref.fa -b bp1.txt \
          --window 10 >t4.fa
EXP4A="$(repeat A 10)"
EXP4B="$(repeat G 10)"
assert_eq "test4 zero-coverage side 1 = ref" "$EXP4A" "$(fasta_seq t4.fa 1)"
assert_eq "test4 zero-coverage side 2 = ref" "$EXP4B" "$(fasta_seq t4.fa 2)"

echo
echo "all tests passed"
