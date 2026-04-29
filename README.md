# breakpointPileup.pl

A standalone consensus-sequence builder for fusion breakpoints. Takes the
read-names list produced by `filterEnds.pl`, plus the original BAM and
reference, and emits **two FASTA records per breakpoint line** — one for
each side of the fusion, in 5'→3' order — each a per-position consensus
over the supporting reads (reference base by default, replaced when a
single non-reference base is seen in more than `--threshold` of the read
bases at that column). Sides whose breakpoint coordinates encode `-`
strand are reverse-complemented before output.

The canonical script lives at `bin/breakpointPileup.pl`. The top-level
`breakpointPileup.pl` is a symlink to it for convenient ad-hoc testing.
A separate copy at
`widgets/fast_fusion_finder/BiodepotFusionFinder/Dockerfiles/breakpointPileup.pl`
is what gets ADDed into the BFF Docker image — keep it in sync with `bin/`
when you change either one.

## Requirements

There is nothing to compile — it is a Perl script that shells out to
`samtools`. You need:

- Perl 5 (any modern distribution; the script uses only `Getopt::Long`,
  which is core)
- `samtools` 1.12 or newer on `$PATH` (1.12 added `samtools view -N`,
  the read-name filter this script depends on)

Verify with:

```bash
perl -v
samtools --version | head -1   # 1.12+
perl -c ./breakpointPileup.pl  # syntax check, prints "syntax OK"
```

## Inputs

1. **`--reads reads.txt`** — newline-separated read names. The expected
   producer is `bin/filterEnds.pl`, whose stdout is exactly this format.
2. **`--bam in.bam`** — the original alignment. Must be coordinate-sorted
   and indexed (`samtools sort` + `samtools index`).
3. **`--ref ref.fa`** — the reference FASTA used for alignment. Must have
   a matching `.fai` (`samtools faidx ref.fa`).
4. **`--breakpointFile breakpoints.txt`** (`-b`) — the same breakpoint file
   `filterEnds.pl` was given. Each non-comment line is one fusion in the
   form `chr:start:end;chr:start:end`. Strand is encoded by coordinate
   order: `start < end` is `+`, `start > end` is `-`. Lines starting with
   `#` are ignored.

Optional:

- **`--guideFile guides.txt`** (`-G`) — overrides the breakpoint-file's
  encoded strand. Format is the same as the input to
  `guidesToBreakpoints.pl`: `group<TAB>chr:pos<TAB>+/-`. Each group's
  halves are matched to breakpoint regions by chromosome+position,
  reverse-complemented per the guide direction, and concatenated in the
  order the entries appear in the file (first listed = 5' partner).
- **`--threshold 0.5`** (`-t`) — fraction of read bases at a column needed
  to override the reference base. Strict majority (`>`, not `>=`).
- **`--minDepth 1`** (`-d`) — columns with fewer reads than this fall
  back to the reference base regardless of what was observed.
- **`--verbose`** (`-v`) — echo the `samtools` invocations to stderr.

## Running it

End-to-end from a sorted+indexed BAM:

```bash
# 1. Get supporting read names for each breakpoint.
bin/filterEnds.pl --first first.bam --last last.bam \
                  -b breakpoints.txt > reads.txt

# 2. Build the fusion FASTA.
./breakpointPileup.pl --reads reads.txt \
                     --bam   all.bam \
                     --ref   ref.fa \
                     -b      breakpoints.txt \
                     > fusion.fa
```

Two FASTA records per non-comment line in `breakpoints.txt` — the 5'
partner first, then the 3' partner. Each header carries that side's
region and resolved strand, e.g.
`>chr15:74318559-74336132+ length=...` followed by
`>chr17:38428464-38512385- length=...`.

## Quick smoke test

```bash
samtools sort -o test.sorted.bam test.bam
samtools index test.sorted.bam
samtools faidx ref.fa

# Minimal breakpoint file, one fusion:
printf 'chr15:74318559:74336132;chr17:38428464:38512385\n' > bp.txt

bin/filterEnds.pl --first test.sorted.bam --last test.sorted.bam \
                  -b bp.txt > reads.txt
./breakpointPileup.pl --reads reads.txt --bam test.sorted.bam \
                     --ref ref.fa -b bp.txt
```

If `reads.txt` is empty, every consensus column will fall back to the
reference base, so the output FASTA is the reference sequence over each
region (RC'd where the breakpoint encodes `-`). That's the cheapest way
to confirm the plumbing without needing real fusion-supporting reads.
