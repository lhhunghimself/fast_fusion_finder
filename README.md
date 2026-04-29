# Fusion breakpoint consensus

Build per-base consensus sequences over the supporting reads for every
fusion breakpoint detected by the BFF pipeline. The two halves of each
fusion are emitted as **separate FASTA records** (5' partner, then 3'
partner), windowed around the junction; columns where a non-reference
base wins a strict majority of the reads are emitted in lowercase.

Two entry points:

- [`bin/breakpointPileup.pl`](bin/breakpointPileup.pl) — the core script.
  Consumes a read-names list, the original BAM, the reference FASTA, and
  a breakpoint file; produces one FASTA record per side and (optionally)
  a per-position TSV of base counts. Use this when you already have
  those four pieces.
- [`bin/processBreakpointDir.sh`](bin/processBreakpointDir.sh) — batch
  wrapper for a `breakpoint_files/`-style directory. For each cluster it
  rebuilds the SAM into a sorted BAM, derives the reads list and the
  breakpoint encoding from `*.nearby.readnamescoords`, and runs the core
  script. Use this when you have BFF cluster output and a reference.

## Layout

```
bin/
├── breakpointPileup.pl       canonical script
├── processBreakpointDir.sh   batch wrapper
├── filterEnds.pl             upstream: read-name filter
└── ...
breakpointPileup.pl           symlink → bin/breakpointPileup.pl
tests/test_breakpointPileup.sh   self-contained fixture test
widgets/.../Dockerfiles/      copies that get baked into biodepot/bff:test
```

`breakpointPileup.pl` lives at `bin/breakpointPileup.pl`; the top-level
file is a symlink for convenience. The Dockerfiles directory holds a
build-context copy of each script (kept in sync manually).

## Requirements

There is nothing to compile.

- `samtools` 1.12 or newer on `$PATH` (1.12 added `samtools view -N`,
  the read-name filter the script depends on)
- `perl` 5 with `Getopt::Long` and `File::Temp` (both core)
- the reference FASTA must have a `.fai` (built lazily by the wrapper
  when missing, but only if the directory is writable)

Verify:

```bash
samtools --version | head -1   # 1.12+
perl -c bin/breakpointPileup.pl
```

## Three ways to run

### 1. Container (zero host setup)

The `biodepot/bff:test` image bundles `samtools`, `perl`, and both
scripts under `/usr/local/bin`.

```bash
# build once
docker build -t biodepot/bff:test \
    widgets/fast_fusion_finder/BiodepotFusionFinder/Dockerfiles

# run on a breakpoint_files-style directory
docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v /path/to/breakpoint_files:/in:ro \
    -v /path/to/reference_dir:/ref:ro \
    -v /path/to/output:/out \
    biodepot/bff:test \
    processBreakpointDir.sh /in /ref/genome.fa /out 100
```

The `-u` flag makes the output files owned by your host user instead of
root. The reference directory needs the `.fa.fai`; create it once on
the host with `samtools faidx genome.fa` if it isn't there.

### 2. Standalone

If `samtools` and `perl` are on the host, just put the repo's `bin/`
on `$PATH` (so the wrapper can find `breakpointPileup.pl`) and run:

```bash
PATH="$PWD/bin:$PATH" \
bin/processBreakpointDir.sh \
    /path/to/breakpoint_files \
    /path/to/genome.fa \
    /path/to/output \
    100
```

### 3. As part of the BFF widget pipeline

`runBFF.sh` (the widget entrypoint) calls `breakpointPileup.pl`
automatically when both `$BAM` and `$REFERENCE` are set:

```bash
WINDOW=100 \
pileup=/data/fusion.fa \
pileupTable=/data/variants.tsv \
filterends=/data/reads.txt \
BAM=/data/all.sorted.bam REFERENCE=/data/genome.fa \
runBFF.sh --first first.bam --last last.bam -b breakpoints.txt
```

If `$BAM` or `$REFERENCE` is empty the pileup step is skipped and
`runBFF.sh` runs `filterEnds.pl` alone, matching the legacy behaviour.

## Inputs

### For `processBreakpointDir.sh`

| Path | Purpose |
|------|---------|
| `INPUT_DIR` | A `breakpoint_files/`-style directory with `*.nearby.sam` (alignment dump, no `@HD`/`@SQ`) and matching `*.nearby.readnamescoords`. |
| `REFERENCE_FASTA` | The genome FASTA used for alignment. Needs a `.fai`. |
| `OUTPUT_DIR` | Created if missing. One `<name>.fa` and one `<name>.tsv` per cluster. |
| `WINDOW` | Optional half-window size in bp (default `100`). Each side gets `WINDOW` bases on the transcript-internal flank of the breakpoint. |

The wrapper auto-detects UCSC vs Ensembl contig naming. If the SAM uses
`chr15`-style names but the reference doesn't, `chr` is stripped on the
fly from `RNAME`, `RNEXT`, and the breakpoint encoding before mpileup.

### For `breakpointPileup.pl`

| Flag | Required | Notes |
|------|----------|-------|
| `--reads reads.txt` | yes | Newline-separated read names — typically the stdout of `filterEnds.pl`. |
| `--bam in.bam` | yes | Coordinate-sorted and indexed. |
| `--ref ref.fa` | yes | With matching `.fai`. |
| `-b breakpoints.txt` | yes | One fusion per line: `chr:start:end;chr:start:end`. Strand encoded by coord order — `start<end` is `+`, `start>end` is `-`. Lines starting with `#` are ignored. |
| `--window N` | no, default 100 | bp per side. |
| `--pileupOut path` | no | Write a per-position TSV (chr, pos, ref, consensus, depth, A, C, G, T) alongside stdout. |
| `--threshold F` | no, default 0.5 | Strict majority needed to call a non-reference base. |
| `--minDepth D` | no, default 1 | Columns below this depth fall back to the reference. |
| `--guideFile guides.txt` (`-G`) | no | `group<TAB>chr:pos<TAB>+/-`. Overrides the breakpoint-file's encoded strand and groups output by guide name. |

## Output

For each non-comment line of the breakpoint file, two FASTA records
are written to stdout (or to a file if you redirect):

```
>15:74024123-74024222+ length=100
GAGAGTCTCCAGGAGTCTTTGACTAAATCTCCAAGCTGGAATGTGAGCTCTGAGCAGCTAGTATAGGATGCCCTTACTGG
GAAAGGGGAGGGAGGCTATG
>17:40346374-40346473+ length=100
GAGAGGCGGGGCCCAGGGCAAACGGTGGATTAGGAGGGGTGGGGAGGTCAGTGCCTTCTTCCTCTGCTTGTCGGAATGCT
GACCAAGATTCTAGGCCATG
```

Headers are `<region><strand> length=<N>`; sides are emitted 5' first
then 3', in the order they appear on the breakpoint line.
Lowercase letters mark columns where a non-reference base won
`> --threshold` of the parsed reads. The strand annotation in the
header reflects what was decoded from the breakpoint file's coord
order (or, with `--guideFile`, the guide direction).

When `--pileupOut` is set, the companion TSV has one row per pileup
column **in genomic order** (so for `-` strand sides the row order
is the reverse of the FASTA letter order):

```
#chr   pos    ref  consensus  depth  A  C  G  T
15     74024124  G  G   2  0  0  2  0
15     74024125  A  a   8  0  0  6  2
...
```

## Conventions worth knowing

- **Strand from coord order.** The upstream BFF generator emits each
  side with `start<end` for `+` and `start>end` for `-`. The script
  trusts that encoding — no fallback to BAM FLAGs. If you supply your
  own breakpoint file, follow the same convention or use `--guideFile`.
- **Breakpoint position per side.** Side 1 (5' partner) → second coord
  of the encoded `chr:start:end`. Side 2 (3' partner) → first coord.
  This matches what `findBreakPoints` in `filterEnds.pl` writes.
- **Window placement.** `5'+` and `3'-` use `[bp-W+1, bp]`; `5'-` and
  `3'+` use `[bp, bp+W-1]`. So each FASTA record is `W` bases of the
  partner gene's transcript-side flank ending at (or starting at) the
  junction.
- **Variant case.** Lowercase = non-reference call. Reverse complement
  preserves case (`tr/ACGTNacgtn/TGCANtgcan/`).

## Smoke test

`tests/test_breakpointPileup.sh` builds a synthetic chr1/chr2 reference
and ten reads, then exercises forward+forward, narrow window,
forward+reverse with RC variants, and zero-coverage fallback:

```bash
bash tests/test_breakpointPileup.sh
```

All assertions should pass. The same suite runs unmodified inside the
container against `/usr/local/bin/breakpointPileup.pl`.
