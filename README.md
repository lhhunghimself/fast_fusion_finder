# Fusion breakpoint consensus

Build per-base consensus sequences over the supporting reads for every
fusion breakpoint detected by the BFF pipeline. There are two output
modes depending on what you want the consensus *for*:

- **Per-side mpileup consensus** (`processBreakpointDir.sh`,
  `breakpointPileup.pl`) — fast, reference-anchored, emits the two
  halves of each fusion as **separate FASTA records** (5' partner, then
  3' partner) windowed around the junction. Columns where a
  non-reference base wins a strict majority are emitted in lowercase.
  Good for variant inspection and quick QC.
- **Polished chimeric consensus** (`polishFusionConsensus.sh`) — builds
  a chimeric reference per fusion (5' window + 3' window in transcript
  order), realigns the supporting reads to it with `minimap2 -ax map-ont`,
  then polishes with `racon`. Output is a single continuous record
  spanning the junction with indels resolved at the read level.
  **Use this for primer design.**

Entry points:

- [`bin/runFusionConsensus.sh`](bin/runFusionConsensus.sh) — **top-level
  driver. Build the Docker image lazily, run either pipeline inside it,
  emit a TSV table.** Accepts either a directory of `*.nearby.sam` files
  (`-d DIR`) or an explicit list of SAM paths. Use this when you just
  want a sequence table.
- [`bin/breakpointPileup.pl`](bin/breakpointPileup.pl) — the per-side
  pileup core script. Consumes a read-names list, the original BAM, the
  reference FASTA, and a breakpoint file; produces one FASTA record per
  side and (optionally) a per-position TSV of base counts.
- [`bin/processBreakpointDir.sh`](bin/processBreakpointDir.sh) — batch
  wrapper around `breakpointPileup.pl` for a `breakpoint_files/`-style
  directory. Two records per cluster, mpileup-derived.
- [`bin/polishFusionConsensus.sh`](bin/polishFusionConsensus.sh) — the
  chimeric-realign + racon-polish wrapper. One polished record per
  cluster, primer-design-grade.

## Layout

```
bin/
├── runFusionConsensus.sh      top-level driver (docker, builds image lazily)
├── breakpointPileup.pl        per-side mpileup core script
├── processBreakpointDir.sh    batch wrapper around breakpointPileup.pl
├── polishFusionConsensus.sh   chimeric-realign + racon polish wrapper
├── extractReadsFasta.sh       names + FASTQ → FASTA
├── filterEnds.pl              upstream: read-name filter
└── ...
breakpointPileup.pl            symlink → bin/breakpointPileup.pl
tests/test_breakpointPileup.sh self-contained fixture test
widgets/.../Dockerfiles/       copies that get baked into biodepot/bff:test
```

`breakpointPileup.pl` lives at `bin/breakpointPileup.pl`; the top-level
file is a symlink for convenience. The Dockerfiles directory holds a
build-context copy of each script (kept in sync manually).

## Requirements

There is nothing to compile.

For the mpileup pipeline (`processBreakpointDir.sh`,
`breakpointPileup.pl`):

- `samtools` 1.12 or newer on `$PATH` (1.12 added `samtools view -N`,
  the read-name filter the script depends on)
- `perl` 5 with `Getopt::Long` and `File::Temp` (both core)
- the reference FASTA must have a `.fai` (built lazily by the wrapper
  when missing, but only if the directory is writable)

For the polish pipeline (`polishFusionConsensus.sh`), additionally:

- `minimap2` 2.x on `$PATH` (the `map-ont` preset is used)
- `racon` 1.5+ on `$PATH`

The `biodepot/bff:test` Docker image ships with all four — easiest way
to get the polish pipeline running without touching your host.

Verify:

```bash
samtools --version | head -1   # 1.12+
perl -c bin/breakpointPileup.pl
```

## Quick start: the top-level driver

For most users, `bin/runFusionConsensus.sh` is the only command needed.
It builds `biodepot/bff:test` if it isn't there, mounts the inputs into
the container, runs the chosen pipeline, and prints a tab-separated
table of cluster sequences. Run with `-h` (or `--help`) to see all
flags and a copy-pasteable example:

```bash
bin/runFusionConsensus.sh --help
```

### Common usages

```bash
# 1. Polish every *.nearby.sam in a breakpoint_files-style directory.
#    Default mode; one primer-grade row per cluster.
bin/runFusionConsensus.sh \
    -d /path/to/breakpoint_files \
    /path/to/genome.fa \
    > consensus.tsv

# 2. Pass SAM files directly (each .nearby.sam must have its
#    matching .nearby.readnamescoords next to it).
bin/runFusionConsensus.sh \
    /path/to/genome.fa \
    /path/to/cluster1.nearby.sam /path/to/cluster2.nearby.sam \
    > consensus.tsv

# 3. Per-side mpileup mode for variant inspection
#    (two rows per cluster, lowercase = non-reference call).
bin/runFusionConsensus.sh -m pileup -w 100 \
    -d /path/to/breakpoint_files \
    /path/to/genome.fa \
    > sides.tsv

# 4. Keep the per-cluster FASTAs around (default: dropped after the
#    table is built). Useful when you want the .fa / .tsv files too.
bin/runFusionConsensus.sh -D out/ \
    -d /path/to/breakpoint_files \
    /path/to/genome.fa \
    > consensus.tsv

# 5. Dry-run: list the clusters that would be processed and exit
#    without invoking docker. Useful for sanity-checking inputs.
bin/runFusionConsensus.sh --dry-run \
    -d /path/to/breakpoint_files \
    /path/to/genome.fa
```

### Worked example (the `bff-test` fixtures)

Against the two 2-read clusters in
`/mnt/pikachu/bff-test/breakpoint_files` and the cellranger Ensembl
FASTA:

```bash
bin/runFusionConsensus.sh \
    -d /mnt/pikachu/bff-test/breakpoint_files \
    -o consensus.tsv \
    /mnt/pikachu/autoindex_110_44/cellranger_ref_cache/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

`consensus.tsv`:

```
cluster                          side      length  sequence
chr15-74024222-chr17-40346374    junction  626     TTGCGCGCCCTGCCTGAGTGTTGGAC...
chr15-74024223-chr17-40346371    junction  621     TGCGCGCCCTGCCTGAGTGTTGGACG...
```

Add `-m pileup -w 50` to the same command to get four rows instead
(5' and 3' for each cluster, with lowercase variants where they
appear).

### Output columns

| Column | Meaning |
|--------|---------|
| `cluster` | Cluster name, parsed from the source filename (`chrA-posA-chrB-posB`). |
| `side` | `junction` in polish mode; `5p` / `3p` in pileup mode. |
| `length` | Sequence length in bp. |
| `sequence` | The consensus, on a single line. Lowercase letters in pileup mode mark non-reference consensus calls; polish mode is plain uppercase. |

The driver only needs `docker` on the host. The image carries
`samtools`, `perl`, `minimap2`, `racon`, and the wrapper scripts —
nothing else has to be installed.

## Three ways to run the underlying pipeline

### 1. Container (zero host setup)

The `biodepot/bff:test` image bundles `samtools`, `perl`, `minimap2`,
`racon`, and the wrappers under `/usr/local/bin`.

```bash
# build once
docker build -t biodepot/bff:test \
    widgets/fast_fusion_finder/BiodepotFusionFinder/Dockerfiles

# Mode A — fast per-side mpileup consensus (variant-aware, two records per cluster)
docker run --rm -u "$(id -u):$(id -g)" \
    -v /path/to/breakpoint_files:/in:ro \
    -v /path/to/reference_dir:/ref:ro \
    -v /path/to/output:/out \
    biodepot/bff:test \
    processBreakpointDir.sh /in /ref/genome.fa /out 100

# Mode B — polished chimeric consensus for primer design (one record per cluster)
docker run --rm -u "$(id -u):$(id -g)" \
    -v /path/to/breakpoint_files:/in:ro \
    -v /path/to/reference_dir:/ref:ro \
    -v /path/to/output:/out \
    biodepot/bff:test \
    polishFusionConsensus.sh /in /ref/genome.fa /out 300
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

### For `polishFusionConsensus.sh`

Same arguments as above, but the default `WINDOW` is **300** because
minimap2 needs a longer chimeric contig to anchor reads cleanly. Per
cluster the wrapper:

1. Builds a chimeric reference by extracting `WINDOW` bases of the 5'
   partner ending at the breakpoint and `WINDOW` bases of the 3'
   partner starting at the breakpoint (RC where the breakpoint encoding
   says `-`), concatenated as one record.
2. Pulls the cluster's supporting reads' sequences from the SAM via
   `samtools fasta` (in their original orientation, FLAG-aware).
3. Realigns those reads to the chimeric ref with
   `minimap2 -ax map-ont --secondary=no`.
4. Polishes with `racon`, writing one FASTA record per cluster to
   `OUTPUT_DIR/<name>.fa`. Racon's header carries `RC:i:<n>`
   (read coverage) and `XC:f:<frac>` (corrected coverage).

Edge cases:

- **Very low depth** (e.g., 2-read clusters): racon still runs but
  the polish is fundamentally limited by what 2 reads can correct.
  Output is technically valid but expect residual nanopore errors.
- **Racon fails outright** (no overlaps): the wrapper falls back to
  the unpolished chimeric draft and prints a warning. Reported in the
  end-of-run summary as `fallback_to_draft=N`.
- **Lowercase variants**: racon discards case, so its output is plain
  uppercase. Use the mpileup pipeline if you need the variant
  annotation.

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
