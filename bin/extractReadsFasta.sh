#!/usr/bin/env bash
# extractReadsFasta.sh — pull reads by name from a FASTQ and emit them as FASTA.
#
# Usage:
#   extractReadsFasta.sh NAMES_FILE FASTQ [FASTQ ...]  > out.fa
#
# - NAMES_FILE is a text file of read names. Only the first whitespace-
#   separated field of each line is consulted, so a TSV (e.g. the
#   .readnamescoords used elsewhere in this repo: "name<TAB>chr:s:e;...")
#   works without preprocessing. Lines starting with '#' are ignored.
# - FASTQ may be plain or gzipped (auto-detected by the '.gz' suffix).
# - Multiple FASTQs can be passed; they're concatenated in the order
#   given. Output records appear in FASTQ order (not names-file order).
# - Names not found in any FASTQ are listed on stderr at the end.
# - Read-name comparison strips the leading '@' and any trailing
#   whitespace+comment from the FASTQ header line, so headers like
#   "@abc123 length=12345 ..." match the bare name "abc123".
#
# Portability: works with gawk or mawk, requires only POSIX awk
# features. No Perl, no Python, no samtools. zcat is needed only for
# gzipped inputs.

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "usage: $0 NAMES_FILE FASTQ [FASTQ ...] > out.fa" >&2
    exit 1
fi

names=$1
shift

[[ -f "$names" ]] || { echo "names file not found: $names" >&2; exit 1; }

# Concatenate every FASTQ (decompressing .gz on the fly) into one stream
# and let awk pluck the wanted records out of it.
{
    for f in "$@"; do
        [[ -f "$f" ]] || { echo "fastq not found: $f" >&2; exit 1; }
        case "$f" in
            *.gz)  zcat -- "$f";;
            *)     cat  -- "$f";;
        esac
    done
} | awk -v namesFile="$names" '
BEGIN {
    nWanted = 0
    while ((getline line < namesFile) > 0) {
        sub(/^[[:space:]]+/, "", line)
        if (line == "" || substr(line, 1, 1) == "#") continue
        n = split(line, parts, /[[:space:]]+/)
        if (parts[1] != "" && !(parts[1] in wanted)) {
            wanted[parts[1]] = 1
            nWanted++
        }
    }
    close(namesFile)
    if (nWanted == 0) {
        print "no read names parsed from " namesFile > "/dev/stderr"
        exit 1
    }
}

# Walk the concatenated FASTQ four lines at a time. Sequence lines never
# start with "@" by spec, but quality lines can; reading exactly four
# lines from a known "@" header keeps us in sync.
{
    if (substr($0, 1, 1) != "@") next
    name = substr($0, 2)
    sub(/[[:space:]].*$/, "", name)
    if ((getline seq)  <= 0) exit
    if ((getline plus) <= 0) exit
    if ((getline qual) <= 0) exit
    if (name in wanted && !(name in seen)) {
        print ">" name
        print seq
        seen[name] = 1
        nFound++
    }
}

END {
    if (nFound + 0 < nWanted) {
        for (n in wanted) {
            if (!(n in seen)) print "missing: " n > "/dev/stderr"
        }
        printf("extracted %d/%d reads\n", nFound + 0, nWanted) > "/dev/stderr"
    }
}
'
