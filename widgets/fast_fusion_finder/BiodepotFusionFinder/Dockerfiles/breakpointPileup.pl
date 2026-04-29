#!/usr/bin/perl

# breakpointPileup.pl
#
# Companion to filterEnds.pl. filterEnds.pl prints a list of read names whose
# split alignments support a fusion breakpoint; this script takes that list
# and emits the two halves of the fusion as separate FASTA records, each
# windowed around the breakpoint.
#
# Per side, the breakpoint position comes from the breakpoint file itself:
# the second coord of side 1 (5' partner) and the first coord of side 2
# (3' partner), matching the convention findBreakPoints uses to write them.
# Each side is then piled up over a window of --window bases on the
# transcript-internal flank of the breakpoint:
#   5' partner: the W bases ending at the breakpoint in transcript order
#   3' partner: the W bases starting at the breakpoint in transcript order
# Output is one FASTA record per side, in line order (5' first, 3' second).
#
# At each pileup column the reference base is used by default; if a single
# non-reference base is seen in more than --threshold fraction of the parsed
# read bases at that column, that base is used instead and emitted in lower
# case. Reverse-complementing for - strand sides preserves case.
#
# Strand per side comes from the start:end coord order in the breakpoint
# file: start<end is +, start>end is - (the upstream majority-rule generator
# only emits correctly stranded matches). --guideFile overrides that.
#
# --pileupOut writes a per-position TSV (chr, pos, ref, consensus, depth,
# A, C, G, T) for every pileup column visited. Genomic order, one row per
# column. The "consensus" column matches what goes into the FASTA, so
# uppercase = ref, lowercase = called variant.
#
# Usage:
#   breakpointPileup.pl --reads reads.txt --bam in.bam --ref ref.fa \
#                       --breakpointFile breakpoints.txt \
#                       [--window 100] [--pileupOut variants.tsv] \
#                       [--guideFile guides.txt] \
#                       [--threshold 0.5] [--minDepth 1] [--verbose]
#
# The BAM must be coordinate-sorted and indexed, and the FASTA must have a
# matching .fai. The breakpoint file uses the same format as filterEnds.pl:
#   chr15:74318559:74336132;chr17:38428464:38512385
# Lines starting with '#' are ignored.

use warnings;
use strict;
use Getopt::Long;
use File::Temp qw(tempdir);

my $readsFile;
my $bamFile;
my $refFile;
my $breakpointFile;
my $guideFile;
my $pileupOut;
my $threshold = 0.5;
my $minDepth  = 1;
my $window    = 100;
my $verbose   = 0;

GetOptions(
    "reads=s"            => \$readsFile,
    "bam=s"              => \$bamFile,
    "ref=s"              => \$refFile,
    "breakpointFile|b=s" => \$breakpointFile,
    "guideFile|G=s"      => \$guideFile,
    "pileupOut=s"        => \$pileupOut,
    "threshold|t=f"      => \$threshold,
    "minDepth|d=i"       => \$minDepth,
    "window|w=i"         => \$window,
    "verbose|v"          => \$verbose,
) or die "Error in command line arguments\n";

defined $readsFile      or die "--reads is required\n";
defined $bamFile        or die "--bam is required\n";
defined $refFile        or die "--ref is required\n";
defined $breakpointFile or die "--breakpointFile is required\n";
$window > 0             or die "--window must be positive\n";

my $pileupFh;
if (defined $pileupOut) {
    open($pileupFh, '>', $pileupOut) or die "can't write $pileupOut: $!\n";
    print $pileupFh "#chr\tpos\tref\tconsensus\tdepth\tA\tC\tG\tT\n";
}

my $tmpDir     = tempdir(CLEANUP => 1);
my $tmpCounter = 0;

my @pairs = readBreakpointPairs($breakpointFile);
@pairs or die "no usable breakpoint lines in $breakpointFile\n";

my %byRegion;
foreach my $pair (@pairs) {
    foreach my $side (@$pair) {
        $side->{region} = windowForSide($side);
        my $region = $side->{region};
        next if exists $byRegion{$region};
        my ($chr, $start, $end) = parseRegion($region);
        my $seq = consensusForRegion($region, $pileupFh);
        $byRegion{$region} = { chr => $chr, start => $start, end => $end, seq => $seq };
    }
}

if (defined $guideFile) {
    emitFusedGroups($guideFile, \%byRegion);
}
else {
    foreach my $pair (@pairs) {
        emitFusedPair($pair, \%byRegion);
    }
}

close($pileupFh) if $pileupFh;

sub emitFusedPair {
    my ($pair, $byRegion) = @_;
    foreach my $side (@$pair) {
        my $info = $byRegion->{ $side->{region} };
        my $seq  = $info->{seq};
        $seq = revcomp($seq) if $side->{dir} eq '-';
        emitFasta("$side->{region}$side->{dir}", $seq);
    }
}

sub consensusForRegion {
    my ($region, $pileupFh) = @_;

    # mpileup -r requires an indexed BAM, so we materialize a per-region
    # filtered BAM under $tmpDir rather than piping samtools view in.
    $tmpCounter++;
    my $regionBam = "$tmpDir/region$tmpCounter.bam";
    my $filterCmd = sprintf(
        "samtools view -b -N %s -o %s %s %s",
        shellQuote($readsFile), shellQuote($regionBam),
        shellQuote($bamFile),   shellQuote($region),
    );
    $verbose and print STDERR "# $filterCmd\n";
    runOrDie($filterCmd);
    runOrDie("samtools index " . shellQuote($regionBam));

    my $pileCmd = sprintf(
        "samtools mpileup -aa -A -B -Q 0 -f %s -r %s %s",
        shellQuote($refFile), shellQuote($region), shellQuote($regionBam),
    );
    $verbose and print STDERR "# $pileCmd\n";

    open(my $fp, "-|", $pileCmd) or die "can't run mpileup: $!";
    my $seq = "";
    while (my $line = <$fp>) {
        chomp $line;
        my ($chr, $pos, $ref, $depth, $bases) = split(/\t/, $line);
        my ($call, $counts) = consensusBase($ref, $depth // 0, $bases // "");
        $seq .= $call;
        if ($pileupFh) {
            print $pileupFh join("\t",
                $chr, $pos, uc($ref // 'N'), $call, $depth // 0,
                $counts->{A}, $counts->{C}, $counts->{G}, $counts->{T},
            ), "\n";
        }
    }
    close($fp);
    return $seq;
}

sub readBreakpointPairs {
    my ($file) = @_;
    open(my $fp, '<', $file) or die "can't open $file: $!\n";
    my @pairs;
    while (my $line = <$fp>) {
        next if substr($line, 0, 1) eq '#';
        chomp $line;
        next if $line =~ /^\s*$/;
        my @halves = split /;/, $line;
        my @sides;
        for (my $i = 0; $i < @halves; $i++) {
            my $half = $halves[$i];
            my ($chr, $a, $b) = split /:/, $half;
            next unless defined $chr && length($chr);
            unless (defined $a && length($a)
                 && defined $b && length($b)) {
                warn "skipping '$half' on '$line': pileup needs chr:start:end\n";
                next;
            }
            my $dir = ($a > $b) ? '-' : '+';
            my $bp  = ($i == 0) ? $b : $a;
            push @sides, { chr => $chr, bp => 0 + $bp, dir => $dir, idx => $i };
        }
        push @pairs, \@sides if @sides;
    }
    close($fp);
    return @pairs;
}

# Genomic window for one side of a fusion pair.
# 5' partner (idx 0) wants the W bases ending at bp in transcript order.
# 3' partner (idx 1) wants the W bases starting at bp in transcript order.
# Combined with strand: + transcript flows up genomically, - transcript
# flows down. The two true cases below pick "going up" vs "going down".
sub windowForSide {
    my ($side) = @_;
    my $bp  = $side->{bp};
    my $dir = $side->{dir};
    my $idx = $side->{idx};
    my ($lo, $hi);
    if (($idx == 0 && $dir eq '+') || ($idx == 1 && $dir eq '-')) {
        $lo = $bp - $window + 1;
        $hi = $bp;
    }
    else {
        $lo = $bp;
        $hi = $bp + $window - 1;
    }
    $lo = 1 if $lo < 1;
    return "$side->{chr}:$lo-$hi";
}

sub parseRegion {
    my ($r) = @_;
    $r =~ /^([^:]+):(\d+)-(\d+)$/ or die "bad region '$r'\n";
    return ($1, $2, $3);
}

sub emitFusedGroups {
    my ($file, $byRegion) = @_;
    open(my $fp, '<', $file) or die "can't open $file: $!\n";
    my %groups;
    my @order;
    while (my $line = <$fp>) {
        next if substr($line, 0, 1) eq '#';
        chomp $line;
        next if $line =~ /^\s*$/;
        my ($group, $loc, $dir) = split /\t/, $line;
        unless (defined $group && defined $loc && defined $dir) {
            warn "skipping malformed guide line: $line\n";
            next;
        }
        $dir =~ s/\s+//g;
        unless ($dir eq '+' || $dir eq '-') {
            warn "skipping guide '$group' $loc: direction must be + or -\n";
            next;
        }
        my ($chr, $pos) = split /:/, $loc;
        unless (defined $chr && defined $pos && $pos =~ /^\d+$/) {
            warn "skipping guide '$group' $loc: expected chr:pos\n";
            next;
        }
        push @order, $group unless exists $groups{$group};
        push @{ $groups{$group} }, { chr => $chr, pos => $pos, dir => $dir };
    }
    close($fp);

    foreach my $group (@order) {
        foreach my $g (@{ $groups{$group} }) {
            my $info = matchRegion($g->{chr}, $g->{pos}, $byRegion);
            unless ($info) {
                warn "no breakpoint region matches guide '$group' $g->{chr}:$g->{pos}\n";
                next;
            }
            my $seq = $info->{seq};
            $seq = revcomp($seq) if $g->{dir} eq '-';
            emitFasta("$group $g->{chr}:$info->{start}-$info->{end}$g->{dir}", $seq);
        }
    }
}

sub matchRegion {
    my ($chr, $pos, $byRegion) = @_;
    my $chrFallback;
    foreach my $r (keys %$byRegion) {
        my $info = $byRegion->{$r};
        next unless $info->{chr} eq $chr;
        $chrFallback ||= $info;
        return $info if $pos >= $info->{start} && $pos <= $info->{end};
    }
    return $chrFallback;
}

sub emitFasta {
    my ($header, $seq) = @_;
    print ">$header length=", length($seq), "\n";
    for (my $i = 0; $i < length($seq); $i += 80) {
        print substr($seq, $i, 80), "\n";
    }
}

sub revcomp {
    my ($s) = @_;
    $s =~ tr/ACGTNacgtn/TGCANtgcan/;
    return scalar reverse $s;
}

# Returns ($base, \%counts). $base is uppercase when it equals the reference
# (or no variant call could be made), lowercase when a non-reference base won
# >threshold of the parsed reads. %counts is A/C/G/T tallies for the column.
#
# mpileup encoding: '.'/',' = match to ref; ACGTN (any case) = mismatch;
# '^X' starts a read (X is mapq, skip it); '$' ends a read; '+N<seq>' /
# '-N<seq>' are insertions/deletions to the reference following this column;
# '*'/'#' is a deletion that covers this column.
sub consensusBase {
    my ($ref, $depth, $bases) = @_;
    $ref = uc($ref // 'N');

    my %count = (A => 0, C => 0, G => 0, T => 0, N => 0);
    parseBases($bases, $ref, \%count) if length($bases);

    if (!$depth || $depth < $minDepth || !length($bases) || $bases eq '*') {
        return ($ref, \%count);
    }

    my $total = $count{A} + $count{C} + $count{G} + $count{T} + $count{N};
    return ($ref, \%count) if $total == 0;

    foreach my $b (qw(A C G T)) {
        next if $b eq $ref;
        return (lc($b), \%count) if $count{$b} / $total > $threshold;
    }
    return ($ref, \%count);
}

sub parseBases {
    my ($bases, $ref, $count) = @_;
    return if $bases eq '*';
    my @c = split(//, $bases);
    my $i = 0;
    while ($i <= $#c) {
        my $ch = $c[$i];
        if ($ch eq '^') { $i += 2; next; }
        if ($ch eq '$') { $i += 1; next; }
        if ($ch eq '+' || $ch eq '-') {
            $i++;
            my $n = "";
            while ($i <= $#c && $c[$i] =~ /\d/) { $n .= $c[$i]; $i++; }
            $i += int($n);
            next;
        }
        if ($ch eq '*' || $ch eq '#') { $i++; next; }
        if ($ch eq '.' || $ch eq ',') { $count->{$ref}++; $i++; next; }
        $count->{uc($ch)}++;
        $i++;
    }
}

sub shellQuote {
    my ($s) = @_;
    $s =~ s/'/'\\''/g;
    return "'$s'";
}

sub runOrDie {
    my ($cmd) = @_;
    system($cmd) == 0 or die "command failed (exit $?): $cmd\n";
}
