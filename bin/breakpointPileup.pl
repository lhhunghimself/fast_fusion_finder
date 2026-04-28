#!/usr/bin/perl

# breakpointPileup.pl
#
# Companion to filterEnds.pl. filterEnds.pl prints a list of read names whose
# split alignments support a fusion breakpoint; this script takes that list and
# emits a consensus sequence over each breakpoint region.
#
# For every region in the breakpoint file we run samtools mpileup over only the
# reads in the read-name list and walk the pileup column by column. At each
# column the reference base is used by default; if a single non-reference base
# is seen in more than --threshold fraction of the parsed read bases at that
# column, that base is used instead.
#
# With --guideFile (same input format as guidesToBreakpoints.pl, namely
# "group<TAB>chr:pos<TAB>+/-"), each fusion group's two halves are matched to
# breakpoint regions by chromosome+position, reverse-complemented when the
# guide direction is '-', and concatenated in the order the entries appear in
# the guides file (first entry of the group = 5' partner). One FASTA record is
# emitted per group instead of per region.
#
# Usage:
#   breakpointPileup.pl --reads reads.txt --bam in.bam --ref ref.fa \
#                       --breakpointFile breakpoints.txt \
#                       [--guideFile guides.txt] \
#                       [--threshold 0.5] [--minDepth 1] [--verbose]
#
# The BAM must be coordinate-sorted and indexed, and the FASTA must have a
# matching .fai. The breakpoint file uses the same format as filterEnds.pl:
#   chr15:74318559:74336132;chr17:38428464:38512385
# Each side becomes one FASTA record (or one half of a fused record under
# --guideFile). Lines starting with '#' are ignored.

use warnings;
use strict;
use Getopt::Long;

my $readsFile;
my $bamFile;
my $refFile;
my $breakpointFile;
my $guideFile;
my $threshold = 0.5;
my $minDepth  = 1;
my $verbose   = 0;

GetOptions(
    "reads=s"            => \$readsFile,
    "bam=s"              => \$bamFile,
    "ref=s"              => \$refFile,
    "breakpointFile|b=s" => \$breakpointFile,
    "guideFile|G=s"      => \$guideFile,
    "threshold|t=f"      => \$threshold,
    "minDepth|d=i"       => \$minDepth,
    "verbose|v"          => \$verbose,
) or die "Error in command line arguments\n";

defined $readsFile      or die "--reads is required\n";
defined $bamFile        or die "--bam is required\n";
defined $refFile        or die "--ref is required\n";
defined $breakpointFile or die "--breakpointFile is required\n";

my @regions = readRegions($breakpointFile);
@regions or die "no usable regions in $breakpointFile\n";

my %byRegion;
foreach my $region (@regions) {
    my ($chr, $start, $end) = parseRegion($region);
    my $seq = consensusForRegion($region);
    $byRegion{$region} = { chr => $chr, start => $start, end => $end, seq => $seq };
}

if (defined $guideFile) {
    emitFusedGroups($guideFile, \%byRegion);
}
else {
    foreach my $region (@regions) {
        emitFasta($region, $byRegion{$region}{seq});
    }
}

sub consensusForRegion {
    my ($region) = @_;
    my $cmd = sprintf(
        "samtools view -b -N %s %s %s | "
      . "samtools mpileup -aa -A -B -Q 0 -f %s -r %s -",
        shellQuote($readsFile), shellQuote($bamFile), shellQuote($region),
        shellQuote($refFile),   shellQuote($region),
    );
    $verbose and print STDERR "# $cmd\n";

    open(my $fp, "-|", $cmd) or die "can't run pileup pipeline: $!";
    my $seq = "";
    while (my $line = <$fp>) {
        chomp $line;
        my ($chr, $pos, $ref, $depth, $bases) = split(/\t/, $line);
        $seq .= consensusBase($ref, $depth // 0, $bases // "");
    }
    close($fp);
    return $seq;
}

sub readRegions {
    my ($file) = @_;
    open(my $fp, '<', $file) or die "can't open $file: $!\n";
    my @out;
    while (my $line = <$fp>) {
        next if substr($line, 0, 1) eq '#';
        chomp $line;
        next if $line =~ /^\s*$/;
        foreach my $half (split /;/, $line) {
            my ($chr, $start, $end) = split /:/, $half;
            next unless defined $chr && length($chr);
            unless (defined $start && length($start)
                 && defined $end   && length($end)) {
                warn "skipping '$half': pileup needs chr:start:end\n";
                next;
            }
            ($start, $end) = ($end, $start) if $start > $end;
            push @out, "$chr:$start-$end";
        }
    }
    close($fp);
    return @out;
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
        my @parts;
        my @anno;
        foreach my $g (@{ $groups{$group} }) {
            my $info = matchRegion($g->{chr}, $g->{pos}, $byRegion);
            unless ($info) {
                warn "no breakpoint region matches guide '$group' $g->{chr}:$g->{pos}\n";
                next;
            }
            my $part = $info->{seq};
            $part = revcomp($part) if $g->{dir} eq '-';
            push @parts, $part;
            push @anno, "$g->{chr}:$info->{start}-$info->{end}$g->{dir}";
        }
        next unless @parts;
        my $fusion = join('', @parts);
        emitFasta("$group " . join(';', @anno), $fusion);
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

# Walks an mpileup bases string, returning the consensus base for the column.
# mpileup encoding: '.'/',' = match to ref; ACGTN (any case) = mismatch;
# '^X' starts a read (X is mapq, skip it); '$' ends a read; '+N<seq>' /
# '-N<seq>' are insertions/deletions to the reference following this column;
# '*'/'#' is a deletion that covers this column.
sub consensusBase {
    my ($ref, $depth, $bases) = @_;
    $ref = uc($ref // 'N');
    return $ref if !$depth || $depth < $minDepth;
    return $ref if !length($bases) || $bases eq '*';

    my %count;
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
        if ($ch eq '.' || $ch eq ',') { $count{$ref}++; $i++; next; }
        $count{uc($ch)}++;
        $i++;
    }

    my $total = 0;
    foreach my $b (qw(A C G T N)) { $total += $count{$b} if $count{$b}; }
    return $ref if $total == 0;

    foreach my $b (qw(A C G T)) {
        next if $b eq $ref;
        return $b if defined $count{$b} && $count{$b} / $total > $threshold;
    }
    return $ref;
}

sub shellQuote {
    my ($s) = @_;
    $s =~ s/'/'\\''/g;
    return "'$s'";
}
