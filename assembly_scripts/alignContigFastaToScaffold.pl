#!/usr/bin/perl
# This script is designed to take a contig fasta file and align the sequence to a scaffold fasta file
# It takes "pairs" of reads from the ends of the contig fasta and aligns them to the scaffold fasta file
# Only proper pairs (correct orientation) are considered alignments
# Output is a bed file showing contig position and orientation on the scaffold

use strict;
use Getopt::Std;

my $usage = "perl $0 -c <contig fasta file> -s <scaffold fasta file> -b <output bed file of positions and orientations>\n";

my %opts;
getopt('csb', \%opts);

unless(defined($opts{c}) && defined($opts{b}) && defined($opts{s})){
    print $usage;
    exit;
}

# Check if the FAI file is there
unless( -s "$opts{c}.fai"){
    print STDERR "Did not find fasta index file, creating one now...\n";
    system("samtools faidx $opts{c}");
}

# Part 1: Identify target sites from the FAI file
my %contigLens;
open(my $IN, "< $opts{c}.fai") || die "could not open contig fasta index file!\n";
while(my $line = <$IN>){
    chomp $line;
    my @s = split(/\t/, $line);
    $contigLens{$s[0]} = $s[1];
}
close $IN;

print STDERR "Read input fasta index file...\n";

# Part 2: generate the temp fa files
open(my $O1, "> temp.contig.1.fa");
open(my $O2, "> temp.contig.2.fa");
foreach my $contig (keys(%contigLens)){
    my $len = $contigLens{$contig};
    my $lastStart = $len - 2000;
    
    # Get beginning seq
    my $startSeq = getSamtoolsSeq($opts{c}, $contig, 1, 2000);
    print {$O1} ">$contig\n$startSeq\n";
    
    # Get end seq
    my $endSeq = getSamtoolsSeq($opts{c}, $lastStart, $len);
    my $revEnd = getRevCompSeq($endSeq);
    print {$O2} ">$contig\n$revEnd\n";
}
close $O1;
close $O2;

print STDERR "Finished generating input temp fa files...\n";

# Part 3: align the contig fa file to the scaffold fasta and check orientation/alignment
open(my $IN, "bwa mem $opts{s} temp.contig.1.fa temp.contig.2.fa |");
my %contigOrder; # {scaffoldChr} ->[] -> [start, end, contigname, orient]
my %prevContigInfo; # {contig} -> [scaffoldchr, start, first/sec, for/rev]
my @unmaps; # [contig, scaffoldchr1, pos, scaffoldchr2, pos]
while(my $line = <$IN>){
    if($line =~ /^@/){next;}
    chomp $line;
    my @segs = split(/\t/, $line);
    
    if($segs[1] & 2048){next;} # Skipping alternative alignments for now
    
    my $isFirst =($segs[1] & 64)? 1 : 0;
    my $isRev = ($segs[1] & 16)? 1 : 0;
    
    if(exists($prevContigInfo{$segs[0]})){
	# We previously saw this one. Let's populate the contig order variable
	my $prevInfo = $prevContigInfo{$segs[0]};
	
	# Check to see if we're in compliance
	my $matched = 0;
	if($segs[2] eq $prevInfo->[0] && $segs[2] ne "*"){
	    # A good start, let's see if the orientation is correct
	    if($isFirst != $prevInfo->[2] 
	       && $isRev != $prevInfo->[3] 
	       && abs($segs[3] - $prevInfo->[1]) > $contigLens{$segs[0]} - 6000){
		# also added a size check to see if the contig fits in the same size profile as what we expect
		
		# The contig is good. Add the data to the pile
		my $orient = ($isFirst && $isRev)? "-" : "+";
		my @pos = ($segs[3], $prevInfo->[1]);
		@pos = sort{$a <=> $b} @pos;
		
		push(@{$contigOrder{$segs[2]}}, [$pos[0] + 0, $pos[1] + 0, $segs[0], $orient . ""]);
		$matched = 1;
	    }
	}
	
	unless($matched){
	    # Some criteria failed. Saving this for the log dump
	    push(@unmaps, [$segs[0] . "", $segs[2] . "", $segs[3] + 0, $prevInfo->[0], $prevInfo->[1]]);
	}
    }else{
	$prevContigInfo{$segs[0]} = [$segs[2] . "", $segs[3] + 0, $isFirst, $isRev];
    }
}

close $IN;

open(my $BED, "> $opts{b}");
open(my $UNMAP, "> $opts{b}.unmapped");

foreach my $scaff (keys(%contigOrder)){
    foreach my $row (sort{$a->[0] <=> $b->[0]} @{$contigOrder{$scaff}}){
	print {$BED} "$scaff\t";
	print {$BED} join("\t", @{$row}) . "\n";
    }
}
close $BED;

# Print out unmapped contigs
foreach my $row (@unmaps){
    print {$UNMAP} join("\t", @{$row}) . "\n";
}
close $UNMAP;

exit;

sub getSamtoolsSeq{
    my ($fasta, $contig, $start, $end) = @_;
    open(my $IN, "samtools $fasta $contig:$start-$end |");
    my $h =<$IN>;
    my $seq = '';
    while(my $line = <$IN>){
	chomp $line;
	$seq .= $line;
    }
    close $IN;
    return $seq;
}

sub getRevCompSeq{
    my ($seq) = @_;
    $seq =~ tr/ACGTacgt/TGCATGCA/;
    return reverse($seq);
}
