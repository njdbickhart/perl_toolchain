#!/usr/bin/perl
# This script processes a UCSC bed file and generates fasta sequence using the "start and end" 
# coordinates present in the bed columns

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -b <UCSC bed> -f <Reference fasta> -e <[optional] only search for these> -o <output>\n";
getopt('bfeo', \%opts);

unless(defined($opts{'b'}) && defined($opts{'f'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my %keepList;
if(defined($opts{'e'})){
	open(my $IN, "< $opts{e}") || die "Could not open input Search list!\n";
	while(my $line = <$IN>){
		chomp $line;
		$keepList{$line} = 1;
	}
	close $IN;
}

# Check if fasta file has been fasta indexed by samtools
unless( -s "$opts{f}.fai"){
	print STDERR "Could not find fasta index file for reference: $opts{f}\n";
	print STDERR "Generating now...\n";
	system("samtools faidx $opts{f}");
}

# Open bed file and process needed entries
open(my $OUT, "> $opts{o}");
open(my $IN, "< $opts{b}") || die "Could not open UCSC bed file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	# Take care of the "version" dot notation that UCSC loves to place on EnsTran IDs
	$segs[1] =~ s/\.\d{1,2}$//; 
	
	if(defined($opts{'e'})){
		if(exists($keepList{$segs[1]})){
			grepSamtoolsFasta($segs[1], $segs[2], $segs[3], $segs[9], $segs[10], $opts{'f'}, $OUT);
		}
	}else{
		grepSamtoolsFasta($segs[1], $segs[2], $segs[3], $segs[9], $segs[10], $opts{'f'}, $OUT);
	}
}
close $IN;
close $OUT;

print STDERR "Output is in: $opts{o}\n";

exit;

sub grepSamtoolsFasta{
	my ($name, $chr, $orient, $scoords, $ecoords, $ref, $OUT) = @_;
	
	# Get rid of trailing comma
	$scoords =~ s/\,$//;
	$ecoords =~ s/\,$//;

	my @starts = split(/,/, $scoords);
	my @ends = split(/,/, $ecoords);
	
	my $exoncount = ($orient eq '-')? scalar(@starts) : 1;
	for(my $x = 0; $x < scalar(@starts); $x++){
		my $ucsc = "$chr:$starts[$x]-$ends[$x]";
		open(my $SAM, "samtools faidx $ref $ucsc | ");
		while(my $line = <$SAM>){
			if($line =~ /^>/){
				# Steal the fasta carrot and replace it
				print $OUT ">$name\_$exoncount\n";
			}else{
				print $OUT $line;
			}
		}
		close $SAM;
		$exoncount += ($orient eq '-')? -1 : 1;
	}
			
			
}