#!/usr/bin/perl
# This is a short script designed to output key alignments of sam files without taking up the entire real estate of the screen
# Also, processes alignments to bed file for fast 

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 (-s <input sam> || -b <input bam>) [-d output to bed file]\n";

getopt('sb', \%opts);

unless(defined($opts{'s'}) || defined($opts{'b'})){
	print $usage;
	exit;
}

my $IN;
if(defined($opts{'s'})){
	open($IN, "< $opts{s}") || die "Could not open input sam file: $opts{s}!\n$usage";
}else{
	open($IN, "samtools view $opts{b} |");
}

while(my $line = <$IN>){
	chomp $line;
	if($line =~ /^@/){next;}
	
	my @segs = split(/\t/, $line);
	my $alnEnd = getAlignLen($segs[5]) + $segs[3];
	my $seqLen = length($segs[9]);
	
	print "$segs[0]\t$segs[1]\t$segs[2]\t$segs[3]\t$alnEnd\t$seqLen\t$segs[4]\n";
}

sub getAlignLen{
	my ($cigar) = @_;
	my $len = 0;
	while($cigar =~ /(\d{1,6})(\D{1})/g){
		my $bp = $1;
		my $code = $2;
		if($code =~ /[M=XDN]/){
			$len += $bp;
		}
	}
	return $len;
}
