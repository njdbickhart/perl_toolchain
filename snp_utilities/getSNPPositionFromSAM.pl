#!/usr/bin/perl
# This script is designed to pull likely variant site information from a sam file. 
# Assumes that the SNP probe sequence was used, and that the 3' flanking sequence was reverse complemented
# Also assumes that there is a ".f" or ".r" suffix on the read name

use strict;

my $usage = "perl $0 <sam file> <output tab delimited file>\n";

chomp(@ARGV);
unless(scalar(@ARGV) == 2){
	print $usage;
	exit;
}

my %snpcontainer;

open(my $IN, "< $ARGV[0]") || die "Could not open sam file!\n$usage";
open(my $OUT, "> $ARGV[1]");
while(my $line = <$IN>){
	chomp $line;
	if($line =~ /^@/){next;}
	
	my @segs = split(/\t/, $line);
	# Count the number of insertions and deletions in the CIGAR string
	my $ins = 0;
	my $del = 0;
	my $hasclipping = 0;
	while($segs[5] =~ /(\d+)([MDISH])/){
		if($2 eq 'D'){
			$del += $1;
		}elsif($2 eq 'I'){
			$ins += $1;
		}elsif($2 eq 'S' || $2 eq 'H'){
			$hasclipping = 1;
		}
	}
	
	my $pos = $segs[3] + length($segs[9]) - $ins + $del;
	my @namesegs = split(/\./, $segs[0]);
	if(!exists($snpcontainer{$namesegs[0]})){
		$snpcontainer{$namesegs[0]} = VariantSite->new('chr' => $segs[2], 'pos' => $pos);
	}else{
		$snpcontainer{$namesegs[0]}->checkConsistency($segs[2], $pos);
	}
	
	if($namesegs[1] eq 'f'){
		if($del > 0 || $ins > 0 || $hasclipping){
			$snpcontainer{$namesegs[0]}->cigarfor(1);
		}
	}else{
		if($del > 0 || $ins > 0 || $hasclipping){
			$snpcontainer{$namesegs[0]}->cigarrev(1);
		}
	}
}
close $IN;

print $OUT "ProbeName\tChr\tPos\tConsistent\tCigarProblem5prime\tCigarProblem3prime\tAltChr\tAltPos\n";
foreach my $name (sort{$a cmp $b} keys(%snpcontainer)){
	print $OUT "$name\t";
	print $OUT $snpcontainer{$name}->chr . "\t";
	print $OUT $snpcontainer{$name}->pos . "\t";
	print $OUT $snpcontainer{$name}->consistent . "\t";
	print $OUT $snpcontainer{$name}->cigarfor . "\t";
	print $OUT $snpcontainer{$name}->cigarrev . "\t";
	
	if($snpcontainer{$name}->altpos != 0){
		print $OUT $snpcontainer{$name}->altchr . "\t";
		print $OUT $snpcontainer{$name}->altpos . "\t";
	}
	
	print $OUT "\n";
}

print STDERR "Data is in $ARGV[1]\n";

exit;

BEGIN{
package VariantSite;
use Mouse;
use namespace::autoclean;

has ['chr', 'altchr'] => (is => 'rw', isa => 'Str', default => '');
has ['pos', 'altpos'] => (is => 'rw', isa => 'Num', default => 0);
has ['consistent', 'cigarfor', 'cigarrev'] => (is => 'rw', isa => 'Bool', default => 0);

sub checkConsistency{
	my ($self, $chr, $pos) = @_;
	if($self->chr eq $chr && $self->pos == $pos){
		$self->consistent(1);
	}else{
		$self->consistent(0);
		$self->altchr($chr);
		$self->altpos($pos);
	}
}

__PACKAGE__->meta->make_immutable;

1; 
}