#!/usr/bin/perl
# This script estimates regions of excessive hard or soft clipping in an area of the genome
# Input is a bed file

use strict;
use Getopt::Std;

my $usage = "perl $0 -i <input bed file with coordinates> -b <input indexed bam for searching> -o <output>\n";
my %opts;
getopt('ibo', \%opts);

unless(defined($opts{'i'}) && defined($opts{'b'})){
	print $usage;
	exit;
}

# Read in bed file and create samtools search for each entry
open(IN, "< $opts{i}") || die "Could not open bed file!\n$usage";
open(OUT, "> $opts{o}");
print OUT "snpname\tclipstart\tclipend\tregstart\tregend\n";
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my $str = "$segs[0]:$segs[1]-$segs[2]";
	my @positions; # []->[s1, s2, e1, e2]
	open(SAM, "samtools view $opts{b} \'$str\' |") || die "Could not index bam file!\n$usage";
	while(my $samstr = <SAM>){
		chomp $samstr;
		my @ssegs = split(/\t/, $samstr);
		my $cigar = Cigar->new();
		$cigar->populateCigar($ssegs[5]);
		
		my($s1, $s2, $e1, $e2) = determineLength($ssegs[3], $cigar);
		push(@positions, [$s1, $s2, $e1, $e2]);
	}
	close SAM;
	
	my @coords = concordance(\@positions, $segs[1], $segs[2], $segs[3]);
	foreach  my $c (@coords){
		print OUT $c;
	}
}
close IN;
close OUT;

exit;

sub concordance{
	my ($array, $start, $end, $name) = @_;
	my $neg = 0;
	my ($s1, $s2, $e1, $e2);
	my $len = $end - $start;
	$s1 = 9999999999;
	$s2 = 0;
	$e1 = 9999999999;
	$e2 = 0;
	
	my @output;
	foreach my $row (@{$array}){
		
		if($row->[0] == 0 && $row->[1] == 0 && $neg == 0 && $row->[2] < ($end - ($len * 0.5)) && $row->[3] > ($start - ($len * 0.5))){
			# There were no hard/softclipped bases and the region overlapped with the predicted start and end position by 50%
			my $str = outStr($row, $name);
			push(@output, $str);
		}
		
		if($row->[0] < $s1){
			$s1 = $row->[0];
		}
		if($row->[1] > $s2){
			$s2 = $row->[1];
		}
		if($row->[2] < $e1){
			$e1 = $row->[2];
		}
		if($row->[3] > $e2){
			$e2 = $row->[3];
		}
	}
	my $str = outStr([$s1, $s2, $e1, $e2], $name);
	push(@output, $str);
	return @output;
}

sub outStr{
	my ($array, $name) = @_;
	my $str = "$name\t" . $array->[0] . "\t" . $array->[1] . "\t" . $array->[2] . "\t" . $array->[3] . "\n";
	return $str;
}

sub determineLength{
	my ($pos, $cigar) = @_;
	my $e1 = $pos;
	my ($s2, $s1, $e2);
	$s2 = 0;
	$e1 = 0;
	$e2 = 0;
	
	# I'm going to make $s1 and $s2 the hard/softclip bed locs
	my $current = $pos;
	for(my $x = 0; $x < $cigar->num(); $x++){
		my $t = $cigar->get_tag($x);
		my $v = $cigar->get_value($x);
		if($t eq 'S' || $t eq 'H'){
			$s1 = $current;
			$s2 = $current + $v;
		}else{
			$e2 = $current;
		}
		$current += $v;
		
	}
	return ($s1, $s2, $e1, $e2);
}

BEGIN{
package Cigar;
use Mouse;
use namespace::autoclean;
use strict;

has 'num' => (is => 'rw', isa => 'Int', default => 0);

has 'values' => (
	traits => ['Array'],
	is => 'rw',
	isa => 'ArrayRef[Int]',
	default => sub{[]},
	handles => {
		'add_value' => 'push',
		'get_value' => 'get',
	}
);

has 'tags' => (
	traits => ['Array'],
	is => 'rw',
	isa => 'ArrayRef[Str]',
	default => sub{[]},
	handles => {
		'add_tag' => 'push',
		'get_tag' => 'get',
	}
);

sub populateCigar{
	my ($self, $cstr) = @_;
	while($cstr =~ /(\d+)(\D{1})/g){
		$self->add_value($1);
		$self->add_tag($2);
		$self->num($self->num + 1);
	}
}

__PACKAGE__->meta->make_immutable;

1;
}