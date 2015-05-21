#!/usr/bin/perl
# This is a modification of my hard clipping check script that is designed to resolve ambiguous liftover mappings
# It is highly specific to my work in Goat right now, but I can expand it in the future

use strict;
use Getopt::Std;

my $usage = "perl $0 -i <input bam with tags in readnames> -o <tabular output>\n";
my %opts;
getopt('io', \%opts);

unless(defined($opts{'i'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

# Read in bed file and create samtools search for each entry
open(IN, "samtools view $opts{i} |") || die "Could not open bam file!\n$usage";
open(OUT, "> $opts{o}");
print OUT "region\tgffID\talign1chr\tflag\tcigar\tstart\tend\talign2chr\tflag\tcigar\tstart\tend\n";
my %store; # {readname} -> [align1, align2]
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	push(@{$store{$segs[0]}}, \@segs);
}
close IN;

foreach my $read (sort { $a cmp $b } keys(%store)){
	my @rows;
	my @temp = sort{$a->[1] <=> $b->[1]} @{$store{$read}};
	
	my @rtags = split(/_/, $read);
	push(@rows, @rtags);
	
	foreach my $aln (@temp){
		my $cigar = Cigar->new();
		$cigar->populateCigar($aln->[5]);
		my $len = $cigar->getLenLargestClip();
		push(@rows, $aln->[2]);
		my $elen = ($aln->[1] & 0x100)? (length($aln->[9]) - $len) + $aln->[3] : $aln->[9] + $aln->[3];
		push(@rows, ($aln->[1], $aln->[5], $aln->[3], $elen));
	}		
	
	print OUT join("\t", @rows); 
	print OUT "\n";
}

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
	$s1 = 0;
	$e2 = 0;
	
	# I'm going to make $s1 and $s2 the hard/softclip bed locs
	my $current = $pos;
	for(my $x = 0; $x < $cigar->num(); $x++){
		my $t = $cigar->get_tag($x);
		my $v = $cigar->get_value($x);
		if($t eq 'S' || $t eq 'H'){
			$s1 = $current;
			$s2 = $current + $v;
			if($x == 0){
				# The first part of the read was bookended by a hard clip
				$e1 = $current + $v;
			}
		}else{
			$e2 = $current + $v;
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
		'Vcount' => 'count',
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

# I am assuming that the largest clip area is at the end of the seq string
sub getLenLargestClip{
	my ($self) = @_;
	my $largest = 0;
	for(my $x = 0; $x < $self->Vcount; $x++){
		if($self->get_tag($x) eq 'H' || $self->get_tag($x) eq 'X'){
			if($self->get_value($x) > $largest){
				$largest = $self->get_value($x);
			}
		}
	}
	return $largest;
}

__PACKAGE__->meta->make_immutable;

1;
}