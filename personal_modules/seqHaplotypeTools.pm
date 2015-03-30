#!/usr/bin/perl
# This is a highly specialized class that deals with the data files that Jana's haplotype program generates
# First, the user must generate a HapSegMap to identify haplotype genomic coordinates
# Next, the user must generate the Individual table to associate haplotypes with individuals sequenced
package HapSegMap;
use Mouse;
use strict;
use namespace::autoclean;

has 'Map' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'getS' => 'get',
		'num_values' => 'count',
		'keyS' => 'keys',
		'putS' => 'set',});

sub loadFromFile(){
	my ($self, $file) = @_;
	open(IN, "< $file") || die "[HapSegMap] Could not open file for haplotype loading!\n";
	while(my $line = <IN>){
		chomp $line;
		$line =~ s/\r//g;
		my @segs = split(/\t/, $line);
		my $seg = HapSeg->new('chr' => "chr" . $segs[1], 'start' => $segs[2], 'end' => $segs[3]);
		$self->Map->putS($segs[0], $seg);
	}
	close IN;
}



__PACKAGE__->meta->make_immutable;

package Individual;
use Mouse;
use strict;
use namespace::autoclean;

has 'HapMap' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'getH' => 'get',
		'putH' => 'set',
		'existsH' => 'exists',
		'numH' => 'count',
		'keyH' => 'keys',});
		
has 'intersections' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', predicate => 'hasIntersections',
	handles => {
		'getI' => 'get',
		'putI' => 'set',
		'keyI' => 'keys',});


sub returnIntersections(){
	my ($self) = @_;
	my @values;
	
	foreach my $h ($self->intersections->keyI){
		foreach my $a (keys %{$self->intersections->getI($h)}){
			my $ans = $self->intersections->getI($h)->{$a};
			if(scalar(@{$ans}) > 1){
				push(@values, Intersects->new('hap' => HapNum->new('hap' => $h, 'allele' => $a),
					'animals' => $ans);
			}
		}
	}
	return @values;
}

sub calcIntersections(){
	my ($self) = @_;
	my %inters;
	foreach my $i ($self->HapMap->keyH){
		my $haps = $self->HapMap->getH($i);
		foreach my $h (@{$haps}){
			push(@{$inters{$h->hap}->{$h->allele}}, $i);
		}
	}
	$self->intersections(\%inters);
}

sub loadFromFile(){
	my ($self, $file) = @_;
	open(IN, "< $file") || die "[Individual] Could not open file for individual data loading!\n";
	while(my $line = <IN>){
		chomp $line;
		$line =~ s/\r//g;
		my @segs = split(/\t/, $line);
		my @haps = split(/,/, $segs[5]);
		my @containing;
		foreach my $h (@haps){
			my @data = split(/\./, $h);
			if(scalar(@data) != 2){next;}
			my $hap = HapNum->new('hap' => $data[0], 'allele' => $data[1]);
			push(@containing, $hap);
		}
		
		foreach my $name (@{$segs[0], $segs[1], $segs[2]}){
			if($name ne "" || length($name) > 2){
				$self->HapMap->putH($name, \@containing);
			}
		}
	}
	close IN;
}
		
__PACKAGE__->meta->make_immutable;

package Intersects;
use Mouse;
use strict;
use namespace::autoclean;

has 'hap' => (is => 'ro', isa => 'HapNum');
has 'animals' => (traits => ['Array'], is => 'ro', isa => 'ArrayRef[Any]', predicate => 'hasIntersects',
	handles => {
		'joinAn' => 'join',
		'num' => 'count',});

__PACKAGE__->meta->make_immutable;

package HapNum;
use Mouse;
use strict;
use namespace::autoclean;

has ['hap', 'allele'] => (is => 'ro', isa => 'Int');

__PACKAGE__->meta->make_immutable;

package HapSeg;
use Mouse;
use strict;
use namespace::autoclean;

has 'chr' => (is => 'ro', isa => 'Str');
has ['start', 'end'] => (is => 'ro', isa => 'Int');

__PACKAGE__->meta->make_immutable;

1;