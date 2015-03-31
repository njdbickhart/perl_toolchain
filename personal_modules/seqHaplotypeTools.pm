#!/usr/bin/perl
# This is a highly specialized class that deals with the data files that Jana's haplotype program generates
# First, the user must generate a HapSegMap to identify haplotype genomic coordinates
# Next, the user must generate the Individual table to associate haplotypes with individuals sequenced
package VCFHetFilter;
use Mouse;
use strict;
use namespace::autoclean;

has ['vcffile', 'outbase'] => (is => 'ro', isa => 'Str', required => 1);
has ['retained', 'possiblehom', 'possiblehet', 'filtered'] => (is => 'rw', isa => 'Str');
has 'vcfhead' => (traits => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub{[]}, handles => {
	'add' => 'push',});

sub startMultHetFilter(){
	my ($self) = @_;
	
	# Generate output file names
	$self->retained($self->outbase . ".fixed.vcf");
	$self->possiblehom($self->outbase . ".homozygous.vcf");
	$self->possiblehet($self->outbase . ".heterozygous.vcf");
	$self->filtered($self->outbase . ".filtered.vcf");
	
	open(IN, "< " . $self->vcffile) || die "[VCFHETFILT] Could not open input vcf file!\n";
	open(FIX, "> " . $self->retained);
	open(HOM, "> " . $self->possiblehom);
	open(HET, "> " . $self->possiblehet);
	open(FIL, "> " . $self->filtered);
	
	my $header = 0;
	my $filtered = 0;
	my $fixed = 0;
	while(my $line = <IN>){
		chomp $line;
		if($line =~ /^#/){
			$self->add($line);
			next;
		}
		
		if($header == 0){
			$header = 1;
			my @head = @{$self->vcfhead};
			print FIX join("\n", @head) . "\n";
			print HOM join("\n", @head) . "\n";
			print HET join("\n", @head) . "\n";
			print FIL join("\n", @head) . "\n";
		}
		
		my @segs = split(/\t/, $line);
		my $num = scalar(@segs);
		
		my $hets = 0;
		my $homs = 0;
		foreach my $gt (@$segs[9 .. $num]){
			my @gsegs = split(/:/, $gt);
			my @csegs = split(/\//, $gsegs[0]);
			if($csegs[0] eq '0' || $csegs[1] eq '0'){
				$hets++;
			}elsif($csegs[0] eq $cesegs[1]){
				$homs++;
			}else{
				$hets++;
			}
		}
		
		if($homs == $num){
			print FIX $line . "\n";
			$fixed++;
		}elsif($hets > 1){
			print HET $line . "\n";
		}elsif($homs > 1){
			print HOM $line . "\n";
		}else{
			print FIL $line . "\n";
			$filtered++;
		}
	}
	close IN;
	close FIX;
	close HET;
	close HOM;
	close FIL;
	
	print "[VCFHETFILT] For file: " . $self->vcffile . " found $fixed fixed homozygotes and $filtered unconfident heterozygote calls\n";
}


__PACKAGE__->meta->make_immutable;

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
	my ($self, $available) = @_;
	my @values;
	my @singles;
	
	foreach my $h ($self->intersections->keyI){
		foreach my $a (keys %{$self->intersections->getI($h)}){
			my $ans = $self->intersections->getI($h)->{$a};
			my $num = scalar(@{$ans});
			if($num > 1){
				my @found;
				foreach my $b (@{$available}){
					foreach my $c (@{$ans}){
						if($b eq $c){
							push(@found, $c);
						}
					}
				}
				if(scalar(@found) > 1){
					push(@values, Intersects->new('hap' => HapNum->new('hap' => $h, 'allele' => $a),
						'animals' => \@found);
				}elsif(scalar(@found) == 1){
					push(@singles, Intersects->new('hap' => HapNum->new('hap' => $h, 'allele' => $a),
						'animals' => \@found);
				}
			}
		}
	}
	return \@values, \@singles;
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

sub printUCSC(){
	my ($self) = @_;
	my $str = $self->chr() . ":";
	$str .= $self->start() . "-" . $self->end();
	return $str;
}

__PACKAGE__->meta->make_immutable;

1;