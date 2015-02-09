#!/usr/bin/perl
# A collection of classes used to manipulate bed files in Perl

package BedCoord;
use strict;
use Mouse;
use namespace::autoclean;

has 'chr' => ( is => 'rw', isa => 'Str',);
has 'start' => ( is => 'rw', isa => 'Int',);
has 'end' => ( is => 'rw', isa => 'Int',);
has 'other' => ( is => 'rw', isa => 'ArrayRef[Any]',);

__PACKAGE__->meta->make_immutable;

package BedContainer;
use strict;
use Mouse;
use kentBinTools;
use namespace::autoclean;

has 'bed' => (is => 'rw', isa => 'HashRef[HashRef[ArrayRef[genomic_coord]]]', predicate => 'has_bed',); # {chr}->{bin}->[0]->genomic_coord
has 'lines' => (is => 'rw', isa => 'Int',);

sub loadFile {
	my ($self, $file) = @_;
	
	my $binner = kentBinTools->new();
	open(IN, "< $file") || die "[perlBed] Could not open $file!\n";
	while (my $line = <IN>){
		chomp  $line;
		$line =~ s/\r//g;
		if($line eq '' || $line =~ /^\s+$/){next;}
		my @segs = split(/\t/, $line);
		if(!(scalar(@segs) >= 4)){
			print STDERR "[perlBed] Error! Line has less than four categories! $line\n";
			exit;
		}
		
		my $segNum = scalar(@segs);
		my $entry = BedCoord->new('chr' => $segs[0], 'start' => $segs[1], 'end' => $segs[2], 'other' => [@$segs[3 .. $segNum]]);
		my $bin = $binner->getbin($segs[1], $segs[2]);
		if(!($self->has_bed())){
			my %h;
			push(@{$h{$segs[0]}->{$bin}}, $entry);
			$self->bed(\%h);
		}else{
			my %h = $self->bed();
			push(@{$h{$segs[0]}->{$bin}}, $entry);
			$self->bed(\%h);
		}
		$self->lines($self->lines() + 1);
	}
	close IN;
	print STDERR "[perlBed] Finished loading bed file: $file\n";
}

# Bool: does the input information intersect with the stored data
sub intersects {
	my ($self, $chr, $start, $end) = @_;
	
	my $binner = kentBinTools->new();
	unless($self->has_bed()){
		print STDERR "[perlBed] Error! Called data_intersects before loading data!\n";
		exit;
	}
	if(!(exists($self->bed()->{$chr}))){
		return 0;
	}
	
	my @bins = $binner->searchbins($start, $end);
	foreach my $b (@bins){
		if(!(exists($self->bed()->{$chr}->{$b}))){
			next;
		}
		foreach my $gc (@{$self->bed()->{$chr}->{$b}}){
			if($binner->overlap($gc->start(), $gc->end(), $start, $end) > 1){
				return 1;
			}
		}
	}
	return 0;
}

__PACKAGE__->meta->make_immutable;

1;