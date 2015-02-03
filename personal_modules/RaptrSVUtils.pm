#!/usr/bin/perl
# This is a collection of utilities designed to process RAPTR-SV files

package SVEntry;
use Mouse;
use Mouse::Util::TypeConstraints;
use strict;
use namespace::autoclean;

has 'line' => (is => 'rw', isa => 'Str', required => 1, clearer => 'remove_str');

has 'chr' => (is => 'rw', isa => 'Str');
has ['outStart', 'inStart', 'inEnd', 'outEnd', 'divSup', 'splitSup', 'unbalSup'] => (is => 'rw', isa => 'Int', default => 0);
has 'totSupp' => (is => 'rw', isa => 'Num');
has 'svType' => (is => 'rw', isa => enum([qw[ DELETION INSERTION TANDEM ]]));


# Gets the raw number of reads supporting this event
sub getRawReadCount{
	my ($self) = @_;
	return $self->divSup + $self->splitSup + $self->unbalSup;
}

# calculates the total support divided by the total number of reads (divet, split and unbal) 
sub getSupportRatio{
	my ($self) = @_;
	my $countreads = $self->getRawReadCount;
	return $self->totSupp / $countreads;
}

# Returns a string that is a simple bed format standard (without the newline!)
sub printBed{
	my ($self) = @_;
	
	if($self->svType eq "TANDEM" || $self->svType eq "INSERTION"){
		my @elements = ($self->chr, $self->outStart, $self->outEnd, $self->svType, $self->totSupp);
		return join("\t", @elements);
	}elsif($self->svType eq "DELETION"){
		my $start = ($self->inStart > $self->inEnd)? $self->inEnd : $self->inStart;
		my $end = ($self->inStart > $self->inEnd)? $self->inStart : $self->inEnd;
		my @elements = ($self->chr, $start, $end, $self->svType, $self->totSupp);
		return join("\t", @elements);
	}
}

before 'chr', 'outStart', 'inStart', 'inEnd', 'outEnd', 'divSup', 'splitSup', 'totSupp', 'svType' => sub {
	my ($self) = @_;
	my @segs = split(/\t/, $self->line());
	
	$self->chr($segs[0]);
	$self->outStart($segs[1]);
	$self->inStart($segs[2]);
	$self->inEnd($segs[3]);
	$self->outEnd($segs[4]);
	$self->svType($segs[5]);
	$self->divSup($segs[6]);
	$self->splitSup($segs[7]);
	$self->unbalSup($segs[8]);
	$self->totSupp($segs[9]);
	
	$self->remove_str();
};

__PACKAGE__->meta->make_immutable;

1;