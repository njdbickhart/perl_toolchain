#!/usr/bin/perl
# This is a collection of handy utilities designed to make memory management and time tracking easier

package BenchMarker;
use Mouse;
use strict;
use namespace::autoclean;
use Benchmark qw(:all) ;

has 'cmd' => (is => 'ro', isa => 'ArrayRef[Str]', required => 1);
has 'time' => (is => 'rw', isa => 'Num');
has 'results' => (
	traits => ['Array'], 
	is => 'rw', 
	isa => 'ArrayRef[Any]',
	default => sub{[]},
	handles => {
		'push_results' => 'push',
		'count_results' => 'count',
	},
);

sub run{
	my ($self) = @_;
	foreach  my $c (@{$self->cmd}){
		my $bobject = timethis(1, system("$self->cmd"));
		$self->push_results($bobject);
	}
}

sub formatOut{
	my ($self) = @_;
	my @results;
	
	for(my $x = 0; $x < $self->count_results; $x++){
		my $cmd = $self->cmd->[$x];
		my $b  = $self->results->[$x];
		push(@results, "Command:\t$cmd\tCPU time:\t" . $b->cpu_a . "\tWall time:\t" . $b->real );
	}
	return \@results;
}

sub totalTimeOut{
	my ($self) = @_;
	my $cpua = 0;
	my $realt = 0;
	
	for(my $x = 0; $x < $self->count_results; $x++){
		#my $cmd = $self->cmd->[$x];
		my $b  = $self->results->[$x];
		$cpua += $b->cpu_a;
		$real += $b->real;
	}
	
	return "Command:\tTOTAL\tCPU time:\t$cpua\tWall time:\t$realt";
}

__PACKAGE__->meta->make_immutable;


package BenchCompare;
use Mouse;
use strict;
use kentBinTools;
use namespace::autoclean;



__PACKAGE__->meta->make_immutable;