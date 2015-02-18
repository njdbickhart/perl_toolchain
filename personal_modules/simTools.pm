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
use PerlBed;
use kentBinTools;
use namespace::autoclean;

has ['RefFile', 'CompFile'] => (is => 'ro', isa => 'CompareFile', required => 1);
has ['RefVars', 'TruePos', 'CompVars', 'CloseCall'] => (is => 'rw', isa => 'Int', default => 0);

#NOTE: must input reference column keys and one of them must be "altAllele"
sub CalculateComp{
	my ($self, $refcolkeys) = @_;
	my $refContainer = $self->RefFile->CreateBedContainer($refcolkeys);
	$self->RefVars($refContainer->lines);
	
	# Read comp file and determine bed coordinates
	my $comp = $self->CompFile->file;
	my $chrcol = $self->CompFile->chrCol;
	my $startcol = $self->CompFile->startCol;
	my $endcol = $self->CompFile->endCol;
	my $altAlleleCol = $self->CompFile->altAlleleCol;
	open(IN, "< $comp") || die "[BENCHCOMP] Could not open comparison file!\n";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		my $retbed = $refContainer->firstIntersect($segs[$chrcol], $segs[$startcol], $segs[$endcol]);
		
		if($retbed != 0){
			# There was an intersection
			my $refalt = $retbed->get_value("altAllele");
			if($refalt eq $segs[$altAlleleCol]){
				$self->TruePos($self->TruePos + 1);
			}else{
				$self->CloseCall($self->CloseCall + 1);
			}
		}
		$self->CompVars($self->CompVars + 1);
	}
	close IN;
	
	print STDERR "[BENCHCOMP] Completed comparison\n";
}

sub CreateStats{
	my ($self, $outfile) = @_;
	my $falsePos = $self->CompVars - $self->TruePos;
	my $truePos = $self->TruePos;
	my $falseNeg = abs($self->TruePos - $self->RefVars);
	my $knownPos = $self->RefVars;
	my $testCalls = $self->CompVars;
	
	my $sens = $truePos / $knownPos;
	my $pres = $truePos / $testCalls;
	
	print "Reference known variants: $knownPos\tTest true positives: $truePos\n";
	print "Test total calls: $testCalls\tTest false positives: $falsePos\tTest false negatives: $falseNeg\n";
	print "Sensitivity: $sens\n";
	print "Precision: $pres\n";
	if(defined($outfile)){
		open(OUT, "> $outfile");
		print OUT "Reference known variants: $knownPos\tTest true positives: $truePos\n";
		print OUT "Test total calls: $testCalls\tTest false positives: $falsePos\tTest false negatives: $falseNeg\n";
		print OUT "Sensitivity: $sens\n";
		print OUT "Precision: $pres\n";
		close OUT;
	}
}

__PACKAGE__->meta->make_immutable;

package CompareFile;
use Mouse;
use strict;
use perlBed;
use namespace::autoclean;

has 'file' => (is => 'ro', isa => 'Str', required => 1);
has ['chrCol', 'startCol', 'endCol', 'altAlleleCol'] => (is => 'ro', isa => 'Int', required => 1);

# Assumes that the file is a typical bed file
sub CreateBedContainer{
	my ($self, $refcolkeys) = @_;
	print STDERR "[COMPARISON] Creating bed container...\n";
	my $comp = BedContainer->new();
	$comp->loadFile($self->file, $refcolkeys);
	return $comp;
}

__PACKAGE__->meta->make_immutable;

1;