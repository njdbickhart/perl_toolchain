#!/usr/bin/perl
# This script is designed to combine the results output of EBSeq into a singular tab delimited file
# Here is the output plan I'm looking to generate:
# GeneID	PPDE	BestPattern	ProbPattern	C1rawcount	...	CNrawcount	pattern	C2vsC1FC	...
# Something similar to what I remember from the CLCBio output


use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -c <condmeans file> -f <base EBSeq output file> -o <output tab file>\n";

getopt('cfo', \%opts);

unless(defined($opts{'c'}) && defined($opts{'f'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

open(COND, "< $opts{c}") || die "Could not open condmeans file!\n";
open(BASE, "< $opts{f}") || die "Could not open base output file!\n";
open(OUT, "> $opts{o}"); 

my $temp = <COND>;
chomp($temp);
$temp =~ s/\"//g;
my @condhead = split(/\t/, $temp);
$temp = <BASE>;
chomp($temp);
$temp =~ s/\"//g;
my @outhead = split(/\t/, $temp);

my $printhead = 1;
while(my $line = <COND>){
	my $base = <BASE>;
	chomp($line, $base);
	$line =~ s/\"//g;
	$base =~ s/\"//g;
	
	my @condsegs = split(/\t/, $line);
	my @outsegs = split(/\t/, $base);
	
	for(my $x = 1; $x < scalar(@condsegs); $x++){
		$condsegs[$x] =~ s/NA/0/g;
	}
	for(my $x = 1; $x < scalar(@outsegs); $x++){
		$outsegs[$x] =~ s/NA/0/g;
	}
	
	my $worker = EBseqCont->new();
	$worker->Format(\@condsegs, \@outsegs, \@condhead, \@outhead);
	
	if($printhead){
		my @headersegs = $worker->GetHeader(\@condhead);
		print OUT join("\t", @headersegs);
		print OUT "\n";
		$printhead = 0;
	}
	
	my @results = $worker->GetArrayOut();
	print OUT join("\t", @results);
	print OUT "\n";
}

close OUT;

exit;

BEGIN{
# This package is just a formatting workhorse
# It takes two arrays from the two Ebseq output files and generates another array
package EBseqCont;
use Mouse;
use namespace::autoclean;

has ['geneid', 'bestpat'] => (is => 'rw', isa => 'Str');
has ['ppde', 'probpat'] => (is => 'rw', isa => 'Num');
has 'expcounts' => (traits => ['Array'], is => 'rw', isa => 'ArrayRef[Any]', default => sub{[]},
	handles => {
		'addexp' => 'push',
		'getexp' => 'elements',
	});
has 'fcdata' => (traits => ['Array'], is => 'rw', isa => 'ArrayRef[FCCalc]', default => sub{[]},
	handles => {
		'addfc' => 'push',
		'getarray' => 'elements',
	});

sub Format{
	my ($self, $condref, $outref, $condhead, $outhead) = @_;
	# condhead and outhead are array refs of the headers of each file
	# outhead: "Pattern1"      "Pattern2"      "Pattern3"      "Pattern4"      "Pattern5"      "MAP"   "PPDE"
	# condhead: "C1"    "C2"    "C3"
	
	$self->geneid($outref->[0]);
	$self->bestpat($outref->[6]);
	$self->ppde($outref->[7]);
	for(my $x = 0; $x < scalar(@{$outhead}); $x++){
		if($outhead->[$x] eq $outref->[6]){
			$self->probpat($outref->[$x + 1]);
			last;
		}
	}
	
	for(my $x = 1; $x < scalar(@{$condref}); $x++){
		$self->addexp($condref->[$x]);
	}
	
	# Pairwise comparison loops -- always going to try to be alternating
	# first comp: c2 vs c1, second c3 vs c1, third c3 vs c2
	for(my $x = 1; $x < scalar(@{$condref}) - 1; $x++){
		for(my $y = $x +1; $y < scalar(@{$condref}); $y++){
			$self->addfc(FCCalc->new(
				'c1c' => $condref->[$x],
				'c2c' => $condref->[$y],
				'c1name' => $condhead->[$x - 1],
				'c2name' => $condhead->[$y - 1])
			);
		}
	}
}

sub GetArrayOut{
	my ($self) = @_;
	# GeneID	PPDE	BestPattern	ProbPattern	C1rawcount	...	CNrawcount	C2vsC1FC	...
	my @results;
	push(@results, ($self->geneid, $self->ppde, $self->bestpat, $self->probpat));
	
	foreach my $e ($self->getexp){
		push(@results, $e);
	}
	
	foreach my $fc ($self->getarray){
		push(@results, $fc->getFC);
	}
	return @results;
}

sub GetHeader{
	my ($self, $condhead) = @_;
	my @outputhead = ("GeneID", "PPDE", "BestPattern", "ProbPattern");
	
	foreach my $c (@{$condhead}){
		push(@outputhead, "$c" . 'NormCount');
	}
	
	foreach my $fc($self->getarray){
		push(@outputhead, $fc->getTitle);
	}
	return @outputhead;
}

__PACKAGE__->meta->make_immutable;

package FCCalc;
use Mouse;
use namespace::autoclean;
# FC calc is always going to be the first condition in the input divided by the second

has ['c1c', 'c2c'] => (is => 'ro', isa => 'Num', required => 1);
has ['c1name', 'c2name'] => (is => 'ro', isa => 'Str', required => 1);

sub getFC{
	my ($self) = @_;
	return ($self->c2c / $self->c1c);
}

sub getTitle{
	my ($self) = @_;
	return $self->c2name . "_vs_" . $self->c1name;
}

__PACKAGE__->meta->make_immutable;
}