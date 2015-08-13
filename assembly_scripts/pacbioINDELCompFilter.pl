#!/usr/bin/perl
# This script processes a samtools-HTSLIB vcf file of Illumina reads against a reference genome
# The goal is to categorize and filter INDELs found within the file for assembly correction

use strict;
use Getopt::Std;

my $usage = "perl $0 -v <input vcf, gzipped> -o <filtered vcf output>\n";
my %opts;
getopt('vo', \%opts);

unless(defined($opts{'v'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $counter = contigINDELCounter->new();
open(my $IN, "gunzip -c $opts{v} | ");
open(my $OUT, "> $opts{o}");
while(my $line = <$IN>){
	chomp $line;
	if($line =~ /^#/){
		print $OUT "$line\n";
	}elsif($line =~ /INDEL/){
		my @segs = split(/\t/, $line);
		
		# Filtering criteria:
		# 1. skip heterozygotes
		# 2. skip INDELs with IMF <= 0.5
		# 3. skip INDELs with IDV <= 2
		
		if($segs[-1] =~ /^0\/.{1}\:/){
			# Has one copy of the reference allele
			next;
		}
		
		# Store all INFO tags into a hash, minus the INDEL tag
		my %info = map {split(/=/, $_, 2)} grep(/=/, split(/;/, $segs[7]));
		if($info{'IMF'} <= 0.5 || $info{'IDV'} <= 2){
			next;
		}
		
		$counter->collect($segs[0], $segs[1], $segs[3], $segs[4]);
		print $OUT "$line\n";
	}
}
close $IN;
close $OUT;

$counter->printOut;
	

exit;

BEGIN{
package contigINDELCounter;
use Mouse;
use namespace::autoclean;

has 'indlens' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'contigs' => 'keys',
		'has_contig' => 'exists',
		'get_contig' => 'get',
		'set_contig' => 'set',
	});
	
has 'contig' => (is => 'rw', isa => 'Str', default => 'NA');
has 'laspos' => (is => 'rw', isa => 'Int', default => 0);
has ['totINDEL', 'simpINS', 'simpDEL', 'simpHomINS', 'simpHomDEL', 'compINS', 'compDEL'] => (is => 'rw', isa => 'Counter', default => sub{Counter->new()});

sub collect{
	my ($self, $contig, $pos, $ref, $alt) = @_;
	if($contig ne $self->contig){
		# Reset contig position count
		$self->laspos(0);
		$self->contig($contig);
		if(!$self->has_contig($contig)){
			$self->set_contig($contig, []);
		}
	}
	
	my $temp = $self->get_contig($contig);
	push(@{$temp}, $pos - $self->laspos);
	$self->laspos($pos);
	
	$self->totINDEL->inc;
	# Preemptive homopolymeric check
	my @refsegs = split(//, $ref);
	my @altsegs = split(//, $alt);
	my $homopoly = ($refsegs[-1] eq $altsegs[-1])? 1 : 0;
	
	# determine type of INDEL
	if(length($alt) - length($ref) == 1){
		# Simple INS
		if($homopoly){
			$self->simpHomINS->inc;
		}else{
			$self->simpINS->inc;
		}
	}elsif(length($ref) - length($alt) == 1){
		# Simple DEL
		if($homopoly){
			$self->simpHomDEL->inc;
		}else{
			$self->simpDEL->inc;
		}
	}else{
		# Complex variant
		if(length($alt) - length($ref) > 0){
			# INS
			$self->compINS->inc;
		}else{
			$self->compDEL->inc;
		}
	}
	$self->set_contig($contig, $temp);
}

sub printOut{
	my ($self) = @_;
	my @totallens;	
	
	print "Simple_INS\t" . $self->simpINS->value . "\n";
	print "Simple_DEL\t" . $self->simpDEL->value . "\n";
	print "Simple_HOM_INS\t" . $self->simpHomINS->value . "\n";
	print "Simple_HOM_DEL\t" . $self->simpHomDEL->value . "\n";
	print "Complex_INS\t" . $self->compINS->value . "\n";
	print "Complex_DEL\t" . $self->compDEL->value . "\n";
	print "Total_INDELS\t" . $self->totINDEL->value . "\n";
	
	print "---\n";
	print "Contig\tAvg_INDEL_Distance\n";
	
	foreach my $ctg (sort {$a cmp $b} $self->contigs){
		my $avg = $self->average($self->get_contig($ctg));
		push(@totallens, @{$self->get_contig($ctg)});
		print "$ctg\t$avg\n";
	}
	
	my $totavg = $self->average(\@totallens);
	print "TOTAL\t$totavg\n";
}
	

sub average{
	my ($self, $array) = @_;
	if(scalar(@{$array}) == 0){
		return 0;
	}
	my $sum = 0;
	foreach my $a (@{$array}){
		$sum += $a;
	}
	
	return ($sum / scalar(@{$array}));
}

__PACKAGE__->meta->make_immutable;

package Counter;
use Mouse;
use namespace::autoclean;

has 'value' => (is => 'rw', isa => 'Int', default => 0);

sub inc{
	my ($self) = @_;
	$self->value($self->value + 1);
}

__PACKAGE__->meta->make_immutable;

}
