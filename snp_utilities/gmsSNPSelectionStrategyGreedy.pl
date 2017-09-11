#!/usr/bin/perl;
# This script selects markers between two target locations based on the algorithm of Matukumali et al.
# input: chr	pos	VCFQS	MAF	Mapqupstream	Mapqdownstream	SNPName

use strict;
use Getopt::Std;

my $usage = "$0 -c <chr to filter> -i <input tab file with gms snps> -s <start> -e <end> -d <minimum SNP distance> -m <number of snps to aim for>\n";
my %opts;
getopt('cisedm', \%opts);

unless(defined($opts{'i'}) && defined($opts{'s'}) && defined($opts{'d'}) && defined($opts{'m'})){
	print $usage;
	exit;
}

my @snps;
open(IN, "< $opts{i}") || die "Could not open input tab file!\n $usage";
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[0] ne $opts{'c'}){next;}
	if($segs[3] < 0.05 || $segs[3] > 0.51){
		next;
	}
	my $snp = snpgms->new('position' => $segs[1], 'qs' => 10**(-1 * $segs[2] / 10), 'freq' => $segs[3], 'gmsup' => ($segs[4] / 100), 'gmsdown' => ($segs[5] / 100));
	push(@snps, $snp);
}

my $selections = greedyRecursiveSelection(\@snps, $opts{'s'}, $opts{'e'}, 1, $opts{'m'}, $opts{'d'}, []);

foreach my $s (sort{$a->position <=> $b->position} @{$selections}){
	print $s->printOut() . "\n";
}

exit;

sub greedyRecursiveSelection{
	my ($snppool, $start, $end, $curselect, $maxselect, $minlength, $selsnps) = @_;
	my $num = @{$snppool};
	if($num <= 0 || $end - $start <= $minlength || $curselect >= $maxselect){
		# termination conditions:
		# 1. no more SNPs in the region
		# 2. The end/start region is smaller than the minimum length between snps
		# 3. We have selected more than enough SNPs for now
		return;
	}
	
	# Calculate all SNP scores
	foreach my $s (@{$snppool}){
		$s->calcScore($start, $end);
	}
	
	# Run through and select the highest snp; add it to pool if it is beyond the minimum length from previous SNP entries
	my $loop = 1;
	my $selection;
	while($loop){
		my $topSNP;
		my $topScore = 0;
		foreach my $s (@{$snppool}){
			if($s->score() > $topScore){
				$topSNP = $s;
				$topScore = $s->score();
			}
		}
		
		if($topScore == 0){
			# We ran out of SNPs and could not bridge the minimum distance; Break
			return;
		}elsif($topSNP->position() - $start < $minlength || $end - $topSNP->position() < $minlength){
			# The SNP is currently within the minimum distance value; set to zero
			$topSNP->score(0);
		}else{
			$selection = $topSNP;
			$loop = 0;
			last;
		}
	}
	
	push(@{$selsnps}, $selection);
	# Create two pools of snps, before and after the selected SNP and recurse the function
	my @spool;
	foreach my $s (@{$snppool}){
		if($s->position == $selection->position){ next;}
		elsif($s->position < $selection->position){
			push(@spool, $s);
		}
	}
	
	my @epool;
	foreach my $s (@{$snppool}){
		if($s->position == $selection->position){ next;}
		elsif($s->position > $selection->position){
			push(@epool, $s);
		}
	}
	
	greedyRecursiveSelection(\@spool, $start, $selection->position, $curselect * 2, $maxselect, $minlength, $selsnps);
	greedyRecursiveSelection(\@epool, $selection->position, $end, $curselect * 2, $maxselect, $minlength, $selsnps);
	
	# if all goes well, return the array with the selected SNPs
	return $selsnps;
}

BEGIN{
package snpgms;
use Mouse;
use namespace::autoclean;

has 'position' => (is => 'ro', isa =>'Int', required => 1);
has 'qs' => (is => 'ro', isa => 'Num', required => 1);
has 'gmsup' => (is => 'ro', isa => 'Num', required => 1);
has 'gmsdown' => (is => 'ro', isa => 'Num', required => 1);
has 'freq' => (is => 'ro', isa => 'Num', required => 1);

has 'score' => (is => 'rw', isa => 'Num');

sub calcScore{
	my ($self, $start, $end) = @_;
	my $len = $end - $start;
	my $score = $self->_max($self->gmsup, $self->gmsdown) * $self->qs * $self->freq * ($len - abs(2 * $self->position() - ($end + $start)));
	$self->score($score);
}

sub printOut{
	my ($self) = @_;
	return $self->position() . "\t" . $self->qs() . "\t" . $self->freq() . "\t" . $self->gmsup() . "\t" . $self->gmsdown();
}

sub _max{
	my ($self, $v1, $v2) = @_;
	return ($v1 > $v2)? $v1 : $v2;
}
__PACKAGE__->meta->make_immutable;
1;
}