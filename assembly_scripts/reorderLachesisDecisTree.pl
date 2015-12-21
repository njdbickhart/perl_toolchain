#!/usr/bin/perl
# This script is designed to take the Lachesis scaffold decision file and turn it into an agp file

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -t <tab decision file> -f <input fasta> -i <input fai> -o <output basename>\n";
getopt('tfio', \%opts);

unless(defined($opts{'i'}) && defined($opts{'f'}) && defined($opts{'t'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $worker = LachWorker->new('fasta' => $opts{'f'});
$worker->getChrLens($opts{'i'});

print STDERR "Working on cluster addition\n";
$worker->generateClusters($opts{'t'});

print STDERR "Generating output\n";
$worker->generateOutput($opts{'o'} . ".fa", $opts{'o'} . ".agp");


exit;

BEGIN{
package LachWorker;
use Mouse;
use namespace::autoclean;
use Clone;

# hash: {chr num} -> [] -> lachClusters
has 'lachclusters' => (is => 'rw', isa => 'HashRef[Any]');
has 'chrlens' => (is => 'rw', isa => 'HashRef[Any]');
has ['fasta'] => (is => 'ro', isa => 'Str', required => 1);

sub getChrLens{
	my($self, $infai) = @_;
	my %hash;
	open(my $IN, "< $infai");
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$hash{$segs[0]} = $segs[1];
	}
	close $IN;
	$self->chrlens(\%hash);
}

sub generateOutput{
	my($self, $outfa, $outagp) = @_;
	
	open(my $AGP, "> $outagp");
	open(my $FA, "> $outfa");
	my $fasta = $self->fasta();
	
	my %data = %{$self->lachclusters};
	my %lens = %{$self->chrlens};
	my $currentpos = 0;
	foreach my $clust (sort{$a <=> $b} keys(%data)){
		$currentpos = 0;
		my $chrlen = $lens{$clust};
		my $seq;
		my @lclust = sort{$a->lowest <=> $b->lowest} @{$data{$clust}};
		for(my $x = 0; $x < scalar(@lclust); $x++){
			# Should be sorted and ready for printout!
			my @ordered = $lclust[$x]->getOrderedArray();
			for(my $y =0; $y < scalar(@ordered); $y++){
				my $o = $ordered[$y];
				my $end = $currentpos + $chrlen;
				my $partnum = $x + $y;
				my $contig = $o->contig;
				print {$AGP} $o->clustid . "\t" . $currentpos . "\t" . $end . "\t";
				print {$AGP} $partnum . "\tF\t" . $o->contig . "\t0\t$chrlen\t" . $o->lo;
				print {$AGP} "\n";
				
				# Taking care of the fasta
				open(my $SEQ, "samtools faidx $fasta $contig:0-$chrlen |");
				my $h = <$SEQ>;
				
				my $tseq;
				while(my $line = <$SEQ>){
					chomp $line;
					$tseq .= $line;
				}
				close $SEQ;
				
				if($lclust[$x]->revorder){
					$tseq = reverse($tseq);
					$tseq =~ tr/ACGT/TGCA/;
				}
				$seq .= $tseq;
				
				if($y < scalar(@ordered) - 1){
					$currentpos = $end; 
					$end += 5;
					# Now to print the gap
					print {$AGP} $o->clustid . "\t" . $currentpos . "\t" . $end . "\tU\t5\tscaffold\tyes\tmap\n";
					$seq .= "nnnnn";
				}
				$currentpos = $end;
			}
			if($x < scalar(@ordered) - 1){
				$currentpos = $end; 
				$end += 5;
				# Printing a gap if we're not at the end
				print {$AGP} $o->clustid . "\t" . $currentpos . "\t" . $end . "\tU\t5\tscaffold\tyes\tmap\n";
				$currenpos = $end;
			}
		}
		$seq =~ s/(.{1,60})/$1\n/gs;
		print {$FA} $seq;
	}
	close $AGP;
	close $FA;
}

sub generateClusters{
	my ($self, $tab) = @_;
	my %final;
	my $lastchr = -1;
	my $lastmember; my $currentclust;
	open(my $IN, "< $tab");
	my $head = <$IN>;
	while(my $line = <$IN>){
		chomp $line;
		$line =~ s/\r//g;
		my @segs = split(/\t/, $line);
		my $remove = ($segs[8] eq "remove")? 1 : 0;
		my $prob = ($segs[11] eq "problem")? 1 : 0;
		$lastmember = LachLine->new('clustid' => $segs[0], 
					'contig' => $segs[1],
					'rhc' => $segs[3],
					'rhp' => $segs[4],
					'lo' => $segs[5],
					'rho' => $segs[6]);
		if($lastchr == -1 && !$remove){			
			$currentclust = LachCluster->new('clustid' => $segs[0]);
			$currentclust->add($lastmember);
			$lastchr = $segs[0];
		}elsif($lastchr ne $segs[0] || $prob){
			if(!remove){
				$currentclust->add($lastmember)
			}
			push(@{$final{$currentclust->clustid}}, clone($currentclust));
			$currentclust = LachCluster->new('clustid' => $segs[0]);
			$lastchr = $segs[0];
		}else{
			if(!remove){
				$currentclust->add($lastmember);
			}
			$lastchr = $segs[0];
		}
	}
	push(@{$final{$currentclust->clustid}}, clone($currentclust));
	close $IN;
}
				
		

__PACKAGE__->meta->make_immutable;

package LachCluster;
use Mouse;
use namespace::autoclean;

has 'clustid' => (is => 'ro', isa => 'Num', required => 1);
has 'members' => (traits => ['Array'], is => 'rw', isa => 'ArrayRef[LachLine]', default => sub{[]},
	handles => {
		'add' => 'push',
		'get' => 'get',
	});
	
has 'lowest' => (is => 'rw', isa => 'Int', builder => '_findlowest', lazy => 1);
has 'revorder' => (is => 'rw', isa => 'Int', builder => '_getorder', lazy => 1);

sub getOrderedArray{
	my ($self) = @_;
	my @values;
	my $reverse = $self->revorder;
	
	my @array = @{$self->members};	
	if($reverse){
		@array = reverse(@array);
	}
	
	for(my $x = 0; $x < scalar(@array); $x++){
		if($reverse){
			$array[$x]->lo(($array->lo eq "+")? "-" : "+");
		}
		push(@values, $array[$x]);
	}
	return @values;
}

sub _getorder{
	my ($self) = @_;
	my @nums; my @orients;
	foreach my $v (@{$self->members}){
		if($v->rhp eq '.'){
			next;
		}else{
			push(@nums, $v->rhp);
			push(@orients, $v->lo . ";" . $v->rho);
		}
	}
	
	if(scalar(@nums) == 1){
		@b = split(";", $orients[0]);
		if($b[0] ne $b[1]){
			return 0; # needs to be reversed
		}else{
			return 1; # in the same order
		}
	}else{
		if($nums[1] - $nums[0] >= 0){
			return 0; # in the right order
		}else{
			return 1; # in the wrong order compared to rh
		}
	}
}

sub _findlowest{
	my ($self) = @_;
	my $min = 500;
	foreach my $v (@{$self->members}){
		if($v->rhp ne "." && $v->rhp < $min){
			$min = $v->rhp;
		}
	}
	return $min;
}

__PACKAGE__->meta->make_immutable;

package LachLine;
use Mouse;
use namespace::autoclean;

has ['clustid'] => (is => 'ro', isa => 'Int');
has ['rhc', 'rhp'] => (is => 'ro', isa => 'Any');
has ['contig', 'lo', 'rho'] => (is => 'rw', isa => 'Str');


__PACKAGE__->meta->make_immutable;

}