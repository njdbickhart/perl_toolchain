#!/usr/bin/perl
# This script is designed to process a divet file and identify putative transchromosomal events

use strict;
use Getopt::Std;
use Class::Struct;
use kentBinTools;

struct('bed' => {
	chr => '%',
});

struct('binbed' => {
	bin => '%',
});

struct('simpbed' => {
	chr1 => '$',
	start1 => '$',
	end1 => '$',
	chr2 => '$',
	start2 => '$',
	end2 => '$',
	support => '$',
});

my $usage = "perl $0 -d <divet file> -o <output bed file>\n";
my %opts;

getopt('do', \%opts);

unless(defined($opts{'d'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $binner = kentBinTools->new();
my $beds = bed->new();

open(IN, "grep \'transchr\' $opts{d} |") || die "Could not open divet file!\n";
my $last = "none";
my @store;
my $linenum = 0;
while (my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[0] ne $last && $last ne "none"){
		# The read name has changed, so there are likely no more alternative mappings
		if(scalar(@store) == 1){
			# There is only one mapping, process this position
			processRead($beds, $store[0], $binner);
		}
		@store = ();
		push(@store, \@segs);
		$last = $segs[0];
	}else{
		push(@store, \@segs);
		$last = $segs[0];
	}
	$linenum++;
	if($linenum % 100000 == 0){
		print STDERR "Finished with line:\t $linenum\r";
	}
}
print STDERR "\n";

close IN;
open(OUT, "> $opts{o}");
foreach my $chr (sort {my ($x) = $a =~ /chr(.+)/;
			my ($y) = $b =~ /chr(.+)/;
			if($x eq "X"){
				$x = 500;
			}elsif($x eq "Y"){
				$x = 501;
			}elsif($x eq "M" || $x eq "MT"){
				$x = 502;
			}
			
			if($y eq "X"){
				$y = 500;
			}elsif($y eq "Y"){
				$y = 501;
			}elsif($y eq "M" || $x eq "MT"){
				$y = 502;
			}
			$x <=> $y} keys(%{$beds->chr})){
			foreach my $bin (sort{$a <=> $b} keys(%{$beds->chr($chr)->bin})){
				foreach my $row (sort {$a->start1 <=> $b->start1} @{$beds->chr($chr)->bin()->{$bin}}){
					print OUT $row->chr1 . "\t" . $row->start1 . "\t" . $row->end1 . "\t" . $row->chr2 . "\t";
					print OUT $row->start2 . "\t" . $row->end2 . "\t" . $row->support . "\n";
				}
			}
}
close OUT;
exit;

sub processRead{
	my ($beds, $segs, $binner) = @_;
	# Always sort by the least chr
	my ($leastchr, $mostchr) = ($segs->[1] lt $segs->[5])? ($segs->[1], $segs->[5]) : ($segs->[5], $segs->[1]);
	my ($leaststart, $moststart) = ($segs->[1] eq $leastchr)? ($segs->[2], $segs->[6]) : ($segs->[6], $segs->[2]);
	my ($leastend, $mostend) = ($segs->[1] eq $leastchr) ? ($segs->[3], $segs->[7]) : ($segs->[7], $segs->[3]);
	my $bin = $binner->getbin($leaststart, $leastend);
	
	if(!exists($beds->chr()->{$leastchr})){
		$beds->chr($leastchr, binbed->new());
	}
	if(!exists($beds->chr($leastchr)->bin()->{$bin})){
		$beds->chr($leastchr)->bin($bin, []);
		push(@{$beds->chr($leastchr)->bin($bin)}, createBed($leastchr, $leaststart, $leastend, $mostchr, $moststart, $mostend));
		return;
	}
	
	my $searchbins = $binner->searchbins($leaststart, $leastend);
	my $found = 0;
	foreach my $b (@{$beds->chr($leastchr)->bin()->{$bin}}){
		if($b->chr2 ne $mostchr){next;}
		my $overlap1 = $binner->overlap($b->start1, $b->end1, $leaststart, $leastend);
		my $overlap2 = $binner->overlap($b->start2, $b->end2, $moststart, $mostend);
		
		if($overlap1 > 0 && $overlap2 > 0){
			# To speed things up, I will not refine the coordinates for now
			$b->support($b->support + 1);
			$found = 1;
			last;
		}
	}
	if(!$found){
		push(@{$beds->chr($leastchr)->bin()->{$bin}}, createBed($leastchr, $leaststart, $leastend, $mostchr, $moststart, $mostend));
	}
}

sub createBed{
	my ($leastchr, $leaststart, $leastend, $mostchr, $moststart, $mostend) = @_;
	return simpbed->new('chr1' => $leastchr, 'start1' => $leaststart, 'end1' => $leastend,
		'chr2' => $mostchr, 'start2' => $moststart, 'end2' => $mostend, 'support' => 0);
}