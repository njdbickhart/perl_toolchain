#!/usr/bin/perl
# This script is designed to find frequent transchr events that are likely to indicate repeats rather than transchr

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
	start => '$',
	end => '$',
	name => '$',
});

my $usage = "perl $0 -d <input file transchr bedpe> -f <Repeats bed file> -o <output bed file>\n";
my %opts;

getopt('dof', \%opts);

unless(defined($opts{'d'}) && defined($opts{'f'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $binner = kentBinTools->new();
my $beds = bed->new();
my @transchr;

# load filters

open(IN, "< $opts{f}") || die "Could not open repeat bed file!\n";
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[3] < 3){
		next;
	}
	my $bin = $binner->getbin($segs[1], $segs[2]);
	if(!exists($beds->chr()->{$segs[0]})){
		$beds->chr($segs[0], binbed->new());
	}
	if(!exists($beds->chr($segs[0])->bin()->{$bin})){
		$beds->chr($segs[0])->bin($bin, []);
	}
	push(@{$beds->chr($segs[0])->bin()->{$bin}}, simpbed->new('start' => $segs[1], 'end' => $segs[2]), 'name' => $segs[3]);
}

close IN;

# Now process transchr file to remove entries that overlap likely repeats
open(IN, "< $opts{d}") || die "Could not open input file!\n";
open(OUT, "> $opts{o}");
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	
	my @searchbins1 = $binner->searchbins($segs[1], $segs[2]);
	my @searchbins2 = $binner->searchbins($segs[4], $segs[5]);
	
	my $found1 = 0;
	my $found2 = 0;
	my %names; 
	foreach my $b (@searchbins1){
		if(!exists($beds->chr($segs[0])->bin()->{$b})){next;}
		foreach my $bed (@{$beds->chr($segs[0])->bin()->{$b}}){
			if($binner->overlap($segs[1], $segs[2], $bed->start(), $bed->end())){
				$found1 = 1;
				$names{$bed->name()} = 1;
			}
		}
	}
	
	foreach my $b (@searchbins2){
		if(!exists($beds->chr($segs[0])->bin()->{$b})){next;}
		foreach my $bed (@{$beds->chr($segs[3])->bin()->{$b}}){
			if($binner->overlap($segs[4], $segs[5], $bed->start(), $bed->end())){
				$found2 = 1;
				$names{$bed->name()} = 1;
			}
		}
	}
	
	if(!(!$found1 && !$found2) && scalar(keys(%names)) <= 2){
		print OUT join("\t", @segs) . "\n";
	}
}

close OUT;
close IN;

exit;