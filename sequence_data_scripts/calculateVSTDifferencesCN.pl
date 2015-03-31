#!/usr/bin/perl
# I'm going to have to calculate my own Vst method
# Damn. Well I have some metrics to go off of here.
# Pop strucure list: first column has the animal name, second column has the population number
# Gene list is the tab-delimited text output from AnnotateUsingGenomicInfo.jar

use strict;
use Getopt::Std;

my %opts;
getopt('cpos', \%opts);

my $usage = "$0 -c <gene list> -p <pop structure list> -o <output file name>
	-c 	gene list from AnnotateUsingGenomicInfo 
	-p	population structure list (animalname(tab)populationname) 
	-o 	output file name 
	-s 	[skip list, optional, comma separated]\n";

unless(defined($opts{'c'}) && defined($opts{'p'}) && defined($opts{'o'})){
	print "Missing mandatory arguments!\n";
	print $usage;
	exit;
}

my %skip;
if(defined($opts{'s'})){
	my @ssegs = split(/,/, $opts{'s'});
	foreach my $s (@ssegs){
		$skip{$s} = 1;
	}
}

my %popstruct;
my %poplookup;
open(POP, "< $opts{p}") || die "Could not open POP gen file!\n$usage";
while(my $line = <POP>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if(exists($skip{$segs[0]})){next;}
	push(@{$popstruct{$segs[1]}}, $segs[0]);
	$poplookup{$segs[0]} = $segs[1];
}
close POP;

open(CN, "< $opts{c}") || die "Could not open CNVR file!\n$usage";
my $h = <CN>;
my @cols = split(/\t/, $h);
my %genoholder; # {pop}->{animal} ->{marker} = value
my %totalholder; # {markername} -> [] = values
my %genenames; # {markername} -> [] = genes
my @markers;
while(my $line = <CN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my $m = "$segs[1]:$segs[2]-$segs[3]";
	push(@markers, "$segs[1]:$segs[2]-$segs[3]");
	push(@{$genenames{$m}}, $segs[0]);
	for(my $x = 6; $x < scalar(@segs); $x++){
		my $an = $cols[$x];
		if(!exists($poplookup{$an})){next;}
		my $pop = $poplookup{$an};
		
		#$genoholder{$pop}->{$an}->{$m} = $segs[$x];
		#push(@{$totalholder{$m}}, $segs[$x]);
		$genoholder{$pop}->{$an}->{$m} = int(($segs[$x] * 10) + 0.5);
		push(@{$totalholder{$m}}, int(($segs[$x] * 10) + 0.5));
		
	}
}
close CN;

open(OUT, "> $opts{o}");
# Working on single genes at the time 
foreach my $m (@markers){
	print STDERR "Working on marker: $m\r";
	my $vst = Vst(\%totalholder, \%genoholder, [$m], \%skip);
	my $genes = join(";", @{$genenames{$m}});
	my ($chr, $start, $end) = $m =~ m/(chr.+):(\d+)-(\d+)/;
	print OUT "$chr\t$start\t$end\t$vst\t$genes\n";
}
close OUT;
print "\n";
exit;

sub Vst{
	my ($totalholder, $genoholder, $markers, $skip) = @_;
	my @total;
	foreach my $m (@{$markers}){
		push(@total, @{$totalholder{$m}});
	}
	
	my $Vt = variance(\@total);
	
	my %popnums;
	my %popvars;
	my $numerator;
	my $denominator;
	foreach my $p (keys(%{$genoholder})){
		foreach my $an (keys %{$genoholder{$p}}){
			foreach my $m (@{$markers}){
				if(exists($genoholder{$p}->{$an}->{$m})){
					if(exists($skip->{$an})){next;}
					$popnums{$p} += 1;
					push(@{$popvars{$p}}, $genoholder{$p}->{$an}->{$m});
				}
			}
		}
		my $pvar = variance($popvars{$p});
		$numerator += $pvar * $popnums{$p}; 
		$denominator += $popnums{$p};
	}
	
	my $Vs = $numerator / $denominator;
	if($Vt == 0){return 0;}
	return (($Vt - $Vs) / $Vt);
}
sub variance{
	my ($data) = @_;
	if(scalar(@{$data}) <= 1){
		return 0;
	}
	
	my $sum = 0;
	foreach my $b (@{$data}){
		$sum += $b;
	}
	my $avg = $sum / (scalar(@{$data}));
	
	my $sumsquares = 0;
	foreach my $b (@{$data}){
		$sumsquares += (($b - $avg) * ($b - $avg));
	}
	return ($sumsquares / (scalar(@{$data}) - 1));
}

#sub variance{
 #   my ($data) = @_;
 #   if (@$data ==1) {
 #       return 0;
 #   }
 #   my $totalD2 = 0;
 #   foreach (@$data) {
 #      $totalD2 += $_*100
 #   }
 #   my $meanD2 = $totalD2 / (scalar @$data);
 #
 #   my $sqtotalD4 = 0;
 #   foreach (@$data) {
 #       $sqtotalD4 += ($_*100 - $meanD2) ** 2;
 #   }
 #   my $varD4 = $sqtotalD4 / (scalar @$data - 1);
 #   return $varD4/10_000; # convert from fixed point to floating point
#}