#!/usr/bin/perl
# This script creates genomic windows, counts CNVRs within them and then identifies windows that are significantly enriched for CNVRs

use strict;
use Getopt::Std;
use kentBinTools;

my %opts;
getopt('cbsiof', \%opts);
my $usage = "perl $0 (-c <CNVR tab file> OR -b <Bed file list>) -s <Sliding window size (bp)> -i <Sliding window increment (bp)> -o <Output bed file> -f <samtools faidx file>\n";

unless((defined($opts{'c'}) || defined($opts{'b'})) && defined($opts{'s'}) && defined($opts{'i'}) && defined($opts{'o'}) && defined($opts{'f'})){
	print $usage;
	exit;
}

# Create windows
my %windows;
#open(TEST, "> $opts{o}.wins");
open(IN, "< $opts{f}") || die "Could not open faidx file!\n";
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	print "Windows for chr: $segs[0]\r";
	for(my $x = 0; $x < $segs[1]; $x += $opts{'i'}){
		my $end;
		if($x + $opts{'s'} > $segs[1]){
			$end = $segs[1];
		}else{
			$end = $x + $opts{'s'};
		}
		my $bin = kentBinTools->getbin($x, $end);
		
		push(@{$windows{$segs[0]}->{$bin}}, cnvWindow->new(
			'chr' => $segs[0],
			'start' => $x,
			'end' => $end));
		#print TEST "$segs[0]\t$x\t$end\n";
	}
}
close IN;
#close TEST;
print "\n";

# Read in CNVRs
if(defined($opts{'c'})){
	print "Reading CNVR tab file...\n";
	open(IN, "< $opts{c}") || die "Could not open CNVR tab file!\n$usage";
	my $head = <IN>;
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		windowOverlapper(\%windows, $segs[1], $segs[2], $segs[3]);
	}
	close IN;
}elsif(defined($opts{'b'})){
	print "Reading Bed file list...\n";
	open(IN, "< $opts{b}") || die "Could not open bed file!\n$usage";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		open(BED, "< $segs[0]");
		while(my $l = <BED>){
			chomp $l;
			my @bsegs = split(/\t/, $l);
			windowOverlapper(\%windows, $bsegs[0], $bsegs[1], $bsegs[2]);
		}
		close BED;
	}
	close IN;
}

# Calculate Average and Stdev
my $total = 0;
my $sum = 0;
foreach my $c (keys(%windows)){
	foreach my $b (keys(%{$windows{$c}})){
		foreach my $w (@{$windows{$c}->{$b}}){
			$total++;
			$sum += $w->count();
		}
	}
}
my $avg = $sum / $total;

my $sumsquares = 0;
foreach my $c (keys(%windows)){
	foreach my $b (keys(%{$windows{$c}})){
		foreach my $w (@{$windows{$c}->{$b}}){
			$sumsquares += ($avg - $w->count())**2;
		}
	}
}
my $stdev = sqrt($sumsquares / $total);
print "Calculated average: $avg\tCalculated Stdev: $stdev\n";

my $thresh = $avg + ($stdev * 3);
# Identify windows that are 3 fold higher than the average and print them out
open(OUT, "> $opts{o}");
my $numprint = 0;
foreach my $c (keys(%windows)){
	foreach my $b (keys(%{$windows{$c}})){
		foreach my $w (@{$windows{$c}->{$b}}){
			if($w->count() >= $thresh){
				print OUT $w->printOut();
				$numprint++;
			}
		}
	}
}
close OUT;
print "Found $numprint windows above the threshold of $thresh. Total windows: $total\n";

exit;

sub standardDeviation {
	my $array_ref = shift(@_);
	my $avg = shift(@_);
	my @numbers = @{$array_ref};
	return 0 unless(scalar(@numbers));
	
	my $total2 = 0;
	foreach my $num (@numbers) {
	$total2 += ($avg-$num)**2;
	}
	my $mean2 = $total2 / (scalar @numbers);
	
	my $std_dev = sqrt($mean2);
	return $std_dev;
}
sub average {
	my $array_r = shift(@_);
	my @numb = @{$array_r};
	if (scalar(@numb) == 0){
		return 0;
	}
	my $total3 = 0;
	foreach my $num1 (@numb) {
	$total3 += $num1;
	}
	my $mean3 = $total3 / (scalar @numb);
	return $mean3;
}

sub windowOverlapper{
	my ($wref, $c, $s, $e) = @_;
	my @bins = kentBinTools->searchbins($s, $e);
	foreach my $b (@bins){
		foreach my $wins (@{$windows{$c}->{$b}}){
			my $ovlp = kentBinTools->overlap($wins->start(), $wins->end(), $s, $e);
			if($ovlp > 0){
				$wins->inc();
			}
		}
	}
	#print "done\n";
}

BEGIN{
	package cnvWindow;
	use Mouse;
	use namespace::autoclean;
	use strict;
	
	has ['start', 'end'] => (
		is => 'ro',
		isa => 'Int',
	);
	
	has 'chr' => (
		is => 'ro',
		isa => 'Str',
	);
	
	has 'count' => (
		traits => ['Counter'],
		is => 'rw',
		isa => 'Int',
		default => 0,
		handles => {
			'inc' => 'inc',
		}
	);
	
	sub printOut{
		my $self = shift;
		return $self->chr() . "\t" . $self->start() . "\t" . $self->end() . "\n";
	}
	__PACKAGE__->meta->make_immutable();
}