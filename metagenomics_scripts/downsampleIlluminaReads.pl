#!/usr/bin/perl
# this script is designed to downsample several Illumina fastq files using a progressive downsampling algorithm
# The goal is to prepare several subset fastq files for MASH profile estimation.
# Downsample thresholds will be in 5% intervals

use strict;
use FileHandle;

my $usage = "perl $0 <fastq file 1> ...\n";

chomp(@ARGV);
unless(scalar(@ARGV) >= 1){
	print $usage;
	exit;
}

open(my $LOG, "> downsampleIlluminaReads.log");

# Calculate read numbers
my $totalReads = 0;
foreach my $file (@ARGV){
	print {$LOG} "Working on file: $file\n";
	print "Working on file: $file\n";
	$totalReads += GetFastqReadCount($file);
}

print {$LOG} "Total reads:\t$totalReads\n";
print "Total reads:\t$totalReads\n";

# Create output streams
my @ratios = ("1.0", "0.95", "0.90", "0.85", "0.80", "0.75", "0.70", "0.65", "0.60", 
	"0.55", "0.50", "0.45", "0.40", "0.35", "0.30", "0.25", "0.20", "0.15", "0.10", "0.05");
my @fhs; my %counts;
foreach my $r (@ratios){
	push(@fhs, FileHandle->new("| mash sketch -k 21 -s 10000 -r -m 2 -o downsampleSketch_$r -"));
}

# Start processing the files
foreach my $file (@ARGV){
	print {$LOG} "Sketching from file: $file\n";
	print "Sketching from file: $file\n";
	open(my $IN, "gunzip -c $file |");
	while(my $read = <$IN>){
		$read .= <$IN>;
		$read .= <$IN>;
		$read .= <$IN>;
		# Random counter out here to reduce attrition from sequential tests
		my $randm = rand();
		recursiveSelection($randm, -1, \@ratios, \@fhs, $read);
		
		for(my $x = 0; $x < scalar(@ratios); $x++){
			if($randm < $ratios[$x]){
				# Read counter update
				$counts{$ratios[$x]} += 1;
			}
		}
		$read = ""; 
	}
	close($IN);
}

# Log the read counters
foreach my $keys (sort {$b <=> $a} keys(%counts)){
	my $trueRatio = $counts{$keys} / $totalReads;
	print {$LOG} "$keys\t$counts{$keys}\t$trueRatio\n";
	print "$keys\t$counts{$keys}\t$trueRatio\n";
}
close $LOG;

exit;

# Returns the smallest ratio that this read passes
sub recursiveSelection{
	my ($randm, $startIdx, $ratios, $fhs, $read) = @_;
	if($startIdx + 1 >= scalar(@{$ratios})){
		return $randm; # We reached the end of the ratio list!
	}
	my $r = $ratios->[$startIdx + 1];
	
	if($randm < $r){
		# Passed the test
		print {$fhs->[$startIdx + 1]} $read;
		return recursiveSelection($randm, $startIdx + 1, $ratios, $fhs, $read);
	}
	return $randm;
}

sub GetFastqReadCount{
	my ($file) = @_;
	open(my $IN, "gunzip -c $file |");
	my $counter = 0;
	while(my $head = <$IN>){
		my $seq = <$IN>;
		my $plus = <$IN>;
		my $qual = <$IN>;
		$counter++;
	}
	seek($IN, 0, 0);
	close $IN;
	return $counter;
}