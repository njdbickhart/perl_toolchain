#!/usr/bin/perl
# Base algorithm adapted from: http://genomics-array.blogspot.com/2011/02/calculating-n50-of-contig-assembly-file.html
# Edit 9/21/2016 -- Added option to calculate NG50 instead of N50 and option to read fai files

use strict;
use Getopt::Std;

my %opts;

my $usage = "perl $0 -f <input fasta file> OR -i <input fai file> 
	-f	Input fasta file
	OR
	-i	Input samtools fai file
	
	[optional]
	-g	Genome size (for NG50 calculation)\n";
getopt('fig', \%opts);

if(!defined($opts{'f'}) && !defined($opts{'i'})){
	print $usage;
	exit;
}

my $IN;
my $fasta = 0;

if(defined($opts{'f'})){
	open($IN, "< $opts{f}") || die "Could not open input fasta file: $opts{f}!\n";
	$fasta = 1;
}elsif(defined($opts{'i'})){
	open($IN, "< $opts{i}") || die "Could not open input fai file: $opts{i}!\n";
}

my $ng50 = 0;
my $totalLength = 0; 
my @arr;
if(defined($opts{'g'})){
	$ng50 = 1;
	$totalLength = $opts{'g'};
}

## Read Fasta File and compute length ###
if(defined($opts{'f'})){
	my $length = 0;	
	while(my $line = <$IN>){
		chomp $line; 
		if($line =~ /^>/){
			push (@arr, $length);
			unless($ng50){
				$totalLength += $length;
			}
			$length = 0;
			next;
		}
		$length += length($line);
	}
}elsif(defined($opts{'i'})){
	# Read FAI file and compute length
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		push(@arr, $segs[1]);
		unless($ng50){
			$totalLength += $segs[1];
		}
	}
}

close($IN);

my @sort = sort {$b <=> $a} @arr;
my $n50; 
my $L50;
foreach my $val(@sort){
	$n50+=$val;
	$L50++;
	if($n50 >= $totalLength/2){
		my $char = ($ng50)? "g" : "";
		print "N$char\50 length:\t$n50\n";
		print "N$char\50 value:\t$val\n";
		print "L$char\50 value:\t$L50\n";
		last; 
	}
}

exit;
