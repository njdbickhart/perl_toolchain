#!/usr/bin/perl
# This script uses the Forks super package to split up bwa runs
# FQ spreadsheets are tab delimited entries where the path to first file and the second file are separated by tabs
# The last column of the spreadsheet is used for the sample output directory

use strict;
use File::Basename;
use Forks::Super;
use Getopt::Std;

# Define the program inputs and a help menu (good practice)
my %opts;
my $usage = "perl $0 [options -- all mandatory]
	-r <reference genome fasta; File must be indexed by the same version of BWA on your path> 
	-i <input fq spreadsheet; format(per line): fastq1 (tab) fastq2 (tab) sample name> 
	-o <base output dir; will create if not present> 
	-n <Max number of threads>\n";

# Retrieve the arguments from the flags
getopt('rions', \%opts);

# If the user requests "help", then print the usage and exit
if($opts{'h'}){
	print $usage;
	exit;
}

# If the user doesn't enter any flags, or fails to enter one, print the usage and exit
unless(defined($opts{'n'}) && defined($opts{'r'}) && (defined($opts{'i'})) && (defined($opts{'o'}))){
	print "Missing mandatory arguments!\n";
	print $usage;
	exit;
}

# Make the output directory if it doesn't already exist
mkdir("$opts{o}") || print "$!\n";

# Check the reference genome and see if it has been indexed yet by samtools. If not, do that indexing
my $ref = $opts{'r'};
my $reffai = "$opts{r}.fai";
unless(-e $reffai){
	print "Generating reference fasta index...\n";
	system("samtools index $ref");
}

# Now, open the fastq spreadsheet and start processing the files with the workhorse subroutine
open(IN, "< $opts{i}") || die "could not open fastq spreadsheet: $opts{i}!\n";
my %counter;
while(my $line = <IN>){
	chomp $line;
	# Split the line into an array based on tab delimiters
	my @segs = split(/\t/, $line);
	# This counter serves to differentiate the read group IDs on the different aligned bam files
	$counter{$segs[-1]} += 1;
	
	# The command that creates the new alignment process
	fork { sub => \&runBWAAligner, args => [$segs[0], $segs[1], $ref, $reffai, $opts{'o'}, $segs[-1], $counter{$segs[-1]}], max_proc => $opts{'n'} };
	
}

waitall;

exit;

# The workhorse subroutine
sub runBWAAligner{
	# Perl needs to read in the arguments to the subroutine
	my ($fq1, $fq2, $ref, $reffai, $outdir, $base, $num) = @_;
	# Take care of the names of the input and output files
	my $fqbase = basename($fq1);
	my ($fqStr) = $fqbase =~ /(.+)\.(fastq|fq).*/;
	my $bwasam = "$outdir/$base/$fqStr.sam";
	my $bwabam = "$outdir/$base/$fqStr.bam";
	my $bwasort = "$outdir/$base/$fqStr.sorted";
	
	# Create a sub directory if needed
	mkdir("$outdir/$base") || print "$!\n";
	
	# Run the BWA mem command and create a bam file from the sam file
	print "bwa mem -R '\@RG\\tID:$base.$num\\tLB:$base\\tPL:ILLUMINA\\tSM:$base\' -v 1 -M $ref $fq1 $fq2 > $bwasam\n";
	system("bwa mem -R '\@RG\\tID:$base.$num\\tLB:$base\\tPL:ILLUMINA\\tSM:$base\' -v 1 -M $ref $fq1 $fq2 > $bwasam");
	system("samtools view -b -o $bwabam -S $bwasam");
	
	# Sort the bam file
	print "samtools sort $bwabam $bwasort\n";
	system("samtools sort $bwabam $bwasort");
	system("samtools index $bwasort.bam");
	
	# If the bam file exists and is not empty, then remove the sam file
	if( -s $bwabam){
		system("rm $bwasam");
	}
	
	# If the sorted bam file exists and is not empty, then remove the bam file
	if( -s "$bwasort.bam"){
		system("rm $bwabam");
	}
}

sub checkReqs{
	my ($usage, @progs) = @_;
	
	my @notfound;
	foreach my $p (@progs){
		my $found = 0;
		foreach my $path (split(/:/, $ENV{PATH})) {
			if( -f "$path/$p") {
				$found = 1;
				last;
			}
		}
		if(!$found){
			push(@notfound, $p);
		}
	}
	if(scalar(@notfound) > 0){
		print "Error! Could not find the following programs on your path:\n";
		print join(', ', @notfound) . "\n";
		print $usage;
		exit;
	}
}