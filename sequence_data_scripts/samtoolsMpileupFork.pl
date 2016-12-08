#!/usr/bin/perl
# This script generates mpileup commands for samtools for use in downstream analysis
# Requires htslib samtools

use strict;
use Forks::Super;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -r <reference genome fasta> -i <input bams, separated by commas> -o <output bam basename> -n <Max number of threads> -s <segments [overrides default chr processing]> -t <generate only slurm scripts>\n";

getopt('rions', \%opts);

if($opts{'h'} || (!defined($opts{'r'}) && !defined($opts{'i'}) && !defined($opts{'o'}) && !defined($opts{'n'}))){
	print $usage;
	exit;
}

unless(defined($opts{'n'}) && defined($opts{'r'}) && defined($opts{'o'}) && (defined($opts{'i'}))){
	print "Missing mandatory arguments!\n";
	print $usage;
	exit;
}

my $ref = $opts{'r'};
my $reffai = "$opts{r}.fai";
unless(-e $reffai){
	print "Generating reference fasta index...\n";
	system("samtools index $ref");
}


my @chrs;
if(defined($opts{'s'})){
	open(IN, "< $opts{s}") || die "Could not open segment file: $opts{s}!\n";
	while(my $line = <IN>){
		chomp $line;
		push(@chrs, $line);
	}
	close IN;
}else{
	open(IN, "< $reffai") || die "Could not open $reffai reference fasta index!\n";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		push(@chrs, $segs[0]);
	}
	close IN;
}

my @bcfs;
foreach my $c (@chrs){	
	my $bcf = "$opts{o}.mpileup.$c.bcf";
	
	push(@bcfs, $bcf);
	
	fork { sub => \&runSamtoolsBCFCaller, args => [$ref, $bcf, $opts{'i'}, $c], max_proc => $opts{'n'} };
}

waitall();

exit;


sub runSamtoolsBCFCaller{
	my ($ref, $bcf, $bamstr, $chr) = @_;
	my @bams = split(/,/, $bamstr);
	my $whitebams = join(' ', @bams);
	
	#samtools mpileup -go genome.canada.raw.bcf -f ../reference/umd_3_1_reference_1000_bull_genomes.fa
	system("samtools mpileup -go $bcf -f $ref -r $chr $whitebams"); 
	
	print "Finished creating bcf on chr: $chr\n";
}
