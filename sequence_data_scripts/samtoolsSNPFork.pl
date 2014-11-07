#!/usr/bin/perl
# This script uses the Forks super package to split up Samtools chromosome processing
# samtools segments are in UCSC format and are newline separated values: (chr):(start)-(end)
# Samtools must be in the system $PATH variable in order to run this program

use strict;
use Forks::Super;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -r <reference genome fasta> -i <input bams, separated by commas> -o <output bam basename> -n <Max number of threads> -s <segments [overrides default chr processing]>\n";

getopt('rions', \%opts);

if($opts{'h'} || (!defined($opts{'r'}) && !defined($opts{'i'}) && !defined($opts{'o'}) && !defined($opts{'n'}))){
	print $usage;
	exit;
}

unless(defined($opts{'n'}) && defined($opts{'r'}) && defined($opts{'o'}) && (defined($opts{'i'}) || defined($opts{'s'}))){
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
	my $bcf = "$opts{o}.samtools.$c.bcf";
	
	push(@bcfs, $bcf);
	
	fork { sub => \&runSamtoolsBCFCaller, args => [$ref, $bcf, $opts{'i'}, $c], max_proc => $opts{'n'} };
}

waitall();

# merge files and call variants
concatBCFtoVCF(\@bcfs, "$opts{o}.samtools.merged.bcf", "$opts{o}.samtools.filtered.vcf");

# Check to see if the vcf exists and is not empty, then delete bcfs
if( -s "$opts{o}.samtools.filtered.vcf"){
	for(my $x = 0; $x < scalar(@bcfs); $x++){
		system("rm $bcfs[$x]");
	}
}



exit;

sub concatBCFtoVCF{
	my ($bcfs, $mergebcf, $vcf) = @_;
	
	my $bcfstr = join(' ', @bcfs);
	system("bcftools cat $bcfstr > $mergebcf");
	system("bcftools view $mergebcf | vcfutils.pl varFilter -D100 > $vcf");
}

sub runSamtoolsBCFCaller{
	my ($ref, $bcf, $bamstr, $chr) = @_;
	my @bams = split(/,/, $bamstr);
	my $whitebams = join(' ', @bams);
	
	system("samtools mpileup -uf $ref -r $chr $whitebams | bcftools view -bvcg - > $bcf"); 
	
	print "Finished creating bcf on chr: $chr\n";
}