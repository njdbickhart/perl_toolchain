#!/usr/bin/perl
# This script is designed to process the aligned tophat output. It runs cufflinks to generate transcript isoforms,
# then it runs cuffcompare to test cufflinks output data, then it runs cuffmerge to merge transcript isoforms for all samples,
# finally it runs cuff-diff to identify statistically significant expression among samples

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 [options]
	-i	input directory list, separated by commas
	-o	output prefix for cufflinks
	-r	reference genome fasta
	-g	reference genome transcript gtf file
	-c	output prefix for cuffcompare
	-d	output file name for cuffdiff
	-t	number of threads to use\n";
	
getopt('iorgcdt', \%opts);

unless(defined($opts{'i'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my @files = split(/,/, $opts{'i'});
print STDERR "Working on " . scalar(@files) . " directories!\n";
my @linksfiles;
my @bams;
my $cuffmerge = "merged_asm/merged.gtf";

# Process each sample
foreach my $d (@files){
	chomp $d;
	print STDERR "Working on directory: $d\n";
	assemblyAndStats($opts{'r'}, $opts{'g'}, $opts{'o'}, $d, $opts{'c'}, $opts{'t'});
	push(@linksfiles, "$d/transcripts.gtf");
	push(@bams, "$d/accepted_hits.bam");
}

#my $cuffcomparestr = join(" ", @linksfiles);
# Cuffcompare stats 
#system("cuffcompare -o $opts{c} -r $opts{g}");

open(OUT, "> cufflinks_files.txt");
foreach my $f (@linksfiles){
	print OUT "$f\n";
}
close OUT;

system("cuffmerge -g $opts{g} -p $opts{t} cufflinks_files.txt");

my $bamstr = join(" ", @bams);
system("cuffdiff -p $opts{t} $cuffmerge $bamstr");

exit;
sub assemblyAndStats{
	my ($refgenome, $refgtf, $outpre, $outdir, $outcompare, $threads) = @_;
	
	system("cufflinks -o $outdir -g $refgtf -p $threads $outdir/accepted_hits.bam");
	
}