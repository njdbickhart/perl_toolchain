#!/usr/bin/perl
# This is a script designed to process a tab file with sequence data fastqs and produce slurm scripts

use strict;
use File::Basename;
use Getopt::Std;
use slurmTools;

my %opts;
my @modules = ("samtools/1.3-20-gd49c73b", "bwa/0.7.13-r1126");
my $usage = "perl $0 -b <base outfolder name> -t <input tab sequence files> -f <input reference fasta file>\n";
getopt('btf', \%opts);

unless(defined($opts{'b'}) && defined($opts{'o'}) && defined($opts{'t'})){
	print $usage;
	exit;
}

my %slurmWorkers; # each sample gets its own worker
mkdir $opts{'b'} || print "$!\n";
my $scriptCounter = 0;

open(my $IN, "< opts{t}") || die "Could not open input tab file!\n";
while(my $line = <$IN>){
	chomp $line; 
	my @segs = split(/\t/, $line);
	
	if(! exists($slurmWorkers{$segs[-1]})){
		$slurmWorkers{$segs[-1]} = slurmTools->new('workDir' => "$opts{b}/$segs[-1]", 
			'scriptDir' => "$opts{b}/$segs[-1]/scripts", 
			'outDir' => "$opts{b}/$segs[-1]/outLog", 
			'errDir' => "$opts{b}/$segs[-1]/errLog",
			'modules' => \@modules,
			'nodes' => 1,
			'tasks' => 5,
			'mem' => 7000);
	}
	
	my $bname = basename($segs[0]);
	my @bsegs = split(/\./, $bname);
	my $uname = $bsegs[0];
	
	my $cmd = "bwa mem -t 4 $opts{f} $segs[0] $segs[1] | samtools view -S -b - | samtools sort -m 2G -o $opts{b}/$segs[-1]/$uname.sorted.bam -T $opts{b}/$segs[-1]/$uname -";
	
	$slurmWorkers{$segs[-1]}->createGenericCmd($cmd);
	$scriptCounter++;
}
close $IN;
my $numSamples = scalar(keys(%slurmWorkers));

print "Generated $scriptCounter alignment scripts for $numSamples samples!\n";

exit;