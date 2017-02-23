#!/usr/bin/perl
# This is a script designed to process a tab file with sequence data fastqs and produce slurm scripts

use strict;
use File::Basename;
use Getopt::Std;
use slurmTools;
use Cwd;

my %opts;
my @modules = ("samtools/1.3-20-gd49c73b", "bwa/0.7.13-r1126");
my $usage = "perl $0 -b <base outfolder name> -t <input tab sequence files> -f <input reference fasta file>\n";
getopt('btf', \%opts);

unless(defined($opts{'b'}) && defined($opts{'f'}) && defined($opts{'t'})){
	print $usage;
	exit;
}

my %slurmWorkers; # each sample gets its own worker
mkdir $opts{'b'} || print "$!\n";
my $scriptCounter = 0;
my $currentDir = cwd();
my $fasta = $opts{'f'};

if( -e "$currentDir/$opts{f}"){
	$fasta = "$currentDir/$opts{f}";
}else{
	print STDERR "Error locating fasta file in current directory! Please check file path!\n$usage\n";
}

open(my $IN, "< $opts{t}") || die "Could not open input tab file!\n";
while(my $line = <$IN>){
	chomp $line; 
	my @segs = split(/\t/, $line);
	
	if(! exists($slurmWorkers{$segs[-1]})){
		$slurmWorkers{$segs[-1]} = slurmTools->new('workDir' => "$currentDir/$opts{b}/$segs[-1]", 
			'scriptDir' => "$currentDir/$opts{b}/$segs[-1]/scripts", 
			'outDir' => "$currentDir/$opts{b}/$segs[-1]/outLog", 
			'errDir' => "$currentDir/$opts{b}/$segs[-1]/errLog",
			'modules' => \@modules,
			'useTime' => 0,
			'nodes' => 1,
			'tasks' => 7,
			'mem' => 9000);
	}
	
	my $bname = basename($segs[0]);
	my @bsegs = split(/\./, $bname);
	my $uname = $bsegs[0];
	
	my $cmd = "bwa mem -t 5 $fasta $segs[0] $segs[1] | samtools view -S -b - | samtools sort -m 2G -o $uname.sorted.bam -T $uname -";
	
	$slurmWorkers{$segs[-1]}->createGenericCmd($cmd);
	$scriptCounter++;
}
close $IN;
my $numSamples = scalar(keys(%slurmWorkers));

print "Generated $scriptCounter alignment scripts for $numSamples samples!\n";

exit;
