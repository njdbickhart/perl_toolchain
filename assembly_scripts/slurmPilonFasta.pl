#!/usr/bin/perl
# This is a script designed to generate a series of scripts for queued pilon processing

use slurmTools;
use Getopt::Std;
use strict;
use File::Basename;
use Cwd;

my %opts;
my @modules = ("samtools", "bwa", "pilon");
my $usage = "perl $0 -f <input frag bam> -g <input genome fasta> -o <output directory> -p <partition>\n";
my $currentDir = cwd();

getopt('fgo', \%opts);

unless(defined($opts{'f'}) && defined($opts{'g'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

# Check fasta index file
my %chrlengths;
if(! -s "$opts{g}.fai"){
	print STDERR "Could not find $opts{g}.fai! Generating now...\n";
	system("module load samtools; samtools faidx $opts{g}");
}

open(my $IN, "< $opts{g}.fai") || die "Could not find $opts{g}.fai!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	$chrlengths{$segs[0]} = $segs[1];
}
close $IN;

# Generate scripts
mkdir("$opts{o}") || print "$!\n";
my $slurmWorker = slurmTools->new('workDir' => "$currentDir", 
			'scriptDir' => "$currentDir/$opts{o}/scripts", 
			'outDir' => "$currentDir/$opts{o}/outLog", 
			'errDir' => "$currentDir/$opts{o}/errLog",
			'modules' => \@modules,
			'useTime' => 0,
			'nodes' => 1,
			'tasks' => 2,
			'mem' => 25000);
if(defined($opts{'p'})){
	$slurmWorker->partition("$opts{p}");
}
foreach my $chr (keys %chrlengths){
	my $mem = int(int($chrlengths{$chr} / 1000000) * 1.15);
	$slurmWorker->mem($mem);
	
	my $cmd = "pilon --genome $opts{g} --frags $opts{f} --output $chr.pilon --outdir $opts{o} --fix bases --targets $chr --verbose --nostrays";
	$slurmWorker->createGenericCmd($cmd, "pilon$chr");
}

$slurmWorker->queueJobs;
print "Queued slurm jobs!\n";

exit;