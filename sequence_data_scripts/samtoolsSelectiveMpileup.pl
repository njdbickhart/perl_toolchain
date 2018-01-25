#!/usr/bin/perl
# This is a script designed to selectively run samtools mpileup on a select region and automate output generation
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5000
#SBATCH --nodes=1
#SBATCH --job-name=SelMpileup
#SBATCH --output=mpileup_%j.out

use strict;
use Getopt::Std;

my $usage = "perl $0 -b <bam file list> -s <samtools segment> -n <output base name> -f <reference genome name>\n";

my %opts;
getopt('bsnf', \%opts);

unless(defined($opts{'b'}) && defined($opts{'s'}) && defined($opts{'n'}) && defined($opts{'f'})){
	print $usage;
	exit;
}

print "module load samtools; samtools mpileup -go $opts{n}.mpileup.bcf -f $opts{f} -r $opts{s} -b $opts{b}\n";
system("module load samtools; samtools mpileup -go $opts{n}.mpileup.bcf -f $opts{f} -r $opts{s} -b $opts{b}");

exit;
