#!/usr/bin/perl
# This script uses the Forks super package to split up mrsfast runs

use strict;
use File::Basename;
use Forks::Super;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -r <reference genome fasta> -i <input fqs, separated by commas> -o <output dir> -n <Max number of threads>\n";

getopt('rions', \%opts);

if($opts{'h'}){
	print $usage;
	exit;
}

unless(defined($opts{'n'}) && defined($opts{'r'}) && (defined($opts{'i'})) && (defined($opts{'o'}))){
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

my @fqs = split(/,/, $opts{'i'});
foreach my $f (@fqs){
	
	fork { sub => \&runMrsfastAligner, args => [$f, $ref, $reffai, $opts{'o'}], max_proc => $opts{'n'} };
}

waitall;

exit;

sub runMrsfastAligner{
	my ($fq, $ref, $reffai, $outdir) = @_;
	my $fqbase = basename($fq);
	my ($fqStr) = $fqbase =~ /(.+)\.fq/;
	my $mrsfastsam = "$outdir/$fqStr.sam";
	my $mrsfastbam = "$outdir/$fqStr.bam";
	
	my $seqcomp = ($fqbase =~ /.+\.gz/)? "--seqcomp" : "";
	print "mrsfast --search $ref --seq $fq $seqcomp -o $mrsfastsam\n";
	system("mrsfast --search $ref --seq $fq $seqcomp -o $mrsfastsam");
	system("samtools view -bt $reffai -o $mrsfastbam -S $mrsfastsam");
	if( -s $mrsfastbam){
		system("rm $mrsfastsam");
	}
}
