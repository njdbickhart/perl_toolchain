#!/usr/bin/perl
# This is a script designed to generate a listing of smaller scripts for routine sequence analysis
# It can "daisychain" multiple forms of analysis recursively
# Simplicity at first, but I can add on later if needed.

use strict;
use Getopt::Std;
use Switch;

my %opts;
my $usage  = "perl $0 -p <type, comma separated> -j <jobid basename> -i <input files, comma separated> -o <output files, comma separated> [slurm options]\n
Types:
	convertBam		Change chr names to capitals and add read groups
	runBCFtools		Process bam file into bcf for downstream calling
	
Slurm commands:
	(n) nodes		Number of nodes to use
	(t) taskspernode	Tasks per node
	(m) mem			Memory to use
	(w) workdir		Working directory
	(h) time		Time to run in minutes
	(s) stdout		Std output directory
	(e) stderr		Std err directory
	(a) modules		Modules to load [comma separated; applies to all scripts]
\n";

getopt('pjiontmwhsea', \%opts);

unless(defined($opts{'p'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

exit;
sub convertBam{
	my ($inputfile, $outputfile) = @_;
	
	my $str = "perl -e 'chomp @ARGV;
		my @basefile = split(/\/, $ARGV[0]);
		my @basename = split(/[\-\.]/, $basefile[-1]);
		
		system(qq{samtools view -h $inputfile | sed \"s/chr/Chr/g\" | samtools view -bS - | samtools addreplacerg -r \"@RG\tID:$basename[0]\tPL:Illumina\tLB:$basename[0]\tSM:$basename[0]\" -o $ARGV[1]});
		' $inputfile $outputfile";
	return $str;
}

sub runBCFtools{
	my ($inputfile, $outputfile) = @_;
	
	# Not yet implemented
}

sub recursiveJobAdd{
	my ($jobstr, $jobarray, $inputfile, $outputarray, $jobsubid) = @_;
	
	my @subjobs = @{$jobarray}; 
	my $outputfile = $outputarray->[$jobsubid];
	my $newoutput; my $tjstr;
	switch($subjobs[$jobsubid]){
		case "convertBam" { $tjstr = convertBam($inputfile, $outputfile)}
		case "runBCFtools" { $tjstr = runBCFtools($inputfile, $outputfile)}
		else {$tjstr = "#empty\n"}
	}
	
	$jobstr .= "$tjstr";
	
	if($jobsubid + 1 >= scalar(@subjobs)){
		return $jobstr;
	}else{
		return recursiveJobAdd($jobstr, $jobarray, $outputfile, $outputarray, $jobsubid + 1);
	}
}

sub getSlurmHeader{
	my ($nodes, $jidbase, $taskspernode, $mem, $workdir, $time, $stdout, $stderr, $modules, $jobid) = @_;
	
	my $batchStr = "#SBATCH --nodes=$nodes\n";
	$batchStr .= "#SBATCH --ntasks-per-node=$taskspernode\n";
	$batchStr .= "#SBATCH --mem=$mem\n";
	unless( -d $workdir){
		mkdir($workdir);
	}
	$batchStr .= "#SBATCH --workdir=$workdir\n";
	$batchStr .= "#SBATCH --time=$time\n";
	unless( -d $stdout){
		mkdir($stdout);
	}
	$batchStr .= "#SBATCH --output=$stdout\/$jidbase.$jobid.out\n";
	unless( -d $stderr){
		mkdir($stderr);
	}
	$batchStr .= "#SBATCH --error=$stderr\/$jidbase.$jobid.err\n\n";
	
	my @mods = split(/,/, $modules);
	
	foreach my $m (@mods){
		$batchStr .= "module load $m\n";
	}
	
	$batchStr .= "cd $workdir\n";
	
	return $batchStr;
}
	
	

