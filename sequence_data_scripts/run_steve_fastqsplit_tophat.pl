#!/usr/bin/perl
# This script runs through a directory that Steve made after demultiplexing sequencing libraries
# Output is placed in the sample name output folder in the base output directory specified
# It assumes that the libraries are paired end

use strict;
use Getopt::Std;

my %opts;
getopt('botslrd', \%opts);
my $usage = "$0 -b <base directory (project name)> -o <output base folder> -r <reference genome base> -t <number of threads for tophat> -l <log file> -s [comma separated sample skip list] -d [dry run flag]\n";
my $threads = 1;

unless(defined($opts{'b'}) && defined($opts{'o'}) && defined($opts{'l'}) && defined($opts{'r'})){
	print $usage;
	exit;
}

if(defined($opts{'t'})){
	$threads = $opts{'t'};
}

my %skip;
if(defined($opts{'s'})){
	my @segs = split(/,/, $opts{'s'});
	foreach my $s (@segs){
		$skip{$s} = 1;
	}
}
mkdir("$opts{o}") || print "MKDIR: $!\n";
open(LOG, "> $opts{l}");

print LOG "Basedir: $opts{b}\nOutputdir: $opts{o}\nThreads: $opts{t}\n";
opendir(BASE, $opts{'b'}) || die "Could not open base folder location!\n$usage";
while(my $file = readdir(BASE)){
	if($file =~ /^\./ || !(-d "$opts{b}/$file")){next;}
	if(exists($skip{$file})){next;}
	my ($sampbase) = $file =~ /Sample_(.+)/;
	mkdir("$opts{o}/$sampbase") || print "Samp MKDIR: $!\n";
	
	opendir(SAMP, "$opts{b}/$file") || die "Could not open sample folder location!\n$usage";
	
	# Steve's file format: H5_AGTCAA_L008_R1_001.fastq.gz
	my %fastqs; # {last segment}->{read segment} = file
	while(my $f = readdir(SAMP)){
		if($f =~ /^\./ || !($f =~ /.*\.fastq\.gz/)){next;}
		my @fsegs = split(/[\.\_]/, $f);
		$fastqs{$fsegs[4]}->{$fsegs[3]} = $f;
	}
	
	# Generate the fastq1 and fastq2 strings
	my ($fastq1, $fastq2);
	foreach my $smt (keys(%fastqs)){
		$fastq1 .= "$opts{b}/$file/" . $fastqs{$smt}->{'R1'} . ",";
		$fastq2 .= "$opts{b}/$file/" . $fastqs{$smt}->{'R2'} . ",";
	}
	# get rid of last comma
	chop($fastq1);
	chop($fastq2);
	
	# Now to run tophat
	my $time = timestamp();
	print LOG "$time -- tophat2 -p $threads -o $opts{o}/$sampbase/ $opts{r} $fastq1 $fastq2\n";
	print "$time -- tophat2 -p $threads -o $opts{o}/$sampbase/ $opts{r} $fastq1 $fastq2\n";
	
	unless($opts{'d'}){
		system("tophat2 -p $threads -o $opts{o}/$sampbase/ $opts{r} $fastq1 $fastq2");
	}
	
	close SAMP;
}
close BASE;
close LOG;
exit;


sub timestamp{
	my @timearray = localtime(time);
	return "$timearray[4]/$timearray[3]/$timearray[5] - $timearray[2]:$timearray[1]:$timearray[0] ";
}