#!/usr/bin/perl
# This is a script designed to process a list of bam files, subsection genomic regions and call SNPs and INDELS
# using the samtools mpileup pipeline
# Adapted from the original script samtoolsMpileupBamsSlurm.pl

use strict;
use File::Basename;
use Sub::Identify;
use Getopt::Std;
use slurmTools;
use Cwd;

my %opts;
my $processChunks = 10000000; # 10 megabase chunks for variant calling
my @modules = ("samtools/1.9", "bcftools/1.9");
my $usage = "perl $0 -b <base outfolder name> -s <OPTIONAL: ref sections for variant calling in UCSC format> -t <input newline separted bam files> -f <input reference fasta file> -m <boolean: generate and queue merger scripts> -p <optional: partition name>\n";
getopt('btfp', \%opts);

unless(defined($opts{'b'}) && defined($opts{'f'}) && defined($opts{'t'})){
	print $usage;
	exit;
}

my %bcfWorkers; # each section gets its own worker
my %vcfWorkers; # each section gets its own worker
my @slurmBcfs; # each section bcf gets added to the list
my @slurmVcfs; # each section vcf gets added to the list
mkdir $opts{'b'} || print "$!\n";
my $scriptCounter = 0;
my $currentDir = cwd();
my $fasta = $opts{'f'};

if( -e "$currentDir/$opts{f}"){
	$fasta = "$currentDir/$opts{f}";
}else{
	print STDERR "Error locating fasta file in current directory! Please check file path!\n$usage\n";
	exit;
}

if(! -e "$currentDir/$opts{f}.fai"){
	print STDERR "Did not find fasta index file! Generating it...\n";
	system("module load samtools; samtools faidx $currentDir/$opts{f}");
	if(! -e "$currentDir/$opts{f}.fai"){
		print STDERR "Error generating fasta index!\n$usage\n";
		exit;
	}
}

# Check to see if there are SAM sections, otherwise generate them
my @sections;
if(defined($opts{'s'})){
	open(my $IN, "< $opts{s}") || die "Could not open sam file segments!\n";
	while(my $line = <$IN>){
		chomp $line;
		push(@sections, $line);
	}
	close $IN;
}else{
	print STDERR "Generating sam file segments using a chunk size of $processChunks bp!\n";
	# Read fasta index and get chr sizes
	open(my $IN, "< $opts{f}.fai");
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		for(my $x = 1; $x <= $segs[1]; $x += $processChunks){
			my $ne = $x + $processChunks - 1;
			if($ne >= $segs[1]){
				$ne = $segs[1];
			}
			push(@sections, "$segs[0]:$x\-$ne");
		}
	}
	close $IN;
}

# Generate bcf generation scripts
my %bcfJobids; # queued bcf jobs by called section ID
my @vcfJobids; # queued vcf jobs
foreach my $section (@sections){
	
	if(! exists($bcfWorkers{$section})){
		$bcfWorkers{$section} = slurmTools->new('workDir' => "$currentDir/$opts{b}/bcf_files", 
			'scriptDir' => "$currentDir/$opts{b}/bcf_files/scripts", 
			'outDir' => "$currentDir/$opts{b}/bcf_files/outLog", 
			'errDir' => "$currentDir/$opts{b}/bcf_files/errLog",
			'modules' => \@modules,
			'useTime' => 0,
			'nodes' => 1,
			'tasks' => 4,
			'mem' => 25000);
	}

	my $uname = "bcf_segment_$section";
	$uname =~ s/[\:-]/_/g;
	
	my $cmd = "bcftools mpileup -Ob -r $section -o $uname.bcf -f $currentDir/$opts{f} --threads 3 -b $currentDir/$opts{t}";
	push(@slurmBcfs, "$uname.bcf");
	
	my $sname = "mpileup_$section";
	$sname =~ s/[\:-]/_/g;
	$bcfWorkers{$section}->createGenericCmd($cmd, $sname);
	$scriptCounter++;
	$bcfWorkers{$section}->queueJobs;
	
	# Now queue up the variant callers
	$vcfWorkers{$section} = slurmTools->new('workDir' => "$currentDir/$opts{b}/vcf_files",
		'scriptDir' => "$currentDir/$opts{b}/vcf_files/scripts",
		'outDir' => "$currentDir/$opts{b}/vcf_files/outLog",
		'errDir' => "$currentDir/$opts{b}/vcf_files/errLog",
		'modules' => \@modules,
		'useTime' => 0,
		'dependencies' => $bcfWorkers{$section}->jobIds,
		'nodes' => 1,
		'tasks' => 4,
		'mem' => 25000);
	
	my $vname = "vcf_segment_$section";
	$vname =~ s/[\:-]/_/g;
	my $vcmd = "bcftools call -vmO z -r $section -o $vname.vcf.gz --threads 3 $uname.bcf";
	
	$vcfWorkers{$section}->createGenericCmd($cmd, "call_$section");
	$vcfWorkers{$section}->queueJobs;
	
	push(@vcfJobids, @{$vcfWorkers{$section}->jobIds});
	push(@slurmVcfs, "$vname.vcf.gz");
}

my $numSamples = scalar(keys(%bcfWorkers));

print "Generated $scriptCounter mpileup and vcf scripts for $numSamples samples!\n";

# Concatenate VCF files into merged file
my $concatWorker = slurmTools->new('workDir' => "$currentDir/$opts{b}/vcf_files", 
			'scriptDir' => "$currentDir/$opts{b}/vcf_files/scripts", 
			'outDir' => "$currentDir/$opts{b}/vcf_files/outLog", 
			'errDir' => "$currentDir/$opts{b}/vcf_files/errLog",
			'modules' => \@modules,
			'dependencies' => \@vcfJobids,
			'useTime' => 0,
			'nodes' => 1,
			'tasks' => 3,
			'mem' => 10000);

# Generate list of vcf files
open(my $OUT, "> $currentDir/$opts{b}/vcf_files/vcf_files.list");
foreach my $s (@slurmVcfs){
	print {$OUT} "$s\n";
}
close $OUT;

my $vcmd = "bcftools concat -a -d all -f vcf_files.list -O z -o merged_all_section.vcf.gz --threads 3";
$concatWorker->createGenericCmd($vcmd, "concat");
$concatWorker->queueJobs;


exit;

sub urlHash{
        my @alphabet = (('a'..'z'), 0..9);
        my %collection;
        for (1..30) {
                my $len = rand(11) + 10;
                my $key = join '', map {$alphabet[rand(@alphabet)]} 1..$len;
                $collection{$key} ? redo : $collection{$key}++;
        }
        my @sets = keys(%collection);
        return $sets[0];
}
