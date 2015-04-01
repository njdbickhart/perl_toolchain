#!/usr/bin/perl
# This is a combination of four scripts: bamMergeUtility.pl, samtoolsSNPFork.pl, getBamStats.pl and 
# bwaFQFork.pl. It also adds in a SNPEff automation for the vcf. It generates a large merged, indexed
# bam, and an annotated VCF of SNPs and INDELs.
# Spreadsheet columns:
# 	fq1\tfq2\tlibraryname\tanimalname

use strict;
use Getopt::Long;
use File::Basename;
use bamTools;
use Forks::Super;
use simpleConfigParser;
use simpleLogger;

my ($configfile, $javaexe, $picardfolder, $outputfolder, $snpeffjar, $spreadsheet, $refgenome, $threads);
my ($runSNPFork, $runSNPEff);
my @requiredConfig = ("java", "picard", "snpeff", "runSNP", "runEFF");

my $usage = "perl $0 --fastqs <spreadsheet to be processed> --output <base output folder>
Arguments:
	--fastqs	A tab-delimited spreadsheet with fastq file locations and read group info
	--output	The base output folder for all bams and vcfs
	--reference	The reference genome fasta file

Optional:
	--config	Override default configuration (searches for ~/.mergedpipeline.cnfg by default)
	--threads	Number of threads to fork (default: 1)
	";

$configfile = "~/.mergedpipeline.cnfg";
$threads = 1;
	
GetOptions("config=s" => \$configfile,
	"fastqs=s" => \$spreadsheet,
	"output=s" => \$outputfolder,
	"reference=s" => \$refgenome, 
	"threads=i" => \$threads)
	|| die "Could not load command line options!\n$usage";
	
if(!defined($spreadsheet) || !defined($outputfolder)){
	print $usage;
	exit;
}

# Parse configuration
my $cfg = simpleConfigParser->new();
$cfg->loadConfigFile($configfile);
$cfg->checkReqKeys(\@requiredConfig);

# Population config settings
$picardfolder = $cfg->get("picard");
$javaexe = $cfg->get("java");
$snpeffjar = $cfg->get("snpeff");
$runSNPEff = $cfg->get("runEFF");
$runSNPFork = $cfg->get("runSNP");

foreach my $prog ("bwa", "samtools", "vcfutils.pl", "bcftools"){
	StaticUtils::checkReqs($prog);
}

if(! -d $outputfolder){
	mkdir($outputfolder) || die "Could not create base output folder: $outputfolder!\n";
}

# Generate log file
my $log = simpleLogger->new('logFileBaseStr' => 'MergedBamPipeline');
$log->OpenLogger($outputfolder);
$log->Info("Start", "Began pipeline run for spreadsheet: $spreadsheet in $outputfolder using $threads threads");

# Check the reference genome and see if it has been indexed yet by samtools. If not, do that indexing
my $reffai = "$refgenome.fai";
unless(-e $reffai){
	$log->Warn("RefIndexer", "Samtools fasta index not found! Generating it now.");
	system("samtools index $refgenome");
}

# Load spreadsheet for processing
open(IN, "< $spreadsheet") || die "could not open fastq spreadsheet: $spreadsheet!\n";
my %counter;
my %bams; # {sample} -> [] -> sortedbam
while(my $line = <IN>){
	chomp $line;
	# Split the line into an array based on tab delimiters
	my @segs = split(/\t/, $line);
	# This counter serves to differentiate the read group IDs on the different aligned bam files
	$counter{$segs[-1]} += 1;
	
	# The command that creates the new alignment process
	fork { sub => \&runBWAAligner, args => [$segs[0], $segs[1], $refgenome, $reffai, $outputfolder, $segs[-1], $segs[-2], $counter{$segs[-1]}, $log], max_proc => $threads };
	push(@{$bams{$segs[-1]}}, "$outputfolder/$segs[-1]/$segs[-1].$counter{$segs[-1]}.bam");
}
close IN;
waitall;

# After generating all of the bams, now do the merger

exit;

# The workhorse subroutine
sub runBWAAligner{
	# Perl needs to read in the arguments to the subroutine
	my ($fq1, $fq2, $ref, $reffai, $outdir, $base, $lib, $num, $log) = @_;
	# Take care of the names of the input and output files
	my $fqbase = basename($fq1);
	my ($fqStr) = $fqbase =~ /(.+)\.(fastq|fq).*/;
	my $bwasam = "$outdir/$base/$base.$num.sam";
	my $bwabam = "$outdir/$base/$base.$num.bam";
	#my $bwasort = "$outdir/$base/$base.$num.sorted";
	
	# Create a sub directory if needed
	mkdir("$outdir/$base") || print "";
	
	# Run the BWA mem command and create a bam file from the sam file
	print "bwa mem -R '\@RG\\tID:$base.$num\\tLB:$lib\\tPL:ILLUMINA\\tSM:$base\' -v 1 -M $ref $fq1 $fq2 > $bwasam\n";
	system("bwa mem -R '\@RG\\tID:$base.$num\\tLB:$lib\\tPL:ILLUMINA\\tSM:$base\' -v 1 -M $ref $fq1 $fq2 > $bwasam");
	my $samreader = SamFileReader->new('log' => $log);
	$samreader->prepSam($bwasam);
	#system("samtools view -b -o $bwabam -S $bwasam");
	
	# Sort the bam file
	#print "samtools sort $bwabam $bwasort\n";
	#system("samtools sort $bwabam $bwasort");
	#system("samtools index $bwasort.bam");
	
	# If the bam file exists and is not empty, then remove the sam file
	if( -s $bwabam){
		system("rm $bwasam");
		$log->Info("BWAALIGNER", "Completed sorted bam. Removing SAM file");
	}
	
	# $bwabam should be the sorted, indexed bam at the end of the process due to my API
}