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
use fastqcParser;
use threads;
use threads::shared;
use threadPool;

my ($configfile, $javaexe, $picardfolder, $outputfolder, $snpeffjar, $spreadsheet, $refgenome, $threads, $fastacoords, $fastqc);
my ($runSNPFork, $runSNPEff, $runFQC);
my @requiredConfig = ("java", "picard", "snpeff", "fastqc", "runSNP", "runEFF", "runFQC");

my $usage = "perl $0 --fastqs <spreadsheet to be processed> --output <base output folder> --reference <reference genome fasta>
Arguments:
	--fastqs	A tab-delimited spreadsheet with fastq file locations and read group info
	--output	The base output folder for all bams and vcfs
	--reference	The reference genome fasta file

Optional:
	--config	Override default configuration (searches for ~/.mergedpipeline.cnfg by default)
	--coords	Reference genome fasta coordinates for threaded SNP calling
	--threads	Number of threads to fork (default: 1)
	";

$configfile = "$ENV{HOME}/.mergedpipeline.cnfg";
$threads = 1;
	
GetOptions("config=s" => \$configfile,
	"fastqs=s" => \$spreadsheet,
	"output=s" => \$outputfolder,
	"reference=s" => \$refgenome,
	"coords=s" => \$fastacoords, 
	"threads=i" => \$threads);
	
if(!defined($spreadsheet) || !defined($outputfolder)){
	print $usage;
	exit;
}

# Parse configuration
my $cfg = simpleConfigParser->new();
$cfg->loadConfigFile($configfile);
$cfg->checkReqKeys(\@requiredConfig);


# Population config settings
$picardfolder = $cfg->getKey("picard");
$fastqc = $cfg->getKey("fastqc");
$javaexe = $cfg->getKey("java");
$snpeffjar = $cfg->getKey("snpeff");
$runFQC = $cfg->getKey("runFQC");
$runSNPEff = $cfg->getKey("runEFF");
$runSNPFork = $cfg->getKey("runSNP");

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

$log->Info("Start", "Using config file: $configfile");
$log->Info("Start", "Exes: picard-$picardfolder fastqc-$fastqc java-$javaexe snpeff-$snpeffjar");
$log->Info("Start", "Run Settings: runFQC-$runFQC runSNPEFF-$runSNPEff runSNPFork-$runSNPFork");

# Check the reference genome and see if it has been indexed yet by samtools. If not, do that indexing
my $reffai = "$refgenome.fai";
unless(-e $reffai){
	$log->Warn("RefIndexer", "Samtools fasta index not found! Generating it now.");
	system("samtools index $refgenome");
}

$log->Info("SpreadSheet", "Starting fastq file processing.");
# Load spreadsheet for processing
open(IN, "< $spreadsheet") || die "could not open fastq spreadsheet: $spreadsheet!\n";
my %counter;
my %bams; # {sample} -> [] -> sortedbam
my @parsers;
my $alignthreads = ($runFQC)? $threads / 3 : $threads; # If we are running fastqc at the same time, then only run one alignment thread per fastqc thread, rounded down
my $fqcThreads = ($threads / 3) * 2;

# Ensure that we don't have fractional threads!
if($alignthreads < 1){
	$alignthreads = 1;
}
if($fqcThreads < 1){
	$fqcThreads = 1;
}
my $fqcPool = threadPool->new('maxthreads' => $fqcThreads);
while(my $line = <IN>){
	chomp $line;
	# Split the line into an array based on tab delimiters
	my @segs = split(/\t/, $line);
	# This counter serves to differentiate the read group IDs on the different aligned bam files
	$counter{$segs[-1]} += 1;
	
	# If we're generating statistics on all of the fastqs, then fork fastqc
	if($runFQC){
		# file => fastq file, sample => sample_name, library => library_name, readNum => first_or_second_read, log => simpleLogger
		my $fqcparser1 = fastqcParser->new('file' => $segs[0], 'sample' => $segs[-1], 'library' => $segs[-2], 'readNum' => 1, 'log' => $log);
		my $fqcparser2 = fastqcParser->new('file' => $segs[1], 'sample' => $segs[-1], 'library' => $segs[-2], 'readNum' => 1, 'log' => $log);
		
		# Running fastqc takes the longest, so we're only going to fork this section
		#fork { sub => \&fastqcWrapper, args => [$fqcparser1, $fastqc], max_proc => $threads };
		#fork { sub => \&fastqcWrapper, args => [$fqcparser2, $fastqc], max_proc => $threads };
		#fastqcWrapper($fqcparser1, $fastqc);
		#fastqcWrapper($fqcparser2, $fastqc);
		
		my $fq1thr = threads->create(sub{$fqcparser1->runFastqc($fastqc)});
		$fqcPool->submit($fq1thr);
		my $fq2thr = threads->create(sub{$fqcparser2->runFastqc($fastqc)});
		$fqcPool->submit($fq2thr);
		push(@parsers, ($fqcparser1, $fqcparser2));
	}	
	
	# The command that creates the new alignment process
	#fork { sub => \&runBWAAligner, 
	#	args => [$segs[0], $segs[1], $refgenome, $reffai, $outputfolder, $segs[-1], $segs[-2], $counter{$segs[-1]}, $log, $javaexe, $picardfolder], 
	#	max_proc => $threads };
	runBWAAligner($segs[0], $segs[1], $refgenome, $reffai, $outputfolder, $segs[-1], $segs[-2], $counter{$segs[-1]}, $log, $javaexe, $picardfolder);
	push(@{$bams{$segs[-1]}}, "$outputfolder/$segs[-1]/$segs[-1].$counter{$segs[-1]}.nodup.bam");
}
close IN;
waitall;
$fqcPool->joinAll();
$log->Info("Spreadsheet", "All alignment threads completed");

# if we ran fastqc, then start parsing the data
if($runFQC){
	$log->Info("Fastqc", "Parsing data from " . (scalar(@parsers)) . " fastqc instances");
	open(OUT, "> $outputfolder/$spreadsheet.fastqc");
	my @headers = $parsers[0]->getHeaderArray();
	print OUT join("\t", @headers);
	print OUT "\n";
	foreach my $p (@parsers){
		$p->parseStats();
		$p->cleanUp();
		my @rows = $p->getOutArray();
		print OUT join("\t", @rows);
		print OUT "\n";
	}
	close OUT;
	$log->Info("Fastqc", "Quality metric data on fastq files located in: $outputfolder/$spreadsheet.fastqc");
}

# After generating all of the bams, now do the merger
$log->Info("Merger", "Beginning bam merger.");
my @finalbams;
my $mergerPool = threadPool->new('maxthreads' => $threads);
foreach my $sample (keys(%bams)){
	my $finalbam = "$outputfolder/$sample/$sample.merged.bam";
	if(scalar(@{$bams{$sample}}) > 1){
		#fork { sub => \&samMerge,
		#	args => [$bams{$sample}, "$outputfolder/$sample", $log, $sample],
		#	max_proc => $threads };
	
		samMerge($bams{$sample}, "$outputfolder/$sample", $log, $sample);
		my $thr = threads->create('samMerge', $bams{$sample}, "$outputfolder/$sample", $log, $sample);
		$mergerPool->submit($thr);
		push(@finalbams, $finalbam);
	}else{
		push(@finalbams, $bams{$sample}->[0]);
	}
}

waitall;
$mergerPool->joinAll();
$log->Info("Merger", "Finished bam merger.");

my $snpPool = threadPool->new('maxthreads' => $threads);
if($runSNPFork){
	$log->Info("SNPFORK", "Beginning Samtools SNP calling methods");
	mkdir("$outputfolder/vcfs") || $log->Warn("SNPFORK", "Could not create vcf directory!");
	my @coords;
	if(defined($fastacoords)){
		open(IN, "< $fastacoords") || die $log->Fatal("SNPFORK", "Could not open fasta coordinate file! $fastacoords");
		while(my $line = <IN>){
			chomp $line;
			push(@coords, $line);
		}
		close IN;
	}else{
		open(IN, "< $reffai") || die $log->Fatal("SNPFORK", "Could not open reference fai file! $reffai");
		while(my $line = <IN>){
			chomp $line;
			my @segs = split(/\t/, $line);
			push(@coords, "$segs[0]:1-$segs[1]");
		}
		close IN;
	}
	
	my $samexe = SamtoolsExecutable->new('log' => $log);
	my @vcfsegs;
	foreach my $c (@coords){
		my $out = "$outputfolder/vcfs/combined.$c.vcf";
		#fork { sub => \&SamtoolsExecutable::GenerateSamtoolsVCF,
		#	args => [$samexe, \@finalbams, $out, $refgenome, $c],
		#	max_proc => $threads };
		my $thr = threads->create(sub{$samexe->GenerateSamtoolsVCF(\@finalbams, $out, $refgenome, $c);});
		$snpPool->submit($thr);
		push(@vcfsegs, $out);
	}
	
	waitall;
	$snpPool->joinAll();
	
	$log->Info("SNPFORK", "Beginning bam merger");
	$samexe->MergeSamtoolsVCF(\@vcfsegs, "$outputfolder/vcf/combined.vcf");
	
	if(-s "$outputfolder/vcf/combined.vcf"){
		$log->Info("SNPFORK", "Removing temporary vcf files");
		foreach my $v (@vcfsegs){
			system("rm $v");
		}
	}
	
	if($runSNPEff){
		# TODO: Implement SNPeff annotation of vcf file
	}
}

$log->Info("End", "Run completed successfully");
$log->Close();

exit;
sub fastqcWrapper{
	my ($fastqcParser, $fastqcexe) = @_;
	
	$fastqcParser->runFastqc($fastqcexe);
}

sub samMerge{
	my ($file_array, $outputdir, $log, $sample) = @_;
	
	my $output = "$outputdir/$sample.merged.bam";
	my $header = headerCreation($file_array, "$outputdir/$sample.header.sam");
	my $filestr = join(" ", @{$file_array});
	$log->Info("sammerge", "samtools merge -h $header $output $filestr");
	system("samtools merge -h $header $output $filestr");
	$log->Info("sammerge", "samtools index $output");
	system("samtools index $output");
	
	if(-s $output){
		$log->Info("sammerge", "Cleaning up bam files...");
		foreach my $f (@{$file_array}){
			system("rm $f");
		}
	}
	
	if(-s $output){
		system("rm $header");
	}
}

sub headerCreation{
	my ($file_array, $header_name) = @_;
	my %order = ('@HD' => 0, '@SQ' => 1, '@RG' => 2, '@PG' => 3, '@CO' => 4);
	my @store; # [0,1,2,3] -> []
	my $begin = 1;
	
	foreach my $f (@{$file_array}){
		unless(-s $f){
			print "Could not find file: $f for samtools header creation!\n";
		}
		
		open(IN, "samtools view -H $f | ");
		while(my $line = <IN>){
			chomp $line;
			my @segs = split(/\t/, $line);
			
			if(!exists($order{$segs[0]})){next;}
			if($begin && $segs[0] ne '@RG'){
				my $num = $order{$segs[0]};
				push(@{$store[$num]}, $line);
			}elsif($segs[0] eq '@RG'){
				push(@{$store[2]}, $line);
			}
		}
		close IN;
		$begin = 0;
	}
	
	open(OUT, "> $header_name");
	for(my $x = 0; $x < scalar(@store); $x++){
		my %unique;
		foreach my $row (@{$store[$x]}){
			$unique{$row} = 1;
		}
		# Print out only unique lines!
		foreach my $row (keys(%unique)){
			print OUT "$row\n";
		}
	}
	close OUT;
}

# The workhorse subroutine
sub runBWAAligner{
	# Perl needs to read in the arguments to the subroutine
	my ($fq1, $fq2, $ref, $reffai, $outdir, $base, $lib, $num, $log, $java, $picarddir) = @_;
	# Take care of the names of the input and output files
	my $fqbase = basename($fq1);
	my ($fqStr) = $fqbase =~ /(.+)\.(fastq|fq).*/;
	my $bwasam = "$outdir/$base/$base.$num.sam";
	my $bwabam = "$outdir/$base/$base.$num.bam";
	my $bwadedupbam = "$outdir/$base/$base.$num.nodup.bam";
	#my $bwasort = "$outdir/$base/$base.$num.sorted";
	
	# Create a sub directory if needed
	mkdir("$outdir/$base") || print "";
	
	# Run the BWA mem command and create a bam file from the sam file
	$log->Info("BWAALIGNER", "bwa mem -R '\@RG\\tID:$base.$num\\tLB:$lib\\tPL:ILLUMINA\\tSM:$base\' -v 1 -M $ref $fq1 $fq2 > $bwasam");
	system("bwa mem -R '\@RG\\tID:$base.$num\\tLB:$lib\\tPL:ILLUMINA\\tSM:$base\' -v 1 -M $ref $fq1 $fq2 > $bwasam");
	my $samreader = SamFileReader->new('log' => $log);
	$samreader->prepSam($bwasam);
	#system("samtools view -b -o $bwabam -S $bwasam");
	
	# Sort the bam file
	#print "samtools sort $bwabam $bwasort\n";
	#system("samtools sort $bwabam $bwasort");
	#system("samtools index $bwasort.bam");
	$log->Info("BWAALIGNER", "$java -jar $picarddir/MarkDuplicates.jar INPUT=$bwabam OUTPUT=$bwadedupbam METRICS_FILE=$bwadedupbam.metrics VALIDATION_STRINGENCY=LENIENT");
	system("$java -jar $picarddir/MarkDuplicates.jar INPUT=$bwabam OUTPUT=$bwadedupbam METRICS_FILE=$bwadedupbam.metrics VALIDATION_STRINGENCY=LENIENT");
	
	# If the bam file exists and is not empty, then remove the sam file
	if( -s $bwabam){
		system("rm $bwasam");
		$log->Info("BWAALIGNER", "Completed sorted bam. Removing SAM file");
	}
	
	# $bwabam should be the sorted, indexed bam at the end of the process due to my API
	# I just need to remove the base bam after marking duplicates
	if( -s $bwadedupbam){
		system("rm $bwabam");
		$log->Info("BWAALIGNER", "Completed noduplicate bam. Removing BAM file");
	}
}
