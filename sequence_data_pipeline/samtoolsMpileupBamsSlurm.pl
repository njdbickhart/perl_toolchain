#!/usr/bin/perl
# This is a script designed to process a list of bam files, subsection genomic regions and call SNPs and INDELS
# using the samtools mpileup pipeline
# TODO: write code logic and test pipeline

use strict;
use File::Basename;
use Sub::Identify;
use Getopt::Std;
use slurmTools;
use Cwd;

my %opts;
my $processChunks = 1000000; # 1 megabase chunks for variant calling
my @modules = ("samtools/1.3-20-gd49c73b", "bwa/0.7.13-r1126");
my $usage = "perl $0 -b <base outfolder name> -s <OPTIONAL: ref sections for variant calling in UCSC format> -t <input newline separted bam files> -f <input reference fasta file> -m <boolean: generate and queue merger scripts>\n";
getopt('btf', \%opts);

unless(defined($opts{'b'}) && defined($opts{'f'}) && defined($opts{'t'})){
	print $usage;
	exit;
}

my %bcfWorkers; # each section gets its own worker
my @slurmBcfs; # each section bcf gets added to the list
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
my @bcfJobids; # list of queued bcf jobs
foreach my $section (@sections){
	chomp $line; 
	my @segs = split(/\t/, $line);
	
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
	
	my $cmd = "bcftools mpileup -Ob -o $uname.bcf -f $currentDir/$opts{f} --threads = 3 -S $opts{t}";
	push(@slurmBcfs, "$uname.bcf");
	
	$bcfWorkers{$section}->createGenericCmd($cmd, "mpileup_$section");
	$scriptCounter++;
	$bcfWorkers{$section}->queueJobs;
	push(@bcfJobids, @{$bcfWorkers{$section}->jobids});
}

my $numSamples = scalar(keys(%bcfWorkers));

print "Generated $scriptCounter mpileup scripts for $numSamples samples!\n";

if(defined($opts{'m'})){
	# Now to generate the dependency merge scripts
	foreach my $samples (keys(%slurmWorkers)){
		# queue up alignment scripts
		$slurmWorkers{$samples}->queueJobs;
		my $jobids = $slurmWorkers{$samples}->jobIds;
		my $jobNum = scalar(@{$jobids});
		
		print "Sample: $samples has $jobNum queued jobs. JobIds: " . join(" ", @{$jobids}) . "... ";
		
		my $merger = slurmTools->new('workDir' => "$currentDir/$opts{b}/$samples", 
			'scriptDir' => "$currentDir/$opts{b}/$samples/scripts", 
			'outDir' => "$currentDir/$opts{b}/$samples/outLog", 
			'errDir' => "$currentDir/$opts{b}/$samples/errLog",
			'modules' => ["samtools/1.3-20-gd49c73b"],
			'dependencies' => $jobids,
			'useTime' => 0,
			'nodes' => 1,
			'tasks' => 7,
			'mem' => 9000);
			
		my @bams = @{$slurmBams{$samples}};
		my $bwhitespace = join(" ", @bams);
		my @cmds;
		foreach my $b (@bams){
			push(@cmds, "samtools index $b");
		}
		
		push(@cmds, "samtools merge -c -p -@ 6 $samples.sorted.merged.bam $bwhitespace");
		push(@cmds, "samtools index $samples.sorted.merged.bam");
		
		$merger->createArrayCmd(\@cmds, "samMerger");
		$merger->queueJobs;
		print "queued merger script with dependencies!\n";
	}

}

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
