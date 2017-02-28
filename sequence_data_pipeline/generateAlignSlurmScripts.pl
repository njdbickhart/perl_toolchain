#!/usr/bin/perl
# This is a script designed to process a tab file with sequence data fastqs and produce slurm scripts

use strict;
use File::Basename;
use Sub::Identify;
use Getopt::Std;
use slurmTools;
use Cwd;

my %opts;
my @modules = ("samtools/1.3-20-gd49c73b", "bwa/0.7.13-r1126");
my $usage = "perl $0 -b <base outfolder name> -t <input tab sequence files> -f <input reference fasta file> -m <boolean: generate and queue merger scripts>\n";
getopt('btf', \%opts);

unless(defined($opts{'b'}) && defined($opts{'f'}) && defined($opts{'t'})){
	print $usage;
	exit;
}

my %slurmWorkers; # each sample gets its own worker
my %slurmBams; # each sample gets a list of sorted bams
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
# Generate alignment scripts
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
			'mem' => 12000);
	}
	
	my $bname = basename($segs[0]);
	my @bsegs = split(/\./, $bname);
	# Adding a hash of file paths to ensure that the file names are unique
	my $uHash = urlHash($segs[0]);
	my $uname = "$bsegs[0].$uHash";
	
	my $cmd = "bwa mem -t 5 -R '\@RG\tID:$segs[-2]\tSM:$segs[-1]\tLIB:$segs[-2]' $fasta $segs[0] $segs[1] | samtools view -S -b - | samtools sort -m 2G -o $uname.sorted.bam -T $uname -";
	push(@{$slurmBams{$segs[-1]}}, "$uname.sorted.bam");
	
	$slurmWorkers{$segs[-1]}->createGenericCmd($cmd, "bwaAlign");
	$scriptCounter++;
}
close $IN;
my $numSamples = scalar(keys(%slurmWorkers));

print "Generated $scriptCounter alignment scripts for $numSamples samples!\n";

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
	my ($url) = @_;
	$url =~ s/\///g; # Remove redundant path separators
	my $hash = 0;
	foreach my $h (split //, $url) {
		$hash = $hash*33 + ord($h);
	}
	return $hash;
}