#!/usr/bin/perl
# This script is designed to run BWA alignment + samtools variant calling
# and to profile the time/memory it takes to process each step

use strict;
use Getopt::Std;
use POSIX;
use File::Slurp qw/slurp/;

my $usage = "perl $0 -f <fastq 1> -r <fastq 2> -l <log file> -g <ref genome fasta> -o <output base name>\n";
my %opts;
getopt('frglo', \%opts);

unless(defined($opts{'f'}) && defined($opts{'r'})){
	print $usage;
	exit;
}

open(my $LOG, "> $opts{l}");

print STDERR "ALIGN\n";
# First, the alignment
runProfileCommand("ALIGN", "bwa mem $opts{f} $opts{r} $opts{o}.sam", $LOG);

print STDERR "BAM\n";
# Now the Bam conversion
runProfileCommand("BAM", "samtools view -bS $opts{o}.sam > $opts{o}.bam", $LOG);

print STDERR "SORT\n";
# Now the sorting of the BAM
runProfileCommand("SORT", "samtools sort -T temp.sam -o $opts{o}.sort.bam $opts{o}.bam", $LOG);

print STDERR "INDEX\n";
# Now, index it
runProfileCommand("INDEX", "samtools index $opts{o}.sort.bam", $LOG);

print STDERR "CALL\n";
# Finally, call all variants
runProfileCommand("CALL", "samtools mpileup -go $opts{o}.bcf -f $opts{g} $opts{o}.sort.bam", $LOG);

print STDERR "FILTER\n";
# Now, filter the variants and generate a vcf
runProfileCommand("FILTER", "bcftools call -vmO z -o $opts{o}.vcf.gz $opts{o}.bcf", $LOG);


exit;

sub runProfileCommand{
	my ($cmdname, $cmd, $log) = @_;
	
	my $start = time();
	my $pid = fork();
	if($pid == 0){
		# child process
		system($cmd);
		exit;
	}else{
		sleep(0.5);
		my $attempts = 0;
		my $mem = getMem($pid);
		while($attempts < 5000){
			my $res = waitpid($pid, WNOHANG);
			my $tmem = getMem($pid);
			if($tmem > $mem){
				$mem = $tmem;
			}
			
			sleep(3);
			if($res == -1){
				print {$log} "Error running command: $cmd\n";
				return;
			}elsif($res){
				my $end = time() - $start;
				print {$log} "$cmdname\t$cmd\t$end\t$mem\n";
				last;
			}
			$attempts++;
		}
	}
}

sub getMem{
	my ($pid) = @_;
    	my $memory = split / /, slurp("/proc/$pid/statm");
    	return $memory;
}