#!/usr/bin/perl
# This is a script designed to process both fasta and fastq files to calculate a histogram of average GC content of reads
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5000
use Getopt::Std;

my $usage = "perl $0 -f <gzipped fastq or fasta file> -o <output gc occurance histogram>\n";

my %opts;
getopt('fo', \%opts);

unless(defined($opts{'f'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

# check file format
my $fastq = 0;
open(my $TEST, "gunzip -c $opts{f} |");
for(my $x = 0; $x < 2; $x++){
	my $temp = <$TEST>;
}
# Get the third line and check to see if it's a fastq
my $third = <$TEST>;
if($third =~ /^+/){
	$fastq = 1;
}
seek($TEST, 0, 0);
close $TEST;

# Now open the file and count GC
my %gcBins;
open(my $IN, "gunzip -c $opts{f} |");
if($fastq){
	while(my $head = <$IN>){
		my $seq = <$IN>;
		chomp $seq;
		my $gc = calcGC($seq);
		$gcBins{$gc} += 1;
		for(my $x = 0; $x < 2; $x++){
			<$IN>; # Get rid of the next two lines
		}
	}
	close $IN;
}else{
	# This is not the best way, but it's fast to code and memory efficient:
	# I'm calculating GC content per fasta line, rather than per fasta entry
	while(my $line = <$IN>){
		if($line =~ /^>/){
			next;
		}
		chomp($line);
		my $gc = calcGC($line);
		$gcBins{$gc} += 1;
	}
	close $IN;
}

# Now to print out the GC bins
open(my $OUT, "> $opts{o}");
print {$OUT} "GC\tCount\n";
foreach my $gc (sort{$a <=> $b} keys(%gcBins)){
	print {$OUT} "$gc\t$gcBins{$gc}\n";
}
close $OUT;


exit;

sub calcGC{
	my ($seq) = @_;
	my $count = ($seq =~ tr/GCgc/GCgc/);
	$count /= length($seq);
	return sprintf("%.3f", $count);
}
