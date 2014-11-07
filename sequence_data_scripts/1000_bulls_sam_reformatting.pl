#!/usr/bin/perl
# This script is designed to process raw sam data in a threaded fashion
# It takes a bam file and formats it according to the 1000 bulls format specification.
# Namely it:
# 	Changes "chr" designations to the uppercase
#	Removes reads that have 3 or more N's in the sequence
#	Removes reads with an average Phred quality score of less than 20
#	Trims low quality bases from the 3' end of a read
#	Removes trimmed reads that are less than 50 bp in length
#	Puts the International ID of the animal in the RG sample tag, and prints out the BAM with the international sample ID in the filename

use Forks::Super;
use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -i <list of input bam file locations> -o <base output directory> -t <number of threads>\n
	Input file list is formatted by two columns tab-delimited
	1. file location
	2. official animal name\n";
	
getops('iot', \%opts);

unless(defined($opts{'i'}) && defined($opts{'o'}) && defined($opts{'t'})){
	print $usage;
	exit;
}

my $threads = $opts{'t'};

my $num = 0;
open(LIST, "< $opts{i}") || die "Could not open list file!\n$usage";
while(my $l = <LIST>){
	chomp $l;
	my @lsegs = split(/\t/, $l);
	fork { sub => \&processBamFile, args => [$lsegs[0], $opts{'o'}, $lsegs[1], $num], max_proc => $threads};
	$num++;
}
close LIST;
waitall();
exit;

sub processBamFile{
	my ($bam, $outbase, $anname, $tnum) = @_;
	open(OUT, "> $outbase/$anname.reformatted.sorted.sam");
	print STDERR "Thread:\t$tnum\tWorking on SAM header\n";
	open(IN, "samtools view -H $bam | ") || die "could not open $bam bam file!\n";
	# Header reformatting first
	while(my $line = <IN>){
		chomp $line;
		if($line =~ /^\@SQ/){
			my @segs = split(/\t/, $line);
			$segs[1] =~ s/c/C/g;
			print OUT join("\t", @segs);
			print OUT "\n";
		}elsif($line =~ /^\@RG/){
			my @segs = split(/\t/, $line);
			foreach my $s (@segs){
				if($s =~ /^SM/){
					my @b = split(/:/, $s);
					$b[1] = $anname;
					$s = join(":", @b);
				}
			}
			print OUT join("\t", @segs);
			print OUT "\n";
		}else{
			print OUT "$line\n";
		}
	}
	close IN;
	
	# Now for the read formatting
	print STDERR "Thread:\t$tnum\tWorking on SAM body\n";
	open(IN, "samtools view $bam |") || die "could not open $bam bam file!\n";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$segs[2] =~ s/c/C/g;
		my $ncount = ($segs[9] =~ tr/N/N/);
		# First check the "N's" and skip reads that don't meet our criteria
		if($ncount >= 3){
			next;
		}else{
			# Now, check the average quality score
			my $avgqs = getAvgQS($segs[10]);
			if($avgqs < 20){
				next;
			}else{
				# OK, we're now at the read trimming stage
				my ($tseq, $tqs) = threePrimeTrimmer($segs[9], $segs[10]);
				if(length($tseq) < 50){
					next; # This read is too small
				}
				if(length($tseq) != length($segs[9])){
					# Now we need to modify the cigar before proceeding
					my $tcigar = cigarChanger($segs[5], length($segs[9]) - length($tseq));
					$segs[4] = $tcigar;
				}
				$segs[9] = $tseq;
				$segs[10] = $tqs;
			}
		}
		print OUT join("\t", @segs);
		print OUT "\n";
	}
	close IN;
	close OUT;
	
	# Now, the samtools to bam conversion
	print STDERR "Thread:\t$tnum\tSAM to BAM conversion\n";
	system("samtools view -bS $outbase/$anname.reformatted.sorted.sam > $outbase/$anname.reformatted.sorted.bam");
	system("rm $outbase/$anname.reformatted.sorted.sam");
	print STDERR "Thread:\t$tnum\tBAM indexing\n";
	system("samtools index $outbase/$anname.reformatted.sorted.bam");
}

sub cigarChanger{
	my ($cigar, $lendiff) = @_;
	my @matches = $cigar =~ /(\d+[MISD])/g;
	for(my $x = scalar(@matches) -1; $x >= 0; $x--){
		# Go through the cigar groups backwards in order to subtract elements
		my $c = $matches[$x] =~ /(\d+)[MISD]/;
		if($c > $lendiff){
			my $h = $matches[$x] =~ /\d+([MISD])/;
			$matches[$x] = ($c - $lendiff) . $h;
			last;
		}else{
			$matches[$x] = "";
			$lendiff -= $c;
		}
	}
	return join("", @matches);
}

sub threePrimeTrimmer{
	my ($seq, $qstr) = @_;
	my $len = 1;
	my @qsbases = split(//, $qstr);
	for(my $x = scalar(@qsbases) - 1; $x >= 0; $x--){
		my $qs = ord($qsbases[$x]) - 33;
		if($qs > 20){
			$len = $x + 1;
			last;
		}
	}
	return (substr($seq, 0, $len), substr($qstr, 0, $len));
}

sub getAvgQS{
	my ($qstr) = @_;
	my $num = length($qstr);
	if($num == 0){return 0;}
	my $count = 0;
	my @bases = split(//, $qstr);
	foreach my $b (@bases){
		$count += ord($b);
	}
	return ($count / $num) - 33;
}