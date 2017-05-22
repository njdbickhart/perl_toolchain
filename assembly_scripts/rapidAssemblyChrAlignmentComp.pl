#!/usr/bin/perl
#SBATCH --nodes=1
#SBATCH --mem=10000
#SBATCH --ntasks-per-node=2
# This is a script that aligns sampled portions of a chromosome from one assembly to another to check for continuity

use strict;
use Getopt::Std;

my $usage = "perl $0 -c <comparison assembly> -j <comparison chromosome> -n <new assembly> -m <new assembly chromosome> -o <output tab file> -s <sampling window [10000 bp]>\n";

my %opts;
getopt('cjnmos', \%opts);

unless(defined($opts{'c'}) && defined($opts{'j'}) && defined($opts{'n'}) && defined($opts{'m'})){
	print $usage;
	exit;
}

my $sampling = (defined($opts{'s'}))? $opts{'s'} : 10000;

# check if assembly FAI files exist
checkFAIFile($opts{'c'});
checkFAIFile($opts{'n'});

# Generate sampling fasta
my $sampFasta = "$opts{j}.sampling.fa";
open(my $IN, "module load samtools; samtools faidx $opts{c} $opts{j} |");
open(my $OUT, "> $sampFasta");
my $head = <$IN>;
my $bp = 0; my $totalbp = 0; 
while(my $line = <$IN>){
	chomp $line;
	if($bp >= $sampling){
		$bp = 0;
		print {$OUT} ">sample.$opts{j}.$totalbp\n$line\n";
	}
	$bp += length($line);
	$totalbp += length($line);
}
close $IN;
close $OUT;

# Extract chromosome to compare
my $newASMComp = "$opts{m}.chr.comp.fa";
system("module load samtools; samtools faidx $opts{n} $opts{m} > $newASMComp");
system("module load bwa; bwa index $newASMComp");
# align and order sampling fasta
my @data; # [chr, start, end, orient, origpos, origend]
my $start = 0;
open(my $IN, "module load bwa; bwa mem $newASMComp $sampFasta |");
while(my $line = <$IN>){
	if($line =~ /^@/){next;}
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[1] & 2048){next;}
	if($segs[2] eq "*"){next;}
	
	my @rnamesegs = split(/\./, $segs[0]);
	
	if($start == 0){
		$start = 1;
		push(@data, [$segs[2], $segs[3], $segs[3], ($segs[1] & 16)? "-" : "+", $rnamesegs[2], $rnamesegs[2]]);
		next;
	}else{
		my $compRow = $data[-1];
		my $orient = ($segs[1] & 16)? "-" : "+";
		if($compRow->[0] eq $segs[2] 
			&& ($compRow->[2] > $segs[3] - ($sampling * 3) && $compRow->[2] < $segs[3] + ($sampling * 3))
			&& $compRow->[3] eq $orient){
			
			$data[-1]->[2] = $segs[3];
			$data[-1]->[5] = $rnamesegs[2];
		}else{
			push(@data, [$segs[2], $segs[3], $segs[3], $orient, $rnamesegs[2], $rnamesegs[2]]);
		}
	}
}
close $IN;

open(my $OUT, "> $opts{o}");
foreach my $row (@data){
	my $len = $row->[5] - $row->[4];
	print {$OUT} $row->[0] . "\t" . $row->[1] . "\t" . $row->[2] . "\t" . $row->[3] . "\t$opts{m}\t" . $row->[4] . "\t" . $row->[5] . "\t$len\n";
}
close $OUT;

exit;

sub checkFAIFile{
	my ($asm) = @_;
	unless( -s "$asm.fai"){
		print STDERR "Generating fasta index file for: $asm\n";
		system("module load samtools; samtools faidx $asm");
	}
}
