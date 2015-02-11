#!/usr/bin/perl
# This script is designed to simulate a series of SNPs and INDELs that are known, true-positives
# Output is formatted in the bed style:
# chr	start	end	SVtype	subtype	size	ref	alt	MAF
# (SNP subtype is simply, "SNP", but INDELs are "INS" or "DEL")

use strict;
use perlBed;

# SNPs will be 4 times more likely in this simulation
our @SVs = ("INDEL", "SNP", "SNP", "SNP", "SNP");
our @INDELType = ("INS", "DEL");
our @Bases = ("A", "G", "C", "T");

chomp(@ARGV);
if(scalar(@ARGV) < 4){
	print "Usage: perl $0 <input repetitive bed> <input genome fasta> <output bed file name> <event rate <per chr>>\n";
	exit;
}

# Load the repeat bed file
my $repeats = BedContainer->new();
$repeats->loadFile($ARGV[0]);

# Read the FAI file and generate a hash of chromosomes and lengths
my @chrs;
my @lens;
open(IN, "< $ARGV[1].fai") || die "Could not open fasta index for $ARGV[1] fasta file! Please index with samtools!\n";
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\s+/, $line);
	push(@chrs, $segs[0]);
	push(@lens, $segs[1]);
}
close IN;

# Now, run the main program for each chromosome
for(my $y = 0; $y < scalar(@chrs); $y++){
	mainRoutine($ARGV[1], $chrs[$y], $lens[$y], $ARGV[3], $repeats, $ARGV[2]);
}

print STDERR "Finished with the bed file!\n";

exit;

sub mainRoutine{
	require perlBed;
	my ($fasta, $chr, $len, $eventRate, $repeats, $output) = @_;
	open(OUT, ">> $output");
	
	my $thresh = $len * $eventRate;
	for(my $x = 1; $x < $len - 50; $x++){
		# Check if this intersects with a repeat -- if so, skip
		if($repeats->intersects($chr, $x, $x + 49)){
			next;
		}
	
		my $test = int(rand($len));
		# Check to see if there will be a variant at this position
		if($test <= $thresh){
			# Now to see what variant will be generated
			my $sv = $SVs[int(rand(5) - 0.1)];
			# MAF can only be above 5%
			my $maf = (rand(0.5) - 0.05) + 0.05;
			if($sv eq "INDEL"){
				my $endpos = $x + 49;
				open(IN, "samtools faidx $fasta $chr:$x-$endpos |") || die "Could not run samtools on $chr:$x-$endpos for $fasta!\n";
				my $h = <IN>;
				my $seq = <IN>;
				chomp($seq);
				my ($altseq, $type, $size) = simulateINDEL($seq);
				close IN;
				
				if($type eq "INS"){
					$seq = "-";
				}else{
					$seq = substr($seq, 0, $size);
				}
				
				print OUT "$chr\t$x\t$endpos\t$sv\t$type\t$size\t$seq\t$altseq\t$maf\n";
			}else{
				# Assuming SNP here
				open(IN, "samtools faidx $fasta $chr:$x-$x |") || die "Could not run samtools on $chr:$x-$x for $fasta!\n";
				my $h = <IN>;
				my $seq = <IN>;
				chomp($seq);
				my ($altseq) = simulateSNP($seq);
				close IN;
				
				print OUT "$chr\t$x\t$x\t$sv\t$sv\t1\t$seq\t$altseq\t$maf\n";
			}
		}
	}
	close OUT;
	print STDERR "Finished with chr: $chr...\r";
}
	
	

# Will insert random sequence or delete sequence from a string provided to the subroutine
# Insertions will be prior to the start of the sequence (start position)
# Deletions will start at the start of the sequence (start position, again)
sub simulateINDEL{
	my ($refseq) = @_;
	my $type = $INDELType[int(rand(2) - 0.1)];
	my $size = int(rand(50) - 0.1) + 1;
	if($type eq "INS"){
		my $insseq = generateRandSeq($size);
		return $insseq, $type, $size;
	}else{
		return "-", $type, $size;
	}
}

sub generateRandSeq{
	my ($size) = @_;
	my $seq;
	for(my $x = 0; $x < $size; $x++){
		$seq .= $Bases[int(rand(5) - 0.1)];
	}
	return $seq;
}

# Does not need to account for the position of the base
sub simulateSNP{
	my ($refb) = @_;
	my $altb = "N";
	while($altb ne $refb && $altb ne "N"){
		$altb = $Bases[int(rand(5) - 0.1)];
	}
	return $altb;
}