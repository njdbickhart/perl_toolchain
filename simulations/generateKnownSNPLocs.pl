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
our %revcomp = ("A" => "T", "T" => "A", "G" => "C", "C" => "G");

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

# Check if output exists, then delete it and start fresh
if(-s $ARGV[2]){
	system("rm $ARGV[2]");
}

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
	
	my $counter = 0;
	# Just to make things fair, removing the repetitive regions length
	my $repeatLength = $repeats->calculateBedLength($chr);
	my $nonreplen = $len - $repeatLength;
	my $thresh = $len / $eventRate;
	for(my $x = 1; $x < $len - 50; $x++){
		# Check if this intersects with a repeat -- if so, skip
		# I believe that the bed coordinate positions are usually one based
		if($repeats->intersects($chr, $x - 1, $x + 51)){
			# Should speed up calculations by moving past repetitive regions
			my $bed = $repeats->firstIntersect($chr, $x - 1, $x + 51);
			$x = $bed->end;
			next;
		}
	
		my $test = int(rand($thresh) - 0.1);
		# Check to see if there will be a variant at this position
		if($test <= 1){
			# Now to see what variant will be generated
			my $sv = $SVs[int(rand(5) - 0.1)];
			# MAF can only be above 5%
			my $maf = (rand(0.45)) + 0.05;
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
				
				$endpos = $x + $size;
				print OUT "$chr\t$x\t$endpos\t$sv\t$type\t$size\t$seq\t$altseq\t$maf\n";
				#if($type eq "INS"){
					$x += $size + 1;
				#}
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
			$counter++;
			if($counter % 100 == 0){
				print STDERR "Generated $counter loci so far...\r";
			}
		}
	}
	close OUT;
	print STDERR "Finished with chr: $chr...\n";
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
		$seq .= $Bases[int(rand(4) - 0.1)];
	}
	return $seq;
}

# Does not need to account for the position of the base
sub simulateSNP{
	my ($refb) = @_;
	my $altb = $Bases[int(rand(4) - 0.1)];
	


	if($altb eq $refb){
		# I was having difficulty with the loop structure losing references
		$altb = $revcomp{uc($refb)};
	}
	return $altb;
}
