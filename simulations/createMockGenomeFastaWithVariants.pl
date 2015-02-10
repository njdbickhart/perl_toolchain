#!/usr/bin/perl
# This script is designed to create two genome fasta files for use with wgsSim
# Heterozygosity/homozygosity of variants will be determined and generated in the subsequent fastas
# To keep thing simple, genome fasta 1 will always have the sole copy of the heterozygous variants from the known list
# I will use WGSIM's INDEL and SNP creation utilities to generate random events
# WGSIM mutation output format:
#	chr	pos	ref	variant	zygosity ("-" = hom, "+" = het)
# I will use WGSIM's '-h' option in order to only generate homozygous mutations in each haploid 'chromsome'

use strict;
use Getopt::Std;
use perlBed;
use Forks::Super;

# The 4th through 9th columns of the variant bed file
my @variantCols = ("SVtype", "subtype", "size", "ref", "alt", "MAF");
my @posFileCols = ("SVtype", "subtype", "size", "ref", "alt", "zygosity", "knownSite");

my %opts;
my $usage = "perl $0 [All arguments are mandatory]
	-v <known variant file> 
	-g <input, base genome fasta file>
	-d <output directory> 
	-o <output basename for files> 
	-m <mutation rate; approximate number of random SNPs/INDELs per base>
	-i <indel rate; out of the mutation rate, the proportion that are INDELs vs SNPs>
	-n <target X cov>
	-l <read length>
	-e <sequencing error rate>
	-b <debug flag (prevents erasure of temp files)>\n";
getopt('vgordminle', \%opts);

my @missing;
foreach my $key ("v", "g", "o", "r", "d", "m", "i", "n", "l", "e"){
	if(!defined($opts{$key})){
		push(@missing, $key);
	}
}

if(scalar(@missing) > 0){
	print "Missing following mandatory arguments! " . join(", ", @missing) . "\n$usage";
	exit;
}

# Define files, make output directory
mkdir($opts{'d'}) || print STDERR "[Directory] $!\n";
my $genomeFa1 = "$opts{d}/$opts{o}.1.fa";
my $genomeFa2 = "$opts{d}/$opts{o}.2.fa";
my $fastq11 = "$opts{d}/$opts{o}.1.1.fq";
my $fastq12 = "$opts{d}/$opts{o}.1.2.fq";
my $fastq21 = "$opts{d}/$opts{o}.2.1.fq";
my $fastq22 = "$opts{d}/$opts{o}.2.2.fq";
my $wgsimLocs1 = "$opts{d}/$opts{o}.wgsim.1.locs";
my $wgsimLocs2 = "$opts{d}/$opts{o}.wgsim.2.locs";
# Container for temporary files to be deleted later
my @temps = ($genomeFa1, $genomeFa2, $fastq11, $fastq12, $fastq21, $fastq22, $wgsimLocs1, $wgsimLocs2);


my $finalFq1 = "$opts{d}/$opts{o}.1.fq";
my $finalFq2 = "$opts{d}/$opts{o}.2.fq";
my $truePosBed = "$opts{d}/$opts{o}.truepos.bed";

# Load known variant locations
#my $known = BedContainer->new();
#$known->loadFile($opts{'v'}, \@variantCols);
my %beds;
my $numknown = 0;
open(IN, "< $opts{v}") || die "Could not open variant locations file! $opts{v}\n";
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	
	# Determine which variants will be represented in this individual
	my $use = rand(1);
	my $hom = rand(1);
	
	# check if random value is below the MAF threshold
	if($use <= $segs[-1]){
		my $zygo = "het";
		# Now, determine zygosity
		if($hom <= 0.33){
			# Simple mendelian here. 1/3 of positive mutations will be homozygous
			$zygo = "hom";
		}
		my %values = (
			"SVtype" => $segs[3],
			"subtype" => $segs[4],
			"size" => $segs[5],
			"ref" => $segs[6],
			"alt" => $segs[7],
			"zygosity" => $zygo,
			"knownSite" => 1);
		$numknown++;
		push(@{$beds{$segs[0]}}, BedCoord->new('chr' => $segs[0], 'start' => $segs[1], 'end' => $segs[2], 'values' => \%values));
	}
}
close IN;


print STDERR "Will use $numknown known variants for this simulation dataset.\n";


# Ensure that the input genome is faidxed
if(-s "$opts{g}.fai"){
	print STDERR "Generating samtools faidx for genome file: $opts{g}...\n";
	system("samtools faidx $opts{g}");
}

# Calculate number of reads from read length and x coverage
my $genomesize = 0;
my %chrsizes;
open(IDX, "< $opts{g}.fai") || die "Could not open fasta index file!\n";
while(my $line = <IDX>){
	chomp $line;
	my @segs = split(/\t/, $line);
	$genomesize += $segs[1];
	$chrsizes{$segs[0]} = $segs[1];
}
close IDX;
my $readnum = int($genomesize / ($opts{'l'} * 2));

# Create both variant genome fastas
my ($hetandhom, $hom) = createBothVariantFastas($opts{g}, $genomeFa1, $genomeFa2, \%chrsizes, \%beds);

# Run two threads of wgsim to generate the fastqs
fork{ sub => \&wgsimWrapper, args => [$genomeFa1, $fastq11, $fastq12, $wgsimLocs1, $opts{'e'}, $opts{'m'}, $opts{'i'}, $opts{'l'}, $readnum]};

fork{ sub => \&wgsimWrapper, args => [$genomeFa2, $fastq21, $fastq22, $wgsimLocs2, $opts{'e'}, $opts{'m'}, $opts{'i'}, $opts{'l'}, $readnum]};

waitall();

# Read the location files and generate bed coords for later printout
# I'm going to be lazy here and not try to resolve coordinate overlaps between the two files
# Damn, I do need to correct the coordinate locations because the mutant fastas have different base counts
my @locfiles = ($wgsimLocs1, $wgsimLocs2);
my @variantSites = ($hetandhom, $hom);
for (my $x = 0; $x < 2; $x++) {
	my $file = $locfiles[$x];
	my $vars = $variantSites[$x];
	open(IN, "< $file") || print STDERR "Could not find wgsim loc file: $file!\n";
	my $varidx = 0;
	my $varoffset = 0; # Number of bases to subtract from the current location to bring parity between the fasta references
	my $curvar = "NA";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if($curvar ne "NA"){			
			if($curvar->chr ne $segs[0]){
				$varidx = 0;
				$varoffset = 0;
			}
		}	
		
		$curvar = $vars->{$segs[0]}->[$varidx];
		
		while($curvar->start <= $segs[1]){
			# Only change the offset if the current variant site is AFTER the known variant location!
			$varoffset = updateOffset($curvar, $varoffset);
			$varidx++;
			$curvar = $vars->{$segs[0]}->[$varidx];
		}
		
		# Determine SV type
		my $SVtype = "SNP";
		my $subtype = "SNP";
		my $size = 1;
		if($segs[2] eq '-'){
			$SVtype = "INDEL";
			$subtype = "INS";
			$size = length($segs[3]);
		}elsif($segs[3] eq '-'){
			$SVtype = "INDEL";
			$subtype = "DEL";
			$size = length($segs[2]);
		}
		
		my %values = (
			"SVtype" => $SVtype,
			"subtype" => $subtype,
			"size" => $size,
			"ref" => $segs[2],
			"alt" => $segs[3],
			"zygosity" => "het",
			"knownSite" => 0);
		
		my $updateStart = $segs[1] - $varoffset;
		if($updateStart < 1){
			print STDERR "Insertion near beginning of the mutant file! First variant position will be wrong for $file!\n";
			$updateStart = 1;
		}
		my $end = $updateStart + $size;
		push(@{$beds{$segs[0]}}, BedCoord->new('chr' => $segs[0], 'start' => $updateStart, 'end' => $end, 'values' => \%values));
	}
	close IN;
}

# Combine the four fastqs into a pair of final fastqs
combineFastqs($fastq11, $fastq21, $finalFq1);
combineFastqs($fastq12, $fastq22, $finalFq2);

# Print out the true and wgsim variant locations
open(OUT, "> $truePosBed");
my $container = BedContainer->new();
foreach my $chr (keys %beds){
	$container->addBedArray($beds{$chr});
}

my @sorted = $container->getOrderedList();
foreach my $b (@sorted){
	print OUT $b->outString(\@posFileCols);
}
close OUT;

unless($opts{'b'}){
	print STDERR "Deleting temp files...\n";
	foreach my $f (@temps){
		system("rm $f");
	}
}

exit;

# Really primitive wrapper
sub combineFastqs{
	my ($fq1, $fq2, $output) = @_;
	system("cat $fq1 $fq2 > $output");
}
	

# Small sub to calculate mutant fasta offset value from wgsim
sub updateOffset{
	my ($var, $offset) = @_;
	if($var->get_value("SVtype") eq "INDEL"){
		if($var->get_value("subtype") eq "INS"){
			return $offset + $var->get_value("size");
		}elsif($var->get_value("subtype") eq "DEL"){
			return $offset - $var->get_value("size");
		}
	}
	return $offset;
}

sub createBothVariantFastas{
	my ($genomefa, $output1fa, $output2fa, $chrlist, $vlist) = @_;
	
	my %hetAndhom;
	my %hom;
	
	# Divide the known variants into two groups
	foreach my $chr (keys %{$vlist}){
		foreach my $bed (@{$vlist->{$chr}}){
			if($bed->get_value("zygosity") eq "hom"){
				push(@{$hom{$chr}}, $bed);
			}
			push(@{$hetAndhom{$chr}}, $bed);
		}
	}
	
	# create fasta 1
	generateVariantFastaLONG($chrlist, $genomefa, $output1fa, \%hetAndhom);
	
	# create fasta 2
	generateVariantFastaLONG($chrlist, $genomefa, $output2fa, \%hom);
	
	return (\%hetAndhom, \%hom);
}

sub generateVariantFastaLONG{
	my ($chrlist, $genomefa, $outputfa, $vlist) = @_;
	open(OUT, "> $outputfa");
	foreach my $chr (keys %{$chrlist}){
		my $chrsize = $chrlist->{$chr};
		print STDERR "Working on generating $chr for $outputfa...\n";
		my @variants = $vlist->{$chr};
		
		# This is a really terrible, brute force way, but it's the quickest to code
		open(IN, "samtools faidx $genomefa $chr:1-$chrsize |");
		my @seq;
		my $h = <IN>;
		while(my $line = <IN>){
			chomp $line;
			my @segs = split(//, $line);
			push(@seq, @segs);
		}
		close IN;
		
		print OUT ">$chr\n";
		my @buffer;
		#my $printout = 0;
		for(my $x = 0; $x < scalar(@seq); $x++){
			if(scalar(@variants) > 0){
				my $curvar = $variants[0];
				if($x < $curvar->start()){
					push(@buffer, $seq[$x]);
				}elsif($x == $curvar->start()){
					# The start position equals our selected variant location
					if($curvar->get_value("SVtype") eq "INDEL"){
						if($curvar->get_value("subtype") eq "INS"){
							#my $s = $printout;
							my $alt = $curvar->get_value("alt");
							my @asegs = split(//, $alt);
							for(my $y = 0; $y < scalar(@asegs); $y++){
								push(@buffer, $asegs[$y]);
								if(scalar(@buffer) >= 60){
									print OUT join('', @buffer) . "\n";
									@buffer = ();
									#$printout += 60;
								}
							}
							#print STDERR $curvar->printCoords() . "\t$s\t$len\n";
							#$x = $curvar->end();
							shift(@variants);
						}elsif($curvar->get_value("subtype") eq "DEL"){
							$x = $curvar->end();
							#print STDERR $curvar->printCoords() . "\t$printout\n";
							shift(@variants);
						}
					}elsif($curvar->get_value("SVtype") eq "SNP"){
						# Really easy substitution
						push(@buffer, $curvar->get_value("alt"));
					}	
				}
				if(scalar(@buffer) >= 60){
					print OUT join('', @buffer) . "\n";
					@buffer = ();
					#$printout += 60;
				}
			}else{
				push(@buffer, $seq[$x]);
				if(scalar(@buffer) >= 60){
					print OUT join('', @buffer) . "\n";
					@buffer = ();
					#$printout += 60;
				}
			}
		}
		if(scalar(@buffer) > 0){
			print OUT join('', @buffer) . "\n";
		}
	}
	close OUT;
}

sub wgsimWrapper{
	my ($genomefasta, $fastq1, $fastq2, $locfile, $error, $mut, $indel, $readlen, $readnum) = @_;
	
	system("wgsim -e $error -N $readnum -1 $readlen -2 $readlen -r $mut -R $indel -h $genomefasta $fastq1 $fastq2 > $locfile");
}