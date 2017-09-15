#!/usr/bin/perl
# This is a script designed to take a readname sorted bam, align the unmapped mate to a new reference and tabulate the mapping information
# Two output files are generated: (1) a ".reads" file that compares "stacked" read alignments and (2) a ".segs" file that has consensus read mappings
# Version 1.0.2: addition of mapq filters.

use strict;
use Getopt::Std;

my $usage = "perl $0 -b <input readname sorted bam> -r <alternate reference genome> -o <output base name> -m <mapping quality lower threshold> -t <bwa thread count>\n";
my %opts;

getopt('bromt', \%opts);

unless(defined($opts{'b'}) && defined($opts{'r'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $bwathreads = (defined($opts{'t'}))? $opts{'t'} : 1;

# Mapping quality threshold default
my $mapq = 30;
if(defined($opts{'m'})){
	$mapq = $opts{'m'};
}

my %data; # {read name} -> [chr1, start1, end1, chr2, start2, end2]
# open bam file and queue up alignments
open(my $IN, "samtools view $opts{b} |");
open(my $OUT, "> temp.fastq") || die "Could not open temporary fastq file!\n";
while(my $read1 = <$IN>){
	chomp $read1; 
	my @r1segs = split(/\t/, $read1);
	
	my $read2 = <$IN>;
	chomp $read2;
	my @r2segs = split(/\t/, $read2);
	
	# flag for which read is unmapped
	my $r1unmapped = ($r1segs[1] & 4 == 4)? 1 : 0;
	
	# End coordinates will be the read length + the mapping position -- it's not exact, but a good approximation for now
	if($r1unmapped){
		# If the anchor is a low quality alignment, skip
		if($r2segs[4] <= $mapq || $r2segs[4] == 255){next;}
		print {$OUT} "\@$r1segs[0]\n$r1segs[9]\n+\n$r1segs[10]\n";
		push(@{$data{$r2segs[0]}}, $r2segs[2], $r2segs[3], $r2segs[3] + length($r2segs[9]));
	}else{
		# If the anchor is a low quality alignment, skip
		if($r1segs[4] <= $mapq || $r1segs[4] == 255){next;}
		print {$OUT} "\@$r2segs[0]\n$r2segs[9]\n+\n$r2segs[10]\n";
		push(@{$data{$r1segs[0]}}, $r1segs[2], $r1segs[3], $r1segs[3] + length($r1segs[9]));
	}
}

close $IN;
close $OUT;

print "Generated unmapped read fastq -- aligning to alternate reference: $opts{r}\n";

# Now map the reads and generate the table
open(my $IN, "bwa mem -t $bwathreads $opts{r} temp.fastq |");
while(my $align = <$IN>){
	chomp $align;
	# Addition below to remove the SAM header flags
	if($align =~ /^@/){next;}
	my @asegs = split(/\t/, $align);
	
	# if secondary alignment is bad, then fill in placeholder data
	if($asegs[4] <= $mapq || $asegs[4] == 255){
		push(@{$data{$asegs[0]}}, "*", 1, 1);
	}
	
	# Naieve merge of read names with prior information
	push(@{$data{$asegs[0]}}, $asegs[2], $asegs[3], $asegs[3] + length($asegs[9]));
}
close $IN;

# Now to reorder the data in the table so that we can generate a consensus
print "Reordering raw data table...\n";
# Sorting by mapped original chr and position and then merging

# Added more sort criteria for condensing segments together
my @sortedData = sort{
	$a->[0] cmp $b->[0] ||
	$a->[1] <=> $b->[1] ||
	$a->[3] cmp $b->[3] || 
	$a->[4] <=> $b->[4]} values(%data);

print "Data resorted, now condensing...\n";
my @condensedData;
push(@condensedData, $sortedData[0]);
for(my $x = 1; $x < scalar(@sortedData); $x++){
	my $row = $sortedData[$x];

	# For each read, if both mapping locations are within 500 bp of each other, merge the coordinates
	if($row->[0] eq $condensedData[-1]->[0] &&
	   $row->[3] eq $condensedData[-1]->[3] &&
	   overlap($row->[1], $row->[2], $condensedData[-1]->[1], $condensedData[-1]->[2]) < 500 &&
	   overlap($row->[4], $row->[5], $condensedData[-1]->[4], $condensedData[-1]->[5]) < 500){
		# The overlap was successful, now merging coordinates
		$condensedData[-1]->[1] = min($row->[1], $condensedData[-1]->[1]);
		$condensedData[-1]->[2] = max($row->[2], $condensedData[-1]->[2]);
		$condensedData[-1]->[4] = min($row->[4], $condensedData[-1]->[4]);
		$condensedData[-1]->[5] = max($row->[5], $condensedData[-1]->[5]);
	}else{

		# Did not overlap, creating a separate entry
		push(@condensedData, $sortedData[$x]);
	}
}

# Printing out the results, first the raw file
open($OUT, "> $opts{o}.reads");
foreach my $key (keys(%data)){
	my @row = @{$data{$key}};
	print {$OUT} "$key\t";
	print {$OUT} join("\t", @row);
	print {$OUT} "\n";
}
close $OUT;

# Now the condensed coords
open($OUT, "> $opts{o}.segs");
foreach my $row (@condensedData){
	print {$OUT} join("\t", @{$row});
	print {$OUT} "\n";
}
close $OUT;

print "Program complete. Read alignment comparisons are in: $opts{o}.reads and consensus coordinates are in: $opts{o}.segs\n";

#system("rm temp.fastq");
exit;

sub overlap{
	my ($s1, $e1, $s2, $e2) = @_;
	my $ovlp = max($s1, $s2) - min($e1, $e2);
	return $ovlp;
}

sub max{
	my ($p1, $p2) = @_;
	my $max = ($p1 > $p2)? $p1 : $p2;
	return $max;
}

sub min{
	my ($p1, $p2) = @_;
	my $min = ($p1 < $p2)? $p1 : $p2;
	return $min;
}