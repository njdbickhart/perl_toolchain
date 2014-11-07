#!/usr/bin/perl
# This is my rewriting of Bin Zhu's script

use FileHandle;
use strict 'vars';
use Getopt::Std;
use List::Util qw[min max];
use vars qw($program $version $description);
$program='cnv_sd_correlation.pl';
$version='0.1';
$description="$program (ver:$version, Author: Bin Zhu) determines the significance of colocolization between CNVs and SD.";
STDOUT->autoflush(1);

use vars qw($opt_c $opt_s $opt_r $opt_l $opt_o $opt_i);


#DEFAULT PARAMETERS#

if ( !defined $ARGV[0] ) {
print "usage: $0 -c [input CNV records] -s [input SD records] -r [number of replicates] -o [output file]
DESCRIPTION
$description
ARGUMENTS
-c input observed CNV records
-s input observed SD records
-r number of replicates for each random CVN selection

OPTIONAL:
-l maximum extended length of SD records on both ends
-i increment of length extension
-o out put file (default: $opt_i.out)
Example: cnv_sd_correlation.pl -c high_confidence_CNV_cow4.txt -s high_confidence_SD_autosome.txt -r 1000 -i 100000 -l 5000000
\n";
exit ();
}


###GET OPTIONS
getopts('c:s:r:l:o:i:');

defined $opt_c || die "Please use -c to input observed CNV file.\n";
defined $opt_s || die "Please use -s to input observed SD file.\n";
defined $opt_r || die "Please use -r to set number of replicates.\n";

$opt_o ||="sd_cnv_correlation.stats";
$opt_l ||= 4000000;
$opt_i ||= 200000;

# Notice: chrlength[N] refers to chrN, chrlength[0] refers to length of chrX, chrlength[30] refers to the length of chrUN
my @chrlength =(88516663, 
		161106243, 
		140800416, 
		127923604, 
		124454208, 
		125847759, 
		122561022, 
		112078216, 
		116942821, 
		108145351, 
		106383598, 
		110171769, 
		85358539, 
		84419198, 
		81345643, 
		84633453, 
		77906053, 
		76506943, 
		66141439, 
		65312493, 
		75796353, 
		69173390, 
		61848140, 
		53376148, 
		65020233, 
		44060403, 
		51750746, 
		48749334, 
		46084206, 
		51998940, 
		283544868
		);  
		
open (CNV, "$opt_c") || die "Can't open CNV file ($opt_c)!\n";
open (SD, "$opt_s") || die "Can't open SD file ($opt_s)!\n";
open (OUT, ">$opt_o") || die "Can't open output file ($opt_o)!\n";

my @c;
my %sd;
my %cnv;
my @randcnv;
my $overlapcount;

my $header = <SD>;
while (my $line = <SD>) {
  		$line =~ s/\r\n/\n/;
  		chomp $line;
  		my @col= split(/\t/, $line);
  		
  		next if $col[1] !~ /\d/; #skip if $col[3] or $col[4] is not digits/numbers
  		next if $col[2] !~ /\d/;
  		next if $col[0] == 0;
  		($col[1],$col[2])=($col[2],$col[1]) if  $col[1]>$col[2]; #sort $col[1] and $col[2]
  		push (@{$sd{$col[0]}}, [$col[1], $col[2], $col[3]]); #put pointers of arrays into an array: ie array of arrays
}

my $header = <CNV>;
while (<CNV>){
               s/\r\n/\n/;
	       chomp;
	       my @col= split /\t/;
	         	
	       next if $col[1] !~ /\d/; #skip if $col[3] or $col[4] is not digits/numbers
	       next if $col[2] !~ /\d/;
	       next if $col[0] == 0;		
	       ($col[1],$col[2])=($col[2],$col[1]) if  $col[1]>$col[2]; #sort $col[1] and $col[2]
  	       push (@{$cnv{$col[0]}}, [$col[1], $col[2], $col[3]]);
}

# Check CNV and SD colocalization
for (my $i = 0; $i <= $opt_l; $i+=$opt_i){

# Extend SDs if not the first iteration
	if ($i > 0){
		foreach my $chr (keys(%sd)){
			foreach my $row (@{$sd{$chr}}){
				$row->[0] = $row->[0] - $opt_i > 0 ? ($row->[0] - $opt_i) : 1;
				$row->[1] = $row->[1] + $opt_i > $chrlength[$chr] ? $chrlength[$chr] : $row->[1] + $opt_i;
			}
		}
	}
# Merge and reread data
	my %mergesd;
	open (TEMP, "> temp");
	foreach my $chr (keys(%sd)){
		foreach my $row (@{$sd{$chr}}){
			print TEMP "chr$chr\t$row->[0]\t$row->[1]\t$row->[2]\n";
		}
	}
	close (TEMP);
	open(BED, "mergeBed -i temp |") || die "Could not open merged bed output\n";
	while (my $line = <BED>){
		chomp $line;
		my @segs = split(/\t/, $line);
		my ($c_v) = $segs[0] =~ /chr(.+)/;
		push(@{$mergesd{$c_v}}, [$segs[1], $segs[2], $segs[3]]);
	} 
	close(BED);
	%sd = %mergesd;
	
	
# Check overlap count for initial CNV and SD
	open(CTEMP, "> temp1");
	open(STEMP, "> temp2");
	for (my $x = 1; $x < 30; $x++){
		my $c_v = "chr$x";
		foreach my $crow (@{$cnv{$x}}){
			print CTEMP "$c_v\t$crow->[0]\t$crow->[1]\t$crow->[2]\tC\n";
		}
		foreach my $srow (@{$sd{$x}}){
			print STEMP "$c_v\t$srow->[0]\t$srow->[1]\t$srow->[2]\tS\n";
		}
	}
	close (CTEMP);
	close (STEMP);
	open (BED, "intersectBed -a temp2 -b temp1 -wa |") || die "Could not open intersect bed output\n";
	my $ocount = 0;
	while (my $line = <BED>){
		$ocount++;
	}
	close BED;
	print OUT "$ocount\t";
# Generate random non-overlapping CNVs and count the overlap
	for (my $k = 0; $k < $opt_r; $k++){
		my $scount = 0;
		my %randcnv = gen_random_cnv(\%cnv);
		open (TEMP, "> temp1");
		foreach my $chr (keys(%randcnv)){
			foreach my $row (@{$randcnv{$chr}}){
				print TEMP "chr$chr\t$row->[0]\t$row->[1]\t$row->[2]\tC\n";
			}
		}
		close (TEMP);
		open (BED, "intersectBed -a temp2 -b temp1 -wa |") || die "Could not open intersect bed sim output\n";
		while (my $line = <BED>){
			$scount++;
		}
		if ($k < $opt_r - 1) { print OUT "$scount\t"; } else { print OUT "$scount\n"; }
		my $knum = $k + 1;
		my $inum = $i/$opt_i;
		print "\e[K";
		print "No. $inum extension and No $knum random CNV maps generated\r";
	}
	
}
print "\n";
system("rm temp temp1 temp2");
close CNV;
close SD;
close OUT;
exit;

sub gen_random_cnv{

   my $cref = shift(@_);
   my %cnv = %{$cref};
   my %randcnv;
   my $randchr;
   my $randpos;
   my $overlap;
   #my %nonoverlap;
   foreach my $chr (keys(%cnv)){
   	foreach my $row (@{$cnv{$chr}}){
   		my $diff = $row->[1] - $row->[0];
   		$overlap = 1;
   		while($overlap){
   			$randchr = &random_int_between(1,29);
       			$randpos = &random_int_between(1, ($chrlength[$randchr]-$diff));
       			my @record = ($randchr, $randpos, $randpos+$diff, "C");
       			if(!exists($randcnv{$randchr})){
       				push(@{$randcnv{$randchr}},[($randpos, $randpos + $diff, "C")]);
       				$overlap = 0;
       			} else {
       				$overlap = overlap_removal(\@record, \%randcnv);
       				if ($overlap == 0){
       					push(@{$randcnv{$randchr}}, [$randpos, $randpos + $diff, "C"]);
       				}
       			}
       		}
       	}
   }

   return %randcnv;

}

sub random_int_between {

    my($min, $max) = @_;
    # Assumes that the two arguments are integers themselves!
    return $min if $min == $max;
    ($min, $max) = ($max, $min)  if  $min > $max;
    return $min + int rand(1 + $max - $min);
    
} 

sub overlap_removal {  

   my ($a_ref, $h_ref) = @_;
   my @test = @{$a_ref};
   my %conf = %{$h_ref};
   my $count = 0;
	foreach my $row (@{$conf{$test[0]}}){
		if ($test[2] > $row->[0] && $test[1] < $row->[1]){
			return 1;
		}
	}
	return 0;
}


sub isempty { return (!defined($_[0]) || $_[0] eq ''); }