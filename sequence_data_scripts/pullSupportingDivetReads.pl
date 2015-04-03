#!/usr/bin/perl
# This is a pipeline script designed to pull read support from a series of divet files
# in order to corroborate known locations that should display read discordancy
# Primarily designed for regions of no read depth coverage

use strict;
use Getopt::Std;
use Forks::Super;

my $usage = "perl $0 -i <underscore separated, comma separated divet file names<underscore separates pairs, commas separate names from files>> 
	-b <underscore separated, comma separated bam files <same scheme>>
	-c <chromosome>
	-s <start>
	-e <end>
	-o <output directory>\n";
	
my %opts;
getopt('ibocse', \%opts);

unless(defined($opts{'i'}) && defined($opts{'b'}) && defined($opts{'o'}) && defined($opts{'c'})){
	print $usage;
	exit;
}

mkdir($opts{'o'}) || print "$!\n";

my @divetpairs = split(/\_/, $opts{'i'});
my @bampairs = split(/\_/, $opts{'b'});

my %divet = map{split(/,/)} @divetpairs;
my %bam = map{split(/,/)} @bampairs;

print STDERR "Beginning divet filtering\n";
my %pickedDivets;
foreach my $animals (keys(%divet)){
	$pickedDivets{$animals} = "$opts{o}/$animals.selected.divet";
	print STDERR "Working on $opts{o}/$animals.selected.divet\n";
	fork { sub => \&pullDivetSupport, args => [$opts{'c'}, $opts{'s'}, $opts{'e'}, "$opts{o}/$animals.selected.divet", $divet{$animals}] };
}

waitall;
print STDERR "Finished with divet filtering\n";

my %breakpoints; # {animal} -> [start, end]
foreach my $animals (keys(%pickedDivets)){
	my ($start, $end) = determineDelBreakpoints($opts{'c'}, $pickedDivets{$animals});
	$breakpoints{$animals} = [$start, $end];
	print STDERR "Breakpoints for $animals : $opts{c} $start $end\n";
}

print STDERR "Pulling reads for fq generation\n";
foreach my $animals (keys(%breakpoints)){
	my $reads = getBreakpointReads($bam{$animals}, $opts{'c'}, $breakpoints{$animals}->[0], $breakpoints{$animals}->[1]);
	print STDERR "Pulling reads for: $animals $bam{$animals}\n";
	fork { sub => \&createReadFQs, args => [$bam{$animals}, $reads, "$opts{o}/$animals.breakpoint.fq"]};
}

waitall;
print STDERR "Finished with pipeline\n";

exit;

sub createReadFQs{
	my ($bam, $reads, $outputfq) = @_;
	if(scalar(keys(%{$reads})) == 0){
		print STDERR "Did not find suitable reads for $bam!\n";
		return;
	}
	open(IN, "samtools view $bam |");
	open(OUT, "> $outputfq");
	while(my $line = <IN>){
		my @segs = split(/\t/, $line);
		if(exists($reads->{$segs[0]})){
			chomp $segs[10];
			print OUT "\@$segs[0]\n$segs[9]\n+\n$segs[10]\n";
		}
	}
	close IN;
	close OUT;
}

sub getBreakpointReads{
	my ($bam, $chr, $start, $end) = @_;
	
	print STDERR "samtools view $bam $chr\:$start\-$end |\n";
	open(IN, "samtools view $bam $chr\:$start\-$end |");
	my %reads;
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if($segs[3] < $end && $segs[3] > $start){
			$reads{$segs[0]} = 1;
		}
	}
	close IN;
	return \%reads;
}

sub determineDelBreakpoints{
	my ($chr, $input) = @_;
	open(IN, "grep 'deletion' $input |");
	my @leftstarts;
	my @rightends;
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if($segs[1] eq $chr && $segs[5] eq $chr){
			my $left1 = ($segs[2] < $segs[6])? $segs[2] : $segs[6];
			my $left2 = ($segs[2] > $segs[6])? $segs[2] : $segs[6];
			my $right1 = ($segs[3] < $segs[7])? $segs[3] : $segs[7];
			my $right2 = ($segs[3] > $segs[7])? $segs[3] : $segs[7];
			push(@leftstarts, ($left1, $left2));
			push(@rightends, ($right1, $right2));
		}
	}
	
	@leftstarts = sort{$a <=> $b} @leftstarts;
	@rightends = sort{$a <=> $b} @rightends;
	return $leftstarts[-1], $rightends[0];
}

sub pullDivetSupport{
	my ($chr, $start, $end, $outfile, $input) = @_;
	my $sstart = $start - 10000;
	my $send = $end + 10000;
	open(IN, "< $input");
	open(OUT, "> $outfile");
	while(my $line =<IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if($segs[1] eq $chr){
			if($segs[2] < $send && $segs[3] > $sstart){
				print OUT "$line\n";
				next;
			}
		}
		if($segs[5] eq $chr){
			if($segs[6] < $send && $segs[7] > $sstart){
				print OUT "$line\n";
				next;
			}
		}
	}
	close IN;
	close OUT;
}
