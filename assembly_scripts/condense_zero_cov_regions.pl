#!/usr/bin/perl
#

my $usage = "perl $0 <input bedtools coverage tab> <output zero coverage bed>\n";

chomp(@ARGV);

open(my $IN, "< $ARGV[0]") || die "$usage";
open(my $OUT, "> $ARGV[1]");
my $chr = "NA"; my $start = 0; my $end = 0; my $inzero = 0;

my $h = <$IN>;
chomp $h; my @first = split(/\t/, $h);
$chr = $first[0];
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[0] ne $chr){
		if($inzero){
			print {$OUT} "$chr\t$start\t$end\n";
			$inzero = 0;
			$start = 0; $end = 0;
		}
		$chr = $segs[0];
	}elsif($inzero && $segs[2] != 0){
		print {$OUT} "$chr\t$start\t$end\n";
                $inzero = 0;
                $start = 0; $end = 0;
	}elsif($segs[2] == 0 && !$inzero){
		$inzero = 1;
		$start = $segs[1];
		$end = $segs[1];
	}elsif($inzero){
		$end = $segs[1];
	}
}
close $IN;

if($inzero){
	print {$OUT} "$chr\t$start\t$end\n";
}
close $OUT;
exit;
