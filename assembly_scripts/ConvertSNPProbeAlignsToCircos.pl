#!/usr/bin/perl
# This is a script designed to process the output of the AlignAndOrderSnpProbes.pl script and turn it into a circos data track

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -s <segs file> -l <contig lengths> -c <conflicts file> -k <karyotype file> -o <output basename>\n";
getopt('scko', \%opts);

unless(defined($opts{'s'}) && defined($opts{'l'}) && defined($opts{'c'}) && defined($opts{'k'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my %chrs; # {chr} = color
# Read in karyotype file and store info
open(my $IN, "< $opts{k}") || die "Could not open karyotype file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\s+/, $line);
	$chrs{$segs[2]} = $segs[6];
}
close $IN;

my %ctglens;
# Read in contig lengths
open(my $IN, "< $opts{l}") || die "Could not open FAI file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\s+/, $line);
	$segs[0] =~ s/[\|\\]/_/g;
	$ctglens{$segs[0]} = $segs[1];
}
close $IN;


my %contigs; #{contig} = [chr, start]
# Process segs file and draw contig positions from start sites
open(my $IN, "< $opts{s}") || die "Could not open segs file!\n";
open(my $CTG, "> $opts{o}.ctg");
open(my $LNK, "> $opts{o}.lnk");
open(my $TXT, "> $opts{o}.ctxt");
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my $ismapped = ($segs[8] > 0)? 1 : 0;
	my $orient = ($segs[6] eq "+")? 1 : -1;
	if(!exists($contigs{$segs[3]}) && $ismapped){
		# easy case -- just plot at the start point to the contig length
		my $end = $ctglens{$segs[3]} + $segs[1];
		print {$CTG} "$segs[0] $segs[1] $end $orient\n";
		$contigs{$segs[3]} = [$segs[0], $segs[1]];
	}elsif(exists($contigs{$segs[3]}) && $ismapped){
		# check if it's on the same chromosome and only plot if it is 1mb distant
		if($contigs{$segs[3]}->[0] eq $segs[0] && $segs[1] - $contigs{$segs[3]}->[1] > 1000000){
			my $start = $contigs{$segs[3]}->[1]; 
			my $end = $start + 100000;
			print {$TXT} "$segs[0] $start $end $segs[3]\n";
			my $lstart = $start + $ctglens{$segs[3]};
			my $estart = $segs[1] + $ctglens{$segs[3]};
			print {$LNK} "$segs[0] $start $lstart $segs[0] $segs[1] $estart color=dgrey\n";
		}
	}
}
close $IN;


# Process conflicts file and add to link data
open(my $IN, "< $opts{c}") || die "Could not open conflicts file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[0] eq "*"){
		next;
	}
	my @fsegs = split(/\:/, $segs[1]);
	if($fsegs[2] < 0.80){
		# at least 20% of the contig maps somewhere else!
		# only plotting the next best hit to lessen confusion
		my @ssegs = split(/\:/, $segs[2]);
		my $ochr = $contigs{$segs[0]}->[0];
		my $start = $contigs{$segs[0]}->[1];
		
		my $schr = ($ochr eq $fsegs[0])? $ssegs[0] : $fsegs[0];
		my $color = ($ochr eq $fsegs[0])? $chrs{$ssegs[0]} : $chrs{$fsegs[0]};
		my $tend = $start + 100000;
		
		my $lstart = $start + $ctglens{$segs[3]};
		my $estart = $segs[1] + $ctglens{$segs[3]};
		# Note: Perhaps I need to run this first so that I can store the location of the second segment?
		# Note: Right now, it's a fixed width plot at the beginning of the secondary chromosome
		print {$TXT} "$ochr $start $tend $segs[0]\n";
		print {$LNK} "$ochr $start $lstart $schr 1 1000000 color=$color\n";
	}
}
close $IN;

close $TXT;
close $LNK;
close $CTG;

# Now attempt to write out a base configuration file
open(my $CONF, "> $opts{o}.conf");
print {$CONF} "<<include colors.conf>>\n";
print {$CONF} "<ideogram>\nshow = yes\nthickness = 20p\nstroke_thickness = 1\nstroke_color = black\n";
print {$CONF} "fill = yes\nradius = 0.90r\nshow_label = yes\nlabel_font = default\nlabel_radius = dims(ideogram,radius_outer) + 100p\n";
print {$CONF} "label_size = 24p\nlabel_parallel = yes\n</ideogram>\n";

print {$CONF} "\n<image>\ndir = .\nfile = $opts{0}.png\npng = yes\nradius = 1500p\nauto_alpha_colors = yes\nauto_alpha_steps = 5\n</image>\n";

print {$CONF} "\nchromosomes_units = 10000000\nchromosomes_display_default = yes\n";

print {$CONF} "\n<plots>\n";

print {$CONF} "\n<plot>\nshow = yes\ntype = histogram\nfile = $opts{o}.ctg\ncolor = black\nthickness = 2\nfill_under = yes\nfill_color = red\n";
print {$CONF} "r0 = 0.80r\nr1 = 0.90r\nmin=-1\nmax=1\n</plot>\n";

print {$CONF} "\n<plot>\nshow = yes\ntype = text\ncolor = black\nfile = $opts{o}.ctxt\nr0 = 0.80r\nr1 = 0.90r\nlabel_size = 20p\n</plot>\n";

print {$CONF} "\n<\plots>\n";

print {$CONF} "\n<links>\n";

print {$CONF} "\n<link>\nfile = $opts{o}.lnk\nradius = 0.80r\nbezier_radius = 0.1r\nthickness = 2\nribbon = yes\n</link>\n";

print {$CONF} "\n<\plots>\n";

