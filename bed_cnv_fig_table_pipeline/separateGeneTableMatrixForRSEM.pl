#!/usr/bin/perl
# This script reads in a gene matrix file for RSEM, then selects only the samples that are needed for the analysis based on user input
# Finally, the script orders the samples based on the column separated conditions added by the user.

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -f input_gene_matrix -t comma_separated_order_ofconditions -s subtype_selected -k key_file -o output_matrix_file\n";

getopt('ftkos', \%opts);

unless(defined($opts{'f'}) && defined($opts{'t'}) && defined($opts{'k'}) && defined($opts{'o'}) && defined($opts{'s'})){
	print $usage;
	exit;
}

my @order = split(/,/, $opts{'t'}); 

# Key file: condition \t subtype \t colname
# Opening the key file first
my %goodkeys;
open(IN, "< $opts{k}") || die "Could not open input key file!\n";
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[1] eq $opts{'s'}){
		$goodkeys{$segs[2]} = $segs[0];
	}
}
close IN;

if(scalar(keys(%goodkeys)) < 1){
	print "Error! Found no keys matching the requested subtype!\n";
	exit;
}

# It's really sloppy, but this is just going to be a brute-force memorization and print out
open(IN, "< $opts{f}") || die "Could not open input file!\n";
my @reserve;

my $h = <IN>;
chomp $h;
my @headers = split(/\t/, $h);
for(my $x = 1; $x < scalar(@headers); $x++){
	if(exists($goodkeys{$headers[$x]})){
		# I'm saving the index of the column for later retrieval
		push(@reserve, $x);
	}
}

my %data; # {condition}->{colname} -> [entries]
my @rows;
# Now to process the rest of the matrix
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	push(@rows, $segs[0]);
	foreach my $idx (@reserve){
		my $col = $headers[$idx];
		my $val = $segs[$idx];
		my $cond = $goodkeys{$col};
		push(@{$data{$cond}->{$col}}, $val);
	}
}
close IN;

# Now for the data sorting. Going to use a presorted hash array here
my %finalorder;
foreach my $o (@order){
	$finalorder{$o} = [sort{$a cmp $b} keys(%{$data{$o}})];
}

open(OUT, "> $opts{o}");
# Printing out the header first
foreach my $o (@order){
	foreach my $c (@{$finalorder{$o}}){
		print OUT "\t$c";
	}
}
print OUT "\n";

# Now the data matrix
for(my $x = 0; $x < scalar(@rows); $x++){
	print OUT "$rows[$x]";
	foreach my $o(@order){
		foreach my $c (@{$finalorder{$o}}){
			print OUT "\t" . $data{$o}->{$c}->[$x];
		}
	}
	print OUT "\n";
}
close OUT;
exit;
