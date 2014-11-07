#!/usr/bin/perl
# This pipeline is designed to take formatted files in a directory (animal name first) and generate major tables 2 and 4
# Files must be named beds with the usual named format (animal_gain/loss_method)
# Input is ARGV with the ls string to be used.

use strict;
#use Statistics::R;

my $input = $ARGV[0];
chomp $input;
# merged_method generation
my $bedtools = "~/bin";
my $outfolder = "merged_method";
my @b_files = `ls $input`;


my %d_animals;
foreach my $f (@b_files){
	my @f_segs = split(/\_/, $f);
	$d_animals{$f_segs[0]} = 1;
}

mkdir ("$outfolder") || print "$!\n";

my @m_files;
my $an_counter = 0;
foreach my $d (sort {$a cmp $b} keys(%d_animals)){
	print "cat $d\n";
	system ("cat $d*.bed | $bedtools/mergeBed -i stdin -nms > $outfolder/$d\_merge_method.bed");
	push(@m_files, "$outfolder/$d\_merge_method.bed");
	$an_counter++;
}
print "Sorting merged method files\n";
my $cat_str;
foreach my $mm (@m_files){
	$cat_str .= "$mm ";
}
my $cat_out = "$outfolder/all_cat.bed";
system("cat $cat_str > $cat_out");
print "Creating major_table_2\n";
my @meth_list = table_2_generate($cat_out, "./major_table_2.tab");
print "Creating major_table_4\n";
my $merge_cat = "$outfolder/all_merged_cat.bed";
system("$bedtools/mergeBed -i $cat_out -nms > $merge_cat");
table_4_generate($merge_cat, "./major_table_4.tab", $an_counter, \@meth_list);

exit;
#Table2 generation
# "inmethod?" will be a "y" if yes and a "-" if no.
# CNVR will be  in the UCSC format
sub table_2_generate (){
my $filename = shift(@_);
my $output = shift(@_);
chomp $filename;
open (IN, "< $filename") || die "Could not open file!\n";

my @animal_table;
my $c = 0;

my %method_possibilites;
my %m_pos;
while (my $line = <IN>){
	chomp $line;
	$line =~ s/\r//g;
	my @seg = split(/\t/, $line);
	my @method = split(/;/, $seg[3]);
	my $annam;
	my %m_array;
	foreach my $meth (@method){
		my @mini_seg = split(/\_/, $meth);
		$annam = $mini_seg[0];
		$m_array{"$mini_seg[1]\_$mini_seg[2]"} = 1;
		$method_possibilites{"$mini_seg[1]\_$mini_seg[2]"} = 1;
		$m_pos{$mini_seg[2]} = 1;
	}
	push (@animal_table, [$annam, $seg[0], $seg[1], $seg[2], \%m_array]);
}
close IN;
my @meth_list = sort{$a cmp $b} keys(%m_pos);
open (OUT, "> $output");
print OUT "animal\tcnvr";
foreach my $m_pos (sort {$a cmp $b} keys(%method_possibilites)){
	print OUT "\t$m_pos";
}
print OUT "\tnum_method\n";
foreach my $l_ref (@animal_table){
	my $num_pos = 0;
	my %rev_hash = %{$l_ref->[4]};
	print OUT "$l_ref->[0]\t$l_ref->[1]\:$l_ref->[2]\-$l_ref->[3]";
	foreach my $n_pos (sort {$a cmp $b} keys(%method_possibilites)){
		if($rev_hash{$n_pos}){
			$num_pos++;
			print OUT "\ty";
		} else {
			print OUT "\t-";
		}
	}
	print OUT "\t$num_pos\n";
}
close OUT;
return @meth_list;
}

#table4 generation
sub table_4_generate (){
my $filename = shift(@_);
my $output = shift(@_);
my $an_counter = shift(@_);
my $m_ref = shift(@_);
chomp $filename;
open (IN, "< $filename") || die "Could not open file!\n";
open (OUT, "> $output");
print OUT "num\tchr\tstart\tend\tucsc\tevent\t\#method\t\#animal\t\%animal\n";
my $counter = 1;
my %meth_venn;
my %an_venn;
my $x = 0;
while (my $line = <IN>){
	chomp $line;
	my @l_seg = split(/\t/, $line);
	my @n_seg = split(/\;/, $l_seg[3]);
	my $m_num = 0;
	my $a_num = 0;
	my $event = "";
	my %method = ();
	my %animal = ();
	my %event = ();
	foreach my $m (@{$m_ref}){
		# Initialize each value at this CNVR position as "false"
		$meth_venn{$m}->[$x] = "FALSE";
		$an_venn{$m}->[$x] = "FALSE";
	}
	foreach my $a (@n_seg){
		my @m_seg = split(/\_/, $a);
		$method{$m_seg[2]} = 1;
		$event{$m_seg[1]} = 1;
		$animal{$m_seg[0]} = 1;
		$meth_venn{$m_seg[2]}->[$x] = "TRUE";
		$an_venn{$m_seg[0]}->[$x] = "TRUE";
	}
	foreach my $k (keys(%method)){
		$m_num += 1;
	}
	foreach my $k (keys(%animal)){
		$a_num += 1;
	}
	if (scalar(keys(%event)) > 1){
		$event = "both";
	} else {
		my @temp = keys(%event);
		$event = $temp[0];
	}
	my $a_perc = ($a_num / $an_counter) * 100;
	print OUT "$counter\t$l_seg[0]\t$l_seg[1]\t$l_seg[2]\t$l_seg[0]\:$l_seg[1]\-$l_seg[2]\t$event\t$m_num\t$a_num\t$a_perc\n";
	$counter++;
	$x++;
}
venn_diagram_creation(\%meth_venn, "method_venn.png");
venn_diagram_creation(\%an_venn, "animal_venn.png");
close IN;
close OUT;
return;
}

sub venn_diagram_creation(){
	require Statistics::R;
	my ($h_ref, $output) = @_;
	my $pwd = `pwd`;
	chomp $pwd;
	if (scalar(keys(%{$h_ref})) < 2 || scalar(keys(%{$h_ref})) > 3){
		print "Cannot create Venn for this dataset\n";
		return;
	}
	my @files;
	my @f_key = sort{$a cmp $b} keys(%{$h_ref});
	my $num_lines = scalar(@{$h_ref->{$f_key[0]}});
	print "Number of CNVRs: $num_lines\n";
	open (OUT, "> rinput.tmp");
	my $header = join(',', @f_key);
	print OUT "$header\n";
	for (my $y = 0; $y < ($num_lines); $y++){
		my $out_string = "";
		foreach my $v (@f_key){
			$out_string .= "$h_ref->{$v}->[$y],";	
		}
		chop $out_string;
		$out_string .= "\n";
		print OUT "$out_string";
		
		
	}
	close OUT;
	# Write an R script to process:
	open (R, "> temp.R");
	print R qq{library(limma)\nv<-read.table("rinput.tmp", sep=',', header=T)\na<-vennCounts(v)\n};
	print R qq{png(file="$output", bg="transparent")\nvennDiagram(a, circle.col= c("blue","green","red"))\ndev.off()\n };
	close R;
	system("R CMD BATCH temp.R");
	print "Venn diagram is in $pwd/$output\n";
	return;
}

=begin
my $R = Statistics::R->new();
	$R->startR();
	$R->send(q/library(limma)/);
	my $rbind = "";
	for (my $x = 0; $x < scalar(@f_key); $x++){
		$R->send(qq{$f_key[$x]<-read.table("$pwd/$files[$x]", header=F)});
		#$R->send(qq{print($f_key[$x])});
		#my $d = $R->read();
		#print "$d\n";
		$rbind .= "$f_key[$x],";
	}
	chop $rbind;
	print $rbind . "\n";
	$R->send(qq/c<-cbind($rbind)/);
	$rbind =~ s/(\w+),*/\"$1\",/g;
	chop $rbind;
	print $rbind . "\n";

	$R->send(q/a<-vennCounts(c)/);
	$R->send(qq/print(a)/);
	my $ret = $R->read();
	print "$ret\n";
	$R->send(qq{pdf(file="$pwd/$output")});
	#$R->send(q{vennDiagram(a)});
	my $stuff = $R->read();
	print "stuff $stuff\n";
	$R->send(q/dev.off()/);
	
	$R->stopR();
=cut