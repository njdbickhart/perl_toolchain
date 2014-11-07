#!/usr/bin/perl
# This script is designed to process named bed files (typical format: animal_gain/loss_method) in order to produce two tables
# two venn diagrams and a bar contribution chart
# The bed files must have the animal name first in the filename


use strict;
use threads;
use Spreadsheet::WriteExcel;

my $bedtools = "~/bin";

my $usage ="\nUsage: $0 \'lscommand\'\n";

my $ls_command = $ARGV[0];
chomp $ls_command;
my @ls_files = `ls $ls_command`;
if (scalar(@ls_files) < 2 || !defined($ARGV[0])){
	print $usage;
	exit;
}

print "Processing: " . scalar(@ls_files) . " files\n";

#create_table_3(@ls_files, $bedtools);
#exit;
my $t_three;
$t_three = threads->create(\&create_table_3, @ls_files, $bedtools);

my %super_hash;
# {method} {animal} [line] [chr, start, end]
my %an_hash;
my %d_animals; 
foreach my $n (@ls_files){
	my @f_segs = split(/\_/, $n);
	$d_animals{$f_segs[0]} = 1;
}

my @tot_animals = sort {$a cmp $b} keys(%d_animals);
print "Total animals: " . scalar(@tot_animals) . "\n";


my $temp_a;
my $animal;
my @a_animal;
my $method;
my %unique_method;
foreach my $file (@ls_files){
	($temp_a, $animal, $method) = put_file_into_hasharray($file);
	push(@{$super_hash{$animal}->{$method}}, @{$temp_a});
	# Now super_hash{animal} contains a hash reference
	# super_hash{animal}->{method}->[chr,start,end,gl]
	$unique_method{$method} = 1;
}
my @meth_types = sort {$a cmp $b} keys(%unique_method);
foreach my $an_pass (keys(%super_hash)){
	foreach my $mth (@meth_types){
		if (exists($super_hash{$an_pass}->{$mth})){
			next;
		} else {
			$super_hash{$an_pass}->{$mth}->[0] = [0,0,0];
		}
	}
}

my ($table_array) = determine_unique(\%super_hash, \@tot_animals);
# should contain a line per animal, with a cnvr per method (or blank if no overlapping cnvrs)
### work on this, print out array!

open (OUT, "> major_table_1.tab");
print OUT "animal";
foreach my $j (@meth_types){
	print OUT "\t$j";
}
print OUT "\n";

foreach my $v (@{$table_array}){
	print OUT $v;
}
close OUT;
$t_three->join();
exit;

#########################
# 	Subroutines	#
#########################

sub put_file_into_hasharray (){
	# filename
	# returns reference of hash
	my ($file) = @_;
	chomp $file;
	print "Putting $file into memory...\n";
	my @r_array;
	my %an_hash;
	my $animal;
	my $method;
	open (IN, "< $file") || die "could not open file $file\n";
	while (my $line = <IN>){
		chomp $line;
		$line =~ s/\r//g;
		my @seg = split(/\t/, $line);
		my @a_seg = split(/\_/, $seg[3]);
		chomp $a_seg[2];
		$animal = $a_seg[0];
		my $val = $a_seg[2];
		$method = $val;
		if ($val eq ""){ print "error reading file\n"; $val = 0;}
		my $gl = 0;
		if($a_seg[1] eq "gain"){
			$gl = 1;
		}
		push(@r_array, [$seg[0], $seg[1], $seg[2], $gl]);
	}
	close IN;
	return (\@r_array, $animal, $method);
}

sub determine_unique (){
	# super_hash_reference, tot_animals
	# returns reference of table to main sequence
	my ($s_tab, $t_an) = @_;
	# All of these are references to hashes per animal
	my @real_table;
	my ($an_array, $use, $meth_types);
	foreach my $an (@{$t_an}){
		print "Starting on $an\n";
		($an_array, $use) = animal_comp($s_tab->{$an}, $an);
		unless ($use == 0){
			push(@real_table, @{$an_array});
		}
	}
	return \@real_table;	
}

sub animal_comp (){
	# should be arrays of values for each dataset per animal
	my ($method_hash, $an) = @_;
	my $temp_coords;
	my %method_coords;
	my @done_array;
	my $good = 0;
	my @meth_types;
	for (my $x = 1; $x < 30; $x++){
		my $chr = "chr$x";
		foreach my $meth (sort {$a cmp $b} keys(%{$method_hash})){
			$temp_coords = chr_into_array($method_hash->{$meth}, $chr);
			$method_coords{$meth} = $temp_coords;
			
		}
		
		my ($temp_array, $use) = interval_comp(\%method_coords, $chr, $an); 
		if ($use == 0){
			next;
		} else{
			push(@done_array, @{$temp_array});
			$good = 1;
		}
	}
	return \@done_array, $good;
		
}

sub interval_comp (){
	# arrays of values, but just with start and end values.
	my ($method_coords, $chr, $an) = @_;
	my $hit = 0; 
	my @temp;
	my @final_array;
	my $t_value;
	my $new_x = 0;
	my $old_x = 0;
	my $use = 0;
	my $ovalue;
	my $value;
	
	print "Cycling through $chr\n";
	for (my $x = 1; $x < 200000000; $x += 2500){
		$old_x = $x + 50000;
		$value = $old_x;
		my $gl = 3;
		foreach my $key (sort {$a cmp $b} keys (%{$method_coords})){
			($t_value, $hit, $new_x, $value, $gl) = splice_and_return($method_coords->{$key}, $x, $chr, $old_x, $hit, $ovalue);
			push(@temp, $t_value);
			if($ovalue != $value && $value != $old_x){$ovalue = $value;}
			if ($new_x > $old_x){$old_x = $new_x;}
		}
		if ($x + 50000 != $old_x){$x = $old_x - 2499;}
		if ($hit > 0){
			$use = 1; 
			my $t_string = "$an\t";
			for (my $y = 0; $y < scalar(@temp); $y++){
				$t_string = $t_string . $temp[$y];
			}
			$t_string = $t_string . "\n";
			#print $t_string;
			push(@final_array, $t_string);
		}
		@temp = ();
		$hit = 0;
	}
	return \@final_array, $use; 
}
sub chr_into_array (){
	# Just to reduce the number of for loops in animal_comp
	my ($a_ref, $chr) = @_;
	my @return;
	
	foreach my $a (@{$a_ref}){
		if ($a->[0] eq $chr){
			
			push(@return, [$a->[1], $a->[2], $a->[3]]);
		}
	}

	return \@return;
}
sub splice_and_return (){
	# return a_ref string hit new_x
	my ($a_ref, $x, $chr, $old_x, $hit, $spl, $gl) = @_;
	my $new_x = $x;
	my $string = "-\t"; 
	my $last_values = $old_x;
	my %type_vals = (
		0 => "loss",
		1 => "gain",
		3 => "undef");
	if (scalar(@{$a_ref}) == 0){
		return ($string, $hit, $new_x, $last_values, $gl);
	}
	for (my $y = 0; $y < scalar(@{$a_ref}); $y++){
		if ($old_x >= $a_ref->[$y]->[0] && $x <= $a_ref->[$y]->[1]){
			my $len = $a_ref->[$y]->[1] - $a_ref->[$y]->[0];
			my $half = $a_ref->[$y]->[1] - ($len / 2);
			if($a_ref->[$y]->[1] == $spl){last;}
			$hit = 1;
			if($gl = 3){
				$gl = $a_ref->[$y]->[2];
			}elsif($gl != $a_ref->[$y]->[2]){
				open (LOG, ">> both_intervals.log");
				print LOG "$chr\:$a_ref->[$y]->[0]\-$a_ref->[$y]->[1] $gl\n";
				close LOG;
			}
			my $type = $type_vals{$gl};
			$string = "$chr\:$a_ref->[$y]->[0]\-$a_ref->[$y]->[1] ($type)\t";
			if($old_x < $half){
				$new_x = $half;
				$last_values = $a_ref->[$y]->[1];
			}
			
			splice(@{$a_ref}, $y, 1);
			last;
		} elsif ($x > $a_ref->[$y]->[0] && $x > $a_ref->[$y]->[1]){
			last;
		}
	}
	return ($string, $hit, $new_x, $last_values, $gl);
}

sub create_table_3 (){
	require Spreadsheet::WriteExcel;
	my $bedtools = pop(@_);
	my @filelist = @_;
	my %cat_list;
	my %m_hash;
	foreach my $file (@filelist){
		chomp $file;
		my @f_segs = split(/\_/, $file);
		push(@{$cat_list{$f_segs[1]}}, $file);
		
	}
	foreach my $meth (sort {$a cmp $b} keys(%cat_list)){
		my $catf = "cat ";
		foreach my $q (@{$cat_list{$meth}}){
			$catf .= "$q ";
		}
		print "table3: $catf | $bedtools/mergeBed -i stdin -nms > temp\n";
		system("$catf | $bedtools/mergeBed -i stdin -nms > temp");
		open(IN, "< temp");
		#print "Putting $file into memory...\n";
		while (my $line = <IN>){
			chomp $line;
			$line =~ s/\r//g;
			my @seg = split(/\t/, $line);
			my @aseg = split(/\;/, $seg[3]);
			my @an_tag;
			my $gl = "a";
			foreach my $tag (@aseg){
				my @t_segs = split(/\_/, $tag);
				push(@an_tag, $t_segs[0]);
				if ($gl eq "a"){
					$gl = $t_segs[1];
				} elsif ($gl ne $t_segs[1]){
					$gl = "both";
				}
			}
			my @push_array = sort{$a cmp $b} (@an_tag);
			push (@{$m_hash{$meth}}, [$seg[0], $seg[1], $seg[2], $gl, @push_array]);
		}
		close IN;
	}
	system("rm temp");
	my $workbook = Spreadsheet::WriteExcel->new('major_table_3.xls');
	#open(PIE, "> pie_chart_info.tab");
	my $counter = 1;
	foreach my $key (sort {$a cmp $b} (keys(%m_hash))){
		my $acounter = 0;
		print "table3 $key\n";
		#print OUT "$key\n";
		my $worksheet = $workbook->add_worksheet($key);
		$worksheet->write(0,0, "CNVR #");
		$worksheet->write(0,1, "chr");
		$worksheet->write(0,2, "start");
		$worksheet->write(0,3, "end");
		$worksheet->write(0,4, "length");
		$worksheet->write(0,5, "g/l");
		$worksheet->write(0,6, "# animals");
		$worksheet->write(0,7, "animals");
		my $row = 1;
		#print OUT "CNVR #\tchr\tstart\tend\t# animals\tanimals\n";
		foreach my $l (@{$m_hash{$key}}){
			my $chr = shift(@{$l});
			my $start = shift(@{$l});
			my $end = shift(@{$l});
			my $g_l = shift(@{$l});
			my $len = $end - $start;
			$acounter = scalar(@{$l});
			$worksheet->write($row, 0, $counter);
			$worksheet->write($row, 1, $chr);
			$worksheet->write($row, 2, $start);
			$worksheet->write($row, 3, $end);
			$worksheet->write($row, 4, $len);
			$worksheet->write($row, 5, $g_l);
			$worksheet->write($row, 6, $acounter);
			#print OUT "$counter\t$chr\t$start\t$end\t$acounter";
			my $col = 7;
			foreach my $b (@{$l}){
				$worksheet->write($row, $col, $b);
				$col++;
				#print OUT "\t$b";
			}
			#print OUT "\n";
			$row++;
			$counter++;
		}
		#print OUT "\n\n";
	}
	
	#close PIE;
}