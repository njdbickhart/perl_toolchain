#!/usr/bin/perl
# This script is designed to take multiple bed files (from ARGV) and compare them
# My Eventual plan is to try to make a plot using the GD library, but I'll be happy with an algorithm that works for a change!
# Files must have unique start names separated by underscores (eg. "btan_ngs_reads.bed");
# Requires intersectBed, mergeBed and subtractBed in the bin or path
# Can only handle 2, 3, 4 and 5 way comparisons. Anything higher is too arduous to try to make

use strict;
use GD::Simple;

my $usage = "$0 \<output\> \<file1\> \<file2\> ...\n\n";
my $outfile = shift(@ARGV);
my @files = @ARGV;
my $f_num = scalar(@files);
if ($f_num < 2 || $f_num > 5){
	print $usage;
	exit;
}

chomp $files[-1];
my %data_store;
my @prefixes;
my %prefix_file;
my %prefix_index;
my %file_prefix;
# I will try to store information using this scheme: $data_store{prefix} = @{counta, countb, countc ... bases}
# generate file prefixes to keep track of them
for (my $x = 0; $x < scalar(@files); $x ++){
	my @pre = split(/\_/, $files[$x]);
	push(@prefixes, $pre[0]);
	$prefix_file{$pre[0]} = $files[$x];
	$prefix_index{$files[$x]} = $x;
	$file_prefix{$files[$x]} = $pre[0];
	my $cy = $x + 1;
	print "File $cy is $files[$x]\n";
}

# Lets try generating all of the permutations from the prefixes and calculating how many keys there are
# permutations should be 2^x - 1 where x is the number of files
my $num_perm = 2^$f_num - 1;
# each element's permutations should be 2^(x - 1), where x is the number of files
my $list_perm = 2^($f_num -1);

my @permutations;
my %unique_perm;

for (my $x = 0; $x < $f_num; $x++){
	# Separate comparison prefix from the others
	my $starter = $prefixes[$x];
	my @others;
	for (my $y = 0; $y < $f_num; $y++){
		if ($y == $x){
			next;
		} else {
			push(@others, $prefixes[$y]);
		}
	}
	# Now, take into account the venn possibilities
	# Single comparisons
	$unique_perm{$starter} = 1;	
	# Two way comparisons
	foreach my $o (@others){
		if (!defined($o) || $o eq ''){next;}
		my @temp_sorter = ($starter);
		push (@temp_sorter, $o);
		@temp_sorter = sort {$a cmp $b} @temp_sorter;
		my $string = join(';', @temp_sorter);
		$unique_perm{$string} = 1;
	}
	if ($f_num == 2){next;}
	# Three way comparisons
	my $three = 0;
	my $rem = 0;
	if ($f_num == 3){ $three = 1;}elsif($f_num == 4){ $three =  3; $rem = 1;}elsif($f_num == 5){$three = 20; $rem = 2;}
	for (my $y = 0; $y < $three; $y++){
		my @temp_sorter = ($starter);
		my @holder;
		for (my $z = 0; $z < $rem; $z++){
			if ($z == 1){
				push(@holder, shift(@others));
			}else{
				push(@holder, pop(@others));
			}
		}
			
		push(@temp_sorter, @others);
		@temp_sorter = sort{$a cmp $b} @temp_sorter;
		my $string = join(';', @temp_sorter);
		$unique_perm{$string} = 1;
		for (my $z = 0; $z < $rem; $z++){
			unshift(@others, $holder[$z]);
		}
		if ($y == 4){
			my $xt = shift(@others);
			push(@others, $xt);
		}
   		if ($y % 5 == 0){
			randomize_array(\@others);
		}
	}
	# Four way comparisons
	my $four = 0;
	my $rem = 0;
	if ($f_num == 4){ $four = 1;}elsif($f_num == 5){$four = 20; $rem = 1;}
	for (my $y = 0; $y < $three; $y++){
		my @temp_sorter = ($starter);
		my @holder;
		for (my $z = 0; $z < $rem; $z++){
			if($z == 1){
				push(@holder, shift(@others));
			} else{
				push(@holder, pop(@others));
			}
		}
		push(@temp_sorter, @others);
		@temp_sorter = sort{$a cmp $b} @temp_sorter;
		my $string = join(';', @temp_sorter);
		$unique_perm{$string} = 1;
		for (my $z = 0; $z < $rem; $z++){
			unshift(@others, $holder[$z]);
		}
		if ($y % 5 == 0){
			randomize_array(\@others);
		}
	}
	# Five way comparisons
	if ($f_num == 5){
		my @temp_sorter = ($starter);
		push(@temp_sorter, @others);
		@temp_sorter = sort{$a cmp $b} @temp_sorter;
		my $string = join(';', @temp_sorter);
		$unique_perm{$string} = 1;
	}
}
# Initialize novel comparison holders in the datastore
# I will sort prefixes lexically
foreach my $p (keys(%unique_perm)){
	for (my $x = 0; $x < scalar(@prefixes) * 2; $x++){
		$data_store{$p}->[$x] = 0;
	}
}

# Now to do the harder multi-way comparisons
my $ne_lengths = 0;
foreach my $pkey (sort {$a cmp $b} keys(%unique_perm)){
	my @temp_holder = @files;
	my @samples = split(/\;/, $pkey);
	# Novel comp
	if (scalar(@samples) == 1){
		my $comp = $prefix_file{$pkey};
		my @novel = grep($_ ne $comp, @temp_holder);
		unshift(@novel, $comp);
		my ($count, $length) = novel_bed_comp(@novel);
		if ($length == undef){$length = 0;}
		my $p = $prefix_file{$pkey};
		$data_store{$pkey}->[$prefix_index{$p}] = $count;
		$data_store{$pkey}->[$prefix_index{$p} + $f_num] = $length;
	}
	# Two way comp
	if (scalar(@samples) == 2){
		my $file1 = $prefix_file{$samples[0]};
		my $file2 = $prefix_file{$samples[1]};
		my @two_way = grep($_ ne $file1 && $_ ne $file2, @temp_holder);
		unshift(@two_way, $file2);
		unshift(@two_way, $file1);
		my $rem = 0;
		if ($f_num == 3){ $rem = 1;}elsif($f_num == 4){ $rem = 2;}
		my ($count, $length) = multi_way_comp($rem, \%file_prefix, $pkey, @two_way);
		if ($length == undef){$length = 0;}
		$data_store{$pkey}->[$prefix_index{$file1}] = $count;
		$data_store{$pkey}->[$prefix_index{$file1} + $f_num] = $length;
		shift(@two_way);
		shift(@two_way);
		unshift(@two_way, $file1);
		unshift(@two_way, $file2);
		my ($count, $length) = multi_way_comp($rem, \%file_prefix, $pkey, @two_way);
		if ($length == undef){$length = 0;}
		#if ($data_store{$pkey}->[$base_holder] != $length){ print "$pkey lengths not equal\n"; $ne_lengths++;}
		$data_store{$pkey}->[$prefix_index{$file2}] = $count;
		$data_store{$pkey}->[$prefix_index{$file2} + $f_num] = $length;
	}
	# Three way comp
	if (scalar(@samples) == 3){
		my $file1 = $prefix_file{$samples[0]};
		my $file2 = $prefix_file{$samples[1]};
		my $file3 = $prefix_file{$samples[2]};
		my @two_way = grep($_ ne $file1 && $_ ne $file2 && $_ ne $file3, @temp_holder);
		unshift(@two_way, $file3);
		unshift(@two_way, $file2);
		unshift(@two_way, $file1);
		my $rem = 0;
		if ($f_num == 4){ $rem = 1;}elsif($f_num == 5){$rem = 2;}
		my ($count, $length) = multi_way_comp($rem, \%file_prefix, $pkey, @two_way);
		if ($length == undef){$length = 0;}
		$data_store{$pkey}->[$prefix_index{$file1}] = $count;
		$data_store{$pkey}->[$prefix_index{$file1} + $f_num] = $length;
		shift(@two_way);
		shift(@two_way);
		shift(@two_way);
		unshift(@two_way, $file3);
		unshift(@two_way, $file1);
		unshift(@two_way, $file2);
		my ($count, $length) = multi_way_comp($rem, \%file_prefix, $pkey, @two_way);
		if ($length == undef){$length = 0;}
		#if ($data_store{$pkey}->[$base_holder] != $length){ print "$pkey lengths not equal\torig: $data_store{$pkey}->[$base_holder]\t new: $length \n"; $ne_lengths++;}
		$data_store{$pkey}->[$prefix_index{$file2}] = $count;
		$data_store{$pkey}->[$prefix_index{$file2}+ $f_num] = $length;
		shift(@two_way);
		shift(@two_way);
		shift(@two_way);
		unshift(@two_way, $file2);
		unshift(@two_way, $file1);
		unshift(@two_way, $file3);
		my ($count, $length) = multi_way_comp($rem, \%file_prefix, $pkey, @two_way);
		if ($length == undef){$length = 0;}
		#if ($data_store{$pkey}->[$base_holder] != $length){ print "$pkey lengths not equal\torig: $data_store{$pkey}->[$base_holder]\t new: $length \n"; $ne_lengths++;}
		$data_store{$pkey}->[$prefix_index{$file3}] = $count;
		$data_store{$pkey}->[$prefix_index{$file3}+ $f_num] = $length;
	}
	# Four way comp
	# Since I will max out at four comparisons for now, lets just make this simple
	if (scalar(@samples) == 4){
		foreach my $sam (@samples){
			my @four_way = grep($_ ne $prefix_file{$sam}, @temp_holder);
			unshift(@four_way, $prefix_file{$sam});	
			my $rem = 0;
			if ($f_num == 5){ #find the value that does not belong in the comparison and remove it
				$rem = 1;
				my $outcast;
				foreach my $p (@prefixes){
					my $found = 0;
					foreach my $s (@samples){
						if ($s eq $p){
							$found++;
						}
					}
					if($found = 0){
						$outcast = $p;
						last;
					}
				}
				my @t_four_way = grep ($_ ne $file_prefix{$outcast}, @four_way);
				@four_way = (@t_four_way, $file_prefix{$outcast}); #place the file to be removed at the end.
			}
			my ($count, $length) = multi_way_comp($rem, \%file_prefix, $pkey, @four_way);
			if ($length == undef){$length = 0;}
			my $p = $prefix_file{$sam};
			$data_store{$pkey}->[$prefix_index{$p}] = $count;
			$data_store{$pkey}->[$prefix_index{$p} + $f_num] = $length;
		}
	}
	# Five way comp
	# This is pretty daring, but lets see if it works
	if (scalar(@samples) == 5){
		foreach my $sam (@samples){
			my @five_way = grep($_ ne $prefix_file{$sam}, @temp_holder);
			unshift(@five_way, $prefix_file{$sam});
			my $rem = 0;
			my ($count, $length) = multi_way_comp($rem, \%file_prefix, $pkey, @five_way);
			if ($length == undef){$length = 0;}
			my $p = $prefix_file{$sam};
			$data_store{$pkey}->[$prefix_index{$p}] = $count;
			$data_store{$pkey}->[$prefix_index{$p} + $f_num] = $length;
		}
	}
}
foreach my $key (sort {$a cmp $b} keys %data_store){
	my $vals = join("\t", @{$data_store{$key}});
	print "$key\t$vals\n";
}
print_out_png_gd(\%data_store, \@prefixes, $outfile);
system("rm crtemp ctemp mtemp rtemp srtemp stemp tcat");
exit;

sub print_out_png_gd (){
	require GD::Simple;
	my ($h_key, $p_array, $out) = @_;
	my @prefixes = sort {$a cmp $b}@{$p_array};
	my $num_comp = scalar(@prefixes);
	my %data_store = %{$h_key};
	if ($num_comp == 2){
		my $img = GD::Simple->new(434, 306);
		$img->bgcolor('white');
		$img->moveTo(145, 155);
		$img->bgcolor(undef);
		$img->fgcolor('red');
		$img->ellipse(275, 275);
		
		$img->moveTo(280, 155);
		$img->bgcolor(undef);
		$img->fgcolor('blue');
		$img->ellipse(275, 275);
		$img->fgcolor('black');
		$img->fontsize(20);
		$img->moveTo(30,30);
		$img->string("$prefixes[0]");
		$img->moveTo(380, 30);
		$img->string("$prefixes[1]");
		foreach my $keys (keys(%data_store)){
			my @vals = @{$data_store{$keys}};
			my @ints = @vals[0,1];
			my @len = @vals[2,3];
			if($keys eq $prefixes[0]){
				my $int_string = join(' , ', @ints);
				$img->moveTo(65, 135);
				$img->string("$int_string");
				@len = sort {$b <=> $a} @len;
				$img->moveTo(65, 150);
				$img->string("$len[0]");
			} elsif ($keys eq $prefixes[1]){
				my $int_string = join(' , ', @ints);
				$img->moveTo(320, 135);
				$img->string("$int_string");
				@len = sort {$b <=> $a} @len;
				$img->moveTo(320, 150);
				$img->string("$len[0]");
			} else {
				my $int_string = join(' , ', @ints);
				$img->moveTo(185, 135);
				$img->string("$int_string");
				@len = sort {$b <=> $a} @len;
				$img->moveTo(185, 150);
				$img->string("$len[0]");
			}
		}
		open(OUT, "> $out");
		print OUT $img->png;		
	}
	if ($num_comp == 3){
		my @d_keys = sort{$a cmp $b} keys(%data_store);
		my @positions = (220, 90, 
			125, 210,
			215, 240,
			340, 210,
			75, 360,
			215, 385,
			400, 360);			
		my $img = GD::Simple->new(500, 500);
		$img->bgcolor('white');
		$img->moveTo(250, 175);
		$img->bgcolor(undef);
		$img->fgcolor('red');
		$img->ellipse(350, 350);
		
		$img->moveTo(175, 325);
		$img->bgcolor(undef);
		$img->fgcolor('blue');
		$img->ellipse(350, 350);
		
		$img->moveTo(325, 325);
		$img->bgcolor(undef);
		$img->fgcolor('green');
		$img->ellipse(350, 350);
		
		$img->fgcolor('black');
		$img->fontsize(20);
		$img->moveTo(120,25);
		$img->string("$prefixes[0]");
		$img->moveTo(30, 478);
		$img->string("$prefixes[1]");
		$img->moveTo(455, 478);
		$img->string("$prefixes[2]");
		my $y = 0;
		for (my $x = 0; $x < scalar(@d_keys); $x++){
			my $xval = $positions[$y];
			my $yval = $positions[$y+1];
			my $key = $d_keys[$x];
			my @vals = @{$data_store{$key}};
			my @ints = @vals[0,1,2];
			my @len = @vals[3,4,5];
			@len = sort{$b <=> $a} @len;
			my $int_string = join(' , ', @ints);
			$img->moveTo($xval, $yval);
			$img->string("$int_string");
			$yval += 21;
			$xval += 7;
			$img->moveTo($xval, $yval);
			my $len_string;
			if ($len[0] > 100000){
				$len_string = sprintf("%.2f Mb", $len[0] / 1000000);
			}elsif($len[0] == 0){
				$len_string = "0 bp";
			}else{
				$len_string = sprintf("%.2f Kb", $len[0] / 1000);
			} 
			$img->string("$len_string");
			$y += 2;
		}
		open(OUT, "> $out");
		print OUT $img->png;
	}
	if ($num_comp == 4){
		my @d_keys = sort{$a cmp $b} keys(%data_store);
		my @positions = (
			445, 45, #blue
			445, 330, #blue red
			445, 225, #blue red green
			310, 225, #blue red green black
			310, 290, #blue red black
			445, 165, #blue green
			310, 165, #blue green black
			310, 115, #blue black
			100, 330, #red
			100, 225, #red green
			220, 225, #red green black
			240, 290, #red black
			100, 165, #green
			220, 165, #green black
			240, 115  #black
			);
		my $img = GD::Simple->new(600, 400);
		$img->bgcolor('white');
		$img->fgcolor('red');
		$img->bgcolor(undef);
		$img->penSize(2,2);
		$img->rectangle(1, 200, 598, 398);
		$img->fgcolor('blue');
		$img->rectangle(300, 1, 597, 397);

		$img->moveTo(300, 200);
		$img->fgcolor('black');
		$img->ellipse(280,280);

		$img->fgcolor('green');
		$img->moveTo(150,200);
		$img->arc(225, 225, 40, 320);
		$img->moveTo(450, 200);
		$img->arc(225, 225, 220, 140);
	
		$img->moveTo(300, 50);
		$img->arc(204, 204, 50, 130);

		$img->moveTo(300, 350);
		$img->arc(204, 204, 230, 310);
		
		$img->moveTo(2, 15);
		$img->fgcolor('blue');
		$img->string("$prefixes[0]");
		$img->moveTo(2, 35);
		$img->fgcolor('red');
		$img->string("$prefixes[1]");
		$img->moveTo(2, 55);
		$img->fgcolor('green');
		$img->string("$prefixes[2]");
		$img->moveTo(2, 75);
		$img->fgcolor('black');
		$img->string("$prefixes[3]");
		my $y = 0;
		for (my $x = 0; $x < scalar(@d_keys); $x++){
			my $xval = $positions[$y];
			my $yval = $positions[$y+1];
			my $key = $d_keys[$x];
			my @vals = @{$data_store{$key}};
			my @ints = @vals[0,1,2,3];
			my @len = @vals[4,5,6,7];
			@len = sort{$b <=> $a} @len;
			my $int_string = join(',', @ints);
			$img->moveTo($xval, $yval);
			$img->string("$int_string");
			$yval += 15;
			$xval += 7;
			$img->moveTo($xval, $yval);
			my $len_string;
			if ($len[0] > 100000){
				$len_string = sprintf("%.2f Mb", $len[0] / 1000000);
			}elsif($len[0] == 0){
				$len_string = "0 bp";
			}else{
				$len_string = sprintf("%.2f Kb", $len[0] / 1000);
			} 
			$img->string("$len_string");
			$y += 2;
		}
		open(OUT, "> $out");
		print OUT $img->png;
	}
	if ($num_comp == 5){
		my @d_keys = sort{$a cmp $b} keys(%data_store);
		my @positions = (
			640, 50, #blue
			640, 550, #blue red
			680, 396, #blue red green
			578, 380, #blue red green black
			495, 320, #blue red green black orange
			670, 320, #blue red green orange
			527, 420, #blue red black
			460, 412, #blue red black orange
			460, 535, #blue red orange
			680, 187, #blue green
			578, 195, #blue green black
			495, 256, #blue green black orange
			670, 256, #blue green orange
			527, 154, #blue black
			460, 179, #blue black orange
			460, 50, #blue orange
			150, 550, #red
			150, 396, #red green
			290, 380, #red green black
			340, 320, #red green black orange
			180, 320, #red green orange
			330, 420, #red black
			396, 412, #red black orange
			396, 535, #red orange
			150, 187, #green
			280, 195, #green black
			340, 256, #green black orange
			180, 256, #green orange
			320, 154,  #black
			396, 179, #black orange
			396, 50 #orange
			);
		my $img = GD::Simple->new(900, 600);
		$img->bgcolor('white');
		$img->fgcolor('red');
		$img->bgcolor(undef);
		$img->penSize(3,3);
		$img->rectangle(1, 300, 897, 597);
		$img->fgcolor('blue');
		$img->rectangle(450, 1, 896, 596);

		$img->penSize(2,2);
		$img->moveTo(450, 300);
		$img->fgcolor('black');
		$img->ellipse(420,420);

		$img->fgcolor('green');
		$img->moveTo(225,300);
		$img->arc(300, 300, 39, 321);
		$img->moveTo(675, 300);
		$img->arc(300, 300, 220, 140);
		$img->moveTo(450, 75);
		$img->arc(342, 342, 50, 130);
		$img->moveTo(450, 525);
		$img->arc(342, 342, 230, 310);

		$img->fgcolor('orange');
		$img->moveTo(450, 103);
		$img->arc(175,175, 130, 50);
		$img->moveTo(240, 300);
		$img->arc(175,175, 40, 320);

		$img->moveTo(450, 497);
		$img->arc(175,175, 310, 230);
		$img->moveTo(660, 300);
		$img->arc(175,175, 220, 140);

		$img->moveTo(300, 152);
		$img->arc(188,188,9, 86);

		$img->moveTo(300, 448);
		$img->arc(188,188,274, 350);

		$img->moveTo(600, 152);
		$img->arc(188,188,94, 170);

		$img->moveTo(600, 448);
		$img->arc(188,188,190, 267);
		$img->moveTo(2, 15);
		$img->fgcolor('blue');
		$img->string("$prefixes[0]");
		$img->moveTo(2, 35);
		$img->fgcolor('red');
		$img->string("$prefixes[1]");
		$img->moveTo(2, 55);
		$img->fgcolor('green');
		$img->string("$prefixes[2]");
		$img->moveTo(2, 75);
		$img->fgcolor('black');
		$img->string("$prefixes[3]");
		$img->moveTo(2, 95);
		$img->fgcolor('orange');
		$img->string("$prefixes[4]");
		$img->fgcolor('black');
		my $y = 0;
		for (my $x = 0; $x < scalar(@d_keys); $x++){
			my $xval = $positions[$y];
			my $yval = $positions[$y+1];
			my $key = $d_keys[$x];
			my @vals = @{$data_store{$key}};
			my @ints = @vals[0,1,2,3,4];
			my @len = @vals[5,6,7,8,9];
			@len = sort{$b <=> $a} @len;
			my $int_string;
			my @temp_int = sort{$b <=> $a} @ints;
			if ($temp_int[-1] == 0 || $temp_int[-1] == undef){
				if ($temp_int[1] == 0){
					$int_string = "$temp_int[0]";
				}else{
					my $hold;
					for (my $y = 0; $y < scalar(@temp_int); $y++){
						if ($temp_int[$y] != 0 && $temp_int[$y] != undef){
							$hold = $temp_int[$y];
						}
						if($temp_int[$y] == 0 || $temp_int[$y] == undef){
							if ($temp_int[0] == $hold){
								$int_string = "$temp_int[0]";
								last;
							}
							$int_string = "$temp_int[0]\-$hold";
							last;
						}
					}
				}
			}elsif ($temp_int[0] == 0){ $int_string = "$temp_int[0]";
			}else{
				$int_string = "$temp_int[0]\-$temp_int[-1]";
			}
			$img->moveTo($xval, $yval);
			$img->string("$int_string");
			$yval += 15;
			$xval -= 1;
			$img->moveTo($xval, $yval);
			my $len_string;
			if ($len[0] > 100000){
				$len_string = sprintf("%.2f Mb", $len[0] / 1000000);
			}elsif($len[0] == 0){
				$len_string = "0 bp";
			}else{
				$len_string = sprintf("%.2f Kb", $len[0] / 1000);
			} 
			$img->string("$len_string");
			$y += 2;
		}
		open(OUT, "> $out");
		print OUT $img->png;
	}
	close OUT;
	
}

sub novel_bed_comp () {
	my $comp = shift(@_);
	my @beds = @_;
	create_temp_file($comp, "ctemp", "comp");
	create_temp_cat_file(@beds, "stemp", 0);
	open(CNT, "intersectBed -a ctemp -b stemp -v |");
	my $count = 0;
	while(my $line = <CNT>){
		$count++;
	}
	close CNT;
	open(SUB, "subtractBed -a ctemp -b stemp |");
	my @novel_lens;
	while(my $line = <SUB>){
		chomp $line;
		my @segs = split(/\t/, $line);
		push(@novel_lens, $segs[2] - $segs[1]);
	}
	close SUB;
	my $len_val = sum(\@novel_lens);
	return($count, $len_val);
}

sub multi_way_comp () {
	# Take bed files in an ordered array, then do paired calculations of each
	# First bed is the comp, number is how many beds to remove from the end of the array
	my $remove = shift(@_);
	my $p_ref = shift(@_);
	my $prefix = shift(@_);
	my %file_prefix = %{$p_ref};
	my @orig = @_;
	my @beds = @orig;
	if ($beds[0] eq undef || $beds[0] eq ''){
		shift(@beds);
	}
	my @hold_prefixes = split(/;/, $prefix); 
	# Just to determine which prefixes I need to keep in the interval detection
	# remove intervals from both the comp and subjects
	my $comp = shift(@beds);
	my @comp_pre = grep($_ ne $file_prefix{$comp}, @hold_prefixes); 
	my $share_int = join(';', @comp_pre); #this should be the shared interval prefix for the bedtools merge section
	if($comp eq '' || !defined($comp)){ $comp = shift(@beds);}
	if ($remove > 1){
		my @rem;
		for (my $x = 0; $x < $remove; $x++){
			my $holder = pop(@beds);
			push (@rem, $holder);
		}
		create_temp_cat_file(@rem, "rtemp", \%file_prefix);
	} elsif ($remove == 0){
		open(OUT, "> rtemp");
		print OUT "chr20\t300000000\t300000001\n";
		close OUT;
	} else {
		my $holder = pop(@beds);
		create_temp_file($holder, "rtemp");
		push(@beds, $holder);
	}
	create_temp_file($comp, "ctemp", $file_prefix{$comp});
	create_temp_cat_file(@beds, "stemp", \%file_prefix);
	remove_bed_intervals("rtemp", "ctemp", "crtemp", 1);
	remove_bed_intervals("rtemp", "stemp", "srtemp", 1);
	
	# Calculate interval overlap
	my $count = 0;
	open(CNT, "intersectBed -a crtemp -b srtemp -wb |");
	#system ("intersectBed -a crtemp -b srtemp -wa > inttemp");
	#open (CNT, "mergeBed -i inttemp |");
	my $int_lengths;
	while (my $line = <CNT>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if ($segs[7] eq $share_int){
			$count++;
			$int_lengths += $segs[2] - $segs[1];
		}
	}
	close CNT;
	
	# Calculate base overlap
=pod
	my @bedmerge_cat;
	my @indiv_rem;
	foreach my $orf (@orig){
		if (!defined($orf)){next;}
		my $ofname = $file_prefix{$orf} . "_temp";
		my $ofnamer = "$ofname\_r";
		create_temp_file($orf, $ofname, $file_prefix{$orf});
		push(@indiv_rem, $ofname);
		remove_bed_intervals("rtemp", $ofname, $ofnamer, 1);
		push(@bedmerge_cat, $ofnamer);
	}
	create_temp_cat_file(@bedmerge_cat, "mtemp", \%file_prefix);
	my @lengths;
	open(INT, "mergeBed -i mtemp -nms |");
	while (my $line = <INT>){
		chomp $line;
		my @segs = split(/\t/, $line);
		my @reformat = split(/\;/, $segs[3]);
		my $ref_pref = join(';', sort {$a cmp $b} @reformat);
		if ($prefix eq $ref_pref){
			push(@lengths, $segs[2] - $segs[1]);
		}
	}
	close INT;
	my @remove = (@bedmerge_cat, @indiv_rem);
	my $r_string = join(' ', @remove);
	#system("rm $r_string");
	my $int_lengths = sum(\@lengths);
=cut
	return($count, $int_lengths);
}	

sub create_temp_file (){
	my ($a_ref, $output, $prefix) = @_;
	open(IN, "< $a_ref");
	open(OUT, "> $output");
	while(my $line = <IN>){
		chomp $line;
		$line =~ s/\r//g;
		my @segs = split(/\t/, $line);
		print OUT "$segs[0]\t$segs[1]\t$segs[2]\t$prefix\n";
	}
	close IN;
	close OUT;
}

sub create_temp_cat_file (){
	my $p_ref = pop(@_);
	my $output = pop(@_);
	my @refs = @_;
	open (OUT, "> tcat");
	my $prefix;
	foreach my $r (@refs){
		if ($p_ref == 0){
			$prefix = 0;
		} else {
			$prefix = $p_ref->{$r};
		}
		open(IN, "< $r");
		while(my $line = <IN>){
			chomp $line;
			$line =~ s/\r//g;
			my @segs = split(/\t/, $line);
			if ($prefix){
				print OUT "$segs[0]\t$segs[1]\t$segs[2]\t$prefix\n";
			} else {
				print OUT "$segs[0]\t$segs[1]\t$segs[2]\t$segs[3]\n";
			}
		}
		close IN;
	}
	close OUT;
	open(BED, "mergeBed -i tcat -nms |");
	open(OUT, "> $output");
	while(my $line = <BED>){
		chomp $line;
		$line =~ s/\r//g;
		my @segs = split(/\t/, $line);
		my @nam = split(/\;/, $segs[3]);
		my %unique;
		foreach my $n (@nam){
			$unique{$n} = 1;
		}
		my $id = join(';', sort{$a cmp $b}(keys %unique));
		if (scalar(@nam) > 1){
			print OUT "$segs[0]\t$segs[1]\t$segs[2]\t$id\n";
		}else{
			print OUT "$segs[0]\t$segs[1]\t$segs[2]\t$segs[3]\n";
		}
	}
	#system("mergeBed -i tcat -nms > $output");
}
sub remove_bed_intervals (){
	# $v indicates if intersectBed -v should be used, if zero, subtractBed will be used instead
	my ($r_file, $s_file, $outfile, $v) = @_;
	if ($v){
		system("intersectBed -a $s_file -b $r_file -v > $outfile");
	} else{
		system("subtractBed -a $s_file -b $r_file > $outfile");
	}
}

sub sum (){
	my $a_ref = shift(@_);
	my $tot = 0;
	foreach my $val (@{$a_ref}){
		$tot += $val;
	}
	return $tot;
}
sub randomize_array (){
	my $array = shift(@_);
	my $i;
	for ($i = @$array; --$i;){
		my $j = int(rand($i + 1));
		next if $i == $j;
		@$array[$i,$j] = @$array[$j,$i];
	}
}