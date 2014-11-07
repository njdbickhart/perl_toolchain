#!/usr/bin/perl
# This script takes a list of genes and a list of copy number files and checks for intersections
# It then prints them as separate sheets on an excel workbook based on the number of gene databases fed in
# -i list of database file intersections (chr, start, end, accession)
# -c CN coordinate lists per animal (files must start with the animal name)
# -v merged cnvs with animal_type designation merged by bedtools mergebed
# Optional: -p creates a table that can be uploaded for integrated pathway analysis
# 7/28/2011: added the ability to create separate gene xls files for ingenuity ipa upload
# 9/27/2012: added binning and changed input file format

use strict;
use Getopt::Std;
use Spreadsheet::WriteExcel;

my $bedtools = "";

my $usage = "$0 -i \<database files\> -c \<copy number\> -v \<cnv list\> -p \<database name for ipa\>
-s \<skip major gene db table\>
\ninputs are text files with the locations of the files\n";

my %opts;
getopt('icvps', \%opts);

unless(defined($opts{i}) && defined($opts{c})){
	print $usage;
	exit;
}
if ($opts{s} && !defined($opts{p})){
	print $usage;
	exit;
}

my %database_coords; #{db}->{chr}->{bin}->[start, end, accession]
my @file_list;
open (IN, "< $opts{i}");
while(my $line = <IN>){
	chomp $line;
	$line =~ s/\r//g;
	my @p_seg = split(/\t/, $line);
	#my @f_seg = split(/\_/, $p_seg[-1]);
	my $db = $p_seg[1];
	if($opts{s}){           
		if($db != $opts{p}){
			print "Skipping $db\n";
			next;
		}
	}
	push(@file_list, $line);
	open (DB, "< $p_seg[0]");
	#$db =~ tr/[a-z]/[A-Z]/;
	while (my $cnv = <DB>){
		chomp $cnv;
		$cnv =~ s/\r//g;
		my @seg = split(/\t/, $cnv);
		#my $bin = getbin($seg[1], $seg[2]);
		push(@{$database_coords{$db}}, [$seg[0], $seg[1], $seg[2], $seg[3]]);
	}
	close DB;
}
print "Number of database files: " . scalar(@file_list) . "\n";
close IN;

my %an_cn_numbers; #{animal}->{chr}->{bin}->[start, end, value]
open (IN, "< $opts{c}");
while (my $line = <IN>){
	chomp $line;
	
	my @p_seg = split(/\t/, $line);
	open (CN, "< $p_seg[0]");
	#my @f_seg = split(/\_/, $p_seg[-1]);
	my $an = $p_seg[1];
	#$an =~ tr/[a-z]/[A-Z]/;
	while (my $cn = <CN>){
		chomp $cn;
		my @seg = split(/\t/, $cn);
		my $bin = getbin($seg[1], $seg[2]);
		push(@{$an_cn_numbers{$an}->{$seg[0]}->{$bin}}, [$seg[1], $seg[2], $seg[3]]);
	}
	close CN;
	print "Loaded $an\n";
}
close IN;

my %an_cnv_coords; #{chr}->{bin}->[row]->[start, end, animalstr, type]
my @cnv_f_list;
open (IN, "< $opts{v}");
while(my $line = <IN>){
	chomp $line;
	#push(@cnv_f_list, $line);
	my @segs = split(/\t/, $line);
	my $bin = getbin($segs[1], $segs[2]);
	my @nams = split(/\;/, $segs[3]);
	my @store;
	my $t = "";
	foreach my $n (@nams){
		my @type = split(/\_/, $n);
		push(@store, $type[0]);
		if($t eq ""){
			$t = $type[1];
		}elsif($t ne $type[1]){
			$t = "both";
		}else{
			$t = $type[1];
		}
	}
	push(@{$an_cnv_coords{$segs[0]}->{$bin}}, [$segs[1], $segs[2], join(";", @store), $t]);
}
print "Loaded CNVs\n";
=pod
# No longer needed since I will do the bedtools merger myself
print "Number of cnv files: " . scalar(@cnv_f_list) . "\n";
close IN;
print "Bedtools merging...\n";
my $cat = join(' ', @cnv_f_list);
system("cat $cat > temp");
open (BED, "mergeBed -i temp -nms |");
my %merge_cnvrs;
my $wc = 0;
while (my $line = <BED>){
	chomp $line;
	my @seg = split(/\t/, $line);
	my @tags = split(/\;/, $seg[3]);
	my $type;
	my $m_string;
	foreach my $t (@tags){
		my @t_segs = split(/\_/, $t);
		if (!defined($type)){
			$type = $t_segs[1];
		} elsif(length($type) > 1 && $type ne $t_segs[1]){
			$type = "both";
		}
		$t_segs[0] =~ tr/[a-z]/[A-Z]/;
		$m_string .= "$t_segs[0];";
	}
	chop $m_string;
	push (@{$merge_cnvrs{$seg[0]}}, [$seg[1], $seg[2], $m_string, $type]);
	$wc++;
}
close BED;
=cut
#print "Merger output lines: $wc\n";

my @animal_list = sort{$a cmp $b} (keys %an_cn_numbers);
my $num = scalar(@animal_list);
unless($opts{s}){
my $workbook = Spreadsheet::WriteExcel->new('gene_list_cn_100bulls.xls');
foreach my $gene_db (sort{$a cmp $b} keys(%database_coords)){
	my $worksheet = $workbook->add_worksheet($gene_db);
	$worksheet->write(0,0,"gene ID");
	$worksheet->write(0,1,"chr");
	$worksheet->write(0,2,"start");
	$worksheet->write(0,3,"end");
	$worksheet->write(0,4,"gene size");
	$worksheet->write(0,5,"covered bases");
	$worksheet->write(0,6,"covered perc");
	#$worksheet->write(0,7,"cnv type");
	my $col = 8;
	foreach my $a (@animal_list){
		$worksheet->write(0,$col,"$a");
		$col++;
	}
	$worksheet->write(0,$col,"animals with cnvs");
	my $row = 1;
	print "Working on $gene_db\n";
	foreach my $db_row (@{$database_coords{$gene_db}}){
		$worksheet->write($row,0,"$db_row->[3]");
		$worksheet->write($row,1,"$db_row->[0]");
		$worksheet->write($row,2,"$db_row->[1]");
		$worksheet->write($row,3,"$db_row->[2]");
		my $length = $db_row->[2] - $db_row->[1];
		$worksheet->write($row,4,"$length");
		print "$db_row->[3]\t$db_row->[0]\t$db_row->[1]\t$db_row->[2]";
		my ($perc, $an_string, $g_l) = percent_covered($db_row->[0], $db_row->[1], $db_row->[2], \%an_cnv_coords);
		my $cover = $length * $perc / 100;
		$worksheet->write($row,5,"$cover");
		$worksheet->write($row,6,"$perc");
		$worksheet->write($row,7,"$g_l");
		print "\t$perc";
		my $d_col = 8;
		foreach my $animal (@animal_list){
			my ($holder, $val) = cross_index_cn($db_row->[0], $db_row->[1], $db_row->[2], $an_cn_numbers{$animal}, $animal);
			$worksheet->write($row,$d_col,"$val");
			$d_col++;
			print "\t$val";
		}
		$worksheet->write($row,$d_col,"$an_string");
		print "\t$an_string\n";
		$row++;
	}
}
}
if(defined($opts{p})){
	print "Working on pathway analysis\n";
	my $used_db = $opts{p};
	unless (exists($database_coords{$used_db})){
		print "Could not find database name!\n";
		exit;
	}
	my $ingenuity = Spreadsheet::WriteExcel->new('ingenuity_upload.xls');
	my $worksheet = $ingenuity->add_worksheet($used_db);
	$worksheet->write(0,0,"ID");
	$worksheet->write(0,1,"chr");
	$worksheet->write(0,2,"start");
	$worksheet->write(0,3,"end");
	$worksheet->write(0,4,"gene size");
	$worksheet->write(0,5,"covered bases");
	$worksheet->write(0,6,"covered perc");
	$worksheet->write(0,7,"cnv type");
	my $col = 8;
	foreach my $a (@animal_list){
		$worksheet->write(0,$col,"$a");
		$col++;
	}
	print "Working on the ingenuity_upload\n";
	my $row = 1;
	foreach my $db_row (@{$database_coords{$used_db}}){
		my $length = $db_row->[2] - $db_row->[1];
		my ($perc, $an_string, $g_l) = percent_covered($db_row->[0], $db_row->[1], $db_row->[2], \%an_cnv_coords);
		my $cover = $length * $perc / 100;
		if($cover == 0){ next;}
		$worksheet->write($row,0,"$db_row->[3]");
		$worksheet->write($row,1,"$db_row->[0]");
		$worksheet->write($row,2,"$db_row->[1]");
		$worksheet->write($row,3,"$db_row->[2]");		
		$worksheet->write($row,4,"$length");
		#print "$db_row->[3]\t$db_row->[0]\t$db_row->[1]\t$db_row->[2]";
		$worksheet->write($row,5,"$cover");
		$worksheet->write($row,6,"$perc");
		$worksheet->write($row,7,"$g_l");
		#print "\t$perc";
		my $d_col = 8;
		foreach my $animal (@animal_list){
			my ($holder, $val) = cross_index_cn($db_row->[0], $db_row->[1], $db_row->[2], $an_cn_numbers{$animal}, $animal);
			my $norm = 0;
			if ($val < 1.5){
				$norm = -1;
			}elsif($val > 2.5){
				$norm = 1;
			}
			$worksheet->write($row,$d_col,"$norm");
			$d_col++;
			#print "\t$val";
		}
		#$worksheet->write($row,$d_col,"$an_string");
		#print "\t$an_string\n";
		$row++;
	}
}

exit;

sub cross_index_cn (){
	# Takes input coordinates and a reference to a copy number array and calculates average copy number for those coords
	my ($chr, $start, $end, $cn_ref, $a_name) = @_;
	#print "$a_name ";
	my @values;
	my $term = 0;
	my @bins = searchbins($start, $end);
	#foreach my $bref (keys%{$cn_ref->{$chr}}){
		foreach my $b (@bins){
			foreach my $l_ref (@{$cn_ref->{$chr}->{$b}}){
		
				if ($l_ref->[1] > $start && $l_ref->[0] < $end){
					push (@values, $l_ref->[2]);
					$term = 1;
				}
			}
		}
	#}
	
	my $avg1;
	my $avg2;
	my $num = scalar(@values);
	if ($num == 0){
		$num = 1;
	}
	foreach my $v (@values){
		$avg1 += $v;
	}
	$avg2 = $avg1 / $num;
	return($a_name, $avg2);

}

sub percent_covered (){
	# Takes input coordinates and an hash ref (chr index). Returns the percentage of the coordinates covered by the array
	# also returns a string with the animals that have the cnvr
	my ($chr, $start, $end, $h_ref) = @_;
	my $len = $end - $start;
	my $fstart;
	my $fend;
	my $occlusion;
	my %u_animals;
	my $type;
	my @bins = searchbins($start, $end);
	#foreach my $bref (keys%{$h_ref->{$chr}}){
		foreach my $b (@bins){
			foreach my $row (@{$h_ref->{$chr}->{$b}}){
				if($end > $row->[0] && $start < $row->[1]){
					if(!defined($type)){$type = $row->[3];}
					elsif($type ne $row->[3]){$type = "both";}
					if($row->[0] < $start){
						$fstart = $start;
					} else {
						$fstart = $row->[0];
					}
					if ($row->[1] > $end){
						$fend = $end;
					} else {
						$fend = $row->[1];
					}
					$occlusion += $fend - $fstart;
					my @temp = split(/\;/, $row->[2]);
					foreach (@temp){
						$u_animals{$_} = 1;
					}
				}
			}
		}
	#}
	my @t_an = sort{$a cmp $b} keys(%u_animals);
	my $an_string = join(';', @t_an);
	
	my $perc;
	if($occlusion >= $len){
		$perc = "100";
	} elsif ($occlusion < $len) {
		if ($len == 0){ $len = 1;}
		$perc = ($occlusion / $len) * 100;
	} else {
		$perc = 0;
	}
	return ($perc, $an_string, $type);
}

sub most {
	my($a, $b) = @_;
	return ($a > $b)? $a : $b;
}
sub least {
	my ($a, $b) = @_;
	return ($a < $b) ? $a : $b;
}
sub overlap {
	my ($s1, $e1, $s2, $e2) = @_;
	return (least($e1, $e2) - most($s1, $s2));
}
sub searchbins{
	my ($start, $end) = @_;
	my $genomiclength = 536870912;
	my %tmph;
	my $href = collectbins($start, $end, 0,0,0,0,1,0,$genomiclength, \%tmph);
	my @vals = keys(%{$href});
	return @vals;
}
sub getbin{
	my ($start, $end) = @_;
	my $genomiclength = 536870912;
	return calcbin($start, $end, 0,0,0,0,1,0,$genomiclength);
}

# private method
sub collectbins{
	my ($start, $end, $binid, $level, $binrowstart, $rowindex, $binrowcount, $genomicpos, $genomiclength, $hset) = @_;
	my $maxlevel = 4;
	my $childrencount = 8;
	$hset->{$binid} = 1;
	if($level >= $maxlevel){
		return $hset;
	}

	my $childlength = $genomiclength / $childrencount;
	my $childbinrowcount = $binrowcount * $childrencount;
	my $childrowbinstart = $binrowstart + $binrowcount;
	my $firstchildindex = $rowindex * $childrencount;
	my $firstchildbin = $childrowbinstart + $firstchildindex;
	for (my $i = 0; $i < $childrencount; ++$i){
		my $childstart = $genomicpos + $i * $childlength;
		if($start > ($childstart + $childlength) || $end < $childstart){
			next;
		}
		collectbins($start, $end, $firstchildbin + $i, $level + 1, $childrowbinstart, $firstchildindex + $i, $childbinrowcount, $childstart, $childlength, $hset);
	}
	return $hset;
}

# private method
sub calcbin{
	my ($start, $end, $binid, $level, $binrowstart, $rowindex, $binrowcount, $genomicpos, $genomiclength) = @_;
	my $maxlevel = 4;
	my $childrencount = 8;
	if ($start >= $genomicpos && $end <= ($genomicpos + $genomiclength)){
		if($level >= $maxlevel){
			return $binid;
		}
		my $childlength = $genomiclength / $childrencount;
		my $childbinrowcount = $binrowcount * $childrencount;
		my $childrowbinstart = $binrowstart + $binrowcount;
		my $firstchildindex = $rowindex * $childrencount;
		my $firstchildbin = $childrowbinstart + $firstchildindex;
		for (my $i = 0; $i < $childrencount; ++$i){
			my $n = calcbin($start, $end, $firstchildbin+$i, $level + 1, $childrowbinstart, $firstchildindex + $i, $childbinrowcount, $genomicpos + $i * $childlength, $childlength);
			if($n != -1){
				return $n;
			}
		}
		return $binid;
	}
	return -1;
}	