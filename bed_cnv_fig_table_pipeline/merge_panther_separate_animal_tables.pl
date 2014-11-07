#!/usr/bin/perl
# I created this script to change the panther script merger output to create separate workbooks per animal
# This should also filter the results by George's specifications (ie, 5 or more category hits and significant p values)

use strict;
use Spreadsheet::WriteExcel;

my $arg = $ARGV[0];
chomp $arg;
my @files = `$arg`;
my $num_files = scalar(@files);
print "Found $num_files files...\n";
my ($file_suffix) = $arg =~ m/ls \*(.+)\..+/;
if ($num_files < 1){
	exit;
}
print "$file_suffix\n";
my %pathway_files;


foreach my $f (@files){
	chomp $f;
	my @f_segs = split(/\_/, $f);
	$pathway_files{$f_segs[1]}->{$f_segs[0]} = $f;
}

my @animals;
my %file_arrays;
my %panther_ref;
foreach my $fkeys (sort {$a cmp $b} keys(%pathway_files)){
	print "Working on file $fkeys\n";
	@animals = sort {$a cmp $b} keys(%{$pathway_files{$fkeys}});
	foreach my $akeys (sort {$a cmp $b} keys(%{$pathway_files{$fkeys}})){
		my @temp = open_file_create_array($pathway_files{$fkeys}->{$akeys});
		foreach my $l_ref (@temp){
			push(@{$file_arrays{$akeys}->{$fkeys}->{$l_ref->[0]}}, $l_ref->[2], $l_ref->[3], $l_ref->[4], $l_ref->[5]);
			$panther_ref{$l_ref->[0]} = $l_ref->[1];
		}
	}
}
print join(' ', @animals) . "\n";

foreach my $akeys (sort {$a cmp $b} keys(%file_arrays)){
	my $workbook = Spreadsheet::WriteExcel->new("$akeys\_$file_suffix.xls");
	print "Making $akeys excel\n";
	foreach my $fkeys (sort {$a cmp $b} keys(%{$file_arrays{$akeys}})){
	
		my $worksheet = $workbook->add_worksheet($fkeys);
		$worksheet->write(1,0, 'pathway/gene class');
		$worksheet->write(1,1, 'reference number');
		$worksheet->write(0,2, $akeys);
		$worksheet->write(1,2, 'num');
		$worksheet->write(1,3, 'expect');
		$worksheet->write(1,4, '+/-');
		$worksheet->write(1,5, 'pvalue');
		my $row = 2;
		foreach my $gkeys (sort{$a cmp $b} keys(%{$file_arrays{$akeys}->{$fkeys}})){
			$worksheet->write($row,0, "$gkeys");
			$worksheet->write($row,1, "$panther_ref{$gkeys}");
			$worksheet->write($row,2, "$file_arrays{$akeys}->{$fkeys}->{$gkeys}->[0]");
			
			$worksheet->write($row,3, "$file_arrays{$akeys}->{$fkeys}->{$gkeys}->[1]");
			
			$worksheet->write($row,4, "$file_arrays{$akeys}->{$fkeys}->{$gkeys}->[2]");
			
			$worksheet->write($row,5, "$file_arrays{$akeys}->{$fkeys}->{$gkeys}->[3]");
			$row++;
		}
	}
	$workbook->close();
}
exit;

sub open_file_create_array (){
	my $file = shift(@_);
	open(IN, "< $file");
	my @return;
	# -> 0 = gene name   1 = panther ref  2= in dataset 3= expected 4= up or down  5 = p value
	while(my $line = <IN>){
		$line =~ s/\r//g;
		if ($line =~ m/^[\s\n]/g){
			next;
		} elsif ($line =~ m/REFLIST/){
			next;
		}
		my @segs = split(/\t/, $line);
		if ($segs[2] < 5 || $segs[5] > 0.05){
			next;	# Filter
		}else{
		push (@return, [@segs]);
		}
	}
	close IN;
	return @return;
}
