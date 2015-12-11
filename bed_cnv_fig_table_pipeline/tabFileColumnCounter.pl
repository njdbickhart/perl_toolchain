#!/usr/bin/perl
# This script is designed to process a tab delimited file to generate some basic statistics on column field count
# It only works on columns that have a limited range of possible entries (like enumerated types)
# 6/12/2015: added a STDIN functionality

use strict;
use Getopt::Std;

my $usage = "perl $0 [options] -f <input file> -c <column number>
	This script is designed to process a tab delimited file to generate some basic statistics on column field count.
REQUIRED:
	-f	Input space/tab delimited files separated by commas OR \"stdin\" for STDIN. Beginning spaces are removed from each line. Trailing spaces are condensed
	-c	The column number (Zero based!) to count features
	
OPTIONAL:
	-o	Output file [Default option: STDOUT]
	-e	Comment prefix. Exclude lines that begin with this character [Default option: exclude nothing]
	-m	Markdown flag. Formats output into table format [Default option: tab delimited]
	-x	Numerical counter statistics (different logic, calculates stdev, avg, median, etc) [Default option: type counter]
\n";

my %opts;
getopt('fcoex', \%opts);

unless(defined($opts{'f'}) && defined($opts{'c'})){
	print $usage;
	exit;
}

my $mbool = (defined($opts{'m'}))? 1 : 0;

# Determine how many files we're working with here.
my @files = split(/,/, $opts{'f'});
if(scalar(@files) == 1){
	# Just a single comparison
	my $worker = ColumnCounter->new('colnum' => $opts{'c'}, 'mkdwn' => $mbool);
	
	if(defined($opts{'e'})){
		$worker->ignore($opts{'e'});
	}
	
	if(defined($opts{'o'})){
		$worker->output($opts{'o'});
	}
	$worker->readFile($opts{'f'});
	$worker->writeOut();
}else{
	# Now we have more to work with
	my $manager = SummaryManager->new('mkdwn' => $mbool);
	
	if(defined($opts{'o'})){
		$manager->output($opts{'o'});
	}
	for (my $x = 0; $x < scalar(@files); $x++){
		my $filename = "File" . ($x + 1);
		
		my $worker = ColumnCounter->new('colnum' => $opts{'c'}, 'mkdwn' => $mbool);
		
		if(defined($opts{'e'})){
			$worker->ignore($opts{'e'});
		}
		
		$worker->readFile($files[$x]);
		$manager->setTable($filename, $worker->createSumTable($filename));
	}
	
	if(defined($opts{'o'})){
		open(OUT, "> $opts{o}");
		for(my $x = 0; $x < scalar(@files); $x++){
			my $filename = "File" . ($x + 1);
			print OUT "$filename\:  $files[$x]\n";
		}
		print OUT "\n";
		close OUT;
	}else{
		for(my $x = 0; $x < scalar(@files); $x++){
			my $filename = "File" . ($x + 1);
			print "$filename\:  $files[$x]\n";
		}
		print "\n";
	}
	
	$manager->PrintResults();
}



exit;

BEGIN{
package ColumnCounter;
use Mouse;
use FileHandle;
use namespace::autoclean;

# Required attributes
has 'colnum' => (is => 'ro', isa => 'Int', required => 1);
has 'mkdwn' => (is => 'ro', isa => 'Bool', required => 1);

has 'output' => (is => 'rw', isa => 'Str', predicate => 'has_output');
has 'ignore' => (is => 'rw', isa => 'Str', predicate => 'has_ignore');
has 'counter' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'Keys' => 'keys',
		'Empty' => 'is_empty',
		'Get' => 'get',
	});
	
sub readFile{
	my ($self, $file) = @_;
	my $fh; # Filehandle
	if($file eq "stdin"){
		$fh = *STDIN;
	}else{
		$fh = FileHandle->new();
		$fh->open("< $file") || die "[ColumnCounter] Could not open input file: $file\n";
	}
	my $col = $self->colnum;
	my %hash;
	my $com;
	if($self->has_ignore){
		$com = $self->ignore;
	}
	
	while(my $line = <$fh>){
		chomp($line);
		$line =~ s/\r//g;
		$line =~ s/^\s+//;
		
		if($self->has_ignore){
			if($line =~ m/^$com/){
				next;
			}
		}
		
		my @segs = split(/\t/, $line);
		$hash{$segs[$col]} += 1;
	}
	close $fh;
	$self->counter(\%hash);
		
}

sub createSumTable{
	my ($self, $name) = @_;
	my $table = CountSummaryTable->new('name' => $name);
	my @sortedkeys = sort{$a cmp $b} $self->Keys;
	foreach my $k (@sortedkeys){
		my $kvalue = $k; 
		if($k eq ''){
			$kvalue = "<Null>";
		}
		my $count = $self->Get($k);

		$table->addEntry($kvalue);
		$table->addValue($count);
	}
	return $table;
}

sub writeOut{
	my ($self) = @_;
	my $fh; 
	if($self->has_output){
		my $out = $self->output;
		open(OUT, "> $out");
		$fh = *OUT;
	}else{
		$fh = *STDOUT;
	}
	
	my @sortedkeys = sort{$a cmp $b} $self->Keys;	
	
	if($self->mkdwn){
		# For tidiness with markdown, we want proper column spacing
		my $collen = 5;
		my $conlen = 5;
		foreach my $k (@sortedkeys){
			if(length($k) > $collen){
				$collen = length($k);
			}
			if(length($self->Get($k)) > $conlen){
				$conlen = length($self->Get($k));
			}
		}
		$collen++;
		$conlen++;
		my $sepstr = '-' x ($collen - 1);
		my $sepcon = '-' x ($conlen - 1);
		print {$fh} sprintf("\|%-*s\|%*s\|\n", $collen, "Entry", $conlen, "Count");
		print {$fh} sprintf("\|\:%s\|%s\:\|\n", $sepstr, $sepcon);
		foreach my $k (@sortedkeys){
			my $kvalue = $k;
			if($k eq ''){
				$kvalue = "<Null>";
			}
			print {$fh} sprintf("\|%-*s\|%*d\|\n", $collen, $kvalue, $conlen, $self->Get($k));
		}
	}else{
		print {$fh} "Entry\tCount\n";
		foreach my $k (@sortedkeys){
			my $kvalue = $k; 
			if($k eq ''){
				$kvalue = "<Null>";
			}
			my $count = $self->Get($k);
			
			print {$fh} "$kvalue\t$count\n";
		}
	}
	if($self->has_output){
		close OUT;
	}
}

__PACKAGE__->meta->make_immutable;

package SummaryManager;
use Mouse;
use namespace::autoclean;

# Name => AbsSummaryTable
has 'tables' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'setTable' => 'set', 
		'getTable' => 'get',
		'tableKeys' => 'keys',
	});
	
has 'output' => (is => 'rw', isa => 'Str', predicate => 'has_output');
has 'mkdwn' => (is => 'ro', isa => 'Bool', required => 1);

sub PrintResults{
	my ($self) = @_;
	
	my $fh; 
	if($self->has_output){
		my $out = $self->output;
		open(OUT, ">> $out");
		$fh = *OUT;
	}else{
		$fh = *STDOUT;
	}
	my @sortedkeys = sort{$a cmp $b} $self->tableKeys;
	
	foreach my $k (@sortedkeys){
		print {$fh} $self->getTable($k)->formatHeader($self->mkdwn);
	}
	if($self->mkdwn){
		print {$fh} "\|\n";
	}else{
		print {$fh} "\n";
	}
	
	if($self->mkdwn){
		foreach my $k (@sortedkeys){
			print {$fh} $self->getTable($k)->formatSep($self->mkdwn);
		}
		print {$fh} "\|\n";
	}
	
	my @output;
	my $maxlen = 0;
	my @blank;
	foreach my $k (@sortedkeys){
		if($self->getTable($k)->numEntry > $maxlen){
			$maxlen = $self->getTable($k)->numEntry;
		}
		push(@output, $self->getTable($k)->formatOutput($self->mkdwn));
		push(@blank, $self->getTable($k)->formatBlank($self->mkdwn));
	}
	
	for( my $j = 0; $j < $maxlen; $j++){
		for( my $i = 0; $i < scalar(@output); $i++){
			if($j >= scalar(@{$output[$i]})){
				print {$fh} $blank[$i];
			}else{
				print {$fh} $output[$i]->[$j];
			}
		}
		if($self->mkdwn){
			print {$fh} "\|\n";
		}else{
			print {$fh} "\n";
		}
	}
	print "\n";
	if($self->has_output){
		close OUT;
	}
}

__PACKAGE__->meta->make_immutable;

package CountSummaryTable;
use Mouse;
use namespace::autoclean;

with 'AbsSummaryTable';

has 'name' => (is => 'ro', isa => 'Str', required => 1);

has 'entryTab' => (traits => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub{[]},
	handles => {
		'addEntry' => 'push',
		'getEntry' => 'get',
		'numEntry' => 'count',
		'allEntry' => 'elements',
	});
has 'valueTab' => (traits => ['Array'], is => 'rw', isa => 'ArrayRef[Int]', default => sub{[]},
	handles => {
		'addValue' => 'push',
		'getValue' => 'get',
		'numValue' => 'count',
		'allValue' => 'elements',
	});

has 'entryLen' => (is => 'rw', isa => 'Int', builder => '_entryLenCount', lazy => 1);
has 'valueLen' => (is => 'rw', isa => 'Int', builder => '_valueLenCount', lazy => 1);

sub formatOutput{
	my ($self, $mrkdown) = @_;
	my @values;
	if($mrkdown){
		for(my $x = 0; $x < $self->numEntry; $x++){
			push(@values, sprintf("\|%-*s\|%*d", $self->entryLen, $self->getEntry($x), $self->valueLen, $self->getValue($x)));
		}
	}else{
		for(my $x = 0; $x < $self->numEntry; $x++){
			push(@values, $self->getEntry($x) . "\t" . $self->getValue($x) . "\t");
		}
	}
	return \@values;
}

sub formatBlank{
	my ($self, $mrkdown) = @_;
	if($mrkdown){
		return sprintf("\|%-*s\|%*s", $self->entryLen, " ", $self->valueLen, " ");
	}else{
		return "\t\t";
	}
}

sub formatHeader{
	my ($self, $mrkdown) = @_;
	my $str;
	if($mrkdown){
		$str = sprintf("\|%-*s\|%*s", $self->entryLen, $self->name, $self->valueLen, "Count");
	}else{
		$str = $self->name . "\t";
		$str .= "Count\t";
	}
	return $str;
}

sub formatSep{
	my ($self, $mrkdown) = @_;
	my $str;
	if($mrkdown){
		my $sepstr = '-' x ($self->entryLen - 1);
		my $sepcon = '-' x ($self->valueLen - 1);
		$str = sprintf("\|\:%s\|%s\:", $sepstr, $sepcon);
	}else{
		$str = "";
	}
	return $str;
}

sub _entryLenCount{
	my ($self) = @_;
	my $collen = length($self->name);
	foreach my $k ($self->allEntry){
		if(length($k) > $collen){
			$collen = length($k);
		}
	}
	$collen++;
	return $collen;
}

sub _valueLenCount{
	my ($self) = @_;
	my $collen = 5;
	foreach my $k ($self->allValue){
		if(length($k) > $collen){
			$collen = length($k);
		}
	}
	$collen++;
	return $collen;
}

__PACKAGE__->meta->make_immutable;

package AbsSummaryTable;
use Mouse::Role;
#use namespace::autoclean;


requires 'formatOutput';
requires 'formatHeader';
requires 'formatSep';
requires 'formatBlank';

}
