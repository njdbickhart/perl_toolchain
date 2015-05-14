#!/usr/bin/perl
# This script is designed to process a tab delimited file to generate some basic statistics on column field count
# It only works on columns that have a limited range of possible entries (like enumerated types)

use strict;
use Getopt::Std;

my $usage = "perl $0 [options] -f <input file> -c <column number>
	This script is designed to process a tab delimited file to generate some basic statistics on column field count.
REQUIRED:
	-f	Input space/tab delimited file. Beginning spaces are removed from each line. Trailing spaces are condensed
	-c	The column number (Zero based!) to count features
	
OPTIONAL:
	-o	Output file [Default option: STDOUT]
	-e	Comment prefix. Exclude lines that begin with this character [Default option: exclude nothing]
	-m	Markdown flag. Formats output into table format [Default option: tab delimited]
	-x	Numerical counter statistics (different logic, calculates stdev, avg, median, etc) [Default option: type counter]
\n";

my %opts;
getopt('fcoemx', \%opts);

unless(defined($opts{'f'}) && defined($opts{'c'})){
	print $usage;
	exit;
}

my $worker = ColumnCounter->new('colnum' => $opts{'c'}, 'mkdwn' => $opts{'m'});

if(defined($opts{'e'})){
	$worker->ignore($opts{'e'});
}

if(defined($opts{'o'})){
	$worker->output($opts{'o'});
}

$worker->readFile($opts{'f'});
$worker->writeOut();

exit;

BEGIN{
package ColumnCounter;
use Mouse;
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
	open(IN, "< $file") || die "[ColumnCounter] Could not open input file: $file\n";
	my $col = $self->colnum;
	my %hash;
	my $com;
	if($self->has_ignore){
		$com = $self->ignore;
	}
	
	while(my $line = <IN>){
		chomp($line);
		$line =~ s/^\s+//;
		
		if($self->has_ignore){
			if($line =~ m/^$com/){
				next;
			}
		}
		
		my @segs = split(/\t/, $line);
		$hash{$segs[$col]} += 1;
	}
	close IN;
	$self->counter(\%hash);
		
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
		my $collen = 0;
		foreach my $k (@sortedkeys){
			if(length($k) > $collen){
				$collen = length($k);
			}
		}
		$collen++;
		my $sepstr = '-' x ($collen - 1);
		print {$fh} sprintf("\|%-*s\|%*s\|\n", $collen, "Entry", $collen, "Count");
		print {$fh} sprintf("\|\:%s\|%s\:\|\n", $sepstr, $sepstr);
		foreach my $k (@sortedkeys){
			my $kvalue = $k;
			if($k eq ''){
				$kvalue = "<Null>";
			}
			print {$fh} sprintf("\|%-*s\|%*d\|\n", $collen, $kvalue, $collen, $self->Get($k));
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

}