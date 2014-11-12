#!/usr/bin/perl -w

use strict;
use POSIX;


my $sum=0;
my $sqr=0;

my @data = ();
my %dataDis = ();
while(<STDIN>)
{
  chomp;
  next if( $_ !~ /^-?\d+\.?\d*$/ );

  $sum=$sum + $_;
  $sqr=$sqr + $_ * $_;

  push @data, $_;

  $dataDis{$_}++   if( defined $dataDis{$_} );
  $dataDis{$_} = 1 if( !defined $dataDis{$_} );
}

@data = sort { $a <=> $b } @data;

my $ct = scalar(@data);
print STDOUT "total\t$ct\n";
print STDOUT "Minimum\t", $data[0],     "\n";
print STDOUT "Maximum\t", $data[$ct-1], "\n";
#print STDOUT "Average\t", ($sum/$ct),   "\n";
print STDOUT "Average\t";
printf("%f", ($sum/$ct));
print STDOUT   "\n";

my $med = ($ct % 2 == 0) ? ($data[$ct/2-1] + $data[$ct/2])/2 : $data[floor($ct/2)];
print STDOUT "Median\t", $med, "\n";

print STDOUT "Standard Deviation\t";
printf("%f", sqrt( ($sqr * $ct - $sum * $sum)/($ct * ($ct -1)) ) ) unless (1 == $ct);
print STDOUT "\n";


# q50 mode is the value where highest occurance across the whole dataset is found
my ($q50, $maxOcc) = (undef, -1);
foreach my $key (keys %dataDis)
  {
    if( $dataDis{$key} > $maxOcc)
      {
	$maxOcc = $dataDis{$key};
	$q50    = $key;
      }
  }
print STDOUT "Mode(Highest Distributed Value)\t$q50\n\n";
