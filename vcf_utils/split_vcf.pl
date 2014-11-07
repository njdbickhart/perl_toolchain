#!/usr/bin/perl

use strict;
use FileHandle;
use Time::HiRes;
use threads;

chomp(@ARGV);

my @fhs;
open(IN, "< $ARGV[0]") || die "Could not open file $ARGV[0]!\n";
print 
"#################################
#                               #
#       SUPER COW SPLITTER      #
#                               #
#   An Animals in Problematic   #
#      Locations Production     #
#################################
";
while(my $line = <IN>){
if($line =~ /#/){
  if($line =~ /#CHROM/){
      my @segs = split(/\s+/, $line);
      for(my $x = 9; $x < scalar(@segs); $x++){
	push(@fhs, FileHandle->new("$segs[$x].vcf", "w"));
      }
      last;
    }else{
      next;
    }
  }
}
close IN;

my $t = threads->new(\&goofyticker);

open(IN, "< $ARGV[0]") || die "Could not open file $ARGV[0]!\n";
while(my $line = <IN>){
  chomp $line;
  $line =~ s/\r//g;
  
  if($line =~ /^#/){
    if($line =~ /^#CHROM/){
      my @segs = split(/\s+/, $line);
      for(my $x = 9; $x < scalar(@segs); $x++){
        my $y = $x - 9;
        my @bsegs = @segs[0..8];
        my $f = $fhs[$y];
        print $f join("\t", @bsegs) . "\t" . $segs[$x] . "\n";
      }
    }else{
      foreach my $fh (@fhs){
        print $fh "$line\n";
      }
    }
    next;
  }
  
  my @segs = split(/\s+/, $line);
  if($segs[6] !~ /\./){
    next;
  }else{
    for(my $x = 9; $x < scalar(@segs); $x++){
      if($segs[$x] =~ /\.\/\./){
        next;
      }else{
        my $y = $x - 9;
        my $f = $fhs[$y];
        my @bsegs = @segs[0..8];
        my $str = join("\t", @bsegs);
        print $f "$str\t$segs[$x]\n";
      }
    }
  }
}
close IN;

$t->join();
exit;

sub goofyticker{
  require Time::HiRes;
  local $| = 1;
  my $sstatement = "WARNING!!! Too many cows input into program! Working on a solution...";
  for(my $x = 0; $x < 8; $x++){
    print $sstatement . "\r";
    Time::HiRes::usleep(550000);
    print "\e[K";
    Time::HiRes::usleep(250000);
  }
  
  my $newstatement = "GRANTING SUPER COW POWER!!!!!!!";
  for (my $y = 0; $y < length($newstatement); $y++){
    print substr($newstatement, $y, 1);
    Time::HiRes::usleep(250000);
  }
  print "\n";
  
  for (my $x = 0; $x < 3; $x++){
    for (my $y = 0; $y < 100; $y++){
      superscrollbar($y, 100);
      Time::HiRes::usleep(1000);
      
    }
    print "\n";
  }
  print "Super cow power successful\n";
}
sub superscrollbar{
  my ( $got, $total, $width, $char ) = @_;
    $width ||= 25;
    $char  ||= '=';
    my $num_width = length $total;
    local $| = 1;
    printf "|%-${width}s| \r", 
        $char x (($width-1)*$got/$total). '>';
}