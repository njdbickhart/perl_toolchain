#!/usr/bin/perl
# Fastqc data parser
# Updated to create a sum total file and to also check all of the files in a directory

use strict;
use Class::Struct;
use Cwd;
use File::Find;
use File::Basename;

struct('stats' => {
	'pass' => '$',
	'file' => '$',
	'totseq' => '$',
	'filtseq' => '$',
	'seqlen' => '$',
	'gc' => '$',
});

struct('pbpquality' => {
	'pass' => '$',
	'mean' => '%',
	'median' => '%',
});

struct('pseqquality' => {
	'pass' => '$',
	'count' => '%',
});

struct('pbpgccontent' => {
	'pass' => '$',
	'gc' => '%',
});

struct('pbpncontent' => {
	'pass' => '$',
	'ns' => '%',
});

struct('indivrecord' => {
	'stats' => '$',
	'bpquality' => '$',
	'seqquality' => '$',
	'gccontent' => '$',
	'ncontent' => '$',
	'seqlens' => '$',
	'seqdups' => '$',
	'overrep' => '$',
	'kmercon' => '$',
});

my $usage = "$0 <base directory> <output file base name>\n";

chomp(@ARGV);
if(scalar(@ARGV) != 2){
	print $usage;
	exit;
}

opendir(DIR, "$ARGV[0]");
open(SUM, "> $ARGV[0]/$ARGV[1]\_fastqc_summary.tab") || die "Could not create output file!\n";
print SUM "Animal\tTotalReads\tFilteredReads\tReadLength\tTotalGCPerc\tBPQualityPass\tFirst25Qual\tSecond25Qual\tThird25Qual\tFourth25Qual\t";
print SUM "ReadQSPass\t#ReadsQS2-11\t#ReadsQS12-21\t#ReadsQS22-31\t#ReadsQS32-40\tBPGCContentPass\tFirst25GC\tSecond25GC\t";
print SUM "Third25GC\tFourth25GC\tBPNContentPass\tFirst25N\tSecond25N\tThird25N\tFourth25N\tSeqDupPerc\tOverRepPass\tKmerPass\n";
while(my $basedir = readdir(DIR)){
	if($basedir =~ /^\./ || !(-d "$ARGV[0]/$basedir")){next;}
	opendir(FDIR, "$ARGV[0]/$basedir");
	my %anhash;
	while(my $fastqc = readdir(FDIR)){
		if($fastqc =~ /^\./ || !(-d "$ARGV[0]/$basedir/$fastqc")){next;}
		my @files;
		
		find(sub{if($File::Find::name =~ m/fastqc_data\.txt/){push(@files, $File::Find::name);}}, "$ARGV[0]/$basedir/$fastqc");

		#print OUT "Fastqfile\tTotalReads\tFilteredReads\tReadLength\tTotalGCPerc\tBPQualityPass\tFirst25Qual\tSecond25Qual\tThird25Qual\tFourth25Qual\t";
		#print OUT "ReadQSPass\t#ReadsQS2-11\t#ReadsQS12-21\t#ReadsQS22-31\t#ReadsQS32-40\tBPGCContentPass\tFirst25GC\tSecond25GC\t";
		#print OUT "Third25GC\tFourth25GC\tBPNContentPass\tFirst25N\tSecond25N\tThird25N\tFourth25N\tSeqDupPerc\tOverRepPass\tKmerPass\n";

		
		my $filenum = 0;
		my $totfiles = scalar(@files);
		if($totfiles == 0){
			print "Found no fastqc_data.txt files in dir!\n$usage";
			exit;
		}
		my @indiv = split(/[\._]/, basename($fastqc));
		foreach my $f (@files){
			chomp $f;
			open(IN, "< $f") || die "Could not open $f!\n";
			my $fq = indivrecord->new();
			my $head = <IN>;
			
			
			
			if(!exists($anhash{$indiv[0]})){
				$anhash{$indiv[0]} = indivrecord->new();
			}
			
			my $store;
			my $inloop = 0;
			while(my $line = <IN>){
				if($line =~ /^>>/ && $line !~ />>END_MODULE/){
					$inloop = 1;
					$store .= $line;
				}elsif($line =~ />>END_MODULE/){
					$inloop = 0;
					my($name, $ret) = process_module_seg($store);
					$store = '';
					
					my $working = $anhash{$indiv[0]};
					
					if($name eq "Per base sequence quality"){
						$fq->bpquality($ret);
						foreach my $k (keys(%{$ret->mean()})){
							if(defined($working->bpquality)){
								$working->bpquality->mean->{$k} = ($working->bpquality->mean()->{$k} + $ret->mean()->{$k}) / 2;
							}else{
								$working->bpquality(pbpquality->new('mean' => {}));
								$working->bpquality->mean->{$k} = $ret->mean->{$k};
							}
							
						}
						if($ret->pass() =~ m/fail/i || !defined($working->bpquality->pass())){
							$working->bpquality->pass($ret->pass());
						}
					}elsif($name eq "Per sequence quality scores"){
						$fq->seqquality($ret);
						foreach my $k (keys(%{$ret->count()})){
							if(defined($working->seqquality)){
								$working->seqquality->count->{$k} = $working->seqquality->count()->{$k} + $ret->count()->{$k};
							}else{
								$working->seqquality(pseqquality->new('count' => {}));
								$working->seqquality->count->{$k} = $ret->count->{$k};
							}
							
						}
						if($ret->pass() =~ m/fail/i || !defined($working->seqquality->pass())){
							$working->seqquality->pass($ret->pass());
						}
					}elsif($name eq "Per base sequence content"){
						$fq->gccontent($ret);
						foreach my $k (keys(%{$ret->gc()})){
							if(defined($working->gccontent)){
								$working->gccontent->gc->{$k} = ($working->gccontent->gc()->{$k} + $ret->gc()->{$k}) / 2;
							}else{
								$working->gccontent(pbpgccontent->new('gc' => {}));
								$working->gccontent->gc->{$k} = $ret->gc()->{$k};
							}							
						}
						if($ret->pass() =~ m/fail/i || !defined($working->gccontent->pass())){
							$working->gccontent->pass($ret->pass());
						}
					}elsif($name eq "Per base GC content" || $name eq "Per sequence GC content"){
					}elsif($name eq "Per base N content"){
						$fq->ncontent($ret);
						foreach my $k (keys(%{$ret->ns()})){
							if(defined($working->ncontent)){
								$working->ncontent->ns->{$k} = ($working->ncontent->ns()->{$k} + $ret->ns()->{$k}) / 2;
							}else{
								$working->ncontent(pbpncontent->new('ns' => {}));
								$working->ncontent->ns->{$k} = $ret->ns()->{$k};
							}
						}
						if($ret->pass() =~ m/fail/i || !defined($working->ncontent->pass())){
							$working->ncontent->pass($ret->pass());
						}
					}elsif($name eq "Sequence Length Distribution"){
						$fq->seqlens($ret);
						if(defined($working->seqlens())){
							$working->seqlens(($working->seqlens() + $ret) / 2);
						}else{
							$working->seqlens($ret);
						}						
					}elsif($name eq "Sequence Duplication Levels"){
						$fq->seqdups($ret);
						if(defined($working->seqdups())){
							$working->seqdups(($working->seqdups() + $ret) / 2);
						}else{
							$working->seqdups($ret);
						}
					}elsif($name eq "Overrepresented sequences"){
						$fq->overrep($ret);
						$working->overrep($ret) if(!defined($working->overrep()) || $ret =~ m/fail/i);
					}elsif($name eq "Kmer Content"){
						$fq->kmercon($ret);
						$working->kmercon($ret) if(!defined($working->kmercon()) || $ret =~ m/fail/i);
					}elsif($name eq "Basic Statistics"){
						$fq->stats($ret);
						if(defined($working->stats)){
							$working->stats->totseq($working->stats->totseq() + $ret->totseq());
							$working->stats->filtseq($working->stats->filtseq() + $ret->filtseq());
							$working->stats->seqlen(($working->stats->seqlen() + $ret->seqlen()) / 2);
							$working->stats->gc(($working->stats->gc() + $ret->gc()) / 2);
							if($ret->pass() =~ m/fail/i || !defined($working->stats->pass())){
								$working->stats->pass($ret->pass());
							}
						}else{
							$working->stats(stats->new());
							$working->stats->totseq($ret->totseq());
							$working->stats->filtseq($ret->filtseq());
							$working->stats->seqlen($ret->seqlen());
							$working->stats->gc($ret->gc());
							$working->stats->pass($ret->pass());
						}
						
					}
									
				}elsif($inloop){
					$store .= $line;
				}
			}
			
			close IN;
			$filenum++;
			print "Out of $totfiles, finished with:\t$filenum\r";
		} # for files 
		
	} # fastqdir
	foreach my $an (sort{$a cmp $b} keys(%anhash)){
		my $fq = $anhash{$an};

		print SUM $an . "\t" . $fq->stats->totseq() . "\t" . $fq->stats->filtseq() . "\t" . $fq->stats->seqlen() . "\t";
		print SUM $fq->stats->gc() . "\t" . $fq->bpquality->pass();
		
		foreach my $keys(sort {$a cmp $b} keys %{$fq->bpquality->mean()}){
			print SUM "\t" . $fq->bpquality->mean->{$keys};
		}
		# foreach my $keys(%{$fq->pbpquality->median()}){
			# print OUT "\t" . $fq->bpquality->median->{$keys};
		# }
		
		print SUM "\t" . $fq->seqquality->pass();
		foreach my $keys(sort {my($x) = $a =~ /(\d+)-(\d+)/; my($y) = $b =~ /(\d+)-(\d+)/; $x <=> $y}keys %{$fq->seqquality->count()}){
			print SUM "\t" . $fq->seqquality->count->{$keys};
		}
		
		print SUM "\t" . $fq->gccontent->pass();
		foreach my $keys(sort {$a cmp $b}keys %{$fq->gccontent->gc()}){
			print SUM "\t" . $fq->gccontent->gc->{$keys};
		}
		
		print SUM "\t" . $fq->ncontent->pass();
		foreach my $keys(sort {$a cmp $b}keys %{$fq->ncontent->ns()}){
			print SUM "\t" . $fq->ncontent->ns->{$keys};
		}
		
		print SUM "\t" . $fq->seqdups() . "\t" . $fq->overrep() . "\t" . $fq->kmercon() . "\n";
		
	} # for an
	close FDIR;
	
} # Basedir
close SUM;
close DIR;
print "\nFinished\n";

exit;

sub process_module_seg{
	my $input = shift;
		
	my @linesegs = split(/\n/, $input);
	my ($name, $pass) = $linesegs[0] =~ m/>>(.+)\s+(.+)$/;
		
	
	if($name eq "Per base sequence quality"){
		my $pbreturn = pbpquality->new();
		$pbreturn->pass($pass);
		my %meanhash;
		my %medianhash;
		for(my $x = 1; $x < scalar(@linesegs); $x++){
			if($linesegs[$x] =~ /^#/){next;}
			my @msegs = split(/\s+/, $linesegs[$x]);
			if($msegs[0] =~ /(\d+)-(\d+)/){
				for(my $y = $1; $y <= $2; $y++){
					$meanhash{$y} = $msegs[1];
					$medianhash{$y} = $msegs[2];
				}
			}else{
				$meanhash{$msegs[0]} = $msegs[1];
				$medianhash{$msegs[0]} = $msegs[2];
			}
		}
		my @values = sort {$b <=> $a} keys %meanhash;
		my %nmean;
		my %nmedian;
		for(my $x = 1; $x < $values[0] ; $x += int($values[0] * .25 + 0.5)){
			my @meantemp;
			my @mediantemp;
			for(my $y = $x; $y <= $values[0]; $y++){
				push(@meantemp, $meanhash{$y});
				push(@mediantemp, $medianhash{$y});
			}
			my $str = "$x-" . ($x + int($values[0] * .25 + 0.5));
			$nmean{$str} = average(\@meantemp);
			$nmedian{$str} = average(\@mediantemp);
		}
		$pbreturn->mean(\%nmean);
		$pbreturn->median(\%nmedian);
		return $name, $pbreturn;
	}elsif($name eq "Per sequence quality scores"){
		my $psqsreturn = pseqquality->new();
		$psqsreturn->pass($pass);
		my %sqc;
		for(my $x = 1; $x < scalar(@linesegs); $x++){
			if($linesegs[$x] =~ /^#/){next;}
			my @msegs = split(/\s+/, $linesegs[$x]);
			$sqc{$msegs[0]} = $msegs[1];
		}
		my @values = sort {$b <=> $a} keys %sqc;
		my %final;
		for(my $x = 2; $x <= 40 +2; $x+=10){
			my @meanvalue;
			for(my $y = $x; $y <= $values[0]; $y++){
				push(@meanvalue, $sqc{$y}) if exists($sqc{$y});
			}
			my $end = ($x + 10 > 40)? 40 : $x + 10;
			my $str = "$x-$end";
			$final{$str} = average(\@meanvalue);
			if($end == 40){last;}
		}
		$psqsreturn->count(\%final);
		return $name, $psqsreturn;
	}elsif($name eq "Per base sequence content"){
		my $pbreturn = pbpgccontent->new();
		$pbreturn->pass($pass);
		my %meanhash;
		for(my $x = 1; $x < scalar(@linesegs); $x++){
			if($linesegs[$x] =~ /^#/){next;}
			my @msegs = split(/\s+/, $linesegs[$x]);
			if($msegs[0] =~ /(\d+)-(\d+)/){
				for(my $y = $1; $y <= $2; $y++){
					$meanhash{$y} = $msegs[1] + $msegs[2];
				}
			}else{
				$meanhash{$msegs[0]} = $msegs[1] + $msegs[2];
			}
		}
		my @values = sort {$b <=> $a} keys %meanhash;
		my %nmean;
		for(my $x = 1; $x < $values[0]; $x += int($values[0] * .25 + 0.5)){
			my @meantemp;
			for(my $y = $x; $y <= $values[0]; $y++){
				push(@meantemp, $meanhash{$y});
			}
			my $str = "$x-" . ($x + 24);
			$nmean{$str} = average(\@meantemp);
		}
		$pbreturn->gc(\%nmean);
		return $name, $pbreturn;
		
	}elsif($name eq "Per base GC content" || $name eq "Per sequence GC content"){
		return $name, 0;
	}elsif($name eq "Per base N content"){
		my $pbreturn = pbpncontent->new();
		$pbreturn->pass($pass);
		my %meanhash;
		for(my $x = 1; $x < scalar(@linesegs); $x++){
			if($linesegs[$x] =~ /^#/){next;}
			my @msegs = split(/\s+/, $linesegs[$x]);
			if($msegs[0] =~ /(\d+)-(\d+)/){
				for(my $y = $1; $y <= $2; $y++){
					$meanhash{$y} = $msegs[1];
				}
			}else{
				$meanhash{$msegs[0]} = $msegs[1];
			}
		}
		my @values = sort {$b <=> $a} keys %meanhash;
		my %nmean;
		for(my $x = 1; $x < $values[0]; $x += int($values[0] * .25 + 0.5)){
			my @meantemp;
			for(my $y = $x; $y <= $values[0]; $y++){
				push(@meantemp, $meanhash{$y});
			}
			my $str = "$x-" . ($x + 24);
			$nmean{$str} = average(\@meantemp);
		}
		$pbreturn->ns(\%nmean);
		return $name, $pbreturn;
	}elsif($name eq "Sequence Length Distribution"){
		my @lengths;
		for(my $x = 1; $x < scalar(@linesegs); $x++){
			if($linesegs[$x] =~ /^#/){next;}
			my @msegs = split(/\s+/, $linesegs[$x]);
			push(@lengths, "$msegs[0]: $msegs[1]");
		}
		return $name, join(";", @lengths);
	}elsif($name eq "Sequence Duplication Levels"){
		foreach my $v (@linesegs){
			if($v =~ /#Total Duplicate Percentage\s+(\d+)/){
				return $name, $1;
			}
		}
	}elsif($name eq "Overrepresented sequences" || $name eq "Kmer Content"){
		return $name, $pass;
	}elsif($name eq "Basic Statistics"){
		my $stats = stats->new();
		$stats->pass($pass);
		for(my $x = 1; $x < scalar(@linesegs); $x++){
			if($linesegs[$x] =~ /^#/){next;}
			my @msegs = split(/\t/, $linesegs[$x]);
			if($msegs[0] eq "Filename"){
				$stats->file($msegs[1]);
			}elsif($msegs[0] eq "Total Sequences"){
				$stats->totseq($msegs[1]);
			}elsif($msegs[0] eq "Filtered Sequences"){
				$stats->filtseq($msegs[1]);
			}elsif($msegs[0] eq "Sequence length"){
				$stats->seqlen($msegs[1]);
			}elsif($msegs[0] eq "%GC"){
				$stats->gc($msegs[1]);
			}
		}
		return $name, $stats;
	}
	return $name, $pass;
}

sub average{
	my $array_r = shift(@_);
	my @numb = @{$array_r};
	if(scalar(@numb) == 0){
		return 0;
	}
	my $total3 = 0;
	foreach my $num1 (@numb) {
	$total3 += $num1;
	}
	my $mean3 = $total3 / (scalar @numb);
	return $mean3;
}