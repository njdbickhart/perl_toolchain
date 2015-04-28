#!/usr/bin/perl
# This module is designed to be an encapsulated interpreter of a fastqc_data.txt summary file
# It will process the file, summarize the statistics and, on demand, provide a suitable tab delimited summary
# The whole goal is to condense the information from fastqc across an entire sample into a spreadsheet-style row

# Required in constructor:
# file => fastq file, sample => sample_name, library => library_name, readNum => first_or_second_read, log => simpleLogger

# Pipeline:
# 	$object->runFastqc($fastqc_exe)
#	$object->parseStats()
#	$object->cleanUp()  <- optional removal of folders and zip file
#	my @headers = $object->getHeaderArray()
#	my @values = $object->getOutArray()

package fastqcParser;
use Mouse;
use namespace::autoclean;
use simpleLogger;
use File::Basename;

my @seHeaders = ("Sample", "Library", "FastqFile", "ReadNum", "TotalReads", "FilteredReads", "ReadLength", "TotalGCPerc", "BPQualityPass", 
	"First25Qual", "Second25Qual", "Third25Qual", "Fourth25Qual", "ReadQSPass", "#ReadsQS2-11", "#ReadsQS12-21", 
	"#ReadsQS22-31", "#ReadsQS32-40", "BPGCContentPass", "First25GC", "Second25GC", "Third25GC", "Fourth25GC", 
	"BPNContentPass", "First25N", "Second25N", "Third25N", "Fourth25N", "SeqDupPerc", "OverRepPass", "KmerPass");

has ['pass', 'overrep', 'kmer'] => (is => 'rw', isa => 'Str', default => 'PASS');
has ['file', 'sample', 'library'] => (is => 'ro', isa => 'Str', required => 1);
has 'log' => (is => 'ro', isa => 'simpleLogger', required => 1);
has 'readNum' => (is => 'ro', isa => 'Num', required => 1);

# Locations of fastqc folder and zip
has ['folder', 'zip'] => (is => 'rw', isa => 'Str');

# Containers for individual stat modules
has 'sampstats' => (is => 'rw', isa => 'sampStats');
has 'pbpquality' => (is => 'rw', isa => 'pbpQuality');
has 'pseqquality' => (is => 'rw', isa => 'pseqQuality');
has 'pbpgccontent' => (is => 'rw', isa => 'pbpGCContent');
has 'pbpncontent' => (is => 'rw', isa => 'pbpNContent');

sub runFastqc{
	my ($self, $fastqc) = @_;
	if(!defined($fastqc)){
		$self->log->Fatal("[FQCPARSE]", "Did not receive fastqc executable!");
	}
	
	$self->log->Info("[FQCPARSE]", "Beginning fastqc runtime on file: " + $self->file);
	
	# Determine file output 
	my ($filename, $dirs, $ext) = fileparse($self->file);
	if($ext eq 'gz'){
		my @fsegs = split(/\./, $filename);
		$self->folder("$dirs/$fsegs[0]" . "_fastqc");
		$self->zip("$dirs/$fsegs[0]" . "_fastqc.zip");
	}else{
		$self->folder("$dirs/$filenmae" . "_fastqc");
		$self->zip("$dirs/$filenmae" . "_fastqc.zip");
	}
	
	system("$fastqc -q " . $self->file);
	$self->log->Info("[FQCPARSE]", "Finished fastqc runtime on file: " + $self->file);
}

sub parseStats{
	my ($self) = @_;
	
	$self->log->Info("[FQCPARSE]", "Parsing fastqc stats on file: " $self->file);
	
	open(IN, "< " . $self->folder . "/fastqc_data.txt") || $self->log->Fatal("[FQCPARSE]", "Error accessing fastqc_data.txt in folder: " . $self->folder . " for file: " . $self->file . "!");
	my $store;
	my $inloop = 0;
	while(my $line = <IN>){
		if($line =~ /^>>/ && $line !~ />>END_MODULE/){
			$inloop = 1;
			$store .= $line;
		}elsif($line =~ />>END_MODULE/){
			$inloop = 0;
			
			# TODO: rework this section to call individual class parsers rather than a blanket sub
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
}


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

__PACKAGE__->meta->make_immutable;

package sampStats;
use Mouse;
use namespace::autoclean;

has ['totseq', 'filtseq', 'seqlen', 'gc'] => (is => 'rw', isa => 'Num', default => 0);

__PACKAGE__->meta->make_immutable;

package pbpQuality;
use Mouse;
use namespace::autoclean;

'pass' => '$',
	'mean' => '%',
	'median' => '%',

__PACKAGE__->meta->make_immutable;

package pseqQuality;
use Mouse;
use namespace::autoclean;

'pass' => '$',
	'count' => '%',

__PACKAGE__->meta->make_immutable;

package pbpGCContent;
use Mouse;
use namespace::autoclean;

'pass' => '$',
	'gc' => '%',

__PACKAGE__->meta->make_immutable;

package pbpNContent;
use Mouse;
use namespace::autoclean;

'pass' => '$',
	'ns' => '%',

__PACKAGE__->meta->make_immutable;

