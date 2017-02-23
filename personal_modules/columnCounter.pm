package columnCounter;
use Mouse;
use FileHandle;
use namespace::autoclean;

# Required attributes
has 'colnum' => (is => 'ro', isa => 'Int', required => 1);
has 'mkdwn' => (is => 'ro', isa => 'Bool', required => 1);
has 'numeric' => (is => 'rw', isa => 'Bool', default => 0);

has 'output' => (is => 'rw', isa => 'Str', predicate => 'has_output');
has 'ignore' => (is => 'rw', isa => 'Str', predicate => 'has_ignore');
has 'counter' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'Keys' => 'keys',
		'Empty' => 'is_empty',
		'Get' => 'get',
		'Set' => 'set',
	});
has 'data' => (traits => ['Array'], is => 'rw', isa => 'ArrayRef[Num]', default => sub{[]},
	handles => {
		'addValue' => 'push',
		'sortOptions' => 'sort',
		'countEntries' => 'count',
	});
has ['sum', 'sumSquares'] => (is => 'rw', isa => 'Num', default => 0);

	
# Reads file to count columns and data
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
	my @data;
	my $sum = 0; my $ss = 0;
	my $com;
	if($self->has_ignore){
		$com = $self->ignore;
	}
	
	while(my $line = <$fh>){
		chomp($line);
		$line =~ s/\r//g;
		$line =~ s/^\s+//;
		if(length($line) < 2 || $line eq ""){
			next;
		}
		
		if($self->has_ignore){
			if($line =~ m/^$com/){
				next;
			}
		}
		
		my @segs = split(/\s{1}/, $line);
		if(scalar(@segs) < $col + 1){next;}
		
		if($self->numeric){
			#Skip non numeric values
			if($segs[$col] !~ /^-?\d+\.?\d*$/){next;}
			push(@data, $segs[$col]);
			$sum += $segs[$col];
			$ss += $segs[$col] * $segs[$col];
		}
		$hash{$segs[$col]} += 1;
	}
	close $fh;
	if(scalar(keys(%hash)) < 1){
		print STDERR "[ColumnCounter] Did not find any entries within the column specified! Please check file contents!\n";
		exit;
	}
	$self->counter(\%hash);
	
	if($self->numeric){
		$self->data(\@data);
		$self->sum($sum);
		$self->sumSquares($ss);
	}
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

sub _numericCount{
	my ($self) = @_;
	my %values;
	
	my %numHash = %{$self->counter};
	my @data = @{$self->sortOptions};
	my $sum = $self->sum;
	my $sqr = $self->sumSquares;
	my $ct = scalar(@data);
	my @sortKeys = sort{$numHash{$b} <=> $numHash{$a}} keys(%numHash);
	
	$values{"Count"} = $ct;
	$values{"Sum"} = $sum;
	$values{"Minimum"} = $data[0];
	$values{"Maximum"} = $data[-1];
	$values{"Average"} = $sum / $ct;
	$values{"Median"} = ($ct % 2 == 0) ? ($data[$ct/2-1] + $data[$ct/2])/2 : $data[floor($ct/2)];
	$values{"Stdev"} = ($ct == 1) ? 0 : sqrt( ($sqr * $ct - $sum * $sum)/($ct * ($ct -1)) );
	$values{"Mode"} = $sortKeys[0];
	
	$self->counter(\%values);
}

# subroutine that writes directly to outfile or stdout
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
	
	my @sortedkeys;
	my $tClass;
	
	if($self->numeric){
		# Process data and organize values in different order
		@sortedkeys = ("Count", "Sum", "Minimum", "Maximum", "Average", "Median", "Stdev", "Mode");
		$self->_numericCount;
		$tClass = "Value";
	}else{
		@sortedkeys = sort{$a cmp $b} $self->Keys;	
		$tClass = "Count";
	}
	
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
		print {$fh} sprintf("\|%-*s\|%*s\|\n", $collen, "Entry", $conlen, "$tClass");
		print {$fh} sprintf("\|\:%s\|%s\:\|\n", $sepstr, $sepcon);
		foreach my $k (@sortedkeys){
			my $kvalue = $k;
			if($k eq ''){
				$kvalue = "<Null>";
			}
			print {$fh} sprintf("\|%-*s\|%*d\|\n", $collen, $kvalue, $conlen, $self->Get($k));
		}
	}else{
		print {$fh} "Entry\t$tClass\n";
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

1;