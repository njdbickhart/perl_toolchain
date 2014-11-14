#!/usr/bin/perl
# A collection of useful tools to manipulate and bin bed files
# Shamelessly copied from James Kent and 

package kentBinTools;

sub new {
	my $class = shift();
	return bless {}, $class;
}

sub most {
	my($self, $a, $b) = @_;
	return ($a > $b)? $a : $b;
}
sub least {
	my ($self, $a, $b) = @_;
	return ($a < $b) ? $a : $b;
}
sub overlap {
	my ($self, $s1, $e1, $s2, $e2) = @_;
	return ($self->least($e1, $e2) - $self->most($s1, $s2));
}
sub searchbins{
	my ($self, $start, $end) = @_;
	my $genomiclength = 536870912;
	my %tmph;
	my $href = _collectbins($start, $end, 0,0,0,0,1,0,$genomiclength, \%tmph);
	my @vals = keys(%{$href});
	return @vals;
}
sub getbin{
	my ($self, $start, $end) = @_;
	$start *= 1;
	$end *= 1;
	my $genomiclength = 536870912;
	return _calcbin($start, $end, 0,0,0,0,1,0,$genomiclength);
}

# private method
sub _collectbins{
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
		_collectbins($start, $end, $firstchildbin + $i, $level + 1, $childrowbinstart, $firstchildindex + $i, $childbinrowcount, $childstart, $childlength, $hset);
	}
	return $hset;
}

# private method
sub _calcbin{
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
			my $n = _calcbin($start, $end, $firstchildbin+$i, $level + 1, $childrowbinstart, $firstchildindex + $i, $childbinrowcount, $genomicpos + $i * $childlength, $childlength);
			if($n != -1){
				return $n;
			}
		}
		return $binid;
	}
	return -1;
}
1;
