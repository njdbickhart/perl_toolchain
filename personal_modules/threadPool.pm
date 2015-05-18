#!/usr/bin/perl
# This is meant to be a java-style thread manager with a fixed thread pool

package threadPool;
use Mouse;
use threads;
use threads::shared;
use namespace::autoclean;

has 'queue' => (traits => ['Array'], is => 'rw', isa => 'ArrayRef[Any]', default => sub{[]}, 
	handles => {
		'enqueue' => 'push',
		'remove' => 'splice',
		'thr' => 'get',
		'num' => 'count',
	});
has 'maxthreads' => (is => 'ro', isa => 'Int', required => 1);

sub submit{
	my ($self, $thread, $d) = @_;
	
	$self->enqueue($thread);
	if($self->num >= $self->maxthreads){
		while(1){
			my @joinable = threads->list(threads::joinable);
			my @keep;
			if(scalar(@joinable) > 0){
				for(my $x = 0; $x < $self->num; $x++){
					if($self->thr($x)->is_joinable()){
						$self->thr($x)->join();
						if(defined($d)){
							print STDERR "THREADPOOL joining thread: " . $self->thr($x) . "\n";
						}
					}else{
						push(@keep, $self->thr($x));
					}
				}
				$self->queue(\@keep);
				last;
			}
			sleep(3);
		}
	}
}

sub joinAll{
	my ($self, $d) = @_;
	if(defined($d)){
		my $num = $self->num;
		print STDERR "THREADPOOL waited for $num threads!\n";
	}
	foreach my $thr (@{$self->queue}){
		$thr->join();
	}
}

__PACKAGE__->meta->make_immutable;

1;
