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
	my ($self, $thread) = @_;
	
	$self->queue($thread);
	if($self->num >= $self->maxthreads){
		while(1){
			my @joinable = $threads->list(threads::joinable);
			my @remove;
			if(scalar(@joinable) > 0){
				for(my $x = 0; $x < $self->num; $x++){
					if($self->thr->is_joinable()){
						$self->thr->join();
						$remove[$x] = 1;
					}else{
						$remove[$x] = 0;
					}
				}
				my @threads = @{$self->queue};
				for(my $x = $self->num - 1; $x >= 0; $x--){
					@threads = splice(@threads, $x, 1);
				}
				$self->queue(\@threads);
				last;
			}
			sleep(3);
		}
	}
}

sub joinAll{
	my ($self) = @_;
	foreach my $thr (@{$self->queue}){
		$thr->join();
	}
}

__PACKAGE__->meta->make_immutable;

1;