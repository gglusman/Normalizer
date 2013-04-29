package BaseClass;

use strict;

our $VERSION = '1.00';

sub new {
	my $package = shift;
	
	my $obj = {};
	bless $obj, $package;
	my $status = $obj->initialize(@_);
	unless ($status eq "fail") {
		$obj->set_param("_initialized", 1);
		return $obj;
	}
}

sub ready {
	my($self) = @_;
	
	delete $self->{'_state'};
}

sub save {
	my($self) = @_;
	
	#override this with class-specific save code
	#save code should end with a call to $self->ready() to signal we're up to date
	$self->ready();
}

sub initialize {
	my($self) = shift @_;
	
	$self->{'index'}->{'params'} = 'params';
	$self->{'index'}->{'paramtime'} = 'paramtime';
	$self->{'index'}->{'data'} = 'data';
	$self->{'index'}->{'output'} = 'output';
}

sub param {
	#getter routine - what is a parameter's value?
	#returns the stored value if called with one parameter
	#returns a matching array of references to values, if called with several parameters
	my($self, @fields) = @_;
	my($v, $field, @values);
	
	my $p = $self->{'index'}->{'params'};
	if (1<scalar @fields) {
		foreach $field (@fields) {
			push @values, $self->{$p}->{$field};
		}
		return @values;
	}
	$field = $fields[0];
	return unless $field;
	return unless exists $self->{$p}->{$field};
	if (exists $self->{$p}->{$field}) {
		$v = $self->{$p}->{$field};
		if (wantarray) {
			if (ref($v) eq "HASH") {
				return %{$v};
			} elsif (ref($v) eq "ARRAY") {
				return @{$v};
			}
			#else..? How should we return a REF/CODE/GLOB ref when not in scalar context?
			#What about object references?
			return ($v);
		}
		return $v;
	} elsif (wantarray) {
		return ();
	}
}

sub paramtime {
	#time getter routine - when was a parameter set?
	#returns the stored value if called with one parameter
	#returns a matching array of references to values, if called with several parameters
	my($self, @fields) = @_;
	my($v, $field, @values);
	
	my $p = $self->{'index'}->{'paramtime'};
	if (1<scalar @fields) {
		foreach $field (@fields) {
			push @values, $self->{$p}->{$field};
		}
		return @values;
	}
	$field = $fields[0];
	return unless $field;
	return unless exists $self->{$p}->{$field};
	if (exists $self->{$p}->{$field}) {
		$v = $self->{$p}->{$field};
		if (wantarray) {
			return ($v);
		}
		return $v;
	} elsif (wantarray) {
		return ();
	}
}

sub set_param {
	#setter routine, which could be better designed
	#deletes the stored value if called with one parameter
	#sets to a scalar value if called with two parameters
	#if the new value is equal to the old, no change happens (and the timestamp isn't updated)
	#sets to a list value if called with more than two parameters - and in that case always updates the timestamp
	my($self, $field, @values) = @_;
	
	return unless $field;
	my $v;
	my $p = $self->{'index'}->{'params'};
	my $pt = $self->{'index'}->{'paramtime'};
	if (@values) {
		if (1==scalar @values) {
			if ($self->{$p}->{$field} ne $values[0]) {
				$self->{$p}->{$field} = $values[0];
				$self->{$pt}->{$field} = time();
				if (substr($field,0,1) ne "_") {
					$self->{'_state'} = "changed $field to \"$values[0]\"";
				}
			}
			return $self->{$p}->{$field};
		} else {
			$self->{$pt}->{$field} = time();
			$self->{'state'} = "changed";
			return $self->{$p}->{$field} = \@values;
		}
	} elsif (exists $self->{$p}->{$field}) {
		delete $self->{$p}->{$field};
		$self->{$pt}->{$field} = time();
		if (substr($field,0,1) ne "_") {
			$self->{'_state'} = "changed (deleted) $field";
		}
	}
}

sub flush {
	my($self, @what) = shift @_;
	
	@what = qw/data output/ unless @what;
	foreach my $what (@what) {
		my $section = $self->{'index'}->{$what};
		next unless $section;
		delete $self->{$section};
	}
}

sub DESTROY {
    my $self = shift;
    my $state = $self->{'_state'};
    return unless $state;
    
    $self->save();
    $state = $self->{'_state'};
    return unless $state;
    return unless $self->param("verbose");
    
   	my $id = $self->param("id") || $self->param("method")." on ".$self->param("target")->param("id");
   	open LOGF, ">>destroy.log";
   	print LOGF join("\t", time(), $self->{'_state'}, $self, $id), "\n";
   	close LOGF;
}


1;
__END__

=head1 NAME

BaseClass - Perl extension embodying an object with basic getter/setter functions.

=head1 SYNOPSIS

	use BaseClass;
  
	#create an object and store some info on it
	$object = new BaseClass("some identifier");
	$object->set_param("answer", 42);
	
	#retrieve a parameter
	print "The answer is: ", $object->param("answer"), "\n";
	
	#when was the parameter set?
	print "The parameter 'answer' was set at ", $object->paramtime("answer"), "\n";
	
	#when was the object created?
	print "The object was created at ", $object->paramtime("_initialized"), "\n";
	
	#erase a parameter
	$object->set_param("answer");
	
	#erase the data to save memory
	$object->flush("data");
	
=head1 DESCRIPTION

This is a Perl module providing a base class that:
1) Has standardized getter/setter functions, including timestamping.
2) Has basic support for states.
	$self->ready(): signal we're up-to-date
	$self->save(): called when exiting if we're not up-to-date
3) Has a flush function for releasing memory.

Methods to override: initialize(), save()

=head2 EXPORT

None by default.



=head1 SEE ALSO

http://repeatmasker.org/FEAST

=head1 AUTHOR

Gustavo Glusman, E<lt>Gustavo@SystemsBiology.orgE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 by Gustavo Glusman

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


=cut
