package PopGenome_test;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use Getopt::Long qw(GetOptionsFromArray);
use YAML::Tiny;
use FindBin '$Bin';
use lib "$Bin/lib";
use PopGenome_Shared;

sub TEST{
	my $args = shift; 
	my @args = @{$args};
	my %opts;

	GetOptionsFromArray(\@args, \%opts, 
		'config=s','steps=i','flag');
	print "Hey, my config is $opts{config}\n";
	$opts{steps}||= 1;
	print "steps: $opts{steps}\n" if (defined $opts{flag});

	TEST2(\%opts);

	sub TEST3{
		print "triggered TEST3 and show opts $opts{steps}\n";
	}
	TEST3();

}

# This is the only way to pass variables to another subroutine 
sub TEST2{
	my $opts = shift;
	my %opts = %{$opts};
	print "triggered TEST2 and show opts $opts{steps}\n";
}
1;

