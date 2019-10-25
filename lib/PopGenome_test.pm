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
		'config=s','steps=i','flag'
	);
	print "Hey, my config is $opts{config}\n";
	$opts{steps}||= 1;
	print "steps: $opts{steps}\n" if (defined $opts{flag});

	TEST2(\%opts);

	my %var;
	my ($a,$b)=qw(1 2);
	if (($a==2)&&!(defined $var{new})){print "oh\n"}
#	print "$a,hi\n" unless ((! defined $var{new}) && ($var{new} ==3));

}

# This is the only way to pass variables to another subroutine 
sub TEST2{
	my $opts = shift;
	my %opts = %{$opts};
	print "$Bin\n";
	print "triggered TEST2 and show opts $opts{steps}\n";
}
1;

