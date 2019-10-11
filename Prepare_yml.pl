#!/usr/bin/perl -w 

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use YAML::Tiny;

my $usage="
    Usage:
    Options:
        -i <input>   input list_file contain 2 col (or 3 col for PE),name and raw_path(must be fastq format)
        -s <seqtype> input sequence type [SE]|PE
        -o <opath>   output file [./result.yml]
        -j <qsub>    qsub the job [Y] or N run on the nodes directly
        -m <mem>     memory needed[2G]
        -t <thread>  thread needed[4]
        -q <queue>   queue [st.q]
        -P <prj>     Project Number [P18Z10200N0119]
        -h|?Help!
    Example:perl $0 -i name_path.list -j Y 
";
my ($input,$seqtype,$outpath,$job,$mem,$thread,$queue,$prj,$help);
GetOptions(
    '-i=s' => \$input,
    '-s=s' => \$seqtype,
    '-o=s' => \$outpath,
    '-j=s' => \$job,
    '-m=s' => \$mem,
    '-t=i' => \$thread,
    '-q=s' => \$queue,
    '-P=s' => \$prj,
    'h|?'  => \$help,
);


if($help or !$input){die "$usage\n";}
$outpath ||= "./result.yml";
$job ||="Y";
$queue ||= "st.q";
$prj ||= "P18Z10200N0119";
$mem ||= "2G";
$thread ||= 4;
$seqtype ||= "SE";

$config = $Bin/lib/template.yml;
my $yaml = YAML::Tiny->read( $config );
my %cfg = %{$yaml->[0]};

unless (exists $cfg{args}{threads}){$cfg{args}{threads}=$thread}
unless (exists $cfg{args}{prj}){$cfg{args}{prj}=$prj}
unless (exists $cfg{args}{queue}){$cfg{args}{queue}=$queue}

unless (exists $cfg{args}{mem}){$cfg{args}{mem}=$mem}

open IN,"$input" or die $!;
my %sample;
while(<Raw>){
    chomp;
    if($seqtype eq "SE"){
        my($name,$path) = (split /\s+/,$_)[0,1];
        $sample{$name}=1;
        $cfg{fqdata}{$name}{rawdata}{'lib1'}{'Flag'} = "SE";
        $cfg{fqdata}{$name}{rawdata}{'lib1'}{'PL'} = "BGISEQ500";
        $cfg{fqdata}{$name}{rawdata}{'lib1'}{'Phred'} = 33;
        $cfg{fqdata}{$name}{rawdata}{'lib1'}{'fq1'} = $path;
        #$raw{$name} = $path;
    }elsif($seqtype eq "PE"){
        my($name,$path1,$path2) = (split /\s+/,$_)[0,1,2];
        $sample{$name}=1;
  		$cfg{fqdata}{$name}{rawdata}{'lib1'}{'Flag'} = "PE";
        $cfg{fqdata}{$name}{rawdata}{'lib1'}{'PL'} = "BGISEQ500";
        $cfg{fqdata}{$name}{rawdata}{'lib1'}{'Phred'} = 33;
        $cfg{fqdata}{$name}{rawdata}{'lib1'}{'fq1'} = $path1;
        $cfg{fqdata}{$name}{rawdata}{'lib1'}{'fq2'} = $path2;
    }
}
close IN;

foreach my $sample (keys %sample){
	$cfg{population}{$name}{'presumed population'}="All";
}

unless ( -e $outpath){
	# create this yaml object
	$yaml = YAML::Tiny->new( \%cfg );
	# Save both documents to a file
	$yaml->write( $outpath );
}