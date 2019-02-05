#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long;
use Cwd qw(getcwd abs_path);
use FindBin '$Bin';
use lib "$Bin/lib";
use lib '/root/miniconda3/lib/site_perl/5.26.2';
use PopGenome;
use YAML::Tiny;

my ($config, $step, $run, $skipsh, $help);

GetOptions (
	"config=s" => \$config,
	"step=s" => \$step, #string
	"run=s" => \$run,
	"skipsh" => \$skipsh,
	"help"  => \$help)   # flag
or die ("Error in command line arguments\n");

#------------------------------------- default parameters ---------------------------------------
#$filter ||= "-l 15 -m 3 -p ATGC,10 -n 1 -z";
$run ||= '';
$skipsh ||= 0;
$step ||= '0:indexing;1:read_filtering,read_mapping,recalibration,variant_calling,combine_calling,variant_filtering;3:phylogeny,admixture;4:homozygosity,roh,ld,slidingwindow,sfs';
#--------------------------------------- get step options ----------------------------------------
#print "$verbose\n";
if ( defined ($help) ) {
	&usage;
	exit;
}	
unless ( defined ($config) ) {
	&usage;
	exit;
}	
#my $yml_file = shift;
my $yaml = YAML::Tiny->read( $config );
my %cfg = %{$yaml->[0]};

`cp -f $config $cfg{args}{outdir}/allcfg.yml` unless ( -e "$cfg{args}{outdir}/allcfg.yml");
my $allcfg = "$cfg{args}{outdir}/allcfg.yml";

my %step;
if(defined $step){
	my @step = split /;/, $step;
	foreach my $sopt(@step){
		chomp $sopt;
		my $sid = $1 if ($sopt =~ /^(\d+)\:/);
		$sopt =~ s/^\d+\://g;
		my @a = split /,/, $sopt;
		foreach (@a){
			my $opt = $1 if(/(\w+)/);
			$step{$sid}{$opt}=1 if (defined $opt);
		}
	}
}

#-------------------------------------- Main steps  --------------------------------------------------
my $main = "$cfg{args}{outdir}/HPGAP.main.sh";
open MH, ">$main"; print MH "#!/bin/sh\ncd $cfg{args}{outdir}\n";
my $udocker_cmd="udocker run ";
	for (my $i=0;$i<@{$cfg{args}{mount}};$i++){
		$udocker_cmd .= "-v $cfg{args}{mount}->[$i]->{hostpath}:$cfg{args}{mount}->[$i]->{dockerpath} ";
	}
	$udocker_cmd .= "--env=\"$cfg{args}{env}\" $cfg{args}{container} /bin/bash -c ";
	
#01.indexing
print MH "$udocker_cmd 'HPGAP.pl --run step0_indexing --config $allcfg'\n" if (exists $step{0}{indexing});
PopGenome::INDEXING($allcfg,$skipsh) if ($run eq 'step0_indexing');

print MH "$udocker_cmd 'HPGAP.pl --run step1_read_filtering --config $allcfg'\n" if (exists $step{1}{read_filtering});
PopGenome::DATA_FILTERING($allcfg,$skipsh) if ($run eq 'step1_read_filtering');

print MH "$udocker_cmd 'HPGAP.pl --run step1_read_mapping --config $allcfg'\n" if (exists $step{1}{read_mapping});
PopGenome::MAPPING($allcfg,$skipsh) if ($run eq 'step1_read_mapping');

print MH "$udocker_cmd 'HPGAP.pl --run step1_recalibration --config $allcfg'\n" if (exists $step{1}{recalibration});
PopGenome::CALIBRATION($allcfg,$skipsh) if ($run eq 'step1_recalibration');

print MH "$udocker_cmd 'HPGAP.pl --run step1_variant_calling --config $allcfg'\n" if (exists $step{1}{variant_calling});
PopGenome::VARIANT_CALLING($allcfg,$skipsh) if ($run eq 'step1_variant_calling');

print MH "$udocker_cmd 'HPGAP.pl --run step1_combine_calling --config $allcfg'\n" if (exists $step{1}{combine_calling});
PopGenome::COMBINE_CALLING($allcfg,$skipsh) if ($run eq 'step1_combine_calling');

print MH "$udocker_cmd 'HPGAP.pl --run step1_variant_filtering --config $allcfg'\n" if (exists $step{1}{variant_filtering});
PopGenome::VARIANT_FILTERING($allcfg,$skipsh) if ($run eq 'step1_variant_filtering');


#print "$udocker_cmd 'HPGAP.pl --run step1_comparison --config $allcfg'\n" if (exists $step{1}{Comparison});
#PopGenome::VARIANT_COMPARISON($allcfg,$skipsh) if ($run eq 'step1_comparison');
#############################################
#											#
#	Require statistics of variants	
#											#
#############################################

#02.SampleFiltering
print MH "$udocker_cmd 'HPGAP.pl --run step2_relatedness --config $allcfg'\n" if (exists $step{2}{relatedness});
PopGenome::RELATEDNESS($allcfg,$skipsh) if ($run eq 'step2_relatedness');

print MH "$udocker_cmd 'HPGAP.pl --run step3_phylogeny --config $allcfg'\n" if (exists $step{3}{phylogeny});
PopGenome::PHYLOGENY($allcfg,$skipsh) if ($run eq 'step3_phylogeny');

print MH "$udocker_cmd 'HPGAP.pl --run step3_admixture --config $allcfg'\n" if (exists $step{3}{admixture});
PopGenome::ADMIXTURE($allcfg,$skipsh) if ($run eq 'step3_admixture');

print MH "$udocker_cmd 'HPGAP.pl --run step4_homozygosity --config $allcfg'\n" if (exists $step{5}{homozygosity});
PopGenome::HOMOZYGOSITY($allcfg,$skipsh) if ($run eq 'step4_homozygosity');

print MH "$udocker_cmd 'HPGAP.pl --run step4_roh --config $allcfg'\n" if (exists $step{4}{roh});
PopGenome::ROH($allcfg,$skipsh) if ($run eq 'step4_roh');

print MH "$udocker_cmd 'HPGAP.pl --run step4_ld --config $allcfg'\n" if (exists $step{4}{ld});
PopGenome::LD($allcfg,$skipsh) if ($run eq 'step4_ld');

#05.IntraPopulation
print MH "$udocker_cmd 'HPGAP.pl --run step4_diversity --config $allcfg'\n" if (exists $step{4}{slidingwindow});
PopGenome::SLIDINGWINDOW($allcfg,$skipsh) if ($run eq 'step4_slidingwindow');

print MH "$udocker_cmd 'HPGAP.pl --run step4_sfs --config $allcfg'\n" if (exists $step{4}{sfs});
PopGenome::SFS($allcfg,$skipsh) if ($run eq 'step4_sfs');

#06.Selection
print MH "$udocker_cmd 'HPGAP.pl --run step6_mktest --config $allcfg'\n" if (exists $step{6}{MKtest});
PopGenome::MKTEST($allcfg,$skipsh) if ($run eq 'step6_mktest');



close MH;
#----------------------------------- usage sub progamm ------------------------------------------
sub usage{
	print STDERR <<USAGE;

Description
	For Popolation genetic analysis in helminths.

Version 
	06 Feb 2019: Version v1.0.0

Author
	Daxi Wang
	email: darcywdx\@gmail.com
	please contact me if you find any bug.
							   	
Usage
	HPGAP.pl --samplelist sample.list --reference reference.fa [-options]
	
	--run <String> use this option choose one of steps below (the option should not be used with --step at the same time)
		step0_indexing
		step1_read_filtering
		step1_read_mapping
		step1_recalibration
		step1_variant_calling
		step1_combine_calling
		step1_variant_filtering
		step2_relatedness
		step3_phylogeny
		step3_admixture
		step4_homozygosity
		step4_roh
		step4_ld
		step4_slidingwindow
		step4_sfs

	--config path to the .yml config file
	--step <String>	specified steps , separated by commas (e.g., "0:A;1:A,B,C,E,F,G").
		0:indexing;
		1:read_filtering,read_mapping,recalibration,variant_calling,combine_calling,variant_filtering;
		3:phylogeny,admixture;
		5:homozygosity,roh,ld,slidingwindow,sfs
	--skipsh use this to skip running bash inside the pipeline
	--help

Note 
	1. No symbolic link to the files outside the mounted volume, which means all the data files themselves should be located within the mounted volume.
	2. For each pair of fastq data, the second colomn (library or flowcell code) should always be unique.
	3. All the paths need to be absolute path

Example 
	To be done
USAGE
}
