package PopGenome_arrange;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use Getopt::Long qw(GetOptionsFromArray);
use YAML::Tiny;
use FindBin '$Bin';
use lib "$Bin/lib";
use PopGenome_Shared;

sub ARRANGE{
	my $args = @_;
	my %opts;
	GetOptionsFromArray($args, \%opts, 'config=i');
	print "Hey, my config is $opts{config}\n";

	my ($config, $step, $skipsh, $help, $run_flag);

	#------------------------------------- default parameters ---------------------------------------
	#$filter ||= "-l 15 -m 3 -p ATGC,10 -n 1 -z";
	$run ||= '';
	$skipsh ||= 0;
	$step ||= '0:indexing;1:read_filtering,read_mapping,read_report,recalibration,variant_calling,combine_calling,variant_filtering;3:phylogeny,admixture;4:homozygosity,roh,ld,slidingwindow,sfs';
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
	#read in and add other settings into the configuration file
	my $yaml = YAML::Tiny->read( $config );
	my %cfg = %{$yaml->[0]};
	unless (exists $cfg{args}{threads}){$cfg{args}{threads}=40}

	my $allcfg = "$cfg{args}{outdir}/allcfg.yml";

	unless ( -e $allcfg){
		# create this yaml object
		$yaml = YAML::Tiny->new( \%cfg );
		# Save both documents to a file
		$yaml->write( $allcfg );
	}

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
	# when the analysis is not specified
		# set the running environment for udocker 
	my $main = "$cfg{args}{outdir}/hpgap.main.sh";
	open MH, ">$main"; print MH "#!/bin/sh\ncd $cfg{args}{outdir}\n";
	my $udocker_cmd="udocker run ";
		for (my $i=0;$i<@{$cfg{args}{mount}};$i++){
			if (exists $cfg{args}{mount}->[$i]->{host_tmp}){
				$udocker_cmd .= "-v $cfg{args}{mount}->[$i]->{host_tmp}:/tmp ";
			}elsif (exists $cfg{args}{mount}->[$i]->{host_path}){
				$udocker_cmd .= "-v $cfg{args}{mount}->[$i]->{host_path}:$cfg{args}{mount}->[$i]->{host_path} ";
			}
		}
		$udocker_cmd .= "--env=\"$cfg{args}{env}\" $cfg{args}{container} /bin/bash -c ";
		
	#01.indexing
	print MH "time $udocker_cmd 'hpgap --run indexing --config $allcfg' && echo 'indexing done'\n" if (exists $step{0}{indexing});

	print MH "time $udocker_cmd 'hpgap --run read_filtering --config $allcfg'  && echo 'read_filtering done'\n" if (exists $step{1}{read_filtering});

	print MH "time $udocker_cmd 'hpgap --run read_mapping --config $allcfg' && echo 'read_mapping done'\n" if (exists $step{1}{read_mapping});
	
	print MH "time $udocker_cmd 'hpgap --run read_report --config $allcfg' && echo 'read_report done'\n" if (exists $step{1}{read_report});

	print MH "time $udocker_cmd 'hpgap --run variant_calling --config $allcfg' && echo 'variant_calling done'\n" if (exists $step{1}{variant_calling});

	print MH "time $udocker_cmd 'hpgap --run combine_calling --config $allcfg' && echo 'combine_calling done'\n" if (exists $step{1}{combine_calling});

	print MH "time $udocker_cmd 'hpgap --run variant_filtering --config $allcfg' && echo 'variant_filtering done'\n" if (exists $step{1}{variant_filtering});


	#print "$udocker_cmd 'hpgap --run step1_comparison --config $allcfg'\n" if (exists $step{1}{Comparison});
	#PopGenome::VARIANT_COMPARISON($allcfg,$skipsh) if ($run eq 'step1_comparison');
	#############################################
	#											#
	#	Require statistics of variants	
	#											#
	#############################################

	#02.SampleFiltering
	print MH "time $udocker_cmd 'hpgap --run relatedness --config $allcfg' && echo 'relatedness done'\n" if (exists $step{2}{relatedness});

	print MH "time $udocker_cmd 'hpgap --run phylogeny --config $allcfg' && echo 'phylogeny done'\n" if (exists $step{3}{phylogeny});

	print MH "time $udocker_cmd 'hpgap --run admixture --config $allcfg' && echo 'admixture done'\n" if (exists $step{3}{admixture});

	print MH "time $udocker_cmd 'hpgap --run homozygosity --config $allcfg' && echo 'homozygosity done'\n" if (exists $step{5}{homozygosity});

	print MH "time $udocker_cmd 'hpgap --run roh --config $allcfg' && echo 'roh done'\n" if (exists $step{4}{roh});

	print MH "time $udocker_cmd 'hpgap --run ld --config $allcfg' && echo 'ld done'\n" if (exists $step{4}{ld});

	#05.IntraPopulation
	print MH "time $udocker_cmd 'hpgap --run slidingwindow --config $allcfg' && echo 'slidingwindow done'\n" if (exists $step{4}{slidingwindow});

	print MH "time $udocker_cmd 'hpgap --run sfs --config $allcfg' && echo 'sfs done'\n" if (exists $step{4}{sfs});

	#06.Selection
	print MH "time $udocker_cmd 'hpgap --run mktest --config $allcfg' && echo 'MKtest done'\n" if (exists $step{6}{MKtest});
	close MH;
}
