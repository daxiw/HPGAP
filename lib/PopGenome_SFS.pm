package PopGenome_SFS;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use Getopt::Long qw(GetOptionsFromArray);
use FindBin '$Bin';
use YAML::Tiny;
use lib "$Bin/lib";
use PopGenome_Shared;

sub Main{
	my $args = shift; 
	my @args = @{$args};
	my %opts;
	my %var;

	GetOptionsFromArray (\@args, \%opts, 
		'config=s',
		'overwrite',
		'allsteps',
		'samplelist=s',
		'vcf=s',
		'threads',
		'SFS',
		'help',
		'skipsh');

	if (defined $opts{allsteps}){
		$opts{variant_filtering} = 1;
		$opts{intersection} = 1;
	}

	%var = %{PopGenome_Shared::CombineCfg("$Bin/lib/parameter.yml",\$opts,"SFS")}

	$var{genome}= PopGenome_Shared::LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});

	if (defined $opts{SFS}){ 
		& SFS (\%var,\%opts);
	}
}

#---------------------------------------05.DemographicHistory------------------------------------------
sub SFS{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};
	my %pop = %{$var{pop}};
	
	my $sub_outpath = "$var{outpath}/SFS/";if ( !-d "$var{outpath}/SFS" ) {make_path "$var{outpath}/SFS" or die "Failed to create path: $var{outpath}/SFS";}
	my $slidingwindowpath = "$cfg{args}{outdir}/05.IntraPopulation/Slidingwindow/";

	my $var{vcf} = $cfg{step1}{variant_filtering}{vcf};

	my $window_size = $cfg{slidingwindow}{windowsize};
	my $scaffold_number_limit = $cfg{slidingwindow}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{slidingwindow}{scaffold_length_cutoff};

	my $discoal = "/root/discoal/discoal";
	my $diploSHIC = "python /root/diploSHIC/diploSHIC.py";

	open CL1, ">$var{shpath}/SFS_s1.list";
	#loop for each available population
	
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 6);
		if ( !-d "$sub_outpath/$pop_name" ) {make_path "$sub_outpath/$pop_name" or die "Failed to create path: $sub_outpath/$pop_name";}
		open OT, ">$sub_outpath/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;
		my $pop_size = $var{ploidy} * $pop{$pop_name}{count};
		open IDSH, ">$var{shpath}/SFS.$pop_name.sh";
		print IDSH "cd $sub_outpath/$pop_name\n";
		### project SFS from the file
		print IDSH "easySFS.py -i $slidingwindowpath/$pop_name.SNP.vcf.gz -p $sub_outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $sub_outpath/$pop_name/total_sfs && mv $sub_outpath/$pop_name/output $sub_outpath/$pop_name/total_sfs \n";
		print IDSH "easySFS.py -i $slidingwindowpath/$pop_name.snpEff.nonsyn.vcf.gz -p $sub_outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $sub_outpath/$pop_name/nonsyn_sfs && mv $sub_outpath/$pop_name/output $sub_outpath/$pop_name/nonsyn_sfs \n";
		print IDSH "easySFS.py -i $slidingwindowpath/$pop_name.snpEff.syn.vcf.gz -p $sub_outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $sub_outpath/$pop_name/syn_sfs  && mv $sub_outpath/$pop_name/output $sub_outpath/$pop_name/syn_sfs \n";
		print IDSH "cp -f $Bin/lib/1PopBot20Mb.est $sub_outpath/$pop_name/ && cp -f $Bin/lib/1PopBot20Mb.tpl $sub_outpath/$pop_name/ && cp -f $sub_outpath/$pop_name/total_sfs/fastsimcoal2/pop1_MAFpop0.obs $sub_outpath/$pop_name/1PopBot20Mb_MAFpop0.obs\n";
		print IDSH "sed -i `s/18/$pop_size/g` $sub_outpath/$pop_name/1PopBot20Mb_MAFpop0.obs\n";
		for (my $i = 0; $i <10; $i++){
			print IDSH "fsc26 -t 1PopBot20Mb.tpl -e 1PopBot20Mb.est -n 10000 -d -M -L 40 -q -0 -c 40 -B 40 --foldedSFS\ncp -r 1PopBot20Mb 1PopBot20Mb.b$i\n";
		}
		print IDSH "cat $sub_outpath/$pop_name/1PopBot20Mb.b*/1PopBot20Mb.bestlhoods | grep -v NCUR |sort -k1,1n > $sub_outpath/$pop_name/EstimatedNe.list\n";
		close IDSH;
		print CL1 "sh $var{shpath}/SFS.$pop_name.sh 1>$var{shpath}/SFS.$pop_name.sh.o 2>$var{shpath}/SFS.$pop_name.sh.e\n";
	}
	close CL1;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_SFS_s1_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/SFS_s1.list` unless (defined $opts{skipsh});
	

	foreach my $pop_name (keys %pop){
		open SH, ">$var{shpath}/Simulation.$pop_name.sh";
		next unless ($pop{$pop_name}{count} > 6);

		open NELT, "$sub_outpath/$pop_name/EstimatedNe.list";
		my @nelt = <NELT>;
		close NELT;

		my $length = @nelt;
		my @a = split /\s+/, $nelt[$length/2];
		print "$a[0]\n";

		my $m_rate=0.000000025;
		my $Bigwindowsize = 11*$cfg{step4}{slidingwindow}{windowsize};
		my $Ne = $a[0];
		my $Nan = $a[1];
		my $Nbot = $a[2];
		my $Tbot = $a[3];
		my $Tan = $a[6];

		my $Pt_min = 4 * $Ne * $m_rate * $Bigwindowsize/3.16227766;
		my $Pt_max = 4 * $Ne * $m_rate * $Bigwindowsize * 3.16227766;
		my $Pre_min = 4 * $Ne * $m_rate * $Bigwindowsize/3.16227766;
		my $Pre_max = 4 * $Ne * $m_rate * $Bigwindowsize * 3.16227766;
		my $Pa_min = 4 * $Ne * 0.01/3.16227766;
		my $Pa_max = 4 * $Ne * 1 * 3.16227766;
		my $T = 10000/$Ne;

		my $step = 1/11;
		my $init = 0.5/11;
		my $pop_size = $var{ploidy} * $pop{$pop_name}{count};

		open SIMUCL, ">$var{shpath}/Simulation.$pop_name.cmd.list";
		open FVECCL, ">$var{shpath}/SimFvec.$pop_name.cmd.list";

		for (my $i=0; $i<11;$i++){
			my $x = $init + $step*$i;

			# simulate hard sweeps
			print SIMUCL "$discoal $pop_size $cfg{step4}{discoal}{hard_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max, " -Pa ", $Pa_min, " ", $Pa_max, " -Pu 0.000000 0.000040 -ws 0 ";
			print SIMUCL "-en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;

			print SIMUCL " -x $x >$sub_outpath/$pop_name/hard_$i.msOut\n";

			if ($var{ploidy} == 2){
				print FVECCL "$diploSHIC fvecSim diploid hard_$i.msOut hard_$i.fvec --totalPhysLen $Bigwindowsize\n";
			}elsif ($var{ploidy} == 1){
				print FVECCL "$diploSHIC fvecSim haploid hard_$i.msOut hard_$i.fvec --totalPhysLen $Bigwindowsize\n";
			}
			# simulate soft sweeps
			print SIMUCL "$diploSHIC $pop_size $cfg{discoal}{soft_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max, " -Pa ", $Pa_min, " ", $Pa_max, " -Pu 0.000000 0.000040 -ws 0 -Pf 0 0.1 ";
			print SIMUCL "-en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;
			print SIMUCL " -x $x >$sub_outpath/$pop_name/soft_$i.msOut\n";

			if ($var{ploidy} == 2){
				print FVECCL "$diploSHIC fvecSim diploid soft_$i.msOut soft_$i.fvec --totalPhysLen $Bigwindowsize\n";
			} elsif ($var{ploidy} == 1){
				print FVECCL "$diploSHIC fvecSim haploid soft_$i.msOut soft_$i.fvec --totalPhysLen $Bigwindowsize\n";
			}
		}
		print SIMUCL "$discoal $pop_size $cfg{step4}{discoal}{neut_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max;
		print SIMUCL " -en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;
		print SIMUCL " >$sub_outpath/$pop_name/neut.msOut\n";
		print FVECCL "$diploSHIC fvecSim diploid neut.msOut neut.fvec --totalPhysLen $Bigwindowsize\n" if ($var{ploidy} == 2);
		print FVECCL "$diploSHIC fvecSim haploid neut.msOut neut.fvec --totalPhysLen $Bigwindowsize\n" if ($var{ploidy} == 1);
		close SIMUCL;

		print SH "cd $sub_outpath/$pop_name/\n";
		print SH "parallel -j $var{threads} < $var{shpath}/Simulation.$pop_name.cmd.list\n";
		print SH "parallel -j $var{threads} < $var{shpath}/SimFvec.$pop_name.cmd.list\n";
		print SH "mkdir -p $sub_outpath/$pop_name/rawFVFiles && mv $sub_outpath/$pop_name/*.fvec rawFVFiles/\n";
		print SH "mkdir -p $sub_outpath/$pop_name/trainingSets\n";
		print SH "$diploSHIC makeTrainingSets rawFVFiles/neut.fvec rawFVFiles/soft rawFVFiles/hard 5 0,1,2,3,4,6,7,8,9,10 trainingSets/\n";
		print SH "mkdir -p $sub_outpath/$pop_name/updatedSets\n";
		print SH "less trainingSets/linkedHard.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$sub_outpath/$pop_name/updatedSets/linkedHard.fvec\n";
		print SH "less trainingSets/hard.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$sub_outpath/$pop_name/updatedSets/hard.fvec\n";
		print SH "less trainingSets/linkedSoft.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$sub_outpath/$pop_name/updatedSets/linkedSoft.fvec\n";
		print SH "less trainingSets/neut.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$sub_outpath/$pop_name/updatedSets/neut.fvec\n";
		print SH "less trainingSets/soft.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$sub_outpath/$pop_name/updatedSets/soft.fvec\n";
		print SH "python /root/diploSHIC/diploSHIC.py train updatedSets/ updatedSets/ bfsModel\n";
		print SH "mkdir -p $sub_outpath/$pop_name/observedFVFiles && cp $slidingwindowpath/*.$pop_name.SNP.fvec $sub_outpath/$pop_name/observedFVFiles/\n";

		my $i = 1;
		open PREDCL, ">$var{shpath}/$pop_name.predict.cmd.list";
		foreach my $id(sort { $var{genome}->{len}{$b} <=> $var{genome}->{len}{$a} } keys %{$var{genome}->{len}}){
			if (($var{genome}->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){
				my $len=$var{genome}->{len}{$id};
				print PREDCL "python /root/diploSHIC/diploSHIC.py predict bfsModel.json bfsModel.weights.hdf5 $sub_outpath/$pop_name/observedFVFiles/$id.$pop_name.SNP.fvec $sub_outpath/$pop_name/observedFVFiles/$id.$pop_name.SNP.preds\n";
				$i++;
			}
		}
		close PREDCL;
		print SH "parallel -j $cfg{args}{threads} < $var{shpath}/$pop_name.predict.cmd.list\n";
		close SH;
		`sh $var{shpath}/Simulation.$pop_name.sh 1>$var{shpath}/Simulation.$pop_name.sh.o 2>$var{shpath}/Simulation.$pop_name.sh.e` unless ($skipsh ==1);
	}
}

1;