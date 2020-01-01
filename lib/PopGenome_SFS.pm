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
		'genome=s',
		'threads',
		'SFS',
		'help',
		'skipsh');

	if (defined $opts{allsteps}){
		$opts{variant_filtering} = 1;
		$opts{intersection} = 1;
	}

	%var = %{PopGenome_Shared::CombineCfg("$Bin/lib/parameter.yml",\%opts,"DemographicHistory")};

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
	
	my $sfs_outpath = "$var{outpath}/SFS/";if ( !-d "$var{outpath}/SFS" ) {make_path "$var{outpath}/SFS" or die "Failed to create path: $var{outpath}/SFS";}
	my $slidingwindow_outpath = "$var{outpath}/../Slidingwindow/";

	if (defined $opts{genome}){
		$var{genome} = PopGenome_Shared::LOADREF($opts{genome});
	}else {
		$var{genome} = PopGenome_Shared::LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	}

	my $window_size = $cfg{sfs}{windowsize};
	my $scaffold_number_limit = $cfg{sfs}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{sfs}{scaffold_length_cutoff};

	my $path = "/ldfssz1/ST_INFECTION/P17Z10200N0536_Echinococcus/USER/wangdaxi/Github/";
	#my $discoal = "/root/discoal/discoal";
	my $discoal = "$path/discoal/discoal";
	#my $diploSHIC = "python /root/diploSHIC/diploSHIC.py";
	my $diploSHIC = "python $path/diploSHIC.py";

	open CL1, ">$var{shpath}/SFS_s1.list";
	#loop for each available population
	
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 6);
		if ( !-d "$sfs_outpath/$pop_name" ) {make_path "$sfs_outpath/$pop_name" or die "Failed to create path: $sfs_outpath/$pop_name";}
		open OT, ">$sfs_outpath/$pop_name.raw.list";
		print OT $pop{$pop_name}{line};
		close OT;

		open IN, "$sfs_outpath/$pop_name.raw.list";
		open OT, ">$sfs_outpath/$pop_name.list";
		while (<IN>){
			chomp;
			print OT "$_\t$pop_name\n";
		}
		close IN;
		close OT;
		
		my $pop_size = $var{ploidy} * $pop{$pop_name}{count};
		open IDSH, ">$var{shpath}/SFS.$pop_name.sh";
		print IDSH "cd $sfs_outpath/$pop_name\n";
		### project SFS from the file
		print IDSH "cp -rf $slidingwindow_outpath/$pop_name.SNP.vcf.gz $sfs_outpath/$pop_name/ && gunzip -f $sfs_outpath/$pop_name/$pop_name.SNP.vcf.gz\n";
		print IDSH "python $Bin/Tools/easySFS.py -i $sfs_outpath/$pop_name/$pop_name.SNP.vcf -p $sfs_outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $sfs_outpath/$pop_name/total_sfs && mv $sfs_outpath/$pop_name/output $sfs_outpath/$pop_name/total_sfs \n";
		if (defined $opts{nonsyn}){
			print IDSH "cp -rf $slidingwindow_outpath/$pop_name.snpEff.nonsyn.vcf.gz $sfs_outpath/$pop_name/ && gunzip $sfs_outpath/$pop_name/$pop_name.snpEff.nonsyn.vcf.gz\n";
			print IDSH "python $Bin/Tools/easySFS.py -i $sfs_outpath/$pop_name/$pop_name.snpEff.nonsyn.vcf -p $sfs_outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $sfs_outpath/$pop_name/nonsyn_sfs && mv $sfs_outpath/$pop_name/output $sfs_outpath/$pop_name/nonsyn_sfs \n";

			print IDSH "cp -rf $slidingwindow_outpath/$pop_name.snpEff.syn.vcf.gz $sfs_outpath/$pop_name/ && gunzip $sfs_outpath/$pop_name/$pop_name.snpEff.syn.vcf.gz\n";
			print IDSH "python $Bin/Tools/easySFS.py -i $sfs_outpath/$pop_name/$pop_name.snpEff.syn.vcf -p $sfs_outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $sfs_outpath/$pop_name/syn_sfs  && mv $sfs_outpath/$pop_name/output $sfs_outpath/$pop_name/syn_sfs \n";
		}
		print IDSH "cp -f $Bin/lib/1PopBot20Mb.est $sfs_outpath/$pop_name/$pop_name.est && cp -f $Bin/lib/1PopBot20Mb.tpl $sfs_outpath/$pop_name/$pop_name.tpl && cp -f $sfs_outpath/$pop_name/total_sfs/fastsimcoal2/$pop_name\_MAFpop0.obs $sfs_outpath/$pop_name/$pop_name\_MAFpop0.obs\n";
		print IDSH "sed -i `s/18/$pop_size/g` $sfs_outpath/$pop_name/$pop_name.tpl\n";
		
		for (my $i = 0; $i <10; $i++){
			print IDSH "rm -rf $pop_name.b$i\n";
			print IDSH "fsc26 -t $pop_name.tpl -e $pop_name.est -n 10000 -d -M -L 40 -q -0 -c 40 -B 40 --foldedSFS\ncp -r $pop_name $pop_name.b$i\n";
		}

		print IDSH "cat $sfs_outpath/$pop_name/$pop_name.b*/$pop_name.bestlhoods | grep -v NCUR |sort -k1,1n > $sfs_outpath/$pop_name/EstimatedNe.list\n";
		close IDSH;
		print CL1 "sh $var{shpath}/SFS.$pop_name.sh 1>$var{shpath}/SFS.$pop_name.sh.o 2>$var{shpath}/SFS.$pop_name.sh.e\n";
	}
	close CL1;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_SFS_s1_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/SFS_s1.list` unless (defined $opts{skipsh});
	
=head
	foreach my $pop_name (keys %pop){
		open SH, ">$var{shpath}/Simulation.$pop_name.sh";
		next unless ($pop{$pop_name}{count} > 6);

		open NELT, "$sfs_outpath/$pop_name/EstimatedNe.list";
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

			print SIMUCL " -x $x >$sfs_outpath/$pop_name/hard_$i.msOut\n";

			if ($var{ploidy} == 2){
				print FVECCL "$diploSHIC fvecSim diploid hard_$i.msOut hard_$i.fvec --totalPhysLen $Bigwindowsize\n";
			}elsif ($var{ploidy} == 1){
				print FVECCL "$diploSHIC fvecSim haploid hard_$i.msOut hard_$i.fvec --totalPhysLen $Bigwindowsize\n";
			}
			# simulate soft sweeps
			print SIMUCL "$diploSHIC $pop_size $cfg{discoal}{soft_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max, " -Pa ", $Pa_min, " ", $Pa_max, " -Pu 0.000000 0.000040 -ws 0 -Pf 0 0.1 ";
			print SIMUCL "-en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;
			print SIMUCL " -x $x >$sfs_outpath/$pop_name/soft_$i.msOut\n";

			if ($var{ploidy} == 2){
				print FVECCL "$diploSHIC fvecSim diploid soft_$i.msOut soft_$i.fvec --totalPhysLen $Bigwindowsize\n";
			} elsif ($var{ploidy} == 1){
				print FVECCL "$diploSHIC fvecSim haploid soft_$i.msOut soft_$i.fvec --totalPhysLen $Bigwindowsize\n";
			}
		}
		print SIMUCL "$discoal $pop_size $cfg{step4}{discoal}{neut_simulation_times} $Bigwindowsize -Pt ", $Pt_min, " ", $Pt_max, " -Pre ", $Pt_min, " ", $Pt_max;
		print SIMUCL " -en ", $Tbot/$Ne, " 0 ", $Nbot/$Ne, " -en ", $Tan/$Ne, " 0 ", $Nan/$Ne;
		print SIMUCL " >$sfs_outpath/$pop_name/neut.msOut\n";
		print FVECCL "$diploSHIC fvecSim diploid neut.msOut neut.fvec --totalPhysLen $Bigwindowsize\n" if ($var{ploidy} == 2);
		print FVECCL "$diploSHIC fvecSim haploid neut.msOut neut.fvec --totalPhysLen $Bigwindowsize\n" if ($var{ploidy} == 1);
		close SIMUCL;

		print SH "cd $sfs_outpath/$pop_name/\n";
		print SH "parallel -j $var{threads} < $var{shpath}/Simulation.$pop_name.cmd.list\n";
		print SH "parallel -j $var{threads} < $var{shpath}/SimFvec.$pop_name.cmd.list\n";
		print SH "mkdir -p $sfs_outpath/$pop_name/rawFVFiles && mv $sfs_outpath/$pop_name/*.fvec rawFVFiles/\n";
		print SH "mkdir -p $sfs_outpath/$pop_name/trainingSets\n";
		print SH "$diploSHIC makeTrainingSets rawFVFiles/neut.fvec rawFVFiles/soft rawFVFiles/hard 5 0,1,2,3,4,6,7,8,9,10 trainingSets/\n";
		print SH "mkdir -p $sfs_outpath/$pop_name/updatedSets\n";
		print SH "less trainingSets/linkedHard.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$sfs_outpath/$pop_name/updatedSets/linkedHard.fvec\n";
		print SH "less trainingSets/hard.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$sfs_outpath/$pop_name/updatedSets/hard.fvec\n";
		print SH "less trainingSets/linkedSoft.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$sfs_outpath/$pop_name/updatedSets/linkedSoft.fvec\n";
		print SH "less trainingSets/neut.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$sfs_outpath/$pop_name/updatedSets/neut.fvec\n";
		print SH "less trainingSets/soft.fvec|perl -ne '\@a = split /\\t/; \$i++;\$l=\@a;if (\$l == 132){print;}' >$sfs_outpath/$pop_name/updatedSets/soft.fvec\n";
		print SH "$diploSHIC train updatedSets/ updatedSets/ bfsModel\n";
		print SH "mkdir -p $sfs_outpath/$pop_name/observedFVFiles && cp $slidingwindow_outpath/*.$pop_name.SNP.fvec $sfs_outpath/$pop_name/observedFVFiles/\n";

		my $i = 1;
		open PREDCL, ">$var{shpath}/$pop_name.predict.cmd.list";
		foreach my $id(sort { $var{genome}->{len}{$b} <=> $var{genome}->{len}{$a} } keys %{$var{genome}->{len}}){
			if (($var{genome}->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){
				my $len=$var{genome}->{len}{$id};
				print PREDCL "python /root/diploSHIC/diploSHIC.py predict bfsModel.json bfsModel.weights.hdf5 $sfs_outpath/$pop_name/observedFVFiles/$id.$pop_name.SNP.fvec $sfs_outpath/$pop_name/observedFVFiles/$id.$pop_name.SNP.preds\n";
				$i++;
			}
		}
		close PREDCL;
		print SH "parallel -j $cfg{args}{threads} < $var{shpath}/$pop_name.predict.cmd.list\n";
		close SH;
		`sh $var{shpath}/Simulation.$pop_name.sh 1>$var{shpath}/Simulation.$pop_name.sh.o 2>$var{shpath}/Simulation.$pop_name.sh.e` unless ($skipsh ==1);
	}

=cut
}

1;