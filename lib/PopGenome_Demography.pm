package PopGenome_Demography;
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
		'smcpp',
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

	if (defined $opts{smcpp}){ 
		& SMCPP (\%var,\%opts);
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
		
		print IDSH "vcftools --gzvcf $var{vcf} --keep $sfs_outpath/$pop_name.raw.list --recode --stdout >$sfs_outpath/$pop_name/$pop_name.SNP.vcf\n";
		# calculate syn and nonsyn diversity values in each sliding window
		if (defined $opts{nonsyn}){
			print IDSH "vcftools --gzvcf $var{outpath}/snpEff.nonsyn.vcf.gz --keep $sfs_outpath/$pop_name.raw.list --recode --stdout >$sfs_outpath/$pop_name/$pop_name.snpEff.nonsyn.vcf\n";
			print IDSH "vcftools --gzvcf $var{outpath}/snpEff.syn.vcf.gz --keep $sfs_outpath/$pop_name.raw.list --recode --stdout >$sfs_outpath/$pop_name/$pop_name.snpEff.syn.vcf\n";
		}

		### project SFS from the file
		#print IDSH "cp -rf $slidingwindow_outpath/$pop_name.SNP.vcf.gz $sfs_outpath/$pop_name/ && gunzip -f $sfs_outpath/$pop_name/$pop_name.SNP.vcf.gz\n";
		print IDSH "python $Bin/Tools/easySFS.py -i $sfs_outpath/$pop_name/$pop_name.SNP.vcf -p $sfs_outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $sfs_outpath/$pop_name/total_sfs && mv $sfs_outpath/$pop_name/output $sfs_outpath/$pop_name/total_sfs \n";
		if (defined $opts{nonsyn}){
			#print IDSH "cp -rf $slidingwindow_outpath/$pop_name.snpEff.nonsyn.vcf.gz $sfs_outpath/$pop_name/ && gunzip $sfs_outpath/$pop_name/$pop_name.snpEff.nonsyn.vcf.gz\n";
			print IDSH "python $Bin/Tools/easySFS.py -i $sfs_outpath/$pop_name/$pop_name.snpEff.nonsyn.vcf -p $sfs_outpath/$pop_name.list --ploidy $cfg{args}{ploidy} --proj $pop_size -f -a && rm -rf $sfs_outpath/$pop_name/nonsyn_sfs && mv $sfs_outpath/$pop_name/output $sfs_outpath/$pop_name/nonsyn_sfs \n";

			#print IDSH "cp -rf $slidingwindow_outpath/$pop_name.snpEff.syn.vcf.gz $sfs_outpath/$pop_name/ && gunzip $sfs_outpath/$pop_name/$pop_name.snpEff.syn.vcf.gz\n";
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
}

sub SMCPP{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};
	my %pop = %{$var{pop}};
	
	my $smcpp_outpath = "$var{outpath}/SMCPP/";
	if ( !-d "$var{outpath}/SMCPP" ) {make_path "$var{outpath}/SMCPP" or die "Failed to create path: $var{outpath}/SMCPP";}

	if (defined $opts{genome}){
		$var{genome} = PopGenome_Shared::LOADREF($opts{genome});
	}else {
		$var{genome} = PopGenome_Shared::LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	}

	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 6);
		if ( !-d "$smcpp_outpath/$pop_name" ) {make_path "$smcpp_outpath/$pop_name" or die "Failed to create path: $smcpp_outpath/$pop_name";}
		if ( !-d "$smcpp_outpath/$pop_name/out" ) {make_path "$smcpp_outpath/$pop_name/out" or die "Failed to create path: $smcpp_outpath/$pop_name/out";}
		if ( !-d "$smcpp_outpath/$pop_name/analysis" ) {make_path "$smcpp_outpath/$pop_name/analysis" or die "Failed to create path: $smcpp_outpath/$pop_name/analysis";}
		open OT, ">$smcpp_outpath/$pop_name.raw.list"; print OT $pop{$pop_name}{line}; close OT;

		open IN, "$smcpp_outpath/$pop_name.raw.list";
		my $sampletext="";
		while (<IN>){
			chomp;
			$sampletext .= $_;
			$sampletext .= ",";
		}
		close IN;
		$sampletext =~ s/\,$//g;

		my %genome = %{$var{genome}};

		open LIST1, ">$var{shpath}/smcpp.$pop_name.cmd1.list";

		foreach my $scaffold (keys %{$genome{len}}){
			next if ($genome{len}{$scaffold} < 1000000);
			print LIST1 "cd $smcpp_outpath/$pop_name && conda activate smcpp && smc++ vcf2smc --length $genome{len}{$scaffold} $var{vcf} out/$scaffold.smc.gz $scaffold $pop_name:$sampletext && conda deactivate\n";
		}
		close LIST1;
		`perl $Bin/lib/qsub.pl -d $var{shpath}/smcpp.$pop_name.cmd1_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/smcpp.$pop_name.cmd1.list` unless (defined $opts{skipsh});

		open SH, ">$var{shpath}/smcpp.$pop_name.s2.sh";
		print SH "cd $smcpp_outpath/$pop_name\n";
		print SH "conda activate smcpp && smc++ cv -o analysis/ 1.25e-8 out/*.smc.gz\n";
		print SH "smc++ plot plot.pdf analysis/model.final.json\n";
		print SH "conda deactivate\n";
		close SH;

		open LIST2,">$var{shpath}/smcpp.$pop_name.cmd2.list";
		print LIST2 "sh $var{shpath}/smcpp.$pop_name.s2.sh 1>$var{shpath}/smcpp.$pop_name.s2.sh.o 2>$var{shpath}/smcpp.$pop_name.s2.sh.e\n";
		close LIST2;

		`perl $Bin/lib/qsub.pl -d $var{shpath}/smcpp.$pop_name.cmd2_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/smcpp.$pop_name.cmd2.list` unless (defined $opts{skipsh});
	}
}

1;