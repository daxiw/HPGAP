package PopGenome_Homozygosity;
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
		'homozygosity',
		'roh',
		'ld',
		'help',
		'skipsh');

	if (defined $opts{allsteps}){
		$opts{variant_filtering} = 1;
		$opts{intersection} = 1;
	}

	%var = %{PopGenome_Shared::CombineCfg("$Bin/lib/parameter.yml",\%opts,"Homozygosity")};

	if (defined $opts{homozygosity}){ & HOMOZYGOSITY (\%var,\%opts);}
	if (defined $opts{roh}){ & ROH (\%var,\%opts);}
	if (defined $opts{ld}){ & LD (\%var,\%opts);}
}

#################################
#			   #
#   	step4_Homozygosity   	#
#			   #
#################################
sub HOMOZYGOSITY{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};
	my %pop = %{$var{pop}};

	if ( !-d "$var{outpath}/Homozygosity" ) {make_path "$var{outpath}/Homozygosity" or die "Failed to create path: $var{outpath}/Homozygosity";}

	open CL, ">$var{shpath}/cmd_Homozygosity.list";
	
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 4);

		open OT, ">$var{outpath}/Homozygosity/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;

		open SH, ">$var{shpath}/Homozygosity.$pop_name.sh";
		print SH "cd $var{outpath}/Homozygosity\n";
		print SH "vcftools --gzvcf $var{vcf} --keep $var{outpath}/Homozygosity/$pop_name.list --het --stdout ",'|perl -ne \'if(/INDV/){print "INDV\tO_HOM\tE_HOM\tN_SITES\tPERCENTAGE\tF\n";}else{@a=split /\t/;$p=$a[1]/$a[3]*100;$p=sprintf("%.2f", $p);print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$p\t$a[4]";}\'', ">$pop_name.het.list\n";
		close SH;
		print CL "sh $var{shpath}/Homozygosity.$pop_name.sh 1>$var{shpath}/Homozygosity.$pop_name.sh.o 2>$var{shpath}/Homozygosity.$pop_name.sh.e \n";
	}	
	
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_Homozygosity_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/cmd_Homozygosity.list` unless (defined $opts{skipsh});

}

sub ROH{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};
	my %pop = %{$var{pop}};
	
	if (defined $opts{genome}){
		$var{genome} = PopGenome_Shared::LOADREF($opts{genome});
	}else {
		$var{genome} = PopGenome_Shared::LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	}

	if ( !-d "$var{outpath}/ROH" ) {make_path "$var{outpath}/ROH" or die "Failed to create path: $var{outpath}/ROH";}

	my $window_size = $cfg{ROH}{windowsize};
	my $scaffold_number_limit = $cfg{ROH}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{ROH}{scaffold_length_cutoff};

	### load reference and map scaffold to number list 
	my $i=1;my $j=0;
	open OT, ">$var{outpath}/ROH/chr_map.list";
	foreach my $id(sort { $var{genome}->{len}{$b} <=> $var{genome}->{len}{$a} } keys %{$var{genome}->{len}}){
		print OT "$id\t$i\n";
		$i++;
		if (($var{genome}->{len}{$id} >= $cfg{ROH}{scaffold_length_cutoff})&&($j < $cfg{ROH}{scaffold_number_limit})){
			$j++;
		}
	}
	close OT;
	$cfg{ROH}{scaffold_number_limit} = $j;

	open CL, ">$var{shpath}/cmd_ROH.list";
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 4);
		open SH, ">$var{shpath}/ROH.$pop_name.sh";

		open OT, ">$var{outpath}/ROH/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;

		print SH "cd $var{outpath}/ROH\n";
		print SH 'bcftools annotate --threads 10 --rename-chrs', " $var{outpath}/ROH/chr_map.list" ," $var{vcf}", '|perl -ne \'if (/#\S+ID=(\d+)/){if($1<=',"$scaffold_number_limit",'){print;}}elsif(/^#/){print;}elsif(/^(\d+)\s+/){if($1<= ',"$scaffold_number_limit",'){print;}}\'', ">$var{outpath}/ROH/$pop_name.pre.vcf\n";

		print SH "vcftools --vcf $var{outpath}/ROH/$pop_name.pre.vcf --keep $var{outpath}/ROH/$pop_name.list --plink --out $pop_name.plink \n";
		print SH "plink --file $pop_name.plink --make-bed --chr-set $scaffold_number_limit no-xy no-mt no-y --out $pop_name.plink\n";
		print SH "plink --bfile $var{outpath}/ROH/$pop_name.plink --chr-set $scaffold_number_limit no-xy no-mt no-y --homozyg --homozyg-window-het 2 --homozyg-kb $window_size --homozyg-gap 50 --homozyg-window-snp 100 --homozyg-window-missing 5 --out $pop_name\n";
		close SH;

		print CL "sh $var{shpath}/ROH.$pop_name.sh 1>$var{shpath}/ROH.$pop_name.sh.o 2>$var{shpath}/ROH.$pop_name.sh.e \n";
	}
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_ROH_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/cmd_ROH.list` unless (defined $opts{skipsh});
}

sub LD{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};
	my %pop = %{$var{pop}};

	if (defined $opts{genome}){
		$var{genome} = PopGenome_Shared::LOADREF($opts{genome});
	}else {
		$var{genome} = PopGenome_Shared::LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	}

	if ( !-d "$var{outpath}/LD" ) {make_path "$var{outpath}/LD" or die "Failed to create path: $var{outpath}/LD";}
	my $ld_outpath = "$var{outpath}/LD";

	my $ori_gzvcf;
	if ($var{ploidy} == 1 ){
		$ori_gzvcf = $var{vcf};
		$var{vcf} = "$var{outpath}/LD/diploid.high_confidence.vcf.gz";
	}

	`cp -f $Bin/lib/LD.R $var{shpath}/LD.R`;
	
	if ($var{ploidy} == 1 ){
		open CL, ">$var{shpath}/cmd_LD_s0.list";
		open SH, ">$var{shpath}/LD_s0.sh";
		print SH <<EOF;
zcat $ori_gzvcf|java -jar /root/jvarkit/dist/vcffilterjdk.jar -e 'return new VariantContextBuilder(variant).genotypes(variant.getGenotypes().stream().map(G->!G.isCalled()?GenotypeBuilder.createMissing(G.getSampleName(),2):G).map(G->G.isCalled() && G.getPloidy()==1?new GenotypeBuilder(G).alleles(Arrays.asList(G.getAllele(0),G.getAllele(0))).make():G).collect(Collectors.toList())).attribute("AC",variant.getGenotypes().stream().mapToInt(G->G.isCalled() && G.getPloidy()==1 && !G.getAllele(0).isReference()?2:G.getAlleles().size()).sum()).attribute("AN",variant.getGenotypes().stream().mapToInt(G->G.isCalled() && G.getPloidy()==1 ?2:G.getAlleles().size()).sum()).make();'|bgzip -c >$var{vcf}
EOF
		close SH;
		close CL;
		print CL "sh $var{shpath}/LD_s0.sh 1>$var{shpath}/LD_s0.sh.o 2>$var{shpath}/LD_s0.sh.e \n";
		`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_LD_s0_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=1 -binding linear:1' -m 100 -r $var{shpath}/cmd_LD_s0.list` unless (defined $opts{skipsh});
	}
	
	open CL, ">$var{shpath}/cmd_LD_s1.list";
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 6);
		open OT, ">$ld_outpath/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;

		#write a script to generate a vcf file for each population
		open SH, ">$var{shpath}/LD_s1_$pop_name.sh";
		print SH "cd $ld_outpath\n";
		print SH "vcftools --gzvcf $var{vcf} --keep $ld_outpath/$pop_name.list --recode --stdout |bgzip -c >$ld_outpath/$pop_name.snp.vcf.gz\n";
		print SH "vcftools --gzvcf $pop_name.snp.vcf.gz --singletons --stdout >$ld_outpath/$pop_name.out.singletons\n";
		print SH "vcftools --gzvcf $ld_outpath/$pop_name.snp.vcf.gz --exclude-positions $ld_outpath/$pop_name.out.singletons --recode --stdout |bgzip -c  >$ld_outpath/$pop_name.snp.noSingle.vcf.gz\n";
		close SH;
		print CL "sh $var{shpath}/LD_s1_$pop_name.sh 1>$var{shpath}/LD_s1_$pop_name.sh.o 2>$var{shpath}/LD_s1_$pop_name.sh.e \n";
	}
	close CL;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_LD_s1_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/cmd_LD_s1.list` unless (defined $opts{skipsh});

	print CLSH ">$var{shpath}/cmd_LD_R.list";
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 6);

		my $i=1;
		open PCL, ">$var{shpath}/LD.$pop_name.cmd.list";

#		open RIN, ">$ld_outpath/$pop_name.R.input";
		open GRIN, ">$ld_outpath/$pop_name.GLD.R.input";
		foreach my $id(sort { $var{genome}->{len}{$b} <=> $var{genome}->{len}{$a} } keys %{$var{genome}->{len}}){
			next unless (($var{genome}->{len}{$id}>=$cfg{ld}{scaffold_length_cutoff})&&($i<=$cfg{ld}{scaffold_number_limit}));
			open IDSH, ">$var{shpath}/$pop_name.$id.LD.sh";
			print IDSH "cd $ld_outpath\n";
			print IDSH "vcftools --gzvcf $pop_name.snp.noSingle.vcf.gz --chr $id --recode --stdout |bgzip -c  >$pop_name.$id.snp.noSingle.vcf.gz\n";
			print IDSH "rm -rf $pop_name.$id.snp.noSingle.beagl*\n";
			#print IDSH "beagle gt=$pop_name.$id.snp.noSingle.vcf.gz out=$pop_name.$id.snp.noSingle.beagle\n";
#				print IDSH "vcftools --gzvcf $pop_name.$id.snp.noSingle.beagle.vcf.gz --hap-r2 --ld-window-bp 1000000 --stdout |grep -v nan > $pop_name.$id.LD_window_1M.list\n";
			print IDSH "vcftools --gzvcf $pop_name.$id.snp.noSingle.vcf.gz --geno-r2 --ld-window-bp 1000000 --stdout |grep -v nan | perl -ne '\@a=split /\\t+/;if (/CHR/){print;}elsif(((\$a[2]-\$a[1])>5000)&&(\$i!= 100)){\$i++;}elsif((\$a[2]-\$a[1])<=5000){print;}elsif(((\$a[2]-\$a[1])>5000 )&& (\$i ==100)){print;\$i=0;}'> $pop_name.$id.GLD_window_1M.list\n";
		#	print IDSH "vcftools --gzvcf $pop_name.$id.snp.noSingle.beagle.vcf.gz --geno-r2 --ld-window-bp 1000000 --stdout |grep -v nan | perl -ne '\@a=split /\\t+/;if (/CHR/){print;}elsif(((\$a[2]-\$a[1])>5000)&&(\$i!= 100)){\$i++;}elsif((\$a[2]-\$a[1])<=5000){print;}elsif(((\$a[2]-\$a[1])>5000 )&& (\$i ==100)){print;\$i=0;}'> $pop_name.$id.GLD_window_1M.list\n";
			print IDSH "perl $Bin/lib/window_LD.pl $pop_name.$id.GLD_window_1M.list $cfg{slidingwindow}{windowsize} ",$var{genome}->{len}{$id}," > $pop_name.$id.GLD_window.stats\n";
#				print IDSH "cat $pop_name.$id.LD_window_1M.list",'|sed 1,1d | awk -F " " \'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($3-$2),$5}\''," >$pop_name.$id.LD_window_1M.summary\n";
			print IDSH "cat $pop_name.$id.GLD_window_1M.list",'|sed 1,1d | awk -F " " \'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($3-$2),$5}\''," >$pop_name.$id.GLD_window_1M.summary\n";
			close IDSH;

			print PCL "sh $var{shpath}/$pop_name.$id.LD.sh 1>$var{shpath}/$pop_name.$id.LD.sh.o 2>$var{shpath}/$pop_name.$id.LD.sh.e\n";
	
#				print RIN "$pop_name.$id.LD_window_1M.summary\n";
			print GRIN "$pop_name.$id.GLD_window_1M.summary\n";
			$i++;
		}

		close PCL;
		`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd.LD.$pop_name.qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/LD.$pop_name.cmd.list` unless (defined $opts{skipsh});
#		close RIN;
		close GRIN;

#		my @p;
#		push @p, $ld_outpath; 
#		push @p, "$pop_name.R.input";
#		push @p, "$pop_name.LD.png";
#		my $p = join (' ',@p);
#		print CLSH "Rscript --vanilla $var{shpath}/LD.R $p\n";

		my @gp;
		push @gp, $ld_outpath; 
		push @gp, "$pop_name.GLD.R.input";
		push @gp, "$pop_name.GLD.png";
		my $gp = join (' ',@gp);
		print CLSH "Rscript --vanilla $var{shpath}/LD.R $gp\n";
	}
	close CLSH;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_LD_R.qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/cmd_LD_R.list` unless (defined $opts{skipsh});
}

1;