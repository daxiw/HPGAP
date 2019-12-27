package PopGenome_Slidingwindow;
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
		'threads',
		'vcf=s',
		'nonsyn',
		'genome=s',
		'slidingwindow',
		'help',
		'skipsh');

	die "please provide the correct configuration file" unless ((defined $opts{config}) && (-e $opts{config}));

	if (defined $opts{allsteps}){
		$opts{slidingwindow} = 1;
	}

	%var = %{PopGenome_Shared::CombineCfg("$Bin/lib/parameter.yml",\%opts,"Slidingwindow")};

	if (defined $opts{slidingwindow}){ &SLIDINGWINDOW (\%var,\%opts);}

}

#################################
#			   #
#   	step4_SlidingWindow 	#
#			   #
#################################
sub SLIDINGWINDOW{

	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %pop = %{$var{pop}};
	
	if (defined $opts{genome}){
		$var{genome} = PopGenome_Shared::LOADREF($opts{genome});
	}else {
		$var{genome} = PopGenome_Shared::LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	}
	#my $path = /root/diploSHIC
	my $path = "/ldfssz1/ST_INFECTION/P17Z10200N0536_Echinococcus/USER/wangdaxi/Github/diploSHIC";

	my $window_size = $cfg{slidingwindow}{windowsize};
	my $scaffold_number_limit = $cfg{slidingwindow}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{slidingwindow}{scaffold_length_cutoff};

	`cp -f $Bin/lib/Diversity.R $var{shpath}/Diversity.R`;
	
	my $gzvcf = $var{vcf};

	#generate a scaffold size list
	open SL, ">$var{outpath}/chr.size.list";
	foreach my $id(sort { $var{genome}->{len}{$b} <=> $var{genome}->{len}{$a} } keys %{$var{genome}->{len}}){
		print SL "$id\t$var{genome}->{len}{$id}\n";
	}close SL;

	#`grep CDS $cfg{slidingwindow}{gff} >$var{outpath}/CDS.gff`;
	#generate Bed files
	my $i = 1;
	foreach my $id(sort { $var{genome}->{len}{$b} <=> $var{genome}->{len}{$a} } keys %{$var{genome}->{len}}){
		if (($var{genome}->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){
			open BED, ">$var{outpath}/$id.bed";
			for (my $j = 0;($j + $window_size) <= $var{genome}->{len}{$id};$j+=$window_size ){
				my $start = $j; my $end = $j + $window_size; 
				print BED "$id\t$start\t$end\n";
			}
			close BED;
		#	`bedtools intersect -a $var{outpath}/$id.bed -b $var{outpath}/CDS.gff -wo >$var{outpath}/$id.intersected.bed`;
		}
	}

	open CL1, ">$var{shpath}/Slidingwindow.cmd1.list";
	open CL2, ">$var{shpath}/Slidingwindow.cmd2.list";
	open CL3, ">$var{shpath}/Slidingwindow.cmd3.list";
	open CL4, ">$var{shpath}/Slidingwindow.cmd4.list";
	open CL5, ">$var{shpath}/Slidingwindow.cmd5.list";

	#loop for each available population
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 6);
		open OT, ">$var{outpath}/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;
		# generate a vcf file for each population
		print CL1 "vcftools --gzvcf $gzvcf --keep $var{outpath}/$pop_name.list --recode --stdout |bgzip -c >$var{outpath}/$pop_name.SNP.vcf.gz\n";
		
		
		# calculate syn and nonsyn diversity values in each sliding window
		if (defined $opts{nonsyn}){
			print CL1 "vcftools --gzvcf $var{outpath}/snpEff.nonsyn.vcf.gz --keep $var{outpath}/$pop_name.list --recode --stdout |bgzip -c >$var{outpath}/$pop_name.snpEff.nonsyn.vcf.gz\n";
			print CL1 "vcftools --gzvcf $var{outpath}/snpEff.syn.vcf.gz --keep $var{outpath}/$pop_name.list --recode --stdout |bgzip -c >$var{outpath}/$pop_name.snpEff.syn.vcf.gz\n";
			print CL2 "vcftools --gzvcf $var{outpath}/$pop_name.snpEff.nonsyn.vcf.gz --window-pi $window_size --stdout >$var{outpath}/$pop_name.snpEff.nonsyn.pi.list\n";
			print CL2 "vcftools --gzvcf $var{outpath}/$pop_name.snpEff.syn.vcf.gz --window-pi $window_size --stdout >$var{outpath}/$pop_name.snpEff.syn.pi.list\n";
			print CL3 "PiSynNonSyn.pl $var{outpath}/$pop_name.snpEff.nonsyn.pi.list $var{outpath}/$pop_name.snpEff.nonsyn.pi.list SynSitesBin.list >$var{outpath}/$pop_name.PiSynNonSyn.list\n";
		}

		# plot diversity
		my $i = 1;
		foreach my $id(sort { $var{genome}->{len}{$b} <=> $var{genome}->{len}{$a} } keys %{$var{genome}->{len}}){
			if (($var{genome}->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){
				
				my $len=$var{genome}->{len}{$id};
				my $bigwin=11*$window_size;
				# get the diploid slidingwindow stats
				print CL4 "python $path/diploSHIC.py fvecVcf diploid $var{outpath}/$pop_name.SNP.vcf.gz $id $len $var{outpath}/$id.$pop_name.SNP.fvec --winSize $bigwin --statFileName $var{outpath}/$id.$pop_name.stat.result 1>$var{shpath}/$id.$pop_name.stat.result.o 2>$var{shpath}/$id.$pop_name.stat.result.e\n" if ($var{ploidy} == 2);
				# get the haploid slidingwindow stats
				print CL4 "python $path/diploSHIC.py fvecVcf haploid $var{outpath}/$pop_name.SNP.vcf.gz $id $len $var{outpath}/$id.$pop_name.SNP.fvec --winSize $bigwin --ancFileName $cfg{ref}{db}{$cfg{ref}{choose}}{path} --statFileName $var{outpath}/$id.$pop_name.stat.result 1>$var{shpath}/$id.$pop_name.stat.result.o 2>$var{shpath}/$id.$pop_name.stat.result.e\n" if ($var{ploidy} == 1);

				# calculate GC content, generate the table
				open GC, ">$var{outpath}/$id.$pop_name.GCstat.result";
				for (my $j = 0;($j + $window_size) <= $var{genome}->{len}{$id};$j+=$window_size ){
					my $GCcount = uc(substr($var{genome}->{seq}{$id},$j,$window_size)) =~ tr/GC//;
					my $start = $j + 1; my $end = $j + $window_size; my $GCcontent = $GCcount/$window_size;
					print GC "$id\t$start\t$end\t$GCcontent\n";
				}
				close GC;

				my @p;
				push @p, $var{outpath}; 
				push @p, "$id.$pop_name.Pi.list";
				push @p, "$id.$pop_name.SNPdensity.list";
				push @p, "$id.$pop_name.TajimaD.list";
				push @p, "$id.$pop_name.diversity.png";
				push @p, "$id.$pop_name.TajimaD.png";
				push @p, $pop{$pop_name}{count};
				push @p, $window_size;
				push @p, $id;
				my $p = join (' ',@p);
				
				open IDSH, ">$var{shpath}/$id.$pop_name.Slidingwindow.sh";
				print IDSH "grep -w \"$id\" $var{outpath}/$pop_name.Pi.list >$var{outpath}/$id.$pop_name.Pi.list\n";
				print IDSH "grep -w \"$id\" $var{outpath}/$pop_name.SNPdensity.list >$var{outpath}/$id.$pop_name.SNPdensity.list\n";
				print IDSH "grep -w \"$id\" $var{outpath}/$pop_name.TajimaD.list >$var{outpath}/$id.$pop_name.TajimaD.list\n";
				print IDSH "Rscript --vanilla $var{shpath}/Diversity.R $p\n";
				close IDSH;
				print CL5 "sh $var{shpath}/$id.$pop_name.Slidingwindow.sh 1>$var{shpath}/$id.$pop_name.Slidingwindow.sh.o 2>$var{shpath}/$id.$pop_name.Slidingwindow.sh.e\n";

				$i++;
			}
		}
	}
	close CL1;
	close CL2;
	close CL3;
	close CL4;
	close CL5;
	
	###
	if (defined $opts{nonsyn}){
		open CL0, ">$var{shpath}/Slidingwindow.cmd0.list";
		open SH, ">$var{shpath}/Slidingwindow.sh";
		print SH "cd $var{outpath}\n";
		### annotate  vcf file
		print SH "snpEff $cfg{slidingwindow}{snpeff_species} $gzvcf |bgzip -c >$var{outpath}/snpEff.vcf.gz\n";
		### extract nonsyn snvs
		print SH "zcat $var{outpath}/snpEff.vcf.gz |SnpSift filter \"(ANN[*].BIOTYPE has 'protein_coding') & ((ANN[*].EFFECT has 'missense_variant') | (ANN[*].EFFECT = 'start_lost')| (ANN[*].EFFECT = 'stop_gained')| (ANN[*].EFFECT = 'stop_lost'))\" |perl -ne 'if(/#/){print;}elsif(/ANN/){s/ANN\\S+/./g;print;}else{print;}'|bgzip -c >$var{outpath}/snpEff.nonsyn.vcf.gz\n";
		### extract syn snvs
		print SH "zcat $var{outpath}/snpEff.vcf.gz |SnpSift filter \"(ANN[*].BIOTYPE has 'protein_coding') & ((ANN[*].EFFECT = 'synonymous_variant') | (ANN[*].EFFECT = 'stop_retained_variant') | (ANN[*].EFFECT = 'start_retained'))\" |perl -ne 'if(/#/){print;}elsif(/ANN/){s/ANN\\S+/./g;print;}else{print;}'|bgzip -c >$var{outpath}/snpEff.syn.vcf.gz\n";
		print SH "Number_syn_nonsyn_sites.pl $cfg{ref}{db}{$cfg{ref}{choose}}{path} $cfg{slidingwindow}{gff} test >$var{outpath}/nonsyn.sites.list\n";
		print SH "SynSitesBin.pl chr.size.list $var{outpath}/nonsyn.sites.list |grep -v NULL|sort -k1,1 -k2,2n >$var{outpath}/SynSitesBin.list\n";

		print CL0 "sh $var{shpath}/Slidingwindow.sh 1>$var{shpath}/Slidingwindow.sh.o 2>$var{shpath}/Slidingwindow.sh.e\n";
		close CLO;
		`perl $Bin/lib/qsub.pl -d $var{shpath}/Slidingwindow.cmd0.list.qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/Slidingwindow.cmd0.list` unless (defined $opts{skipsh});
	}

	`perl $Bin/lib/qsub.pl -d $var{shpath}/Slidingwindow.cmd1.list.qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/Slidingwindow.cmd1.list` unless (defined $opts{skipsh});
	
	if (defined $opts{nonsyn}){
		`perl $Bin/lib/qsub.pl -d $var{shpath}/Slidingwindow.cmd2.list.qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/Slidingwindow.cmd2.list` unless (defined $opts{skipsh});
		`perl $Bin/lib/qsub.pl -d $var{shpath}/Slidingwindow.cmd3.list.qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/Slidingwindow.cmd3.list` unless (defined $opts{skipsh});
	}

	`perl $Bin/lib/qsub.pl -d $var{shpath}/Slidingwindow.cmd4.list.qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/Slidingwindow.cmd4.list` unless (defined $opts{skipsh});
	#`perl $Bin/lib/qsub.pl -d $var{shpath}/Slidingwindow.cmd5.list.qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/Slidingwindow.cmd5.list` unless (defined $opts{skipsh});

	close SH;
=head

#This is a comment.  Note the blank lines.

	foreach my $pop_name (keys %pop){
		open ALLSTAT, ">$var{outpath}/$pop_name.allstat.txt";
		next unless ($pop{$pop_name}{count} > 6);
		my $i = 1;
		foreach my $id(sort { $var{genome}->{len}{$b} <=> $var{genome}->{len}{$a} } keys %{$var{genome}->{len}}){
			if (($var{genome}->{len}{$id}>=$scaffold_length_cutoff)&&($i<=$scaffold_number_limit)){

				open IDALLSTAT, ">$var{outpath}/$id.$pop_name.allstat.txt";
				my %stats;
				open STAT, "$var{outpath}/$id.$pop_name.stat.result";
				while(<STAT>){
					next if (/chrom/);
					chomp;
					my @a = split /\s+/;
					$stats{$a[1]}{diversity}="$a[3]\t$a[4]\t$a[5]";
					#print "$a[0]\t$a[1]\t$a[2]\t$stats{$a[1]}{diversity}\n";
				}
				close STAT;

				open STAT, "$var{outpath}/$id.intersected.bed";
				while(<STAT>){
					next if (/chrom/);
					chomp;
					my @a = split /\s+/;
					my $start = $a[1]+1;
					my $p = $a[-1]/$window_size;
					if (exists $stats{$start}{p_coding}){
						$stats{$start}{p_coding}+= $p; 
					}else{
						$stats{$start}{p_coding}=$p;
					}
					#print "$a[0]\t$a[1]\t$a[2]\t$stats{$a[1]}{diversity}\n";
				}
				close STAT;

				open STAT, "$var{outpath}/../LD/$pop_name.$id.GLD_window.stats";
				while(<STAT>){
					chomp;
					my @a = split /\t+/;
					$stats{$a[1]}{r2}="$a[3]";
				}
				close STAT;

				open STAT, "$var{outpath}/$id.$pop_name.GCstat.result";
				while(<STAT>){
					chomp;
					my @a = split /\s+/;
					$stats{$a[1]}{r2} = "nan" unless (exists $stats{$a[1]}{r2});
					$stats{$a[1]}{p_coding} = 0 unless (exists $stats{$a[1]}{p_coding});
					if (exists $stats{$a[1]}{diversity}){
						print IDALLSTAT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$stats{$a[1]}{p_coding}\t$stats{$a[1]}{diversity}\t$stats{$a[1]}{r2}\n";
						print ALLSTAT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$stats{$a[1]}{p_coding}\t$stats{$a[1]}{diversity}\t$stats{$a[1]}{r2}\n";
					}
				}
				close STAT;
				close IDALLSTAT;
			}
			$i++;
		}
		close ALLSTAT;

		my %stats;
		open STAT, "$var{outpath}/$pop_name.PiSynNonSyn.list";
		while(<STAT>){
			next if (/CHROM/);
			chomp;
			my @a = split /\t+/;
			$stats{$a[0]}{$a[1]}{pnps}="$a[3]\t$a[4]\t$a[5]";
		}
		close STAT;

		open ALLSTAT, ">$var{outpath}/final.$pop_name.allstat.txt";
		print ALLSTAT "chrom\tstart\tend\tgc_content\tcoding_percentage\tpi\ttheta\ttajimsD\tr2\tps\tpn\tpnps\n";
		open STAT, "$var{outpath}/$pop_name.allstat.txt";
		while(<STAT>){
			chomp;
			my @a = split /\s+/;
			$stats{$a[0]}{$a[1]}{pnps} = "nan\tnan\tnan" unless (exists $stats{$a[0]}{$a[1]}{pnps});
			print ALLSTAT "$_\t$stats{$a[0]}{$a[1]}{pnps}\n";
		}
		close STAT;
		close ALLSTAT;
		`Rscript --vanilla $Bin/lib/Correlation.R $var{outpath} $var{outpath}/final.$pop_name.allstat.txt`;
	}
=cut
}

1;