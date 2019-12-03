package PopGenome_Variant_Filtering;
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
		'gatk_filtering',
		'freebayes_filtering',
		'intersection',
		'advanced_filtering',
		'help',
		'skipsh');

	die "please provide the correct configuration file" unless ((defined $opts{config}) && (-e $opts{config}));

	if (defined $opts{allsteps}){
		$opts{variant_filtering} = 1;
		$opts{intersection} = 1;
	}

	my $yaml = YAML::Tiny->read( $opts{config} );
	my %cfg = %{$yaml->[0]};
	my %samplelist = %{$cfg{fqdata}};
	
	$var{samplelist}=\%samplelist;
	$var{cfg}=\%cfg;

	#set the number of threads
	if (defined $opts{threads}){
		$var{threads} = $opts{threads};
	}elsif(defined $cfg{args}{threads}){
		$var{threads} = $cfg{args}{threads};
	}else{
		$var{threads} = 4;
	}


	$var{outpath} = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$cfg{ref}{choose}/Filtered_variants";
	if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 

	$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$cfg{ref}{choose}/Filtered_variants";
	if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}

	die "please add genome path into configuration file" unless (defined $cfg{ref}{db}{$cfg{ref}{choose}}{path});
	$var{reference} = $cfg{ref}{db}{$cfg{ref}{choose}}{path};
	die "$var{reference} does not exists" unless (-e $var{reference});

	if (defined $opts{variant_filtering}){ &VariantFiltering (\%var,\%opts);}
	if (defined $opts{freebayes_filtering}){ &FreebayesBasicFiltering (\%var,\%opts);}
	if (defined $opts{advanced_filtering}){ &AdvancedFiltering (\%var,\%opts);}
#	if (defined $opts{intersection}){ &ReadReport (\%var,\%opts);}
}

#################################
#			   #
#    Step 1g Variant filtering  #
#			   #
#################################
sub SampleFiltering {
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};


	
}

sub GATKBasicFiltering {
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	open SH, ">$var{shpath}/gatk_basic_filtering.sh";
	print SH "cd $var{outpath}\n";
	print SH "gatk SelectVariants \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "	-V $var{outpath}/../JointCalling/JointCalling.HC.vcf.gz \\\n";
	print SH "	--select-type-to-include SNP \\\n";
	print SH "	-O $var{outpath}/../JointCalling/JointCalling_raw_snps1st.vcf && echo \"** GVCF JointCalling/JointCalling_raw_snps1st done\" && \\\n";

	print SH "gatk VariantFiltration \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "	-V $var{outpath}/../JointCalling/JointCalling_raw_snps1st.vcf \\\n";
	print SH "	--filter-expression \"$cfg{step1}{variant_filtering}{snp}\" \\\n";
	print SH "	--filter-name \"my_snp_filter\" \\\n";
	print SH "	-O $var{outpath}/../JointCalling/JointCalling_filtered_snps1st.vcf && echo \"** GVCF JointCalling/JointCalling_raw_snps1st done\" \n";

	print SH "gatk SelectVariants \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "	-V $var{outpath}/../JointCalling/JointCalling.HC.vcf.gz \\\n";
	print SH "	--select-type-to-include INDEL \\\n";
	print SH "	--maxIndelSize 60 \\\n";
	print SH "	-O $var{outpath}/../JointCalling/JointCalling_raw_indels1st.vcf && echo \"** GVCF JointCalling/JointCalling_raw_snps1st done\" && \\\n";

	print SH "gatk VariantFiltration \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "	-V $var{outpath}/../JointCalling/JointCalling_raw_indels1st.vcf \\\n";
	print SH "	--filter-expression \"$cfg{step1}{variant_filtering}{indel}\" \\\n";
	print SH "	--filter-name \"my_indel_filter\" \\\n";
	print SH "	-O $var{outpath}/../JointCalling/JointCalling_filtered_indels1st.vcf && echo \"** JointCalling/JointCalling_raw_snps1st done\" \n";
	
	close SH;

	#switch on the bash running
	open CL, ">$var{shpath}/cmd_gatk_basic_filtering.list";
	print CL "sh $var{shpath}/gatk_basic_filtering.sh 1>$var{shpath}/gatk_basic_filtering.sh.o 2>$var{shpath}/gatk_basic_filtering.sh.e\n";
	close CL;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/gatk_basic_filtering_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/cmd_gatk_basic_filtering.list` unless (defined $opts{skipsh});
}

sub FreebayesBasicFiltering {
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	open SH, ">$var{shpath}/freebayes_basic_filtering.sh";

	print SH "#!/bin/sh\ncd $var{outpath}/\n";
	print SH "vcftools --gzvcf $var{outpath}/../FreebayesCalling/freebayes_joint_calling.vcf.gz --missing-indv\n";
	print SH "awk \'\$5 > 0.2\' out.imiss | cut -f1 > lowDP.indv\n";
	print SH "vcftools --gzvcf $var{outpath}/../FreebayesCalling/freebayes_joint_calling.vcf.gz --max-missing 0.8 --max-alleles 2 --minQ 30 --remove-filtered-all --remove lowDP.indv --recode --recode-INFO-all --stdout | bgzip -c > freebayes_joint_calling_filtered1.vcf.gz";

	close SH;

	open CL, ">$var{shpath}/cmd_freebayes_basic_filtering.list";
	print CL "sh $var{shpath}/freebayes_basic_filtering.sh 1>$var{shpath}/freebayes_basic_filtering.sh.o 2>$var{shpath}/freebayes_basic_filtering.sh.e\n";
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/freebayes_basic_filtering_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/cmd_freebayes_basic_filtering.list` unless (defined $opts{skipsh});
}

sub AdvancedFiltering {
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	#print SH "vcftools --gzvcf freebayes_joint_calling_filtered1.vcf.gz --max-alleles 2 --remove-filtered-all --recode --recode-INFO-all --stdout \| vcfsnps\|bgzip -c > $var{outpath}/basic_snp.vcf.gz \n";
	`cp -f $var{outpath}/freebayes_joint_calling_filtered1.vcf.gz $var{outpath}/basic_snp.vcf.gz`;
	my $genome=PopGenome_Shared::LOADREF($var{reference});

	###filter SNP by depth if needed
	my $sample_size = keys %{$cfg{fqdata}};
	open (IN, "zcat $var{outpath}/basic_snp.vcf.gz|") or die $!;
	my @dp;
	while (<IN>){
    	next if (/#/);
    	if (/DP=([\d\.]+)/){push @dp, $1;}
	}
	my $p50 = int(@dp * .5);
	my $avgdp50 = 3*$dp[$p50]/$sample_size;
	close IN; 
	
	my $chrcmd="";
	###select SNPs by chromosomes if needed
	if (defined $cfg{step1}{variant_filtering}{chr}){
		my @a = split /,/, $cfg{step1}{variant_filtering}{chr};
		foreach (@a){
			$chrcmd .= " --chr $_";
		}
	}

	open CL, ">$var{shpath}/cmd_variant_filtering_s1.list";
	open SH, ">$var{shpath}/variant_filtering_s1.sh";
	#default SNV sites 1: PASS + max-meanDP: 3*$avgdp50 + max-missing: 0.8 + biallelic + selected chr (optional)
	print SH "vcftools --gzvcf $var{outpath}/basic_snp.vcf.gz $chrcmd --min-meanDP 5 --max-meanDP $avgdp50 --max-missing 0.8 --max-alleles 2 --remove-filtered-all --recode --recode-INFO-all --stdout \| bgzip -c > $var{outpath}/basic_snp_dp.vcf.gz \n";
	#switch on the bash running
	print CL "sh $var{shpath}/variant_filtering_s1.sh 1>$var{shpath}/variant_filtering_s1.sh.o 2>$var{shpath}/variant_filtering_s1.sh.e\n";
	close CL;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/variant_filtering_s1_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/cmd_variant_filtering_s1.list` unless (defined $opts{skipsh});

	$cfg{variant_filtering}{basic_vcf}="$var{outpath}/basic_snp_dp.vcf.gz";
	### prepare the setting for low LD prunning
	$cfg{variant_filtering}{scaffold_number_limit} = 95 unless (defined $cfg{variant_filtering}{scaffold_number_limit});
	$cfg{variant_filtering}{scaffold_length_cutoff} = 0 unless (defined $cfg{variant_filtering}{scaffold_length_cutoff});
	#set ld prunning cut off
	$cfg{variant_filtering}{ldwindowsize} = 10 unless (defined $cfg{variant_filtering}{ldwindowsize});
	$cfg{variant_filtering}{ldwindowstep} = 5 unless (defined $cfg{variant_filtering}{ldwindowstep});
	$cfg{variant_filtering}{ldcutoff} = 0.3 unless (defined $cfg{variant_filtering}{ldcutoff});

	### load reference and map scaffold to number list 
	my $i=1;my $j=0;
	open OT, ">$var{outpath}/chr_map.list";
	foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
		print OT "$id\t$i\n";
		$i++;
		if (($genome->{len}{$id} >= $cfg{variant_filtering}{scaffold_length_cutoff})&&($j < $cfg{variant_filtering}{scaffold_number_limit})){
			$j++;
		}
	}
	close OT;

	$cfg{variant_filtering}{scaffold_number_limit} = $j;
	my $scaffold_number_limit = $cfg{variant_filtering}{scaffold_number_limit};
	print "scaffold_number_limit: $scaffold_number_limit\n";
	
	#Based on default SNV sites + no singletons + no missing data + minDP 2 + minQ 30  (high quality SNPs)
	open CL, ">$var{shpath}/cmd_variant_filtering_s2.list";
	open SH, ">$var{shpath}/variant_filtering_s2.sh";
	print SH "#!/bin/sh\ncd $var{outpath}\n";
	print SH "vcftools --gzvcf $cfg{variant_filtering}{basic_vcf} --singletons --stdout >$var{outpath}/singletons.list\n";
	print SH "vcftools --gzvcf $cfg{variant_filtering}{basic_vcf} --exclude-positions $var{outpath}/singletons.list --max-missing 1 --max-alleles 2 --recode --recode-INFO-all --stdout |bgzip -c >$var{outpath}/nosingle_snp_dp.vcf.gz\n";
	print SH "tabix -f $var{outpath}/nosingle_snp_dp.vcf.gz\n";
	
	#Based on high quality SNV sites + low LD
	print SH 'bcftools annotate --threads $var{threads} --rename-chrs', " $var{outpath}/chr_map.list" ," $var{outpath}/nosingle_snp_dp.vcf.gz", '|perl -ne \'if (/#\S+ID=(\d+)/){if($1<=',"$scaffold_number_limit",'){print;}}elsif(/^#/){print;}elsif(/^(\d+)\s+/){if($1<= ',"$scaffold_number_limit",'){print;}}\'|vcftools --vcf - --plink',"\n";
	
	print SH "plink --file out --make-bed --chr-set $scaffold_number_limit no-xy no-mt no-y\n";
	print SH "plink --bfile plink --indep-pairwise $cfg{variant_filtering}{ldwindowsize} $cfg{variant_filtering}{ldwindowstep} $cfg{variant_filtering}{ldcutoff} --chr-set $scaffold_number_limit no-xy no-mt no-y\n";
	print SH "plink --bfile plink --extract plink.prune.in --make-bed --out nosingle_snp_dp_prunned --chr-set $scaffold_number_limit no-xy no-mt no-y\n";
	print SH "plink --bfile nosingle_snp_dp_prunned --chr-set $scaffold_number_limit no-xy no-mt no-y --recode vcf --out nosingle_snp_dp_pre\n";
	print SH "cat $var{outpath}/nosingle_snp_dp_pre.vcf |perl -ne \'print unless (/CHROM/);if(/CHROM/){s/_\\S+//g;print;}\' | bgzip -c >$var{outpath}/nosingle_snp_dp_prunned.vcf.gz\n";
	print SH "ls $var{outpath}/nosingle_snp_dp_prunned.vcf.gz && echo \"** finish variant_filtering_s2 **\" > $var{shpath}/variant_filtering_s2.finished.txt";
	close SH;
	print CL "sh $var{shpath}/variant_filtering_s2.sh 1>$var{shpath}/variant_filtering_s2.sh.o 2>$var{shpath}/variant_filtering_s2.sh.e\n";
	close CL;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/variant_filtering_s2_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/cmd_variant_filtering_s2.list` unless (defined $opts{skipsh});

	$cfg{variant_filtering}{plink_data}="$var{outpath}/nosingle_snp_dp_prunned";
	$cfg{variant_filtering}{high_confidence_vcf}="$var{outpath}/nosingle_snp_dp.vcf.gz";
	$cfg{variant_filtering}{lowld_high_confidence_vcf}="$var{outpath}/nosingle_snp_dp_prunned.vcf.gz";
	
	### format a report
	my %report;
	if ($cfg{variant_filtering}{basic_vcf} =~ /.gz$/) { 
		open(IN, "gunzip -c $cfg{variant_filtering}{basic_vcf} |") || die "can’t open pipe to $cfg{variant_filtering}{basic_vcf}";
	} else { 
		open(IN, $cfg{variant_filtering}{basic_vcf}) || die "can’t open $cfg{variant_filtering}{basic_vcf}";
	}

	$report{snv1}{number}=0;
	while (<IN>){
		next if (/^#/);
		$report{snv1}{number} ++;
	}
	close IN;

	$report{snv1}{singletons} = `wc -l $var{outpath}/singletons.list | cut -d" " -f1`;
	$report{snv1}{singletons} = $report{snv1}{singletons} - 1;

	if ($cfg{step1}{variant_filtering}{high_confidence_vcf} =~ /.gz$/) {
		open(IN, "gunzip -c $cfg{step1}{variant_filtering}{high_confidence_vcf} |") || die "can’t open pipe to $cfg{step1}{high_confidence_vcf}{vcf}";
	}
	else {
		open(IN, $cfg{step1}{variant_filtering}{high_confidence_vcf}) || die "can’t open $cfg{step1}{variant_filtering}{high_confidence_vcf}";
	}
	$report{snv2}{number}=0;
	while (<IN>){
		next if (/^#/);
		$report{snv2}{number} ++;
	}
	close IN;

	if ($cfg{variant_filtering}{lowld_high_confidence_vcf} =~ /.gz$/) {
		open(IN, "gunzip -c $cfg{step1}{variant_filtering}{lowld_high_confidence_vcf} |") || die "can’t open pipe to $cfg{step1}{variant_filtering}{lowld_high_confidence_vcf}";
	}
	else {
		open(IN, $cfg{variant_filtering}{lowld_high_confidence_vcf}) || die "can’t open $cfg{step1}{variant_filtering}{lowld_high_confidence_vcf}";
	}

	$report{snv3}{number}=0;
	while (<IN>){
		next if (/^#/);
		$report{snv3}{number} ++;
	}
	close IN;

	####
	open OT, ">$var{outpath}/snv.summary.txt";
	print OT "SNV set\tSNV size\tSingleton size\n";
	print OT "PASS.SNP.DP.vcf.gz\t$report{snv1}{number}\t$report{snv1}{singletons}\n";
	print OT "high_confidence.vcf.gz\t$report{snv2}{number}\t0\n";
	print OT "high_confidence_prunned.vcf.gz\t$report{snv3}{number}\t0\n";
	close OT; 

}

1;