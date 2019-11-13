package PopGenome_Variant_Calling;
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
		'threads=i',
		'individual_variant_calling',
		'joint_calling',
		'freebayes_calling',
		'recalibration_mode',
		'help',
		'skipsh');

	die "please provide the correct configuration file" unless ((defined $opts{config}) && (-e $opts{config}));

	if (defined $opts{allsteps}){
		$opts{individual_variant_calling} = 1;
		$opts{joint_calling} = 1;
	}

	my $yaml = YAML::Tiny->read( $opts{config} );
	my %cfg = %{$yaml->[0]};
	my %samplelist = %{$cfg{fqdata}};
	
	$var{samplelist}=\%samplelist;
	$var{cfg}=\%cfg;

	#set ploidy
	if (defined $cfg{args}{ploidy}){
		$var{ploidy} = $cfg{args}{ploidy};
	}else{
		$var{ploidy} = 2;
	}

	#set the number of threads
	if (defined $opts{threads}){
		$var{threads} = $opts{threads};
	}elsif(defined $cfg{args}{threads}){
		$var{threads} = $cfg{args}{threads};
	}else{
		$var{threads} = 4;
	}

	#set running mode
	if (defined $opts{recalibration_mode}){
		$var{variant_calling_mode} = "slow";
	} else {
		$var{variant_calling_mode} = "fast";
	}

	die "please add reference genome path into configuration file" unless (defined $cfg{ref}{db}{$cfg{ref}{choose}}{path});
	$var{reference} = $cfg{ref}{db}{$cfg{ref}{choose}}{path};
	die "$var{reference} does not exists" unless (-e $var{reference});

	$var{outpath} = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$cfg{ref}{choose}"; 
	if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 
	$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$cfg{ref}{choose}";
	if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}

	if (defined $opts{freebayes_calling}){ &FreebayesCalling (\%var,\%opts);}

	if (defined $opts{individual_variant_calling}){ &IndividualVariantCalling (\%var,\%opts);}

	if (defined $opts{joint_calling}){ &JointCalling (\%var,\%opts);}

	#### estimate phylogeny of mt genomes ###

}

sub IndividualVariantCalling {
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	open CL, ">$var{shpath}/cmd_variant_calling.list";

	foreach my $sample (keys %samplelist){
		my $sample_outpath="$var{outpath}/$sample"; 
		if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}

		### skip finished samples
		if ((-e "$sample_outpath/$sample.HC.gvcf.gz") && (!defined $opts{overwrite})){
			$samplelist{$sample}{finish_flag} = "finished";
			next;
		}else{
			`rm -f $var{shpath}/$sample.variant_calling.finished.txt`;
		}

		open SH, ">$var{shpath}/$sample.variant_calling.sh";

		print SH "#!/bin/sh\ncd $sample_outpath\n";


		if ($var{variant_calling_mode} eq 'fast'){

			print SH "gatk HaplotypeCaller \\\n";
			print SH "	--emit-ref-confidence GVCF \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-ploidy $var{ploidy} \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
			print SH "	-O $sample.HC.gvcf.gz && echo \"** GVCF ${sample}.HC.g.vcf.gz done\" && \\\n";
		  	print SH "rm -f $sample.sorted.markdup.BQSR2nd.bam && echo \"** variant calling done **\" > $var{shpath}/$sample.variant_calling.finished.txt\n";
		}
		else{ #recalibration mode
			# HaplotypeCaller
		  	print SH "gatk HaplotypeCaller \\\n";
		  	print SH "	-R $var{reference} \\\n";
		  	print SH "	-ploidy $var{ploidy} \\\n";
		 	print SH "	-I $sample.sorted.markdup.bam \\\n";
		  	print SH "	-O $sample.HC.g1st.vcf.gz && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";
		  	###########SNP extraction and filtering#######
		  	# SelectVariants
		  	print SH "gatk SelectVariants \\\n";
		  	print SH "	-R $var{reference} \\\n";
		  	print SH "	-V $sample.HC.g1st.vcf.gz \\\n";
		  	print SH "	--select-type-to-include SNP \\\n";
			print SH "	-O $sample\_raw_snps1st.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";
			# VariantFiltration
			print SH "gatk VariantFiltration \\\n";
		  	print SH "	-R $var{reference} \\\n";
		  	print SH "	-V $sample\_raw_snps1st.vcf \\\n";
		  	print SH "	--filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" \\\n";
			print SH "	--filter-name \"my_snp_filter\" \\\n";
			print SH "	-O $sample\_filtered_snps1st.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

			###########INDEL extraction and filtering#######
			print SH "gatk SelectVariants \\\n";
		  	print SH "	-R $var{reference} \\\n";
		  	print SH "	-V $sample.HC.g1st.vcf.gz \\\n";
		  	print SH "	--select-type-to-include INDEL \\\n";
			print SH "	-O $sample\_raw_indels1st.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

			print SH "gatk VariantFiltration \\\n";
		  	print SH "	-R $var{reference} \\\n";
		  	print SH "	-V $sample\_raw_indels1st.vcf \\\n";
		  	print SH "	--filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" \\\n";
			print SH "	--filter-name \"my_indel_filter\" \\\n";
			print SH "	-O $sample\_filtered_indels1st.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

			print SH "bgzip -f $sample\_filtered_snps1st.vcf\n";
			print SH "tabix -f $sample\_filtered_snps1st.vcf.gz \n";
			print SH "bgzip -f $sample\_filtered_indels1st.vcf\n";
			print SH "tabix -f $sample\_filtered_indels1st.vcf.gz \n";

			print SH "gatk BaseRecalibrator \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
		  	print SH "	-O $sample.sorted.markdup.recal_data.table \\\n";
			print SH "	--known-sites $sample\_filtered_snps1st.vcf.gz \\\n";
			print SH "	--known-sites $sample\_filtered_indels1st.vcf.gz && echo \"** $sample.sorted.markdup.recal_data.table done **\" \n";

			print SH "gatk ApplyBQSR \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
			print SH "	-O $sample.sorted.markdup.BQSR.bam \\\n";
			print SH "	--bqsr-recal-file $sample.sorted.markdup.recal_data.table && \\\n";
			print SH "samtools index $sample.sorted.markdup.BQSR.bam && echo \"** $sample.sorted.markdup.BQSR.bam index done **\" \n";

			print SH "gatk BaseRecalibrator \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.BQSR.bam  \\\n";
		  	print SH "	-O $sample.sorted.markdup.recal_data1st_after.table \\\n";
			print SH "	--known-sites $sample\_filtered_snps1st.vcf.gz \\\n";
			print SH "	--known-sites $sample\_filtered_indels1st.vcf.gz && echo \"** $sample.sorted.markdup.recal_data.table done **\" \n";

		  	# HaplotypeCaller
		  	print SH "gatk HaplotypeCaller \\\n";
		  	print SH "	-R $var{reference} \\\n";
		  	print SH "	-ploidy $var{ploidy} \\\n";
		  	print SH "	-I $sample.sorted.markdup.BQSR.bam \\\n";
		  	print SH "	-O $sample.HC.g2nd.vcf.gz && echo \"** GVCF $sample.HC.g.vcf.gz done\" && \\\n";
		  	print SH "rm -f $sample.sorted.markdup.BQSR.bam\n";
		  	
		  	###########SNP extraction and filtering#######
		  	# SelectVariants
		  	print SH "gatk SelectVariants \\\n";
		  	print SH "	-R $var{reference} \\\n";
		  	print SH "	-V $sample.HC.g2nd.vcf.gz \\\n";
		  	print SH "	--select-type-to-include SNP \\\n";
			print SH "	-O $sample\_raw_snps2nd.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";
			# VariantFiltration
			print SH "gatk VariantFiltration \\\n";
		  	print SH "	-R $var{reference} \\\n";
		  	print SH "	-V $sample\_raw_snps2nd.vcf \\\n";
		  	print SH "	--filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" \\\n";
			print SH "	--filter-name \"my_snp_filter\" \\\n";
			print SH "	-O $sample\_filtered_snps2nd.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

			###########INDEL extraction and filtering#######
			print SH "gatk SelectVariants \\\n";
		  	print SH "	-R $var{reference} \\\n";
		  	print SH "	-V $sample.HC.g2nd.vcf.gz \\\n";
		  	print SH "	--select-type-to-include INDEL \\\n";
			print SH "	-O $sample\_raw_indels2nd.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

			print SH "gatk VariantFiltration \\\n";
		  	print SH "	-R $var{reference} \\\n";
		  	print SH "	-V $sample\_raw_indels2nd.vcf \\\n";
		  	print SH "	--filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" \\\n";
			print SH "	--filter-name \"my_indel_filter\" \\\n";
			print SH "	-O $sample\_filtered_indels2nd.vcf && echo \"** GVCF $sample.HC.g.vcf.gz done\" \n";

			print SH "bgzip -f $sample\_filtered_snps2nd.vcf\n";
			print SH "tabix -f $sample\_filtered_snps2nd.vcf.gz \n";
			print SH "bgzip -f $sample\_filtered_indels2nd.vcf\n";
			print SH "tabix -f $sample\_filtered_indels2nd.vcf.gz \n";

			print SH "gatk BaseRecalibrator \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
		  	print SH "	-O $sample.sorted.markdup.recal_data2nd.table \\\n";
			print SH "	--known-sites $sample\_filtered_snps2nd.vcf.gz \\\n";
			print SH "	--known-sites $sample\_filtered_indels2nd.vcf.gz && echo \"** $sample.sorted.markdup.recal_data.table done **\" \n";

			print SH "gatk ApplyBQSR \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
			print SH "	-O $sample.sorted.markdup.BQSR2nd.bam \\\n";
			print SH "	--bqsr-recal-file $sample.sorted.markdup.recal_data2nd.table && \\\n";
			print SH "samtools index $sample.sorted.markdup.BQSR2nd.bam && echo \"** $sample.sorted.markdup.BQSR2nd.bam index done **\" \n";

			print SH "gatk BaseRecalibrator \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-I $sample.sorted.markdup.BQSR2nd.bam \\\n";
		  	print SH "	-O $sample.sorted.markdup.recal_data2nd_after.table \\\n";
			print SH "	--known-sites $sample\_filtered_snps2nd.vcf.gz \\\n";
			print SH "	--known-sites $sample\_filtered_indels2nd.vcf.gz && echo \"** $sample.sorted.markdup.recal_data.table done **\" \n";

			print SH "gatk HaplotypeCaller \\\n";
			print SH "	--emit-ref-confidence GVCF \\\n";
			print SH "	-R $var{reference} \\\n";
			print SH "	-ploidy $var{ploidy} \\\n";
			print SH "	-I $sample.sorted.markdup.BQSR2nd.bam \\\n";
			print SH "	-O $sample.HC.gvcf.gz && echo \"** GVCF ${sample}.HC.g.vcf.gz done\" && \\\n";
		  	print SH "rm -f $sample.sorted.markdup.BQSR2nd.bam\n";

		  	print SH "gatk AnalyzeCovariates \\\n";
		  	print SH "	--before-report-file $sample.sorted.markdup.recal_data.table \\\n";
		  	print SH "	--after-report-file $sample.sorted.markdup.recal_data1st_after.table \\\n";
		  	print SH "	--plots-report-file $sample.recalQC.1st.pdf\n";

		  	print SH "gatk AnalyzeCovariates \\\n";
		  	print SH "	--before-report-file $sample.sorted.markdup.recal_data2nd.table \\\n";
		  	print SH "	--after-report-file $sample.sorted.markdup.recal_data2nd_after.table \\\n";
		  	print SH "	--plots-report-file $sample.recalQC.2nd.pdf && echo \"** variant calling done **\" > $var{shpath}/$sample.variant_calling.finished.txt\n";
	  	}
		close SH;
		print CL "sh $var{shpath}/$sample.variant_calling.sh 1>$var{shpath}/$sample.variant_calling.sh.o 2>$var{shpath}/$sample.variant_calling.sh.e\n";
	}
	close CL;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_variant_calling_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=10G,num_proc=1 -binding linear:1' -m 100 $var{shpath}/cmd_variant_calling.list` unless (defined $opts{skipsh});

	while(1){
		sleep(10);
		my $flag_finish = 1; 
		foreach my $sample (keys %samplelist){
			next if (defined $samplelist{$sample}{finish_flag});
			if(-e "$var{shpath}/$sample.variant_calling.finished.txt"){
				next;
			}else{
				$flag_finish = 0;
			}
		}
		my $datestring = localtime();
		print "waiting for sample.variant_calling to be done at $datestring\n";
		last if($flag_finish == 1);
	}
}

sub JointCalling{

	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	if ( !-d "$var{outpath}/JointCalling/" ) {
		make_path "$var{outpath}/JointCalling/" or die "Failed to create path: $var{outpath}/JointCalling/";
	}

	open CL, ">$var{shpath}/cmd_joint_calling.list";

	open SH, ">$var{shpath}/joint_calling.sh";
	## Joint genotyping
	## First, merge all the gvcf results, then perform GenotypeGVCFs
	my $sample_gvcfs = "";
	foreach my $sample (keys %samplelist){
		$sample_gvcfs .= "	-V $var{outpath}/$sample/$sample.HC.gvcf.gz \\\n";
	}
	print SH "#!/bin/sh\ncd $var{outpath}\n";

	print SH "gatk CombineGVCFs \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "$sample_gvcfs";
	print SH "	-O $var{outpath}/JointCalling/JointCalling.HC.g.vcf.gz && echo \"** JointCalling.HC.g.vcf.gz done ** \"\n";

	print SH "gatk GenotypeGVCFs \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "	-ploidy $var{ploidy} \\\n";
	print SH "	-V $var{outpath}/JointCalling/JointCalling.HC.g.vcf.gz \\\n";
	print SH "	-O $var{outpath}/JointCalling/JointCalling.HC.vcf.gz && echo \"** JointCalling.HC.vcf.gz done ** \"\n";

	close SH;
	print CL "sh $var{shpath}/joint_calling.sh 1>$var{shpath}/joint_calling.sh.o 2>$var{shpath}/joint_calling.sh.e\n";
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_joint_calling_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=10G,num_proc=1 -binding linear:1' -m 100 $var{shpath}/cmd_joint_calling.list` unless (defined $opts{skipsh});
}

sub FreebayesCalling {

	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	if ( !-d "$var{outpath}/FreebayesCalling/" ) {
		make_path "$var{outpath}/FreebayesCalling/" or die "Failed to create path: $var{outpath}/FreebayesCalling/";	
	}

	if ( !-d "$var{outpath}/FreebayesCalling/SplitScaffolds/" ) {
		make_path "$var{outpath}/FreebayesCalling/SplitScaffolds/" or die "Failed to create path: $var{outpath}/FreebayesCalling/SplitScaffolds/";
	}

	if ( !-d "$var{shpath}/FreebayesCalling_SplitScaffolds/" ) {
		make_path "$var{shpath}/FreebayesCalling_SplitScaffolds/" or die "Failed to create path: $var{shpath}/FreebayesCalling_SplitScaffolds/";
	}

	open BAMLIST, ">$var{outpath}/FreebayesCalling/bam.list";
	foreach my $sample (keys %samplelist){
		print BAMLIST "$var{outpath}/$sample/$sample.sorted.markdup.bam\n";
	}
	close BAMLIST;

	open REF, "$var{reference}";
	my %genome;
	while (<REF>){
		chomp;
		if (/\>(\S+)/){
			$genome{$1}=1;
		}
	}

	open CL, ">$var{shpath}/cmd_freebayes_calling.list";

	foreach my $id (keys %genome){
		open SH, ">$var{shpath}/FreebayesCalling_SplitScaffolds/freebayes_calling_$id.sh";
		print SH "freebayes -f $var{reference} -L $var{outpath}/FreebayesCalling/bam.list -p $var{ploidy} -r $id >$var{outpath}/FreebayesCalling/SplitScaffolds/freebayes_joint_calling_$id.vcf";
		close SH;
		print CL "sh $var{shpath}/FreebayesCalling_SplitScaffolds/freebayes_calling_$id.sh 1>$var{shpath}/FreebayesCalling_SplitScaffolds/freebayes_calling_$id.sh.o 2>$var{shpath}/FreebayesCalling_SplitScaffolds/freebayes_calling_$id.sh.e\n";
	}

	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_freebayes_calling_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=2G,num_proc=1 -binding linear:1' -m 100 $var{shpath}/cmd_freebayes_calling.list` unless (defined $opts{skipsh});
	
}

1;