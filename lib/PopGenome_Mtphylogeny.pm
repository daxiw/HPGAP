package PopGenome_Mtphylogeny;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use Getopt::Long qw(GetOptionsFromArray);
use FindBin '$Bin';
use YAML::Tiny;
use lib "$Bin/lib";
use PopGenome_Shared;

############################
#			  			   #
#    MTPHYLOGENY           #
#			               #
############################
sub Main{
	my $args = shift; 
	my @args = @{$args};
	my %opts;
	my %var;

	GetOptionsFromArray (\@args, \%opts, 
		'config=s',
		'overwrite',
		'allsteps',
		'mt_genome_mapping',
		'mt_genome_variant_calling',
		'mt_genome_phylogeny',
		'help',
		'skipsh');

	die "please provide the correct configuration file" unless ((defined $opts{config}) && (-e $opts{config}));

	if (defined $opts{allsteps}){
		$opts{mt_genome_mapping} = 1;
		$opts{mt_genome_variant_calling = 1;
		$opts{mt_genome_phylogeny} = 1;
	}

	my $yaml = YAML::Tiny->read( $opts{config} );
	my %cfg = %{$yaml->[0]};
	my %samplelist = %{$cfg{fqdata}};
	

	#set ploidy
	$var{ploidy} = 1;

	#set the number of threads
	if (defined $opts{threads}}){
		$var{threads} = $opts{threads};
	}elsif(defined $cfg{args}{threads}){
		$var{threads} = $cfg{args}{threads};
	}else{
		$var{threads} = 4;
	}

	die "please add mt genome information into the configuration file" unless (defined $cfg{mtref}{db});
	
	foreach my $temp_ref(keys %{$cfg{mtref}{db}}){
		$var{outpath} = "$cfg{args}{outdir}/01.QualityControl/mt_phylogeny.$temp_ref"; 
		if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 
		$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/mt_phylogeny.$temp_ref";
		if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}

		die "please add mt genome path into configuration file" unless (defined $cfg{mtref}{db}{$temp_ref}{path});
		my $var{reference} = $cfg{mtref}{db}{$temp_ref}{path};
		die "$var{reference} does not exists" unless (-e $var{reference});

		if (defined $opts{mt_genome_mapping}){ &MtGenomeMapping (\%var,\%opts);}

		if (defined $opts{mt_genome_variant_calling}){ &MtGenomeVariantCalling (\%var,\%opts);}

		if (defined $opts{mt_genome_phylogeny}){ &ReadReport (\%var,\%opts);}
		#### estimate phylogeny of mt genomes ###
		
	}	

}

sub MtGenomeMapping {
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	open CL, ">$var{shpath}/cmd_mt_genome_mapping.list";
	foreach my $sample (keys %samplelist){
		my $sample_outpath="$var{outpath}/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}

		open SH, ">$var{shpath}/$sample.mt_genome_mapping.sh";
		#if(-e "$var{shpath}/$sample.mt_genome_mapping.finished.txt"){`rm $var{shpath}/$sample.mt_genome_mapping.finished.txt`;}	

		print SH "#!/bin/sh\ncd $sample_outpath\n";

		my ($n_lib, $n_done, $run_flag, $written_flag) = qw (0 0 0 0);

		foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
			$n_lib ++;
			if (-e "$sample_outpath/$lib\_filt.bamstat.txt"){ $n_done++; next unless (defined $opts{overwrite});}

			if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "PE"){
				print SH "bwa mem $var{reference} ../../read_filtering/$sample/$lib\_1.filt.fq.gz ../../read_filtering/$sample/$lib\_2.filt.fq.gz -t $var{threads} -R \"\@RG\\tID:$lib\\tSM:$sample\\tLB:$lib\\tPL:$samplelist{$sample}{rawdata}{$lib}{PL}\"\| samtools view -bS -@ $var{threads} -F 4 - -o $lib\_filt.bam && \\\n";
			}
			elsif($samplelist{$sample}{rawdata}{$lib}{Flag} eq "SE"){
				print SH "bwa mem $var{reference} ../../read_filtering/$sample/$lib\_1.filt.fq.gz -t $var{threads} -R \"\@RG\\tID:$lib\\tSM:$sample\\tLB:$lib\\tPL:$samplelist{$sample}{rawdata}{$lib}{PL}\"\| samtools view -bS -@ 10 -F 4 - -o $lib\_filt.bam && \\\n";
			}
			#summarise statistics for each library bam file 
			print SH "samtools stats -@ $var{threads} $lib\_filt.bam 1>$lib\_filt.bamstat.txt 2>$lib\_filt.bamstat.txt.e && echo \"** $lib\_filt.bamstat.txt done **\"\n";
			#then sort each library bam file 
			print SH "samtools sort -@ $var{threads} $lib\_filt.bam -o $lib\_filt.sort.bam --output-fmt BAM && \\\n";
			#then remove the unsorted bam file
			print SH "rm -f $lib\_filt.bam\n";
			$written_flag ++;		
		}

		$run_flag = 1 if (!(defined $opts{overwrite}) && (-e "$sample_outpath/bam.stats.txt"));

		#when there is only one library/lane for each sample
		if (($n_lib == 1) && ($run_flag == 0)){
			foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
				print SH "mv $lib\_filt.sort.bam $sample.sorted.bam\n";
			}
			print SH "samtools stats -@ $var{threads} $sample.sorted.bam 1>bam.stats.txt 2>bam.stats.txt.e && echo \"** finish mt_genome_mapping **\" > $var{shpath}/$sample.mt_genome_mapping.finished.txt\n";
			$written_flag ++;
		}

		#when there is more than one library/lane for each sample
		elsif (($n_lib > 1) && ($run_flag == 0)){
			print SH "samtools merge -nr -@ $var{threads} $sample.sorted.bam *_filt.sort.bam && echo \"** $sample.sort.bam done **\" && rm -f *_filt.sort.bam\n";
			print SH "samtools stats -@ $var{threads} $sample.sorted.bam 1>bam.stats.txt 2>bam.stats.txt.e && echo \"** finish mt_genome_mapping **\" > $var{shpath}/$sample.mt_genome_mapping.finished.txt\n";
			$written_flag ++;
		}
		elsif (($n_lib > 1) && ($run_flag > 0) && ($n_done!=$n_lib)){
			print SH "mv $sample.sorted.bam $sample.previous.sorted.bam && \\\n";
			print SH "samtools merge -nr -@ $var{threads} $sample.sorted.bam $sample.previous.sorted.bam *_filt.sort.bam && echo \"** $sample.sort.bam done **\" && rm -f *_filt.sort.bam\n";
			print SH "samtools stats -@ $var{threads} $sample.sorted.bam 1>bam.stats.txt 2>bam.stats.txt.e && echo \"** finish mt_genome_mapping **\" > $var{shpath}/$sample.mt_genome_mapping.finished.txt\n";
			$written_flag ++;
		}

		# calculate read depth for each base
		if (($run_flag == 0)||($n_done!=$n_lib)){
			print SH "bedtools genomecov -d -ibam $sample.sorted.bam > $sample.genomecov\n";
			$written_flag ++;
		}

		close SH;
		print CL "sh $var{shpath}/$sample.mt_genome_mapping.sh 1>$var{shpath}/$sample.mt_genome_mapping.sh.o 2>$var{shpath}/$sample.mt_genome_mapping.sh.e \n" if ($written_flag > 0);
	}
	close CL;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_mt_genome_mapping_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=2G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/cmd_mt_genome_mapping.list` unless (defined $opts{skipsh});

	my $flag_finish = 0;
	while(1){
		sleep(20);
		
		my $sample_number = keys %samplelist;

		foreach my $sample (keys %samplelist){

			if(-e "$var{shpath}/$sample.mt_genome_mapping.finished.txt"){$flag_finish +=1;}

		}
		my $datestring = localtime();
		print "waiting for mt_genome_mapping to be done [$flag_finish out of $sample_number samples are finished] at $datestring\n";

		last if($flag_finish == $sample_number);
	}
}

sub MtGenomeVariantCalling{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	## check whether the fai and dict files of reference exist 
	my $ref_name=$var{reference};
	if ($ref_name =~ /\.fasta$/){
		$ref_name =~ s/\.fasta$//; 
	}elsif ($ref_name =~ /\.fa$/){
		$ref_name =~ s/\.fa$//;
	}elsif ($ref_name =~ /\.fas$/){
		$ref_name =~ s/\.fas$//;
	}

	if ( !-e "$var{reference}.fai" ) {
	  		`samtools faidx $var{reference}`;
	}
	if ( !-e "$ref_name.dict" ) {
	  		`picard CreateSequenceDictionary R=$var{reference} O=$ref_name.dict`;
	}

	###############

	open CL, ">$var{shpath}/cmd_mt_genome_variant_calling.list";
	foreach my $sample (keys %samplelist){
		my $sample_outpath="$var{outpath}/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}

		if(-e "$var{shpath}/$sample.mt_genome_variant_calling.finished.txt"){`rm -f $var{shpath}/$sample.mt_genome_variant_calling.finished.txt`;}	

		if (-e "$var{outpath}/$sample/$sample.genomecov"){
			open IN, "$var{outpath}/$sample/$sample.genomecov";
			my $n_base = 0;
			my $n_sum = 0;
			my %depth;
			$depth{1000}=0;
			$depth{100}=0;
			$depth{50}=0;
			$depth{10}=0;
			$depth{1}=0;
			$depth{0}=0;

			while (<IN>){
				$n_base ++;
				my @a = split /\t/;
				$n_sum += $a[2];
				if ($a[2]>=1000){
					$depth{1000}++;
				}elsif($a[2]>=100){
					$depth{100}++;
				}elsif($a[2]>=50){
					$depth{50}++;
				}elsif($a[2]>=10){
					$depth{10}++;
				}elsif($a[2]>=10){
					$depth{1}++;
				}elsif($a[2]==0){
					$depth{0}++;
				}
			}
			close IN;

			open OT, ">$var{outpath}/$sample/$sample.genomecov.summary.txt";
			print OT "SampleID\t0\t1-9\t10-49\t50-99\t100-999\t>1000\tn_bases\tn_sum\n";
			print OT $sample, "\t";
			print OT $depth{0}, "\t";
			print OT $depth{1}, "\t";
			print OT $depth{10}, "\t";
			print OT $depth{50}, "\t";
			print OT $depth{100}, "\t";
			print OT $depth{1000}, "\t";
			print OT $n_base, "\t";
			print OT $n_sum/$n_base, "\n";
			close OT;
		}

		open SH, ">$var{shpath}/$sample.mt_genome_variant_calling.sh";

		print SH "#!/bin/sh\ncd $sample_outpath\n";

		# MarkDuplicates
		print SH "gatk MarkDuplicates \\\n";
	  	print SH "	--INPUT $sample.sorted.bam \\\n";
	  	print SH "	--OUTPUT $sample.sorted.markdup.bam \\\n";
	  	print SH "	--METRICS_FILE $sample.sorted.markdup_metrics.txt && \\\n";
	  	print SH "rm -f $sample.sorted.bam && \\\n";
	  	print SH "samtools index $sample.sorted.markdup.bam && echo \"** $sample.sorted.markdup.bam index done **\" \n";

	  	print SH "gatk HaplotypeCaller \\\n";
		print SH "	--emit-ref-confidence GVCF \\\n";
		print SH "	-R $var{reference} \\\n";
		print SH "	-ploidy $var{ploidy} \\\n";
		print SH "	-I $sample.sorted.markdup.bam \\\n";
		print SH "	-O $sample.HC.gvcf.gz && echo \"** GVCF ${sample}.HC.g.vcf.gz done\" && \\\n"; 
		print SH "echo \"** finish mt_genome_variant_calling **\" > $var{shpath}/$sample.mt_genome_variant_calling.finished.txt \n";

		close SH;
		print CL "sh $var{shpath}/$sample.mt_genome_variant_calling.sh 1>$var{shpath}/$sample.mt_genome_variant_calling.sh.o 2>$var{shpath}/$sample.mt_genome_variant_calling.sh.e \n";
	}
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_mt_genome_variant_calling_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=1 -binding linear:1' -m 100 -r $var{shpath}/cmd_mt_genome_variant_calling.list` unless (defined $opts{skipsh});

	$flag_finish = 0;
	while(1){
		sleep(20);
		
		my $sample_number = keys %samplelist;
		foreach my $sample (keys %samplelist){
			if(-e "$var{shpath}/$sample.mt_genome_variant_calling.finished.txt"){$flag_finish +=1;}
		}

		my $datestring = localtime();
		print "waiting for mt_genome_variant_calling to be done [$flag_finish out of $sample_number samples are finished] at $datestring\n";
		last if($flag_finish == $sample_number);
	}

	#### joint variant calling on mt genomes ###
	if ( !-d "$var{outpath}/Joint_calling" ) {
		make_path "$var{outpath}/Joint_calling" or die "Failed to create path: $var{outpath}/Joint_calling";
	} 		
	if(-e "$var{shpath}/mt_genome_joint_calling.finished.txt"){`rm -f $var{shpath}/mt_genome_joint_calling.finished.txt`;}

	open CL, ">$var{shpath}/cmd_mt_genome_joint_calling.list";

	print SH "#!/bin/sh\ncd $var{outpath}/Joint_calling\n";
	open SH, ">$var{shpath}/mt_genome_joint_calling.sh";
	## Joint genotyping
	## First, merge all the gvcf results, then perform GenotypeGVCFs
	my $sample_gvcfs = "";
	foreach my $sample (keys %samplelist){
		$sample_gvcfs .= "	-V $var{outpath}/$sample/$sample.HC.gvcf.gz \\\n";
	}

	print SH "gatk CombineGVCFs \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "$sample_gvcfs";
	print SH "	-O $var{outpath}/Joint_calling/Joint.HC.g.vcf.gz && echo \"** Joint.HC.g.vcf.gz done ** \"\n";

	print SH "gatk GenotypeGVCFs \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "	-ploidy $var{ploidy} \\\n";
	print SH "	-V $var{outpath}/Joint_calling/Joint.HC.g.vcf.gz \\\n";
	print SH "	-O $var{outpath}/Joint_calling/Joint.HC.vcf.gz && echo \"** Joint.HC.vcf.gz **\"\n";

	print SH "gatk SelectVariants \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "	-V $var{outpath}/Joint_calling/Joint.HC.vcf.gz \\\n";
	print SH "	--select-type-to-include SNP \\\n";
	print SH "	-O $var{outpath}/Joint_calling/Joint_raw_snps1st.vcf && echo \"** Joint_raw_snps1st.vcf done\" && \\\n";

	print SH "gatk VariantFiltration \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "	-V $var{outpath}/Joint_calling/Joint_raw_snps1st.vcf \\\n";
	print SH "	--filter-expression \"$cfg{step1}{variant_filtering}{snp}\" \\\n";
	print SH "	--filter-name \"my_snp_filter\" \\\n";
	print SH "	-O $var{outpath}/Joint_calling/Joint_filtered_snps1st.vcf && echo \"** Joint_filtered_snps1st done\" \n";

	print SH "gatk SelectVariants \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "	-V $var{outpath}/Joint.HC.vcf.gz \\\n";
	print SH "	--select-type-to-include INDEL \\\n";
	print SH "	-O $var{outpath}/Joint_calling/Joint_raw_indels1st.vcf && echo \"** Joint_raw_indels1st.vcf done\" && \\\n";

	print SH "gatk VariantFiltration \\\n";
	print SH "	-R $var{reference} \\\n";
	print SH "	-V $var{outpath}/Joint_calling/Joint_raw_indels1st.vcf \\\n";
	print SH "	--filter-expression \"$cfg{step1}{variant_filtering}{indel}\" \\\n";
	print SH "	--filter-name \"my_indel_filter\" \\\n";
	print SH "	-O $var{outpath}/Joint_calling/Joint_filtered_indels1st.vcf && echo \"** finish mt_genome_joint_calling **\" > $var{shpath}/mt_genome_joint_calling.finished.txt \n";

	close SH;
	print CL "sh $var{shpath}/mt_genome_joint_calling.sh 1>$var{shpath}/mt_genome_joint_calling.sh.o 2>$var{shpath}/mt_genome_joint_calling.sh.e \n";
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_mt_genome_joint_calling_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=1 -binding linear:1' -m 100 -r $var{shpath}/cmd_mt_genome_joint_calling.list` unless ($opts{skipsh} ==1);


	$flag_finish = 0;
	while(1){
		sleep(20);
		my $datestring = localtime();
		print "waiting for mt_genome_joint_calling to be done at $datestring\n";
		if(-e "$var{shpath}/mt_genome_joint_calling.finished.txt"){last;}
	}
}

sub MtGenomePhylogeny{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};

	if ( !-d "$var{outpath}/Mt_genome_phylogeny" ) {make_path "$var{outpath}/Mt_genome_phylogeny" or die "Failed to create path: $var{outpath}/Mt_genome_phylogeny";}

	open CL, ">$var{shpath}/cmd_mt_genome_phylogeny.list";

	print SH "#!/bin/sh\ncd $var{outpath}/Mt_genome_phylogeny\n";
	open SH, ">$var{shpath}/mt_genome_phylogeny.sh";

	print SH "cd $var{outpath}/Mt_genome_phylogeny\n";
	print SH "rm -rf $var{outpath}/Mt_genome_phylogeny/RAxML_*\n";
	
	print SH "zcat $var{outpath}/Joint_calling/Joint.HC.g.vcf.gz|vcf-to-tab >$var{outpath}/Mt_genome_phylogeny/mt_genome.tab\n";
	print SH "vcf_tab_to_fasta_alignment.pl -i $var{outpath}/Mt_genome_phylogeny/mt_genome.tab > $var{outpath}/mt_genome.fasta\n";

	print SH "seqmagick convert $var{outpath}/Mt_genome_phylogeny/mt_genome.fasta $var{outpath}/Mt_genome_phylogeny/mt_genome.phy\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -s $var{outpath}/Mt_genome_phylogeny/mt_genome.phy -n trees -T 24 -# 20 -p 12345\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -s $var{outpath}/Mt_genome_phylogeny/mt_genome.phy -n boots -T 24 -# 100 -p 23456 -b 23456\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -p 12345 -f b -t RAxML_bestTree.trees -T 2 -z RAxML_bootstrap.boots -n consensus\n";
	#print SH "sumtrees.py --percentages --min-clade-freq=0.50 --target=RAxML_bestTree.trees --output=result2.tre RAxML_bootstrap.boots";

	close SH;
	print CL "sh $var{shpath}/mt_genome_phylogeny.sh 1>$var{shpath}/mt_genome_phylogeny.sh.o 2>$var{shpath}/mt_genome_phylogeny.sh.e \n";
	close CL;

	#`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_mt_genome_phylogeny_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=6 -binding linear:1' -m 100 -r $var{shpath}/cmd_mt_genome_phylogeny.list` unless (defined $opts{skipsh});
}

1;