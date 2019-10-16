package PopGenome_Mtphylogeny;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use FindBin '$Bin';
use YAML::Tiny;
#use Bio::SeqIO;
use lib "$Bin/lib";
#use lib '/root/miniconda3/lib/site_perl/5.26.2';
use PopGenome_Shared;

############################
#			  			   #
#    MTPHYLOGENY           #
#			               #
############################
sub MTPHYLOGENY{
	my ($yml_file,$skipsh) = @_;
	my $yaml = $yml_file;
	#my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my %samplelist = %{$cfg{fqdata}};

	die "please add mt genome information into configuration file" unless (defined $cfg{mtref}{db});
	
	foreach my $temp_ref(keys %{$cfg{mtref}{db}}){
		my $outpath = "$cfg{args}{outdir}/01.QualityControl/mt_genome_mapping.$temp_ref"; 
		if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";} 
		my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/mt_genome_mapping.$temp_ref";
		if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

		die "please add mt genome path into configuration file" unless (defined $cfg{ref}{db}{$temp_ref}{path});
		my $reference = $cfg{mtref}{db}{$temp_ref}{path};
		die "$reference does not exists" unless (-e $reference);

		open CL, ">$shpath/cmd_mt_genome_mapping.list";
		foreach my $sample (keys %samplelist){
			my $sample_outpath="$outpath/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}

			open SH, ">$shpath/$sample.mt_genome_mapping.sh";
			if(-e "$shpath/$sample.mt_genome_mapping.01_mapping.finished.txt"){`rm $shpath/$sample.mt_genome_mapping.01_mapping.finished.txt`;}	

			print SH "#!/bin/sh\ncd $sample_outpath\n";
			foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
				if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "PE"){
					print SH "bwa mem $reference ../../read_filtering/$sample/$lib\_1.filt.fq.gz ../../read_filtering/$sample/$lib\_2.filt.fq.gz -t $cfg{args}{threads} -R \"\@RG\\tID:$lib\\tSM:$sample\\tLB:$lib\\tPL:$samplelist{$sample}{rawdata}{$lib}{PL}\"\| samtools view -bS -@ $cfg{args}{threads} -F 4 - -o $lib\_filt.bam && \\\n";
				}
				elsif($samplelist{$sample}{rawdata}{$lib}{Flag} eq "SE"){
					print SH "bwa mem $reference ../../read_filtering/$sample/$lib\_1.filt.fq.gz -t $cfg{args}{threads} -R \"\@RG\\tID:$lib\\tSM:$sample\\tLB:$lib\\tPL:$samplelist{$sample}{rawdata}{$lib}{PL}\"\| samtools view -bS -@ 10 -F 4 - -o $lib\_filt.bam && \\\n";
				}
				#summarise statistics for each library bam file 
				print SH "samtools stats -@ $cfg{args}{threads} $lib\_filt.bam 1>$lib\_filt.bamstat.txt 2>$lib\_filt.bamstat.txt.e && echo \"** $lib\_filt.bamstat.txt done **\"\n";
				#then sort each library bam file 
				print SH "samtools sort -@ $cfg{args}{threads} $lib\_filt.bam -o $lib\_filt.sort.bam --output-fmt BAM && \\\n";
				#then remove the unsorted bam file
				print SH "rm -f $lib\_filt.bam\n";			
			}

			#when there is only one library/lane for each sample
			if (keys %{$samplelist{$sample}{rawdata}} == 1){
				foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
					print SH "mv $lib\_filt.sort.bam $sample.sorted.bam\n";}
				print SH "samtools stats -@ $cfg{args}{threads} $sample.sorted.bam 1>bam.stats.txt 2>bam.stats.txt.e && echo \"** finish mt_genome_mapping **\" > $shpath/$sample.mt_genome_mapping.finished.txt\n";
			}

			#when there is more than one library/lane for each sample
			elsif (keys %{$samplelist{$sample}{rawdata}} > 1){
				print SH "samtools merge -nr -@ $cfg{args}{threads} $sample.sorted.bam *_filt.sort.bam && echo \"** $sample.sort.bam done **\" && rm -f *_filt.sort.bam\n";
				#print SH "samtools sort -@ $cfg{args}{threads} $sample.bam -o $sample.sorted.bam --output-fmt BAM && echo \"** $sample.sorted.bam done **\" && rm -f $sample.bam\n";
				print SH "samtools stats -@ $cfg{args}{threads} $sample.sorted.bam 1>bam.stats.txt 2>bam.stats.txt.e && echo \"** finish mt_genome_mapping **\" > $shpath/$sample.mt_genome_mapping.finished.txt\n";
			}
			close SH;
			print CL "sh $shpath/$sample.mt_genome_mapping.sh 1>$shpath/$sample.mt_genome_mapping.sh.o 2>$shpath/$sample.mt_genome_mapping.sh.e \n";
		}
		close CL;
		my $threads = $cfg{args}{threads};
		`perl $Bin/lib/qsub.pl -d $shpath/cmd_mt_genome_mapping_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=$cfg{args}{mem},num_proc=$threads -binding linear:1' -m 100 -r $shpath/cmd_mt_genome_mapping.list` unless ($skipsh ==1);

		my $flag_finish = 0;
		while(1){
			sleep(10);
			my $sample_number = keys %samplelist;
			foreach my $sample (keys %samplelist){
				if(-e "$shpath/$sample.mt_genome_mapping.finished.txt"){$flag_finish +=1;}
			}
			last if($flag_finish == $sample_number);
		}
		open CL, ">$shpath/cmd_mt_genome_variant_calling.list";
		foreach my $sample (keys %samplelist){
			my $sample_outpath="$outpath/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}
			open SH, ">$shpath/$sample.mt_genome_variant_calling.sh";

			print SH "#!/bin/sh\ncd $sample_outpath\n";

			# MarkDuplicates
			print SH "gatk MarkDuplicates \\\n";
		  	print SH "	--INPUT $sample.sorted.bam \\\n";
		  	print SH "	--OUTPUT $sample.sorted.markdup.bam \\\n";
		  	print SH "	--METRICS_FILE $sample.sorted.markdup_metrics.txt && \\\n";
		  	print SH "rm -f $sample.sorted.bam && \\\n";
		  	print SH "echo \"** $sample.sorted.markdup.bam done **\" \n";
		  	print SH "samtools index $sample.sorted.markdup.bam && echo \"** $sample.sorted.markdup.bam index done **\" \n";

		  	print SH "gatk HaplotypeCaller \\\n";
			print SH "	--emit-ref-confidence GVCF \\\n";
			print SH "	-R $reference \\\n";
			print SH "	-ploidy 1 \\\n";
			print SH "	-I $sample.sorted.markdup.bam \\\n";
			print SH "	-O $sample.HC.gvcf.gz && echo \"** GVCF ${sample}.HC.g.vcf.gz done\" && echo \"** finish mt_genome_variant_calling*\" > $shpath/$sample.mt_genome_variant_calling.finished.txt **\" \n";

			close SH;
			print CL "sh $shpath/$sample.mt_genome_variant_calling.sh 1>$shpath/$sample.mt_genome_variant_calling.sh.o 2>$shpath/$sample.mt_genome_variant_calling.sh.e \n";
		}
		
		close CL;
		my $threads = $cfg{args}{threads};
		`perl $Bin/lib/qsub.pl -d $shpath/cmd_mt_genome_variant_calling_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=1 -binding linear:1' -m 100 -r $shpath/cmd_mt_genome_variant_calling.list` unless ($skipsh ==1);
	
		my $flag_finish = 0;
		while(1){
			sleep(10);
			my $sample_number = keys %samplelist;
			foreach my $sample (keys %samplelist){
				if(-e "$shpath/$sample.mt_genome_variant_calling.finished.txt"){$flag_finish +=1;}
			}
			last if($flag_finish == $sample_number);
		}


		if ( !-d "$outpath/Joint_calling" ) {make_path "$outpath/Joint_calling" or die "Failed to create path: $outpath/Joint_calling";} 
		
		open CL, ">$shpath/cmd_mt_genome_joint_calling.list";

		print SH "#!/bin/sh\ncd $outpath/Joint_calling\n";
		open SH, ">$shpath/mt_genome_joint_calling.sh";
		## Joint genotyping
		## First, merge all the gvcf results, then perform GenotypeGVCFs
		my $sample_gvcfs = "";
		foreach my $sample (keys %samplelist){
			$sample_gvcfs .= "	-V $outpath/$sample/$sample.HC.gvcf.gz \\\n";
		}

		print SH "gatk CombineGVCFs \\\n";
		print SH "	-R $reference \\\n";
		print SH "$sample_gvcfs";
		print SH "	-O $outpath/Joint_calling/Joint.HC.g.vcf.gz && echo \"** Combined.HC.g.vcf.gz done ** \"\n";

		print SH "gatk GenotypeGVCFs \\\n";
		print SH "	-R $reference \\\n";
		print SH "	-ploidy $cfg{args}{ploidy} \\\n";
		print SH "	-V $outpath/Joint_calling/Joint.HC.g.vcf.gz \\\n";
		print SH "	-O $outpath/Joint_calling/Joint.HC.vcf.gz && echo \"** Combined.HC.vcf.gz done ** \"\n";

		close SH;
		print CL "sh $shpath/mt_genome_joint_calling.sh 1>$shpath/mt_genome_joint_calling.sh.o 2>$shpath/mt_genome_joint_calling.sh.e \n";
		close CL;

		`perl $Bin/lib/qsub.pl -d $shpath/cmd_mt_genome_joint_calling_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=1 -binding linear:1' -m 100 -r $shpath/cmd_mt_genome_joint_calling.list` unless ($skipsh ==1);

	}	

	# create this yaml object
	#$yaml = YAML::Tiny->new( \%cfg );
	# Save both documents to a file
	#$yml_file = s/\.yml$//g;
	#$yml_file = $yml_file."_mapping_modified.yml";
	#$yaml->write( $yml_file );
}

1;