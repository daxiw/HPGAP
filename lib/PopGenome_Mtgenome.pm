package PopGenome_Mtgenome;
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
		'threads=s',
		'samplelist=s',
		'downsize=s',
		'mt_genome_mapping',
		'mt_genome_variant_calling',
		'mt_genome_phylogeny',
		'help',
		'skipsh');

	die "please provide the correct configuration file" unless ((defined $opts{config}) && (-e $opts{config}));

	if (defined $opts{allsteps}){
		$opts{mt_genome_mapping} = 1;
		$opts{mt_genome_variant_calling} = 1;
		$opts{mt_genome_phylogeny} = 1;
	}

	my $yaml = YAML::Tiny->read( $opts{config} );
	my %cfg = %{$yaml->[0]};
	my %samplelist_ori = %{$cfg{fqdata}};
	my %samplelist = %samplelist_ori;

	if (defined $opts{samplelist}){
		my %selected_sample;
		open IN, $opts{samplelist};
		while (<IN>){
			/(\S+)/;
			$selected_sample{$1}=1;
		}
		close IN;
		foreach my $id (keys %samplelist){
			unless (exists $selected_sample{$id}){
				delete $samplelist{$id};
			}
		}
	}

	$var{samplelist}=\%samplelist;
	$var{cfg}=\%cfg;

	#set ploidy
	$var{ploidy} = 2;

	#set the number of threads
	if (defined $opts{threads}){
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
		$var{reference} = $cfg{mtref}{db}{$temp_ref}{path};

		die "$var{reference} does not exists" unless (-e $var{reference});

		open IN, $var{reference};
		$var{temp_ref}{length}=0;
		while (<IN>){
			if(/\>/){next;}
			elsif(/\w+/){
				$var{temp_ref}{length} += length;
			}
		}
		close IN;

		if (defined $opts{mt_genome_mapping}){ & MtGenomeMapping (\%var,\%opts);}
		if (defined $opts{mt_genome_variant_calling}){ & MtGenomeVariantCalling (\%var,\%opts);}
		if (defined $opts{mt_genome_phylogeny}){ & MtGenomePhylogeny (\%var,\%opts);}		
	}	
}

sub MtGenomeMapping {
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	open CL, ">$var{shpath}/cmd_mt_genome_mapping.list";
	my $sample_number = 0;
	foreach my $sample (keys %samplelist){
		my $sample_outpath="$var{outpath}/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}

		### skip finished samples 
		$samplelist{$sample}{finish_flag}="finished";
		foreach my $readgroup (keys %{$samplelist{$sample}{cleandata}}){
			if ((-e "$sample_outpath/$readgroup\_filt.bamstat.txt") && (!defined $opts{overwrite})){
				$samplelist{$sample}{cleandata}{$readgroup}{finish_flag}="finished";
			}
			else{
				$samplelist{$sample}{finish_flag}="unfinished";
				$samplelist{$sample}{cleandata}{$readgroup}{finish_flag}="unfinished";
			#	`rm -f $var{shpath}/$sample.variant_calling.finished.txt`;
			}
		}
		next if($samplelist{$sample}{finish_flag} eq "finished");

		open SH, ">$var{shpath}/$sample.mt_genome_mapping.sh";		
		print SH "#!/bin/sh\ncd $sample_outpath\n";
		foreach my $readgroup (keys %{$samplelist{$sample}{cleandata}}){
			if($samplelist{$sample}{cleandata}{$readgroup}{Flag} eq "PE"){
				print SH "bwa mem $var{reference} $samplelist{$sample}{cleandata}{$readgroup}{fq1} $samplelist{$sample}{cleandata}{$readgroup}{fq2} -t $var{threads} -R \"\@RG\\tID:$readgroup\\tSM:$sample\\tLB:$readgroup\\tPL:$samplelist{$sample}{cleandata}{$readgroup}{PL}\"\| samtools view -bS -@ $cfg{args}{threads} -F 4 - -o $readgroup\_filt.bam && \\\n";
			}
			elsif($samplelist{$sample}{cleandata}{$readgroup}{Flag} eq "SE"){
				print SH "bwa mem $var{reference} $samplelist{$sample}{cleandata}{$readgroup}{fq1} -t $var{threads} -R \"\@RG\\tID:$readgroup\\tSM:$sample\\tLB:$readgroup\\tPL:$samplelist{$sample}{cleandata}{$readgroup}{PL}\"\| samtools view -bS -@ 10 -F 4 - -o $readgroup\_filt.bam && \\\n";
			}
			#summarise statistics for each library bam file 
			print SH "samtools stats -@ $var{threads} $readgroup\_filt.bam 1>$readgroup\_filt.bamstat.txt 2>$readgroup\_filt.bamstat.txt.e && echo \"** $readgroup\_filt.bamstat.txt done **\"\n";
			#then sort each library bam file 
			print SH "samtools sort -@ $var{threads} $readgroup\_filt.bam -o $readgroup\_filt.sort.bam --output-fmt BAM && \\\n";
			#then remove the unsorted bam file
			print SH "rm -f $readgroup\_filt.bam\n";
		}

		#when there is only one library/lane for each sample
		if (keys %{$samplelist{$sample}{cleandata}} == 1){
			foreach my $readgroup (keys %{$samplelist{$sample}{cleandata}}){
				print SH "mv $readgroup\_filt.sort.bam $sample.sorted.bam\n";
			}
		}

		#when there is more than one library/lane for each sample
		if (keys %{$samplelist{$sample}{cleandata}} > 1){
			print SH "samtools merge -f -@ $var{threads} $sample.sorted.bam *_filt.sort.bam && echo \"** $sample.sorted.bam done **\" && rm -f *_filt.sort.bam\n";
		}

		if (defined $opts{downsize}){
			my $required_coverage = $opts{downsize}*$var{temp_ref}{length};
			print SH "samtools stats -@ $var{threads} $sample.sorted.bam 1>$sample.sorted.bam.stats.txt 2>$sample.sorted.bam.stats.error\n";
			print SH "a=\`cat $sample.sorted.bam.stats.txt|perl -ne \'if(/cigar\\):\\s+(\\d+)/){print $required_coverage/\$1;}\'\`\n";
			print SH "if [ $a -le 1 ]\n";
			print SH "then\n";
			print SH "	mv $sample.sorted.bam $sample.sorted.ori.bam && \\\n";
			print SH "	samtools view -s \$a $sample.sorted.ori.bam -o $sample.sorted.bam\n";
			print SH "fi\n";
		}

		print SH "bedtools genomecov -d -ibam $sample.sorted.bam > $sample.genomecov\n";

		print SH "gatk MarkDuplicates \\\n";
	  	print SH "	--INPUT $sample.sorted.bam \\\n";
	  	print SH "	--OUTPUT $sample.sorted.markdup.bam \\\n";
	  	print SH "	--METRICS_FILE $sample.sorted.markdup_metrics.txt && \\\n";
	  	print SH "rm -f $sample.sorted.bam && \\\n";
	  	print SH "echo \"** $sample.sorted.markdup.bam done **\" \n";
	  	print SH "samtools index $sample.sorted.markdup.bam && \\\n";
	  	print SH "echo \"** $sample.sorted.markdup.bam index done **\" \n";

	  	print SH "samtools stats -@ $var{threads} $sample.sorted.markdup.bam 1>$sample.bam.stats.txt 2>$sample.bam.stats.error && echo \"** bam.stats.txt done **\" > $var{shpath}/$sample.readmapping.finished.txt\n";

		close SH;
		print CL "sh $var{shpath}/$sample.mt_genome_mapping.sh 1>$var{shpath}/$sample.mt_genome_mapping.sh.o 2>$var{shpath}/$sample.mt_genome_mapping.sh.e \n";
	}
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_mt_genome_mapping_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=2G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/cmd_mt_genome_mapping.list` unless (defined $opts{skipsh});
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

	open BAMLIST, ">$var{outpath}/FreebayesCalling/bam.list";
	if ( !-d "$var{outpath}/FreebayesCalling") {make_path "$var{outpath}/FreebayesCalling" or die "Failed to create path: $var{outpath}/FreebayesCalling";}
	foreach my $sample (keys %samplelist){
		if ( !-d "$var{outpath}/$sample") {make_path "$var{outpath}/$sample" or die "Failed to create path: $var{outpath}/$sample";}
		if ( -e "$var{outpath}/$sample/$sample.genomecov"){
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
				if ($a[2]>=1000){ $depth{1000}++; }
				elsif($a[2]>=100){ $depth{100}++; }
				elsif($a[2]>=50){ $depth{50}++; }
				elsif($a[2]>=10){ $depth{10}++; }
				elsif($a[2]>=10){ $depth{1}++; }
				elsif($a[2]==0){ $depth{0}++; }
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
		print BAMLIST "$var{outpath}/$sample/$sample.sorted.markdup.bam\n";
	}
	close BAMLIST;

	#### freebayes variant calling on mt genomes ###
	open CL, ">$var{shpath}/cmd_mt_genome_freebayes.list";
	open SH, ">$var{shpath}/mt_genome_freebayes.sh";

	print SH "#!/bin/sh\ncd $var{outpath}/FreebayesCalling\n";
	print SH "freebayes -f $var{reference} -L $var{outpath}/FreebayesCalling/bam.list -p $var{ploidy} --standard-filters | vcfsnps >$var{outpath}/FreebayesCalling/freebayes_joint_calling.vcf";
	close SH;
	
	print CL "sh $var{shpath}/mt_genome_freebayes.sh 1>$var{shpath}/mt_genome_freebayes.sh.o 2>$var{shpath}/mt_genome_freebayes.sh.e \n";
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_mt_genome_freebayes_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=1 -binding linear:1' -m 100 -r $var{shpath}/cmd_mt_genome_freebayes.list` unless (defined $opts{skipsh});

}

sub MtGenomePhylogeny{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	if ( !-d "$var{outpath}/Mt_genome_phylogeny" ) {make_path "$var{outpath}/Mt_genome_phylogeny" or die "Failed to create path: $var{outpath}/Mt_genome_phylogeny";}

	my $i = 0;my %name;
	open OT, ">$var{outpath}/Mt_genome_phylogeny/name_map.list";
	foreach my $sample (keys %samplelist){
		print OT "$sample ", "NS$i", "E\n";
		my $newid = "NS$i"."E";
		$name{$newid} = $sample;
		$i++;
	}
	close OT;

	open CL, ">$var{shpath}/cmd_mt_genome_phylogeny.list";

	open SH, ">$var{shpath}/mt_genome_phylogeny.sh";
	print SH "#!/bin/sh\ncd $var{outpath}/Mt_genome_phylogeny\n";

	print SH "bcftools reheader --samples $var{outpath}/Mt_genome_phylogeny/name_map.list -o $var{outpath}/FreebayesCalling/freebayes_joint_calling_rename.vcf $var{outpath}/FreebayesCalling/freebayes_joint_calling.vcf\n";
	print SH "rm -rf $var{outpath}/Mt_genome_phylogeny/RAxML_*\n";

	print SH "cat $var{outpath}/FreebayesCalling/freebayes_joint_calling_rename.vcf|$Bin/Tools/vcf-to-tab >$var{outpath}/Mt_genome_phylogeny/mt_genome.tab\n";
	print SH "$Bin/Tools/vcf_tab_to_fasta_alignment.pl -i $var{outpath}/Mt_genome_phylogeny/mt_genome.tab > $var{outpath}/Mt_genome_phylogeny/mt_genome.fasta\n";

	print SH "$Bin/Tools/fasta-to-phylip --input-fasta $var{outpath}/Mt_genome_phylogeny/mt_genome.fasta --output-phy $var{outpath}/Mt_genome_phylogeny/mt_genome.phy\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -s $var{outpath}/Mt_genome_phylogeny/mt_genome.phy -n trees -T $var{threads} -# 20 -p 12345\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -s $var{outpath}/Mt_genome_phylogeny/mt_genome.phy -n boots -T $var{threads} -# 100 -p 23456 -b 23456\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -p 12345 -f b -t RAxML_bestTree.trees -T $var{threads} -z RAxML_bootstrap.boots -n consensus\n";
	print SH "sumtrees.py --percentages --min-clade-freq=0.50 --target=RAxML_bestTree.trees --output=result2.tre RAxML_bootstrap.boots";

	close SH;
	print CL "sh $var{shpath}/mt_genome_phylogeny.sh 1>$var{shpath}/mt_genome_phylogeny.sh.o 2>$var{shpath}/mt_genome_phylogeny.sh.e \n";
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_mt_genome_phylogeny_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/cmd_mt_genome_phylogeny.list` unless (defined $opts{skipsh});

	open my $fh, '<', "$var{outpath}/Mt_genome_phylogeny/result2.tre" or die "error opening $var{outpath}/Mt_genome_phylogeny/result2.tre: $!";
	my $data = do { local $/; <$fh> };

	foreach my $newid (keys %name){
		$data =~ s/($newid)/($name{$newid})/g;
	}
	close $fh;

	open OT, ">$var{outpath}/Mt_genome_phylogeny/result2.final.tre";
	print OT $data;
	close OT;
	
}

1;