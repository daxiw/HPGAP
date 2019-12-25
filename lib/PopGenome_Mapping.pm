package PopGenome_Mapping;
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
		'outcfg=s',
		'threads=s',
		'downsize=s',
		'read_mapping',
		'mapping_report',
		'reference_selection',
		'help',
		'skipsh');

	die "please provide the correct configuration file" unless ((defined $opts{config}) && (-e $opts{config}));
	my %cfg = %{YAML::Tiny->read( $opts{config} )->[0]};

	if (defined $opts{allsteps}){
		$opts{read_mapping} = 1;
		$opts{mapping_report} = 1;
		$opts{reference_selection} = 1;
	}

	my %mapping;
	foreach my $temp_ref(keys %{$cfg{ref}{db}}){

		my %var = %{PopGenome_Shared::CombineCfg("$Bin/lib/parameter.yml",\%opts, "01.QualityControl/read_mapping.$temp_ref")};
		
		die "please add genome path into configuration file" unless (defined $cfg{ref}{db}{$temp_ref}{path});
		$var{reference} = $cfg{ref}{db}{$temp_ref}{path};
		die "$var{reference} does not exists" unless (-e $var{reference});

		open IN, $var{reference};
		$var{temp_ref}{name}=$temp_ref;
		$var{temp_ref}{length}=0;
		while (<IN>){
			if(/\>/){next;}
			elsif(/\w+/){
				$var{temp_ref}{length} += length;
			}
		}
		close IN;

		$var{mem} = $var{threads}."G";

		if (defined $opts{read_mapping}){ $mapping{$temp_ref} = & ReadMapping (\%var,\%opts);}
		if (defined $opts{mapping_report}){ & MappingReport (\%var,\%opts);}
	}

	if (defined $opts{reference_selection}){ & SelectReference (\%var,\%opts);}
}

############################
#			   #
#    Step 1b Mapping       #
#			   #
############################
sub ReadMapping {
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};
	
	open CL, ">$var{shpath}/cmd_readmapping.list";
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

		open SH, ">$var{shpath}/$sample.readmapping.sh";		
		print SH "#!/bin/sh\ncd $sample_outpath\n";
		foreach my $readgroup (keys %{$samplelist{$sample}{cleandata}}){

			next if ($samplelist{$sample}{cleandata}{$readgroup}{finish_flag} eq "finished");

			if($samplelist{$sample}{cleandata}{$readgroup}{Flag} eq "PE"){
				print SH "bwa mem $var{reference} $samplelist{$sample}{cleandata}{$readgroup}{fq1} $samplelist{$sample}{cleandata}{$readgroup}{fq2} -t $var{threads} -R \"\@RG\\tID:$readgroup\\tSM:$sample\\tLB:$readgroup\\tPL:$samplelist{$sample}{cleandata}{$readgroup}{PL}\"\| samtools view -bS -@ $var{threads} -F 4 - -o $readgroup\_filt.bam && \\\n";
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
		print CL "sh $var{shpath}/$sample.readmapping.sh 1>$var{shpath}/$sample.readmapping.sh.o 2>$var{shpath}/$sample.readmapping.sh.e \n";
	}
	close CL;

	`perl $Bin/lib/qsub.pl -r -d $var{shpath}/cmd_readmapping_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=$var{mem},num_proc=$var{threads} -binding linear:1' -m 100 $var{shpath}/cmd_readmapping.list` unless (defined $opts{skipsh});

	return (\%samplelist);

}

sub MappingReport {
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	my $report_outpath="$var{outpath}/Report"; 
	if ( !-d $report_outpath ) {make_path $report_outpath or die "Failed to create path: $report_outpath";}
	my $report_sample_outpath="$var{outpath}/Report/Samples"; 
	if ( !-d $report_sample_outpath ) {make_path $report_sample_outpath or die "Failed to create path: $report_outpath";}

	open FOT, "$var{outpath}/../read_filtering/Report/sample_quality_summary.xls";
	my %sample_summary;
	while (<FOT>) {
		chomp;
		next if (/ampleID/);
		my @a = split /\t/;
		$sample_summary{$a[0]}{"raw_data"} = $a[2];
		$sample_summary{$a[0]}{"clean_data"} = $a[4];
		print $a[0],"\t",$sample_summary{$a[0]}{"clean_data"},"\n";
		if ($a[2] > 0) {
			$sample_summary{$a[0]}{"percentage_of_cleandata"} = $a[4]/$a[2];
		} else {
			$sample_summary{$a[0]}{"percentage_of_cleandata"} = 0;
		}
	}
	close FOT;

	open SOT, ">$var{outpath}/Report/sample_mapping_summary.xls";

	print SOT "sampleID","\t";
	print SOT "reference","\t";
	print SOT "percentage_of_cleandata","\t";
	print SOT "clean_data","\t";
	print SOT "mean_covreage","\t";
	print SOT "mapping_rate","\t";
	print SOT "covered_regions","\t";
	print SOT "insert_size","\n";

	foreach my $sample (keys %samplelist){
		my $sample_report_outpath="$var{outpath}/Report/Samples/$sample";
		if ( !-d $sample_report_outpath ) {make_path $sample_report_outpath or die "Failed to create path: $sample_report_outpath";}
		`cp $var{outpath}/$sample/$sample.bam.stats.txt $sample_report_outpath`;
		open IN, "$sample_report_outpath/$sample.bam.stats.txt";

		$sample_summary{$sample}{"covered_regions"} = 0;

		while (<IN>){
			if (/SN\s+bases\s+mapped\s+\(cigar\)\:\s+(\d+)/){
				$sample_summary{$sample}{"mapped_bases"}=$1;
				$sample_summary{$sample}{"mean_covreage"}=$sample_summary{$sample}{"mapped_bases"}/$var{temp_ref}{length};
				$sample_summary{$sample}{"mapping_rate"}=$sample_summary{$sample}{"mapped_bases"}/$sample_summary{$sample}{"clean_data"};
			}

			if (/SN\s+insert\s+size\s+average\:\s+(\d+)/){
				$sample_summary{$sample}{"insert_size"} = $1;
			}

			if (/SN\s+insert\s+size\s+standard\s+deviation\:\s+(\d+)/){
				$sample_summary{$sample}{"insert_size_std"} = $1;
			}

			if(/COV\s+\S+\]\s+(\d+)\s+(\d+)/){
				if($1>10){
					$sample_summary{$sample}{"covered_regions"} += $2;
				}
			}
		}

		close IN;
		print SOT "$sample\t";
		print SOT "$var{temp_ref}{name}\t";
		print SOT $sample_summary{$sample}{"percentage_of_cleandata"},"\t";
		print SOT $sample_summary{$sample}{"clean_data"},"\t";
		print SOT $sample_summary{$sample}{"mean_covreage"},"\t";
		print SOT $sample_summary{$sample}{"mapping_rate"},"\t";
		print SOT $sample_summary{$sample}{"covered_regions"}, "\t";
		print SOT $sample_summary{$sample}{"insert_size"}, "\n";
		#SN      insert size average:    486.4
		#SN      insert size standard deviation: 1454.8

	}
	close SOT;

	open SOT, ">$var{outpath}/Report/sample_mapping_summary_formatted.xls";

	print SOT "Sample ID","\t";
	print SOT "Reference","\t";
	print SOT "Percentage of Cleandata (%)","\t";
	print SOT "Clean_data (Gb)","\t";
	print SOT "Mean Covreage","\t";
	print SOT "Mapping Rate (%)","\t";
	print SOT "Covered Regions","\t";
	print SOT "Insert Size","\n";

	foreach my $sample (keys %samplelist){
		print SOT "$sample\t";
		print SOT "$var{temp_ref}{name}\t";

		print SOT sprintf("%.2f",100*$sample_summary{$sample}{"percentage_of_cleandata"}),"\t";
		print SOT sprintf("%.2f",$sample_summary{$sample}{"clean_data"}/1000000000),"\t";
		print SOT sprintf("%.2f",$sample_summary{$sample}{"mean_covreage"}),"\t";
		print SOT sprintf("%.2f",100*$sample_summary{$sample}{"mapping_rate"}),"\t";
		print SOT $sample_summary{$sample}{"covered_regions"}, "\t";
		print SOT $sample_summary{$sample}{"insert_size"}, "\n";
	}
	close SOT;
}

1;