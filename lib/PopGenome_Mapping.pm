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

	if (defined $opts{allsteps}){
		$opts{read_mapping} = 1;
		$opts{mapping_report} = 1;
		$opts{reference_selection} = 1;
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

	$var{mem} = $var{threads}."G";
	if (defined $opts{reference_selection}){ & SelectReference (\%var,\%opts); last;}

	my %mapping;
	foreach my $temp_ref(keys %{$cfg{ref}{db}}){
		$var{outpath} = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$temp_ref"; 
		if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 
		$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$temp_ref";
		if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}

		die "please add mt genome path into configuration file" unless (defined $cfg{ref}{db}{$temp_ref}{path});
		$var{reference} = $cfg{ref}{db}{$temp_ref}{path};
		die "$var{reference} does not exists" unless (-e $var{reference});

		open IN, $var{reference};
		$var{$temp_ref}{length}=0;
		while (<IN>){
			if(/\>/){next;}
			elsif(/\w+/){
				$var{$temp_ref}{length} += length;
			}
		}
		close IN;

		if (defined $opts{read_mapping}){ $mapping{$temp_ref} = & ReadMapping (\%var,\%opts);}

		#### estimate phylogeny of mt genomes ###		
	}

	# wait for the read mapping 
	while(defined $opts{read_mapping}){
		sleep(10);
		my $flag_finish = 1; 

		foreach my $temp_ref(keys %{$cfg{ref}{db}}){
			$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$temp_ref";
			#####

			my %subsamplelist = %{ $mapping{$temp_ref} };

			foreach my $sample (keys %samplelist){
				next if ($subsamplelist{$sample}{finish_flag} eq "finished");
				if(-e "$var{shpath}/$sample.readmapping.finished.txt"){
					next;
				}else{
					$flag_finish = 0;
				}
			}
			#####
		}
		my $datestring = localtime();
		print "waiting for readmapping to be done at $datestring\n";
		last if($flag_finish == 1);
	}

	foreach my $temp_ref(keys %{$cfg{ref}{db}}){
		$var{outpath} = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$temp_ref"; 
		if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 
		$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$temp_ref";
		if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}

		if (defined $opts{mapping_report}){ & MappingReport (\%var,\%opts);}
	}

	#$var{outfig} = $opts{config};
	#$var{outfig} =~ s/\.yml|\.yaml/_mapping_done\.yml/g;
	#$opts{outcfg} ||= $var{outfig};
	# create this yaml object
    #$yaml = YAML::Tiny->new( $var{newcfg} );
    # Save both documents to a file
    #$yaml->write( $opts{outcfg} );
#    print "$opts{outpath}\n";
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
			my $required_coverage = $opts{downsize}*$var{$temp_ref}{length};
			print SH "samtools stats -@ $var{threads} $sample.sorted.bam 1>$sample.sorted.bam.stats.txt 2>$sample.sorted.bam.stats.error\n";
			print SH "a=\`cat $sample.sorted.bam.stats.txt|perl -ne \'if(/cigar\\):\\s+(\\d+)/){\$b=$required_coverage/\$1;if(\$b<1){print \$b}else{print 1}}\'\`\n";
			print SH "mv $sample.sorted.bam $sample.sorted.ori.bam";
			print SH "samtools view -s \$a $sample.sorted.ori.bam -o $sample.sorted.bam\n";
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

1;