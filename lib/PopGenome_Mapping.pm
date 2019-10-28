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
		'threads',
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

	if (defined $opts{reference_selection}){ & SelectReference (\%var,\%opts); last;}

	foreach my $temp_ref(keys %{$cfg{ref}{db}}){
		$var{outpath} = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$temp_ref"; 
		if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 
		$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$temp_ref";
		if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}

		die "please add mt genome path into configuration file" unless (defined $cfg{ref}{db}{$temp_ref}{path});
		$var{reference} = $cfg{ref}{db}{$temp_ref}{path};
		die "$var{reference} does not exists" unless (-e $var{reference});

		if (defined $opts{read_mapping}){ & ReadMapping (\%var,\%opts);}

		if (defined $opts{mapping_report}){ & MappingReport (\%var,\%opts);}
	
		#### estimate phylogeny of mt genomes ###		
	}
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

		open SH, ">$shpath/$sample.readmapping.sh";		
		print SH "#!/bin/sh\ncd $sample_outpath\n";
		foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
			if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "PE"){
				print SH "bwa mem $reference ../../read_filtering/$sample/$lib\_1.filt.fq.gz ../../read_filtering/$sample/$lib\_2.filt.fq.gz -t $var{threads} -R \"\@RG\\tID:$lib\\tSM:$sample\\tLB:$lib\\tPL:$samplelist{$sample}{rawdata}{$lib}{PL}\"\| samtools view -bS -@ $cfg{args}{threads} -F 4 - -o $lib\_filt.bam && \\\n";
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
			print SH "samtools stats -@ $var{threads} $sample.sorted.bam 1>bam.stats.txt 2>bam.stats.txt.e && echo \"** bam.stats.txt done **\"\n";
		}

		#when there is more than one library/lane for each sample
		if (keys %{$samplelist{$sample}{rawdata}} > 1){
			print SH "samtools merge -nr -@ $cfg{args}{threads} $sample.sorted.bam *_filt.sort.bam && echo \"** $sample.sort.bam done **\" && rm -f *_filt.sort.bam\n";
			#print SH "samtools sort -@ $cfg{args}{threads} $sample.bam -o $sample.sorted.bam --output-fmt BAM && echo \"** $sample.sorted.bam done **\" && rm -f $sample.bam\n";
			print SH "samtools stats -@ $var{threads} $sample.sorted.bam 1>bam.stats.txt 2>bam.stats.txt.e && echo \"** bam.stats.txt done **\"\n";
		}
		close SH;
		print CL "sh $shpath/$sample.readmapping.sh 1>$shpath/$sample.readmapping.sh.o 2>$shpath/$sample.readmapping.sh.e \n";
	}
	close CL;
	my $threads = $cfg{args}{threads};
	`perl $Bin/lib/qsub.pl -d $shpath/cmd_readmapping_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=$cfg{args}{mem},num_proc=$var{threads} -binding linear:1' -m 100 -r $shpath/cmd_readmapping.list` unless (defined $opts{skipsh});
}

1;