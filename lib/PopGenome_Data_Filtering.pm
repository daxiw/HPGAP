package PopGenome_Data_Filtering;
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
#			   			   #
#   	Data filtering     #
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
		'filter',
		'report',
		'help',
		'skipsh=i');

	$opts{config}||="allcfg.yml";
	$opts{skipsh}||= 0;

	if (defined $opts{allsteps}){
		$opts{filter} = 1;
		$opts{report} = 1;
	}

	my $yaml = YAML::Tiny->read( $opts{config} );
	my %cfg = %{$yaml->[0]};
	my %samplelist = %{$cfg{fqdata}};

	$var{outpath} = "$cfg{args}{outdir}/01.QualityControl/read_filtering/"; 
	if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 
	$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_filtering/";
	if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}
	$var{samplelist}=\%samplelist;
	$var{cfg}=\%cfg;

	if (defined $opts{filter}){ &DataFiltering (\%var,\%opts);}

	if (defined $opts{report}){ &ReadReport (\%var,\%opts);}
}

sub DataFiltering{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	open CL, ">$var{shpath}/cmd_read_filtering.list";
	foreach my $sample (keys %samplelist){
		my $sample_outpath="$var{outpath}/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}

		open SH, ">$var{shpath}/$sample.read_filtering.sh";
		print SH "#!/bin/sh\ncd $sample_outpath\n";
		foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
			
			my $read;
			if ($samplelist{$sample}{rawdata}{$lib}{fq1} =~ /gz$/){
				$read = `gunzip -c $samplelist{$sample}{rawdata}{$lib}{fq1}|head -n 2|tail -n 1`;
			}else{
				$read = `cat $samplelist{$sample}{rawdata}{$lib}{fq1}|head -n 2|tail -n 1`;
			}
			my @temp = split //, $read;
			$samplelist{$sample}{rawdata}{$lib}{Length} = @temp;
			$samplelist{$sample}{rawdata}{$lib}{Length} = int($samplelist{$sample}{rawdata}{$lib}{Length}*0.5);
			
			if (-e "$sample_outpath/$lib\_1.filt.fq.gz"){ print SH "#" unless (defined $opts{overwrite});}
			if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "PE"){
				print SH "fastp -i $samplelist{$sample}{rawdata}{$lib}{fq1} -I $samplelist{$sample}{rawdata}{$lib}{fq2} -o $lib\_1.filt.fq.gz -O $lib\_2.filt.fq.gz --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --detect_adapter_for_pe --disable_trim_poly_g -q 20 -u 30 -n 2 --length_required $samplelist{$sample}{rawdata}{$lib}{Length} -w 4 -j $lib.fastp.json -h $lib\_1.fastp.html -R \"$sample $lib fastp report\" && echo \"** finish mt_genome_mapping **\" > $var{shpath}/$sample.$lib.read_filtering.finished.txt\n";
			}
			if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "SE"){
				print SH "fastp -i $samplelist{$sample}{rawdata}{$lib}{fq1} -o $lib\_1.filt.fq.gz --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --detect_adapter_for_pe --disable_trim_poly_g -q 20 -u 30 -n 2 --length_required $samplelist{$sample}{rawdata}{$lib}{Length} -w 4 -j $lib.fastp.json -h $lib\_1.fastp.html -R \"$sample $lib fastp report\" && echo \"** finish mt_genome_mapping **\" > $var{shpath}/$sample.$lib.read_filtering.finished.txt\n";
			}

		}
		close SH;
		print CL "sh $var{shpath}/$sample.step1a.sh 1>$var{shpath}/$sample.step1a.sh.o 2>$var{shpath}/$sample.step1a.sh.e\n";
	}
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_read_filtering_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=2G,num_proc=4 -binding linear:1' -m 100 -r $var{shpath}/cmd_read_filtering.list` unless ($skipsh ==1);

	my $flag_finish = 0;
	my $sample_number = 0;
	while(1){
		sleep(10);
		foreach my $sample (keys %samplelist){
			foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
				$sample_number ++;
				if(-e "$var{shpath}/$sample.$lib.read_filtering.finished.txt"){$flag_finish +=1;}
			}
		}
		my $datestring = localtime();
		print "waiting for read_filtering to be done [$flag_finish out of $sample_number samples are finished] at $datestring\n";
		last if($flag_finish == $sample_number);
	}
}

sub ReadReport{

	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};

	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	my $report_outpath="$var{outpath}/Report"; 
	if ( !-d $report_outpath ) {make_path $report_outpath or die "Failed to create path: $report_outpath";}
	my $report_sample_outpath="$var{outpath}/Report/Samples"; 
	if ( !-d $report_sample_outpath ) {make_path $report_sample_outpath or die "Failed to create path: $report_outpath";}

	open OT, ">$var{outpath}/Report/read_quality_summary.xls";
    
    print OT "sampleID","\t";
    print OT "before_filtering_reads","\t";
    print OT "before_filtering_bases","\t";
    print OT "after_filtering_reads","\t";
    print OT "after_filtering_bases","\t";
    print OT "q20_rate","\t";
    print OT "q30_rate","\t";
    print OT "gc_content","\t";
    print OT "duplication_rate","\t";
    print OT "adapter_trimmed_reads","\n";

	foreach my $sample (keys %samplelist){	
		my $sample_report_outpath="$var{outpath}/Report/Samples/$sample"; 
		if ( !-d $sample_report_outpath ) {make_path $sample_report_outpath or die "Failed to create path: $sample_report_outpath";}

		# copy filtering statistics from read_filtering

		`cp $var{outpath}/$sample/*json $sample_report_outpath`;

		my $json;
		{
		  local $/; #Enable 'slurp' mode
		  open my $fh, "<", "$sample_report_outpath/$sample.fastp.json";
		  $json = <$fh>;
		  close $fh;
		}

		my $data = decode_json($json);

		print OT "$sample\t";
		print OT $data->{'summary'}->{'before_filtering'}->{'total_reads'}, "\t";
		print OT $data->{'summary'}->{'before_filtering'}->{'total_bases'}, "\t";
		print OT $data->{'summary'}->{'after_filtering'}->{'total_reads'}, "\t";
		print OT $data->{'summary'}->{'after_filtering'}->{'total_bases'}, "\t";
		print OT $data->{'summary'}->{'after_filtering'}->{'q20_rate'}, "\t";
		print OT $data->{'summary'}->{'after_filtering'}->{'q30_rate'}, "\t";
		print OT $data->{'summary'}->{'after_filtering'}->{'gc_content'}, "\t";
		print OT $data->{'duplication'}->{'rate'}, "\t";
        print OT $data->{'adapter_cutting'}->{'adapter_trimmed_reads'}, "\n";
	}
	close OT;
}

1;