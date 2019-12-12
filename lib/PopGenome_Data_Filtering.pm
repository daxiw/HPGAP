package PopGenome_Data_Filtering;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use Getopt::Long qw(GetOptionsFromArray);
use FindBin '$Bin';
use YAML::Tiny;
use JSON;
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

	GetOptionsFromArray (\@args, \%opts, 
		'config=s',
		'overwrite',
		'allsteps',
		'threads=s',
		'filter',
		'report',
		'updateconfig',
		'help',
		'skipsh');

	die "please provide the correct configuration file" unless ((defined $opts{config}) && (-e $opts{config}));

	if (defined $opts{allsteps}){
		$opts{filter} = 1;
		$opts{report} = 1;
	}

	my %var = %{PopGenome_Shared::CombineCfg("$Bin/lib/parameter.yml",\%opts, "01.QualityControl/read_filtering")};

	if (defined $opts{filter}){ 
		&DataFiltering (\%var,\%opts);
		&WriteCfg (\%var,\%opts);
	}

	if (defined $opts{report}){ &ReadReport (\%var,\%opts);}
	if (defined $opts{updateconfig}){ &WriteCfg (\%var,\%opts);}
}

sub DataFiltering{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist_ori = %{$var{samplelist}};
	my %samplelist = %samplelist_ori;

	open CL, ">$var{shpath}/cmd_read_filtering.list";
	foreach my $sample (keys %samplelist){
		my $sample_outpath="$var{outpath}/$sample"; 
		if ( !-d $sample_outpath ) {
			make_path $sample_outpath or die "Failed to create path: $sample_outpath";
		}

		### skip finished samples 
		$samplelist{$sample}{finish_flag}="finished";
		foreach my $readgroup (keys %{$samplelist{$sample}{rawdata}}){
			if ((-e "$sample_outpath/$readgroup\_1.filt.fq.gz") && (!defined $opts{overwrite})){
				$samplelist{$sample}{rawdata}{$readgroup}{finish_flag}="finished";
			}
			else{
				$samplelist{$sample}{finish_flag}="unfinished";
				$samplelist{$sample}{rawdata}{$readgroup}{finish_flag}="unfinished";
			}
		}
		next if($samplelist{$sample}{finish_flag} eq "finished");

		open SH, ">$var{shpath}/$sample.read_filtering.sh";
		print SH "#!/bin/sh\ncd $sample_outpath\n";

		foreach my $readgroup (keys %{$samplelist{$sample}{rawdata}}){
			
			next if ($samplelist{$sample}{rawdata}{$readgroup}{finish_flag} eq "finished");

			#estimate read length
			my $read;
			if ($samplelist{$sample}{rawdata}{$readgroup}{fq1} =~ /gz$/){
				$read = `gunzip -c $samplelist{$sample}{rawdata}{$readgroup}{fq1}|head -n 2|tail -n 1`;
			}else{
				$read = `cat $samplelist{$sample}{rawdata}{$readgroup}{fq1}|head -n 2|tail -n 1`;
			}
			my @temp = split //, $read;
			$samplelist{$sample}{rawdata}{$readgroup}{Length} = @temp;
			$samplelist{$sample}{rawdata}{$readgroup}{Length} = int($samplelist{$sample}{rawdata}{$readgroup}{Length}*0.7);
			
			if($samplelist{$sample}{rawdata}{$readgroup}{Flag} eq "PE"){
				print SH "fastp -i $samplelist{$sample}{rawdata}{$readgroup}{fq1} -I $samplelist{$sample}{rawdata}{$readgroup}{fq2} -o $readgroup\_1.filt.fq.gz -O $readgroup\_2.filt.fq.gz --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --detect_adapter_for_pe --disable_trim_poly_g -q 20 -u 30 -n 2 --length_required $samplelist{$sample}{rawdata}{$readgroup}{Length} -w $var{threads} -j $readgroup.fastp.json -h $readgroup\_1.fastp.html -R \"$sample $readgroup fastp report\" && echo \"** finish mt_genome_mapping **\" > $var{shpath}/$sample.$readgroup.read_filtering.finished.txt\n";
			}
			if($samplelist{$sample}{rawdata}{$readgroup}{Flag} eq "SE"){
				print SH "fastp -i $samplelist{$sample}{rawdata}{$readgroup}{fq1} -o $readgroup\_1.filt.fq.gz --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --detect_adapter_for_pe --disable_trim_poly_g -q 20 -u 30 -n 2 --length_required $samplelist{$sample}{rawdata}{$readgroup}{Length} -w $var{threads} -j $readgroup.fastp.json -h $readgroup\_1.fastp.html -R \"$sample $readgroup fastp report\" && echo \"** finish mt_genome_mapping **\" > $var{shpath}/$sample.$readgroup.read_filtering.finished.txt\n";
			}
		}
		close SH;
		print CL "sh $var{shpath}/$sample.read_filtering.sh 1>$var{shpath}/$sample.read_filtering.sh.o 2>$var{shpath}/$sample.read_filtering.sh.e\n";
	}
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_read_filtering_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=2G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/cmd_read_filtering.list` unless (defined $opts{skipsh});
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

	open OT, ">$var{outpath}/Report/read_group_quality_summary.xls";
    
    print OT "sampleID","\t";
    print OT "readgroup","\t";
    print OT "before_filtering_reads","\t";
    print OT "before_filtering_bases","\t";
    print OT "after_filtering_reads","\t";
    print OT "after_filtering_bases","\t";
    print OT "q20_rate","\t";
    print OT "q30_rate","\t";
    print OT "gc_content","\t";
    print OT "duplication_rate","\t";
    print OT "adapter_trimmed_reads","\n";

    open SOT, ">$var{outpath}/Report/sample_quality_summary.xls";
    
    print SOT "sampleID","\t";
    print SOT "before_filtering_reads","\t";
    print SOT "before_filtering_bases","\t";
    print SOT "after_filtering_reads","\t";
    print SOT "after_filtering_bases","\n";

	foreach my $sample (keys %samplelist){
		print "$sample\n";
		my $sample_report_outpath="$var{outpath}/Report/Samples/$sample";
		print "$sample_report_outpath/$sample.fastp.json\n";
		if ( !-d $sample_report_outpath ) {make_path $sample_report_outpath or die "Failed to create path: $sample_report_outpath";}

		# copy filtering statistics from read_filtering

		`cp $var{outpath}/$sample/*json $sample_report_outpath`;

		my $before_filtering_reads = 0;
		my $before_filtering_bases = 0;
		my $after_filtering_reads = 0;
		my $after_filtering_bases = 0;

		foreach my $readgroup (keys %{$samplelist{$sample}{rawdata}}){
			my $json;
			{
			  local $/; #Enable 'slurp' mode
			  open JS, "$sample_report_outpath/$readgroup.fastp.json";
			  $json = <JS>;
			  close JS;
			}

			my $data = JSON::decode_json($json);

			print OT "$sample\t";
			print OT "$readgroup\t";
			print OT $data->{'summary'}->{'before_filtering'}->{'total_reads'}, "\t";
			print OT $data->{'summary'}->{'before_filtering'}->{'total_bases'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'total_reads'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'total_bases'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'q20_rate'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'q30_rate'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'gc_content'}, "\t";
			print OT $data->{'duplication'}->{'rate'}, "\t";
	        print OT $data->{'adapter_cutting'}->{'adapter_trimmed_reads'}, "\n";

	        $before_filtering_reads += $data->{'summary'}->{'before_filtering'}->{'total_reads'};
	        $before_filtering_bases += $data->{'summary'}->{'before_filtering'}->{'total_bases'};
	        $after_filtering_reads += $data->{'summary'}->{'after_filtering'}->{'total_reads'};
	        $after_filtering_bases += $data->{'summary'}->{'after_filtering'}->{'total_bases'};
    	}
    	print SOT "$sample\t";
    	print SOT "$before_filtering_reads\t";
    	print SOT "$before_filtering_bases\t";
    	print SOT "$after_filtering_reads\t";
    	print SOT "$after_filtering_bases\n";
	}
	close OT;
	close SOT;
}

sub WriteCfg{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	foreach my $sample (keys %samplelist){
		my $sample_outpath="$var{outpath}/$sample"; 
		if ( !-d $sample_outpath ) {
			make_path $sample_outpath or die "Failed to create path: $sample_outpath";
		}

		foreach my $readgroup (keys %{$samplelist{$sample}{rawdata}}){
			if($samplelist{$sample}{rawdata}{$readgroup}{Flag} eq "PE"){
				if (-e "$sample_outpath/$readgroup\_1.filt.fq.gz"){
					$cfg{fqdata}{$sample}{cleandata}{$readgroup}{fq1}="$sample_outpath/$readgroup\_1.filt.fq.gz";
				}
				if (-e "$sample_outpath/$readgroup\_2.filt.fq.gz"){
					$cfg{fqdata}{$sample}{cleandata}{$readgroup}{fq2}="$sample_outpath/$readgroup\_2.filt.fq.gz";
					$cfg{fqdata}{$sample}{cleandata}{$readgroup}{'Flag'} = $cfg{fqdata}{$sample}{rawdata}{$readgroup}{'Flag'};
        			$cfg{fqdata}{$sample}{cleandata}{$readgroup}{'PL'} = $cfg{fqdata}{$sample}{rawdata}{$readgroup}{'PL'};
       				$cfg{fqdata}{$sample}{cleandata}{$readgroup}{'Phred'} = $cfg{fqdata}{$sample}{rawdata}{$readgroup}{'Phred'};
				}
			}
			if($samplelist{$sample}{rawdata}{$readgroup}{Flag} eq "SE"){
				#$readgroup\_1.filt.fq.gz
				if (-e "$sample_outpath/$readgroup\_1.filt.fq.gz"){
					$cfg{fqdata}{$sample}{cleandata}{$readgroup}{fq1}="$sample_outpath/$readgroup\_1.filt.fq.gz";
					$cfg{fqdata}{$sample}{cleandata}{$readgroup}{'Flag'} = $cfg{fqdata}{$sample}{rawdata}{$readgroup}{'Flag'};
        			$cfg{fqdata}{$sample}{cleandata}{$readgroup}{'PL'} = $cfg{fqdata}{$sample}{rawdata}{$readgroup}{'PL'};
       				$cfg{fqdata}{$sample}{cleandata}{$readgroup}{'Phred'} = $cfg{fqdata}{$sample}{rawdata}{$readgroup}{'Phred'};
				}
			}
		}
	}

    my $yaml = YAML::Tiny->new( \%cfg );
    # Save both documents to a file
    $yaml->write( "$cfg{args}{outdir}/.db.yml" );
	# print "$opts{outpath}\n";
}

1;