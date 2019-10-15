package PopGenome_Read_Report;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use FindBin '$Bin';
use YAML::Tiny;
use JSON;
use lib "$Bin/lib";
use PopGenome_Shared;

#################################
#			   #
#    Step 1c Report summary	    #
#			   #
#################################
sub READ_REPORT{
	my ($yml_file,$skipsh) = @_;
	my $yaml = YAML::Tiny->read( $yml_file );
	my %cfg = %{$yaml->[0]};
	my %samplelist = %{$cfg{fqdata}};

	my %rs;
	my %rs_ref;
	foreach my $temp_ref(keys %{$cfg{ref}{db}}){
		my $outpath = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$temp_ref"; 
		if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";} 
		my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$temp_ref";
		if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

		my $report_outpath="$outpath/Report"; if ( !-d $report_outpath ) {make_path $report_outpath or die "Failed to create path: $report_outpath";}

		my $report_sample_outpath="$outpath/Report/Samples"; if ( !-d $report_sample_outpath ) {make_path $report_sample_outpath or die "Failed to create path: $report_outpath";}

		open OT, ">$outpath/Report/read_quality_summary.xls";

		foreach my $sample (keys %samplelist){
			
			my $sample_report_outpath="$outpath/Report/Samples/$sample"; if ( !-d $sample_report_outpath ) {make_path $sample_report_outpath or die "Failed to create path: $sample_report_outpath";}

			# copy filtering statistics from read_filtering

			`cp $cfg{args}{outdir}/01.QualityControl/read_filtering/$sample/*json $sample_report_outpath`;

			my $json;
			{
			  local $/; #Enable 'slurp' mode
			  open my $fh, "<", "$sample.fastp.json";
			  $json = <$fh>;
			  close $fh;
			}

			my $data = decode_json($json);

			print "$sample\t";
			print OT $data->{'summary'}->{'before_filtering'}->{'total_reads'}, "\t";
			print OT $data->{'summary'}->{'before_filtering'}->{'total_bases'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'total_reads'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'total_bases'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'q20_rate'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'q30_rate'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'gc_content'}, "\t";
			print OT $data->{'summary'}->{'after_filtering'}->{'total_reads'}, "\t";

		}
		close OT;

	}
}

1;