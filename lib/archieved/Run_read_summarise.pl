use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use FindBin '$Bin';
use YAML::Tiny;
use Getopt::Long;

my ($yml_file,$skipsh);
GetOptions(
    '-i=s' => \$yml_file,
    '-s=s' => \$skipsh,
);

if(!$yml_file){die;}
$skipsh ||= 0;

my $yaml = YAML::Tiny->read( $yml_file );
my %cfg = %{$yaml->[0]};
my %samplelist = %{$cfg{fqdata}};

foreach my $temp_ref(keys %{$cfg{ref}{db}}){
	my $outpath = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$temp_ref"; 
	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";} 
	my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$temp_ref";
	if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}
	
    my $report_outpath="$outpath/Report"; if ( !-d $report_outpath ) {make_path $report_outpath or die "Failed to create path: $report_outpath";}
    my $report_sample_outpath="$outpath/Report/Samples"; if ( !-d $report_sample_outpath ) {make_path $report_sample_outpath or die "Failed to create path: $report_outpath";}

	foreach my $sample (keys %samplelist){
        my $sample_report_outpath="$outpath/Report/Samples/$sample"; if ( !-d $sample_report_outpath ) {make_path $sample_report_outpath or die "Failed to create path: $sample_report_outpath";} 
        
        `cat $output/$sample/bam.stats.txt|grep ^COV | cut -f 2- >$sample_report_outpath/COV.stat.txt`;
        `cat $output/$sample/bam.stats.txt|grep ^IS | cut -f 2- >$sample_report_outpath/IS.stat.txt`;
        `cat $output/$sample/bam.stats.txt|grep ^SN | cut -f 2-|sed "s/ //g;" >$sample_report_outpath/SN.stat.txt`;
        
        open IN, "$output/$sample/bam.stats.txt";
            
        close IN;

        #foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
        #    if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "PE"){
        #	    next;
        #	}
        #	elsif($samplelist{$sample}{rawdata}{$lib}{Flag} eq "SE"){
        #	    next;
        #	}
        #}
	}
}
