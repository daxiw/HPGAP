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

my $outpath = "$cfg{args}{outdir}/01.QualityControl/centrifuge"; 
if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";} 
my $shpath = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/centrifuge";
if ( !-d $shpath ) {make_path $shpath or die "Failed to create path: $shpath";}

open CL, ">$shpath/centrifuge.list";
foreach my $sample (keys %samplelist){
	my $sample_outpath="$outpath/$sample"; if ( !-d $sample_outpath ) {make_path $sample_outpath or die "Failed to create path: $sample_outpath";}

	open SH, ">$shpath/$sample.sh";		
	print SH "#!/bin/sh\ncd $sample_outpath\n";
	foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
		if($samplelist{$sample}{rawdata}{$lib}{Flag} eq "PE"){
			print SH "zcat ../../read_filtering/$sample/$lib\_1.filt.fq.gz | head -n 1000000 >$lib\_head1M_1.fq && \\\n";
			print SH "zcat ../../read_filtering/$sample/$lib\_2.filt.fq.gz | head -n 1000000 >$lib\_head1M_2.fq \n";
		}
		elsif($samplelist{$sample}{rawdata}{$lib}{Flag} eq "SE"){
			print SH "zcat ../../read_filtering/$sample/$lib\_1.filt.fq.gz | head -n 1000000 >$lib\_head1M_1.fq \n";
		}

	}

	foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
	#when there is only one library/lane for each sample
		if ($samplelist{$sample}{rawdata}{$lib}{Flag} eq "PE"){
			#foreach my $lib (keys %{$samplelist{$sample}{rawdata}}){
			print SH "/ldfssz1/ST_INFECTION/P17Z10200N0536_Echinococcus/USER/wangdaxi/Programs/centrifuge/centrifuge-1.0.3-beta/centrifuge -x /ldfssz1/ST_INFECTION/P17Z10200N0536_Echinococcus/USER/wangdaxi/Programs/centrifuge/index/p_compressed -1 $lib\_head1M_1.fq -2 $lib\_head1M_2.fq --report-file $lib\_centrifuge.report -S $lib\_centrifuge.classification\n";
		}

		#when there is more than one library/lane for each sample
		elsif ($samplelist{$sample}{rawdata}{$lib}{Flag} eq "SE"){
			print SH "/ldfssz1/ST_INFECTION/P17Z10200N0536_Echinococcus/USER/wangdaxi/Programs/centrifuge/centrifuge-1.0.3-beta/centrifuge -x /ldfssz1/ST_INFECTION/P17Z10200N0536_Echinococcus/USER/wangdaxi/Programs/centrifuge/index/p_compressed -1 $lib\_head1M_1.fq --report-file $lib\_centrifuge.report -S $lib\_centrifuge.classification\n";
		}
	}

	close SH;
	print CL "sh $shpath/$sample.sh 1>$shpath/$sample.sh.o 2>$shpath/$sample.sh.e \n";
}
close CL;

my $threads = $cfg{args}{threads};
`perl $Bin/qsub.pl -d $shpath/centrifuge_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=$cfg{args}{mem},num_proc=$threads -binding linear:1' -m 100 -r $shpath/centrifuge.list` unless ($skipsh ==1);

