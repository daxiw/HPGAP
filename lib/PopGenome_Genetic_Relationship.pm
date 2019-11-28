package PopGenome_Genetic_Relationship;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use Getopt::Long qw(GetOptionsFromArray);
use FindBin '$Bin';
use YAML::Tiny;
use lib "$Bin/lib";


sub Main{
	my $args = shift; 
	my @args = @{$args};
	my %opts;
	my %var;

	GetOptionsFromArray (\@args, \%opts, 
		'config=s',
		'overwrite',
		'allsteps',
		'samplelist=s',
		'threads',
		'admixture',
		'phylogeny',
		'pca',
		'help',
		'skipsh');

	die "please provide the correct configuration file" unless ((defined $opts{config}) && (-e $opts{config}));

	if (defined $opts{allsteps}){
		$opts{variant_filtering} = 1;
		$opts{intersection} = 1;
	}

	my $yaml = YAML::Tiny->read( $opts{config} );
	my %cfg = %{$yaml->[0]};	
	$var{cfg}=\%cfg;

	my %samplelist_ori = %{$cfg{fqdata}};
	my %samplelist = %samplelist_ori;
	$var{samplelist}=\%samplelist;

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

	#set the number of threads
	if (defined $opts{threads}){
		$var{threads} = $opts{threads};
	}elsif(defined $cfg{args}{threads}){
		$var{threads} = $cfg{args}{threads};
	}else{
		$var{threads} = 4;
	}

	$var{outpath} = "$cfg{args}{outdir}/GeneticRelationship/";
	if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 

	$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/GeneticRelationship/";
	if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}

	die "please add genome path into configuration file" unless (defined $cfg{ref}{db}{$cfg{ref}{choose}}{path});
	$var{reference} = $cfg{ref}{db}{$cfg{ref}{choose}}{path};
	die "$var{reference} does not exists" unless (-e $var{reference});

	if (defined $opts{admixture}){ &ADMIXTURE (\%var,\%opts);}
	if (defined $opts{phylogeny}){ &PHYLOGENY (\%var,\%opts);}
	if (defined $opts{pca}){ &PCA (\%var,\%opts);}
}

#################################
#			   #
#   	GeneticRelationships_admixture    		#
#			   #
#################################
sub ADMIXTURE{

	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	my $outpath = "$var{outpath}/Admixture/";
	my $plink_data;
	if ( !-d $outpath ) {make_path $outpath or die "Failed to create path: $outpath";}
	if (defined $cfg{variant_filtering}{plink_data}){
		$plink_data = $cfg{variant_filtering}{plink_data};
	}else {
		$plink_data = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$cfg{ref}{choose}/Filtered_variants/nosingle_snp_dp_prunned";
	}

	open CL, ">$var{shpath}/cmd_admixture_s1.list";
	open SH, ">$var{shpath}/admixture.sh";
	print SH "cd $outpath\n";
	print SH "cp -fr $plink_data* ./\n";
	
	for (my $k=0; $k<9;$k++){
		print CL "admixture --cv $plink_data.bed $k 1>$var{shpath}/$k.admixture.o 2>$var{shpath}/$k.admixture.e\n";
	}
	close CL;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/admixture_s1_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/cmd_admixture_s1.list` unless (defined $opts{skipsh});

	# generate the sample list file (for the sample order in the graph)
	`cp -f $Bin/lib/admixture.R $var{shpath}/admixture.R`;
	
	open SL, ">$outpath/sample.list";
	foreach (keys %samplelist){
		print SL "$_\n";
	}
	close SL;

	# generate the arguments for Rscript
	#undef @pp if (defined @pp);
	my @pp = ();
	push @pp, $outpath; 
	push @pp, "$outpath/sample.list";
	push @pp, "K";
	my $p = join (' ',@pp);

	`cat $var{shpath}/*admixture.o|grep \"CV error\" >$outpath/CV.error.txt`;

	open CL, ">$var{shpath}/cmd_admixture_s2.list";
	open SH, ">$var{shpath}/admixture_s2.sh";
	print SH "Rscript --vanilla $var{shpath}/admixture.R $p\n";
	close SH;
	close CL;
	
	`perl $Bin/lib/qsub.pl -d $var{shpath}/admixture_s2_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/cmd_admixture_s2.list` unless (defined $opts{skipsh});
}

1;