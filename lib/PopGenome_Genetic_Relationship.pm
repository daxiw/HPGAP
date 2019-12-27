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

	%var = %{PopGenome_Shared::CombineCfg("$Bin/lib/parameter.yml",\%opts,"GeneticRelationship")};

	if (defined $opts{admixture}){ &ADMIXTURE (\%var,\%opts);}
	if (defined $opts{phylogeny}){ &PHYLOGENY (\%var,\%opts);}
	if (defined $opts{pca}){ &PCA (\%var,\%opts);}
}

sub ADMIXTURE{

	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	my $admixture_outpath = "$var{outpath}/Admixture/";
	if ( !-d "$var{outpath}/Admixture" ) {make_path "$var{outpath}/Admixture" or die "Failed to create path: $var{outpath}/Admixture";}
	
	my $plink_data;
	if ( !-d $admixture_outpath ) {make_path $admixture_outpath or die "Failed to create path: $admixture_outpath";}
	if (defined $cfg{variant_filtering}{plink_data}){
		$plink_data = $cfg{variant_filtering}{plink_data};
	}else {
		$plink_data = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$cfg{ref}{choose}/Filtered_variants/nosingle_snp_dp_prunned";
	}

	open CL, ">$var{shpath}/cmd_admixture_s1.list";
	open SH, ">$var{shpath}/admixture.sh";
	print SH "#!/bin/sh\ncd $var{outpath}/Admixture/\n";
	print SH "cp -fr $plink_data* ./\n";
	
	for (my $k=0; $k<9;$k++){
		print CL "admixture --cv $plink_data.bed $k 1>$var{shpath}/$k.admixture.o 2>$var{shpath}/$k.admixture.e\n";
	}
	close CL;
	`perl $Bin/lib/qsub.pl -d $var{shpath}/admixture_s1_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/cmd_admixture_s1.list` unless (defined $opts{skipsh});

	# generate the sample list file (for the sample order in the graph)
	`cp -f $Bin/lib/admixture.R $var{shpath}/admixture.R`;
	
	open SL, ">$admixture_outpath/sample.list";
	foreach (keys %samplelist){
		print SL "$_\n";
	}
	close SL;

	# generate the arguments for Rscript
	#undef @pp if (defined @pp);
	my @pp = ();
	push @pp, $admixture_outpath; 
	push @pp, "$admixture_outpath/sample.list";
	push @pp, "K";
	my $p = join (' ',@pp);

	`cat $var{shpath}/*admixture.o|grep \"CV error\" >$admixture_outpath/CV.error.txt`;

	open CL, ">$var{shpath}/cmd_admixture_s2.list";
	open SH, ">$var{shpath}/admixture_s2.sh";
	print SH "Rscript --vanilla $var{shpath}/admixture.R $p\n";
	close SH;
	close CL;
	
	`perl $Bin/lib/qsub.pl -d $var{shpath}/admixture_s2_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=1G,num_proc=2 -binding linear:1' -m 100 -r $var{shpath}/cmd_admixture_s2.list` unless (defined $opts{skipsh});
}

sub PHYLOGENY{

	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	if ( !-d "$var{outpath}/Phylogeny" ) {make_path "$var{outpath}/Phylogeny" or die "Failed to create path: $var{outpath}/Phylogeny";}

	my $i = 0;my %name;
	open OT, ">$var{outpath}/Phylogeny/name_map.list";
	foreach my $sample (keys %samplelist){
		print OT "$sample ", "NS$i", "E\n";
		my $newid = "NS$i"."E";
		$name{$newid} = $sample;
		$i++;
	}
	close OT;

	open CL, ">$var{shpath}/cmd_Phylogeny.list";

	open SH, ">$var{shpath}/Phylogeny.sh";
	print SH "#!/bin/sh\ncd $var{outpath}/Phylogeny\n";

	print SH "bcftools reheader --samples $var{outpath}/Phylogeny/name_map.list -o $var{outpath}/FreebayesCalling/freebayes_joint_calling_rename.vcf $var{outpath}/FreebayesCalling/freebayes_joint_calling.vcf\n";
	print SH "rm -rf $var{outpath}/Phylogeny/RAxML_*\n";

	print SH "cat $var{outpath}/FreebayesCalling/freebayes_joint_calling_rename.vcf|$Bin/Tools/vcf-to-tab >$var{outpath}/Phylogeny/popgenome.tab\n";
	print SH "$Bin/Tools/vcf_tab_to_fasta_alignment.pl -i $var{outpath}/Phylogeny/popgenome.tab > $var{outpath}/Phylogeny/popgenome.fasta\n";

	print SH "$Bin/Tools/fasta-to-phylip --input-fasta $var{outpath}/Phylogeny/popgenome.fasta --output-phy $var{outpath}/Phylogeny/popgenome.phy\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -s $var{outpath}/Phylogeny/popgenome.phy -n trees -T $var{threads} -# 20 -p 12345\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -s $var{outpath}/Phylogeny/popgenome.phy -n boots -T $var{threads} -# 100 -p 23456 -b 23456\n";
	print SH "raxmlHPC-PTHREADS -m GTRGAMMA -p 12345 -f b -t RAxML_bestTree.trees -T $var{threads} -z RAxML_bootstrap.boots -n consensus\n";
	print SH "sumtrees.py --percentages --min-clade-freq=0.50 --target=RAxML_bestTree.trees --output=result2.tre RAxML_bootstrap.boots";

	close SH;
	print CL "sh $var{shpath}/Phylogeny.sh 1>$var{shpath}/Phylogeny.sh.o 2>$var{shpath}/Phylogeny.sh.e \n";
	close CL;

	`perl $Bin/lib/qsub.pl -d $var{shpath}/cmd_Phylogeny_qsub -q $cfg{args}{queue} -P $cfg{args}{prj} -l 'vf=4G,num_proc=$var{threads} -binding linear:1' -m 100 -r $var{shpath}/cmd_Phylogeny.list` unless (defined $opts{skipsh});

	open my $fh, '<', "$var{outpath}/Phylogeny/result2.tre" or die "error opening $var{outpath}/Phylogeny/result2.tre: $!";
	my $data = do { local $/; <$fh> };

	foreach my $newid (keys %name){
		$data =~ s/($newid)/($name{$newid})/g;
	}
	close $fh;

	open OT, ">$var{outpath}/Phylogeny/result2.final.tre";
	print OT $data;
	close OT;
}

1;