package PopGenome_Shared;
use File::Basename;
use File::Path qw(make_path);
use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/lib";
use YAML::Tiny;

#----------------------------------- other subroutines ---------------------------------------------
sub LOADREF{
	my ($reference) = @_;
	my %genome;
	open REF, "$reference";
	my $id;
	while (<REF>){
		chomp;
		if (/\>(\S+)/){
			$id = $1;
			$genome{seq}{$id}="";
		}
		elsif(/(\w+)/){
			$genome{seq}{$id}.=$1;
		}
	}
	close REF;
	$genome{sumlen}=0;
	foreach my $id(keys %{$genome{seq}}){
		$genome{len}{$id}=length($genome{seq}{$id});
		$genome{sumlen} += $genome{len}{$id};
	}

	print "LOADREF finished\n";
	return \%genome;
}

sub CombineCfg{
	my ($par,$opts,$folder) = @_;
	my %par = %{YAML::Tiny->read( $par )->[0]};
	my %opts = %{$opts};
	my %cfg = %par;
	my %var;

	die "please provide the correct configuration file" unless ((defined $opts{config}) && (-e $opts{config}));
	my %user_cfg = %{YAML::Tiny->read( $opts{config} )->[0]};
	my %db_cfg = %{YAML::Tiny->read( "$user_cfg{args}{outdir}/.db.yml" )->[0]} if (-e "$user_cfg{args}{outdir}/.db.yml");

	if (defined $user_cfg{'fqdata'}){
		foreach my $sample(keys %{$user_cfg{'fqdata'}}){
			if (defined $user_cfg{'fqdata'}{$sample}{rawdata}){
				$cfg{'fqdata'}{$sample}{rawdata} = $user_cfg{'fqdata'}{$sample}{rawdata};
			}if (defined $user_cfg{'fqdata'}{$sample}{cleandata}){
				$cfg{'fqdata'}{$sample}{cleandata} = $user_cfg{'fqdata'}{$sample}{cleandata};
			}elsif (defined $db_cfg{'fqdata'}{$sample}{cleandata}){
				$cfg{'fqdata'}{$sample}{cleandata} = $db_cfg{'fqdata'}{$sample}{cleandata};
			}
		}
	}

	if (defined $user_cfg{'population'}){
		$cfg{'population'} = $user_cfg{'population'};
	}elsif(defined $db_cfg{'population'}){
		$cfg{'population'} = $db_cfg{'population'};
	}

	if (defined $user_cfg{'ref'}{'choose'}){
		$cfg{'ref'}{'choose'} = $user_cfg{'ref'}{'choose'};
	}elsif(defined $db_cfg{'ref'}{'choose'}){
		$cfg{'ref'}{'choose'} = $db_cfg{'ref'}{'choose'};
	}

	foreach my $key_L1(keys %par){
		foreach my $key_L2(keys %{$par{$key_L1}}){
			if (defined $user_cfg{$key_L1}{$key_L2}){
				$cfg{$key_L1}{$key_L2} = $user_cfg{$key_L1}{$key_L2};
			}
		}
	}

	$var{cfg} = \%cfg;

	if (defined $opts{threads}){
		$var{threads} = $opts{threads};
	}elsif(defined $cfg{args}{threads}){
		$var{threads} = $cfg{args}{threads};
	}else{
		$var{threads} = 1;
	}

	# select samples if sample list is provided
	my %samplelist_ori = %{$cfg{fqdata}};
	my %samplelist = %samplelist_ori;
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
	$var{samplelist}=\%samplelist;

	# use defined vcf if provided
	if (defined $opts{vcf}){
		$var{vcf} = $opts{vcf};
	}elsif (defined $cfg{variant_filtering}{high_confidence_vcf}){
		$var{vcf} = $cfg{variant_filtering}{high_confidence_vcf};
	}

	if (defined $cfg{population}){
		my %pop;
		foreach my $sample (keys %{$cfg{population}}){
			next unless (defined $var{samplelist});
			unless (defined $pop{$cfg{population}{$sample}{'presumed_population'}}{line}){
				$pop{$cfg{population}{$sample}{'presumed_population'}}{count} = 0;
				$pop{$cfg{population}{$sample}{'presumed_population'}}{line} = "";
			}
			$pop{$cfg{population}{$sample}{'presumed_population'}}{count}++;
			$pop{$cfg{population}{$sample}{'presumed_population'}}{line} .= "$sample\n";
		}
		$var{pop}=\%pop;
	}

	if (defined $opts{ploidy}){
		$var{ploidy} = $opts{ploidy};
	}elsif (defined $cfg{args}{ploidy}){
		$var{ploidy} = $cfg{args}{ploidy};
	}

	die "fill in folder" unless (defined $folder);
	if  ($folder eq "Filtered_variants"){
		$var{outpath} = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$cfg{ref}{choose}/$folder/";
		if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 

		$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$cfg{ref}{choose}/$folder/";
		if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}

	}elsif  ($folder eq "VariantCalling"){
		$var{outpath} = "$cfg{args}{outdir}/01.QualityControl/read_mapping.$cfg{ref}{choose}/";
		if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 

		$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/01.QualityControl/read_mapping.$cfg{ref}{choose}/";
		if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}

	}elsif ($folder ne "NULL"){

		$var{outpath} = "$cfg{args}{outdir}/$folder/";
		print '$var{outpath}',"\t",$var{outpath},"\n"; #test_flag
		if ( !-d $var{outpath} ) {make_path $var{outpath} or die "Failed to create path: $var{outpath}";} 

		$var{shpath} = "$cfg{args}{outdir}/PipelineScripts/$folder/";
		if ( !-d $var{shpath} ) {make_path $var{shpath} or die "Failed to create path: $var{shpath}";}

#		die "please add genome path into configuration file" unless (defined $cfg{ref}{db}{$cfg{ref}{choose}}{path});
	}

	if (defined $cfg{ref}{db}{$cfg{ref}{choose}}{path}){
		$var{reference} = $cfg{ref}{db}{$cfg{ref}{choose}}{path};
		die "$var{reference} does not exists" unless (-e $var{reference});
	}

	# create this yaml object
    my $yaml = YAML::Tiny->new( \%cfg );
    # Save both documents to a file
    $yaml->write( "$cfg{args}{outdir}/.combined.yml" );
	# print "$opts{outpath}\n";

	return \%var;	
}

1;