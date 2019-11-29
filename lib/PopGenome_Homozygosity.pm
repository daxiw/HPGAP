package PopGenome_Homozygosity;
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
		'homozygosity',
		'roh',
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

	if (defined $opts{homozygosity}){ &HOMOZYGOSITY (\%var,\%opts);}
	if (defined $opts{roh}){ &ROH (\%var,\%opts);}
}

#################################
#			   #
#   	step4_Homozygosity   	#
#			   #
#################################
sub HOMOZYGOSITY{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	my $genome = PopGenome_Shared::LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});
	my $gzvcf = $cfg{variant_filtering}{high_confidence_vcf};

	if ( !-d "$var{outpath}/Homozygosity" ) {make_path "$var{outpath}/Homozygosity" or die "Failed to create path: $var{outpath}/Homozygosity";}

	my %pop;
	foreach my $sample (keys %{$cfg{population}}){
		unless (exists $pop{$cfg{population}{$sample}{'presumed population'}}{line}){
			$pop{$cfg{population}{$sample}{'presumed population'}}{count} = 0;
			$pop{$cfg{population}{$sample}{'presumed population'}}{line} = "";
		}
		$pop{$cfg{population}{$sample}{'presumed population'}}{count}++;
		$pop{$cfg{population}{$sample}{'presumed population'}}{line} .= "$sample\n";
	}

	open SH, ">$shpath/Homozygosity.sh";
	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 4);
		open OT, ">"$var{outpath}/Homozygosity"/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;

		print SH "cd "$var{outpath}/Homozygosity";vcftools --gzvcf $gzvcf --keep "$var{outpath}/Homozygosity"/$pop_name.list --het --stdout ",'|perl -ne \'if(/INDV/){print "INDV\tO_HOM\tE_HOM\tN_SITES\tPERCENTAGE\tF\n";}else{@a=split /\t/;$p=$a[1]/$a[3]*100;$p=sprintf("%.2f", $p);print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$p\t$a[4]";}\'', ">$pop_name.het.list\n";
	}	
	close SH;

	`sh $shpath/Homozygosity.sh 1>$shpath/Homozygosity.sh.o 2>$shpath/Homozygosity.sh.e` unless ($skipsh ==1);

}

sub ROH{
	my ($var,$opts) = @_;
	my %opts = %{$opts};
	my %var = %{$var};
	my %cfg = %{$var{cfg}};
	my %samplelist = %{$var{samplelist}};

	my $gzvcf = $cfg{step1}{variant_filtering}{high_confidence_vcf};
	my $genome = PopGenome_Shared::LOADREF($cfg{ref}{db}{$cfg{ref}{choose}}{path});

	if ( !-d "$var{outpath}/ROH" ) {make_path "$var{outpath}/ROH" or die "Failed to create path: $var{outpath}/ROH";}
	if ( !-d $shpath ) { make_path $shpath or die "Failed to create path: $shpath";}

	$cfg{step4}{ROH}{windowsize} = 100 unless (defined $cfg{step4}{ROH}{windowsize});
	$cfg{step4}{ROH}{scaffold_number_limit} = 95 unless (defined $cfg{step4}{ROH}{scaffold_number_limit});
	$cfg{step4}{ROH}{scaffold_length_cutoff} = 1000000 unless (defined $cfg{step4}{ROH}{scaffold_length_cutoff});
	my $window_size = $cfg{step4}{ROH}{windowsize};
	my $scaffold_number_limit = $cfg{step4}{ROH}{scaffold_number_limit};
	my $scaffold_length_cutoff = $cfg{step4}{ROH}{scaffold_length_cutoff};

	### load reference and map scaffold to number list 
	my $i=1;my $j=0;
	open OT, ">"$var{outpath}/Homozygosity"/chr_map.list";
	foreach my $id(sort { $genome->{len}{$b} <=> $genome->{len}{$a} } keys %{$genome->{len}}){
		print OT "$id\t$i\n";
		$i++;
		if (($genome->{len}{$id} >= $cfg{step4}{ROH}{scaffold_length_cutoff})&&($j < $cfg{step4}{ROH}{scaffold_number_limit})){
			$j++;
		}
	}
	close OT;

	$cfg{step4}{ROH}{scaffold_number_limit} = $j;

	my %pop;
	foreach my $sample (keys %{$cfg{population}}){
		unless (exists $pop{$cfg{population}{$sample}{'presumed population'}}{line}){
			$pop{$cfg{population}{$sample}{'presumed population'}}{count} = 0;
			$pop{$cfg{population}{$sample}{'presumed population'}}{line} = "";
		}
		$pop{$cfg{population}{$sample}{'presumed population'}}{count}++;
		$pop{$cfg{population}{$sample}{'presumed population'}}{line} .= "$sample\n";
	}

	foreach my $pop_name (keys %pop){
		next unless ($pop{$pop_name}{count} > 4);
		open SH, ">$shpath/ROH.sh";
		open OT, ">"$var{outpath}/Homozygosity"/$pop_name.list";
		print OT $pop{$pop_name}{line};
		close OT;
		print SH "cd "$var{outpath}/Homozygosity"\n";
		print SH 'bcftools annotate --threads 10 --rename-chrs', " "$var{outpath}/Homozygosity"/chr_map.list" ," $gzvcf", '|perl -ne \'if (/#\S+ID=(\d+)/){if($1<=',"$scaffold_number_limit",'){print;}}elsif(/^#/){print;}elsif(/^(\d+)\s+/){if($1<= ',"$scaffold_number_limit",'){print;}}\'', ">"$var{outpath}/Homozygosity"/$pop_name.pre.vcf\n"; 
		print SH "vcftools --vcf "$var{outpath}/Homozygosity"/$pop_name.pre.vcf --keep "$var{outpath}/Homozygosity"/$pop_name.list --plink --out $pop_name.plink \n";
		print SH "plink --file $pop_name.plink --make-bed --chr-set $scaffold_number_limit no-xy no-mt no-y --out $pop_name.plink\n";
		print SH "plink --bfile "$var{outpath}/Homozygosity"/$pop_name.plink --chr-set $scaffold_number_limit no-xy no-mt no-y --homozyg --homozyg-window-het 2 --homozyg-kb $window_size --homozyg-gap 50 --homozyg-window-snp 100 --homozyg-window-missing 5 --out $pop_name\n";
		close SH; 
	}
	`sh $shpath/ROH.sh 1>$shpath/ROH.sh.o 2>$shpath/ROH.sh.e` unless ($skipsh ==1);
}

1;