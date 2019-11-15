#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Cwd;

my $pwd = getcwd();

my($qsub_opt,@allJobs,$qsubDir,$shell);

use vars qw($opt_d $opt_l $opt_q $opt_N $opt_P $opt_n $opt_b $opt_m $opt_s $opt_r $opt_h);
getopts("d:l:q:N:P:n:b:m:s:rh");

#print "d:l:q:N:n:m:s:S:h\n";
#print "$opt_d\t$opt_l\t$opt_q\t$opt_N\t$opt_n\t$opt_m\t$opt_s\t$opt_S\t$opt_h\n";

if($opt_h or @ARGV == 0){
    &usage();
    exit;
}

$shell = shift;
my $shell_name = (split /\//,$shell)[-1];
$qsubDir = $opt_d || (split /\//,$shell)[-1]."_qsub";
`rm -rf $qsubDir` if(-e $qsubDir);
`mkdir $qsubDir`;
`rm $shell.log` if(-e "$shell.log");
`rm $shell.error` if(-e "$shell.error");
`rm $shell.finished` if(-e "$shell.finished");

$opt_l = $opt_l || "vf=1G";
$qsub_opt = "qsub -cwd -S /bin/bash -l $opt_l ";
$qsub_opt .= "-q $opt_q " if($opt_q);
$qsub_opt .= "-P $opt_P " if($opt_P);
$qsub_opt .= "-l h=$opt_n" if($opt_n);
$opt_N = $opt_N || "work";

my $lines = $opt_b || 1;
my $maxJob = $opt_m || 30;
my $sleepTime = $opt_s || 120;
my $max_try = 3;
$max_try = 1 if(!$opt_r);


my $split_number;
open IS,$shell or die "can\'t open shell.sh: $shell\n";
while(<IS>){
    chomp;
    $split_number++;
    my $num = 1;
    open OUTS,">$qsubDir/$opt_N\_$split_number.sh" or die "can\'t open split shell: $qsubDir/$opt_N\_$split_number.sh\n";
    print OUTS $_;
    while($num < $lines){
        $num++;
        last if(eof(IS));
        chomp(my $command = <IS>);
        print OUTS "\n$command";
    }
    print OUTS " && echo this-work-is-complete\n";
    close OUTS;
    push @allJobs,"$qsubDir/$opt_N\_$split_number.sh";
}
close IS;
my $suspend_flag = 1;

&qsub_and_wait();

sub qsub_and_wait{
    chomp(my $user = `whoami`);

    my(%runJob,%error,@wait);
    my $sub_num = $maxJob > $split_number ? $split_number : $maxJob;
    @wait = (1..$split_number);

    my $qnum = 0;
    while(@wait and $qnum < $sub_num){
        my $i = shift @wait;
        print "$qsub_opt -o $qsubDir/$opt_N\_$i.sh.o -e $qsubDir/$opt_N\_$i.sh.e -N $opt_N\_$i\_$shell_name $qsubDir/$opt_N\_$i.sh\n";
        chomp(my $qmess = `$qsub_opt -o $qsubDir/$opt_N\_$i.sh.o -e $qsubDir/$opt_N\_$i.sh.e -N $opt_N\_$i\_$shell_name $qsubDir/$opt_N\_$i.sh`);
        if($qmess =~ /^[Yy]our\sjob\s(\d+)\s\(\".*\"\)\shas\sbeen\ssubmitted.?$/){
            $runJob{$1} = "$qsubDir/$opt_N\_$i.sh";
            $qnum++;
        }else{
            unshift @wait,$i;
        }
    }

    while(@wait or keys %runJob){
        sleep($sleepTime);
        &check_job($user,\%error,\@wait,\%runJob);
        $qnum = keys %runJob;
#        print "run job numbers: $qnum\nunqsub jobs: ",join "\t",@wait,"\n";
        while(@wait and $qnum < $sub_num){
            my $i = shift @wait;
            print "$qsub_opt -o $qsubDir/$opt_N\_$i.sh.o -e $qsubDir/$opt_N\_$i.sh.e -N $opt_N\_$i\_$shell_name $qsubDir/$opt_N\_$i.sh\n";
            chomp(my $qmess = `$qsub_opt -o $qsubDir/$opt_N\_$i.sh.o -e $qsubDir/$opt_N\_$i.sh.e -N $opt_N\_$i\_$shell_name $qsubDir/$opt_N\_$i.sh`);
            if($qmess =~ /^[Yy]our\sjob\s(\d+)\s\(\".*\"\)\shas\sbeen\ssubmitted.?$/){
                $runJob{$1} = "$qsubDir/$opt_N\_$i.sh";
                $qnum++;
            }else{
                unshift @wait,$i;
            }
        }
    }

    open OUTL,">>$shell.log" or die "can\'t open shell.log\n";
    if(keys %error){
        print OUTL "There are some job can't run finish, check the shell and qsub again\n";
        for(sort {$a cmp $b} keys %error){
            print OUTL "$_\n";
        }
    }else{
        print OUTL "All jobs are finished correctly\n";
    }
    close OUTL;
}

sub check_job{
    my($userName,$error,$wait,$run) = @_;
    my %dead;
    &dead_nodes(\%dead);
    my %running;
    my $qsub_stat = `qstat -xml -u $userName`;
    while($qsub_stat =~ /<JB_job_number>(\d+?)<\/JB_job_number>.*?
            <JB_name>(.+?)<\/JB_name>.*?
            <state>(.+?)<\/state>.*?
            <queue_name>(.*?)<\/queue_name>
            /gxs){
        my ($jbnum, $jbname, $jbstat, $jbqueue) = ($1, $2, $3, $4);
        if($jbname =~ /$opt_N\_(\d+)/){
            my $num = $1;
            my $split_shell = $$run{$jbnum};
            if($jbstat eq "Eqw" or $jbstat eq "T" or ($jbqueue =~ /^.+@(.+)\.local$/ and exists $dead{$1})){
                $$error{$split_shell}++;
                `qdel $jbnum`;
                `echo $split_shell has not finished! >>$shell.error`;
                if($$error{$split_shell} < $max_try){
                    `rm $split_shell.[oe]`;
                    unshift @$wait,$num;
                    `echo $split_shell has been reqsub >>$shell.error`;
                }
                delete $$run{$jbnum};
            }
            $running{$jbnum} = undef;
        }
    }

    foreach my $id (sort {$a <=> $b} keys %$run){
        my $split_shell = $$run{$id};
        if(!exists $running{$id}){
            delete $$run{$id};
            chomp(my $log = `tail -1 $split_shell.o`);
            if($log eq "this-work-is-complete"){
                delete($$error{$split_shell});
                `echo $split_shell is finished! >> $shell.finished`;
            }else{
                `echo $split_shell has not finished! >>$shell.error`;
                $$error{$split_shell}++;
                if($$error{$split_shell} < $max_try){
                    `rm $split_shell.[oe]`;
                    my $num = $1 if($split_shell =~ /$opt_N\_(\d+)\.sh/);
                    unshift @$wait,$num;
                    `echo $split_shell has been reqsub >>$shell.error`;
                }
            }
        }
    }
}

sub waitJobs{
    my $userName = shift;
    my $qsub_stat = `qstat -xml -u $userName`;
    while($qsub_stat =~ /<JB_name>(.+?)<\/JB_name>.*?
            <state>(.+?)<\/state>/gxs){
        my $jbname = $1;
        my $jbstat = $2;
#       print "$1\t$2\n";
        return $jbname if($jbname =~ /$opt_N\_\d+\_$shell_name$/ and $jbstat eq "qw");
    }
    return "NO";
}

sub dead_nodes{
    my $dead = shift;
    chomp(my @nodeMess = `qhost`);
    shift @nodeMess for(1..3);
    foreach(@nodeMess){
        my @temp = split;
        my $node_name = $temp[0];
        $dead->{$node_name} = undef if($temp[3]=~/-/ || $temp[5]=~/-/ || $temp[7]=~/-/ || $temp[4]=~/-/ || $temp[6]=~/-/);
    }
}

sub usage{
    print <<EOD
usage: perl $0 [options] shell.sh
    Options:
        -d  qsub script and log dir, default ./shell.sh_qsub/
        -l  the qsub -l option argument: vf=xxG[,p=xx,...] (default vf=1G)
        -q  queue list, default all availabile queues
        -N  set the prefix tag for qsubed jobs, default work
        -P  project_name, default not
        -n  compute node, default all availabile nodes
        -b  set number of lines to form a job, default 1
        -m  set the maximum number of jobs to throw out, default 30
        -s  set interval time of checking by qstat, default 120 seconds
        -r  mark to reqsub the job which was finished error, max reqsub 10 times, default not
        -h  show this help
EOD
}
