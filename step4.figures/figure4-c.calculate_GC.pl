#!/usr/bin/perl -w
use strict;
use POSIX;
use Getopt::Long;
use File::Basename;
use threads;
use threads::shared;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path);
use Statistics::Descriptive;

#=========================================================================#
my (%Par, $command, %command);
my %func = (
  seq => \&SEQ,
  fq  => \&FQ,
  fa  => \&FA,
  bam => \&BAM,
  subseqfa => \&SUBSEQFA,
  report   => \&REPORT,
  subseqfq => \&SUBSEQFQ,
  subseqbam=> \&SUBSEQBAM,
);
&SetDefault();
&ParameterReading();
&{$func{$Par{inputtype}}};
system("rm -fr $Par{outputdir}/tmp.*");
exit(0);

#=========================================================================#
sub BAM {
  &PrintInfo2ERR("Start calculate GC content of BAM: $Par{input}");
  `mkdir -p $Par{outputdir}/tmp.roi $Par{outputdir}/tmp.sam $Par{outputdir}/tmp.gc`;
  if( $Par{region} ){
    `$Par{samtools} view $Par{input} $Par{region} > $Par{outputdir}/tmp.sam/splitROIaa`;
  }else{
    &PrintInfo2ERR("build tmp directory for ROI");
    $Par{regionlist} ? `less $Par{regionlist} | awk 'BEGIN{OFS=\"\\t\"}{a=\$2-1; print \$1,a,\$3}' > $Par{outputdir}/tmp.roi/roi.bed` :
                       `ln -s $Par{input}.fai  $Par{outputdir}/tmp.roi/roi.bed`;
    my $line_roi = `wc -l $Par{outputdir}/tmp.roi/roi.bed | awk '{print \$1}'`;
    chomp $line_roi;
    my $line_split = $line_roi > $Par{thread} ? int($line_roi / $Par{thread}) : 1;
    `$Par{split} --lines=$line_split $Par{outputdir}/tmp.roi/roi.bed $Par{outputdir}/tmp.roi/splitROI `;
    foreach (`ls $Par{outputdir}/tmp.roi | grep splitROI`){
      chomp;
      push @{$command{0}}, "for i in `cat $Par{outputdir}/tmp.roi/$_ | awk '{print \$1\":\"\$2\"-\"\$3}'`;do $Par{samtools} view $Par{input} \$i ;done > $Par{outputdir}/tmp.sam/$_ \n";
    }
  }
  foreach (`ls $Par{outputdir}/tmp.roi | grep splitROI`){
    chomp;
    push @{$command{1}}, "$0 -i $Par{outputdir}/tmp.sam/$_ -it subseqbam -o $Par{outputdir}/tmp.gc/$_ \n";
  }
  push @{$command{2}}, "$0 -i $Par{outputdir}/tmp.gc -it report -o $Par{outputdir}/ \n";
  &runcommand();
  &PrintInfo2ERR("GC content calcualtation done!");
}

sub SUBSEQBAM {
  open I, $Par{input};
  open O,">$Par{outputdir}";
  while(<I>){
    chomp;
    my @a = split;
    my $gc = calcgc($a[9]);
    my $len = length($a[9]);
    print O "$a[0]\t$len\t$gc\n";
  }
  close I;
  close O;
} 

sub FQ {
  &PrintInfo2ERR("Start calculate GC content of FastQ: $Par{input}");
  `mkdir -p $Par{outputdir}/tmp.fq $Par{outputdir}/tmp.seq $Par{outputdir}/tmp.gc`;
  # split FQ
  my $line = $Par{input} =~ /gz$/ ? `gzip -cd $Par{input} | wc -l | awk '{print \$1}'` :
                                    `wc -l $Par{input} | awk '{print \$1}'`;
  chomp $line;
  $line /= 4;
  my $line_split = $line > $Par{thread} ? int($line / $Par{thread}) : 1;
  $line_split *= 4;
  $Par{input} =~ /gz$/ ? `gzip -cd $Par{input} > $Par{outputdir}/tmp.fq/input.fq` :
                         `ln -s $Par{input} $Par{outputdir}/tmp.fq/input.fq`;
  `$Par{split} --lines=$line_split $Par{outputdir}/tmp.fq/input.fq $Par{outputdir}/tmp.seq/splitFQ `;
  
  # calculate GC and Qual for each split-FQ
  foreach (`ls $Par{outputdir}/tmp.seq | grep splitFQ`){
    chomp;
    push @{$command{0}}, "$0 -i $Par{outputdir}/tmp.seq/$_ -it subseqfq -o $Par{outputdir}/tmp.gc/$_ \n";
  }

  # collect GC results
  push @{$command{1}}, "$0 -i $Par{outputdir}/tmp.gc -it report -o $Par{outputdir}/ \n";
  &runcommand();
  &PrintInfo2ERR("GC content calcualtation done!");
}

sub SUBSEQFQ {
  open I, $Par{input};
  open O,">$Par{outputdir}";
  while(my $id = <I>){
    chomp $id;
    my $seq = <I>;
    chomp $seq;
    <I>;
    my $qua = <I>;
    chomp $qua;
    $id = (split /\s+/, $id)[0];
    my $gc = calcgc($seq);
    my $len = length($seq);
    my $qual = calcqual($qua) if $Par{qualcalc} == 1;
    print O "$id\t$len\t$gc";
    $Par{qualcalc} == 1 ? print O "\t$qual\n" : print O "\n";
  }
  close I;
  close O;
}

sub FA {
  # extract a single sequence
  if( $Par{region} ){
    my $seq = `$Par{samtools} faidx $Par{input} $Par{region} | grep -v '^>' | awk '{printf "%s",\$0}'`;
    chomp $seq;
    #&PrintInfo2ERR("The Seqeunce of region is \n\t>>>$seq>>>");
    print STDOUT calcgc($seq)."%\n";
  }else{
    &PrintInfo2ERR("Start calculate GC content of FastA: $Par{input}");
    &PrintInfo2ERR("build tmp directory for ROI");
    `mkdir -p $Par{outputdir}/tmp.roi $Par{outputdir}/tmp.seq $Par{outputdir}/tmp.gc`;
    # region of interset = BED ? BED : FA.fai
    $Par{regionlist} ? `less $Par{regionlist} | awk 'BEGIN{OFS=\"\\t\"}{a=\$2-1; print \$1,a,\$3}' > $Par{outputdir}/tmp.roi/roi.bed` :
                       `ln -s $Par{input}.fai  $Par{outputdir}/tmp.roi/roi.bed`;
    my $line_roi = `wc -l $Par{outputdir}/tmp.roi/roi.bed | awk '{print \$1}'`;
    chomp $line_roi;
    my $line_split = $line_roi > $Par{thread} ? int($line_roi / $Par{thread}) : 1;
    # split ROI according to thread, we cannot process 1 line per thread, because >1k files in one directory is unacceptable
    `$Par{split} --lines=$line_split $Par{outputdir}/tmp.roi/roi.bed $Par{outputdir}/tmp.roi/splitROI `;
    
    # calculate GC for each ROI file
    foreach (`ls $Par{outputdir}/tmp.roi | grep splitROI`){
      chomp;
      # extract sequence of ROI, there may be many ROI in one file, so using seqtk subseq
      push @{$command{0}}, "$Par{seqtk} subseq $Par{input} $Par{outputdir}/tmp.roi/$_ > $Par{outputdir}/tmp.seq/$_  \n";
      # calcualte GC for each sequence of ROI
      push @{$command{1}}, "$0 -i $Par{outputdir}/tmp.seq/$_ -it subseqfa -o $Par{outputdir}/tmp.gc/$_ \n";
    }
    # collect GC results
    push @{$command{2}}, "$0 -i $Par{outputdir}/tmp.gc -it report -o $Par{outputdir}/ \n";
    &runcommand();
    &PrintInfo2ERR("GC content calcualtation done!");
  }
}

sub REPORT {
  my %stat;
  my (@len, @gc, @gcc, @qual);
  open O,">$Par{outputdir}/GCcontent.detail.xls";
  print O "region\tlength\tGC";
  $Par{inputtype} =~ /fq|bam/ ? print O "\tQual\n" : print O "\n";
  foreach my $gcf (`ls $Par{input}/`){
    chomp $gcf;
    open I, "$Par{input}/$gcf";
    while(<I>){
      chomp;
      my @a = split;
      push @gc, $a[2];
      push @len, $a[1];
      push @gcc, $a[1] * $a[2] / 100;
      push @qual, $a[3] if @a == 4;
      print O "$_\n";
      my $k = int($a[2] / $Par{bin});
      $stat{$k}[0] += 1;
      $stat{$k}[1] += $a[1];
    }
    close I;
  }
  close O;

  my @gcstat = arraystat(\@gc);
  my @lenstat = arraystat(\@len);
  my @gccstat = arraystat(\@gcc);
  my $avegc = sprintf("%.2f", $gccstat[1] / $lenstat[1] * 100);
  open O,">$Par{outputdir}/GCcontent.brief.xls";
  print O "Task\t$Par{task}
SequenceCount\t$lenstat[0]
SequenceLength\t$lenstat[1]
SequenceMeanLegnth\t$lenstat[4]
AverageGC(%)\t$avegc
MeanGC(%)\t$gcstat[4]
MedianGC(%)\t$gcstat[5]
MinGC(%)\t$gcstat[2]
MaxGC(%)\t$gcstat[3]\n";
  for(my $i = 0; $i < int(100 / $Par{bin}); $i++){
    my $s = $Par{bin} * $i;
    my $e = $s + $Par{bin};
    my $v = defined $stat{$i}[1] ? $stat{$i}[1] : "0";
    print O "$s-$e\t$v\n";
  }
  close O;

  `$Par{Rscript} $FindBin::RealBin/GCdistribution.R $Par{outputdir}/GCcontent.detail.xls $Par{outputdir}/GCcontent.pdf $avegc $gcstat[4]`;
}

sub SUBSEQFA {
  my ($chr, $seq, $len, $gc);
  $seq = "";
  open I, $Par{input};
  open O,">$Par{outputdir}";
  while(<I>){
    chomp;
    if(/^>/){
      #$chr = (split /\s+/)[0];
      #$chr =~ s/^>//;
	  if($chr =~ /\S+/){
        &PrintInfo2ERR("calculate GC content of $chr");
		#print "A\$chr\t$seq\n";
		$gc = calcgc($seq);
	    $len = length($seq);
        $seq = "";
        print O "$chr\t$len\t$gc\n";
	  }
      $chr = (split /\s+/)[0];
      $chr =~ s/^>//;
      next;
    }
    $seq .= $_;
    #&PrintInfo2ERR("calculate GC content of $chr");
    #$gc = calcgc($_);
    #$len = length($_);
    #print O "$chr\t$len\t$gc\n";
  }
  close I;
  close O;
}

sub SEQ {
  print STDOUT calcgc($Par{input})."%\n";
}

sub calcgc{
  my $seq = $_[0];
  my $count = 0;
  $count++ while $seq =~ m/[$Par{letter}]/ig;
  #print "CG=$count\n";
  my $gc = sprintf("%.2f", $count / length($seq) * 100);
  #print "len=".length($seq)."\n";
  return $gc;
}

sub calcqual {
  my @qual;
  my $add = $Par{qualsys} == 1 ? 33 : 64;
  foreach (split //, $_[0]){
    chomp;
    push @qual, ord($_) - $add;
  }
  my @out = arraystat(\@qual);
  return $out[4];
}

sub arraystat{
  my ($array) = (@_);
  my $stat = Statistics::Descriptive::Full->new();
  $stat->add_data(@$array);
  my $count = $stat->count();
  my $length = $stat->sum();
  my $min = $stat->min();
  my $max = $stat->max();
  my $mean = sprintf("%.2f", $stat->mean());
  my $median = $stat->median();
  return ($count, $length, $min, $max, $mean, $median);
}

#=========================================================================#
# basic monitor subs
#=========================================================================#
sub RunInParallel{
  my ($command, $index, $maxThread) = (@_);

  my $threadid = 0;
  my $commandcount = @$command;
  &PrintInfo2ERR("Start $index of $commandcount tasks with $maxThread in parallel ...");

  while(@{$command} > 0){
    if(scalar(threads->list(threads::all)) < $maxThread){
      my $oneCommand = shift @{$command};
      $threadid++;
      threads->new(\&RunCommand, $oneCommand, $index, $threadid);
    }else{
      sleep($Par{"sleeptime"});
    }
    foreach my $thread (threads->list(threads::all)){
      $thread->join() if $thread->is_joinable();
    }
  }

  while(scalar(threads->list(threads::all)) > 0){
    foreach my $thread (threads->list(threads::all)){
      $thread->join() if $thread->is_joinable();
    }
  }
  &PrintInfo2ERR("Finish $index inParallel");
}

sub RunCommand{
  my ($oneCommand, $index, $threadid) = (@_);
  &PrintInfo2ERR("Start $index inParallel $threadid ...");
  chomp $oneCommand;
  print "-----> $index $threadid ::: $oneCommand\n" if $Par{showcom} eq "Y";
  system("$oneCommand") if $Par{run} eq "Y";
  &PrintInfo2ERR("Finish $index inParallel $threadid ");
}

sub PrintInfo2ERR {
  my ($msg) = (@_);
  my $logTime =  strftime("[%Y-%m-%d %H:%M:%S]", localtime());
  print STDOUT "$logTime INFO: $msg\n";
}

sub runcommand{
  foreach my $kid (sort {$a <=> $b} keys %command){
    &RunInParallel(\@{$command{$kid}}, "step $kid", $Par{thread});
  }
}

sub ParameterReading{
  GetOptions(
    "o=s"           => \$Par{outputdir},
    "i=s"           => \$Par{input},
    "it=s"          => \$Par{inputtype},
    "rl=s"          => \$Par{regionlist},
    "r=s"           => \$Par{region},
    "task=s"        => \$Par{task},
    "h|help"        => \$Par{help},
    "t=i"           => \$Par{thread},
    "l=s"           => \$Par{letter},
    "b=s"           => \$Par{bin},
    "qc=s"          => \$Par{qualcalc},
    "qs=s"          => \$Par{qualsys},

    "run=s"         => \$Par{run},
    "sleeptime=i"   => \$Par{sleeptime},
    "showcom=s"     => \$Par{showcom},

    "seqtk=s"       => \$Par{seqtk},
    "samtools=s"    => \$Par{samtools},
    "split=s"       => \$Par{split},
  );

  my $usage = "
  Program: calgc (calculate GC content)
  Version: $Par{version}
  Usage  : calgc [options]

  Options:

      -i    *   input
      -it       input format, one of seq,fa,fq,bam [ seq ]
      -o        output directory, useless when input type is 'seq' [ ./ ]
      -b        bin size for GC output [ 5 ]
      -r        target region in format: chr:start-end, useful when input type is fa|bam
      -rl       target region list in BED format, one region per line of three columns:  chr  start   end
      -l        letter for calculate [ CG ]
      -task     task name [ calculateGC ]
      -t        thread [ 10 ]
      -h|help   help
  ";

  die "$usage \n" if defined $Par{help} || !defined $Par{input};
  $Par{outputdir} = abs_path($Par{outputdir});
  $Par{input}     = abs_path($Par{input}) if $Par{inputtype} ne "seq";
  $Par{regionlist}= abs_path($Par{regionlist}) if $Par{regionlist};
};

sub SetDefault{
  %Par = (
    "outputdir"  => "./",
    "task"       => "CalculateGC",
    "version"    => "1.0",
    "inputtype"  => "seq",
    "bin"        => 5,
    "letter"     => "CG",
    "qualcalc"   => 0,
    "qualsys"    => 1,
    "thread"     => 10,
    "run"        => "Y",
    "sleeptime"  => 30,
    "showcom"    => "N",

    "seqtk"      => "$FindBin::RealBin/seqtk",
    "samtools"   => "$FindBin::RealBin/samtools",
    "split"      => "$FindBin::RealBin/split",
    "Rscript"    => "$FindBin::RealBin/Rscript",
  );
};

