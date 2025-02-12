#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;

die "perl $0 <indir> <out>"unless @ARGV==2;

my %softtype;
my (%tn, @tn);
open IN, "softlist.xls";
while(<IN>){
  chomp;
  my @a = split;
  $a[1] = uc($a[1]);
  $softtype{$a[1]} = $a[0];
	push @tn, $a[1];
	$tn{$a[1]} = $a[2];
}
close IN;

my @samples;
open OT,">$ARGV[1].size_dis.ggplot.xls";
print OT "Software\tAlgorithm\tSample\tPlatform\tType\tSize\n";
my (%count, %length);
foreach my $file (`ls $ARGV[0]/*.bed`){
  chomp $file;
  next if $file =~ /T7/;
  my $base = `basename $file`;
  chomp $base;
  my @base = split /\./, $base;
  my $soft = $base[0];
  my $sam = $base[1];
	push @samples, $sam;

  my $k = uc($soft);
  my $method = $softtype{$k};
  my $pl = $sam =~ /GI/ ? "DNBSEQ" : "Illumina";

  open IN, $file;
  while(<IN>){
    chomp;
    my @a = split;
    my $len = $a[2] - $a[1] + 1;
    $a[3] = $a[3] eq "MEI" ? "INS" : $a[3];
    next if $a[3] !~ /DEL|DUP|INS|INV|TRA/;
    $count{$k}{$a[3]}{$sam}++;
    $length{$k}{$a[3]}{$sam}+=$len;

    print OT "$tn{$k}\t$method\t$sam\t$pl\t$a[3]\t$len\n";
  }
  close IN;

}
close OT;

build_stat_per_dataset("number", \%count);
build_stat_per_dataset("length", \%length);

sub build_stat_per_dataset {
	my ($flag, $hash) = (@_);

	foreach my $sam (@samples){
	  open OT,">$ARGV[1].$sam.$flag.txt";
	  for(my $i = 0; $i < @tn; $i++){
		  my $k = $i + 1;
		  $k = "hs$k";

		  next if !$$hash{$tn[$i]}{DEL}{$sam} && !$$hash{$tn[$i]}{DUP}{$sam} && !$$hash{$tn[$i]}{INS}{$sam} && !$$hash{$tn[$i]}{INV}{$sam} && !$$hash{$tn[$i]}{TRA}{$sam};

		  my $del = $$hash{$tn[$i]}{DEL}{$sam} ? $$hash{$tn[$i]}{DEL}{$sam} : 0;
		  my $dup = $$hash{$tn[$i]}{DUP}{$sam} ? $$hash{$tn[$i]}{DUP}{$sam} : 0;
		  my $ins = $$hash{$tn[$i]}{INS}{$sam} ? $$hash{$tn[$i]}{INS}{$sam} : 0;
		  my $inv = $$hash{$tn[$i]}{INV}{$sam} ? $$hash{$tn[$i]}{INV}{$sam} : 0;
		  my $tra = $$hash{$tn[$i]}{TRA}{$sam} ? $$hash{$tn[$i]}{TRA}{$sam} : 0;
		  $tra = 0 if $flag eq "length";

		  print OT "$k\t0\t1\t$del,$dup,$ins,$inv,$tra\n";
	  }
	  close OT;
	}
}


open OT,">$ARGV[1]";
print OT "Software\tAlgorithm\tPlatform\tType\tKey\tMean\tSD\tSEM\n";
foreach my $soft (sort keys %count){
  my $k = uc($soft);
  my $method = $softtype{$k};

  foreach my $svt (sort keys %{$count{$soft}}){
    my (@mgic, @illc, @mgil, @illl) = () x 4;

    foreach my $sam (sort keys %{$count{$soft}{$svt}}){
      if($sam =~ /GI/){
        push @mgic, $count{$soft}{$svt}{$sam};
        push @mgil, $length{$soft}{$svt}{$sam} if $length{$soft}{$svt}{$sam};
      }else{
        push @illc, $count{$soft}{$svt}{$sam};
        push @illl, $length{$soft}{$svt}{$sam} if $length{$soft}{$svt}{$sam};
      }
    }

    my @out = (0) x 12;

    ($out[0], $out[1],  $out[2])  = array_stat(\@mgic);
    ($out[3], $out[4],  $out[5])  = array_stat(\@mgil) if @mgil > 0;
    ($out[6], $out[7],  $out[8])  = array_stat(\@illc);
    ($out[9], $out[10], $out[11]) = array_stat(\@illl) if @illl > 0;

    print OT "$tn{$soft}\t$method\tDNBSEQ\t$svt\tNumber\t$out[0]\t$out[1]\t$out[2]\n";
    print OT "$tn{$soft}\t$method\tDNBSEQ\t$svt\tLength\t$out[3]\t$out[4]\t$out[5]\n";
    print OT "$tn{$soft}\t$method\tIllumina\t$svt\tNumber\t$out[6]\t$out[7]\t$out[8]\n";
    print OT "$tn{$soft}\t$method\tIllumina\t$svt\tLength\t$out[9]\t$out[10]\t$out[11]\n";
  }
}
close OT;

sub array_stat{
  my ($array) = (@_);
print "A\t@$array\n";
  my $stat    = Statistics::Descriptive::Full->new();
  $stat->add_data(@$array);
  my $count   = $stat->count();
  my $mean    = sprintf("%.2f", $stat->mean());
  my $sd      = sprintf("%.2f", $stat->standard_deviation());
  my $sem     = sprintf("%.2f", $sd / sqrt($count));
  return ($mean, $sd, $sem);
}

