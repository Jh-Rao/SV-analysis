#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;

die "perl $0 <indir> <out>"unless @ARGV==2;

my %softtype;
open IN, "softlist.xls";
while(<IN>){
  chomp;
  my @a = split;
  $a[1] = uc($a[1]);
  $softtype{$a[1]} = $a[0];
}
close IN;

my (%count, %length);
foreach my $file (`ls $ARGV[0]/*.bed`){
  chomp $file;
  my $base = `basename $file`;
  chomp $base;
  my @base = split /\./, $base;
  my $soft = $base[0];
  my $sam = $base[1];

  open IN, $file;
  while(<IN>){
    chomp;
    my @a = split;
    my $len = $a[2] - $a[1] + 1;
    $a[3] = $a[3] eq "MEI" ? "INS" : $a[3];
    next if $a[3] !~ /DEL|DUP|INS|INV|TRA/;
    $count{$soft}{$a[3]}{$sam}++;
    $length{$soft}{$a[3]}{$sam}+=$len;
  }
  close IN;

}

open OT,">$ARGV[1]";
print OT "Software\tAlgorithm\tKey\tType\tDNBSEQ\tIllumina\n";
foreach my $soft (sort keys %count){
  my $k = uc($soft);
  my $method = $softtype{$k};

  foreach my $svt (sort keys %{$count{$soft}}){
    my (@mgic, @illc, @mgil, @illl) = () x 4;
    my ($mc, $ml, $ic, $il) = 0 x 4;

    foreach my $sam (sort keys %{$count{$soft}{$svt}}){
      if($sam =~ /GI/){
        push @mgic, $count{$soft}{$svt}{$sam};
        push @mgil, $length{$soft}{$svt}{$sam} if defined $length{$soft}{$svt}{$sam};
      }else{
        push @illc, $count{$soft}{$svt}{$sam};
        push @illl, $length{$soft}{$svt}{$sam} if defined $length{$soft}{$svt}{$sam};
      }
    }

    $mc = array_stat(\@mgic) if @mgic > 0;
    $ml = array_stat(\@mgil) if @mgil > 0;
    if($soft eq "inGAP-sv" && $svt eq "INS"){
      $ic = 0;
      $il = 0;
    }else{
      $ic = array_stat(\@illc) if @illc > 0;
      $il = array_stat(\@illl) if @illl > 0;
    }

    print OT "$soft\t$method\tNumber\t$svt\t$mc\t$ic\n";
    print OT "$soft\t$method\tLength\t$svt\t$ml\t$il\n" if @mgil > 0 && $svt ne "TRA";
  }
}
close OT;

sub array_stat{
  my ($array) = (@_);
  my $stat    = Statistics::Descriptive::Full->new();
  $stat->add_data(@$array);
  my $count   = $stat->count();
  my $mean    = sprintf("%.4f", $stat->mean());
  my $median  = sprintf("%.4f", $stat->median());
  my $sum     = $stat->sum();
  my $sd      = $stat->standard_deviation();
  return $mean;
}

