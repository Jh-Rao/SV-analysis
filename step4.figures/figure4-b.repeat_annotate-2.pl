#!/usr/bin/perl -w
use strict;

die "perl $0 <rpdir> <homerdir> <pl>"unless @ARGV==3;

my %repeat;
foreach my $svt (qw/DEL DUP INS INV/){
  foreach my $anno (qw/rmsk sva trf/){
	open IN, "$ARGV[0]/$ARGV[2].$anno.$svt.out";
	while(<IN>){
	  chomp;
	  my @a = split;
      next if @a == 3;
	  my $s = $a[1] + 1;
	  my $k = "$a[0]:$s-$a[2]";
	  $repeat{$k} .= "$a[3].";
	}
	close IN;
  }
}

open IN, "$ARGV[1]/$ARGV[2].output.bed";
open OT,">$ARGV[0]/$ARGV[2].bed";
while(<IN>){
  chomp;
  my @a = split /\t/;
  next if /^PeakID/;
  my $svt = $a[0];
  $svt =~ s/\d+$//;
  my $k = "$a[1]:$a[2]-$a[3]";
  my $s = $a[2] - 1;
  my $flag = $a[8] =~ /n|Simple_repeat/ ? "STR" :
             $a[8] =~ /Alu$/ ? "Alu" :
             $a[8] =~ /Low_complexity/ ? "Complex" :
             $a[8] =~ /LINE|L1/ ? "L1" :
             $a[8] =~ /SVA/ ? "SVA" :
             $a[8] =~ /HERV/ ? "HERV" :
             $a[8] =~ /LTR/ ? "LTR" : 
             $a[8] =~ /UTR|CpG|exon|Intergenic|intron|non-coding|TSS|Unknown/ ? "NoRepeat" : "Other";
  $repeat{$k} .= "$flag.";

  $flag = "NoRepeat";
  foreach my $t (qw/TRF STR Alu L1 SVA HERV LTR Complex Other/){
    if($repeat{$k} =~ /$t/){
      $flag = $t;
      last;
    }
  }
  print OT "$a[1]\t$s\t$a[3]\t$svt\t$flag\n";
}
close OT;
close IN;
