#!/usr/bin/perl -w
use strict;
die "perl $0 <in> <out> <sv>"unless @ARGV==3;

my $bedtools = "miniconda3/bin/bedtools";

my %trf;my %rp;

overlapf("$ARGV[0]/hg19.diff.bed", \%rp);

sub overlapf {
  my ($f, $h) = (@_);
  open IN, "$bedtools intersect -a $ARGV[2] -b $f -wao |";
  while(<IN>){
  chomp;
  my @a = split;
  next if $a[5] eq ".";
  my $k = "$a[0]:$a[1]:$a[2]:$a[3]";
  my $t = $a[8];
  $$h{$k}{$t} += $a[-1];
  }
  close IN;
}

open IN, "$ARGV[2]";
open OT,">$ARGV[1]";
while(<IN>){
  chomp;
  my @a = split;
  my $k = "$a[0]:$a[1]:$a[2]:$a[3]";
  my $len = $a[2] - $a[1] + 1;
  my $flag = "Easy";
  foreach my $t (sort keys %{$rp{$k}}){
    next if $rp{$k}{$t} / $len < 0.5;
    $flag = $t if $flag eq "Easy";
  }
  print OT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$flag\n";
}
close IN;
close OT;


