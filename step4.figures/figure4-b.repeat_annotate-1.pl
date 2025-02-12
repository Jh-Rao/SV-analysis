#!/usr/bin/perl -w
use strict;
die "perl $0 <sv> <anno> <out> <svt>"unless @ARGV==4;

my $bedtools = "miniconda3/bin/bedtools";

my %rp;
my $size = 50;
my $frate = 0.3;

$ARGV[3] =~ /INS/ ? overlapins($ARGV[1], \%rp) : overlapf("$ARGV[1]", \%rp);

sub overlapins{
  my ($f, $h) = (@_);
  open IN, "$bedtools closest -a $ARGV[0] -b $f -d |";
  while(<IN>){
    chomp;
    my @a = split;
    next if $a[5] eq ".";
    next if $a[-1] > $size;
    my $k = "$a[0]:$a[1]:$a[2]:$a[3]";
    my $t = $a[8] =~ /trf/ ? "TRF" :
          $a[8] =~ /Alu/ ? "Alu" :
          $a[8] =~ /LINE.L1/  ? "L1"  :
          $a[8] =~ /HERV/  ? "HERV"  :
          $a[8] =~ /SVA/ ? "SVA" :
          $a[8] =~ /LTR/ ? "LTR" :
          $a[8] =~ /Low_complexity/ ? "Complex" :
          $a[8] =~ /Simple_repeat/  ? "STR" : "Other";
    $$h{$k}{$t} = $a[-1];
  }
  close IN;
}

sub overlapf {
  my ($f, $h) = (@_);
  open IN, "$bedtools intersect -a $ARGV[0] -b $f -wao |";
  while(<IN>){
  chomp;
  my @a = split;
  next if $a[5] eq ".";
  my $k = "$a[0]:$a[1]:$a[2]:$a[3]";
  my $t = $a[8] =~ /trf/ ? "TRF" :
          $a[8] =~ /Alu/ ? "Alu" :
          $a[8] =~ /LINE.L1/  ? "L1"  :
          $a[8] =~ /HERV/  ? "HERV"  :
          $a[8] =~ /SVA/ ? "SVA" :
          $a[8] =~ /LTR/ ? "LTR" :
          $a[8] =~ /Low_complexity/ ? "Complex" :
          $a[8] =~ /Simple_repeat/  ? "STR" : "Other";
  $$h{$k}{$t} += $a[-1];
  }
  close IN;
}

open IN, "$ARGV[0]";
open OT,">$ARGV[2]";
while(<IN>){
  chomp;
  my @a = split;
  my $k = "$a[0]:$a[1]:$a[2]:$a[3]";
  my $len = $a[2] - $a[1] + 1;
  my @flag;
  foreach my $t (sort keys %{$rp{$k}}){
    next if $rp{$k}{$t} / $len < 0.5 && $ARGV[3] ne "INS";
    push @flag, $t;
  }
  print OT "$a[0]\t$a[1]\t$a[2]\t".(join ".",@flag)."\n";
}
close IN;
close OT;


