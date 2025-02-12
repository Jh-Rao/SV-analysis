#!/usr/bin/perl -w
use strict;

die "perl $0 <in1> <in2> <name2> <svtype> <method: u|n> <out>"unless @ARGV==6;
my ($file1, $file2, $name2, $svtype, $method, $out) = @ARGV;
# input contains 5 columns

my $ins_sd = 200;
my $rate = 0.5;
my (%sv1, %flag1, %sv2, %flag2, %ovlap1, %ovlap2);

# read SV
&readSV($file1, \%sv1, \%flag1);
&readSV($file2, \%sv2, \%flag2);

# overlap
foreach my $chr (sort keys %sv1){
  next if !$sv2{$chr};
  for(my $i = $flag1{$chr}; $i < @{$sv1{$chr}}; $i++){
    for(my $j = $flag2{$chr}; $j < @{$sv2{$chr}}; $j++){
      if($svtype =~ /INS/){
        if($sv1{$chr}[$i][1] < $sv2{$chr}[$j][1] - $ins_sd){ last; }
        if($sv1{$chr}[$i][1] > $sv2{$chr}[$j][1] + $ins_sd){ $flag2{$chr} = $j; next; }
        my $k1 = "$sv1{$chr}[$i][0]:$sv1{$chr}[$i][1]:$sv1{$chr}[$i][2]";
        my $k2 = "$sv2{$chr}[$j][0]:$sv2{$chr}[$j][1]:$sv2{$chr}[$j][2]";
        $ovlap1{$k1} = 1;
        $ovlap2{$k2} = 1;
      }else{
        if($sv1{$chr}[$i][2] < $sv2{$chr}[$j][1]){ last; }
        if($sv1{$chr}[$i][1] > $sv2{$chr}[$j][2]){ $flag2{$chr} = $j; next; }
        my $s = $sv1{$chr}[$i][1] > $sv2{$chr}[$j][1] ? $sv1{$chr}[$i][1] : $sv2{$chr}[$j][1];
        my $e = $sv1{$chr}[$i][2] < $sv2{$chr}[$j][2] ? $sv1{$chr}[$i][2] : $sv2{$chr}[$j][2];
        my $ovlen = $e - $s + 1;
        my $r1 = $ovlen / ($sv1{$chr}[$i][2] - $sv1{$chr}[$i][1] + 1);
        my $r2 = $ovlen / ($sv2{$chr}[$j][2] - $sv2{$chr}[$j][1] + 1);
        next if $r1 < $rate || $r2 < $rate;
        my $k1 = "$sv1{$chr}[$i][0]:$sv1{$chr}[$i][1]:$sv1{$chr}[$i][2]";
        my $k2 = "$sv2{$chr}[$j][0]:$sv2{$chr}[$j][1]:$sv2{$chr}[$j][2]";
        $ovlap1{$k1} = 1;
        $ovlap2{$k2} = 1;
      }
    }
  }
}

# output
open OT,">$out";
&outputSV($file1, \%sv1, \%ovlap1, \*OT, 1);
&outputSV($file2, \%sv2, \%ovlap2, \*OT, 2) if $method eq "u";
close OT;

sub outputSV{
  my ($f, $sv, $ovlap, $fh, $flag) = (@_);

  open IN, $f;
  while(<IN>){
    chomp;
    my @a = split;
    next if $a[3] ne $svtype;
    my $k = "$a[0]:$a[1]:$a[2]";
    $$ovlap{$k} ||= 0;
    next if $method eq "n" && !$$ovlap{$k};
    next if $method eq "u" && $$ovlap{$k} && $flag == 2;
    $a[4] .= ":$name2" if $$ovlap{$k} && $flag == 1;
    $a[4] = $name2 if $flag == 2;
    print OT (join "\t", @a)."\n";
  }
  close IN;

=head
  foreach my $chr (sort keys %$sv){
    for(my $i = 0; $i < @{$$sv{$chr}}; $i++){
      my $k = "$$sv{$chr}[$i][0]:$$sv{$chr}[$i][1]:$$sv{$chr}[$i][2]";
      next if $method eq "n" && !$$ovlap{$k};
      next if $method eq "u" && $$ovlap{$k} && $flag == 2;
      $$sv{$chr}[$i][4] .= ":$name2" if $$ovlap{$k} && $flag == 1;
      print OT (join "\t", @{$$sv{$chr}[$i]})."\n";
    }
  }
=cut

}

sub readSV {
  my ($file, $hash1, $hash2) = (@_);
  open IN, $file;
  while(<IN>){
	chomp;
	my @a = split;
    next if $a[3] ne $svtype;
	push @{$$hash1{$a[0]}}, [@a];
    $$hash2{$a[0]} = 0;
  }
  close IN;
}

