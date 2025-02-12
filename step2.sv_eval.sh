#===========================================================#
# Below are the commands for evaluating the SVs detected in BASIL-ANISE across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=basil_anise
ntool=BASIL-ANISE

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/basil.vcf \
    > $dir/$i.vcf

  perl $epl $dir/$i.vcf
done

#===========================================================#
# Below are the commands for evaluating the SVs detected in BICseq2 across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=BICseq2
ntool=BICseq2

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/chr*.cnv.txt` \
    > $i.vcf

  perl $epl $i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in BreakDancer across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=breakdancer
ntool=BreakDancer

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/chr*.cnv.txt` \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in BreakSeek across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=breakseek
on=$otool
ntool=BreakSeek

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/*INDEL_list.txt` \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in BreakSeq2 across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=breakseq2
ntool=BreakSeq2

for i in `ls $dir/../data/alignment`
do
  gzip -cd $dir/$i/breakseq.vcf.gz \
     > $dir/$i/breakseq.vcf
  perl $cpl \
    -t $ntool \
    $dir/$i/breakseq.vcf \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in CNVnator across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=cnvnator
ntool=CNVnator

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/cnv.txt \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in Control-FREEC across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=freec
ntool=Control-FREEC

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/chr*/*CNVs` \
    > $i.vcf

  perl $epl $i.vcf
done


#===========================================================#
# Below are the commands for evaluating the SVs detected in DELLY across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=delly
ntool=DELLY

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/*bcf` \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in DINUMT across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=dinumt
ntool=DINUMT

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/dinumt.vcf \
    > $dir/$i.vcf

  perl $epl \
    -st INS \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in ERDS across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=erds
ntool=ERDS

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/chr*.raw.vcf` \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in FermiKit across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=fermikit
ntool=FermiKit

for i in `ls $dir/../data/alignment`
do
  gzip -cd $dir/$i/prefix.sv.vcf.gz \
    > $dir/$i/prefix.sv.vcf
  perl $cpl \
    -t $ntool \
    $dir/$i/prefix.sv.vcf \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in GASVpro across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=gasvpro
ntool=GASVpro

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/chr*/*gasv.in.clusters.GASVPro.clusters.pruned.clusters` \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in GenomeSTRiP2 across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=genomestrip
ntool=GenomeSTRiP2

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/svdiscovery.dels.vcf \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in GRIDSS across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=gridss
ntool=GRIDSS

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/$i.sv.vcf \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done

#===========================================================#
# Below are the commands for evaluating the SVs detected in Hydra-sv across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=hydrasv
ntool=Hydra-sv

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/all.${i}.sv.final \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in iCopyDAV across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=icopydav
ntool=iCopyDAV

dir=`pwd`

  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/chr*/*_tvm_pCNVD.bed` \
    > $i.vcf

  perl $epl $i.vcf
  done

#===========================================================#
# Below are the commands for evaluating the SVs detected in IndelMINER across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=indelminer
ntool=IndelMINER

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  gzip -cd $dir/$i/$otool.final.vcf.gz \
    | grep -v ^# \
    | perl -e 'while(<>){
                  chomp;@a=split;
		  @b=split /\;|\=/, $a[7];
		  $a[4] =~ s/[\<\>]//g;
		  $a[7] =~ s/RSS/READS/;
		  next if $b[1] eq "INS";
		  next if $b[1] ne "INS" && $b[3] < 50;
		  next if $b[3] > 2000000;
		  print "$a[0]\t$a[1]\t$a[4]\t.\t.\t.\tPASS\t$a[7]\n";
               }' \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done


#===========================================================#
# Below are the commands for evaluating the SVs detected in ITIS across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=itis
ntool=ITIS

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/*/*.filtered.bed` \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done

#===========================================================#
# Below are the commands for evaluating the SVs detected in laSV across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=lasv
ntool=laSV

dir=`pwd`

  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    $dir/$i/${i}.SVs.vcf \
    > $i.vcf

  perl $epl $i.vcf
  done

#===========================================================#
# Below are the commands for evaluating the SVs detected in Lumpy across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=lumpy
ntool=Lumpy

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/${i}.vcf \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done

#===========================================================#
# Below are the commands for evaluating the SVs detected in Manta across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=manta
ntool=Manta

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/results/variants/diploidSV.vcf.gz \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in MATCHCLIP across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=matchclip
ntool=MATCHCLIP

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/${i}.out \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in Meerkat across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=meerkat
ntool=Meerkat

dir=`pwd`

  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    $dir/$i/${i}.variants \
    > $i.vcf

  perl $epl $i.vcf
  done


#===========================================================#
# Below are the commands for evaluating the SVs detected in MELT across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=melt
ntool=MELT

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/*final_comp.vcf` \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in MetaSV across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=metasv
ntool=MetaSV

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/variants.vcf.gz \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in MindTheGap across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=mindthegap
ntool=MindTheGap

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/s2.find.insertions.fasta \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in PennCNV-Seq across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=penncnvseq
ntool=PennCNV-Seq

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/*/$i.*.rawcnv` \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done

#===========================================================#
# Below are the commands for evaluating the SVs detected in Pindel across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=pindel
ntool=Pindel

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/pindel_D \
    $dir/$i/pindel_TD \
    $dir/$i/pindel_SI \
    $dir/$i/pindel_INV \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
#===========================================================#
# Below are the commands for evaluating the SVs detected in Popins across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=popins
ntool=Popins

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/*/*/insertions.vcf` \
    > $i.vcf

  perl $epl $i.vcf
done


#===========================================================#
# Below are the commands for evaluating the SVs detected in PRISM across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=prism
ntool=PRISM

dir=`pwd`

  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/*/PRISM_output/del_* |grep -v RO$ ` \
    `ls $dir/$i/*/PRISM_output/ins_* ` \
    $dir/$i/*/PRISM_output/dup \
    $dir/$i/*/PRISM_output/inv \
    > $i.vcf

  perl $epl $i.vcf
  done

#===========================================================#
# Below are the commands for evaluating the SVs detected in RetroSeq across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=retroseq
ntool=RetroSeq

dir=`pwd`

  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    $dir/$i/retroseq.vcf \
    > $i.vcf

  perl $epl $i.vcf
  done

#===========================================================#
# Below are the commands for evaluating the SVs detected in Socrates across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=socrates
ntool=Socrates

step=1
dir=`pwd`

if [ $step == "1" ];then
  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/results_Socrates_paired_*txt` \
    > $i.vcf

  perl $epl $i.vcf
  done
fi

#===========================================================#
# Below are the commands for evaluating the SVs detected in SoftSV across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=softsv
ntool=SoftSV

step=1
dir=`pwd`

if [ $step == "1" ];then
  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/*/*.txt` \
    > $i.vcf

  perl $epl $i.vcf
  done
fi

#===========================================================#
# Below are the commands for evaluating the SVs detected in Sprites across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=sprites
ntool=Sprites

step=1
dir=`pwd`

if [ $step == "1" ];then
  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/*/out` \
    > $i.vcf

  perl $epl $i.vcf
  done
fi

#===========================================================#
# Below are the commands for evaluating the SVs detected in SvABA across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=svaba
ntool=SvABA

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/$i.svaba.unfiltered.sv.vcf \
    > $i.vcf

  perl $epl $i.vcf
done

#===========================================================#
# Below are the commands for evaluating the SVs detected in SVelter across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=svelter
ntool=SVelter

step=1
dir=`pwd`

if [ $step == "1" ];then
  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    $dir/$i/${i}.vcf \
    > $i.vcf

  perl $epl $i.vcf
  done
fi

#===========================================================#
# Below are the commands for evaluating the SVs detected in TEMP across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=temp
ntool=TEMP

step=1
dir=`pwd`

if [ $step == "1" ];then
  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    $dir/$i/$i.insertion.refined.bp.summary \
    > $i.vcf

  perl $epl $i.vcf
  done
fi

#===========================================================#
# Below are the commands for evaluating the SVs detected in TIDDIT across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=tiddit
ntool=TIDDIT

step=1
dir=`pwd`

if [ $step == "1" ];then
  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    $dir/$i/${i}.vcf \
    > $i.vcf

  perl $epl $i.vcf
  done
fi

#===========================================================#
# Below are the commands for evaluating the SVs detected in Ulysses across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

otool=ulysses
ntool=Ulysses

step=1
dir=`pwd`

if [ $step == "1" ];then
  cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
  epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl
  for i in `ls $dir/../data/alignment`
  do
  perl $cpl \
    -t $ntool \
    `ls $dir/$i/*/*.vcf` \
    > $i.vcf

  perl $epl $i.vcf
  done
fi

#===========================================================#
# Below are the commands for evaluating the SVs detected in Wham across 10 datasets 
#  using the benchmark of NA12878 and evaluation method proposed by Shunichi et al. 
# The commands sequentially iterate through the 10 datasets to perform SV evaluation.
#===========================================================#

dir=`pwd`
cpl=software/EvalSVcallers/scripts/convert_SV_callers_vcf.pl
epl=software/EvalSVcallers/scripts/evaluate_SV_callers.pl

otool=wham
ntool=Wham

for i in `ls $dir/../data/alignment`
do
  perl $cpl \
    -t $ntool \
    $dir/$i/wham${i}.vcf \
    > $dir/$i.vcf

  perl $epl \
    $dir/$i.vcf
done
