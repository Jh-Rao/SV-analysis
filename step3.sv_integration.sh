dir=`pwd`

ref=hg19/hg19.fa

mkdir -p $dir/svset_soft
#===========================================================#
# Construct SV sets for different software tools. 
# Select a software, traverse the SVs detected by the software across all DNBSEQ datasets, 
#  and extract the SVs consistently identified in all DNBSEQ datasets as the software-specific SV set.
#===========================================================#
pl=DNBSEQ
length=${#dnblist[@]}
for soft in ${softlist[@]}
do
  echo --- --- ---- --- $soft
  for svt in DEL DUP INS INV
  do
    echo --- --- $svt
    let pre=length-1
    let star=pre-1
    less $dir/RSS3_bed/$soft.${dnblist[pre]}.bed | awk -v P=${dnblist[pre]} 'BEGIN{OFS="\t"}{print $0,P}' > $dir/tmp.S1.bed
    sed -i 's/MEI/INS/g' $dir/tmp.S1.bed
    for ((i=$star; i>=0; i--))
    do
      if [[ $soft == "BreakSeek" ]] && [[ ${dnblist[i]} == "MGISEQ-2000_PE100" ]];then
        echo pass
      elif [[ $soft == "Popins" ]] && [[ ${dnblist[i]} == "MGISEQ-2000_PE100" ]];then
        echo pass
      elif [[ $soft == "BreakSeek" ]] && [[ ${dnblist[i]} == "BGISEQ-500_PE100" ]];then
        echo pass
      elif [[ $soft == "Popins" ]] && [[ ${dnblist[i]} == "BGISEQ-500_PE100" ]];then
        echo pass
      else
		echo --- ${dnblist[i]}
		cp $dir/RSS3_bed/$soft.${dnblist[i]}.bed $dir/tmp.S2.bed
		sed -i 's/MEI/INS/g' $dir/tmp.S2.bed
		perl $dir/merge_2in1.pl $dir/tmp.S1.bed $dir/tmp.S2.bed ${dnblist[i]} $svt n $dir/tmp.S3.bed
		mv $dir/tmp.S3.bed $dir/tmp.S1.bed
      fi
    done
    cp $dir/tmp.S1.bed $dir/tmp.$soft.$svt.bed
  done
  cat $dir/tmp.$soft.*.bed | sort -k 1,1 -k 2,2n > $dir/svset_soft/$soft.$pl.bed
  rm -f $dir/tmp.*
done

#===========================================================#
# Construct SV sets for different software tools.
# Select a software, traverse the SVs detected by the software across all Illumina datasets,
#  and extract the SVs consistently identified in all Illumina datasets as the software-specific SV set.
#===========================================================#
pl=Illumina
for soft in ${softlist[@]}
do
  echo --- --- ---- --- $soft
  for svt in DEL DUP INS INV
  do
    echo --- --- $svt
    less $dir/RSS3_bed/$soft.${illlist[0]}.bed | awk -v P=${illlist[0]} 'BEGIN{OFS="\t"}{print $0,P}' > $dir/tmp.S1.bed
    sed -i 's/MEI/INS/g' $dir/tmp.S1.bed
    for ((i=1; i<${#illlist[@]}; i++))
    do
      echo --- ${illlist[i]}
      cp $dir/RSS3_bed/$soft.${illlist[i]}.bed $dir/tmp.S2.bed
      sed -i 's/MEI/INS/g' $dir/tmp.S2.bed
      perl $dir/merge_2in1.pl $dir/tmp.S1.bed $dir/tmp.S2.bed ${illlist[i]} $svt n $dir/tmp.S3.bed
      mv $dir/tmp.S3.bed $dir/tmp.S1.bed
    done
    cp $dir/tmp.S1.bed $dir/tmp.$soft.$svt.bed
  done
  cat $dir/tmp.$soft.*.bed | sort -k 1,1 -k 2,2n > $dir/svset_soft/$soft.$pl.bed
  rm -f $dir/tmp.*
done

#===========================================================#
# Construct platform-specific SV sets.
# Traverse the software-specific SV sets across all software tools for the DNBSEQ or Illumina platforms,
#  and select the SVs detected by two or more software tools as the platform-specific SV set.
#===========================================================#
mkdir -p $dir/svset_pl
for pl in DNBSEQ Illumina
do
  echo --- --- --- $pl
  less $dir/svset_soft/${softlist[0]}.$pl.bed |awk -v P=${softlist[0]} 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,P}' > $dir/tmp.F1.$pl.bed
  for ((i=1; i<${#softlist[@]}; i++))
  do
    echo --- --- ${softlist[i]}
    for svt in DEL DUP INS INV
    do
      echo --- $svt
      less $dir/svset_soft/${softlist[i]}.$pl.bed |awk -v P=${softlist[i]} 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,P}' > $dir/tmp.F2.$pl.bed
      perl $dir/merge_2in1.pl $dir/tmp.F1.$pl.bed $dir/tmp.F2.$pl.bed ${softlist[i]} $svt u $dir/tmp.F3.$pl.$svt.bed
    done
    cat $dir/tmp.F3.$pl.*.bed | sort -k 1,1 -k 2,2n > $dir/tmp.F1.$pl.bed
    rm -f $dir/tmp.F3.$pl.*.bed $dir/tmp.F2.$pl.bed
  done
  cp $dir/tmp.F1.$pl.bed $dir/svset_pl/$pl.bed
  less $dir/svset_pl/$pl.bed | grep ':' |awk -v P=$pl 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,P}' > $dir/svset_pl/$pl.2soft.bed
done

for pl in DNBSEQ Illumina
do
  for svt in DEL DUP INS INV
  do
    less $dir/svset_pl/$pl.2soft.bed | grep $svt > $dir/svset_pl/split/$pl.$svt.bed
  done
done

