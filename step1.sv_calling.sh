#===========================================================#
# Below are the commands for SV detection on 10 datasets using BASIL-ANISE. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`

asgenseng=software/AS-GENSENG
anise_basil=software/anise_basil/build/bin
runlog=$1
mem=20
cpu=10

for sam in `ls ${dir}/../data/alignment`
do
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}

  date
  ${anise_basil}/basil \\
    -ir ${ref} \\
    -im ${bam} \\
    -ov ${dir}/${sam}/basil.vcf \\
    --filter-max-coverage 105 --filter-min-aln-quality 3
  date
  ${anise_basil}/anise \\
    -ir ${ref} \\
    -im ${bam} \\
    -iv ${dir}/${sam}/basil.vcf \\
    -of ${dir}/${sam}/basil.INS.fa \\
    --num-threads 6 --read-mapping-error-rate 2 --overlapper-max-error-rate 3
  date
  " > t1.${sam}.sh
  if [ $runlog == "1" ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
      id=`echo $qsub | cut -d " " -f 3`
      perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using BICseq2. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`

bicseq2=software/BICseq2
mem=10
cpu=10
runlog=$1

for i in {0..10}
do
  sam=${samplelist[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    mkdir -p ${dir}/${sam}
    date
    samtools view ${bam} \\
      | awk '\$6 !~ /H/' \\
      | perl software/samtools-1.3/misc/samtools.pl unique - \\
      | cut -f 4 \\
      > ${dir}/${sam}/chr${chr}.seq
    echo -e \"ChromName\\tfaFile\\tMapFile\\treadPosFile\\tbinFileNorm\\nchr${chr}\\t${ref}\\t${mappability}\\t${dir}/${sam}/chr${chr}.seq\\t${dir}/${sam}/chr${chr}.norm.bin\" \\
      > ${dir}/${sam}/chr${chr}.config
    date
    perl ${bicseq2}/NBICseq-norm_v0.2.4/NBICseq-norm.pl \\
      -l ${readlength[i]} -s ${insertsize[i]} \\
      --tmp ${dir}/${sam}/chr${chr}.tmp \\
      ${dir}/${sam}/chr${chr}.config \\
      ${dir}/${sam}/chr${chr}.norm.out
    rm -fr ${dir}/${sam}/chr${chr}.tmp
    date
    echo -e \"ChromName\\tbinFileNorm\\nchr${chr}\\t${dir}/${sam}/chr${chr}.norm.bin\" \\
      > ${dir}/${sam}/chr${chr}.config.seq.txt
    perl ${bicseq2}/NBICseq-seg_v0.7.2/NBICseq-seg.pl \\
      --bootstrap ${dir}/${sam}/chr${chr}.config.seq.txt \\
      ${dir}/${sam}/chr${chr}.cnv.txt
    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
            qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
            id=`echo $qsub | cut -d " " -f 3`
            perl ${dir}/../usage.pl ${id} > time.${id}.log &
      fi

  done
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using BreakDancer. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`

bicseq2=software/BICseq2

for i in {0..10}
do
  sam=${samplelist[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/ref.chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    mkdir -p ${dir}/${sam}
    date
    perl software/breakdancer/perl/bam2cfg.pl \\
      ${bam} \\
      > ${dir}/${sam}/chr${chr}.config.txt
    date
    software/breakdancer/bin/bin/breakdancer-max \\
      -y 20 -x 200 -r 3 -m 10000000 \\
      ${dir}/${sam}/chr${chr}.config.txt \\
      > ${dir}/${sam}/chr${chr}.cnv.txt
    date
    " > t1.${sam}.chr${chr}.sh

       id=`qsub -clear -cwd -binding linear:1 -l num_proc=1 -l vf=6G t1.${sam}.chr${chr}.sh | awk '{print $3}'`
       perl ${dir}/../usage.pl ${id} > time.${sam}.${chr}.${id}.log &
  done
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using BreakSeek. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
breakseek=software/breakseek/breakseek.py
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  isd=${insertsizesd[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    mkdir -p ${dir}/${sam}
    cd ${dir}/${sam}
    
    date
    samtools sort -n \\
      ${bam} \\
      | samtools view - \\
      > ${dir}/${sam}/chr${chr}.sn.sam
    python ${breakseek} \\
      -f ${dir}/${sam}/chr${chr}.sn.sam \\
      -r ${ref} \\
      -o ${dir}/${sam}/${chr} \\
      -m ${is} -s ${isd} -q ${rl}
    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using BreakSeq2. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

breakseq2=software/anaconda2/bin/run_breakseq2.py
bplib_gff=software/breakseq2/breakseq2_bplib_20150129.hg19.gff
bwa=Software/bwa-0.7.13/bwa
samtools=software/samtools-1.3/bin/samtools
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}

  date
  ${breakseq2} \\
    --bams ${bam} \\
    --sample ${sam} \\
    --work ${dir}/${sam} \\
    --nthreads ${cpu} \\
    --bwa ${bwa} \\
    --samtools ${samtools} \\
    --reference ${ref} \\
    --bplib_gff ${bplib_gff}

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using CNVnator. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

cnvnator=software/CNVnator-0.3.3/cnvnator
mem=13
cpu=10
chr_list='chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX'

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  refdir=hg19/chr_dir
  echo "
  export LD_LIBRARY_PATH=software/pcre-8.41/lib/:\$LD_LIBRARY_PATH
  source software/root-6.12.06/bin/thisroot.sh
  mkdir -p ${dir}/${sam}

  date
  ${cnvnator} \\
    -root ${dir}/${sam}/out.root \\
    -chrom ${chr_list} \\
    -tree ${bam} \\
    -unique

  date
  ${cnvnator} \\
    -root ${dir}/${sam}/out.root \\
    -chrom ${chr_list} \\
    -his ${bin[i]} \\
    -d ${refdir}

  date
  ${cnvnator} \\
    -root ${dir}/${sam}/out.root \\
    -chrom ${chr_list} \\
    -stat ${bin[i]}

  date
  ${cnvnator} \\
    -root ${dir}/${sam}/out.root \\
    -chrom ${chr_list} \\
    -partition ${bin[i]}

  date
  ${cnvnator} \\
    -root ${dir}/${sam}/out.root \\
    -chrom ${chr_list} \\
    -call ${bin[i]} \\
    > ${dir}/${sam}/cnv.txt

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using Control-FREEC. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

freecdir=software/FREEC
freec=${freecdir}/src/freec
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappabilityfile=${freecdir}/data/out100m2_hg19.gem
    snpfile=${freecdir}/data/hg19_snp142.SingleDiNucl.1based.txt
    chrlen=${ref}.fai
    echo "
    mkdir -p ${dir}/${sam}/chr${chr}

    date
    echo -e \"
[general]
chrLenFile=${chrlen}
chrFiles=${dir}/../data/refsplit
gemMappabilityFile=${mappabilityfile}
ploidy=2
step=1000
window=50000
outputDir=${dir}/${sam}/chr${chr}
sex=XX
breakPointType=4
noisyData=FALSE
breakPointThreshold=1.5
maxThreads=${cpu}
samtools=software/samtools-1.3/bin/samtools
[sample]
mateFile=${bam}
inputFormat=BAM
mateOrientation=FR
[BAF]
SNPfile=${snpfile}
makePileup=${freecdir}/data/hg19_snp142.SingleDiNucl.1based.bed
fastaFile=hg19.fa
minimalCoveragePerPosition=5
    \" > ${dir}/${sam}/chr${chr}/chr${chr}.config

    ${freec} -conf ${dir}/${sam}/chr${chr}/chr${chr}.config

    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using DELLY. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
delly=software/delly-0.7.2/src/delly
dellyet=software/delly/excludeTemplates/human.hg19.excl.tsv
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}
  date
  $delly \\
    -g $ref \\
    -o $dir/$sam/del.bcf \\
    -t DEL \\
    $bam \\
    -x $dellyet
  date

  $delly -g $ref -o $dir/$sam/dup.bcf -t DUP $bam -x $dellyet

  $delly -g $ref -o $dir/$sam/inv.bcf -t INV $bam -x $dellyet
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using DINUMT. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
dinumt=software/dinumt-master
mem=20
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}
  date
  perl $dinumt/dinumt.pl \\
    --mask_filename=$dinumt/refNumts.bed \\
    --input_filename=$bam \\
    --reference=$ref \\
    --output_filename=$dir/$sam/dinumt.vcf \\
    --min_reads_cluster=1 \\
    --prefix=NUMT \\
    --len_cluster_include=800 \\
    --len_cluster_link=1600 \\
    --max_read_cov=150 \\
    --ucsc

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using ERDS. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
edrs=software/ERDS/erds_tcag/src/erds_pipeline.pl
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    snv=$dir/../data/snvsplit/$sam.chr$chr.vcf.gz
    echo "
    mkdir -p ${dir}/${sam}
    date
    perl $edrs \\
      -b $bam \\
      -r $ref \\
      -v $snv \\
      -o $dir/$sam/chr$chr \\
      --sd b37 \\
      --large 1000000 \\
      --small 50
    date

    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../zoo.usage.pl ${id} > time.${id}.log &
    fi
  done
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using FermiKit. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

fermi=fermikit/fermi.kit
mem=90
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  fq1=$dir/../data/sequence/$sam/${sam}_1.fq.gz
  fq2=$dir/../data/sequence/$sam/${sam}_2.fq.gz
  echo "
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  $fermi/fermi2.pl unitig -s2g -l$rl -t$cpu -p prefix \\
    \"$fermi/seqtk mergepe $fq1 $fq2 | $fermi/trimadap-mt -p$cpu\" > $dir/$sam/prefix.mak

  date
  make -f $dir/$sam/prefix.mak

  date
  $fermi/run-calling -o $dir/$sam/prefix -t $cpu $ref $dir/$sam/prefix.mag.gz | sh

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using GASVpro. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
gasv=software/gasv
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    mkdir -p ${dir}/${sam}/chr${chr}
    cd $dir/$sam/chr${chr}
    ln -s $bam
    ln -s $gasv

    date
    cp $gasv/bin/GASVPro-HQ.sh $dir/$sam/chr$chr

    date
    sed -i 's/BAMFILE=/BAMFILE=${sam}.chr${chr}.bam/g' $dir/$sam/chr$chr/GASVPro-HQ.sh
    sed -i 's/GASVDIR=/GASVDIR=gasv/g' $dir/$sam/chr$chr/GASVPro-HQ.sh

    date
    sh GASVPro-HQ.sh

    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done

#===========================================================#
# Below are the commands for SV detection on 10 datasets using GenomeSTRiP2. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
svdir=software/svtoolkit
mem=200
cpu=20
ref=sv_genomestrip/fa/hg19.fa
dict=sv_genomestrip/fa/hg19.dict

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  mkdir -p $dir/$sam
  cd $dir/$sam
  echo "
  export SV_DIR=$svdir
  classpath=\"\$SV_DIR/lib/SVToolkit.jar:\$SV_DIR/lib/gatk/GenomeAnalysisTK.jar:\$SV_DIR/lib/gatk/Queue.jar\"
  cd $dir/$sam

  java -jar picard/build/libs/picard.jar AddOrReplaceReadGroups \\
    I=$bam \\
    O=$dir/$sam/$sam.newRG.bam \\
    RGID=$sam \\
    RGLB=$sam \\
    RGPL=COMPLETE \\
    RGSM=$sam \\
    RGPU=unit1

  samtools view -H $dir/$sam/$sam.newRG.bam | grep -v random | grep -v chrUn | grep -v hap > $dir/$sam/$sam.newRG.sam
  samtools view $dir/$sam/$sam.newRG.bam >> $dir/$sam/$sam.newRG.sam
  less $dir/$sam/$sam.newRG.sam | samtools view -bhS - > $dir/$sam/$sam.newRG.bam
  rm $dir/$sam/$sam.newRG.sam
  samtools index $dir/$sam/$sam.newRG.bam

  echo -e \"$sam\\tFemale\" > $dir/$sam/gendermap.txt

  date
  java -Xmx${mem}g -Xms${mem}g \\
    -cp \${classpath} \\
    org.broadinstitute.gatk.queue.QCommandLine \\
    -S \${SV_DIR}/qscript/SVPreprocess.q \\
    -S \${SV_DIR}/qscript/SVQScript.q \\
    -cp \${classpath} \\
    -gatk \${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \\
    -configFile \${SV_DIR}/conf/genstrip_parameters.txt \\
    -R $ref \\
    -I $dir/$sam/$sam.newRG.bam \\
    -md $dir/$sam/out_meta \\
    -ploidyMapFile \$SV_DIR/hg19_ploidy.map \\
    -jobLogDir $dir/$sam/logDir \\
    -jobRunner ParallelShell \\
    -gatkJobRunner ParallelShell \\
    -run

  date
  java -Xmx${mem}g -Xms${mem}g \\
    -cp \${classpath} \\
    org.broadinstitute.gatk.queue.QCommandLine \\
    -S \${SV_DIR}/qscript/SVDiscovery.q \\
    -S \${SV_DIR}/qscript/SVQScript.q \\
    -cp \${classpath} \\
    -gatk \${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \\
    -configFile \${SV_DIR}/conf/genstrip_parameters.txt \\
    -R $ref \\
    -I $dir/$sam/$sam.newRG.bam \\
    -md $dir/$sam/out_meta \\
    -runDirectory $dir/$sam/discovery \\
    -jobLogDir $dir/$sam/logs \\
    -O $dir/$sam/svdiscovery.dels.vcf \\
    -genderMapFile $dir/$sam/gendermap.txt \\
    -jobRunner ParallelShell  \\
    -gatkJobRunner ParallelShell  \\
    -minimumSize 30 \\
    -maximumSize 2000000 \\
    -debug true \\
    -run
  date

  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi

  cd ..

done



#===========================================================#
# Below are the commands for SV detection on 10 datasets using GRIDSS. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

gridss=software/gridss-2.10.2
mem=30
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  bed=$gridss/wgEncodeDacMapabilityConsensusExcludable.bed
  echo "
  mkdir -p ${dir}/${sam}/ref
  for i in hg19.dict hg19.fa hg19.fa.2bit hg19.fa.amb hg19.fa.ann hg19.fa.bwt hg19.fa.dict hg19.fa.fai hg19.fa.img hg19.fa.pac hg19.fa.sa
  do
    ln -s hg19/\$i $dir/$sam/ref
  done

  date
  sh $gridss/gridss.sh \\
    -r $dir/$sam/ref/hg19.fa \\
    -o $dir/$sam/$sam.sv.vcf \\
    -a $dir/$sam/$sam.gridss.assembly.bam \\
    -j $gridss/gridss-2.10.2-gridss-jar-with-dependencies.jar \\
    -w $dir/$sam \\
    -b $bed \\
    $bam

  date
  Rscript $dir/simple-event-annotation.R \\
    $dir/$sam/$sam.sv.vcf \\
    $dir/$sam/$sam.sv.anno.vcf

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using Hydra-sv. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

hydra=software/Hydra
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  export PATH=software/anaconda2/bin:$hydra/bin:$hydra/scripts:\$PATH
  set ulimit -f 16384
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  echo -e \"$sam\\t$bam\" > $dir/$sam/config.stub.txt

  sh $hydra/hydra-multi.sh run -t $cpu -p 100 -o $sam $dir/$sam/config.stub.txt

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using iCopyDAV. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
icopydav=software/icopydav-master
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    map=$icopydav/mapp/chr$chr.dat
    gc=$icopydav/gc/chr$chr.gc
    chrl=$icopydav/chrlen/chr$chr.gen
    echo "
    export LD_LIBRARY_PATH=software/openmpi-4.1.0/lib:\$LD_LIBRARY_PATH
    mkdir -p ${dir}/${sam}/chr${chr}
    cd $dir/$sam/chr$chr
    
    chrlen=`cat $ref.fai | awk '{print \$2}'`
    echo -e \"minSize=100\\ngenomeSize=\$chrlen\\npercCNLoss=0.05\\npercCNGain=0.05\\nfdr=0.01\\noverDispersion=3\\nploidy=2\" > chr$chr.config.txt
    ln -s $bam
    ln -s $bam.bai

    date
    $icopydav/calOptBinSize -c chr$chr.config.txt -i ${sam}.chr${chr}.bam

    date
    $icopydav/prepareData -m $map -g $gc --genome_file $chrl -o chr$chr --win 1000

    date
    $icopydav/pretreatment -i ${sam}.chr${chr}.bam -o chr$chr -z chr${chr}_1000.bin --mapfile chr${chr}_1000.map --gcfile chr${chr}_1000.gc

    date
    $icopydav/runSegmentation -o chr${chr} -t

    date
    $icopydav/callCNV -o chr${chr} -z chr${chr}_1000.bin --hg19

    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done

#===========================================================#
# Below are the commands for SV detection on 10 datasets using IndelMINER. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
indelminer=software/indelMINER-0.2/src
mem=100
cpu=5

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  sd=${insertsizesd[i]}
  let min=${is}-$sd*3
  let max=${is}+$sd*3
  dep=${sampledepthlist[i]}
  rg=${readgroup[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    mkdir -p ${dir}/${sam}/chr${chr}
    cd ${dir}/${sam}/chr${chr}

    date
    cp $dir/is/$sam.cf indelminer.config
    echo -e \"RC\\tchr$chr\\t$dep\" >> indelminer.config

    samtools view -H $bam > $dir/$sam/chr${chr}.sam
    samtools view $bam | awk '\$6 !~ /[NHP]/' >> $dir/$sam/chr${chr}.sam
    samtools view -bhS $dir/$sam/chr${chr}.sam -o $dir/$sam/chr${chr}.bam
    samtools index $dir/$sam/chr${chr}.bam

    $indelminer/indelminer \\
      $ref \\
      -i $dir/$sam/chr$chr/indelminer.config \\
      sample=$dir/$sam/chr${chr}.bam \\
      -e 3 \\
      -c chr$chr \\
      > $dir/$sam/chr${chr}/indelminer.chr${chr}.flt.vcf

    date

    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done

#===========================================================#
# Below are the commands for SV detection on 10 datasets using ITIS. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
itis=software/ITIS
mem=20
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  fq1=$dir/../data/sequence/$sam/${sam}_1.fq.gz
  fq2=$dir/../data/sequence/$sam/${sam}_2.fq.gz
  for t in ALU LINE1 SVA HERVK
  do
    echo "
    export PATH=software/samtools-0.1.19:\$PATH
    mkdir -p ${dir}/${sam}/$t
    cd ${dir}/${sam}/$t

    date
    $itis/itis.pl \\
      -g $ref \\
      -t $itis/ref/$t.fa \\
      -N $t \\
      -l $is -e Y -c 10,3,3 \\
      -1 $fq1 \\
      -2 $fq2 \\
      -D $dir/$sam/$t

    date
    " > t1.${sam}.$t.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.$t.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
	done
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using laSV. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
lasv=software/laSV
mem=200
cpu=20

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  echo "
  export PATH=software/samtools-0.1.19:\$PATH

  mkdir -p ${dir}/${sam}
  cd $dir/$sam
  
  date
  $lasv/scripts/run_laSV.sh \\
    -I $dir/../data/sequence/$sam/$sam \\
    -D $lasv \\
    -R $lasv/ref \\
    -G hg19 \\
    -f $is \\
    -l $rl \\
    -k 31 \\
    -t $cpu \\
    -o $dir/$sam

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi

  cd ..
done

#===========================================================#
# Below are the commands for SV detection on 10 datasets using Lumpy. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
lumpy=software/lumpy-sv
svtyper=software/anaconda2/bin/svtyper
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  sd=${insertsizesd[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  samtools view -uF 0x0002 $bam \\
   | samtools view -uF 0x100  - \\
   | samtools view -uF 0x0004 - \\
   | samtools view -uF 0x0008 - \\
   | samtools view -bF 0x0400 - \\
   | samtools sort -o $dir/$sam/$sam.discordant.sort.bam -
  samtools index $dir/$sam/$sam.discordant.sort.bam

  samtools view -h $bam \\
   | $lumpy/scripts/extractSplitReads_BwaMem -i stdin \\
   | samtools view -bS - \\
   | samtools sort -o $dir/$sam/$sam.sr.sort.bam -
  samtools index $dir/$sam/$sam.sr.sort.bam

  samtools view $bam \\
   | tail -n+100000 \\
   | python $lumpy/scripts/pairend_distro.py \\
     -r 100 -X 4 -N 10000 \\
     -o $dir/$sam/$sam.histo

  date
  $lumpy/bin/lumpy \\
    -mw 4 \\
    -tt 0.0 \\
    -pe bam_file:$dir/$sam/$sam.discordant.sort.bam,histo_file:$dir/$sam/$sam.histo,mean:$is,stdev:$sd,read_length:$rl,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,min_mapping_threshold:20,id:1 \\
    -sr bam_file:$dir/$sam/$sam.sr.sort.bam,back_distance:20,weight:1,min_mapping_threshold:20,id:2 \\
    > $dir/$sam/$sam.vcf

  $svtyper \\
    -i $dir/$sam/$sam.vcf \\
    -B $bam \\
    -l $dir/$sam/bam.json \\
    > $dir/$sam/$sam.gt.vcf
  date

  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done



#===========================================================#
# Below are the commands for SV detection on 10 datasets using Manta. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
manta=software/manta-1.6.0/build/bin
manta=software/manta-1.6.0.centos6_x86_64/bin
mem=200
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  #export PATH=software/anaconda2/bin:\$PATH
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  software/anaconda2/bin/python \\
    $manta/configManta.py \\
    --bam $bam \\
    --referenceFasta $ref \\
    --runDir $dir/$sam

  date
  software/anaconda2/bin/python \\
    $dir/$sam/runWorkflow.py \\
    -m local -j $cpu -g $mem

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using MATCHCLIP. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
matchclip=software/matchclips2
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  export PATH=$matchclip:\$PATH
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  matchclips -b $bam -f $ref -t $cpu -o $dir/$sam/$sam.out

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using Meerkat. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

samclip=software/samclip/samclip
meerkat=software/Meerkat/scripts
mem=20
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment_bwaaln/${sam}/${sam}.sort.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  echo "
  export LD_LIBRARY_PATH=software/Meerkat/src/mybamtools/lib/:\$LD_LIBRARY_PATH
  export PATH=software/BLAST/blast-2.2.26/bin:\$PATH
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  ln -s $bam $sam.bam
  ln -s $bam.bai $sam.bam.bai

  for f in directory.index.dir directory.index.pag hg19.dict hg19.fa hg19.fa.2bit hg19.fa.amb hg19.fa.ann hg19.fa.bwt hg19.fa.dict hg19.fa.fai hg19.fa.gridsscache hg19.fa.img hg19.fa.index hg19.fa.pac hg19.fa.ploidy hg19.fa.sa rmsk_hg19.txt
  do
    cp hg19/\$f ./
  done

  date
  $meerkat/pre_process.pl \\
    -b $sam.bam \\
    -I hg19.fa \\
    -A hg19.fa.fai \\
    -k 200 -f 0 -l 0 -t $cpu \\
    -W software/fermikit/fermi.kit \\
    -S software/samtools-0.1.19 

  date
  perl $meerkat/meerkat.pl \\
    -b $sam.bam \\
    -u 1 -a 0 -p 3 -o 1 -q 2 -z 1000000 -l 0 -t $cpu \\
    -F $dir/$sam \\
    -B software/BLAST/blast-2.2.26/bin \\
    -W software/fermikit/fermi.kit \\
    -S software/samtools-0.1.19
 
  date
  perl $meerkat/mechanism.pl \\
    -b $sam.bam \\
    -R rmsk_hg19.txt

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi

  cd ..
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using MELT. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
melt=software/MELTv2.2.2/MELT.jar
db=software/MELTv2.2.2/me_refs/hg19_ucsc
genebed=software/MELTv2.2.2/add_bed_files/1KGP_Hg19/hg19.genes.bed
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}
  date
  java -Xmx${mem}G -jar $melt Preprocess -bamfile $bam -h $ref

  date
  for type in ALU SVA LINE1 HERVK chrM
  do
    java -Xmx${mem}G -jar $melt Single \\
    -bamfile $bam -h $ref -n $genebed -t $db/\${type}_MELT.zip \\
    -w $dir/$sam -r $rl -e $is -d 40000000 -c ${sampledepthlist[i]}
  done

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using MetaSV. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

metasv=software/metasv-0.5.2
gapbed=$metasv/metasv/resources/hg19.gaps.bed
spadesdir=software/SPAdes-3.15.3-Linux/bin/spades.py

mem=200
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  isd=${insertsizesd[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  pindel=$dir/../Pindel/$sam
  bk=$dir/../BreakDancer/$sam/breakdancer.all.txt
  bkseq=$dir/../BreakSeq2/$sam/breakseq.gff
  cnvnator=$dir/../CNVnator/$sam/cnv.txt
  echo "
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date

  $metasv/build/scripts-2.7/run_metasv.py \\
    --reference $ref \\
    --bam $bam \\
    --filter_gaps --gaps $gapbed \\
    --outdir $dir/$sam \\
    --sample $sam \\
    --num_threads $cpu \\
    --isize_mean $is \\
    --isize_sd $isd \\
    --min_support_ins 3 \\
    --disable_assembly \\
    --spades $spadesdir \\
    --pindel_native $pindel/pindel_D $pindel/pindel_TD $pindel/pindel_INV $pindel/pindel_SI \\
    --breakdancer_native $bk \\
    --breakseq_native $bkseq \\
    --cnvnator_native $cnvnator

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using MindTheGap. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

mindthegap=software/MindTheGap-v2.2.2-bin-Linux/bin/MindTheGap
mem=40
cpu=10
maxm=40000

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  fq1=$dir/../data/sequence/$sam/${sam}_1.fq.gz
  fq2=$dir/../data/sequence/$sam/${sam}_2.fq.gz
  echo "
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  $mindthegap find -in $fq1,$fq2 -ref $ref -out $dir/$sam/s1.find -nb-cores $cpu -max-memory $maxm

  date
  $mindthegap fill -in $fq1,$fq2 -bkpt $dir/$sam/s1.find.breakpoints -out $dir/$sam/s2.find -nb-cores $cpu -max-memory $maxm

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using PennCNV-Seq. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

penn1=software/PennCNV-1.0.5
penns=software/PennCNV-Seq
mem=30
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    export PATH=$penn1:$penns:software/perl-5.14.2/bin:\$PATH
    mkdir -p ${dir}/${sam}/$chr
    cd $dir/$sam/$chr

    date
    penncnv-seq_example.sh $penn1 $penns/reference hg19 EAS $ref $bam $chr

    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done

#===========================================================#
# Below are the commands for SV detection on 10 datasets using Pindel. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

pindel=software/pindel

mem=30
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  echo -e \"$bam\\t$is\\t$sam\" > pindel.cfg

  $pindel/pindel \\
    -f $ref \\
    -i pindel.cfg \\
    -o $sam \\
    -Y $ref.ploidy \\
    -c All -x 2 -M 3 -v 100 -d 50 -E 0.92 -w 10 -t $cpu --MIN_DD_MAP_DISTANCE 3000 -g

 $pindel/pindel2vcf \\
   -f $ref \\
   -P $sam \\
   -R hg19 \\
   -d 20210702 \\
   -v $sam.vcf \\
   --min_size 50 \\
   --only_balanced_samples

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using Popins. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

popins=software/popins/popins
mem=10
cpu=1

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    mkdir -p ${dir}/${sam}/$chr
    cd $dir/$sam/$chr
    ln -s $ref     genome.fa
    ln -s $ref.fai genome.fa.fai

    date
    $popins assemble -t $cpu -s $sam.$chr $bam
    $popins merge -p ./
    $popins contigmap $sam.$chr
    $popins place-refalign   --readLength $rl
    $popins place-splitalign --readLength $rl $sam.$chr
    $popins place-finish

    less insertions.vcf | awk 'BEGIN{OFS=\"\\t\"}{if(\$1 == \"#CHROM\") \$9 = \"FORMAT\"; if(\$1 !~ /#/ ) \$9 =\".\"; print \$0;}' > insertions_edit.vcf
    mv insertions_edit.vcf insertions.vcf

    $popins genotype $sam.$chr

    less $sam.$chr/insertions.vcf | grep ^# > $sam.$chr/i.vcf
    less $sam.$chr/insertions.vcf | grep -v ^# | awk 'BEGIN{OFS=\"\\t\"}{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$11}' >> $sam.$chr/i.vcf
    mv $sam.$chr/i.vcf $sam.$chr/insertions.vcf

    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using PRISM. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
prism=software/PRISM_1_1_6
mem=5
cpu=1

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    export PRISM_PATH=$prism
    mkdir -p ${dir}/${sam}/$chr
    cd $dir/$sam/$chr

    date
    samtools sort -n $bam | samtools view - -o $sam.$chr.sn.sam

    $prism/toolkit/run_PRISM.sh -m $is -e ${insertsizesd[i]} -p 3 -l $rl -r $ref -i $sam.$chr.sn.sam

    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using RetroSeq. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
te=software/RetroSeq-master/data/telist
retro=software/RetroSeq-master/bin/retroseq.pl

mem=50
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  export PATH=software/samtools-0.1.19:\$PATH
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  $retro -discover -bam $bam -refTEs $te -output retroseq.out -q 20

  date
  $retro -call -bam $bam -input retroseq.out -ref $ref -output retroseq.vcf -hets -filter $te -reads 3 -q 20

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done



#===========================================================#
# Below are the commands for SV detection on 10 datasets using Socrates. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bb=miniconda3/bin/bowtie2-build
jar=software/socrates-1.13.1/socrates-1.13.1-jar-with-dependencies.jar
mem=50
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  ln -s $ref 
  $bb hg19.fa hg19.fa

  date
  java -Xms${mem}G -Xmx${mem}G -jar $jar -p 97 -t $cpu hg19.fa $bam 

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using SoftSV. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
softsv=software/SoftSV_1.4.2/SoftSV
mem=10
cpu=1

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    mkdir -p ${dir}/${sam}/$chr
    cd $dir/$sam/$chr

    date
    $softsv -i $bam -r chr$chr -o $dir/$sam/$chr

    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using Sprites. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
sprites=software/sprites/build/sprites

mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  isd=${insertsizesd[i]}
  rl=${readlength[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../data/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../data/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    mkdir -p ${dir}/${sam}/$chr
    cd $dir/$sam/$chr

    date
    $sprites \\
     -r $ref \\
     -o $dir/$sam/$chr/out \\
     $bam

    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done

#===========================================================#
# Below are the commands for SV detection on 10 datasets using SvABA. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
svaba=software/svaba/bin/svaba
mem=30
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  $svaba \\
    run \\
    -t $bam \\
    -G $ref \\
    -a $sam \\
    -p $cpu \\
    --germline

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using SVelter. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

py=software/miniconda3/bin/python
svdir=software/svelter
svelter=$svdir/build/scripts-3.7/svelter.py
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  export PATH=software/anaconda3/bin/:\$PATH
  mkdir -p ${dir}/${sam}
  cd $dir/$sam
  ln -s $bam $sam.bam
  ln -s $bam.bai $sam.bam.bai

  #date
  $svelter Setup --reference $ref --workdir ./ --support $svdir/Support/hg19 --ref-index $svdir/Support/ref-index/hg19

  date
  $svelter NullModel --sample $sam.bam --workdir ./

  date
  $svelter BPSearch --sample $sam.bam --workdir ./

  for i in \`ls BreakPoints.$sam.bam/ | grep -E '_hap|_gl'\`;do rm BreakPoints.$sam.bam/\$i;done
  less BreakPoints.$sam.bam/$sam.SPCff6.CluCff12.AlignCff0.2.chromLNs |grep -v '_gl' | grep -v '_hap' | grep -v 'chrY' |grep -v 'chrM' > a
  mv a BreakPoints.$sam.bam/$sam.SPCff6.CluCff12.AlignCff0.2.chromLNs

  date
  $svelter BPIntegrate --sample $sam.bam --workdir ./

  date
  $svelter SVPredict --sample $sam.bam --workdir ./ --bp-file bp_files.$sam.bam/$sam.txt

  date
  $svelter SVIntegrate --workdir ./ --prefix $sam --input-path bp_files.$sam.bam

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using TEMP. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

temp=software/TEMP
mem=30
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  isd=${insertsizesd[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment_bwaaln/${sam}/${sam}.sort.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  ln -s $bam $sam.sorted.bam
  ln -s $bam.bai $sam.sorted.bam.bai

  date
  sh $temp/scripts/TEMP_Insertion.sh \\
    -i $sam.sorted.bam \\
    -s $temp/scripts \\
    -m 3 \\
    -r $temp/database/te.fa \\
    -t $temp/database/hg19_rmsk.bed \\
    -c $cpu

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using TIDDIT. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

bicseq2=software/BICseq2
ti=software/TIDDIT
python=miniconda3/bin/python
tipy=software/TIDDIT/TIDDIT.py
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  isd=${insertsizesd[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  export LD_LIBRARY_PATH=$ti/build/lib/bamtools/src/api/:\$LD_LIBRARY_PATH
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  $python $tipy --sv --bam $bam -o $sam --ref $ref
  #$ti/bin/TIDDIT --sv -b $bam -o $sam

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


#===========================================================#
# Below are the commands for SV detection on 10 datasets using Ulysses. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

uly=software/ulysses-ulysses-v1.0
mem=20
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  isd=${insertsizesd[i]}
  rl=${readlength[i]}
  for chr in {1..22} X
  do
    bam=${dir}/../20210408.input_split/bamsplit/${sam}.chr${chr}.bam
    ref=${dir}/../20210408.input_split/refsplit/chr${chr}.fa
    mappability=${bicseq2}/nonNregion/hg19.nonN.chr${chr}
    echo "
    mkdir -p ${dir}/${sam}/$chr
    cd $dir/$sam/$chr
    ln -s $bam
    ln -s $bam.bai

    date
    $uly/ReadBAM.py -out uly.$chr -stats $sam.$chr.stats.txt -p uly_conf.$chr $sam.chr$chr.bam

    date
    $uly/Ulysses.py -p uly_conf.$chr -typesv DEL

    $uly/Ulysses.py -p uly_conf.$chr -typesv DUP
    
    $uly/Ulysses.py -p uly_conf.$chr -typesv INV
    date
    " > t1.${sam}.chr${chr}.sh

    if [ $runlog -eq 1 ];then
      qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.chr${chr}.sh`
        id=`echo $qsub | cut -d " " -f 3`
        perl ${dir}/../usage.pl ${id} > time.${id}.log &
    fi
  done
done

#===========================================================#
# Below are the commands for SV detection on 10 datasets using Wham. 
# The commands sequentially iterate through the 10 datasets, 
#  construct specific scripts for each dataset, 
#  and finally execute the scripts while recording the time and storage consumption of the analysis.
#===========================================================#

dir=`pwd`
runlog=$1

wham=software/wham/bin/whamg
mem=10
cpu=10

for i in {0..10}
do
  sam=${samplelist[i]}
  is=${insertsize[i]}
  isd=${insertsizesd[i]}
  rl=${readlength[i]}
  bam=${dir}/../data/alignment/${sam}/${sam}.sorted.bam
  vcf=${dir}/../data/snvindel/${sam}.vcf.gz
  ref=hg19/hg19.fa
  echo "
  mkdir -p ${dir}/${sam}
  cd $dir/$sam

  date
  $wham -f $bam -a $ref -x $cpu -c  chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX > wham$sam.vcf
  #$wham -f $bam -a $ref -x $cpu -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y > wham$sam.vcf

  date
  " > t1.${sam}.sh

  if [ $runlog == "1" ];then
    qsub=`qsub -clear -cwd -binding linear:${cpu} -l num_proc=${cpu} -l vf=${mem}G t1.${sam}.sh`
    id=`echo $qsub | cut -d " " -f 3`
    perl ${dir}/../usage.pl ${id} > time.${id}.log &
  fi
done


