#/usr/bin/bash
##### Whole Exome Sequence Variant Calling Pipeline #####
echo "Welcome to Exome variant calling Pipeline"
NUMCPUS=30
SOFTWAREDIR="/home/cml/Downloads/gatk-4.2.6.1"
BASEDIRINDEX="/home/cml/humandb/index-GRCh38.p13/HG38-chr"
DATABASEDIR="/home/cml/Downloads/snpEff/data/GRCh38.p13"
BEDINDEX="$BASEDIRINDEX/Twist_ComprehensiveExome_targets_hg38.bed"
GENOMEIDX1="$BASEDIRINDEX/human_ref38_patch13.fna"
GENOMEIDX2="$DATABASEDIR/dbSNP155_GRCh38.chr.vcf"
BASEDIRDATA="$PWD"
reads1=($BASEDIRDATA/*1P.fastq.gz)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1P.fastq.gz/_2P.fastq.gz}")
if [ $reads1 -a $reads2 ]
then
echo "***Following SAMPLES are submitted for Germline variant calling***"
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    sample="${sample%_*}"
    echo "$sample"
    done
else
echo "ERROR: $reads1 $reads2 program not found, please edit the configuration script."
  exit 1
fi
### GATK haplotype caller ###
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    sample="${sample%_*}"
    stime=`date +"%Y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $sample"

    echo [$stime] "   * Alignment of reads to genome (bwa)"
    bwa mem -t $NUMCPUS ${GENOMEIDX1}   \
               ${BASEDIRDATA}/${reads1[$i]} \
              ${BASEDIRDATA}/${reads2[$i]} \
                -o ${BASEDIRDATA}/${sample}_bwa.sam
                
    echo [$stime] "   * Convertion of SAM to BAM and sorting"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" SortSam \
               --INPUT ${BASEDIRDATA}/${sample}_bwa.sam     \
               --OUTPUT ${BASEDIRDATA}/${sample}_bwa_sort.bam --SORT_ORDER coordinate
               
    echo [$stime] "   * Remove SAM"
    rm ${BASEDIRDATA}/${sample}_bwa.sam
    
    echo [$stime] "   * BuildBamIndex"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" BuildBamIndex --INPUT ${BASEDIRDATA}/${sample}_bwa_sort.bam
    
     $SOFTWAREDIR/gatk --java-options "-Xmx200G" CollectAlignmentSummaryMetrics --REFERENCE_SEQUENCE ${GENOMEIDX1} --INPUT ${BASEDIRDATA}/${sample}_bwa_sort.bam --OUTPUT ${BASEDIRDATA}/${sample}_alignment_metrics.txt --VALIDATION_STRINGENCY LENIENT
    
    echo [$stime] "   * AddOrReplaceReadGroups : first"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" AddOrReplaceReadGroups --INPUT ${BASEDIRDATA}/${sample}_bwa_sort.bam --OUTPUT ${BASEDIRDATA}/${sample}_rg.bam --RGLB lib1 --RGPL illumina --RGPU NONE --RGSM BG
    
    echo [$stime] "   * MarkDuplicates"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" MarkDuplicates --INPUT ${BASEDIRDATA}/${sample}_rg.bam --OUTPUT ${BASEDIRDATA}/${sample}_md.bam --METRICS_FILE metrics.txt -AS true --VALIDATION_STRINGENCY LENIENT
    
    echo [$stime] "   * Remove rg.bam"
    rm $BASEDIRDATA/${sample}_rg.bam
    rm $BASEDIRDATA/${sample}_rg.bai
    
    echo [$stime] "   * BuildBamIndex"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" BuildBamIndex --INPUT ${BASEDIRDATA}/${sample}_md.bam
    
    echo [$stime] "   * BaseRecalibrator"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" BaseRecalibrator -R ${GENOMEIDX1} -I ${BASEDIRDATA}/${sample}_md.bam --known-sites ${GENOMEIDX2} -O ${BASEDIRDATA}/${sample}_recal_data.table

    echo [$stime] "   * ApplyBQSR"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" ApplyBQSR -R ${GENOMEIDX1} -I ${BASEDIRDATA}/${sample}_md.bam -bqsr ${BASEDIRDATA}/${sample}_recal_data.table -O ${BASEDIRDATA}/${sample}_rr.bam

    echo [$stime] "   * Remove md.bam"
    rm ${BASEDIRDATA}/${sample}_md.bam
    rm ${BASEDIRDATA}/${sample}_md.bai
done

for file in {*N*_rr.bam,*B*_rr.bam; do
    echo "$file"
    sample1=${file}
    echo $sample1
    sample=${sample1/_rr.bam/}
    echo $sample
    stime=`date +"%Y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $sample"
    echo [$stime] "   * AddOrReplaceReadGroups2"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" AddOrReplaceReadGroups --INPUT $BASEDIRDATA/${sample}_rr.bam --OUTPUT $BASEDIRDATA/${sample}_sm.bam --RGLB lib1 --RGPL illumina --RGPU NONE --RGSM NORMAL_"${sample}"
    
    echo [$stime] "   * BuildBamIndex"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" BuildBamIndex --INPUT $BASEDIRDATA/${sample}_sm.bam
done

for file in {*T*_rr.bam,*P*_rr.bam}; do
    echo "$file"
    sample1=${file}
    echo $sample1
    sample=${sample1/_rr.bam/}
    echo $sample
    stime=`date +"%Y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $sample"
    echo [$stime] "   * AddOrReplaceReadGroups2"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" AddOrReplaceReadGroups --INPUT $BASEDIRDATA/${sample}_rr.bam --OUTPUT $BASEDIRDATA/${sample}_sm.bam --RGLB lib1 --RGPL illumina --RGPU NONE --RGSM TUMOR_"${sample}"
    
    echo [$stime] "   * BuildBamIndex"
    $SOFTWAREDIR/gatk --java-options "-Xmx200G" BuildBamIndex --INPUT $BASEDIRDATA/${sample}_sm.bam
    	
    echo [$stime] "   * MoveBamIndex"
done
