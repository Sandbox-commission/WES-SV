#/usr/bin/bash
##### Whole Exome Sequence Variant Calling Pipeline #####
echo "Welcome to Exome variant calling Pipeline : Part 2"
NUMCPUS=30
SOFTWAREDIR="/home/cml/Downloads/gatk-4.2.6.1"
BASEDIRINDEX="/home/cml/humandb/index-GRCh38.p13/HG38-chr"
DATABASEDIR="/home/cml/Downloads/snpEff/data/GRCh38.p13"
BEDINDEX="$BASEDIRINDEX/Twist_ComprehensiveExome_targets_hg38.bed"
GENOMEIDX1="$BASEDIRINDEX/human_ref38_patch13.fna"
GENOMEIDX2="$DATABASEDIR/dbSNP155_GRCh38.chr.vcf"
GENOMEIDX3="$BASEDIRINDEX/Twist_ComprehensiveExome_targets_hg38.interval_list"
GENOMEIDX4="$BASEDIRINDEX/somatic-hg38_af-only-gnomad.hg38.chr.vcf"
BASEDIRDATA="$PWD"
for file in *N*_sm.bam; do
    echo "$file"
    sample1=${file}
    echo $sample1
    sample=${sample1/_sm.bam/}
    echo $sample
    stime=`date +"%Y-%m-%d %H:%M:%S"`
    echo "[$stime] Processing sample: $sample"
        
    echo [$stime] "   * Generate VCF for Normal Sample: $sample"
    $SOFTWAREDIR/gatk --java-options "-Xmx24G" Mutect2 -R ${GENOMEIDX1} -I $BASEDIRDATA/${sample}_sm.bam --max-mnp-distance 0 -O $BASEDIRDATA/${sample}_sm.vcf.gz
done    


echo [$stime] "   *  Creation Panel of Normals"
$SOFTWAREDIR/gatk --java-options "-Xmx24G" GenomicsDBImport -R $GENOMEIDX1 -L $GENOMEIDX3 --genomicsdb-workspace-path pon_db $(ls $BASEDIRDATA/*_sm.vcf.gz | xargs -I {} echo -n "-V {} ") --merge-input-intervals true --reader-threads 30
$SOFTWAREDIR/gatk CreateSomaticPanelOfNormals -R $GENOMEIDX1 --germline-resource $GENOMEIDX4 -V gendb://pon_db -O PANEL_OF_NORMAL.vcf.gz

mkdir mutect-vcf
chmod 777 mutect-vcf
MUTECTVCFDIR="$BASEDIRDATA/mutect-vcf"
PONIDX="$BASEDIRDATA/PANEL_OF_NORMAL.vcf.gz"

for file1 in *T*_sm.bam; do
    file2="${file1/T/N}"
    # Extract sample name
    samplet=${file1}
    samplet=${samplet/_sm.bam}
    samplen=${file2}
    samplen=${samplen/_sm.bam}
    sample=${file1/T_*_sm.bam}
    echo "[$stime] $sample is being processed"
    $SOFTWAREDIR/gatk --java-options "-Xmx24G" Mutect2 --native-pair-hmm-threads 60 -R $GENOMEIDX1 -I $BASEDIRDATA/${file1} -I $BASEDIRDATA/${file2} -tumor TUMOR_${samplet} -normal NORMAL_${samplen} --germline-resource $GENOMEIDX4 -pon $PONIDX --genotype-pon-sites true --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O $MUTECTVCFDIR/${sample}_mutect.vcf
    echo "[$stime] $sample sample is processed; GOOD LUCK with annotation"
done
