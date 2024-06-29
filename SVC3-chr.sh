#Variant filter and annotation
SOFTWAREDIR="/home/cml/Downloads/gatk-4.2.6.1"
DATABASEDIR="/home/cml/Downloads/snpEff/data/GRCh38.p13"
BASEDIRINDEX="/home/cml/humandb/index-GRCh38.p13/HG38.p13-chr"
GENOMEIDX1="$BASEDIRINDEX/human_ref38_patch13.fna"
MUTECTVCFDIR="$PWD"
SNPEFFDIR="/home/cml/Downloads/snpEff"
ANNOVRDIR="/home/cml/Downloads/annovar"
ANNOVRDB="/home/cml/Downloads/annovar/hg38"
VCF2MAFDIR="/home/cml/Downloads/vcf2maf-1.6.21"
for file in *_mutect.vcf; do
sample1=${file}
echo $sample1
sample=${sample1/_mutect.vcf/}
echo $sample
mkdir ${sample}_CRC
chmod 777 ${sample}_CRC
mkdir ${sample}_CRC/ANNOVAR
chmod 777 ${sample}_CRC/ANNOVAR
mkdir MAF
chmod 777 MAF
$SOFTWAREDIR/gatk FilterMutectCalls -R $GENOMEIDX1 -V $MUTECTVCFDIR/${sample}_mutect.vcf -O $MUTECTVCFDIR/${sample}_CRC/${sample}_filtered.vcf
grep "#" $MUTECTVCFDIR/${sample}_CRC/${sample}_filtered.vcf >> $MUTECTVCFDIR/${sample}_CRC/${sample}_header.vcf
awk '$7=="PASS"' $MUTECTVCFDIR/${sample}_CRC/${sample}_filtered.vcf > $MUTECTVCFDIR/${sample}_CRC/${sample}_PASS_wh.vcf
cat $MUTECTVCFDIR/${sample}_CRC/${sample}_header.vcf $MUTECTVCFDIR/${sample}_CRC/${sample}_PASS_wh.vcf > $MUTECTVCFDIR/${sample}_CRC/${sample}_PASS.vcf
java -Xmx4G -jar $SNPEFFDIR/snpEff.jar -v GRCh38.p13 $MUTECTVCFDIR/${sample}_CRC/${sample}_PASS.vcf > $MUTECTVCFDIR/${sample}_CRC/${sample}_ann.vcf
java -jar $SNPEFFDIR/SnpSift.jar annotate $DATABASEDIR/hg38_1000GENOMES-phase_3.chr.vcf $MUTECTVCFDIR/${sample}_CRC/${sample}_ann.vcf > $MUTECTVCFDIR/${sample}_CRC/${sample}_1000G.vcf
cat $MUTECTVCFDIR/${sample}_CRC/${sample}_1000G.vcf | java -jar $SNPEFFDIR/SnpSift.jar filter "(GEN[*].AD[1] >= 10)" > $MUTECTVCFDIR/${sample}_CRC/${sample}_hf.vcf
cat $MUTECTVCFDIR/${sample}_CRC/${sample}_hf.vcf   |  java -Xmx64g -jar $SNPEFFDIR/SnpSift.jar filter "((! exists MAF) | (MAF <= 0.01))" > $MUTECTVCFDIR/${sample}_CRC/${sample}_MAF01.vcf
java -jar $SNPEFFDIR/SnpSift.jar annotate $DATABASEDIR/dbSNP155_GRCh38.chr.vcf $MUTECTVCFDIR/${sample}_CRC/${sample}_MAF01.vcf > $MUTECTVCFDIR/${sample}_CRC/${sample}_dbsnp.vcf
java -jar $SNPEFFDIR/SnpSift.jar annotate $DATABASEDIR/Cosmic95.chr.vcf $MUTECTVCFDIR/${sample}_CRC/${sample}_dbsnp.vcf > $MUTECTVCFDIR/${sample}_CRC/${sample}_c95.vcf
sed '/#/d' $MUTECTVCFDIR/${sample}_CRC/${sample}_c95.vcf > $MUTECTVCFDIR/${sample}_CRC/${sample}_cosmic_nh.vcf
cat $MUTECTVCFDIR/${sample}_CRC/${sample}_header.vcf $MUTECTVCFDIR/${sample}_CRC/${sample}_cosmic_nh.vcf > $MUTECTVCFDIR/${sample}_CRC/${sample}_cosmic.vcf

perl $ANNOVRDIR/table_annovar.pl $MUTECTVCFDIR/${sample}_CRC/${sample}_cosmic.vcf $ANNOVRDB -buildver hg38 -out $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno -remove -protocol ensGene,cytoBand,exac03,esp6500siv2_all,gnomad211_exome,icgc28,intervar_20180118,nci60,regsnpintron,revel,tmcsnpdb,clinvar_20221231,dbnsfp42c,dbscsnv11,avsnp150 -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish 
echo " Annotation using ANNOVAR is complete for ${sample}"
echo " Let's convert the annovar chr id to vep numeric id"
sed -i 's/chrMT/MT/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chrY/Y/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chrX/X/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr22/22/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr21/21/g'  $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr20/20/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr19/19/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr18/18/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr17/17/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr16/16/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr15/15/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr14/14/g'  $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr13/13/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr12/12/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr11/11/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr10/10/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr9/9/g'  $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr8/8/g'  $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr7/7/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr6/6/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr5/5/g' $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr4/4/g'  $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr3/3/g'  $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr2/2/g'  $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
sed -i 's/chr1/1/g'  $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf
echo " annovar chr id to vep numeric id conversion complete for ${sample}"

echo " Let's convert the annotated vcf format of ${sample} into maf format"
perl $VCF2MAFDIR/vcf2maf.pl --input-vcf $MUTECTVCFDIR/${sample}_CRC/ANNOVAR/${sample}_myanno.hg38_multianno.vcf --output-maf $MUTECTVCFDIR/MAF/${sample}_vep.maf
echo " VCF to MAF format conversion complete for ${sample}"
done
