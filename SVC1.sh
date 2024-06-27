echo "Welcome to Exome structural variant calling Pipeline: "Stage 01""
NUMCPUS=60
SOFTWAREDIR="/pathto/gatk"
BASEDIRINDEX="/pathto/HG38-NC"
GENOMEDICT="$BASEDIRINDEX/human_ref38_patch13.dict"
GENOMEIDX1="$BASEDIRINDEX/human_ref38_patch13.fna"
GENOMEIDX2="$BASEDIRINDEX/Twist_ComprehensiveExome_targets_hg38.bed"
GENOMEIDX3="$BASEDIRINDEX/Twist_ComprehensiveExome_targets_hg38.interval_list"
GENOMEIDX4="$BASEDIRINDEX/Twist_ComprehensiveExome_targets_hg38.preprocessed.interval_list"
GENOMEIDX5="$BASEDIRINDEX/Twist_ComprehensiveExome_targets_hg38.annotated.tsv"
GENOMEIDX6="$BASEDIRINDEX/contig_ploidy_priors.tsv"
BAMDIR="/pathto/bam-dir"
READCOUNTDIR="/pathto/readcounts-dir"
mkdir cnv
chmod 777 cnv
OUTPUTDIR="/pathto/cnv"
BASEDIRDATA="$PWD"
$SOFTWAREDIR/gatk PreprocessIntervals -L $GENOMEIDX2 R $GENOMEIDX1 -imr OVERLAPPING_ONLY --bin-length 0 --padding 250 -O Twist_ComprehensiveExome_targets_hg38.preprocessed.interval_list
$SOFTWAREDIR/gatk AnnotateIntervals -L $ENOMEIDX4 -R $GENOMEIDX1 -imr OVERLAPPING_ONLY -O Twist_ComprehensiveExome_targets_hg38.annotated.tsv

for file in *_sm.bam; do
    sample1=${file}
    sample=${sample1/_sm.bam/}
    echo "[$stime] $sample is being processed"
    $SOFTWAREDIR/gatk CollectReadCounts -L $GENOMEIDX4 -R $GENOMEIDX1 -imr OVERLAPPING_ONLY -I $BASEDIRDATA/${sample}_sm.bam --format TSV -O $READCOUNTDIR/${sample}.tsv
    done

$SOFTWAREDIR/gatk FilterIntervals -L $GENOMEIDX4 --annotated-intervals $GENOMEIDX5 $(ls $READCOUNTDIR/*.tsv | xargs -I {} echo -n "-I {} ") -imr OVERLAPPING_ONLY -O $OUTPUTDIR/cohort.gc.filtered.interval_list

GENOMEIDX7="$OUTPUTDIR/cohort.gc.filtered.interval_list"

$SOFTWAREDIR/gatk DetermineGermlineContigPloidy -L $GENOMEIDX7 --interval-merging-rule OVERLAPPING_ONLY $(ls $READCOUNTDIR/*.tsv | xargs -I {} echo -n "-I {} ") --contig-ploidy-priors $GENOMEIDX5 --output . --output-prefix ploidy --verbosity DEBUG

$SOFTWAREDIR/gatk DetermineGermlineContigPloidy -L $GENOMEIDX7 --interval-merging-rule OVERLAPPING_ONLY $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --contig-ploidy-priors $GENOMEIDX7 --output . --output-prefix ploidy --verbosity DEBUG


$SOFTWAREDIR/gatk IntervalListTools --INPUT cohort.gc.filtered.interval_list --SUBDIVISION_MODE INTERVAL_COUNT --SCATTER_CONTENT 5000 --OUTPUT scatter

$SOFTWAREDIR/gatk GermlineCNVCaller --run-mode COHORT -L scatter/01_scattered.interval_list --annotated-intervals /pathto/Twist_ComprehensiveExome_targets_hg38.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ./ploidy-calls $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --output 01_sc --output-prefix 01_sc_cohort --verbosity DEBUG

$SOFTWAREDIR/gatk GermlineCNVCaller --run-mode COHORT -L scatter/02_scattered.interval_list --annotated-intervals /pathto/Twist_ComprehensiveExome_targets_hg38.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ./ploidy-calls $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --output 02_sc --output-prefix 02_sc_cohort --verbosity DEBUG

$SOFTWAREDIR/gatk GermlineCNVCaller --run-mode COHORT -L scatter/03_scattered.interval_list --annotated-intervals /pathto/Twist_ComprehensiveExome_targets_hg38.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ./ploidy-calls $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --output 03_sc --output-prefix 03_sc_cohort --verbosity DEBUG

$SOFTWAREDIR/gatk GermlineCNVCaller --run-mode COHORT -L scatter/04_scattered.interval_list --annotated-intervals /pathto/Twist_ComprehensiveExome_targets_hg38.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ./ploidy-calls $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --output 04_sc --output-prefix 04_sc_cohort --verbosity DEBUG

$SOFTWAREDIR/gatk GermlineCNVCaller --run-mode COHORT -L scatter/05_scattered.interval_list --annotated-intervals /pathto/Twist_ComprehensiveExome_targets_hg38.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ./ploidy-calls $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --output 05_sc --output-prefix 05_sc_cohort --verbosity DEBUG

$SOFTWAREDIR/gatk GermlineCNVCaller --run-mode COHORT -L scatter/06_scattered.interval_list --annotated-intervals /pathto/Twist_ComprehensiveExome_targets_hg38.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ./ploidy-calls $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --output 06_sc --output-prefix 06_sc_cohort --verbosity DEBUG

$SOFTWAREDIR/gatk GermlineCNVCaller --run-mode COHORT -L scatter/07_scattered.interval_list --annotated-intervals /pathto/Twist_ComprehensiveExome_targets_hg38.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ./ploidy-calls $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --output 07_sc --output-prefix 07_sc_cohort --verbosity DEBUG

$SOFTWAREDIR/gatk GermlineCNVCaller --run-mode COHORT -L scatter/08_scattered.interval_list --annotated-intervals /pathto/Twist_ComprehensiveExome_targets_hg38.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ./ploidy-calls $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --output 08_sc --output-prefix 08_sc_cohort --verbosity DEBUG

$SOFTWAREDIR/gatk GermlineCNVCaller --run-mode COHORT -L scatter/09_scattered.interval_list --annotated-intervals /pathto/Twist_ComprehensiveExome_targets_hg38.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ./ploidy-calls $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --output 09_sc --output-prefix 09_sc_cohort --verbosity DEBUG

$SOFTWAREDIR/gatk GermlineCNVCaller --run-mode COHORT -L scatter/10_scattered.interval_list --annotated-intervals /pathto/Twist_ComprehensiveExome_targets_hg38.annotated.tsv --interval-merging-rule OVERLAPPING_ONLY --contig-ploidy-calls ./ploidy-calls $(ls ./*.tsv | xargs -I {} echo -n "-I {} ") --output 10_sc --output-prefix 10_sc_cohort --verbosity DEBUG


$SOFTWAREDIR/gatk PostprocessGermlineCNVCalls --model-shard-path 01_sc/01_sc_cohort-model --model-shard-path 02_sc/02_sc_cohort-model --model-shard-path 03_sc/03_sc_cohort-model --model-shard-path 04_sc/04_sc_cohort-model --model-shard-path 05_sc/05_sc_cohort-model --model-shard-path 06_sc/06_sc_cohort-model --model-shard-path 07_sc/07_sc_cohort-model --model-shard-path 08_sc/08_sc_cohort-model --model-shard-path 09_sc/09_sc_cohort-model --model-shard-path 10_sc/10_sc_cohort-model --calls-shard-path 01_sc/01_sc_cohort-calls --calls-shard-path 02_sc/02_sc_cohort-calls --calls-shard-path 03_sc/03_sc_cohort-calls --calls-shard-path 04_sc/04_sc_cohort-calls --calls-shard-path 05_sc/05_sc_cohort-calls --calls-shard-path 06_sc/06_sc_cohort-calls --calls-shard-path 07_sc/07_sc_cohort-calls --calls-shard-path 08_sc/08_sc_cohort-calls --calls-shard-path 09_sc/09_sc_cohort-calls --calls-shard-path 10_sc/10_sc_cohort-calls --contig-ploidy-calls ploidy-calls --sample-index 0 --output-genotyped-intervals gcnv/sample_genotyped-intervals.vcf.gz --output-genotyped-segments gcnv/sample_segments-intervals.vcf.gz --output-denoised-copy-ratios gcnv/sample_denoised-copy-ratios.csv --sequence-dictionary /pathto/human_ref38_patch13.dict



# Get the list of directories
model_dirs=$(ls -d *_sc/*_sc_cohort-model)
calls_dirs=$(ls -d *_sc/*_sc_cohort-calls)
# Create the model paths string
model_paths=$(echo $model_dirs | xargs -n1 echo "--model-shard-path ")
calls_paths=$(echo $calls_dirs | xargs -n1 echo "--calls-shard-path ")

# Read file pairs from a text file
file_pairs=$(cat file_pairs.txt)

# Iterate over each line in the file pairs
for files in $file_pairs; do
        # Split the pair into two files
        IFS=',' read -r file1 file2 <<< "$files"

        # Extract sample names
        sampleid=${file1}
        echo $sampleid
        sampleindex=${file2}
        echo $sampleindex
	$SOFTWAREDIR/gatk PostprocessGermlineCNVCalls $model_paths $calls_paths --contig-ploidy-calls ploidy-calls --sample-index ${sampleindex} --output-genotyped-intervals gcnv/${sampleid}_genotyped-intervals.vcf.gz --output-genotyped-segments gcnv/${sampleid}_segments-intervals.vcf.gz --output-denoised-copy-ratios gcnv/${sampleid}_denoised-copy-ratios.csv --sequence-dictionary $DICTIONARY
done
