#Latest release via PyPI
#To install the latest release via PyPI using pip:

pip install sequenza-utils

#Development version
#To use the test suite for the package is necessary to install also bwa Using the latest development version directly from the git repository:

git clone https://bitbucket.org/sequenzatools/sequenza-utils
cd sequenza-utils
python setup.py test
python setup.py install

Generate GC reference file
The GC content source required to generate seqz files must to be in the wiggle track format (WIG). In order to generate the wig file from any fasta file use the gc_wiggle program.

sequenza-utils gc_wiggle --fasta human_ref38_patch13.fna -w 50 -o human_ref38_patch13.gc50base.wig.gz



run_samtools() { sample="$1"; samtools mpileup -f /pathto/human_ref38_patch13.fna "${sample}_bwa_sort.bam" > "${sample}".pileup; }; export -f run_samtools; ls *_bwa_sort.bam | sed 's/_bwa_sort\.bam$//' | xargs -n 1 -P 40 -I {} bash -c 'run_samtools "$@"' _ {} 2>&1 | tee samtools.txt

echo "Welcome to CNV analysis Pipeline: "Stage 02""
BASEDIRINDEX="/pathto/HG38dir"
GENOMEIDX1="$BASEDIRINDEX/human_ref38_patch13.fna"
GENOMEIDX2="$BASEDIRINDEX/human_ref38_patch13.gc50base.wig.gz"
BASEDIRDATA="/pathtobam"
SEQUENZADIR="$PWD"

# Read file pairs from a text file
file_pairs=$(cat file_pairs.txt)

# Iterate over each line in the file pairs
for files in $file_pairs; do
	(
	# Split the pair into two files
	IFS=',' read -r file1 file2 <<< "$files"

	# Extract sample names
	samplen=${file1}
	samplen=${samplen/_sm.bam}
	echo $samplen
	samplet=${file2}
	samplet=${samplet/_sm.bam}
	echo $samplet
	
	#make directory for samples
	mkdir $samplet
	chmod 777 $samplet
	sequenza-utils bam2seqz -p -n ${BASEDIRDATA}/${file1} -t ${BASEDIRDATA}/${file2} --fasta ${GENOMEIDX1} -gc ${GENOMEIDX2} -C NC_000001.11 NC_000002.12 NC_000003.12 NC_000004.12 NC_000005.10 NC_000006.12 NC_000007.14 NC_000008.11 NC_000009.12 NC_000010.11 NC_000011.10 NC_000012.12 NC_000013.11 NC_000014.9 NC_000015.10 NC_000016.10 NC_000017.11 NC_000018.10 NC_000019.10 NC_000020.11 NC_000021.9 NC_000022.11 NC_000023.11 NC_000024.10 --parallel 24 -o ${BASEDIRDATA}/${samplet}/${samplet}.seqz.gz
done

#get result folder-wise
find . -type d | parallel -j 60 'cd {} && for file in *.seqz.gz; do sample=$(basename "$file" .seqz.gz); sequenza-utils seqz_binning -w 50 --seqz "$file" -o "${sample}.small.seqz.gz"; done'

#for all file present in one folder
find . -type f -name "*.seqz.gz" | parallel -j 24 'sample=$(basename {} .seqz.gz); sequenza-utils seqz_binning -w 50 --seqz {} -o ${sample}.small.seqz.gz'

seq_bin() { sample="$1"; sequenza-utils seqz_binning -w 50 --seqz "${sample}.seqz.gz" -o "${sample}.small.seqz.gz"; }; export -f seq_bin; ls ./*.seqz.gz | sed 's/.seqz\.gz$//' | xargs -n 1 -P 18 -I {} bash -c 'seq_bin "$@"' _ {} 2>&1 | tee seq_bin.txt

#Unzip ${sample}.small.seqz.gz
sudo parallel -j60 --eta --bar --progress gunzip -9 ::: */*

#Convert NC-id to chr-id
run_sed() { sample="$1"; echo " Let's convert the NC id to chr id"; sed -i -e 's/NC_000024.10/chrY/g' -e 's/NC_000023.11/chrX/g' -e 's/NC_000022.11/chr22/g' -e 's/NC_000021.9/chr21/g' -e 's/NC_000020.11/chr20/g' -e 's/NC_000019.10/chr19/g' -e 's/NC_000018.10/chr18/g' -e 's/NC_000017.11/chr17/g' -e 's/NC_000016.10/chr16/g' -e 's/NC_000015.10/chr15/g' -e 's/NC_000014.9/chr14/g' -e 's/NC_000013.11/chr13/g' -e 's/NC_000012.12/chr12/g' -e 's/NC_000011.10/chr11/g' -e 's/NC_000010.11/chr10/g' -e 's/NC_000009.12/chr9/g' -e 's/NC_000008.11/chr8/g' -e 's/NC_000007.14/chr7/g' -e 's/NC_000006.12/chr6/g' -e 's/NC_000005.10/chr5/g' -e 's/NC_000004.12/chr4/g' -e 's/NC_000003.12/chr3/g' -e 's/NC_000002.12/chr2/g' -e 's/NC_000001.11/chr1/g' "${sample}.small.seqz"; }; export -f run_sed; ls ./*.small.seqz | sed 's/.small\.seqz$//' | xargs -n 1 -P 30 -I {} bash -c 'run_sed "$@"' _ {} 2>&1 | tee nc2chr.txt

#bgzip seqz files
sudo parallel -j30 --eta --bar --progress bgzip ::: *

#Run tabix
run_tabix() { sample="$1"; tabix -f -s 1 -b 2 -e 2 -S 1 "${sample}_CRC.seqz.gz"; }; export -f run_tabix; ls ./*_CRC.seqz.gz | sed 's/_CRC.seqz\.gz$//' | xargs -n 1 -P 15 -I {} bash -c 'run_tabix "$@"' _ {} 2>&1 | tee tabix.txt

#Merge Chr2-24
for file in *CRC/*.seqz; do 
    sample=$(basename "$file" | sed "s/_NC.*//"); 
    cat ${sample}/${sample}_NC_000002.12.small.seqz \
        ${sample}/${sample}_NC_000003.12.small.seqz \
        ${sample}/${sample}_NC_000004.12.small.seqz \
        ${sample}/${sample}_NC_000005.10.small.seqz \
        ${sample}/${sample}_NC_000006.12.small.seqz \
        ${sample}/${sample}_NC_000007.14.small.seqz \
        ${sample}/${sample}_NC_000008.11.small.seqz \
        ${sample}/${sample}_NC_000009.12.small.seqz \
        ${sample}/${sample}_NC_000010.11.small.seqz \
        ${sample}/${sample}_NC_000011.10.small.seqz \
        ${sample}/${sample}_NC_000012.12.small.seqz \
        ${sample}/${sample}_NC_000013.11.small.seqz \
        ${sample}/${sample}_NC_000014.9.small.seqz \
        ${sample}/${sample}_NC_000015.10.small.seqz \
        ${sample}/${sample}_NC_000016.10.small.seqz \
        ${sample}/${sample}_NC_000017.11.small.seqz \
        ${sample}/${sample}_NC_000018.10.small.seqz \
        ${sample}/${sample}_NC_000019.10.small.seqz \
        ${sample}/${sample}_NC_000020.11.small.seqz \
        ${sample}/${sample}_NC_000021.9.small.seqz \
        ${sample}/${sample}_NC_000022.11.small.seqz \
        ${sample}/${sample}_NC_000023.11.small.seqz \
        ${sample}/${sample}_NC_000024.10.small.seqz | \
    gawk '{if (NR!=1 && $1 != "chromosome") {print $0}}' > ${sample}/${sample}_merged.small.seqz; 
done

#Merge Chr1 with Chr2-24
for file in *_merged.small.seqz; do 
    sample=${$file/_merged.small.seqz}; 
    cat ${sample}_NC_000001.11.small.seqz ${sample}_merged.small.seqz > ${sample}.small.seqz
done

for file in *seqz.gz; do sample=$(basename "$file" | sed "s/_NC.*//"); mkdir $sample/; done
for file in *seqz.gz; do sample=$(basename "$file" | sed "s/_NC.*//"); mv ${file} $sample/; done
for subdir in */; do subdir_name=$(basename "$subdir"); output_file="${subdir_name}_merged.small.seqz.gz"; (cd "$subdir" && zcat $(ls -v *.small.seqz.gz) | gawk '{if (NR!=1 && $1 != "chromosome") {print $0}}' | bgzip > "$output_file"); echo "Merged seqz files in $subdir into $output_file"; done
for subdir in */; do subdir_name=$(basename "$subdir"); output_file="${subdir_name}_merged.small.seqz"; (cd "$subdir" && cat $(ls -v *.small.seqz) | gawk '{if (NR!=1 && $1 != "chromosome") {print $0}}' > "$output_file"); echo "Merged seqz files in $subdir into $output_file"; done

cnvkit.py antitarget Twist_ComprehensiveExome_targets_hg38.bed -g access-10kb.hg38.bed -o antitarget-twist-hg38.bed



run_tabix() { sample="$1"; tabix -f -s 1 -b 2 -e 2 -S 1 "${sample}_CRC.small.seqz.gz"; }; export -f run_tabix; ls ./*_CRC.small.seqz.gz | sed 's/_CRC.small.seqz\.gz$//' | xargs -n 1 -P 15 -I {} bash -c 'run_tabix "$@"' _ {} 2>&1 | tee tabix.txt



