#Latest release via PyPI
#To install the latest release via PyPI using pip:

pip install sequenza-utils

#Development version
#To use the test suite for the package is necessary to install also bwa Using the latest development version directly from the git repository:

git clone https://bitbucket.org/sequenzatools/sequenza-utils
cd sequenza-utils
python setup.py test
python setup.py install
sequenza-utils -h
usage: sequenza-utils [-h] [-v]
                      {bam2seqz,gc_wiggle,pileup2acgt,seqz_binning,seqz_merge,snp2seqz}
                      ...

Sequenza Utils is a collection of tools primarily design to convert bam, pileup and vcf files to seqz files, the format used in the sequenza R package

positional arguments:
    bam2seqz        Process a paired set of BAM/pileup files (tumor and
                    matching normal), and GC-content genome-wide
                    information, to extract the common positions withA and
                    B alleles frequencies
    gc_wiggle       Given a fasta file and a window size it computes the GC
                    percentage across the sequences, and returns a file in
                    the UCSC wiggle format.
    pileup2acgt     Parse the format from the samtools mpileup command, and
                    report the occurrence of the 4 nucleotides in each
                    position.
    seqz_binning    Perform the binning of the seqz file to reduce file
                    size and memory requirement for the analysis.
    seqz_merge      Merge two seqz files
    snp2seqz        Parse VCFs and other variant and coverage formats to
                    produce seqz files

optional arguments:
  -h, --help        show this help message and exit
  -v, --verbose     Show all logging information

This is version 3.0.0 - Favero Francesco - 9 May 2019


sequenza-utils bam2seqz -h
error: the following arguments are required: -n/--normal, -t/--tumor, -gc
usage: sequenza-utils bam2seqz [-p] -n NORMAL -t TUMOR -gc GC [-F FASTA] [-o OUT] [-n2 NORMAL2] [-C CHR [CHR ...]] [--parallel NPROC] [-S SAMTOOLS] [-T TABIX] [-q QLIMIT] [-f QFORMAT] [-N N] [--hom HOM]
                               [--het HET] [--het_f HET_F]

Input/Output:
  Input and output files.

  -p, --pileup          Use pileups as input files instead of BAMs.
  -n NORMAL, --normal NORMAL
                        Name of the BAM/pileup file from the reference/normal sample
  -t TUMOR, --tumor TUMOR
                        Name of the BAM/pileup file from the tumor sample
  -gc GC                The GC-content wiggle file
  -F FASTA, --fasta FASTA
                        The reference FASTA file used to generate the intermediate pileup. Required when input are BAM
  -o OUT, --output OUT  Name of the output file. To use gzip compression name the file ending in .gz. Default STDOUT.
  -n2 NORMAL2, --normal2 NORMAL2
                        Optional BAM/pileup file used to compute the depth.normal and depth-ratio, instead of using the normal BAM.

Genotype:
  Options regarding the genotype filtering.

  --hom HOM             Threshold to select homozygous positions. Default 0.9.
  --het HET             Threshold to select heterozygous positions. Default 0.25.
  --het_f HET_F         Threshold of frequency in the forward strand to trust heterozygous calls. Default -0.2 (Disabled, effective with values >= 0).

Subset indexed files:
  Option regarding samtools and bam indexes.

  -C CHR [CHR ...], --chromosome CHR [CHR ...]
                        Argument to restrict the input/output to a chromosome or a chromosome region. Coordinate format is Name:pos.start-pos.end, eg: chr17:7565097-7590856, for a particular region; eg:
                        chr17, for the entire chromosome. Chromosome names can checked in the BAM/pileup files and are depending on the FASTA reference used for alignment. Default behavior is to not
                        selecting any chromosome.
  --parallel NPROC      Defines the number of chromosomes to run in parallel. The output will be divided in multiple files, one for each chromosome. The file name will be composed by the output argument
                        (used as prefix) and a chromosome name given by the chromosome argument list. This imply that both output and chromosome argument need to be set correctly.
  -S SAMTOOLS, --samtools SAMTOOLS
                        Path of samtools exec file to access the indexes and compute the pileups. Default "samtools"
  -T TABIX, --tabix TABIX
                        Path of the tabix binary. Default "tabix"

Quality and Format:
  Options that change the quality threshold and format.

  -q QLIMIT, --qlimit QLIMIT
                        Minimum nucleotide quality score for inclusion in the counts. Default 20.
  -f QFORMAT, --qformat QFORMAT
                        Quality format, options are "sanger" or "illumina". This will add an offset of 33 or 64 respectively to the qlimit value. Default "sanger".
  -N N                  Threshold to filter positions by the sum of read depth of the two samples. Default 20.




sequenza-utils gc_wiggle -h
error: the following arguments are required: -f/--fasta
usage: sequenza-utils gc_wiggle -f FASTA [-o OUT] [-w WINDOW]

optional arguments:
  -f FASTA, --fasta FASTA
                        the fasta file. It can be a file name or "-" to use STDIN
  -o OUT                Output file "-" for STDOUT
  -w WINDOW             The window size to calculate the GC-content percentage
  
  
  
  
sequenza-utils pileup2acgt -h
error: the following arguments are required: -p/--mpileup
usage: sequenza-utils pileup2acgt -p MPILEUP [-o OUTPUT] [-n N] [-q QLIMIT] [--no-end] [--no-start] [-f QFORMAT]

Output:
  Arguments that involve the output destination.

  -p MPILEUP, --mpileup MPILEUP
                        Name of the input mpileup (SAMtools) file. If the filename ends in .gz it will be opened in gzip mode. If the file name is - it will be read from STDIN.
  -o OUTPUT, --output OUTPUT
                        Name of the output file. To use gzip compression name the file ending in .gz. Default STDOUT.

Filtering and Format:
  Arguments that apply some filter to process the mpileup.

  -n N                  The minimum required read depth on a position to test for mutation.
  -q QLIMIT, --qlimit QLIMIT
                        Minimum nucleotide quality score filter.
  --no-end              Discard the base located at the end of the read
  --no-start            Discard the base located at the start of the read
  -f QFORMAT, --qformat QFORMAT
                        Quality format, options are "sanger" or "illumina". This will add an offset of 33 or 64 respectively to the qlimit value.

  
  
sequenza-utils seqz_binning -h
error: the following arguments are required: -s/--seqz
usage: sequenza-utils seqz_binning -s SEQZ [-w WINDOW] [-o OUT] [-T TABIX]

optional arguments:
  -s SEQZ, --seqz SEQZ  A seqz file.
  -w WINDOW, --window WINDOW
                        Window size used for binning the original seqz file. Default is 50.
  -o OUT                Output file "-" for STDOUT
  -T TABIX, --tabix TABIX
                        Path of the tabix binary. Default "tabix"






sequenza-utils seqz_merge -h
error: the following arguments are required: -1/--seqz1, -2/--seqz2
usage: sequenza-utils seqz_merge [-o OUTPUT] -1 S1 -2 S2 [-T TABIX]

Output:
  Output arguments

  -o OUTPUT, --output OUTPUT
                        Output file. For gzip compressed output name the file ending in .gz. Default STDOUT

Input:
  Input files

  -1 S1, --seqz1 S1     First input file. If both input files contain the same line, the information in the first file will be used
  -2 S2, --seqz2 S2     Second input file

Programs:
  Option to call and control externa programs

  -T TABIX, --tabix TABIX
                        Path of the tabix binary. Default "tabix"



-S, --samples-file

sequenza-utils snp2seqz -h
error: the following arguments are required: -v/--vcf, -gc
usage: sequenza-utils snp2seqz [-o OUTPUT] -v VCF -gc GC [--vcf-depth VCF_DEPTH_TAG] [--vcf-samples {n/t,t/n}] [--vcf-alleles VCF_ALLELES_TAG] [--preset {caveman,mutect,mpileup,strelka2_som}] [--hom HOM]
                               [--het HET] [--het_f HET_F] [-N N] [-T TABIX]

Output:
  Output arguments

  -o OUTPUT, --output OUTPUT
                        Output file. For gzip compressed output name the file ending in .gz. Default STDOUT

Input:
  Input files

  -v VCF, --vcf VCF     VCF input file
  -gc GC                The GC-content wiggle file

VCF:
  Parsing option for the VCF file

  --vcf-depth VCF_DEPTH_TAG
                        Column separated VCF tags in the format column for the read depth for the normal and for the tumor. Default "DP:DP"
  --vcf-samples {n/t,t/n}
                        Order of the normal and tumor sample in the VCF column, choices are "n/t" or "t/n". Default "n/t"
  --vcf-alleles VCF_ALLELES_TAG
                        Column separated VCF tags in the format column for the alleles depth for the normal and for the tumor. Default "AD:AD"
  --preset {caveman,mutect,mpileup,strelka2_som}
                        Preset set of options to parse VCF from popular variant callers

Genotype:
  Genotype filtering options

  --hom HOM             Threshold to select homozygous positions. Default 0.9
  --het HET             Threshold to select heterozygous positions. Default 0.25.
  --het_f HET_F         Threshold of frequency in the forward strand to trust heterozygous calls. Default -0.2 (Disabled, effective with values >= 0).

Programs:
  Option to call and control externa programs

  -T TABIX, --tabix TABIX
                        Path of the tabix binary. Default "tabix"

Filters:
  Filter output file by various parameters

  -N N                  Threshold to filter positions by the sum of read depth of the two samples. Default 20.

