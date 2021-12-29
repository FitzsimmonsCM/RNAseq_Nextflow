// Author: @FitzsimmonsCM
// This nextflow (nf) workflow takes RNAseq data from the NCI CCR genomics core and performs adaptor trimming, mapping, duplicate removal, and counting
// This workflow is designed to be run with Illumina TruSeq RNAseq data that has already been base-called and demultiplexed
// Additionally, this workflow is designed to work in a cluster agnostic manner by using docker / singularity containers.

// params
params.fastq_input = "$projectDir/demo_data/"

// channels
raw_fastq_pair_dir = Channel.fromPath(params.fastq_input)

//singularity
//singularity.enabled = true

// remove illumina adaptors
// 2 INPUTS (PAIRED-END sequencing reads); 1 OUTPUT (paired, trimmed reads--returned as 1 tuple)
process cut_adapters {

    container 'dceoy/cutadapt:latest'

    input:
      file fastq_dir from raw_fastq_pair_dir

    output:
      tuple file('R1.trimm.fastq.gz'), file('R2.trimm.fastq.gz') into trim_adaptors

    """
    cutadapt --nextseq-trim=20 --trim-n -m 20 \
      --cores=$task.cpus \
      -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      -o R1.trimm.fastq.gz -p R2.trimm.fastq.gz \
      $fastq_dir/*R1*.fastq.gz \
      $fastq_dir/*R2*.fastq.gz
    """
}

// map to the genome with STAR
// 1 INPUT (tuple with paired-end read data); 1 OUTPUT (1 mapped BAM file)
/*process STAR_map {

container: STAR??

input:
  tuple file('R1.trimm.fastq.gz'), file('R2.trimm.fastq.gz') from trim_adaptors

output:
  file 'mapped.sorted.bam' into genome_mapped

"""
module load samtools/1.11         || fail "could not load samtools module"
module load STAR/2.7.6a         || fail "could not load STAR module"

cd /data/BatistaLab_NGS/UOK_manuscript/RNAseq/HFM7WDRXY_2021.04_repeatRNAseq_data/mapped || fail "no such directory"
mkdir -p map_gencodeV35
GENOME=/fdb/STAR_indices/2.7.6a/GENCODE/Gencode_human/release_35/genes-100


STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $GENOME \
    --readFilesCommand zcat \
    --sjdbOverhang 100 \
    --readFilesIn $INPUT \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix map_gencodeV35/$OUTNAME \
    –-outFilterType BySJout \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --outFilterMatchNminOverLread 0.66 \
    --quantMode TranscriptomeSAM \
    --peOverlapNbasesMin 10 \
    --alignEndsProtrude 10 ConcordantPair

"""


}


// Return a BAM file sorted by coordinates
// 1 INPUT (unsorted BAM file), 1 OUTPUT (1 BAM file Sorted By Coordinates)
// PicardTools requires sorted files or it will raise an error
*/
