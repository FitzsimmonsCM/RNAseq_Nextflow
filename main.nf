// Author: @FitzsimmonsCM
// This nextflow (nf) workflow takes RNAseq data from the NCI CCR genomics core and performs adaptor trimming, mapping, duplicate removal, and counting
// This workflow is designed to be run with Illumina TruSeq RNAseq data that has already been base-called and demultiplexed
// Additionally, this workflow is designed to work in a cluster agnostic manner by using docker / singularity containers.

// params
params.fastq_input = "$projectDir/demo_data/"

// channels
raw_fastq_pair_dir = Channel.fromPath(params.fastq_input)

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

// map to the human genome with STAR
// 1 INPUT (tuple with paired-end read data); 1 OUTPUT (1 mapped BAM file)
process STAR_map {

container: 'quay.io/biocontainers/star:2.7.6a--0'

input:
  tuple file('R1.trimm.fastq.gz'), file('R2.trimm.fastq.gz') from trim_adaptors

output:
  file 'mapped.sorted.bam' into genome_mapped

"""

GENOME=/fdb/STAR_indices/2.7.6a/GENCODE/Gencode_human/release_35/genes-100
STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \       # non-default
    --genomeDir $GENOME \                     # required
    --readFilesCommand zcat \                 # non-default
    --sjdbOverhang 100 \                      # default
    --readFilesIn $INPUT \                    # required
    --outSAMtype BAM SortedByCoordinate \     # non-default
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \ # non-default
    --outReadsUnmapped Fastx \                # default
    --outFileNamePrefix map_gencodeV35/$OUTNAME \ # required
    –-outFilterType BySJout \                 # non-default
    --outFilterMultimapNmax 10 \              # default
    --outFilterMismatchNmax 999 \             # non-default
    --outFilterMismatchNoverLmax 0.04 \       # non-default
    --alignIntronMin 20 \                     # non-default
    --alignIntronMax 1000000 \                # non-default
    --alignMatesGapMax 1000000 \              # non-default
    --alignSJoverhangMin 8 \                  # non-default
    --alignSJDBoverhangMin 1 \                # non-default
    --sjdbScore 1 \                           # non-default
    --outFilterMatchNminOverLread 0.66 \      # default
    --peOverlapNbasesMin 10 \                 # non-default
    --alignEndsProtrude 10 ConcordantPair     # non-default

"""


}


// Return a BAM file sorted by coordinates
// 1 INPUT (unsorted BAM file), 1 OUTPUT (1 BAM file Sorted By Coordinates)
// PicardTools requires sorted files or it will raise an error
process samtools_sort {

  container: 'quay.io/repository/biocontainers/samtools'

  input:
    file 'mapped.sorted.bam' from genome_mapped

  output:
    file 'mapped.sorted.bam' from samtools_sort

  """
// put something here

  """
}

  // Return a sorted BAM file with duplicates removed
  // 1 INPUT (sorted BAM file), 1 OUTPUT (1 sorted BAM file with duplicates removed)

  process picard_dedup {

    container: 'broadinstitute/picard:2.25.0'

    input:
      file 'mapped.sorted.bam' from samtools_sort

    output:
      file 'mapped.sorted.dedup.bam' from samtools_sort

"""
      java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
      I=$INPUT \
      O=$OUTPUT\
      ASSUME_SORTED=true \
      METRICS_FILE=$OUTNAME.mapped.metric.csv \
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
      VALIDATION_STRINGENCY=LENIENT \
      REMOVE_DUPLICATES=TRUE

"""
}

 // Return a count file to be input into DEseq2 for differential expression analysis
 // 1 INPUT (1 deduplicated BAM file), 1 OUTPUT (1 count file)
process htseq_count {

container:

input:
  file 'mapped.sorted.dedup.bam' from samtools_sort

  output:
    file 'BC.count' from htseq_count
}

"""
module load htseq
htseq-count -m intersection-nonempty -s reverse -f bam -r pos -t exon -i gene_id $INPUT > $OUTPUT

"""
