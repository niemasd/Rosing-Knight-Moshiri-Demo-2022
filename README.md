These plugin descriptions are only for the demo. They may translate to long-term use as well, but they may need to be modified slightly (notes are left throughout regarding these potential improvements for long-term use). Much of this was adapted from [Tianqi's notes](https://git.ucsd.edu/tiz014/fpga-to-qiita-demo/snippets/33). Some notations I'll use in this document:

* `{...}` = file
* `[...]` = CLI option

# DeepVariant Pipeline
## Plugin 1: Map Reads to Human Genome
### CPU Commands
```bash
bwa index {REF_GENOME_FASTA} # should probably be done before the plugin is executed
bwa mem -t [THREADS] -k 21 {REF_GENOME_FASTA} {INPUT_READS} > {OUTPUT_SAM}

# everything below this is the same as FPGA
samtools view -@ [THREADS] -b -S {OUTPUT_SAM} > {OUTPUT_BAM}
samtools sort -@ [THREADS] {OUTPUT_BAM} -o {OUTPUT_SORTED_BAM}
samtools index -@ [THREADS] {OUTPUT_SORTED_BAM}
```

### FPGA Commands
```bash
fpga_alg_dv -tables {TABLES_FOLDER} -r {INPUT_READS} -o {OUTPUT_SAM}

# everything below this is the same as CPU
samtools view -@ [THREADS] -b -S {OUTPUT_SAM} > {OUTPUT_BAM}
samtools sort -@ [THREADS] {OUTPUT_BAM} -o {OUTPUT_SORTED_BAM}
samtools index -@ [THREADS] {OUTPUT_SORTED_BAM}
```

### Descriptions of Files and CLI Options
The only output files we care about are `{OUTPUT_SORTED_BAM}` and the index file produced by `samtools index`. All other intermediate files can be destroyed when the plugin finishes executing (and commands can be piped together if desired, but perhaps not for the demo to be able to get step-by-step runtimes if desired).

* `{REF_GENOME_FASTA}` = The reference genome FASTA (`ucsc_hg19.chr21.fa` for the demo)
* `[THREADS]` = number of threads (e.g. 16)
* `{INPUT_READS}` = Input FASTQ file (`BGISEQ_PE100_NA12878.sorted.chr21.fastq` for the demo)
  * For the demo, the FASTQ is assumed to already be trimmed
  * For long-term use, we would need to add a trimming step at the beginning of this plugin
* `{OUTPUT_SAM}` = The initial read mappings as a SAM file
* `{OUTPUT_BAM}` = The conversion of `{OUTPUT_SAM}` to the BAM format
* `{OUTPUT_SORTED_BAM}` = The result of sorting `{OUTPUT_BAM}`
