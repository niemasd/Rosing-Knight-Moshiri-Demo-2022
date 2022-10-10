These plugin descriptions are only for the demo. They may translate to long-term use as well, but they may need to be modified slightly (notes are left throughout regarding these potential improvements for long-term use). Much of this was adapted from [Tianqi's notes](https://git.ucsd.edu/tiz014/fpga-to-qiita-demo/snippets/33). Some notations I'll use in this document:

* `{...}` = file
* `[...]` = CLI option

# DeepVariant Pipeline
In this pipeline, we start with human reads (assumed to already be trimmed in the demo dataset), and we (1) map to the human reference genome, and then (2) call variants using DeepVariant.

## Plugin 1: Map Reads to Human Reference Genome
### CPU Commands
```bash
bwa index {REF_GENOME_FASTA} # should probably be done before the plugin is executed
bwa mem -t [THREADS] -k 21 {REF_GENOME_FASTA} {READS} > output.sam

# everything below this is the same as FPGA
samtools view -@ [THREADS] -b -S output.sam > output.bam
samtools sort -@ [THREADS] output.bam -o {SORTED_BAM}
samtools index -@ [THREADS] {SORTED_BAM}
```

### FPGA Commands
```bash
fpga_alg_dv -tables {TABLES_FOLDER} -r {READS} -o output.sam

# everything below this is the same as CPU
samtools view -@ [THREADS] -b -S output.sam > output.bam
samtools sort -@ [THREADS] output.bam -o {SORTED_BAM}
samtools index -@ [THREADS] {SORTED_BAM}
```

### Descriptions of Files and CLI Options
The only output files we care about are `{SORTED_BAM}` and the index file produced by `samtools index`. All other intermediate files can be destroyed when the plugin finishes executing (and commands can be piped together if desired, but perhaps not for the demo to be able to get step-by-step runtimes if desired).

* `{REF_GENOME_FASTA}` = The reference genome FASTA (`ucsc_hg19.chr21.fa` for the demo)
* `{TABLES_FOLDER}` = I'm not too sure what this is; Tianqi should clarify
  * My best guess is that this is some preprocessed "index" of the reference genome used by the FPGA mapper?
* `[THREADS]` = Number of threads (e.g. 16)
* `{READS}` = Input FASTQ file (`BGISEQ_PE100_NA12878.sorted.chr21.fastq` for the demo)
  * For the demo, the FASTQ is assumed to already be trimmed
  * For long-term use, we would need to add a trimming step at the beginning of this plugin
* `{SORTED_BAM}` = The sorted output BAM (which has been indexed)

## Plugin 2: Call Variants
### CPU Commands
```bash
sudo docker pull google/deepvariant:[DV_VERSION] # should probably be done before the plugin is executed
sudo docker run --name deepvariant0 -v {INPUT_DIR}:/input -v {OUTPUT_DIR}:/output -it google/deepvariant:[DV_VERSION] /bin/bash
sudo docker exec deepvariant0 /opt/deepvariant/bin/run_deepvariant \
  --model_type="WGS" --ref=/input/{REF_GENOME_FASTA} \
  --reads=/input/{SORTED_BAM} \
  --output_vcf=/output/{VCF} \
  --output_gvcf=/output/{GVCF} \
  --num_shards=[SHARDS] \
  --logging_dir=/output/{DV_LOG} \
  --dry_run=false
```

### FPGA Commands
None: we have not accelerated variant calling yet.

### Descriptions of Files and CLI Options
The only output files we care about are `{SORTED_BAM}` and the index file produced by `samtools index`. All other intermediate files can be destroyed when the plugin finishes executing (and commands can be piped together if desired, but perhaps not for the demo to be able to get step-by-step runtimes if desired).

* `{REF_GENOME_FASTA}` = The reference genome FASTA (`ucsc_hg19.chr21.fa` for the demo)
* `[DV_VERSION]` = DeepVariant version (`1.3.0` for the demo)
* `{INPUT_DIR}` = Directory containing the input files (`{REF_GENOME_FASTA}` and `{SORTED_BAM}`)
* `{OUTPUT_DIR}` = Directory to write the output files
* `{SORTED_BAM}` = The output of Plugin 1 (just the filename, not the full path)
* `{VCF}` = Output VCF file (just the filename, not the full path)
* `{GVCF}` = Output GVCF file (just the filename, not the full path)
* `[SHARDS]` = Number of shards to run in parallel (seemingly number of threads minus 1? e.g. 15)
* `{DV_LOG}` = Output DeepVariant log file (just the filename, not the full path)

# UniFrac Pipeline
In this pipeline, we accelerate host-filtering as well as the to-database read mapping prior to calling Woltka. According to the [Qiita Recommended Default Workflows](https://qiita.ucsd.edu/workflows/), we would just feed our data directly into the Qiita Woltka plugin, and the FPGA read mapper would swap in for whatever read mapper is already being used in the Qiita Woltka plugin (with all else remaining the same).

## Plugin 1: Host Filtering
I'm not sure if our input reads are already trimmed; would be good to check with Tianqi. We've accelerated quality trimming as well, so it would make sense to have an FPGA replacement for the `fastp-minimap2` plugin that uses the FPGA trimmer fed into the FPGA mapper. Tianqi and Tajana would know. The only information I was given was about the actual host filtering (i.e., mapping to human genome), so that's the only step I'll provide commands for here.

### CPU Commands
None: use the existing `fastp-minimap2` host-filtering plugin on Qiita.

### FPGA Commands
```bash
fpga_holst_filter -tables {TABLES_FOLDER} -r {READS} -o {FILTERED_READS}
```
* The spelling `fpga_holst_filter` ("ho**l**st" rather than "host") was how it was given to me by Tianqi; would be good to double check the spelling

### Descriptions of Files and CLI Options
* `{READS}` = Input FASTQ file
  * Double check with Tianqi that the demo reads are already trimmed
* `{TABLES_FOLDER}` = I'm not too sure what this is; Tianqi should clarify
  * My best guess is that this is some preprocessed "index" of the reference genome used by the FPGA mapper?
* `{FILTERED_READS}` = Host-filtered reads

## Plugin 2: Woltka
As mentioned, it doesn't seem like any new plugins need to be written here: it seems like the only change that needs to be made for the demo is that the existing Qiita Woltka plugin (which seems to map non-ribosomal reads to the database before calling Woltka) would need to be modified to add the FPGA read mapper as an option. Thus, in my CPU Commands and FPGA Commands sections below, I'll only mention the read-mapping step, as all other commands should remain exactly as-is within the Qiita Woltka plugin.

### CPU Commands
```bash
bowtie2 -p [THREADS]  --very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" -x db -U {READS} | pigz -p [THREADS] -c > {SAM_GZ}

# REST OF QIITA WOLTKA PLUGIN
```

### FPGA Commands
```bash
fpga_alg_mb -k 16 -tables {TABLES_FOLDER} -r {READS} -o {SAM}
pigz -p [THREADS] {SAM}

# REST OF QIITA WOLTKA PLUGIN
```

### Descriptions of Files and CLI Options
* `[THREADS]` = Number of threads (e.g. 16)
* `{READS}` = Input FASTQ file
* `{TABLES_FOLDER}` = I'm not too sure what this is; Tianqi should clarify
  * My best guess is that this is some preprocessed "index" of the reference genome used by the FPGA mapper?
* `{SAM}` = The output of the FPGA read mapper
* `{SAM_GZ}` = The compressed read mapping output, which will be fed into Woltka

# Per-Sample COVID-19 Pipeline
## Plugin 1: Map and Trim Reads
Map reads (amplicon sequencing) to SARS-CoV-2 reference genome and primer+quality trim.

### CPU Commands
```bash
minimap2 -t [THREADS] -a -x sr {REF_GENOME_FASTA} {READS_R1} {READS_R2} | samtools sort -@ [THREADS] -o untrimmed.bam
ivar trim -x 5 -e -i untrimmed.bam -b {PRIMER_BED} -p trimmed
samtools sort -@ [THREADS] -o {TRIMMED_SORTED_BAM} trimmed.bam
```

### FPGA Commands
I didn't see the per-sample read-mapping + trimming commands from Tianqi for the per-sample viral demo. Tianqi/Behnam should know the correct commands.

### Descriptions of Files and CLI Options
The only output file we care about is `{TRIMMED_SORTED_BAM}`. All other intermediate files can be deleted when the plugin finishes.
* `[THREADS]` = Number of threads (e.g. 16)
* `{REF_GENOME_FASTA}` = The reference genome FASTA ([`NC_045512.2`](https://github.com/niemasd/SD-COVID-Sequencing/blob/main/reference_genome/NC_045512.2.fas) for the demo)
* `{READS_R1}` = The R1 reads FASTQ
* `{READS_R2}` = The R2 reads FASTQ
* `{PRIMER_BED}` = Primer BED file ([`sarscov2_v2_primers.bed`](https://github.com/niemasd/SD-COVID-Sequencing/blob/main/primers/swift/sarscov2_v2_primers.bed) for the demo)
* `{TRIMMED_SORTED_BAM}` = The output trimmed sorted BAM

## Plugin 2: Call Consensus Genome
Call a consensus genome sequence from the trimmed BAM.

### CPU Commands
```bash
samtools mpileup -A -aa -d 0 -Q 0 --reference {REF_GENOME_FASTA} {TRIMMED_SORTED_BAM} | ivar consensus -m 10 -t 0.5 -n N -p {CONSENSUS_FASTA}
```

### FPGA Commands
None: we have not accelerated consensus sequence calling yet.

### Descriptions of Files and CLI Options
The only output file we care about is `{CONSENSUS_FASTA}`. All other intermediate files can be deleted when the plugin finishes.
* `{REF_GENOME_FASTA}` = The reference genome FASTA ([`NC_045512.2`](https://github.com/niemasd/SD-COVID-Sequencing/blob/main/reference_genome/NC_045512.2.fas) for the demo)
* `{TRIMMED_SORTED_BAM}` = The output of Plugin 1
* `{CONSENSUS_FASTA}` = The output consensus genome sequence for this sample
