These plugin descriptions are only for the demo. They may translate to long-term use as well, but they may need to be modified slightly (notes are left throughout regarding these potential improvements for long-term use). Much of this was adapted from [Tianqi's notes](https://git.ucsd.edu/tiz014/fpga-to-qiita-demo/snippets/33). Some notations I'll use in this document:

* `{...}` = file
* `[...]` = CLI option

# DeepVariant Pipeline
In this pipeline, we start with human reads (assumed to already be trimmed in the demo dataset), and we (1) map to the human reference genome, and then (2) call variants using DeepVariant.

## Input Files
* `ucsc_hg19.chr21.fa`: Chromosome 21 of the [hg19](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13) human reference genome
* `BGISEQ_PE100_NA12878.sorted.chr21.fataq`: DeepVariant’s training tutorial (short reads: HG001 sequenced by BGISEQ-500 from individual NA12878)
  * This is actually a FASTQ file; the `.fataq` extension is a typo

## Plugin 1: Map Reads to Human Reference Genome
### CPU Commands
```bash
bwa index {REF_GENOME_FASTA} # should probably be done before the plugin is executed
bwa mem -t [THREADS] -k 21 {REF_GENOME_FASTA} {READS} > output.sam

# everything below this is the same as FPGA
samtools sort -@ [THREADS] -o {SORTED_BAM} output.sam
samtools index -@ [THREADS] {SORTED_BAM}
```

### FPGA Commands
```bash
fpga_alg_dv -tables {TABLES_FOLDER} -r {READS} -o output.sam

# everything below this is the same as CPU
samtools sort -@ [THREADS] -o {SORTED_BAM} output.sam
samtools index -@ [THREADS] {SORTED_BAM}
```

### Descriptions of Files and CLI Options
The only output files we care about are `{SORTED_BAM}` and the index file produced by `samtools index`. All other intermediate files can be destroyed when the plugin finishes executing (and commands can be piped together if desired, but perhaps not for the demo to be able to get step-by-step runtimes if desired).

* `{REF_GENOME_FASTA}` = The reference genome FASTA (`ucsc_hg19.chr21.fa` for the demo)
* `{TABLES_FOLDER}` = Index files (similar to Minimap2's `.mmi` files) used by the FPGA mapper
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

# Woltka/UniFrac Pipeline
In this pipeline, we accelerate host-filtering as well as the to-database read mapping prior to calling Woltka. See [Qiita Recommended Default Workflows](https://qiita.ucsd.edu/workflows/) for details of full pipeline.

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
* `{TABLES_FOLDER}` = Index files (similar to Minimap2's `.mmi` files) used by the FPGA mapper
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
* `{TABLES_FOLDER}` = Index files (similar to Minimap2's `.mmi` files) used by the FPGA mapper
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
```bash
fpga_alg_covid -tables {TABLES_FOLDER} -r <(cat {READS_R1} {READS_R2}) -o untrimmed.sam
fpga_ivar_trim -q 20 -w 4 -m 30 -x 5 -tables ./FOLDER_OF_TABLES/ -r untrimmed.sam -o trimmed.sam # this seems incorrect; input should be the SAM/BAM file
samtools sort -@ [THREADS] -o {TRIMMED_SORTED_BAM} trimmed.sam
```

### Descriptions of Files and CLI Options
The only output file we care about is `{TRIMMED_SORTED_BAM}`. All other intermediate files can be deleted when the plugin finishes.
* `[THREADS]` = Number of threads (e.g. 16)
* `{REF_GENOME_FASTA}` = The reference genome FASTA ([`NC_045512.2`](https://github.com/niemasd/SD-COVID-Sequencing/blob/main/reference_genome/NC_045512.2.fas) for the demo)
* `{READS_R1}` = The R1 reads FASTQ
* `{READS_R2}` = The R2 reads FASTQ
* `{TABLES_FOLDER}` = Index files (similar to Minimap2's `.mmi` files) used by the FPGA mapper
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

# Multi-Sample COVID-19 Pipeline
## Plugin 1: Multiple Sequence Alignment
The input to this plugin is the output of Per-Sample COVID-19 Plugin 2 from many samples (i.e., many consensus genome sequences), and this plugin concatenates them and performs reference-guided Multiple Sequence Alignment.

### CPU Commands
```bash
cat {FASTA_1} {FASTA_2} {FASTA_N} > unaligned.fas
python3 ViralMSA.py -s unaligned.fas -r SARS-CoV-2 -e {EMAIL_ADDRESS} -o viralmsa_out -t [THREADS]
python3 trim_msa.py -i viralmsa_out/*.aln -s 100 -e 50 -o {TRIMMED_MSA}
```

### FPGA Commands
None: we have not accelerated Multiple Sequence Alignment (it's plenty fast).

### Descriptions of Files and CLI Options
The only output file we care about is `{TRIMMED_MSA}`. All other intermediate files can be deleted when the plugin finishes.
* `{FASTA_1}`, `{FASTA_2}`, ..., `{FASTA_N}` = The consensus genome sequences to align (outputs from multiple calls to Per-Sample COVID-19 Pipeline Plugin 2)
* `{EMAIL_ADDRESS}` = User's email address (the reference genome is pulled from NCBI using BioPython)
* `[THREADS]` = Number of threads (e.g. 16)
* `{TRIMMED_MSA}` = The trimmed MSA to use for downstream analyses

## Plugin 2: Pairwise Distances
This plugin computes all pairwise TN93 distances from the given trimmed MSA.

### CPU Commands
```bash
cat {TRIMMED_MSA} | tn93 -t 1 -l 1 > {TN93_CSV}
```

### GPU Commands
```bash
fpga_tn93 -input {TRIMMED_MSA} -o {TN93_CSV}
```

### Descriptions of Files and CLI Options
* `{TRIMMED_MSA}` = The trimmed MSA output by Plugin 1
* `{TN93_CSV}` = The output TN93 pairwise distances

# `Dockerfile` List for Programs
To aid with installing the various tools needed for these plugins, I've created `Dockerfile`s for the mentioned tools, which contain specific commands needed to compile + install the tools.

* **BWA:** https://github.com/Niema-Docker/bwa/blob/main/Dockerfile
* **Samtools:** https://github.com/Niema-Docker/samtools/blob/main/Dockerfile
* **Bowtie2:** https://github.com/Niema-Docker/bowtie2/blob/main/Dockerfile
* **Pigz:** `sudo apt-get install -y pigz`
* **Minimap2:** https://github.com/Niema-Docker/minimap2/blob/main/Dockerfile
* **iVar:** https://github.com/Niema-Docker/ivar/blob/main/Dockerfile
* **ViralMSA:** https://github.com/niemasd/ViralMSA#installation
* **trim_msa.py:** https://github.com/niemasd/tools/blob/master/trim_msa.py
* **tn93:** https://github.com/Niema-Docker/tn93/blob/main/Dockerfile
