# EasySci Pipeline

Processing pipeline for EasySci, adapted from Junyue Cao's lab (see their [repository](https://github.com/JunyueCaoLab/EasySci) and [paper](https://www.nature.com/articles/s41588-023-01572-y)). This pipeline handles the complete workflow from raw FASTQ files to gene expression count matrices. 

**Note**: code is under active development, some features may not be fully implemented or may change.

## Overview

This pipeline processes the sequencing data through multiple stages:
1. **Conda environment setup**
2. **Barcode extraction and demultiplexing**
3. **Adapter trimming**
4. **Read alignment**
5. **BAM filtering**
6. **BAM tagging with cell barcodes and UMIs**
7. **UMI deduplication**
8. **SLAM-seq nascent read filtering** (SLAM-seq experiments only)
9. **Transcriptome index creation**
10. **Read counting and gene assignment**
11. **Count matrix merging across samples**

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/denis-torre/easysci.git
cd easysci
```

### 2. Create Conda Environment
Create the conda environment from the provided YAML file. This will install all required dependencies including Python, BioPython, pysam, HTSeq, STAR, samtools, and other tools:

```bash
conda env create -f easysci-env.yaml
```

This will create an environment named `easysci_env` with all necessary dependencies.

### 3. Activate the Environment
```bash
conda activate easysci_env
```

### 4. Make the Wrapper Script Executable
```bash
chmod +x easysci.sh
```

### 5. (Optional) Add to PATH for System-wide Access
```bash
ln -s $(pwd)/easysci.sh ~/.local/bin/easysci
```

## Pipeline Steps

### Step 1: Barcode Extraction and Demultiplexing

Extract cell barcodes and UMIs from sequencing reads and generate barcoded FASTQ files. 

**Note**: FASTQ files used at this stage should already be demultiplexed using the P7 index barcodes.

**Input:**
- Three FASTQ files from paired-end sequencing:
    - Ligation barcode read (usually the "R1.fastq.gz" file)
    - UMI + cell barcode read (usually the "R2.fastq.gz" file)
    - Biological RNA read (usually the "R3.fastq.gz" file)
- EasySci barcode reference files

**Command:**
```bash
python easysci.py barcode \
  --ligation_barcode_read_file /path/to/ligation_R1.fastq.gz \
  --umi_cell_barcode_read_file /path/to/cell_barcode_R2.fastq.gz \
  --biological_read_file /path/to/biological_R3.fastq.gz \
  --ligation_barcode_file barcodes/ligation_barcodes.tsv \
  --RT_barcode_file barcodes/RT_barcodes.tsv \
  --randomN_barcode_file barcodes/RT_randomN_barcodes.txt \
  --output_basename /path/to/output/sample1
```

**Output:**
- `sample1.R1.fastq.gz` - Barcoded R1 reads
- `sample1.R2.fastq.gz` - Barcoded R2 reads

### Step 2: Adapter Trimming

Trim adapters and poly-A tails from the barcoded reads using Trim Galore.

**Command:**
```bash
trim_galore \
  --paired \
  /path/to/output/sample1.R1.fastq.gz \
  /path/to/output/sample1.R2.fastq.gz \
  -a2 AAAAAAAA \
  --three_prime_clip_R1 1 \
  --stringency 3 \
  --cores 4 \
  -o /path/to/trimmed/
```

**Parameters:**
- `--paired`: Process paired-end reads
- `-a2 AAAAAAAA`: Trim poly-A sequences from R2 (biological read)
- `--three_prime_clip_R1 1`: Remove 1 base from 3' end of R1
- `--stringency 3`: Adapter overlap stringency
- `--cores 4`: Number of cores for parallel processing

**Output:**
- `sample1.R1_val_1.fq.gz` - Trimmed R1 reads
- `sample1.R2_val_2.fq.gz` - Trimmed R2 reads
- Trimming reports

### Step 3: Alignment

Align the trimmed FASTQ files to the reference genome using STAR.

**Command:**
```bash
STAR --runThreadN 8 \
  --genomeDir /path/to/STAR_index \
  --readFilesIn /path/to/trimmed/sample1.R2_val_2.fq.gz /path/to/trimmed/sample1.R1_val_1.fq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix /path/to/aligned/sample1_ \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM MD
```

> **Note (SLAM-seq only):** The `MD` attribute in `--outSAMattributes` is required for the `slam` subcommand. It encodes reference bases at each aligned position, enabling T→C detection without an external reference FASTA. If you have a BAM aligned without `MD`, add it retroactively: `samtools calmd -b input.bam reference.fa > input_with_md.bam`

**Output:**
- `sample1_Aligned.sortedByCoord.out.bam` - Aligned BAM file

### Step 4: BAM Filtering

Filter BAM file to retain only high-quality, properly paired reads.

**Command:**
```bash
samtools view -bh -q 30 -f 2 -F 1804 -@ 32 \
  -o /path/to/filtered/sample1_filtered.bam \
  /path/to/aligned/sample1_Aligned.sortedByCoord.out.bam
```

**Parameters:**
- `-bh`: Output BAM format with header
- `-q 30`: Keep reads with mapping quality ≥ 30
- `-f 2`: Keep properly paired reads
- `-F 1804`: Exclude unmapped, not primary alignment, fails QC, and duplicate reads
- `-@ 32`: Use 32 threads

**Output:**
- `sample1_filtered.bam` - Filtered BAM file with high-quality reads

### Step 5: BAM Tagging

Add cell barcode (CB) and UMI (UR) tags to the BAM file from the read names. This is required for UMI deduplication and counting.

**Command:**
```bash
easysci tag \
  --input_file /path/to/filtered/sample1_filtered.bam \
  --output_file /path/to/tagged/sample1_tagged.bam \
  --num_processes 4
```

**Output:**
- `sample1_tagged.bam` - BAM file with CB and UR tags added

### Step 6: UMI Deduplication

Deduplicate reads based on UMI tags per cell to remove PCR duplicates.

**Command:**
```bash
umi_tools dedup \
  --extract-umi-method tag \
  --cell-tag=CB \
  --umi-tag=UR \
  --per-cell \
  --unpaired-reads discard \
  --chimeric-pairs use \
  --paired \
  -I /path/to/tagged/sample1_tagged.bam \
  -S /path/to/dedup/sample1_dedup_temp.bam && \
samtools sort -n /path/to/dedup/sample1_dedup_temp.bam \
  -o /path/to/dedup/sample1_dedup.bam -@ 32 && \
rm /path/to/dedup/sample1_dedup_temp.bam
```

**Parameters:**
- `--extract-umi-method tag`: Extract UMI from BAM tags
- `--cell-tag=CB`: Cell barcode tag name
- `--umi-tag=UR`: UMI tag name
- `--per-cell`: Deduplicate separately for each cell
- `--unpaired-reads discard`: Discard reads without pairs
- `--chimeric-pairs use`: Use chimeric read pairs
- `--paired`: Process paired-end reads
- Final step sorts by read name and removes temporary file

**Output:**
- `sample1_dedup.bam` - Deduplicated BAM file sorted by name

### Step 8: SLAM-seq Nascent Read Filtering (SLAM-seq experiments only)

Identify and extract reads carrying the T→C substitution signature introduced by s⁴U metabolic labeling. This step operates on the deduplicated, name-sorted BAM from Step 7 and produces a nascent-only BAM that can be passed to `count` (Step 10) separately from the total-RNA BAM.

**Requires:** the `MD` BAM tag must be present (see Step 3 note).

**Algorithm:** For each read with FLAG ∈ {83, 99, 147, 163}, the tool counts (1) all non-SNP mismatches above the base-quality threshold, and (2) T→C mismatches on + strand reads (FLAG 99/163) or A→G mismatches on − strand reads (FLAG 83/147 — the reverse-complement representation of T→C). If T→C mismatches / total mismatches ≥ `--min_tc_ratio`, the read's QNAME is marked nascent. In a second pass, all reads (both mates of a pair) with a nascent QNAME are written to the output BAM.

**Input:**
- `sample1_dedup.bam` — deduplicated, name-sorted BAM from Step 7
- (Optional) background SNP VCF to suppress germline variants

**Command:**
```bash
easysci slam \
  --input_bam /path/to/dedup/sample1_dedup.bam \
  --output_bam /path/to/nascent/sample1_nascent.bam \
  --snp_vcf /path/to/ref_SNPs.vcf \
  --min_base_quality 45 \
  --min_tc_ratio 0.3
```

**Parameters:**
- `--input_bam`: deduplicated, name-sorted BAM (required)
- `--output_bam`: output path for the nascent BAM (required)
- `--snp_vcf`: tab-separated SNP file with columns CHROM, POS (1-based), REF, ALT; lines starting with `#` are skipped. VarScan output (with `Chrom/Position/Ref/Var` header) is accepted as-is. Standard VCF: pre-process with `bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n'` (optional)
- `--min_base_quality`: Phred quality floor (default 45, matching the Cao lab reference)
- `--min_tc_ratio`: T→C fraction threshold (default 0.3, matching the Cao lab reference)
- `--num_processes`: reserved for future parallelism; currently 1 process only

**Output:**
- `sample1_nascent.bam` — nascent-only BAM, name-sorted, ready for `easysci count`
- `sample1_nascent-slam-summary.tsv` — per-category read counts (total processed, skipped by flag/MD/quality, no-mismatch reads, below-threshold reads, nascent reads classified, reads written)

### Step 9: Create Transcriptome Index

Generate a transcriptome annotation index for efficient read counting (only needs to be done once per transcriptome).

**Command:**
```bash
easysci index \
    --input_gtf /path/to/genome_annotation.gtf \
    --output_file /path/to/indices/transcriptome_index.pickle
```

**Output:**
- `transcriptome_index.pickle` - Transcriptome annotation index

### Step 10: Count Reads

Count reads and assign them to genes, generating a cell-by-gene count matrix.

**Command:**
```bash
easysci count \
  --input_bam /path/to/dedup/sample1_dedup.bam \
  --index_file /path/to/indices/transcriptome_index.pickle \
  --output_dir /path/to/counts/ \
  --sample_name sample1 \
  --multigene_reads closest_TES \
  --primer_type shortdT \
  --library_strandedness reverse \
  --count_introns \
  --collapse_by_gene_name
```

**Parameters:**
- `--multigene_reads`: How to handle reads mapping to multiple genes:
  - `discard`: Discard reads mapping to multiple genes (default) ✅ **Implemented**
  - `closest_TES`: Assign to gene with closest transcription end site (recommended for dT libraries) ✅ **Implemented**
  - `uniform`: Distribute uniformly across genes ❌ **Not yet implemented**
  - `prop_unique`: Proportional to unique read counts ❌ **Not yet implemented**
  - `EM`: Expectation-maximization algorithm ❌ **Not yet implemented**
- `--primer_type`: `shortdT` (oligo-dT) or `randomN` (random hexamer)
- `--library_strandedness`: `reverse` (dUTP) or `none` (unstranded)
- `--count_introns`: Include intronic reads (useful for nuclear RNA)
- `--collapse_by_gene_name`: Collapse transcript-level counts to gene-level

**Output:**
Two TSV (tab-separated) files in the specified output directory:
- `{sample_name}-read_assignment.tsv` - Per-read assignments with columns: barcode, umi, gene_id, mapped_feature, mapping
- `{sample_name}-unique_umi_count.tsv` - UMI counts with columns: barcode, gene_id (or gene_name if collapsed), unique_umi_count

**Note**: Output format will be optimized in future versions.

### Step 11: Merge Count Matrices (Optional)

Merge count matrices from multiple samples or replicates into a single sparse matrix.

**Command:**
```bash
easysci merge \
  --count_files /path/to/counts/sample1-unique_umi_count.tsv,/path/to/counts/sample2-unique_umi_count.tsv \
  --output_basename /path/to/merged/combined_
```

**Note:** The output_basename should end with an underscore or path separator as filenames will be appended to it.

**Output:**
Four files in Matrix Market format suitable for downstream analysis:
- `combined_gene_counts.mtx` - Sparse count matrix in Matrix Market format
- `combined_cell_annotation.tsv` - Cell metadata (barcode, genes_detected, total_counts)
- `combined_gene_annotation.tsv` - Gene metadata (gene_id/gene_name, cells_expressed)
- `combined_barcode_rank_plot.png` - Quality control visualization showing UMI distribution

## Complete Example Workflow

```bash
# Set up directories
mkdir -p output/{barcoded,trimmed,aligned,filtered,tagged,dedup,counts,indices}

# Step 1: Barcode extraction
easysci barcode \
  --ligation_barcode_read_file raw/ligation_R1.fastq.gz \
  --umi_cell_barcode_read_file raw/cell_R2.fastq.gz \
  --biological_read_file raw/bio_R3.fastq.gz \
  --ligation_barcode_file barcodes/ligation_barcodes.tsv \
  --RT_barcode_file barcodes/RT_barcodes.tsv \
  --randomN_barcode_file barcodes/RT_randomN_barcodes.txt \
  --output_basename output/barcoded/sample1

# Step 2: Adapter trimming
trim_galore \
  --paired \
  output/barcoded/sample1.R1.fastq.gz \
  output/barcoded/sample1.R2.fastq.gz \
  -a2 AAAAAAAA \
  --three_prime_clip_R1 1 \
  --stringency 3 \
  --cores 4 \
  -o output/trimmed/

# Step 3: Alignment with STAR
STAR --runThreadN 8 \
  --genomeDir /path/to/STAR_index \
  --readFilesIn output/trimmed/sample1.R2_val_2.fq.gz output/trimmed/sample1.R1_val_1.fq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix output/aligned/sample1_ \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM MD

# Step 4: Filter BAM file
samtools view -bh -q 30 -f 2 -F 1804 -@ 32 \
  -o output/filtered/sample1_filtered.bam \
  output/aligned/sample1_Aligned.sortedByCoord.out.bam

# Step 5: Add tags
easysci tag \
  --input_file output/filtered/sample1_filtered.bam \
  --output_file output/tagged/sample1_tagged.bam \
  --num_processes 4

# Step 6: Deduplicate UMIs
umi_tools dedup \
  --extract-umi-method tag \
  --cell-tag=CB \
  --umi-tag=UR \
  --per-cell \
  --unpaired-reads discard \
  --chimeric-pairs use \
  --paired \
  -I output/tagged/sample1_tagged.bam \
  -S output/dedup/sample1_dedup_temp.bam && \
samtools sort -n output/dedup/sample1_dedup_temp.bam \
  -o output/dedup/sample1_dedup.bam -@ 32 && \
rm output/dedup/sample1_dedup_temp.bam

# Step 8: Filter nascent reads (SLAM-seq only)
easysci slam \
  --input_bam output/dedup/sample1_dedup.bam \
  --output_bam output/nascent/sample1_nascent.bam \
  --snp_vcf /path/to/ref_SNPs.vcf

# Step 9: Create index (once per genome)
easysci index \
  --input_gtf /path/to/genes.gtf \
  --output_file output/indices/transcriptome_index.pickle

# Step 8: Count reads
easysci count \
  --input_bam output/dedup/sample1_dedup.bam \
  --index_file output/indices/transcriptome_index.pickle \
  --output_dir output/counts/ \
  --sample_name sample1 \
  --multigene_reads closest_TES \
  --primer_type shortdT \
  --library_strandedness reverse \
  --count_introns \
  --collapse_by_gene_name
```

## Citation

If you use this pipeline, please cite the original EasySci publication: [A global view of aging and Alzheimer’s pathogenesis-associated cell population dynamics and molecular signatures in human and mouse brains](https://www.nature.com/articles/s41588-023-01572-y).

## Support

For issues and questions, please open an issue on the GitHub repository.