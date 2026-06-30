#############################################
########## 1. Import Libraries
#############################################
import os
import sys
import argparse

# Custom
sys.path.append(os.path.join(os.path.dirname(__file__), "scripts"))
from scripts import barcode, tag, index, count, merge, slam, gRNA, barcode_correction
# from scripts import count, index

#############################################
########## 2. Argument parser
#############################################

def run():
    
    # Initialize parsers
    parser = argparse.ArgumentParser(description='Main script with subcommands')
    subparsers = parser.add_subparsers(dest='command', help='Available subcommands')

    # 1. Subparser for barcoding
    barcode_parser = subparsers.add_parser('barcode', help='Barcode EasySci reads for paired-end sequencing.')
    barcode_parser.add_argument("--ligation_barcode_read_file", type=str, required=True, help="File path for the ligation barcode read file.")
    barcode_parser.add_argument("--umi_cell_barcode_read_file", type=str, required=True, help="File path for the UMI and cell barcode + biological read file.")
    barcode_parser.add_argument("--biological_read_file", type=str, required=True, help="File path for the biological read file.")
    barcode_parser.add_argument("--output_basename", type=str, required=True, help="Basename for output files; barcoded R1 and R2 outputs will be generated based on this name.")
    barcode_parser.add_argument("--ligation_barcode_file", type=str, required=True, help="File path containing ligation barcodes.")
    barcode_parser.add_argument("--RT_barcode_file", type=str, required=True, help="File path containing RT barcodes.")
    barcode_parser.add_argument("--randomN_barcode_file", type=str, required=True, help="File path containing randomN barcodes.")

    # 2. Subparser for BAM tagging
    tagging_parser = subparsers.add_parser('tag', help="Add CB and UB tags from read name to a BAM file.")
    tagging_parser.add_argument("--input_file", type=str, required=True, help="Input BAM file.")
    tagging_parser.add_argument("--output_file", type=str, required=True, help="Output BAM file.")
    tagging_parser.add_argument("--num_processes", type=int, default=1, help="Number of processes to use.")    

    # 3. Subparser for index creation
    index_parser = subparsers.add_parser('index', help='Create a genome index for rapid counting.')
    index_parser.add_argument("--input_gtf", type=str, required=True, help="Input GTF file.")
    index_parser.add_argument("--output_file", type=str, required=True, help="Output file.")

    # 4. Subparser for counting
    count_parser = subparsers.add_parser('count', help='Count reads from a BAM file.')
    count_parser.add_argument("--input_bam", type=str, required=True, help="Input BAM file.")
    count_parser.add_argument("--index_file", type=str, required=True, help="EasySci index file. If unstranded, use the unstranded index.")
    count_parser.add_argument("--output_dir", type=str, required=True, help="Output basename.")
    count_parser.add_argument("--sample_name", type=str, required=True, help="Sample name for file and cell naming.")
    count_parser.add_argument("--multigene_reads", type=str, required=False, default='discard', choices=['discard', 'closest_TES', 'uniform', 'prop_unique', 'EM'], help="How to handle reads that map to multiple genes.")
    count_parser.add_argument("--primer_type", type=str, required=False, default='shortdT', help="Primer type used for library preparation (shortdT or randomN).")
    count_parser.add_argument("--library_strandedness", type=str, required=False, default='reverse', choices=['reverse', 'none'], help="Strandedness of the library.")
    count_parser.add_argument("--read_subset", type=int, required=False, default=None, help="Number of reads to subset for testing.")
    count_parser.add_argument("--count_introns", action='store_true', help="Count introns.")
    count_parser.add_argument("--collapse_by_gene_name", action='store_true', help="Collapse by gene name.")

    # 5. Subparser for index creation
    merge_parser = subparsers.add_parser('merge', help='Create a genome index for rapid counting.')
    merge_parser.add_argument("--count_files", type=str, required=True, help="Input count files.")
    merge_parser.add_argument("--output_basename", type=str, required=True, help="Output basename.")

    # 6. Subparser for SLAM-seq nascent read filtering
    slam_parser = subparsers.add_parser('slam', help='Filter nascent reads from a SLAM-seq BAM using the T→C substitution signature.')
    slam_parser.add_argument("--input_bam", type=str, required=True, help="Input deduplicated, name-sorted BAM file (output of umi_tools dedup + samtools sort -n).")
    slam_parser.add_argument("--output_bam", type=str, required=True, help="Output BAM file containing only nascent (T→C-labeled) reads.")
    slam_parser.add_argument("--snp_vcf", type=str, required=False, default=None, help="Tab-separated file of background SNPs to exclude. Columns: CHROM(0) POS(1, 1-based) REF(2) ALT(3). Lines starting with '#' are skipped. VarScan format accepted directly; standard VCF: pre-process with 'bcftools query -f \"%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n\"'.")
    slam_parser.add_argument("--min_base_quality", type=int, required=False, default=45, help="Minimum Phred base quality for a position to be scored (default: 45).")
    slam_parser.add_argument("--min_tc_ratio", type=float, required=False, default=0.3, help="Minimum T→C conversion rate (tc_mismatches / total_T_positions_above_quality) to classify a read as nascent (default: 0.3).")
    slam_parser.add_argument("--min_tc_count", type=int, required=False, default=1, help="Minimum number of T→C mismatches required to classify a read as nascent (default: 1).")
    slam_parser.add_argument("--num_processes", type=int, required=False, default=1, help="Number of processes (default: 1; parallelism not yet implemented).")

    # 7. Subparser for gRNA counting
    grna_parser = subparsers.add_parser('gRNA', help='Count sgRNA UMIs per cell from PerturbSci-Kinetics FASTQ files.')
    grna_parser.add_argument("--r1", type=str, required=True, help="R1 FASTQ (UMI + RT barcode + constant region 1).")
    grna_parser.add_argument("--r2", type=str, required=True, help="R2 FASTQ (inner i7 barcode + constant region 2 + sgRNA sequence).")
    grna_parser.add_argument("--r3", type=str, required=True, help="R3 FASTQ (ligation barcode).")
    grna_parser.add_argument("--ligation_barcode_file", type=str, required=True, help="TSV of ligation barcodes (barcode_1bp_substitution / original_barcode).")
    grna_parser.add_argument("--RT_barcode_file", type=str, required=True, help="TSV of RT barcodes (barcode_1bp_substitution / original_barcode).")
    grna_parser.add_argument("--inner_i7_barcode_file", type=str, required=True, help="TSV of inner i7 barcodes (barcode_1bp_substitution / original_barcode).")
    grna_parser.add_argument("--sgrna_barcode_file", type=str, required=True, help="TSV of sgRNA sequences (barcode_1bp_substitution / original_barcode).")
    grna_parser.add_argument("--sgrna_annotation_file", type=str, required=False, default=None, help="Optional TSV mapping sgRNA sequences to names (gRNA_seq / names). If omitted, raw sequences are used as identifiers.")
    grna_parser.add_argument("--output_dir", type=str, required=True, help="Output directory.")
    grna_parser.add_argument("--sample_name", type=str, required=True, help="Sample name prefix for output files and cell IDs.")
    grna_parser.add_argument("--min_umi_threshold", type=int, required=False, default=1, help="Minimum total sgRNA UMIs per cell to retain (default: 1).")
    grna_parser.add_argument("--n_reads", type=int, required=False, default=None, help="If set, process only the first N reads (useful for testing).")

    # 8. Subparser for gRNA sequence discovery
    grna_discover_parser = subparsers.add_parser('gRNA-discover', help='Single-pass discovery of every sgRNA-position sequence passing the barcode/constant-region filter chain — no whitelisting, no UMI counting.')
    grna_discover_parser.add_argument("--r1", type=str, required=True, help="R1 FASTQ (UMI + RT barcode + constant region 1).")
    grna_discover_parser.add_argument("--r2", type=str, required=True, help="R2 FASTQ (inner i7 barcode + constant region 2 + sgRNA sequence).")
    grna_discover_parser.add_argument("--r3", type=str, required=True, help="R3 FASTQ (ligation barcode).")
    grna_discover_parser.add_argument("--ligation_barcode_file", type=str, required=True, help="TSV of ligation barcodes (barcode_1bp_substitution / original_barcode).")
    grna_discover_parser.add_argument("--RT_barcode_file", type=str, required=True, help="TSV of RT barcodes (barcode_1bp_substitution / original_barcode).")
    grna_discover_parser.add_argument("--inner_i7_barcode_file", type=str, required=True, help="TSV of inner i7 barcodes (barcode_1bp_substitution / original_barcode).")
    grna_discover_parser.add_argument("--output_dir", type=str, required=True, help="Output directory.")
    grna_discover_parser.add_argument("--sample_name", type=str, required=True, help="Sample name prefix for output files.")
    grna_discover_parser.add_argument("--n_reads", type=int, required=False, default=None, help="If set, process only the first N reads (useful for testing).")

    # 9. Subparser for building a 1bp-mismatch correction dict from a barcode whitelist
    correction_parser = subparsers.add_parser('build-correction-dict', help='Build a barcode_1bp_substitution/original_barcode correction TSV from a plain barcode whitelist (for use as --sgrna_barcode_file, --ligation_barcode_file, etc).')
    correction_parser.add_argument("--input_file", type=str, required=True, help="TSV containing a column of whitelist barcode/sgRNA sequences.")
    correction_parser.add_argument("--sequence_column", type=str, required=True, help="Name of the column in --input_file containing the sequences to expand.")
    correction_parser.add_argument("--output_file", type=str, required=True, help="Output TSV path (columns: barcode_1bp_substitution, original_barcode).")

    # Get args
    args = parser.parse_args()

    if args.command == 'barcode':
        barcode.main(args)
    elif args.command == 'tag':
        tag.main(args)
    elif args.command == 'index':
        index.main(args)
    elif args.command == 'count':
        count.main(args)
    elif args.command == 'merge':
        merge.main(args)
    elif args.command == 'slam':
        slam.main(args)
    elif args.command == 'gRNA':
        gRNA.main(args)
    elif args.command == 'gRNA-discover':
        gRNA.main_discover(args)
    elif args.command == 'build-correction-dict':
        barcode_correction.main(args)
    else:
        parser.print_help()

#############################################
########## 4. Run
#############################################
    
if __name__ == '__main__':
    run()
    
    

# # 2. Subparser for deduplication
# deduplication_parser = subparsers.add_parser('deduplicate', help='Remove duplicates from a BAM file.')
# deduplication_parser.add_argument("--input_bam", type=str, required=True, help="Input BAM file.")
# deduplication_parser.add_argument("--output_bam", type=str, required=True, help="Output BAM file.")
# deduplication_parser.add_argument("--version", type=str, required=True, help="Version of deduplication algorithm.")


# # 3. Subparser for creating an index
# index_parser = subparsers.add_parser('index', help='Create an EasySci index for quantification.')
# index_parser.add_argument("--gtf_file", type=str, required=True, help="GTF file.")
# index_parser.add_argument("--output_dir", type=str, required=True, help="Output directory.")

# # 4. Subparser for counting
# count_parser = subparsers.add_parser('count', help='Count reads from a BAM file.')
# count_parser.add_argument("--input_bam", type=str, required=True, help="Input BAM file.")
# # count_parser.add_argument("--index_dir", type=str, required=True, help="Directory of EasySci index.")
# # count_parser.add_argument("--randomN_barcode_file", type=str, required=True, help="File path containing randomN barcodes.")
# count_parser.add_argument("--output_dir", type=str, required=True, help="Output directory.")
# count_parser.add_argument("--sample_name", type=str, required=False, help="Sample name.")


# # 5. Subparser for count merging
# merge_parser = subparsers.add_parser('merge', help='Merge counts from multiple replicates.')
# merge_parser.add_argument("--input_dirs", type=str, required=True, help="Input count files.")
# merge_parser.add_argument("--output_dir", type=str, required=True, help="Output directory.")
# merge_parser.add_argument("--RT_matching_file", type=str, required=True, help="File containing RT barcode matching information.")
# merge_parser.add_argument("--index_dir", type=str, required=True, help="Directory of EasySci index.")


# if args.command == 'barcode':
#     barcode.main(args)
# elif args.command == 'deduplicate':
#     deduplicate.main(args)
# elif args.command == 'index':
#     index.main(args)
# elif args.command == 'count':
#     count.main(args)
# elif args.command == 'merge':
#     merge.main(args)
# else:
#     parser.print_help()
