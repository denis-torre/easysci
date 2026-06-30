#############################################
########## 1. Import Libraries
#############################################
import os
import time
from datetime import timedelta
import pysam

#############################################
########## 2. Constants
#############################################

VALID_FLAGS  = {83, 99, 147, 163}
PLUS_FLAGS   = {99, 163}   # plus-strand reads: T→C is the s⁴U signature
MINUS_FLAGS  = {83, 147}   # minus-strand reads: A→G is the RC representation of T→C

#############################################
########## 3. Processing functions
#############################################

def load_snps(snp_vcf):
    """Load background SNPs into a set of (chrom, pos_1based, ref, alt) tuples.

    Accepts VarScan format (Chrom/Position/Ref/Var header, no '#') and
    4-column TSV produced by: bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n'
    Lines starting with '#' and non-numeric POS are skipped.
    """
    snp_set = set()
    with open(snp_vcf, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            snp_set.add((fields[0], pos, fields[2].upper(), fields[3].upper()))
    return snp_set


def filter_nascent_reads(input_bam, output_bam, snp_vcf, min_base_quality, min_tc_ratio, min_tc_count, num_processes):

    # Parallelism stub
    if num_processes > 1:
        print(f"Warning: --num_processes {num_processes} requested but parallelism is not yet "
              f"implemented; running with 1 process.", flush=True)

    ### 1. Load SNPs
    snp_set = set()
    if snp_vcf:
        print(f"Loading SNPs from {snp_vcf}...", flush=True)
        snp_set = load_snps(snp_vcf)
        print(f"Loaded {len(snp_set):,} SNPs.", flush=True)

    ### 2. Create output directory
    os.makedirs(os.path.dirname(os.path.abspath(output_bam)), exist_ok=True)

    ### 3. Initialise counters
    counters = {
        'total_reads_processed':    0,
        'skipped_wrong_flag':       0,
        'skipped_no_md_tag':        0,
        'skipped_no_qualities':     0,
        'no_mismatches':            0,
        'tc_ratio_below_threshold': 0,
        'nascent_reads_classified': 0,
        'total_reads_written':      0,
    }
    nascent_qnames = set()

    ### 4. Pass 1 — classify reads, collect nascent QNAMEs
    print("Pass 1: Scanning reads for T→C substitution signature...", flush=True)
    start_time = time.time()
    input_bam_fh = pysam.AlignmentFile(input_bam, "rb")
    try:
        for i, read in enumerate(input_bam_fh):
            counters['total_reads_processed'] += 1

            if i % 500000 == 0 and i > 0:
                elapsed = str(timedelta(seconds=int(time.time() - start_time)))
                print(f"Processing read {i:,}... ({elapsed} elapsed since previous cycle)", flush=True)
                start_time = time.time()

            # FLAG filter — only properly paired reads in canonical orientations
            if read.flag not in VALID_FLAGS:
                counters['skipped_wrong_flag'] += 1
                continue

            # MD tag guard — get_aligned_pairs(with_seq=True) requires MD
            if not read.has_tag('MD'):
                counters['skipped_no_md_tag'] += 1
                continue

            # Base quality guard
            if read.query_qualities is None:
                counters['skipped_no_qualities'] += 1
                continue

            chrom = read.reference_name
            total_t_positions = 0
            tc_mismatches = 0

            for query_pos, ref_pos, ref_base in read.get_aligned_pairs(
                    matches_only=True, with_seq=True):

                # Base quality filter (strictly greater than threshold)
                if read.query_qualities[query_pos] <= min_base_quality:
                    continue

                ref_up = ref_base.upper()
                qry_up = read.query_sequence[query_pos].upper()

                # Count T positions above quality threshold (strand-aware), matches and mismatches alike
                if (read.flag in PLUS_FLAGS and ref_up == 'T') or \
                   (read.flag in MINUS_FLAGS and ref_up == 'A'):
                    total_t_positions += 1

                # Skip exact matches
                if qry_up == ref_up:
                    continue

                # SNP filter (ref_pos is 0-based; VCF positions are 1-based)
                if snp_set and (chrom, ref_pos + 1, ref_up, qry_up) in snp_set:
                    continue

                # T→C signature (strand-aware)
                if read.flag in PLUS_FLAGS:
                    if ref_up == 'T' and qry_up == 'C':
                        tc_mismatches += 1
                else:  # MINUS_FLAGS
                    if ref_up == 'A' and qry_up == 'G':
                        tc_mismatches += 1

            # Reads with no T positions above threshold cannot have a T→C ratio
            if total_t_positions == 0:
                counters['no_mismatches'] += 1
                continue

            tc_ratio = tc_mismatches / total_t_positions
            if tc_mismatches >= min_tc_count and tc_ratio >= min_tc_ratio:
                nascent_qnames.add(read.query_name)
                counters['nascent_reads_classified'] += 1
            else:
                counters['tc_ratio_below_threshold'] += 1

    finally:
        input_bam_fh.close()

    print(f"Pass 1 complete. {len(nascent_qnames):,} unique nascent QNAMEs identified.", flush=True)

    ### 5. Pass 2 — write nascent reads to output BAM
    print("Pass 2: Writing nascent reads to output BAM...", flush=True)
    start_time = time.time()
    n_written = 0
    input_bam_fh2 = pysam.AlignmentFile(input_bam, "rb")
    output_bam_fh = pysam.AlignmentFile(output_bam, "wb", template=input_bam_fh2)
    try:
        for i, read in enumerate(input_bam_fh2):
            if i % 1000000 == 0 and i > 0:
                elapsed = str(timedelta(seconds=int(time.time() - start_time)))
                print(f"Writing read {i:,}... ({elapsed} elapsed since previous cycle)", flush=True)
                start_time = time.time()

            if read.query_name in nascent_qnames:
                output_bam_fh.write(read)
                n_written += 1
    finally:
        input_bam_fh2.close()
        output_bam_fh.close()

    counters['total_reads_written'] = n_written
    print(f"Pass 2 complete. Wrote {n_written:,} reads to {output_bam}.", flush=True)

    ### 6. Write summary file
    summary_file = output_bam.replace('.bam', '-slam-summary.tsv')
    with open(summary_file, 'w') as f:
        f.write("read_class\tcount\n")
        for key, val in counters.items():
            f.write(f"{key}\t{val}\n")
    print(f"Summary written to {summary_file}", flush=True)
    print("Done!", flush=True)

#############################################
########## 4. Main
#############################################

def main(args):
    filter_nascent_reads(
        input_bam       = args.input_bam,
        output_bam      = args.output_bam,
        snp_vcf         = args.snp_vcf,
        min_base_quality= args.min_base_quality,
        min_tc_ratio    = args.min_tc_ratio,
        min_tc_count    = args.min_tc_count,
        num_processes   = args.num_processes,
    )

#############################################
########## 5. Run
#############################################

if __name__ == "__main__":
    main(args)
