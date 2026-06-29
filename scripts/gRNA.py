#############################################
########## 1. Import Libraries
#############################################
import gzip
import os
import time
from datetime import timedelta
import pandas as pd
import argparse
from Bio import SeqIO

#############################################
########## 2. Processing functions
#############################################

def _mismatch_variants(seq):
    variants = set()
    for pos in range(len(seq)):
        for base in 'ACGTN':
            variants.add(seq[:pos] + base + seq[pos+1:])
    return variants

CONSTANT1_VARIANTS = _mismatch_variants("CAAGTTGATA")  # R1[18:28]
CONSTANT2_VARIANTS = _mismatch_variants("ATCTTGTGGA")  # R2[10:20]


def count_grna_reads(r1_file, r2_file, r3_file, ligation_barcode_file, RT_barcode_file,
                     inner_i7_barcode_file, sgrna_barcode_file, sgrna_annotation_file,
                     output_dir, sample_name, min_umi_threshold=10):

    ### 1. Load barcode dicts
    print("Loading barcodes...", flush=True)

    def _load_barcode_dict(tsv_file):
        df = pd.read_table(tsv_file)
        return {row['barcode_1bp_substitution']: row['original_barcode'] for _, row in df.iterrows()}

    ligation_barcodes = _load_barcode_dict(ligation_barcode_file)
    RT_barcodes = _load_barcode_dict(RT_barcode_file)
    inner_i7_barcodes = _load_barcode_dict(inner_i7_barcode_file)
    sgrna_barcodes = _load_barcode_dict(sgrna_barcode_file)

    annot_df = pd.read_table(sgrna_annotation_file)
    sgrna_annotation = dict(zip(annot_df['gRNA_seq'], annot_df['names']))

    ### 2. Initialize
    os.makedirs(output_dir, exist_ok=True)

    umi_dict = {}
    read_counter = {
        'filtered_constant_region': 0,
        'filtered_ligation_barcode': 0,
        'filtered_inner_i7': 0,
        'filtered_RT_barcode': 0,
        'filtered_sgrna': 0,
        'passed': 0,
    }
    start_time = time.time()

    ### 3. Process reads
    print("Processing FASTQ files...", flush=True)
    handles = [gzip.open(f, "rt") for f in [r1_file, r2_file, r3_file]]
    try:
        readers = [SeqIO.parse(h, "fastq") for h in handles]
        for read_count, (r1, r2, r3) in enumerate(zip(*readers), 1):

            if read_count % 500000 == 0:
                elapsed = str(timedelta(seconds=int(time.time() - start_time)))
                print(f"Processing read {read_count:,}... ({elapsed} elapsed)", flush=True)
                start_time = time.time()

            # Constant region filter
            if str(r1[18:28].seq) not in CONSTANT1_VARIANTS or str(r2[10:20].seq) not in CONSTANT2_VARIANTS:
                read_counter['filtered_constant_region'] += 1
                continue

            # Ligation barcode (R3[0:10])
            lig = ligation_barcodes.get(str(r3[0:10].seq))
            if not lig:
                read_counter['filtered_ligation_barcode'] += 1
                continue

            # Inner i7 barcode (R2[0:10])
            i7 = inner_i7_barcodes.get(str(r2[0:10].seq))
            if not i7:
                read_counter['filtered_inner_i7'] += 1
                continue

            # RT barcode (R1[8:18])
            rt = RT_barcodes.get(str(r1[8:18].seq))
            if not rt:
                read_counter['filtered_RT_barcode'] += 1
                continue

            # sgRNA sequence (R2[35:55])
            sgrna = sgrna_barcodes.get(str(r2[35:55].seq))
            if not sgrna:
                read_counter['filtered_sgrna'] += 1
                continue

            # Build cell ID and record UMI
            cell_id = f"{sample_name}{i7}.{lig}{rt}"
            umi = str(r1[0:8].seq)
            sgrna_name = sgrna_annotation.get(sgrna, sgrna)
            umi_dict.setdefault(cell_id, {}).setdefault(sgrna_name, set()).add(umi)
            read_counter['passed'] += 1

    finally:
        for h in handles:
            h.close()

    ### 4. Filter cells and write count TSV
    print("Writing output files...", flush=True)

    rows = []
    for cell_id, sgrna_dict in umi_dict.items():
        total_umis = sum(len(s) for s in sgrna_dict.values())
        if total_umis < min_umi_threshold:
            continue
        for sgrna_name, umi_set in sgrna_dict.items():
            rows.append((cell_id, sgrna_name, len(umi_set)))

    count_dataframe = pd.DataFrame(rows, columns=['barcode', 'sgrna_name', 'unique_umi_count'])
    count_file = os.path.join(output_dir, f"{sample_name}-sgrna_count.tsv")
    count_dataframe.to_csv(count_file, sep='\t', index=False)

    ### 5. Write summary TSV
    read_counter['total_reads_processed'] = sum(read_counter.values())
    summary_dataframe = pd.Series(read_counter).reset_index()
    summary_dataframe.columns = ['read_class', 'count']
    summary_dataframe['percent'] = summary_dataframe['count'] / read_counter['total_reads_processed'] * 100
    summary_file = os.path.join(output_dir, f"{sample_name}-sgrna_summary.tsv")
    summary_dataframe.to_csv(summary_file, sep='\t', index=False)

    print(f"Done. {read_counter['passed']:,} reads passed all filters. "
          f"Output: {count_file}", flush=True)

#############################################
########## 3. Main
#############################################

def main(args):
    count_grna_reads(
        r1_file=args.r1,
        r2_file=args.r2,
        r3_file=args.r3,
        ligation_barcode_file=args.ligation_barcode_file,
        RT_barcode_file=args.RT_barcode_file,
        inner_i7_barcode_file=args.inner_i7_barcode_file,
        sgrna_barcode_file=args.sgrna_barcode_file,
        sgrna_annotation_file=args.sgrna_annotation_file,
        output_dir=args.output_dir,
        sample_name=args.sample_name,
        min_umi_threshold=args.min_umi_threshold,
    )

#############################################
########## 4. Run
#############################################
if __name__ == "__main__":
    main(args)
