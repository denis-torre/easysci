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


def _load_barcode_dict(tsv_file):
    df = pd.read_table(tsv_file)
    return {row['barcode_1bp_substitution']: row['original_barcode'] for _, row in df.iterrows()}


def _load_inner_i7_barcode_dict(tsv_file):
    """Load inner i7 barcodes from sample_id/barcode format TSV, building 1bp correction dict on the fly."""
    df = pd.read_table(tsv_file)
    return _build_correction_dict(df['barcode'].tolist())


def _build_correction_dict(whitelist):
    correction = {}
    for seq in whitelist:
        for variant in _mismatch_variants(seq):
            if variant in correction:
                correction[variant] = None  # ambiguous — two whitelist seqs within distance 2
            else:
                correction[variant] = seq
    return {k: v for k, v in correction.items() if v is not None}


def _load_barcodes(ligation_barcode_file, RT_barcode_file, inner_i7_barcode_file):
    ligation_barcodes = _load_barcode_dict(ligation_barcode_file)
    RT_barcodes = _load_barcode_dict(RT_barcode_file)
    inner_i7_barcodes = _load_inner_i7_barcode_dict(inner_i7_barcode_file)
    return ligation_barcodes, RT_barcodes, inner_i7_barcodes


def _count_umi_loop(r1_file, r2_file, r3_file, sgrna_lookup, ligation_barcodes,
                    RT_barcodes, inner_i7_barcodes, sgrna_annotation, sample_name,
                    n_reads, progress_label):
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

    handles = [gzip.open(f, "rt") for f in [r1_file, r2_file, r3_file]]
    try:
        readers = [SeqIO.parse(h, "fastq") for h in handles]
        for read_count, (r1, r2, r3) in enumerate(zip(*readers), 1):

            if n_reads and read_count > n_reads:
                break

            if read_count % 500000 == 0:
                elapsed = str(timedelta(seconds=int(time.time() - start_time)))
                print(f"{progress_label} {read_count:,}... ({elapsed} elapsed)", flush=True)
                start_time = time.time()

            if str(r1[18:28].seq) not in CONSTANT1_VARIANTS or str(r2[10:20].seq) not in CONSTANT2_VARIANTS:
                read_counter['filtered_constant_region'] += 1
                continue

            lig = ligation_barcodes.get(str(r3[0:10].seq))
            if not lig:
                read_counter['filtered_ligation_barcode'] += 1
                continue

            i7 = inner_i7_barcodes.get(str(r2[0:10].seq))
            if not i7:
                read_counter['filtered_inner_i7'] += 1
                continue

            rt = RT_barcodes.get(str(r1[8:18].seq))
            if not rt:
                read_counter['filtered_RT_barcode'] += 1
                continue

            sgrna = sgrna_lookup.get(str(r2[35:55].seq))
            if not sgrna:
                read_counter['filtered_sgrna'] += 1
                continue

            cell_id = f"{sample_name}{i7}.{lig}{rt}"
            umi = str(r1[0:8].seq)
            sgrna_name = sgrna_annotation.get(sgrna, sgrna)
            umi_dict.setdefault(cell_id, {}).setdefault(sgrna_name, set()).add(umi)
            read_counter['passed'] += 1

    finally:
        for h in handles:
            h.close()

    return umi_dict, read_counter


def _write_outputs(umi_dict, read_counter, min_umi_threshold, output_dir, sample_name):
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

    read_counter['total_reads_processed'] = sum(read_counter.values())
    summary_dataframe = pd.Series(read_counter).reset_index()
    summary_dataframe.columns = ['read_class', 'count']
    summary_dataframe['percent'] = summary_dataframe['count'] / read_counter['total_reads_processed'] * 100
    summary_file = os.path.join(output_dir, f"{sample_name}-sgrna_summary.tsv")
    summary_dataframe.to_csv(summary_file, sep='\t', index=False)

    print(f"Done. {read_counter['passed']:,} reads passed all filters. "
          f"Output: {count_file}", flush=True)


def count_grna_reads(r1_file, r2_file, r3_file, ligation_barcode_file, RT_barcode_file,
                     inner_i7_barcode_file, sgrna_barcode_file, sgrna_annotation_file,
                     output_dir, sample_name, min_umi_threshold=10, n_reads=None):

    print("Loading barcodes...", flush=True)
    ligation_barcodes, RT_barcodes, inner_i7_barcodes = _load_barcodes(
        ligation_barcode_file, RT_barcode_file, inner_i7_barcode_file)
    sgrna_barcodes = _load_barcode_dict(sgrna_barcode_file)
    annot_df = pd.read_table(sgrna_annotation_file)
    sgrna_annotation = dict(zip(annot_df['gRNA_seq'], annot_df['names']))

    os.makedirs(output_dir, exist_ok=True)

    print("Processing FASTQ files...", flush=True)
    umi_dict, read_counter = _count_umi_loop(
        r1_file, r2_file, r3_file,
        sgrna_lookup=sgrna_barcodes,
        ligation_barcodes=ligation_barcodes,
        RT_barcodes=RT_barcodes,
        inner_i7_barcodes=inner_i7_barcodes,
        sgrna_annotation=sgrna_annotation,
        sample_name=sample_name,
        n_reads=n_reads,
        progress_label="Processing read",
    )

    _write_outputs(umi_dict, read_counter, min_umi_threshold, output_dir, sample_name)


def count_grna_reads_autodiscovery(r1_file, r2_file, r3_file, ligation_barcode_file,
                                    RT_barcode_file, inner_i7_barcode_file,
                                    output_dir, sample_name,
                                    sgrna_annotation_file=None,
                                    min_umi_threshold=10, min_sgrna_count=10,
                                    n_reads=None):

    print("Loading barcodes...", flush=True)
    ligation_barcodes, RT_barcodes, inner_i7_barcodes = _load_barcodes(
        ligation_barcode_file, RT_barcode_file, inner_i7_barcode_file)

    sgrna_annotation = {}
    if sgrna_annotation_file:
        annot_df = pd.read_table(sgrna_annotation_file)
        sgrna_annotation = dict(zip(annot_df['gRNA_seq'], annot_df['names']))

    os.makedirs(output_dir, exist_ok=True)

    ### Pass 1 — discover sgRNA sequences
    print("Pass 1: discovering sgRNA sequences...", flush=True)
    sgrna_seq_counts = {}
    start_time = time.time()

    handles = [gzip.open(f, "rt") for f in [r1_file, r2_file, r3_file]]
    try:
        readers = [SeqIO.parse(h, "fastq") for h in handles]
        for read_count, (r1, r2, r3) in enumerate(zip(*readers), 1):

            if n_reads and read_count > n_reads:
                break

            if read_count % 500000 == 0:
                elapsed = str(timedelta(seconds=int(time.time() - start_time)))
                print(f"  Pass 1: read {read_count:,}... ({elapsed} elapsed)", flush=True)
                start_time = time.time()

            if str(r1[18:28].seq) not in CONSTANT1_VARIANTS or str(r2[10:20].seq) not in CONSTANT2_VARIANTS:
                continue
            if not ligation_barcodes.get(str(r3[0:10].seq)):
                continue
            if not inner_i7_barcodes.get(str(r2[0:10].seq)):
                continue
            if not RT_barcodes.get(str(r1[8:18].seq)):
                continue

            sgrna_seq = str(r2[35:55].seq)
            sgrna_seq_counts[sgrna_seq] = sgrna_seq_counts.get(sgrna_seq, 0) + 1
    finally:
        for h in handles:
            h.close()

    whitelist = {seq for seq, n in sgrna_seq_counts.items() if n >= min_sgrna_count}
    print(f"  {len(sgrna_seq_counts):,} unique sequences; "
          f"{len(whitelist):,} whitelisted (>= {min_sgrna_count} reads).", flush=True)

    whitelist_rows = [
        {'sgrna_seq': seq, 'read_count': sgrna_seq_counts[seq],
         'sgrna_name': sgrna_annotation.get(seq, seq)}
        for seq in sorted(whitelist, key=lambda s: -sgrna_seq_counts[s])
    ]
    whitelist_file = os.path.join(output_dir, f"{sample_name}-sgrna_whitelist.tsv")
    pd.DataFrame(whitelist_rows).to_csv(whitelist_file, sep='\t', index=False)

    sgrna_correction = _build_correction_dict(whitelist)

    ### Pass 2 — count UMIs using discovered whitelist
    print("Pass 2: counting UMIs...", flush=True)
    umi_dict, read_counter = _count_umi_loop(
        r1_file, r2_file, r3_file,
        sgrna_lookup=sgrna_correction,
        ligation_barcodes=ligation_barcodes,
        RT_barcodes=RT_barcodes,
        inner_i7_barcodes=inner_i7_barcodes,
        sgrna_annotation=sgrna_annotation,
        sample_name=sample_name,
        n_reads=n_reads,
        progress_label="  Pass 2: read",
    )

    _write_outputs(umi_dict, read_counter, min_umi_threshold, output_dir, sample_name)

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
        n_reads=getattr(args, 'n_reads', None),
    )

def main_autodiscovery(args):
    count_grna_reads_autodiscovery(
        r1_file=args.r1,
        r2_file=args.r2,
        r3_file=args.r3,
        ligation_barcode_file=args.ligation_barcode_file,
        RT_barcode_file=args.RT_barcode_file,
        inner_i7_barcode_file=args.inner_i7_barcode_file,
        sgrna_annotation_file=getattr(args, 'sgrna_annotation_file', None),
        output_dir=args.output_dir,
        sample_name=args.sample_name,
        min_umi_threshold=args.min_umi_threshold,
        min_sgrna_count=args.min_sgrna_count,
        n_reads=getattr(args, 'n_reads', None),
    )

#############################################
########## 4. Run
#############################################
if __name__ == "__main__":
    main(args)
