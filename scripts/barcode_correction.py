#############################################
########## 1. Import Libraries
#############################################
import pandas as pd

#############################################
########## 2. Processing functions
#############################################

def _mismatch_variants(seq):
    variants = set()
    for pos in range(len(seq)):
        for base in 'ACGTN':
            variants.add(seq[:pos] + base + seq[pos+1:])
    return variants


def _build_correction_dict(whitelist):
    correction = {}
    for seq in whitelist:
        for variant in _mismatch_variants(seq):
            if variant in correction:
                correction[variant] = None  # ambiguous — two whitelist seqs within distance 2
            else:
                correction[variant] = seq
    return {k: v for k, v in correction.items() if v is not None}


def build_correction_dict_file(input_file, sequence_column, output_file):
    whitelist = pd.read_table(input_file)[sequence_column].drop_duplicates().tolist()
    correction = _build_correction_dict(whitelist)

    correction_dataframe = pd.DataFrame(
        sorted(correction.items()),
        columns=['barcode_1bp_substitution', 'original_barcode'],
    )
    correction_dataframe.to_csv(output_file, sep='\t', index=False)

    print(f"{len(whitelist):,} whitelist sequences -> {len(correction_dataframe):,} corrected variants "
          f"(some variants discarded as ambiguous between two whitelist sequences). "
          f"Output: {output_file}", flush=True)

#############################################
########## 3. Main
#############################################

def main(args):
    build_correction_dict_file(
        input_file=args.input_file,
        sequence_column=args.sequence_column,
        output_file=args.output_file,
    )
