#############################################
########## 1. Import Libraries
#############################################
import sys
import gzip
from functools import partial
import pickle
import re
from multiprocessing import Pool
from Bio.Seq import Seq
import argparse
import os
import time
from datetime import timedelta
import pandas as pd
import json
from Bio import SeqIO

#############################################
########## 2. Processing function
#############################################

def barcode_reads_paired_end(ligation_barcode_read_file, umi_cell_barcode_read_file, biological_read_file, output_basename, ligation_barcode_file, RT_barcode_file, randomN_barcode_file):
    
    ### 1. Read barcodes
    print("Loading barcodes...", flush=True)
    
    # Get the ligation barcode list
    ligation_barcode_dataframe = pd.read_table(ligation_barcode_file)
    ligation_barcodes = {rowData['barcode_1bp_substitution']: rowData['original_barcode'] for index, rowData in ligation_barcode_dataframe.iterrows()}
    
    # Get the RT barcode list
    RT_barcode_dataframe = pd.read_table(RT_barcode_file)
    RT_barcodes = {rowData['barcode_1bp_substitution']: rowData['original_barcode'] for index, rowData in RT_barcode_dataframe.iterrows()}

    # Get randomN barcodes
    with open(randomN_barcode_file, "rt") as barcodes:
        randomN_barcodes = [line.strip() for line in barcodes]
    
    ### 2. Initialize output files
    # Create outdir
    os.makedirs(os.path.dirname(output_basename), exist_ok=True)
    
    # Output files
    R1_output = gzip.open(output_basename + ".R1.fastq.gz", 'wt', encoding='utf-8')
    R2_output = gzip.open(output_basename + ".R2.fastq.gz", 'wt', encoding='utf-8')
        
    ### 3. Read files
    # Counter
    start_time = time.time()
    read_counter = {
        'filtered_no_ligation_barcode': 0,
        'filtered_no_RT_barcode': 0,
        'filtered_too_short': 0,
        'passed_all_filters': 0
    }
    
    # Input files as variables
    files = [ligation_barcode_read_file, umi_cell_barcode_read_file, biological_read_file]
    
    # Open all gzipped files as handles
    handles = [gzip.open(file, "rt", encoding="utf-8") for file in files]
    try:
        
        # Use SeqIO to parse each handle
        readers = [SeqIO.parse(handle, "fastq") for handle in handles]
        
        # Loop through reads
        print("Processing FASTQ files...", flush=True)
        for read_count, records in enumerate(zip(*readers), start=1):                

            # Counter and timer
            if read_count % 500000 == 0:
                formatted_time = str(timedelta(seconds=int(time.time() - start_time)))
                print(f"Processing line {read_count:,}... ({formatted_time} elapsed since previous cycle)", flush=True)
                start_time = time.time()
                
            # Split into separate variables
            ligation_barcode_read, umi_cell_barcode_read, biological_read = records

            # Extract expected ligation barcode from read
            ligation_barcode_sequence = ligation_barcode_read.seq[:10]
              
            # Extract corrected ligation barcode
            corrected_ligation_barcode = ligation_barcodes.get(ligation_barcode_sequence)
            
            # If the ligation barcode is not in the list of expected barcodes, skip the read
            if not corrected_ligation_barcode:
                # print(f"Read {read_count}: No ligation barcode found.", flush=True)
                read_counter['filtered_no_ligation_barcode'] += 1
                continue
            
            # If the ligation barcode is in the list of expected barcodes, continue processing
            else:
                # print(f"Read {read_count}: Ligation barcode found: {corrected_ligation_barcode}", flush=True)
                
                # Extract expected RT barcode from read
                RT_barcode_sequence = umi_cell_barcode_read.seq[8:18]
                
                # Extract corrected ligation barcode
                corrected_RT_barcode = RT_barcodes.get(RT_barcode_sequence)

                #  If the RT barcode is not in the list of expected barcodes, skip the read
                if not corrected_RT_barcode:
                    # print(f"Read {read_count}: No RT barcode found.", flush=True)
                    read_counter['filtered_no_RT_barcode'] += 1
                    continue
                
                # If the RT barcode is in the list of expected barcodes, continue processing
                else:
                    # print(f"Read {read_count}: RT barcode found: {corrected_RT_barcode}", flush=True)

                    # Extract UMI from read
                    UMI_sequence = umi_cell_barcode_read.seq[:8]
                    
                    ### Trim UMI and barcode
                    # Trim the first 18 bases from the UMI cell barcode read (UMI=8bp at the beginning, then 10bp cell barcode)
                    R1_out_record = umi_cell_barcode_read[18:]
                    
                    ### Remove poly-T region if RT barcode is not in randomN barcodes
                    # If the RT barcode is not in the list of randomN barcodes, trim the leading oligo-dT region
                    if corrected_RT_barcode not in randomN_barcodes:
                        # if the read contains more than 5 Ts at the beginning, find the first non-T base location index
                        match = re.match(f"^T{{5,}}", str(R1_out_record.seq))
                        if match:
                            R1_out_record = R1_out_record[match.end():]
                    
                    ### Remove anything after the Illumina adapter sequence, including the adapter sequence
                    # Find the Illumina adapter sequence in the read, allowing for 1 mismatch
                    illumina_adapter_sequence = "CTGTCTCTTATACACAT"
                    max_mismatches = 1
                    pattern = f"({illumina_adapter_sequence}){{s<={max_mismatches}}}"
                    R1_end_location = re.search(illumina_adapter_sequence, str(R1_out_record.seq))
                    
                    # If the adapter sequence is found in the read, trim the read beyond that location
                    if R1_end_location:
                        R1_out_record = R1_out_record[:R1_end_location.start()]
                        
                    ### Remove anything after the reverse complement RT barcode sequence in the second read
                    # Search for the reverse complement of the RT barcode in the biological read
                    R2_end_location = re.search(str(Seq(corrected_RT_barcode).reverse_complement()), str(biological_read.seq))

                    # If the RT barcode is found in the biological read, trim the read beyond that location. This indicates a short RNA fragment.
                    if R2_end_location:
                        R2_out_record = biological_read[:R2_end_location.start()]
                        
                    # If the RT barcode is not found in the biological read, keep the entire read
                    else:
                        R2_out_record = biological_read
                        
                    ### Add barcode and UMI to read names
                    # Add the corrected ligation barcode, corrected RT barcode, and UMI to the read ID
                    R1_out_record.id = corrected_ligation_barcode + corrected_RT_barcode + ',' + UMI_sequence + ',' + umi_cell_barcode_read.id
                    R2_out_record.id = corrected_ligation_barcode + corrected_RT_barcode + ',' + UMI_sequence + ',' + biological_read.id
                    
                    # Make sure trimmed sequences are still long enough, and write them to output
                    if len(R1_out_record.seq) > 20 and len(R2_out_record.seq) > 20:
                        read_counter['passed_all_filters'] += 1
                        
                        # Write to output
                        SeqIO.write(R1_out_record, R1_output, "fastq")
                        SeqIO.write(R2_out_record, R2_output, "fastq")
                        
                    else:
                        read_counter['filtered_too_short'] += 1
            
            # # Stop
            # max_reads = 10
            # if read_count >= max_reads:
            #     print(f"Stopped after {max_reads} reads.")
            #     break
    finally:
        
        # Ensure all files are closed
        for handle in handles:
            handle.close()
            
        # Close output files
        R1_output.close()
        R2_output.close()
        
    # Summary file
    read_counter['total_reads_processed'] = sum(read_counter.values())
    summary_file = output_basename + "-summary.tsv"
    summary_dataframe = pd.Series(read_counter).to_frame().reset_index().rename(columns={'index': 'read_class', 0: 'count'})
    summary_dataframe['percent'] = summary_dataframe['count'] / read_counter['total_reads_processed'] * 100
    summary_dataframe.to_csv(summary_file, index=False, sep='\t')

def barcode_reads_paired_end_v0(ligation_barcode_read_file, umi_cell_barcode_read_file, biological_read_file, output_basename, ligation_barcode_file, RT_barcode_file, randomN_barcode_file):
    
    ### 1. Read barcodes
    print("Loading barcodes...", flush=True)
    
    # generate the ligation barcode list
    ligation_barcode_dataframe = pd.read_table(ligation_barcode_file)
    ligation_barcodes = {rowData['barcode_1bp_substitution']: rowData['original_barcode'] for index, rowData in ligation_barcode_dataframe.iterrows()}
    
    # generate the RT barcode list:
    RT_barcode_dataframe = pd.read_table(RT_barcode_file)
    RT_barcodes = {rowData['barcode_1bp_substitution']: rowData['original_barcode'] for index, rowData in RT_barcode_dataframe.iterrows()}

    # load randomN barcodes
    with open(randomN_barcode_file, "rt") as barcodes:
        randomN_barcodes = [line.strip() for line in barcodes]
        
    ### 2. Read files (incorrect matching is intentional because somehow the naming changes between the files and the code - check original for confirmation)
    print("Reading FASTQ files...", flush=True)
    R1_input = gzip.open(ligation_barcode_read_file, 'rt', encoding='utf-8')
    R2_input = gzip.open(biological_read_file, 'rt', encoding='utf-8')
    I5_input = gzip.open(umi_cell_barcode_read_file, 'rt', encoding='utf-8')
    
    # Create outdir
    os.makedirs(os.path.dirname(output_basename), exist_ok=True)
    
    # Output files
    R1_output = gzip.open(output_basename + ".R1.fastq.gz", 'wt', encoding='utf-8')
    R2_output = gzip.open(output_basename + ".R2.fastq.gz", 'wt', encoding='utf-8')
    
    # Read line 1 of 4 (header) from input files (done before the loop to emulate a do-while loop; we want it to execute at least once. It is done again at the end of the loop)
    R1_in1 = R1_input.readline()
    R2_in1 = R2_input.readline()
    I5_in1 = I5_input.readline()
    
    # Set counts to 0
    read_counter = {
        'filtered_no_ligation_barcode': 0,
        'filtered_no_RT_barcode': 0,
        'filtered_too_short': 0,
        'passed_all_filters': 0
    }
    
    ### 3. Read through files
    print("Processing FASTQ files...", flush=True)
    
    # Start timer (for progress updates)
    start_time = time.time()
    
    # Read through files
    while (R1_in1):
        
        # # Counter and timer
        # if total_line % 1000000 == 0:
        #     formatted_time = str(timedelta(seconds=int(time.time() - start_time)))
        #     print(f"Processing line {total_line:,}... ({formatted_time} elapsed since previous cycle)", flush=True)
        #     start_time = time.time()
            
        # Read line 2 of 4 (sequence) from input files
        R1_in2 = R1_input.readline()
        R2_in2 = R2_input.readline() 
        I5_in2 = I5_input.readline()
        
        # Read line 3 of 4 (+ separator) from input files
        R1_in3 = R1_input.readline()
        R2_in3 = R2_input.readline() 
        I5_in3 = I5_input.readline() # not used
        
        # Read line 4 of 4 (quality score) from input files
        R1_in4 = R1_input.readline()
        R2_in4 = R2_input.readline() 
        I5_in4 = I5_input.readline() # not used
        
        # Extract expected ligation barcode from read
        ligation_barcode_sequence = I5_in2[:10]
        
        # Extract corrected ligation barcode
        corrected_ligation_barcode = ligation_barcodes.get(ligation_barcode_sequence)
        
        # If the ligation barcode is not in the list of expected barcodes, skip the read
        if not corrected_ligation_barcode:
            read_counter['filtered_no_ligation_barcode'] += 1
            continue
        
        # If the ligation barcode is in the list of expected barcodes, continue processing
        else:
            
            # Extract expected RT barcode from read
            RT_barcode_sequence = R1_in2[8:18]
            
            # Extract corrected ligation barcode
            corrected_RT_barcode = RT_barcodes.get(RT_barcode_sequence)

            #  If the RT barcode is not in the list of expected barcodes, skip the read
            if not corrected_RT_barcode:
                read_counter['filtered_no_RT_barcode'] += 1
                continue
            
            # If the RT barcode is in the list of expected barcodes, continue processing
            else:

                # Extract UMI from read
                UMI_sequence = R1_in2[:8]
               
                # Create header lines for output files including barcodes and UMIs
                R1_out1 = '@' + corrected_ligation_barcode + corrected_RT_barcode + ',' + UMI_sequence + ',' + R1_in1[1:]
                R2_out1 = '@' + corrected_ligation_barcode + corrected_RT_barcode + ',' + UMI_sequence + ',' + R2_in1[1:]
                
                # Create sequence and quality lines for R1
                if corrected_RT_barcode in randomN_barcodes:
                    R1_out2 = R1_in2[18:]
                    R1_out4 = R1_in4[18:]
                else:
                    R1_out2 = R1_in2[33:] # longer to trim oligo-dT region
                    R1_out4 = R1_in4[33:]
                    
                # Trim ends of sequence and quality lines for R1
                R1_end_seq = "CTGTCTCTTATACACAT" # Illumina adapter sequence
                R1_end_loc = re.search(R1_end_seq, R1_out2)
                if R1_end_loc is not None:
                    R1_out2 = R1_out2[:R1_end_loc.start()] + "\n"
                    R1_out4 = R1_out4[:R1_end_loc.start()] + "\n"
                    
                # Create and sequence and quality lines for R2              
                R2_out2 = R2_in2
                R2_out4 = R2_in4
                    
                # Trim ends of sequence and quality lines for R2                 
                RT_bar_seq = Seq(RT_barcode_sequence)
                R2_end_seq = str(RT_bar_seq.reverse_complement())
                
                R2_end_loc = re.search(R2_end_seq, R2_out2)
                if R2_end_loc is not None:
                    R2_out2 = R2_out2[:R2_end_loc.start()] + "\n"
                    R2_out4 = R2_out4[:R2_end_loc.start()] + "\n" 
                
                # Pass third separator line to output
                R1_out3 = R1_in3
                R2_out3 = R2_in3
                
                # Make sure trimmed sequences are still long enough, and write them to output
                if len(R1_out2) > 20 and len(R2_out2) > 20:
                    read_counter['passed_all_filters'] += 1
                    
                    R1_output.write(R1_out1)
                    R1_output.write(R1_out2)
                    R1_output.write(R1_out3)
                    R1_output.write(R1_out4)
                    
                    R2_output.write(R2_out1)
                    R2_output.write(R2_out2)
                    R2_output.write(R2_out3)
                    R2_output.write(R2_out4)
                    
                else:
                    read_counter['filtered_too_short'] += 1
                    
                
        # Read line 1 of 4 (header) from input files for next cycle
        R1_in1 = R1_input.readline()
        R2_in1 = R2_input.readline()
        I5_in1 = I5_input.readline()
        
    # Update counter - adding all values together
    read_counter['total_reads_processed'] = sum(read_counter.values())
    
    # Close files
    R1_input.close()
    R2_input.close()
    I5_input.close()
    
    R1_output.close()
    R2_output.close()
    
    # Summary JSON
    summary_file = output_basename + "-summary.json"
    with open(summary_file, 'w') as f:
        f.write(json.dumps(read_counter, indent=4))

#############################################
########## 3. Main
#############################################

# Define main function
def main(args):
    
    # Run
    barcode_reads_paired_end(
        ligation_barcode_read_file = args.ligation_barcode_read_file,
        umi_cell_barcode_read_file = args.umi_cell_barcode_read_file,
        biological_read_file = args.biological_read_file,
        output_basename = args.output_basename,
        ligation_barcode_file = args.ligation_barcode_file,
        RT_barcode_file = args.RT_barcode_file,
        randomN_barcode_file = args.randomN_barcode_file
    )

#############################################
########## 4. Run
#############################################
# Run main function
if __name__ == "__main__":
    main(args)


# from Bio import SeqIO

# def split_fastq_stream(input_file, n_subsets, output_prefix):
#     output_files = [open(f"{output_prefix}_subset_{i+1}.fastq", "w") for i in range(n_subsets)]
#     writers = [SeqIO.SeqIO.write for _ in output_files]  # Create writers

#     for i, record in enumerate(SeqIO.parse(input_file, "fastq")):
#         subset_index = i % n_subsets  # Round-robin assignment
#         SeqIO.write(record, output_files[subset_index], "fastq")

#     for f in output_files:
#         f.close()

# # Example usage
# split_fastq_stream("input.fastq", 3, "output")


# summary_data = {
#         "Read 1": ligation_barcode_read_file,
#         "Read 2": umi_cell_barcode_read_file,
#         "Read 3": biological_read_file,
#         "Ligation Barcode File": ligation_barcode_file,
#         "RT Barcode File": RT_barcode_file,
#         "RandomN Barcode File": randomN_barcode_file,
#         "R1 Output": output_basename + ".R1.fastq.gz",
#         "R2 Output": output_basename + ".R2.fastq.gz",
#         "total_reads": total_line,
#         "Filtered Lines": filtered_line,
#         "Filtering Percentage": filtered_line / total_line * 100
#     }


# def barcode_reads_paired_end_v0(ligation_barcode_read_file, umi_cell_barcode_read_file, biological_read_file, output_basename, ligation_barcode_file, RT_barcode_file, randomN_barcode_file):
#     ''' original code
#     # open the read files
#     # Note: the "R2" fastq file is actually I5, and the "R3" fastq file is actually R2
#     R1_input = gzip.open(input_folder + "/" + sample + ".R1.fastq.gz", 'rt', encoding='utf-8')
#     R2_input = gzip.open(input_folder + "/" + sample + ".R3.fastq.gz", 'rt', encoding='utf-8')
#     I5_input = gzip.open(input_folder + "/" + sample + ".R2.fastq.gz", 'rt', encoding='utf-8')
    
#     # open the output files
#     R1_output = gzip.open(output_folder + "/" + sample + ".R1.fastq.gz", 'wt', encoding='utf-8')
#     R2_output = gzip.open(output_folder + "/" + sample + ".R2.fastq.gz", 'wt',encoding='utf-8')
#     '''
    
#     ### 1. Read barcodes
#     print("Loading barcodes...", flush=True)
    
#     # generate the ligation barcode list
#     ligation_barcode_dataframe = pd.read_table(ligation_barcode_file)
#     ligation_barcodes = {rowData['barcode_1bp_substitution']: rowData['original_barcode'] for index, rowData in ligation_barcode_dataframe.iterrows()}
    
#     # generate the RT barcode list:
#     RT_barcode_dataframe = pd.read_table(RT_barcode_file)
#     RT_barcodes = {rowData['barcode_1bp_substitution']: rowData['original_barcode'] for index, rowData in RT_barcode_dataframe.iterrows()}

#     # load randomN barcodes
#     with open(randomN_barcode_file, "rt") as barcodes:
#         randomN_barcodes = [line.strip() for line in barcodes]
        
#     ### 2. Read files (incorrect matching is intentional because somehow the naming changes between the files and the code - check original for confirmation)
#     print("Reading FASTQ files...", flush=True)
#     R1_input = gzip.open(ligation_barcode_read_file, 'rt', encoding='utf-8')
#     R2_input = gzip.open(biological_read_file, 'rt', encoding='utf-8')
#     I5_input = gzip.open(umi_cell_barcode_read_file, 'rt', encoding='utf-8')
    
#     # Create outdir
#     os.makedirs(os.path.dirname(output_basename), exist_ok=True)
    
#     # Output files
#     R1_output = gzip.open(output_basename + ".R1.fastq.gz", 'wt', encoding='utf-8')
#     R2_output = gzip.open(output_basename + ".R2.fastq.gz", 'wt', encoding='utf-8')
    
#     # Read line 1 of 4 (header) from input files (done before the loop to emulate a do-while loop; we want it to execute at least once. It is done again at the end of the loop)
#     R1_in1 = R1_input.readline()
#     R2_in1 = R2_input.readline()
#     I5_in1 = I5_input.readline()
    
#     # Set counts to 0
#     total_line = 0
#     filtered_line = 0
    
#     ### 3. Read through files
#     print("Processing FASTQ files...", flush=True)
    
#     # Start timer (for progress updates)
#     start_time = time.time()
    
#     # Read through files
#     while (R1_in1):
#         total_line += 1
        
#         # Counter and timer
#         if total_line % 1000000 == 0:
#             formatted_time = str(timedelta(seconds=int(time.time() - start_time)))
#             print(f"Processing line {total_line:,}... ({formatted_time} elapsed since previous cycle)", flush=True)
#             start_time = time.time()
            
#         # Read line 2 of 4 (sequence) from input files
#         R1_in2 = R1_input.readline()
#         R2_in2 = R2_input.readline() 
#         I5_in2 = I5_input.readline()
        
#         # Read line 3 of 4 (+ separator) from input files
#         R1_in3 = R1_input.readline()
#         R2_in3 = R2_input.readline() 
#         I5_in3 = I5_input.readline() # not used
        
#         # Read line 4 of 4 (quality score) from input files
#         R1_in4 = R1_input.readline()
#         R2_in4 = R2_input.readline() 
#         I5_in4 = I5_input.readline() # not used
        
#         # first check if the ligation barcode matches an expected barcode and correct sequence to closest barcode
#         ligation_barcode_sequence = I5_in2[:10]
#         if ligation_barcode_sequence in ligation_barcodes:
#             corrected_ligation_barcode = ligation_barcodes[ligation_barcode_sequence]
            
#             # check if the RT barcode matches an expected barcode and correct sequence to closest barcode
#             RT_barcode_sequence = R1_in2[8:18]
#             if RT_barcode_sequence in RT_barcodes:

#                 RT_barcode = RT_barcodes[RT_barcode_sequence]
#                 UMI = R1_in2[:8]
               
#                 # Create header lines for output files including barcodes and UMIs
#                 R1_out1 = '@' + corrected_ligation_barcode + RT_barcode + ',' + UMI + ',' + R1_in1[1:]
#                 R2_out1 = '@' + corrected_ligation_barcode + RT_barcode + ',' + UMI + ','+ R2_in1[1:]
                
#                 # Create sequence and quality lines for R1
#                 if RT_barcode in randomN_barcodes:
#                     R1_out2 = R1_in2[18:]
#                     R1_out4 = R1_in4[18:]
#                 else:
#                     R1_out2 = R1_in2[33:] # longer to trim oligo-dT region
#                     R1_out4 = R1_in4[33:]
                    
#                 # Trim ends of sequence and quality lines for R1
#                 R1_end_seq = "CTGTCTCTTATACACAT"
                
#                 R1_end_loc = re.search(R1_end_seq, R1_out2)
#                 if R1_end_loc is not None:
#                     R1_out2 = R1_out2[:R1_end_loc.start()] + "\n"
#                     R1_out4 = R1_out4[:R1_end_loc.start()] + "\n"
                    
#                 # Create and sequence and quality lines for R2              
#                 R2_out2 = R2_in2
#                 R2_out4 = R2_in4
                    
#                 # Trim ends of sequence and quality lines for R2                 
#                 RT_bar_seq = Seq(RT_barcode_sequence)
#                 R2_end_seq = str(RT_bar_seq.reverse_complement())
                
#                 R2_end_loc = re.search(R2_end_seq, R2_out2)
#                 if R2_end_loc is not None:
#                     R2_out2 = R2_out2[:R2_end_loc.start()] + "\n"
#                     R2_out4 = R2_out4[:R2_end_loc.start()] + "\n" 
                
#                 # Pass third separator line to output
#                 R1_out3 = R1_in3
#                 R2_out3 = R2_in3
                
#                 # Make sure trimmed sequences are still long enough, and write them to output
#                 if len(R1_out2) > 20 and len(R2_out2) > 20:
#                     filtered_line += 1
                    
#                     R1_output.write(R1_out1)
#                     R1_output.write(R1_out2)
#                     R1_output.write(R1_out3)
#                     R1_output.write(R1_out4)
                    
#                     R2_output.write(R2_out1)
#                     R2_output.write(R2_out2)
#                     R2_output.write(R2_out3)
#                     R2_output.write(R2_out4)
                
#         # Read line 1 of 4 (header) from input files for next cycle
#         R1_in1 = R1_input.readline()
#         R2_in1 = R2_input.readline()
#         I5_in1 = I5_input.readline()
    
#     # Close files
#     R1_input.close()
#     R2_input.close()
#     I5_input.close()
    
#     R1_output.close()
#     R2_output.close()
    
#     # Summary JSON
#     summary_data = {
#         "Read 1": ligation_barcode_read_file,
#         "Read 2": umi_cell_barcode_read_file,
#         "Read 3": biological_read_file,
#         "Ligation Barcode File": ligation_barcode_file,
#         "RT Barcode File": RT_barcode_file,
#         "RandomN Barcode File": randomN_barcode_file,
#         "R1 Output": output_basename + ".R1.fastq.gz",
#         "R2 Output": output_basename + ".R2.fastq.gz",
#         "Total Lines": total_line,
#         "Filtered Lines": filtered_line,
#         "Filtering Percentage": filtered_line / total_line * 100
#     }
#     summary_file = output_basename + "-summary.json"
#     with open(summary_file, 'w') as f:
#         f.write(json.dumps(summary_data, indent=4))

# # Function to split a FASTQ file into subsets
# def split_fastq_stream(input_file, n_subsets, split_directory):
#     # Ensure the split directory exists
#     os.makedirs(split_directory, exist_ok=True)

#     # Create output file paths
#     output_files = [
#         os.path.join(split_directory, os.path.basename(input_file)) + f"-subset_{i+1}.fastq.gz"
#         for i in range(n_subsets)
#     ]
#     output_openfiles = [gzip.open(f, 'wt', encoding='utf-8') for f in output_files]

#     try:
#         # Stream through the input FASTQ file and split it
#         with gzip.open(input_file, "rt", encoding='utf-8') as handle:
#             for i, record in enumerate(SeqIO.parse(handle, "fastq")):
#                 subset_index = i % n_subsets  # Round-robin assignment
#                 SeqIO.write(record, output_openfiles[subset_index], "fastq")
#     finally:
#         # Close all output files
#         for f in output_openfiles:
#             f.close()
            
#     return output_files

# # Function for barcode reads with paired-end and parallel processing
# def barcode_reads_paired_end_parallel(ligation_barcode_read_file, umi_cell_barcode_read_file, biological_read_file, output_basename, ligation_barcode_file, RT_barcode_file, randomN_barcode_file, n_subsets=5):
    
#     # Create directory for split files
#     split_directory = os.path.join(os.path.dirname(output_basename), "split")
#     os.makedirs(split_directory, exist_ok=True)
    
#     # Split R1, R2, and R3 into subsets
#     print("Splitting R1...", flush=True)
#     R1_files = split_fastq_stream(ligation_barcode_read_file, n_subsets, split_directory)
#     print("Splitting R2...", flush=True)
#     R2_files = split_fastq_stream(umi_cell_barcode_read_file, n_subsets, split_directory)
#     print("Splitting R3...", flush=True)
#     R3_files = split_fastq_stream(biological_read_file, n_subsets, split_directory)
    
#     # Output split files (example)
#     print(f"Split files saved in: {split_directory}", flush=True)
#     print(f"R1 files: {R1_files}", flush=True)
#     print(f"R2 files: {R2_files}", flush=True)
#     print(f"R3 files: {R3_files}", flush=True)
    