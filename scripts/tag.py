#############################################
########## 1. Import Libraries
#############################################
import pysam
import os
import time
from datetime import timedelta
import multiprocessing as mp

#############################################
########## 2. Processing function
#############################################

########## 1. Consecutive

def tag_bam(input_file, output_file):
    
    # Read BAM
    print("Reading BAM...", flush=True)
    input_bam = pysam.AlignmentFile(input_file, "rb")
    
    # Initialize output file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    output_bam = pysam.AlignmentFile(output_file, 'wb', template=input_bam)
    
    # Initialize time
    start_time = time.time()

    # Process each read in the BAM file sequentially
    print('Looping through reads...', flush=True)
    for i, read in enumerate(input_bam):
        
        # Update counter and report if it is a multiple of 100000
        if i % 1000000 == 0:
            formatted_time = str(timedelta(seconds=int(time.time() - start_time)))
            print(f"Processing line {i:,}... ({formatted_time} elapsed since previous cycle)", flush=True)
            start_time = time.time()
            
        # Extract barcode and UMI from the read name
        barcode, umi, _ = read.query_name.split(",")
        
        # Add barcode (CB) and UMI (UB) tags
        read.set_tag("CB", barcode)
        read.set_tag("UR", umi)
        
        # Write read
        output_bam.write(read)

    # Close BAM files
    input_bam.close()
    output_bam.close()

    # Update
    print("Done!", flush=True)
    
########## 2. Parallel
def process_reads(reads):
    processed_reads = []
    for read in reads:
        try:
            barcode, umi, _ = read.query_name.split(",")
            read.set_tag("CB", barcode)
            read.set_tag("UB", umi)
            processed_reads.append(read)
        except ValueError:
            # Handle reads without expected format
            pass
    return processed_reads

def process_chunk(chunk_data):
    
    # Unpack chunk information
    input_file, output_file, reference, chunk_index = chunk_data
    input_bam = pysam.AlignmentFile(input_file, "rb")
    output_bam = pysam.AlignmentFile(output_file, "wb", template=input_bam)
    
    # Process each read in the reference chunk
    for read in input_bam.fetch(reference=reference):
        try:
            # Extract barcode and UMI from the read name
            name_parts = read.query_name.split(",")
            
            # Add barcode (CB) and UMI (UB) tags
            read.set_tag("CB", name_parts[0])
            read.set_tag("UB", name_parts[1])
            
            # Write to output BAM
            output_bam.write(read)
        except IndexError:
            continue

    input_bam.close()
    output_bam.close()
    return chunk_index, output_file

def parallel_add_tags(input_file, output_file, num_processes):
    
    # Read BAM
    print("Reading BAM...", flush=True)
    input_bam = pysam.AlignmentFile(input_file, "rb")
    
    # Divide BAM by reference sequences (chromosomes/contigs)
    print("Dividing BAM by reference...", flush=True)
    references = input_bam.references
    chunks = []
    for i, reference in enumerate(references):
        temp_output = f"{output_file}.chunk{i}.bam"
        chunks.append((input_file, temp_output, reference, i))
    
    input_bam.close()
    
    # Parallel processing
    print("Adding tags...", flush=True)
    with mp.Pool(num_processes) as pool:
        results = pool.map(process_chunk, chunks)
    
    # Sort results by chunk index and merge
    print("Merging results...", flush=True)
    sorted_results = sorted(results, key=lambda x: x[0])
    with pysam.AlignmentFile(output_file, "wb", template=pysam.AlignmentFile(input_file, "rb")) as final_bam:
        for _, temp_output in sorted_results:
            with pysam.AlignmentFile(temp_output, "rb") as temp_bam:
                for read in temp_bam:
                    final_bam.write(read)
            os.remove(temp_output)
    
    print("Done!", flush=True)

#############################################
########## 3. Main
#############################################
    
# Define main function
def main(args):
    
    if args.num_processes > 1:
        
        parallel_add_tags(
            input_file = args.input_file,
            output_file = args.output_file,
            num_processes = args.num_processes
        )
        
    else:
        
        tag_bam(
            input_file = args.input_file,
            output_file = args.output_file
        )


#############################################
########## 4. Run
#############################################
# Run main function
if __name__ == "__main__":
    main(args)