#############################################
########## 1. Import Libraries
#############################################
import collections
import numpy as np
import pandas as pd
import multiprocessing
import HTSeq
from functools import partial
import sys
import os
import time
from datetime import timedelta
import pickle

#############################################
########## 2. Supporting functions
#############################################

### Function to get fragment strandedness
# Only support reverse and unstranded for now, need to upgrade for forward but will likely not be necessary
def get_fragment_strand(read1, read2, library_strandedness):
    
    # Get conformation of fragment pair
    fragment_conformation = f'R1:{read1.iv.strand},R2:{read2.iv.strand}'
    
    # Assign fragment strandedness based on experiment strandedness and fragment conformation
    if library_strandedness == 'reverse':
        if fragment_conformation == 'R1:-,R2:+':
            fragment_strandedness = '+'
        elif fragment_conformation == 'R1:+,R2:-':
            fragment_strandedness = '-'
        else:
            fragment_strandedness = '.'
    elif library_strandedness == 'forward':
        raise NotImplementedError("Forward strandedness not implemented yet.")
    elif library_strandedness == 'none':
        fragment_strandedness = '.'
    else:
        raise ValueError(f"Experiment strandedness {strandedness} not recognized.")
    
    return fragment_strandedness

### Function to resolve multi-gene mapping reads (called in assign_read_pair_to_gene)
# Given a read pair, a set of genes, and further parameters, tries to assign the best gene
# find_closest_TES(read1, read2, fragment_strandedness, gtf_coordinate_index, candidate_genes=read_pair_assignment.split(','))
def find_closest_TES(read1, read2, fragment_strandedness, gtf_coordinate_index, candidate_genes):

    # Print
    # print(f"Resolving multi-gene reads for {read1.read.name}: {candidate_genes} READ1 chr{read1.iv.chrom}:{read1.iv.start}-{read1.iv.end} READ2 chr{read2.iv.chrom}:{read2.iv.start}-{read2.iv.end}")
    # Initialize list of TES distances
    tes_distances = []

    # Get fragment end closest to TES
    if fragment_strandedness=='+':
        fragment_end_closest_to_tes=max(read1.iv.end, read2.iv.end)
    elif fragment_strandedness=='-':
        fragment_end_closest_to_tes=min(read1.iv.start, read2.iv.start)
        
    ### Get TES distances using exon distance only - ignore introns! unless the read is intronic in which case, I will need to add intron distance, I suppose

    # Loop through candidate genes
    for candidate_gene in candidate_genes:
        
        #  Loop through TESs of candidate gene
        for tes in gtf_coordinate_index['TES'][candidate_gene]:
            
            # Get distance between fragment end and TES
            tes_coordinate = tes[1]
            distance_to_tes = tes_coordinate-fragment_end_closest_to_tes if fragment_strandedness == '+' else fragment_end_closest_to_tes-tes_coordinate
            tes_distances.append([candidate_gene, tes[3], distance_to_tes])

    # Get gene ID with closest read - using a low negative cutoff because in some cases the read end is just beyond the TES annotation in Ensembl
    tes_distance_dataframe = pd.DataFrame(tes_distances, columns=['gene_id', 'transcript_id', 'distance_to_tes']).sort_values('distance_to_tes')
    
    # Filter
    tes_distance_dataframe_filtered = tes_distance_dataframe.query('distance_to_tes>-50')
    
    # Check length of dataframe
    if len(tes_distance_dataframe_filtered) == 0:
        
        # If no TES found, assign ambiguous
        print(f"No TES found for {candidate_genes} READ1 chr{read1.iv.chrom}:{read1.iv.start}-{read1.iv.end} READ2 chr{read2.iv.chrom}:{read2.iv.start}-{read2.iv.end}. All TESs:")
        print(tes_distance_dataframe)
        assigned_gene_id = '_ambiguous'
    else:
        
        # If at least one TES found, assign the closest TES
        assigned_gene_id = tes_distance_dataframe_filtered.iloc[0]['gene_id']
        # print(f"Resolved multi-gene read for {read1.read.name}: {candidate_genes} READ1 chr{read1.iv.chrom}:{read1.iv.start}-{read1.iv.end} READ2 chr{read2.iv.chrom}:{read2.iv.start}-{read2.iv.end} to {assigned_gene_id} with distance {tes_distance_dataframe.iloc[0]['distance_to_tes']}")            

    # Return
    return assigned_gene_id

### Get overlaping genes
def get_overlapping_genes(read1, read2, fragment_strandedness, gtf_features):
    
    ### Part 1. Get overlapping genes
    # Create read dictionary - for looping through reads
    reads = {'read1': read1, 'read2': read2}

    # Initialize gene IDs - exon/intron to get feature-specific counting
    gene_ids = {x: set() for x in reads.keys()}
    
    # Loop through read pair
    for read_name, read in reads.items():

        # Loop through CIGAR portions of read
        for cigop in read.cigar:
            
            # Set strandedness of fragment - needed to adjust for reverse strandedness of read pairs
            cigop.ref_iv.strand = fragment_strandedness
            
            # Only include M type CIGAR alignments (match or mismatch - aka aligned to the genome, as opposed to N-introns, I-insertions, D-deletions and S-soft clipping)
            if cigop.type == 'M':
                
                try:

                    # Loop through read coordinates
                    for iv, val in gtf_features[cigop.ref_iv].steps():

                        # Find overlaps for the read and features, and add to the set (union operation)
                        gene_ids[read_name] |= val
                        
                except KeyError:
                    pass
                
    # Get union and intersection
    overlapping_genes = {
        'union': list(set.union(*gene_ids.values())),
        'intersection': list(set.intersection(*gene_ids.values()))
    }
    
    ### Part 2. Assign gene ID
    # Multigene binary
    multigene_read = False

    # First, check if there is a single gene in the union set, and if so, assign that
    if len(overlapping_genes['union']) == 1:
        read_pair_assignment = overlapping_genes['union'][0]
        
    # If there is more than one gene in the union set, try different approaches based on the intersection set
    elif len(overlapping_genes['union']) > 1:
        
        # If there is a single gene in the intersection set, assign that
        if len(overlapping_genes['intersection']) == 1:
            read_pair_assignment = overlapping_genes['intersection'][0]
            
        # If there is more than one gene in the intersection set, assign all
        elif len(overlapping_genes['intersection']) > 1:
                read_pair_assignment = ','.join(overlapping_genes['intersection'])
                multigene_read = True
        
        # If there is no gene in the exon intersection set, assign ambiguous (this means that each read in the pair is mapping to different genes)
        elif len(overlapping_genes['intersection']) == 0:
            read_pair_assignment = '_ambiguous'

    # If there is no gene in the union set, assign no feature overlap
    elif len(overlapping_genes['union']) == 0:
        read_pair_assignment = '_no_feature_overlap'
        
    # Return
    return read_pair_assignment, multigene_read
               
                
### Function to assign read pair to gene ID
# Given a read pair, a GTF index and further parameters, assigns a read pair to a gene
def assign_read_pair_to_gene(read1, read2, gtf_coordinate_index, count_introns, primer_type, library_strandedness, multigene_reads):

    ### 0. Get initial info    
    # Get fragment strandedness
    fragment_strandedness = get_fragment_strand(read1, read2, library_strandedness)
    
    ### 1. First, try exon match
    # Get overlapping genes
    search_space = 'exon'
    read_pair_assignment, multigene_read = get_overlapping_genes(read1, read2, fragment_strandedness, gtf_features=gtf_coordinate_index[search_space])
    
    # If there is more than one matching gene in the intersection set, try to resolve or label as multimapping
    if multigene_read:
        if multigene_reads == 'discard':
            read_pair_assignment = '_multiple'
        elif multigene_reads == 'closest_TES':
            read_pair_assignment = find_closest_TES(read1, read2, fragment_strandedness, gtf_coordinate_index, candidate_genes=read_pair_assignment.split(','))
            multigene_read = False

    ### 2. Second, expand to intron match
    # Expand to intron search, if needed
    if read_pair_assignment == '_no_feature_overlap' and count_introns:
        
        # Search full transcript overlap
        search_space = 'transcript'
        read_pair_assignment, multigene_read = get_overlapping_genes(read1, read2, fragment_strandedness, gtf_features=gtf_coordinate_index[search_space])
           
        # If there is more than one matching gene in the intersection set, try to resolve or label as multimapping
        if multigene_read:
            if multigene_reads == 'discard':
                read_pair_assignment = '_multiple'
            elif multigene_reads == 'closest_TES':
                read_pair_assignment = find_closest_TES(read1, read2, fragment_strandedness, gtf_coordinate_index, candidate_genes=read_pair_assignment.split(','))
                multigene_read = False
                
    ### 3. Set mapped feature, for annotation purpose (not used in counting)
    mapped_feature = search_space if read_pair_assignment != '_no_feature_overlap' else 'none'
     
    # Return
    return read_pair_assignment, multigene_read, mapped_feature


def resolve_multigene_reads(read_assignment_dataframe, multigene_reads, collapse_by_gene_name):
    raise NotImplementedError(f"Multigene reads resolution method {multigene_reads} not implemented yet.")
        
#############################################
########## 3. Processing function
#############################################

def easysci_count(input_bam, index_file, output_dir, library_strandedness, count_introns, primer_type, multigene_reads, read_subset, collapse_by_gene_name, sample_name):
    
    ### Read BAM
    # Open reader
    print("Reading BAM file...")
    bam_reader = HTSeq.BAM_Reader(input_bam)

    ### Read index
    # Read index
    print('Reading index file...')
    with open(index_file, 'rb') as f:
        easysci_index = pickle.load(f)
        
    # Extract stranded or unstranded index
    if library_strandedness == 'none':
        gtf_coordinate_index = easysci_index['unstranded']
    elif library_strandedness in ['reverse']: # forward may be eventually added, but since I don't have data available at time of development, skipping
        gtf_coordinate_index = easysci_index['stranded']
    else:
        raise ValueError(f"Experiment strandedness {library_strandedness} not recognized.")

    ### Initialize variables
    # Counters
    i = 0
    read_assignments = []
    start_time = time.time()

    ### Loop through BAM file
    # Loop through paired bundles
    print('Looping through reads...', flush=True)
    for bundle in HTSeq.pair_SAM_alignments(bam_reader, bundle=True):

        # Counter
        i += 1
        if read_subset and i > read_subset:
            break

        # Update counter and report
        if i % 500000 == 0:
            formatted_time = str(timedelta(seconds=int(time.time() - start_time)))
            print(f"Processing line {i:,}... ({formatted_time} elapsed since previous cycle)", flush=True)
            start_time = time.time()

        # Skip multi-mapping reads
        if len(bundle) != 1:
            continue

        # Assign reads
        read1, read2 = bundle[0]

        # Skip lines (unpaired reads, if present)
        if read1 is None or read2 is None or read1.read.name != read2.read.name:
            continue
        
        # Get gene ID assignment
        # Read pair assignnment will be a single string, either a gene ID, a comma-separated gene IDs (if multiple gene adjustment is EM, prop_unique, uniform), or a special string indicating no feature overlap, ambiguous, or multiple (if multimapping is discarded)
        # Closest_TES automatically assigns genes during assignment step, returning just one gene
        # Mapped feature is either exon or transcript, based on the search space used to assign the read pair
        read_pair_assignment, multigene_read, mapped_feature = assign_read_pair_to_gene(read1, read2, gtf_coordinate_index, count_introns=count_introns, primer_type=primer_type, library_strandedness=library_strandedness, multigene_reads=multigene_reads)
       
        # Add to result
        barcode, umi, _ = read1.read.name.split(',')
        read_assignments.append((barcode, umi, read_pair_assignment, mapped_feature, f'{read1.read.name} READ1 chr{read1.iv.chrom}:{read1.iv.start}-{read1.iv.end} READ2 chr{read2.iv.chrom}:{read2.iv.start}-{read2.iv.end}'))

    # Close BAM
    bam_reader.close()

    # Concatenate
    read_assignment_dataframe = pd.DataFrame(read_assignments, columns=['barcode', 'umi', 'gene_id', 'mapped_feature', 'mapping'])
    
    # Prepend sample name to barcode
    read_assignment_dataframe['barcode'] = sample_name + '.' + read_assignment_dataframe['barcode']
    
    # If multi-gene reads are to be resolved by approaches that require the counting to be completed, pass to the resolve function
    if multigene_reads in ['EM', 'prop_unique', 'uniform']:
        read_count_dataframe = resolve_multigene_reads(read_assignment_dataframe, multigene_reads, collapse_by_gene_name)
            
    # Otherwise, count unique UMIs per cell and gene, since the disambiguation has already been resolved above
    elif multigene_reads in ['discard', 'closest_TES']:
        
        # Collapse by gene name
        if collapse_by_gene_name:
            read_count_dataframe = read_assignment_dataframe.merge(easysci_index['gene_names'], on='gene_id').groupby(['barcode', 'gene_name'])['umi'].nunique().reset_index(name='unique_umi_count')
            
        # Collapse by gene ID
        else:
            read_count_dataframe = read_assignment_dataframe.groupby(['barcode', 'gene_id'])['umi'].nunique().reset_index(name='unique_umi_count')
 
    # Create directory
    os.makedirs(output_dir, exist_ok=True)
    assignment_outfile = os.path.join(output_dir, sample_name+'-read_assignment.tsv')
    count_outfile = os.path.join(output_dir, sample_name+'-unique_umi_count.tsv')
    
    # Write output
    print(f'Writing output to {assignment_outfile} ...')
    read_assignment_dataframe.to_csv(assignment_outfile, index=False, sep='\t')
    # read_assignment_dataframe.sort_values('compatible_genes', ascending=False).to_csv(assignment_outfile, index=False, sep='\t')
    
    print(f'Writing output to {count_outfile} ...')
    read_count_dataframe.to_csv(count_outfile, index=False, sep='\t')
    print('Done!')

#############################################
########## 4. Main
#############################################

# Define main function
def main(args):
    
    # Run
    easysci_count(
        input_bam=args.input_bam,
        index_file=args.index_file,
        output_dir=args.output_dir,
        sample_name=args.sample_name,
        library_strandedness=args.library_strandedness,
        count_introns=args.count_introns,
        primer_type=args.primer_type,
        multigene_reads=args.multigene_reads,
        read_subset=args.read_subset,
        collapse_by_gene_name=args.collapse_by_gene_name
    )

#############################################
########## 5. Run
#############################################
# Run main function
if __name__ == "__main__":
    main(args)