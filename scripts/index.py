#############################################
########## 1. Import Libraries
#############################################
import os
import argparse
import HTSeq
import pickle
import pandas as pd

#############################################
########## 2. Processing function
#############################################

# def EasySci_index_v2(input_gtf, output_file):
    
#     # Initialize reader
#     gtf_reader = HTSeq.GFF_Reader(input_gtf)
    
#     # Define strandedness
#     strandedness = ['stranded', 'unstranded']
    
#     # Initialize
#     easysci_index = {}
#     gene_names = {}
    
#     # Print
#     print(f'Creating {stranded} index...', flush=True)
        
#     # Initialize coordinate array for read-gene matching
#     feature_types = ('transcript', 'exon', 'TES')
#     gtf_coordinate_index = {x: HTSeq.GenomicArrayOfSets("auto", stranded=stranded) if x != 'TES' else {} for x in feature_types}
    
#     # Loop through features
#     for feature in gtf_reader:

#         # Check if the feature is either transcript or exon
#         if feature.type in ('transcript', 'exon'):

#             # Add feature coordinates and match with gene ID
#             gtf_coordinate_index[feature.type][feature.iv] += feature.attr['gene_id']

#             # Get endpoints for transcripts
#             if feature.type=='transcript':
                
#                 # Get transcript information
#                 strand = feature.iv.strand
#                 gene_id = feature.attr["gene_id"]

#                 # Get TES based on strand
#                 if strand == "+":  # TES is the end coordinate on the + strand
#                     tes = feature.iv.end
#                 elif strand == "-":  # TES is the start coordinate on the - strand
#                     tes = feature.iv.start
#                 else:
#                     continue  # Skip if no valid strand is defined

#                 # Add gene ID to dictionary keys
#                 if gene_id not in gtf_coordinate_index['TES'].keys():
#                     gtf_coordinate_index['TES'][gene_id] = []
        
#                 # Add coordinates
#                 tes_strand = strand if stranded else '.'
#                 gtf_coordinate_index['TES'][gene_id].append((feature.iv.chrom, tes, tes_strand, feature.attr['transcript_id']))
#         elif feature.type == 'gene' and stranded=='stranded': # save gene IDs and names once to keep for gene_name collapsing later, if specified
#             gene_names[feature.attr['gene_id']] = feature.attr.get('gene_name', feature.attr['gene_id'])
#             # gene_names[feature.attr['gene_id']] = feature.attr['gene_name']
            
#     # Add to the index
#     easysci_index[stranded] = gtf_coordinate_index
    
#     ### Gene names dataframe
#     # Get gene names
#     gene_dataframe = pd.DataFrame(gene_names.items(), columns=['gene_id', 'gene_name'])
    
#     # If the gene name is empty or None, replace with gene ID, using list comprehension
#     gene_dataframe['gene_name'] = [x if x else y for x, y in zip(gene_dataframe['gene_name'], gene_dataframe['gene_id'])]
        
#     # Add gene names to the index, dropping duplicate rows
#     easysci_index['gene_names'] = gene_dataframe.drop_duplicates()
    
#     # Make directory if it doesn't exist
#     os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
#     # Save to file
#     with open(output_file, 'wb') as f:
#         pickle.dump(easysci_index, f)

#     # Print a message to indicate that the index has been created
#     print('Index created successfully!')

def EasySci_index(input_gtf, output_file):
    
    # Initialize reader
    gtf_reader = HTSeq.GFF_Reader(input_gtf)
    
    # Initialize
    easysci_index = {}
    gene_names = {}
    
    # Define strandedness
    strandedness = ['stranded', 'unstranded']
    
    # Do both stranded and unstranded
    for stranded in strandedness:
        
        # Print
        print(f'Creating {stranded} index...', flush=True)
        
        # Initialize coordinate array for read-gene matching
        feature_types = ('transcript', 'exon', 'TES')
        gtf_coordinate_index = {x: HTSeq.GenomicArrayOfSets("auto", stranded=stranded) if x != 'TES' else {} for x in feature_types}
        
        # Loop through features
        for feature in gtf_reader:

            # Check if the feature is either transcript or exon
            if feature.type in ('transcript', 'exon'):

                # Add feature coordinates and match with gene ID
                gtf_coordinate_index[feature.type][feature.iv] += feature.attr['gene_id']

                # Get endpoints for transcripts
                if feature.type=='transcript':
                    
                    # Get transcript information
                    strand = feature.iv.strand
                    gene_id = feature.attr["gene_id"]

                    # Get TES based on strand
                    if strand == "+":  # TES is the end coordinate on the + strand
                        tes = feature.iv.end
                    elif strand == "-":  # TES is the start coordinate on the - strand
                        tes = feature.iv.start
                    else:
                        continue  # Skip if no valid strand is defined

                    # Add gene ID to dictionary keys
                    if gene_id not in gtf_coordinate_index['TES'].keys():
                        gtf_coordinate_index['TES'][gene_id] = []
            
                    # Add coordinates
                    tes_strand = strand if stranded else '.'
                    gtf_coordinate_index['TES'][gene_id].append((feature.iv.chrom, tes, tes_strand, feature.attr['transcript_id']))
            elif feature.type == 'gene' and stranded=='stranded': # save gene IDs and names once to keep for gene_name collapsing later, if specified
                gene_names[feature.attr['gene_id']] = feature.attr.get('gene_name', feature.attr['gene_id'])
                # gene_names[feature.attr['gene_id']] = feature.attr['gene_name']
                
        # Add to the index
        easysci_index[stranded] = gtf_coordinate_index
        
    ### Gene names dataframe
    # Get gene names
    gene_dataframe = pd.DataFrame(gene_names.items(), columns=['gene_id', 'gene_name'])
    
    # If the gene name is empty or None, replace with gene ID, using list comprehension
    gene_dataframe['gene_name'] = [x if x else y for x, y in zip(gene_dataframe['gene_name'], gene_dataframe['gene_id'])]
        
    # Add gene names to the index, dropping duplicate rows
    easysci_index['gene_names'] = gene_dataframe.drop_duplicates()
    
    # Make directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save to file
    with open(output_file, 'wb') as f:
        pickle.dump(easysci_index, f)

    # Print a message to indicate that the index has been created
    print('Index created successfully!')
              
#############################################
########## 3. Main
#############################################

# Define main function
def main(args):
    
    # Run
    EasySci_index(
        input_gtf = args.input_gtf,
        output_file = args.output_file
    )

#############################################
########## 4. Run
#############################################
# Run main function
if __name__ == "__main__":
    main(args)
