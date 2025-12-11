# Import
import glob
import pickle
import pandas as pd

# Get files
barcode_files = glob.glob('/home/torred1/packages/easysci-dt/barcode_files/*.pickle2')

# Loop
for barcode_file in barcode_files:
    
    # Load
    with open(barcode_file, 'rb') as f:
        barcode = pickle.load(f)
        
    # Convert to dataframe
    barcode_dataframe = pd.Series(barcode).to_frame().reset_index().rename(columns={'index': 'barcode_1bp_substitution', 0: 'original_barcode'})
    
    # Get outfile
    outfile = barcode_file.replace('.pickle2', '.tsv')
    
    # Write
    barcode_dataframe.to_csv(outfile, sep='\t', index=False)
