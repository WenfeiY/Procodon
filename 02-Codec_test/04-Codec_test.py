"""
Test codec, calculate GC bias and write into json files
"""

from Procodon.seqloader import Seqloader
from Procodon.dnaencoder import Dnaencoder
from Procodon.utils import count_GC_diff
import json
import multiprocessing
import os
import sys

# Read species prefix
prefix = sys.argv[1]
species_type = sys.argv[2]
csv_filename = prefix + '.codon_frequency.csv'

#Load CDS and protein sequences
seqload = Seqloader(seq_dir = prefix)

#Load codecs and binary test data
with open(prefix + '.codec_pass_qua_dict.json', 'r') as f:
    codec_dict = json.loads(f.read())

with open('../' + prefix + '.binary_test_data_out.json', 'r') as f:
    test_data_dict = json.loads(f.read())

# Process each gene
def process_gene(
    gene_data
):
    # Get gene name, codec (encoding dictionary) and codon usage path
    gene, encoding_dict, csv_filename = gene_data
    
    # Initialize encoding program
    encode = Dnaencoder(encoding_dict, csv_filename, species_type = species_type)
    cds_seq = seqload.gene_seq_dict[gene]['CDS']
    prot_seq = seqload.gene_seq_dict[gene]['Protein']
    test_data = test_data_dict[gene]
    
    # Recode genes and calculate GC bias, store in dictionary
    gene_result = {}
    for prop_1 in test_data:
        gene_result[prop_1] = []
        for binary_string in test_data[prop_1]:
            
            # Recode DNA
            recoded_seq, _, _, _ = encode.binary_to_DNA(cds_seq, prot_seq, binary_string)
            
            # Calculate GC bias
            gc_bias = count_GC_diff(recoded_seq, cds_seq)
            gene_result[prop_1].append(gc_bias)
    
    return gene, gene_result


if __name__ == '__main__':
    
    if not os.path.exists(prefix + '.Codec.test/'):
        os.mkdir(prefix + '.Codec.test/')
    if not os.path.exists(prefix + '.GC_bias/'):
        os.mkdir(prefix + '.GC_bias/')
    
    i = 1
    for num in codec_dict:
        encoding_dict = codec_dict[num]
        recoded_gc_bias_dict = {}
        
        # Prepare tasks
        tasks = [
            (gene, encoding_dict, csv_filename)
            for gene in seqload.gene_seq_dict
        ]
        
        # Process by multiple thread pool
        with multiprocessing.Pool(processes=os.cpu_count() - 5) as pool:
            results = pool.map(process_gene, tasks)
            
            # Collect result
            for gene, gene_result in results:
                recoded_gc_bias_dict[gene] = gene_result
        
        with open(prefix + '.Codec.test/' + num + '.codec_test.gc_bias.json', 'w') as hd:
            json.dump(recoded_gc_bias_dict, hd)
        
        os.system("python 01-Summarizing_GC_bias.py " + prefix + ".Codec.test/" + num + ".codec_test.gc_bias.json " + prefix + ".GC_bias/" + num + '.gc_bias.tsv &')
        if i % 50 == 0:
            print(f"{i} codecs complete")
