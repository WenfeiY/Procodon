"""
Test run time for decoding different number of genes
"""

from Procodon.seqloader import Seqloader
from Procodon.decoder import Dnadecoder
import time
import sys

#  Load Protein and CDS sequences
seqload = Seqloader(seq_dir = 'S.cerevisiae')

#  Convert keys to gene list
gene_list = list(seqload.gene_seq_dict.keys())

def count_decoding_time(
    num: int
):
    """
    Test the time for decoding genes with given number

    Parameter:
        num: gene number, int
    
    Output (print):
        num: gene number
        elapsed_time: decoding time
    """
    #  Time start
    start_time = time.perf_counter()
    #  Extract the subset from full CDS list
    sub_gene_list = gene_list[:num]
    cds_list = seqload.get_cds_seqs(sub_gene_list)
    #  Initialize Dnadecoder
    gene_decode = Dnadecoder()
    #  Decode CDSs
    complete_binary_string = gene_decode.decode_cds_list(cds_list)
    #  Write binary string into file
    with open('tmp', 'w') as hd:
        hd.write(complete_binary_string)
    #  Time end
    end_time = time.perf_counter()
    elapsed_time = (end_time - start_time)
    print(num)
    print(elapsed_time)

if __name__ == '__main__':
    #  Get gene number parameter
    gene_num = int(sys.argv[1])
    count_decoding_time(gene_num)
