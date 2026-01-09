"""
Encode Simple notes into sequences
"""

from Procodon.seqloader import Seqloader
from Procodon.codec.codec import Codecloader
from Procodon.dnaencoder import Dnaencoder
from Procodon.fileloader import Fileloader
from Procodon.utils import count_mutation
from Bio.SeqUtils import gc_fraction

#  Load CDSs and protein sequences
seqload= Seqloader(seq_dir = 'S.cerevisiae')

#  Load codec
codecload=Codecloader()
codecload.load_codec_from_file('codec_112417.json')

#  Initialize Dnaencoder
encode = Dnaencoder(codecload.encoding_dict, 'S.cerevisiae.codon_frequency.csv')

#  Data and gene for encoding
data = '1 1 5 5 6 6 5\n4 4 3 3 2 2 1\n5 5 4 4 3 3 2\n5 5 4 4 3 3 2\n1 1 5 5 6 6 5\n4 4 3 3 2 2 1'  #  Simple notes from "Twinkle Twinkle Little Star"
gene_list = ['YAL003W', 'YAL032C', 'YAL041W']  #  Can be replaced to other genes

# Initialize Fileloader
fileload = Fileloader()
#  Convert text to binary data
binary_data = fileload.text_to_bin(data)

#  Recode sequences
rest_bin = binary_data
with open('Three_genes.recoded_info.tsv', 'w') as hd:
    for gene in gene_list:
        recoded_seq, _, _, rest_bin = encode.binary_to_DNA(seqload.gene_seq_dict[gene]['CDS'], seqload.gene_seq_dict[gene]['Protein'], rest_bin)
        if rest_bin == '':
            break
        hd.write(gene + '\t' +
                 recoded_seq + '\t' + 
                 str(count_mutation(recoded_seq, seqload.gene_seq_dict[gene]['CDS'])) + '\t' + 
                 str(round(gc_fraction(seqload.gene_seq_dict[gene]['CDS']) * 100, 2)) + '\t' +
                 str(round(gc_fraction(recoded_seq) * 100, 2)) + '\n'
                 )
