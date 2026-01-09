"""
Encode text data into single sequence
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
data = 'God does not play dice.'  #  Can be replaced to other texts
gene = 'YAL003W'  #  Can be replaced to other genes

# Initialize Fileloader
fileload = Fileloader()
#  Convert text to binary data
binary_data = fileload.text_to_bin(data)

#  Recode sequence
recoded_seq, _, _, _ = encode.binary_to_DNA(seqload.gene_seq_dict[gene]['CDS'], seqload.gene_seq_dict[gene]['Protein'], binary_data)

with open(gene + '.recoded_info.tsv', 'w') as hd:
    hd.write(gene + '\t' +
             recoded_seq + '\t' + 
             str(count_mutation(recoded_seq, seqload.gene_seq_dict[gene]['CDS'])) + '\t' + 
             str(round(gc_fraction(seqload.gene_seq_dict[gene]['CDS']) * 100, 2)) + '\t' +
             str(round(gc_fraction(recoded_seq) * 100, 2)) + '\n'
             )
