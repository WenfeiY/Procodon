"""
Encode text data into single sequence
"""

import random
from Procodon.seqloader import Seqloader
from Procodon.codec.codec import Codecloader
from Procodon.dnaencoder import Dnaencoder
from Procodon.utils import count_mutation
from Procodon.utils import count_available_space
from Bio.SeqUtils import gc_fraction

# Dictionary for converting note number (Based on C major) to 3-bit binary
note_number_encoding_dict = {'0': '000',
                             '1': '001',
                             '2': '010',
                             '3': '011',
                             '4': '100',
                             '5': '101',
                             '6': '110',
                             '7': '111'
                             }

# Dictionary for converting octave change to 3-bit binary
octave_change_encoding_dict = {'-3': ['000'],
                               '-2': ['001'],
                               '-1': ['010'],
                               '0': ['011', '100'],
                               '1': ['101'],
                               '2': ['110'],
                               '3': ['111']
                               }

# Dictionary for converting accidental change (half-tone) to 1-bit binary
accidental_encoding_dict = {'0': '0',
                            '#': '1'
                            }

# Dictionary for converting duration (1 represent one whole note) to 3-bit binary
duration_encoding_dict = {'1/8': '000',
                          '1/4': '001',
                          '1/2': '010',
                          '3/4': '011',
                          '1': '100',
                          '3/2': '101',
                          '2': '110',
                          '3': '111'
                          }

def num_to_binary(
    num: str
):
    """
    Convert note to 10-bit binary string
    
    Parameter:
        num: note info in string format, e.g. '5,-1,0,1'
    
    Return:
        note_binary: 10-bit binary representing 1 note
    """
    
    #  Split note info
    note_number = num.split(',')[0]
    octave_change = num.split(',')[1]
    half_tone = num.split(',')[2]
    duration = num.split(',')[3]
    #  Convert info of each part to binary
    bin1 = note_number_encoding_dict[note_number]
    bin2 = random.choice(octave_change_encoding_dict[octave_change])
    bin3 = accidental_encoding_dict[half_tone]
    bin4 = duration_encoding_dict[duration]
    
    return bin1 + bin2 + bin3 + bin4

def num_list_to_binary(
    num_list: list
):
    """
    Convert note info list to string
    
    Parameter:
        num_list: list containing note infos, list, e.g. ['5,-1,0,1', ...]
    
    Return:
        binary_string: binary string generated from num_list
    """
    binary_string = ''
    
    for num in num_list:
        binary_string += num_to_binary(num)
    
    return binary_string

#  Load CDSs and protein sequences
seqload= Seqloader(seq_dir = 'S.cerevisiae')

#  Load codec
codecload=Codecloader()
codecload.load_codec_from_file('codec_112417.json')

#  Initialize Dnaencoder
encode = Dnaencoder(codecload.encoding_dict, 'S.cerevisiae.codon_frequency.csv')

#  File containing notes
music_file = 'The_blue_danube.YIL036WC.txt'
gene = 'YIL036W'  #  Can be replaced to other genes

#  Read music data
with open(music_file, 'r') as f:
    num_list = f.read().splitlines()

#  Convert notes to binary data
binary_data = num_list_to_binary(num_list)

#  Add melody end in binary data
avil_space = count_available_space(seqload.gene_seq_dict[gene]['Protein'])
if avil_space - len(binary_data) >= 10:
    binary_string = binary_data + '0000000000'
else:
    binary_string = binary_data

recoded_seq, _, _, _ = encode.binary_to_DNA(seqload.gene_seq_dict[gene]['CDS'], seqload.gene_seq_dict[gene]['Protein'], binary_string)

with open(gene + '.recoded_info.tsv', 'w') as hd:
    hd.write(gene + '\t' +
             recoded_seq + '\t' + 
             str(count_mutation(recoded_seq, seqload.gene_seq_dict[gene]['CDS'])) + '\t' + 
             str(round(gc_fraction(seqload.gene_seq_dict[gene]['CDS']) * 100, 2)) + '\t' +
             str(round(gc_fraction(recoded_seq) * 100, 2)) + '\n'
             )
