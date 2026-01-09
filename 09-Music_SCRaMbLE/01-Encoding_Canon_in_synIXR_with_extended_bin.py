"""
Encode Canon in synIXR with extended binaries (150 bit)
"""

import random
from Procodon.seqloader import Seqloader
from Procodon.codec.codec import Codecloader
from Procodon.dnaencoder import Dnaencoder
from Procodon.utils import count_available_space
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
import json
import os

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

def int_to_binary(
    n: int, 
    bits: int = 15
) -> str:
    """
    Convert integer to binary string with given length
    
    Parameters:
        n: integer, int
        bits: length of binary string, int
    
    Return:
        binary_string: result binary string, str
    """
    #  Calculate length of mask
    mask = (1 << bits) - 1
    
    n_mod = n & mask
    
    #  Generate binary string
    return format(n_mod, f'0{bits}b')

def generate_random_num_0_1(
    b: str, 
    bits: int = 15
):
    """
    Generate extended binary based on given binary bit and extended length
    
    Parameters:
        b: 1-bit binary, str
        bits: length of extended binary string, int
    
    Return:
        extended_binary: extended binary string, str
    """
    if b == '0':
        bot = 0
        top = 2** (bits - 1) - 1
    else:
        bot = 2** (bits - 1)
        top = 2** bits - 1
    
    random_number = random.randint(bot, top)
    extended_binary = int_to_binary(random_number, bits = bits)
    return extended_binary

def extend_binary(
    binary_string: str, 
    bits: int = 15
):
    """
    Generate extended binaries based on given binary string and extended length, join to form extended binary string
    
    Parameters:
        binary_string: provided binary string, str
        bits: length of extended binary string, int
    
    Return:
        extended_binary: extended binary string, str
    """
    
    extended_binary_string = ''
    
    for b in binary_string:
        extended_binary_string += generate_random_num_0_1(b, bits = bits)
    return extended_binary_string



#  Load CDSs and protein sequences
seqload= Seqloader(seq_dir = 'S.cerevisiae')

#  Load codec
codecload=Codecloader()
codecload.load_codec_from_file('codec_112417.json')

#  Initialize Dnaencoder
encode = Dnaencoder(codecload.encoding_dict, 'S.cerevisiae.codon_frequency.csv')
#  Load genes
with open('synIXR_chrR_gene_list.txt', 'r') as f:
    gene_list = f.read().splitlines()

#  Read notes
with open('Canon.synIXR.150_bit_per_note.txt', 'r') as f:
    num_list = f.read().splitlines()
#  Convert notes to binary string
binary_string = num_list_to_binary(num_list)

#  Extend the decoded binary string
extended_binary_string = extend_binary(binary_string)


#  Encode extended binary string into genes in synIXR. Only complete notes are encoded.
idx = 0
with open('Recoded_synIXR_Canon_gene_150bit.fa', 'w') as hd:
    for gene in gene_list:
        data_len = int(count_available_space(seqload.gene_seq_dict[gene]['Protein'])/150) * 150
        if data_len == 0:
            continue
        binary = extended_binary_string[idx : idx + data_len]
        recoded_seq, _, _, _ = encode.binary_to_DNA(seqload.gene_seq_dict[gene]['CDS'], seqload.gene_seq_dict[gene]['Protein'], binary)
        hd.write('>' + gene + '\n' + recoded_seq + '\n')
        idx += data_len
