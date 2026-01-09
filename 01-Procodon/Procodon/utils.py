"""
Utilities for Procodon module
"""
import math
from PIL import Image
import numpy as np 
from Bio.SeqUtils import gc_fraction
import itertools

def split_CDS(
    cds_seq: str
):
    
    """
    Split CDS sequence into 3-latter codons, return codon list
    
    Parameter:
        cds_seq: CDS sequence, str
    
    Return:
        codon_list: list containing 3-latter codons, list[str]
    """
    
    cds_seq = cds_seq.upper()
    
    if not len(cds_seq) % 3 == 0:
        print('The CDS sequence length is not a multiple of 3, please check it!')
        return
    else:
        return [cds_seq[i:i+3] for i in range(0, len(cds_seq), 3)]

def get_empty_codon_dict():
    
    """
    Generate empty_codon_dict
    
    Return:
        empty_codon_dict:
            {
                'AAA': 0,
                'AAT': 0,
                'AAC': 0,
                ...
            }
    """
    
    empty_codon_dict = {}
    
    for a in ['A', 'T', 'C', 'G']:
        for b in ['A', 'T', 'C', 'G']:
            for c in ['A', 'T', 'C', 'G']:
                empty_codon_dict[a + b + c] = 0
    
    return empty_codon_dict

def get_empty_select_dict():

    """
    Generate empty_select_dict
    
    Return:
        empty_select_dict: dictionary for codon selection
            {
                'F': {'0': [], '1': []},
                'L': {'0': [], '1': []},
                'I': {'0': [], '1': []},
                ...
            }
    """
    
    empty_select_dict = {
        'F': {'0': [], '1': []},
        'L': {'0': [], '1': []},
        'I': {'0': [], '1': []},
        'V': {'0': [], '1': []},
        'S': {'0': [], '1': []},
        'P': {'0': [], '1': []},
        'T': {'0': [], '1': []},
        'A': {'0': [], '1': []},
        'Y': {'0': [], '1': []},
        'H': {'0': [], '1': []},
        'Q': {'0': [], '1': []},
        'N': {'0': [], '1': []},
        'K': {'0': [], '1': []},
        'D': {'0': [], '1': []},
        'E': {'0': [], '1': []},
        'C': {'0': [], '1': []},
        'R': {'0': [], '1': []},
        'G': {'0': [], '1': []}
        }
    return empty_select_dict

def int_len(
    number: int
) -> int:
    
    """
    Calculate the length of the given integer
    
    Parameter:
        number: given integer, int
    
    Return:
        number_len: length of the given integer, int
    """
    
    return int(math.log10(number)) + 1

def text_to_bin(
    text,
    encoding: str = 'utf-8'
) -> str:
    
    """
    Convert given text to binary string
    
    Parameters:
        text: given text for encoding, str
        encoding: encoding format for computer, default UTF-8
    
    Return:
        binary_string: binary string for encoding, str, e.g. '01011010101...'
    """
    
    bytes_data = text.encode(encoding)
    binary_string = ''.join(format(byte, '08b') for byte in bytes_data)
    
    return binary_string

def binary_array_to_text(
    binary_array
):
    
    """
    Transform binary array data into text. Stop at error position.
    
    Parameter:
        binary_array: bytearray data, bytearray
    
    Return:
        text: decoded text, str
    """
    
    pos = len(binary_array)
    while pos >= 0:
        try:
            text = binary_array[:pos].decode('utf-8', errors='strict')
            return text
        except UnicodeDecodeError as e:
            pos = e.start
            if pos == 0:
                return ''


def normalize_png(
    in_file_path: str,
    out_file_path: str,
    cut_off: int
):
    
    """
    Transform grey PNG image into black-and-white PNG image with cut-off to define black or white (0~255), write into file
    
    Parameters:
        in_file_path: PNG input file path, str
        out_file_path: normalized PNG output file path, str
        cut_off: integer for defining black or white (0~255), int
    
    Return:
        normalized image file
    """
    
    image = Image.open(in_file_path)
    image_array = np.array(image)
    pixels = image_array.tolist()
    coords = []
    for row in range(len(pixels)):
        for col in range(len(pixels[0])):
            val = pixels[row][col]
            if val < cut_off:
                coords.append((col, row))
     
    # width and height of PNG image
    width, height = len(pixels[0]), len(pixels)
     
    # create new gray image
    image = Image.new('L', (width, height), color='white')

    for x, y in coords:
        image.putpixel((x, y), 0)
     
    # save image
    image.save(out_file_path, lossless=True)

def write_bin_to_file(
    binary_array,
    out_file_path: str
):
    
    """
    Write byte array data into file
    
    Parameters:
        binary_array: bytearray data, bytearray
        out_file_path: binary array data output file path, str
    
    Return:
        output file with binary data
    """
    
    with open(out_file_path, 'wb') as output_file:
        # write binary array data into file
        output_file.write(binary_array)

def count_identity(
    seq1: str,
    seq2: str
) -> float:
    
    """
    Calculate the identity of 2 DNA sequences by comparing these sequences with same length, base-by-base
    
    Parameters:
        seq1: DNA sequence 1, str
        seq2: DNA sequence 2, str
    
    Return:
        identity: sequence identity between 2 sequences, float
    """
    
    if len(seq1) == len(seq2):
        seq1 = seq1.upper()
        seq2 = seq2.upper()
        count = 0
        for i in range(len(seq1)):
            if seq1[i] == seq2[i]:
                count += 1
    else:
        raise ValueError(
            '2 sequences must be qeual in length!'
        )
    return count / len(seq1)

def write_fa(
    output_path: str, 
    seq_id: str, 
    seq: str
):
    """
    Write sequence into file in FASTA format

    Parameters:
        output_path: file path for writing, str
        seq_id: sequence ID, str
        seq: sequence, str
    """
    
    with open(output_path, 'w') as hd:
        hd.write('>' + seq_id + '\n' + seq + '\n')

def write_multi_fa(
    output_path: str, 
    seq_id: str, 
    seq_list: list
):
    """
    Write sequence into file in FASTA format

    Parameters:
        output_path: file path for writing, str
        seq_id: sequence ID, str
        seq_list: list containing sequences in str
    """

    seq_num = 1
    with open(output_path, 'w') as hd:
        for seq in seq_list:
            hd.write('>' + seq_id + '_' + str(seq_num) + '\n' + seq + '\n')
            seq_num += 1

def count_available_space(
    prot_seq: str
) -> int:
    
    """
    calculate number of amino acids available for encoding

    Parameter:
        prot_seq: protein sequence, str
    
    Return:
        count: available encoding space (bit), int
    """
    
    count = 0
    for aa in prot_seq:
        if not aa == '*' and not aa == 'M' and not aa == 'W':
            count += 1
    return count

def count_GC_diff(
    query: str, 
    subject: str
):
    """
    Calculate the difference between 2 DNA sequences

    Parameters:
        query: query DNA sequence, str
        subject: subject DNA sequence, str
    
    Return:
        delta_GC: GC deviation, float
    """

    GC1 = round(gc_fraction(query) * 100, 2)
    GC2 = round(gc_fraction(subject) * 100, 2)
    return GC1 - GC2


def binary_string_to_byte_array(
    binary_string: str
):
    """
    Convert binary string to byte array

    Parameter:
        binary_string: binary string, str
    
    Return:
        byte_array: byte array, bytearray
    """

    length = int(len(binary_string) / 8)
    byte_array = bytearray()
    for i in range(length):
        byte_chunk = binary_string[i * 8 : (i + 1) * 8]
        byte_value = int(byte_chunk, 2)
        byte_array.append(byte_value)
    
    return byte_array

def reverse_bin(
    binary_string: str
):
    """
    Convert binary string to 0/1 switched binary string

    Parameter:
        binary_string: binary string, str
    
    Return:
        rev_binary_string: 0/1 switched binary string, str
    """
    
    switch_dict = {'0':'1','1':'0'}
    rev_binary_string = ''
    for b in binary_string:
        rev_binary_string += switch_dict[b]
    
    return rev_binary_string

def count_mutation(
    recode_seq: str, 
    cds_seq: str
):
    """
    Calculate mutation rate of 2 given sequences

    Parameters:
        recode_seq: recoded CDS sequence, str
        cds_seq: original CDS sequence, str
    
    Return:
        mutation_rate: persentage of mutation rate, float
    """
    
    m = 0
    for i in range(len(recode_seq)):
        if not recode_seq[i].upper() == cds_seq[i].upper():
            m += 1
    return round(m / len(recode_seq) * 100, 2)

def binary_to_text(
    binary_string: str
):
    """
    Convert binary string to text

    Parameter:
        binary_string: binary string, str
    
    Return:
        text: text data, str
    """
    
    chunk_size = 8
    #  Convert binary string to byte list
    bytes_list =  [int(''.join(chunk), 2) for chunk in itertools.zip_longest(*[iter(binary_string)]*chunk_size, fillvalue='')]
    act_bytes_list = bytes_list
    byte_array = bytearray(act_bytes_list)
    #  Convert byte list to text. Stop when meet error.
    pos = len(byte_array)
    while pos >= 0:
        try:
            text = byte_array[:pos].decode('utf-8', errors='strict')
            return text
        except UnicodeDecodeError as e:
            pos = e.start
            if pos == 0:
                return ''
