"""
Decode melofies from SCRaMbLEd strains with extended binaries (150 bit)
"""

from Procodon.decoder import Dnadecoder
from Bio import SeqIO
import os
import shutil
from musicpy import *
from fractions import Fraction

# Dictionary for converting 3-bit binary to note number (Based on C major)
note_number_decoding_dict = {'000': '0',  # Rest note
                             '001': '1',  # 'C'
                             '010': '2',  # 'D'
                             '011': '3',  # 'E'
                             '100': '4',  # 'F'
                             '101': '5',  # 'G'
                             '110': '6',  # 'A'
                             '111': '7'  # 'B'
                             }

# Dictionary for converting 3-bit binary to octave change
octave_change_decoding_dict = {'000': '-3',
                               '001': '-2',
                               '010': '-1',
                               '011': '0',
                               '100': '0',
                               '101': '1',
                               '110': '2',
                               '111': '3'
                               }

# Dictionary for converting 1-bit binary to accidental change (half-tone)
accidental_decoding_dict = {'0': '0',
                            '1': '#'
                            }

# Dictionary for converting 3-bit binary to duration (1 represent one whole note)
duration_decoding_dict = {'000': '1/8',
                          '001': '1/4',
                          '010': '1/2',
                          '011': '3/4',
                          '100': '1',
                          '101': '3/2',
                          '110': '2',
                          '111': '3'
                          }

# Dictionary for converting note number to note name (Based on C major)
c_major_note_num_dict = {'1': 'C',
                         '2': 'D',
                         '3': 'E',
                         '4': 'F',
                         '5': 'G',
                         '6': 'A',
                         '7': 'B'
                         }

def split_note_binary_string(
    binary_string: str
):
    """
    Split binary string to 10-bit chunks
    
    Parameter:
        binary_string: decoded binary string, str
    
    Return:
        chunk_list: list containing 10-bit chunks
    """
    return [binary_string[i * 10:i * 10 + 10] for i in range(int(len(binary_string)/10))]

def note_bin_to_num(
    note_bin: str
):
    """
    Convert note binary to note info

    Parameter:
        note_bin: 10-bit binary for 1 note, str
    
    Return:
        note_info: list containing note number, octave change, accidental change and duration
    """
    
    #  Split binaries
    note_number_bin = note_bin[0:3]
    octave_change_bin = note_bin[3:6]
    half_tone_bin = note_bin[6]
    duration_bin = note_bin[7:10]
    
    #  Convert data by dictionary
    note_number = note_number_decoding_dict[note_number_bin]
    octave_change = octave_change_decoding_dict[octave_change_bin]
    half_tone = accidental_decoding_dict[half_tone_bin]
    duration = duration_decoding_dict[duration_bin]
    
    return ','.join([note_number, octave_change, half_tone, duration])


def note_binary_string_to_num_list(
    binary_string: str
):
    """
    Split binary string and convert note binaries to note infos

    Parameter:
        binary_string: provided binary string, str
    
    Return:
        num_list: list containing note infos
    """
    
    #  Split binary string
    note_bin_list = split_note_binary_string(binary_string)
    #  Convert note binaries to note infos
    num_list = []
    for note_bin in note_bin_list:
        num = note_bin_to_num(note_bin)
        if num.split(',')[0] == '0' and not num.split(',')[1] == '0':
            return num_list
        else:
            num_list.append(note_bin_to_num(note_bin))
    
    return num_list

def generate_C_major_chord(
    num_list: list,
    speed: float
):
    '''
    Input:
        num_list: list containing note info, ['5,-1,1']
    
    '''
    note_list = []
    dur_list = []
    interval_list = []
    if num_list[0].split(',')[0] == '0':
        start_time = float(Fraction(num_list[0].split(',')[3]))*speed
        num_list = num_list[1:]
    else:
        start_time = 0
    for num in num_list:
        
        alphabet_number = num.split(',')[0]
        half_tone = num.split(',')[2]
        octave_change = int(num.split(',')[1])
        duration = speed * float(Fraction(num.split(',')[3]))
        if alphabet_number == '0':
            if interval_list:
                interval_list[-1] += duration
        else:
            dur_list.append(duration)
            interval_list.append(duration)
            alphabet = c_major_note_num_dict[alphabet_number]
            if half_tone == '#':
                alphabet += '#'
            octave = 4 + octave_change
            alphabet += str(octave)
            note_list.append(alphabet)

    if note_list:
        full_chord = chord(note_list, duration = dur_list, interval = interval_list, start_time = start_time)
    return full_chord


#  Default parameters
codec_path = 'Codec_112417.json'
chunk_size = 150
bin_size = 15
recoded_seq_path = 'Recoded_synIXR_Canon_gene_150bit.fa'
synIXR_structure = 'synIXR.LoxPsym_unit_structure.txt'
SCRaMbLEd_structure_dir = 'synIXR.SCRaMbLE.structure'
SCRaMbLEd_gene_sort_dir = 'synIXR.SCRaMbLE.gene.sort'
SCRaMbLEd_note_info_dir = 'synIXR.SCRaMbLE.note_info'
SCRaMbLEd_melody_dir = 'synIXR.SCRaMbLE.melody'

#  Read recoded sequences
recoded_seq_dict = {}
for record in SeqIO.parse(recoded_seq_path, 'fasta'):
    recoded_seq_dict[record.id] = str(record.seq)

#  Decode CDS
decoded_bin_dict = {}
dec = Dnadecoder(codec_path = codec_path)
for gene in recoded_seq_dict:
    gene_binary = dec.decode_cds(recoded_seq_dict[gene])
    decoded_bin_dict[gene] = gene_binary

#  Compress binary
bin_dict = {}
for gene in decoded_bin_dict:
    gene_bin = ''
    b = decoded_bin_dict[gene]
    notes = [b[i:i+chunk_size] for i in range(0, len(b), chunk_size) if i + chunk_size <= len(b)]
    if notes:
        for note in notes:
            bins = [note[i:i+bin_size] for i in range(0, len(note), bin_size) if i + bin_size <= len(note)]
            for bs in bins:
                if bs[0] == '1':
                    gene_bin += '1'
                else:
                    gene_bin += '0'
        bin_dict[gene] = gene_bin

#  Parse chromosome structure
lu_gene_dict = {}
with open(synIXR_structure, 'r') as f:
    lines = f.read().splitlines()
    for i in range(len(lines)):
        lu_gene_dict[str(i + 1)] = lines[i].split(' ')
        lu_gene_dict['-' + str(i + 1)] = lines[i].split(' ')[::-1]

#  Parse SCRaMbLEd structure, extract gene order and direction
for structure_file in os.listdir(SCRaMbLEd_structure_dir):
    with open(SCRaMbLEd_structure_dir + '/' + structure_file, 'r') as f:
        lus = f.read().splitlines()
    gene_list = []
    for lu in lus:
        if lu[0] == '-':
            for gene in lu_gene_dict[lu]:
                if gene:
                    gene_list.append('-' + gene)
        else:
            for gene in lu_gene_dict[lu]:
                if gene:
                    gene_list.append(gene)

    with open(SCRaMbLEd_gene_sort_dir + '/' + structure_file.split('.')[0] + '.lu_gene.txt', 'w') as hd:
        hd.write('\n'.join(gene_list) + '\n')

#  Read gene list
SCRaMbLEd_gene_dict = {}

for gene_list_file in os.listdir(SCRaMbLEd_gene_sort_dir):
    with open(SCRaMbLEd_gene_sort_dir + '/' + gene_list_file, 'r') as f:
        gene_list = f.read().splitlines()
        SCRaMbLEd_gene_dict[gene_list_file.split('.')[0]] = gene_list

#  Convert to note infos
for strain in SCRaMbLEd_gene_dict:
    num_list = []
    for gene in SCRaMbLEd_gene_dict[strain]:
        if gene[0] == '-':
            if gene[1:] in bin_dict:
                num_list += note_binary_string_to_num_list(bin_dict[gene[1:]])[::-1]
        else:
            if gene in bin_dict:
                num_list += note_binary_string_to_num_list(bin_dict[gene])
    with open(SCRaMbLEd_note_info_dir + '/' + strain, 'w') as hd:
        hd.write('\n'.join(num_list) + '\n')

#  Generate melodies
for note_info_file in os.listdir(SCRaMbLEd_note_info_dir):
    with open(SCRaMbLEd_note_info_dir + '/' + note_info_file, 'r') as f:
        num_list = f.read().splitlines()
        full_chord = generate_C_major_chord(Codec_strain_num_list_dict[strain], 1/2)
        play(full_chord)

        shutil.move('temp.mid', SCRaMbLEd_melody_dir + '/' + strain + '.mid')
