"""
Decode binaries from CDS extracted from synthetic yeast genome and generate music library
"""

import os
import json
from Procodon.decoder import Dnadecoder
from musicpy import *
from fractions import Fraction
import shutil

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

# Split long binary string to 10-bit binary chunks for decoding notes
def split_note_binary_string(
    binary_string: str
):
    """
    Parameters:
        binary_string: binary in string format, e.g. '1010101010....'

    Returns
        note_bin_list: list containing 10-bit binary strings split from input binary string
    """
    
    return [binary_string[i * 10:i * 10 + 10] for i in range(int(len(binary_string)/10))]

# Convert 10-bit binary chunk to note information
def note_bin_to_num(
    note_bin: str
):
    '''
    Parameters:
        note_bin: 10-bit binary in string format, e.g. '1010101010'
    
    Return:
        num: note info in string format, e.g. '5,-1,0,1'
    '''
    
    note_number_bin = note_bin[0:3]
    octave_change_bin = note_bin[3:6]
    half_tone_bin = note_bin[6]
    duration_bin = note_bin[7:10]

    note_number = note_number_decoding_dict[note_number_bin]
    octave_change = octave_change_decoding_dict[octave_change_bin]
    half_tone = accidental_decoding_dict[half_tone_bin]
    duration = duration_decoding_dict[duration_bin]
    
    return ','.join([note_number, octave_change, half_tone, duration])

# Split long binary string & convert to note information
def note_binary_string_to_num_list(
    binary_string: str
):
    """
    Parameters:
        binary_string: binary in string format, e.g. '1010101010....'
    
    Return:
        num_list: list containingf notes info in string format, e.g. ['5,-1,0,1', ...]
    """
    
    note_bin_list = split_note_binary_string(binary_string)
    num_list = []
    for note_bin in note_bin_list:
        num = note_bin_to_num(note_bin)
        if num.split(',')[0] == '0':
            num = '0,0,0,' + num.split(',')[-1]
        num_list.append(num)
    return num_list

def generate_C_major_chord(
    num_list: list,
    speed: float
):
    '''
    Parameters:
        num_list: list containing note info, ['5,-1,0,1', ...]
        speed: playing speed of melodies
    
    Return:
        full_chord: melody object
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

#  Create output directories
if not os.path.exists('synthetic_chr_num_list'):
    os.mkdir('synthetic_chr_num_list')
if not os.path.exists('synthetic_chr_mid'):
    os.mkdir('synthetic_chr_mid')

# Load CDS sequences from json file
with open('syn_chr_CDS_seq_dict.json', 'r') as f:
    CDS_seq_dict = json.loads(f.read())

# Decode sequences by available codecs from yeast
for codec_file in os.listdir('S.cerevisiae.codec/'):
    
    # Initialize Dnadecoder with given codec
    dec = Dnadecoder(codec_path = 'S.cerevisiae.codec/' + codec_file)
    
    # Decode & concatenate binaries
    codec_bin = ''
    for gene in CDS_seq_dict:
        b = dec.decode_cds(CDS_seq_dict[gene])
        codec_bin += b
    
    # Convert binary string to note info list
    num_list = note_binary_string_to_num_list(codec_bin)
    
    # Write note info list into files
    with open('synthetic_chr_num_list/' + codec_file.split('.')[0] + '.num_list.txt', 'w') as hd:
        hd.write('\n'.join(num_list) + '\n')
    
    # Play melody & change file names
    full_chord = generate_C_major_chord(num_list, 1/2)
    play(full_chord)
    shutil.move('temp.mid', 'synthetic_chr_mid/' + codec_file.split('.')[0] + '.mid')
    