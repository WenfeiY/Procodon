"""
Generating potential available codecs and write into files for dowmstream codec test.
"""

import json
from Procodon.codec.codec import Codecgenerator

def cal_codec_num_and_write(
    sp: str, 
    codec_0_path: str, 
    codon_freq_path: str, 
    aa_freq_path: str, 
    qua_dict_path: str, 
    GBS_path: str
):
    """
    Filter potential available codecs of species;
    Print available codec number;
    Write the first, first-quarter, middle, three-quarter and last GBS-pass codecs sorted by GBS;
    Write GBS of filtered codecs.
    
    Parameters:
        sp: species prefix
        codec_0_path: file path of the original codec
        codon_freq_path: file path of codon frequencies in Comma-Separated Values (CSV) format
        aa_freq_path: file path of the total number of each amino acid accross all proteins
        qua_dict_path: output files path of the first, first-quarter, middle, three-quarter and last GBS-pass codecs sorted by GBS
        GBS_path: output files path of GBS
    """  
    
    #  Initialize Codecgenerator
    codec_gen = Codecgenerator(codec_0_path, codon_freq_path, aa_freq_path)
    
    #  Filter potential available codecs
    combination_pass, prop_pass = codec_gen.filter_GSC(1)
    
    #  Sort codecs by GBS
    combination_pass_dict = {}
    for comb_info in combination_pass:
        combination_pass_dict[comb_info[0]] = comb_info[1]
    prop_pass_sorted = sorted(prop_pass, key=lambda x: abs(x[1] - 0.5))
    
    #  Write GBS into file
    with open(GBS_path, 'w') as hd:
        for prop in prop_pass_sorted:
            hd.write(str(prop[0]) + '\t' + str(abs(prop[1] - 0.5)) + '\n')
    
    #  Output species and number of available codecs
    print(sp)
    print(len(prop_pass))
    
    #  Extract the first, first-quarter, middle, three-quarter and last GBS-pass codecs sorted by GBS
    codec_sort = [1]
    for x in (0.25, 0.5, 0.75, 1):
        codec_sort.append(int(len(prop_pass)*x))
    
    qua_dict = {}
    for codec_pos in codec_sort:
        codec_num = prop_pass_sorted[codec_pos - 1][0]
        print(str(codec_pos) + '\t' + str(codec_num))
        
        codec = codec_gen.generate_codec(combination_pass_dict[codec_num])
        qua_dict[codec_num] = codec
    
    #  Write the first, first-quarter, middle, three-quarter and last GBS-pass codecs into file
    with open(qua_dict_path, 'w') as hd:
        json.dump(qua_dict, hd)


for sp in ['M.vulcanius', 'M.formicicum', 'P.syntrophicum', 'T.sedimenti', 'E.coli', 'S.collinus', 'S.cerevisiae', 'P.patens', 'H.sapiens']:
    codec_0_path = sp + '.codec_0.json'
    codon_freq_path = sp + '.codon_frequency.csv'
    aa_freq_path = sp + '.aa_frequency_dict.json'
    qua_dict_path = sp + '.codec_pass_qua_dict.json'
    GBS_path = sp + '.GBS.txt'
    cal_codec_num_and_write(sp, codec_0_path, codon_freq_path, aa_freq_path, qua_dict_path, GBS_path)
