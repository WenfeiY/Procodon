"""
Process codecs
"""
import json
from itertools import product
import pandas as pd
from Bio.SeqUtils import gc_fraction

def generate_decoding_dict(
    encoding_dict: dict
) -> dict:
    
    """
    Generate decoding dictionary from encoding dictionary
    
    Parameter:
        encoding_dict: dictionary contianing amino acid - codon - 0 or 1 encoding information
            e.g.
            {'F': {'0': ['TTT'], '1': ['TTC']},
             'L': {'0': ['CTA', 'CTC', 'CTT', 'CTG'], '1': ['TTG', 'TTA']},
             'I': {'0': ['ATA', 'ATC'], '1': ['ATT']},
             ...
             }
            
    Return:
        decoding_dict:dictionary contianing codon - 0 or 1 decoding information
            e.g.
            {'TTT': '0',
             'TTC': '1',
             'CTA': '0',
             'CTC': '0',
             'CTT': '0',
             ...
             }
    """
    
    decoding_dict = {}
    
    for aa in encoding_dict:
        for num in encoding_dict[aa]:
            for codon in encoding_dict[aa][num]:
                decoding_dict[codon] = num
    
    return decoding_dict

class Codonusageloader:
    def __init__(
        self,
        codon_table_path
    ):
        
        """
        Load amino acid - codon - relative frequency information from csv file
        File format: https://github.com/Edinburgh-Genome-Foundry/python_codon_tables/tree/master/python_codon_tables/codon_usage_data/tables
        
        Parameter:
            codon_table_path: file path of codon usage, str
        
        Returns:
            codon_frequency_dict: dictionary containing codon relative frequency, dict
                e.g.
                {'TAA': np.float64(0.47),
                 'TAG': np.float64(0.23),
                 'TGA': np.float64(0.3),
                 'GCA': np.float64(0.29),
                 'GCC': np.float64(0.22),
                 'GCG': np.float64(0.11),
                 'GCT': np.float64(0.38),
                 'TGC': np.float64(0.37),
                 'TGT': np.float64(0.63),
                 'GAC': np.float64(0.35),
                 ...
                 }
            codon_aa_dict: dictionary containing codon - amino acid information, dict
                e.g.
                {'TAA': '*',
                 'TAG': '*',
                 'TGA': '*',
                 'GCA': 'A',
                 'GCC': 'A',
                 'GCG': 'A',
                 'GCT': 'A',
                 'TGC': 'C',
                 'TGT': 'C',
                 'GAC': 'D',
                 'GAT': 'D',
                 ...
                 }
            aa_codon_usage_dict: dictionary containing amino acid - codon - relative frequency, dict
                e.g.
                {'A': [['GCA', np.float64(0.29)],
                  ['GCC', np.float64(0.22)],
                  ['GCG', np.float64(0.11)],
                  ['GCT', np.float64(0.38)]],
                 'C': [['TGC', np.float64(0.37)], ['TGT', np.float64(0.63)]],
                 'D': [['GAC', np.float64(0.35)], ['GAT', np.float64(0.65)]],
                 'E': [['GAA', np.float64(0.7)], ['GAG', np.float64(0.3)]],
                 'F': [['TTC', np.float64(0.41)], ['TTT', np.float64(0.59)]],
                 'G': [['GGA', np.float64(0.22)],
                  ['GGC', np.float64(0.19)],
                  ['GGG', np.float64(0.12)],
                  ['GGT', np.float64(0.47)]],
                 ...
                 }
        """
        
        codon_aa_dict = {}
        self.codon_table_path = codon_table_path
        codon_table = pd.read_csv(self.codon_table_path)
        
        codon_frequency_dict = {}
        aa_codon_usage_dict = {}
        
        for idx in codon_table.index:
            aa = codon_table.loc[idx, 'amino_acid']
            codon = codon_table.loc[idx, 'codon'].upper().replace('U', 'T')
            usage = codon_table.loc[idx, 'relative_frequency']
            codon_frequency_dict[codon] = usage
            codon_aa_dict[codon] = aa
            if not aa in aa_codon_usage_dict:
                aa_codon_usage_dict[aa] = [[codon, usage]]
            else:
                aa_codon_usage_dict[aa].append([codon, usage])
        
        self.codon_frequency_dict = codon_frequency_dict
        self.codon_aa_dict = codon_aa_dict
        self.aa_codon_usage_dict = aa_codon_usage_dict


class Codecloader:
    
    """
    Load codecs from file or genrate codec by codec number
    """
    
    def __init__(
        self
    ):
        
        self.codec_tem_path = None
        self.codec_num = None

    @staticmethod
    def generate_comb(
        codec_tem_dict: dict
    ):
        
        """
        Parse all encoding combinations from a codec template
        
        Parameter:
            codec_tem_path: json file containing template encoding dictionary, str
        
        Return:
            all combinations of the codec template encoding 0/1
        """
        
        comb_dict = {}
        for aa in codec_tem_dict:
            comb_dict[aa] = {'0': 0, '1': 0}        
        
        for para in product(*codec_tem_dict.values()):
            yield {k:v for k, v in zip(codec_tem_dict.keys(), para)}
    
    def load_codec_from_dict(
        self,
        encoding_dict: dict
    ):
        
        """
        Save encoding_dict and generate decoding_dict
        
        Parameter:
            encoding_dict: dictionary for protein encoding:
                e.g.
                {'F': {'0': ['TTT'], '1': ['TTC']},
                 'L': {'0': ['CTA', 'CTC', 'CTT', 'CTG'], '1': ['TTG', 'TTA']},
                 'I': {'0': ['ATA', 'ATC'], '1': ['ATT']},
                 ...
                 }
        
        Return:
            self.encoding_dict: dictionary for protein encoding:
                e.g.
                {'F': {'0': ['TTT'], '1': ['TTC']},
                 'L': {'0': ['CTA', 'CTC', 'CTT', 'CTG'], '1': ['TTG', 'TTA']},
                 'I': {'0': ['ATA', 'ATC'], '1': ['ATT']},
                 ...
                 }
            self.decoding_dict: dictionary contianing codon - 0 or 1 decoding information
                e.g.
                {'TTT': '0',
                 'TTC': '1',
                 'CTA': '0',
                 'CTC': '0',
                 'CTT': '0',
                 ...
                 }
        """
        
        self.encoding_dict = encoding_dict
        self.decoding_dict = generate_decoding_dict(self.encoding_dict)
    
    def load_codec_from_file(
        self,
        codec_path: str
    ):
        
        """
        Read encoding_dict from json file
        
        Parameter:
            codec_path: json file containing encoding dictionary
        
        Returns:
            self.encoding_dict: dictionary for protein encoding:
                e.g.
                {'F': {'0': ['TTT'], '1': ['TTC']},
                 'L': {'0': ['CTA', 'CTC', 'CTT', 'CTG'], '1': ['TTG', 'TTA']},
                 'I': {'0': ['ATA', 'ATC'], '1': ['ATT']},
                 ...
                 }
            self.decoding_dict: dictionary contianing codon - 0 or 1 decoding information
                e.g.
                {'TTT': '0',
                 'TTC': '1',
                 'CTA': '0',
                 'CTC': '0',
                 'CTT': '0',
                 ...
                 }
        """
        
        with open(codec_path, 'r') as f:
            encoding_dict = json.loads(f.read())
        
        self.encoding_dict = encoding_dict
        self.decoding_dict = generate_decoding_dict(self.encoding_dict)
    
    def generate_codec(
        self,
        codec_tem_path: str,
        codec_num: str
    ):
        
        """
        Read codec template from json file, generate combinations and get codec by codec number
        
        Parameter:
            codec_tem_path: json file path containing codec template, str
            codec_num: codec ID, str
        
        Return:
            self.encoding_dict: dictionary for protein encoding:
                e.g.
                {'F': {'0': ['TTT'], '1': ['TTC']},
                 'L': {'0': ['CTA', 'CTC', 'CTT', 'CTG'], '1': ['TTG', 'TTA']},
                 'I': {'0': ['ATA', 'ATC'], '1': ['ATT']},
                 ...
                 }
            self.decoding_dict: dictionary contianing codon - 0 or 1 decoding information
                e.g.
                {'TTT': '0',
                 'TTC': '1',
                 'CTA': '0',
                 'CTC': '0',
                 'CTT': '0',
                 ...
                 }
        """
        
        self.codec_tem_path = codec_tem_path
        self.codec_num = codec_num
        
        with open(codec_tem_path, 'r') as f:
            codec_tem_dict = json.loads(f.read())
        
        self.codec_combs = list(self.generate_comb(codec_tem_dict))
        
        encoding_dict = {}
        for aa in self.codec_combs[codec_num]:
            if self.codec_combs[codec_num][aa] == '0':
                encoding_dict[aa] = codec_tem_dict[aa]
            else:
                encoding_dict[aa] = {'0': codec_tem_dict[aa]['1'], '1': codec_tem_dict[aa]['0']}
        
        self.encoding_dict = encoding_dict
        
        self.decoding_dict = generate_decoding_dict(self.encoding_dict)

class Codecgenerator:
    
    """
    Load codec template, generate combinations, calculate GC bias score (GBS), filter codecs and write codecs into file
    """
    
    def __init__(
        self,
        codec_tem_path: str,
        codon_table_path: str,
        aa_frequency_path: str
    ):
        
        """
        Load codec template, load codon & amino acid frequency and calculate sum of GSC
        
        Parameter:
            codec_tem_path: json file containing template encoding dictionary, str
            codon_table_path: file path of codon usage, str
                File format: https://github.com/Edinburgh-Genome-Foundry/python_codon_tables/tree/master/python_codon_tables/codon_usage_data/tables
            aa_frequency_path: json file containing amino acid frequency dictionary, str
                (generate by custom python script)
        
        Return:
            self.codec_tem_path: codec_tem_path
            self.codon_table_path: codon_table_path
            self.codec_tem_dict: template dictionary for protein encoding:
                e.g.
                {'F': {'0': ['TTT'], '1': ['TTC']},
                 'L': {'0': ['CTA', 'CTC', 'CTT', 'CTG'], '1': ['TTG', 'TTA']},
                 'I': {'0': ['ATA', 'ATC'], '1': ['ATT']},
                 ...
                 }
            self.aa_frequency_dict: dictionary containing amino acid - frequency (aa number count in all protein sequences)
                e.g.
                {
                 "F": 129516,
                 "L": 277988,
                 "I": 191677,
                 ...
                 }
            self.codon_frequency_dict: dictionary containing codon relative frequency, dict
                e.g.
                {'TAA': np.float64(0.47),
                 'TAG': np.float64(0.23),
                 'TGA': np.float64(0.3),
                 'GCA': np.float64(0.29),
                 'GCC': np.float64(0.22),
                 'GCG': np.float64(0.11),
                 'GCT': np.float64(0.38),
                 'TGC': np.float64(0.37),
                 'TGT': np.float64(0.63),
                 'GAC': np.float64(0.35),
                 ...
                 }
            self.aa_GC_score_dict: dictionary containing GSC for each amino acid encoding group, dict
                e.g.
                {'F': {'0': np.float64(0.0), '1': np.float64(43172.0)},
                 'L': {'0': np.float64(47144.2),
                  '1': np.float64(128464.2)},
                 'I': {'0': np.float64(0.0), '1': np.float64(31343.4)},
                 'V': {'0': np.float64(71973.8), '1': np.float64(81321.0)},
                 ....
                 }
            self.aa_GC_score_sum: sum of all GSC for each amino acid encoding group, float
        """
        
        self.codec_tem_path = codec_tem_path
        self.codon_table_path = codon_table_path
        
        with open(codec_tem_path, 'r') as f:
            self.codec_tem_dict = json.loads(f.read())
        
        with open(aa_frequency_path, 'r') as f:
            self.aa_frequency_dict = json.loads(f.read())
        
        codon_table = pd.read_csv(self.codon_table_path)
        
        codon_frequency_dict = {}
        for idx in codon_table.index:
            codon = codon_table.loc[idx, 'codon'].upper().replace('U', 'T')
            usage = codon_table.loc[idx, 'relative_frequency']
            codon_frequency_dict[codon] = usage
        
        self.codon_frequency_dict = codon_frequency_dict
        
        aa_GC_score_dict = {}
        for aa in self.codec_tem_dict:
            aa_GC_score_dict[aa] = {'0': 0, '1': 0}
            aa_frequency = self.aa_frequency_dict[aa]
            #  Calculate GS-0 and GS-1 of each AA
            for num in ['0', '1']:
                bin_GC_score = 0
                codons = self.codec_tem_dict[aa][num]
                codon_frequencies = []
                for codon in codons:
                    codon_frequencies.append(self.codon_frequency_dict[codon])
                for codon in codons:
                    bin_GC_score += aa_frequency * self.codon_frequency_dict[codon] / sum(codon_frequencies) * gc_fraction(codon)
                aa_GC_score_dict[aa][num] = bin_GC_score
        
        #  Calculate sum of GC Score
        aa_GC_score_sum = 0
        for aa in aa_GC_score_dict:
            for num in ['0', '1']:
                aa_GC_score_sum += aa_GC_score_dict[aa][num]
        
        self.aa_GC_score_dict = aa_GC_score_dict
        self.aa_GC_score_sum = aa_GC_score_sum
    
    @staticmethod
    def generate_comb(
        codec_tem_dict: dict
    ):
        
        """
        Parse all encoding combinations from a codec template
        
        Parameter:
            codec_tem_path: json file containing template encoding dictionary, str
        
        Return:
            all combinations of the codec template encoding 0/1
        """
        
        comb_dict = {}
        for aa in codec_tem_dict:
            comb_dict[aa] = {'0': 0, '1': 0}        
        
        for para in product(*codec_tem_dict.values()):
            yield {k:v for k, v in zip(codec_tem_dict.keys(), para)}
    
    def filter_GBS(
        self,
        cut_off: float
    ):
        
        """
        Filter codec combinations by percentage difference between GS encoding 0 and 1 with given cut-off
        
        Parameter:
            cut_off: cut-off percentage difference between GBS encoding 0 and 1, float
        
        Return:
            combination_pass: filtered encoding combinations, list
            prop_pass: codec ID and GS-0 info, list
        """
        
        codec_combs = list(self.generate_comb(self.codec_tem_dict))
        
        combination_pass = []
        prop_pass = []
        
        for i in range(len(codec_combs)):
            combination = codec_combs[i]
            aa_GC_score = 0
            for aa in combination:
                aa_GC_score += self.aa_GC_score_dict[aa][combination[aa]]

            #  Filter codec by GBS
            prop = aa_GC_score / self.aa_GC_score_sum
            if 0.5 - (cut_off/100/2) <= prop <= 0.5 + (cut_off/100/2):
                combination_pass.append([i, combination])
                prop_pass.append([i, abs(prop)])
        return combination_pass, prop_pass
    
    def generate_codec(
        self,
        combination
    ):
        
        """
        Generate encoding dict from encoding combination
        
        Parameter:
            combination: dictionary containing group combination, dict
                e.g.
                {
                 "F": "0",
                 "L": "0",
                 "I": "0",
                 ...
                 }
        
        Return:
            generated_dict: dictionary contianing generated encoding dictionary, dict
                e.g.
                {'F': {'0': ['TTT'], '1': ['TTC']},
                 'L': {'0': ['CTA', 'CTC', 'CTT', 'CTG'], '1': ['TTG', 'TTA']},
                 'I': {'0': ['ATA', 'ATC'], '1': ['ATT']},
                 ...
                 }
        """
        
        generated_dict = {}
        for aa in combination:
            if combination[aa] == '0':
                generated_dict[aa] = self.codec_tem_dict[aa]
            else:
                generated_dict[aa] = {'0': self.codec_tem_dict[aa]['1'], '1': self.codec_tem_dict[aa]['0']}
        
        return generated_dict
    
    def generate_codecs(
        self,
        combinations
    ):
        
        """
        Generate encoding dicts from encoding combinations
        
        Parameter:
            combinations: dictionary containing combinations, dict
                e.g.
                {
                "0": {
                    "F": "0",
                    "L": "0",
                    "I": "0",
                    ...
                    }
                "1": {
                    "F": "0",
                    "L": "0",
                    "I": "0",
                    ...
                    }
                }
        
        Return:
            generated_dicts: dictionary contianing generated encoding dictionarys, dict
                e.g.
                {'0':{
                    'F': {'0': ['TTT'], '1': ['TTC']},
                    'L': {'0': ['CTA', 'CTC', 'CTT', 'CTG'], '1': ['TTG', 'TTA']},
                    'I': {'0': ['ATA', 'ATC'], '1': ['ATT']},
                    ...
                    },
                ...
                }
        """
        
        generated_dicts = {}
        
        for combination_info in combinations:
            index, combination = combination_info
            generated_dict = {}
            for aa in combination:
                if combination[aa] == '0':
                    generated_dict[aa] = self.codec_tem_dict[aa]
                else:
                    generated_dict[aa] = {'0': self.codec_tem_dict[aa]['1'], '1': self.codec_tem_dict[aa]['0']}
            
            generated_dicts[index] = generated_dict
            
        return generated_dicts
    
    def get_filtered_codecs(
        self,
        cut_off: float
    ):
        
        """
        Filter codec combinations and generate encoding_dicts
        
        Parameter:
            cut_off: cut-off percentage difference between GSC encoding 0 and 1, float
        
        Return:
            encoding_dict_pass: dictionary containing filtered encoding_dicts, dict
                e.g.
                {'0':{
                    'F': {'0': ['TTT'], '1': ['TTC']},
                    'L': {'0': ['CTA', 'CTC', 'CTT', 'CTG'], '1': ['TTG', 'TTA']},
                    'I': {'0': ['ATA', 'ATC'], '1': ['ATT']},
                    ...
                    },
                ...
                }
            prop_pass: codec ID and GS-0 info, list
        """
        
        combination_pass, prop_pass = self.filter_GBS(cut_off)
        encoding_dict_pass = self.generate_codecs(combination_pass)
        
        return encoding_dict_pass, prop_pass
    
    def get_and_write_filtered_codecs(
        self,
        codec_out_path: str,
        cut_off = 1
    ):
        
        """
        Filter codec combinations, generate encoding_dicts and write into json file
        
        Parameter:
            codec_out_path: Return json file path of filtered codecs, str 
            cut_off: cut-off percentage difference between GSC encoding 0 and 1, float
        Return:
            prop_pass: codec ID and GS-0 info, list
        """
        
        combination_pass, prop_pass = self.filter_GBS(cut_off)
        encoding_dict_pass = self.generate_codecs(combination_pass)
        
        self.combination_pass = combination_pass
        self.encoding_dict_pass = encoding_dict_pass
        
        with open(codec_out_path, 'w') as hd:
            json.dump(encoding_dict_pass, hd)
        return prop_pass
"""
