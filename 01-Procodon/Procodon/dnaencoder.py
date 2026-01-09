"""
DNA encoder module for protein storage
"""
import random
from .utils import split_CDS, get_empty_codon_dict, get_empty_select_dict, count_available_space
from .codec.codec import Codonusageloader

class Dnaencoder:
    
    """
    load codec and codon usage information
    generate recoded CDS sequences with WT CDS, binary string and codec
    The start codons of some prokaryotic CDS are not ATG. 
    To avoid systematic error, while species_type == 'pr', the start codon will not be checked.
    """
    
    def __init__(
        self, 
        encoding_dict: dict,
        codon_table_path: str,
        species_type = 'eu'
    ):
        
        """
        Parameters:
            encoding_dict: dictionary for protein encoding:
                {
                    "F": {
                        "0": [
                            "TTT"
                        ],
                        "1": [
                            "TTC"
                        ]
                    },
                    "L": {
                        "0": [
                            "CTA",
                            "CTC",
                            "CTT",
                            "CTG"
                        ],
                        "1": [
                            "TTG",
                            "TTA"
                        ]
                    },
                ...
                }
            codon_table_path: file path of codon usage, download from https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/blob/master/codon_usage_data/tables/ 
            species_type: 'eu' (eukaryotic) or 'pr' (prokaryotic), str
        
        Output:
            codon_frequency_dict: {codon: relative_frequency, ...}
            codon_aa_dict: {codon: relevant_amino_acid, ...}
            aa_codon_usage_dict: {aa1: [[codon1, relative_frequency_1], [codon2, relative_frequency_2]], ...}
        """
        
        self.encoding_dict = encoding_dict
        culoader = Codonusageloader(codon_table_path)
        
        self.codon_frequency_dict = culoader.codon_frequency_dict
        self.codon_aa_dict = culoader.codon_aa_dict
        self.aa_codon_usage_dict = culoader.aa_codon_usage_dict
        
        species_types = ['eu', 'pr']
        if species_type not in species_types:
            raise ValueError(
                f'Invalid species type {species_type}. Should be one of: '
                f'{", ".join(species_types)}.'
            )
        
    def get_empty_codon_usage_dict(
        self
    ) -> dict:
        
        """
        Generate empty codon usage dictionary
        
        Return: 
            empty_encode_codon_dict:
                dictionary structure:
                    {amino_acid_1: {
                        '0': [[codon_1, codon_2], [0, 0]],
                        '1': [[codon_3, codon_4], [0, 0]]
                        },
                    amino_acid_2: {
                        '0': [[codon_5], [0]],
                        '1': [[codon_6], [0]]
                        }
                    ...
                    }
        """
        
        empty_encode_codon_dict = {}
        
        for aa in self.encoding_dict:
            empty_encode_codon_dict[aa] = {
                '0': [self.encoding_dict[aa]['0'],[0] * len(self.encoding_dict[aa]['0'])],
                '1': [self.encoding_dict[aa]['1'], [0] * len(self.encoding_dict[aa]['1'])]
            }
        return empty_encode_codon_dict
    
    def codon2aa(
        self,
        codon: str
    ) -> str:
        
        """
        Convert codon to 1-latter amino acid
        
        Parameter:
            codon: 3-latter codons, str
        
        Return:
            1-latter amino acid, str
        """
        
        return self.codon_aa_dict[codon]

    def codon_list_to_aa_list(
        self,
        codon_list: list
    ) -> list:
        
        """
        Convert codon list to 1-latter amino acids list
        
        Parameter:
            codon_list: list containing 3-latter codons, list[str]
        
        Return:
            aa_list: list containing translated 1-latter amino acid, list[str]
        """
        
        aa_list = []
        for codon in codon_list:
            aa_list.append(self.codon2aa(codon))
        return aa_list
    
    def cds_to_protein_seq(
        self,
        cds_seq: str
    ) -> str:
        
        """
        Convert CDS to protein sequence
        
        Parameter:
            cds_seq: CDS, str
        
        Return:
            protein sequence, str
        """
        
        codon_list = split_CDS(cds_seq)
        aa_list = self.codon_list_to_aa_list(codon_list)
        return ''.join(aa_list)
    
    def cds_match_protein(
        self,
        protein_seq: str,
        cds_seq: str,
        codon_aa_dict: dict
    ):  
        
        """
        Check whether the CDS can be translated to the given protein sequence with codon_aa_dict
        The start codons of some prokaryotic CDS are not ATG. To avoid systematic error, while species_type == 'pr', the start codon will not be checked
        
        Parameters:
            protein_seq: protein sequence, str
            cds_seq: CDS sequence, str
            codon_aa_dict: dictionary containing AAs with relevant codons as keys {codon: relevant_amino_acid, ...}
        
        Output:
            protein sequence, str
        """
        
        protein_seq = protein_seq.upper()
        cds_seq = cds_seq.upper()
        if not len(protein_seq) == len(cds_seq)/3:
            if len(protein_seq) == len(cds_seq)/3 - 1:
                protein_seq += '*'
            else:
                raise ValueError(
                    'CDS and protein length not match!'
                )
        translated_prot_seq = self.cds_to_protein_seq(self, cds_seq)
        if translated_prot_seq == protein_seq:
            return True
        elif self.species_type == 'pr':
            if 'M' + translated_prot_seq[1:] == protein_seq:
                return True
        raise ValueError(
            'CDS and protein sequence not match!'
        )
    
    @staticmethod
    def random_choose_codon(
        codon_select_list: list
    ):
        
        """
        Randomly choose codon from a codon list, and return the selected codon with the rest codon list
        
        Parameter:
            codon_select_list: list containing determined number of each codon encoding 1 aa as 0 or 1, list, e.g. ['ACT', 'ACT', 'ACT', 'ACG']
        
        Return:
            selected_codon: selected codon for recoding, str
            rest_codon_list: codon_select_list minus selected codon, list[str]
        """
        
        codon_num = random.choice(list(range(len(codon_select_list))))
        selected_codon = codon_select_list[codon_num]
        return selected_codon, codon_select_list[:codon_num] + codon_select_list[codon_num + 1:]
    
    def count_codon_frequency(
        self,
        codon_list: list
    ):
        
        """
        Count each codon number in the given codon_list
        
        Parameter:
            codon_list: list containing 3-latter codons, list[str]
        
        Return:
            codon_usage_dict: dictionary containing amino acid, codon and codon number used in the given codon_list 
                dictionary structure:
                    {amino_acid_1: {
                        '0': [[codon_1, codon_2], [1, 5]],
                        '1': [[codon_3, codon_4], [3, 6]]
                        },
                    amino_acid_2: {
                        '0': [[codon_5], [5]],
                        '1': [[codon_6], [9]]
                        }
                    ...
                    }
        """
        
        codon_dict = get_empty_codon_dict()
        codon_usage_dict = self.get_empty_codon_usage_dict()
        
        for codon in codon_list:
            codon_dict[codon] += 1
        
        for aa in codon_usage_dict:
            for num in codon_usage_dict[aa]:
                for i in range(len(codon_usage_dict[aa][num][0])):
                    codon = codon_usage_dict[aa][num][0][i]
                    codon_usage_dict[aa][num][1][i] = codon_dict[codon]
        return codon_usage_dict
    
    @staticmethod
    def count_aa_bin_num(
        aa_list: list,
        binary_string: str
    ):
        
        """
        Count binary code number for each amino acid
        
        Parameters:
            aa_list: list containing 1-latter amino acid of encoding sequence, list[str]
            binary_string: binary string for encoding, str, e.g. '01010111011...'
        
        Return:
            aa_bin_count_dict: dictionary containing binary code number of each amino acid
                e.g.
                {
                    'F': {'0': 4, '1': 3},
                    'L': {'0': 6, '1': 9},
                    'I': {'0': 1, '1': 2},
                    ...
                }
        """
        
        aa_bin_count_dict = {
            'F': {'0': 0, '1': 0},
            'L': {'0': 0, '1': 0},
            'I': {'0': 0, '1': 0},
            'V': {'0': 0, '1': 0},
            'S': {'0': 0, '1': 0},
            'P': {'0': 0, '1': 0},
            'T': {'0': 0, '1': 0},
            'A': {'0': 0, '1': 0},
            'Y': {'0': 0, '1': 0},
            'H': {'0': 0, '1': 0},
            'Q': {'0': 0, '1': 0},
            'N': {'0': 0, '1': 0},
            'K': {'0': 0, '1': 0},
            'D': {'0': 0, '1': 0},
            'E': {'0': 0, '1': 0},
            'C': {'0': 0, '1': 0},
            'R': {'0': 0, '1': 0},
            'G': {'0': 0, '1': 0}
        }
        j = 0
        for i in range(len(aa_list)):
            if aa_list[i] in aa_bin_count_dict:
                if binary_string[j] == '0':
                    aa_bin_count_dict[aa_list[i]]['0'] += 1
                else:
                    aa_bin_count_dict[aa_list[i]]['1'] += 1
                j += 1
            if j == len(binary_string):
                break
        return aa_bin_count_dict

    def build_select_list(
        self,
        codon_sub_list: list, 
        codon_frequency_list: list, 
        aa_bin_count: int
    ):
        
        """
        Generate codon select list for current encoding sequence
        
        Parameters:
            codon_sub_list: list containing codons for 1 amino acid, 1 binary code, list[str], e.g. ["ATA", "ATC"] (encoding 0)
            codon_frequency_list: list containing relevant codons frequency of codon_sub_list, list[int], e.g. [3, 5]
            aa_bin_count: amino acid number for encoding the designate binary code, int
        
        Return:
            codon_select_list: list containing determined number of each codon encoding 1 aa as 0 or 1, list, e.g. ['ACT', 'ACT', 'ACT', 'ACG']
        """
        
        if len(codon_sub_list) == 1:
            return codon_sub_list * aa_bin_count
        codon_select_list = []
        sum_frequency = sum(codon_frequency_list)
        if sum_frequency == 0:
            codon_frequency_list = []
            for curr_codon in codon_sub_list:
                for aa in self.aa_codon_usage_dict:
                    for codon_term in self.aa_codon_usage_dict[aa]:
                        if codon_term[0] == curr_codon:
                            codon_frequency_list.append(codon_term[1])
            sum_frequency = sum(codon_frequency_list)
        for i in range(len(codon_sub_list) - 1):
            codon_select_list += round(codon_frequency_list[i] / sum_frequency * aa_bin_count) * [codon_sub_list[i]]
        codon_select_list += (aa_bin_count - len(codon_select_list)) * [codon_sub_list[-1]]
        return codon_select_list

    def build_select_dict(
        self,
        codon_usage_dict: dict,
        aa_list: list,
        binary_string: str
    ):
        
        """
        Generate every codon select list for current encoding sequence, stored in dictionary
        
        Parameters:
            codon_usage_dict: dictionary containing amino acid, codon and codon number used in the given codon_list 
                dictionary structure:
                    {amino_acid_1: {
                        '0': [[codon_1, codon_2], [1, 5]],
                        '1': [[codon_3, codon_4], [3, 6]]
                        },
                    amino_acid_2: {
                        '0': [[codon_5], [5]],
                        '1': [[codon_6], [9]]
                        }
                    ...
                    }
            aa_list: list containing translated 1-latter amino acid, list[str]
            binary_string: binary string for encoding, str, e.g. '01010111011...'
        
        Return:
            select_dict: dictionary containing all codon_select_list for encoding
                e.g.
                {
                    'F': {'0': ['TTT', 'TTT', 'TTT'], '1': ['TTC']},
                    'L': {'0': ['CTA', 'CTA', 'CTT'], '1': ['TTA', 'TTA']},
                    'I': {'0': ['ATA'], '1': ['ATT', 'ATT', 'ATT']},
                    ...
                }
        """

        select_dict = get_empty_select_dict()
        aa_bin_count_dict = self.count_aa_bin_num(aa_list, binary_string)
        for aa in select_dict:
            for num in ['0', '1']:
                select_dict[aa][num] = self.build_select_list(codon_usage_dict[aa][num][0], codon_usage_dict[aa][num][1], aa_bin_count_dict[aa][num])
        return select_dict
    
    def cal_encode_edge(
        self,
        prot_seq: str,
        binary_string: str
    ) -> int:
        
        """
        Determine the position of protein used for encoding and return the rest binary string if not used
        
        Parameters:
            prot_seq: protein sequence, str
            binary_string: binary string for encoding, str, e.g. '01010111011...'
        
        Return:
            encode_edge: the end index of amino acid used + 1, int
            rest_bin: unused binary_string, str
        """
        
        encode_space = count_available_space(prot_seq)
        count = 0
        for i in range(len(prot_seq)):
            if not prot_seq[i] in ['M', 'W', '*']:
                count += 1
                if count == len(binary_string):
                    break
        if encode_space >= len(binary_string):
            rest_bin = ''
        else:
            rest_bin = binary_string[count:]
        return i + 1, rest_bin
    
    def binary_to_DNA(
        self,
        cds_seq: str,
        prot_seq: str,
        binary_string: str
    ):
        
        """
        Encode given binary string into protein, return recoded CDS, encode_edge, used binary string and unused binary_string
        
        Parameters:
            cds_seq: CDS sequence, str
            prot_seq: protein sequence, str
            binary_string: binary string for encoding, str, e.g. '01010111011...'
        Return:
            recoded_seq: recoded CDS sequence, str
            encode_edge: the end index of amino acid used + 1, int
            encoded_bin: used binary string, str
            rest_bin: unused binary_string
        """
        
        prot_seq = prot_seq.upper()
        codon_list = split_CDS(cds_seq)
        aa_list = [aa for aa in prot_seq]
        encode_edge, rest_bin = self.cal_encode_edge(prot_seq, binary_string)
        codon_usage_dict = self.count_codon_frequency(codon_list[:encode_edge])
        select_dict = self.build_select_dict(codon_usage_dict, aa_list, binary_string)
        recoded_seq = ''
        j = 0
        for i in range(len(aa_list)):
            if aa_list[i] in self.encoding_dict:
                selected_codon, select_dict[aa_list[i]][binary_string[j]] = self.random_choose_codon(select_dict[aa_list[i]][binary_string[j]])
                recoded_seq += selected_codon
                j += 1
            else:
                recoded_seq += codon_list[i]
            if j == len(binary_string):
                break
        encoded_bin = binary_string[:j]
        for x in range(encode_edge, len(aa_list)):
            recoded_seq += codon_list[x]
        return recoded_seq, encode_edge, encoded_bin, rest_bin
