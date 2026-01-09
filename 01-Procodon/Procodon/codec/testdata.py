# -*- coding: utf-8 -*-
"""
Generate test data for finding available codecs
"""

import json
import random
from ..utils import count_available_space



class Testdata:
    
    '''
    Generate test binary data for each CDS/protein with different 0/1 proportion and write into json files
    default:
        '1' proportion: [0, 0.2, 0.4, 0.6, 0.8, 1]
        for each gene and proportion, randomly generate 20 binary strings
    '''
    
    def __init__(
        self,
        test_seq_dict: dict
    ):
        
        '''
        load CDS/protein seqs for test
        
        Input:
            test_seq_dict: dictionary with gene as key and CDS/protein seq as value
            e.g.:
                {
                    'gene1':
                        {
                            'CDS': '...',
                            'protein': '...'
                        },
                    'gene2':
                        {
                            'CDS': '...',
                            'protein': '...'
                        },
                }
        '''
        
        self.test_seq_dict = test_seq_dict
    
    @staticmethod
    def generate_random_bin(
        bin_len: int, 
        bin_count: int, 
        prop_1: float
    ) -> list:
        
        '''
        Generate test binary string for by length, number and 0/1 proportion
        
        Input:
            bin_len: length of binary string, int
            bin_count: number of binary string to be generated, int
            prop_1: '1' proportion in binary string, float
            
        Return:
            bin_list: list containing generated binary strings
        '''
        
        count_1 = int(bin_len * prop_1)
        count_0 = bin_len - count_1
        
        bin_list = []
        
        for i in range(bin_count):
            binary_string_list = ['0'] * count_0 + ['1'] * count_1
            random.shuffle(binary_string_list)
            bin_list.append(''.join(binary_string_list))
        
        return bin_list
    
    @staticmethod
    def write(
        gene_bin_dict: dict,
        gene_bin_dict_out: str
    ):
        '''
        Write binary test data into json file
        
        Input:
            gene_bin_dict: dictionary with gene as key and binary string seq as value
            e.g.:
                {
                    'gene1':
                        {
                            0: ['0101010....',
                                '1010101...,',
                                ...
                                ],
                            0.1: ['0101010....',
                                  '1010101...,',
                                  ...
                                  ],
                            ...
                        },
                    'gene2':
                        {
                            0: ['0101010....',
                                '1010101...,',
                                ...
                                ],
                            0.1: ['0101010....',
                                  '1010101...,',
                                  ...
                                  ],
                            ...
                        },
                    ...
                }
        '''
        
        with open(gene_bin_dict_out, 'w') as hd:
            json.dump(gene_bin_dict, hd)
    
    def prep_test_data(
        self,
        bin_count: int,
        prop_1_list: list[float],
        gene_bin_dict_out: str
    ):
        
        '''
        Generate test binary data for each gene by available coding length, number and 0/1 proportion.
        Write test binary data into json file.
        
        Input:
            bin_count: number of binary string to be generated, int
            prop_1_list: List containing '1' proportion in binary string, list[float]
            gene_bin_dict_out: output path for generated test binary data
        '''
        
        gene_bin_dict = {}
        
        i = 1
        for gene in self.test_seq_dict:
            gene_bin_dict[gene] = {}
            prot_seq = self.test_seq_dict[gene]['Protein']
            encode_space = count_available_space(prot_seq)
            for prop_1 in prop_1_list:
                if prop_1 == 0 or prop_1 == 1:
                    gene_bin_dict[gene][prop_1] = self.generate_random_bin(encode_space, 1, prop_1)
                else:
                    gene_bin_dict[gene][prop_1] = self.generate_random_bin(encode_space, bin_count, prop_1)
            if i % 100 == 0:
                print(str(i) + ' genes done')
            i += 1
        
        print('Writing into ' + gene_bin_dict_out + ' ...')
        self.write(gene_bin_dict, gene_bin_dict_out)
