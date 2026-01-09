"""
Decoder for Procodon
"""
from .utils import split_CDS, binary_array_to_text, write_bin_to_file
import json
import os

empty_bytes = bytearray(b'')

class Keydecoder(object):
    def __init__(
        self, 
        encoding_key: str,
        protein_name_abbr_path: str,
        protein_list_path: str
    ):
        """
        Read encoding key and decompress information, including codec, gene & index and storage mode
        
        Input:
            encoding_key: the storage key: str
            protein_name_abbr_path: file containing all protein name, abbreviation and gene type (simple represent no intron, intron represent gene with intron)
                e.g. file
                
        """        
        self.encoding_key = encoding_key
        self.protein_name_abbr_path = protein_name_abbr_path
        self.protein_list_path = protein_list_path
        
        self.codec_num = int(self.encoding_key.split(';')[0])
        self.protein_indexs = self.encoding_key.split(';')[1]
        self.storage_mode = self.encoding_key.split(';')[2]  #'S': single strain; 'M': mixed strains
        
        protein_abbr_dict = {}
        protein_type_dict = {}
        protein_abbr_list = []
        with open(protein_name_abbr_path, 'r') as f:
            for line in f:
                protein_abbr_dict[line.strip().split('\t')[1]] = line.strip().split('\t')[0]
                protein_type_dict[line.strip().split('\t')[0]] = line.strip().split('\t')[2]
                protein_abbr_list.append(line.strip().split('\t')[1])
        
        self.protein_abbr_dict = protein_abbr_dict
        self.protein_type_dict = protein_type_dict
        self.protein_abbr_list = protein_abbr_list
    
    #  Extract subset protein list with the start and end protein
    def continuous_abbr_indexs(
        self,
        start_abbr: str,
        end_abbr: str
    ) -> list:
    
        start_index = self.protein_abbr_list.index(start_abbr)
        end_index = self.protein_abbr_list.index(end_abbr)
    
        if start_index < end_index:
            return self.protein_abbr_list[start_index : end_index + 1]
        else:
            return self.protein_abbr_list[end_index : start_index + 1][::-1]
    
    #  Extract the abbreviation list of protein from key
    def decompress_abbr_indexs(
        self
    ):
        
        key_abbr_list = []
        blocks = self.protein_indexs.split('--')
        if len(blocks) == 1:
            key_abbr_list = blocks[0].split('-')
        else:
            for i in range(len(blocks) - 1):
                cur_block = blocks[i]
                if i == 0:
                    key_abbr_list += cur_block.split('-')[:-1]
                else:
                    key_abbr_list += cur_block.split('-')[1:-1]
                next_block = blocks[i + 1]
                cont_list = self.continuous_abbr_indexs(
                                cur_block.split('-')[-1], 
                                next_block.split('-')[0], 
                            )
                key_abbr_list += cont_list
            key_abbr_list += blocks[-1].split('-')[1:]
        self.key_abbr_list = key_abbr_list
    
    #  Change abbreviation to full name
    def abbr_to_full(
        self
    ):
        
        key_full_list = []
        for abbr in self.key_abbr_list:
            key_full_list.append(self.protein_abbr_dict[abbr])
        self.key_full_list = key_full_list
    
    #  Write the extracted protein names into files
    def write_protein_list(
        self
    ):
        with open(self.protein_list_path, 'w') as hd:
            for protein in self.key_full_list:
                hd.write(protein + '\t' + self.protein_type_dict[protein] + '\n')
    
    #  Decompress the information from key
    def decompress_key(
        self
    ):
        
        self.decompress_abbr_indexs()
        self.abbr_to_full()
        self.write_protein_list()
        
        print(self.codec_num)
        print(self.storage_mode)


class Dnadecoder:

    """
    Load codec
    Decode binary string from given CDSs
    """
    
    def __init__(
        self, 
        codec_path: str = 'Codec_112417.json', #  The codec used in this work
    ):
        """
        Convert CDSs into binary string

        Parameter:
            codec_path: json file path containing codec
        """
        #  Load codec
        with open(codec_path, 'r') as f:
            self.prot_encoding_dict = json.loads(f.read())
        self.codon_decode_dict = {}
        for aa in self.prot_encoding_dict:
            for num in self.prot_encoding_dict[aa]:
                for codon in self.prot_encoding_dict[aa][num]:
                    self.codon_decode_dict[codon] = num
    
    def decode_cds(
        self,
        cds_seq: str
    ) -> str:
        """
        Convert CDS into binary string by loaded codec

        Parameter:
            cds_seq: CDS in string
        
        Return:
            binary_string: decoded binary string
        """
        binary_string = ''
        #  Split CDS into list containing 3-bp codons
        codon_list = split_CDS(cds_seq)
        
        #  Decode CDS
        for i in range(0, len(codon_list)):
            if codon_list[i] in self.codon_decode_dict:
                binary_string += self.codon_decode_dict[codon_list[i]]
        
        return binary_string
    
    def decode_cds_list(
        self,
        cds_seq_list: list
    ) -> str:
        """
        Convert multiple CDSs into binary string by loaded codec

        Parameter:
            cds_seq_list: list containing CDSs in string
        
        Return:
            binary_string: decoded and joined binary string
        """
        complete_binary_string = ''
        
        #  Decode and join binaries
        for cds_seq in cds_seq_list:
            complete_binary_string += self.decode_cds(cds_seq)
        
        return complete_binary_string


class Pngdecoder():
    def __init__(
        self,
        out_dir: str = 'Decode_files'
    ):
        self.out_dir = out_dir
    
    @staticmethod
    def check_PNG_head_bin(
        binary_array,
        max_file_name_len: int = 100
    ):
        """
        Find PNG head in binary data and separate file name

        Parameters:
            binary_array: byte array object, converted from binary string
            max_file_name_len: maximum length of file name, preventing dubious data
        """
        #  Find standard binary on PNG head, with restrited length of file name
        position = binary_array[:max_file_name_len].find(b'\x89PNG\r\n\x1a\n')
        #  No PNG head was found
        if position == -1:
            return empty_bytes, empty_bytes
        
        # Extract file name before image data
        file_name_data = binary_array[:position]
        file_name = binary_array_to_text(file_name_data)
        
        return file_name, binary_array[position:]
    
    @staticmethod
    def check_PNG_tail_bin(
        binary_array,
    ):
        """
        Find PNG tail in binary data and separate PNG data

        Parameter:
            binary_array: byte array object, converted from binary string
        """
        position = binary_array.find(b'\x00\x00\x00\x00IEND\xaeB`\x82')
        
        if position == -1:
            return empty_bytes, empty_bytes
        
        return binary_array[:position + 12], binary_array[position + 12:]

    def decode_PNG(
        self,
        binary_string: str,
        max_file_name_len: int = 100
    ):
        """
        Extract PNG data from binary string and write into files

        Parameters:
            binary_string: given binary in string
            max_file_name_len: maximum length of file name, preventing dubious data
        """
        #  create directory if not exists
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        #  Convert binary string to binary array
        binary_array = bytearray(int(binary_string[i:i+8], 2) for i in range(0, len(binary_string), 8))
        self.binary_array = binary_array
        
        while binary_array:
            #  Extract filenames and PNG data
            file_name, binary_array = self.check_PNG_head_bin(binary_array, max_file_name_len)
            png_binary_array, binary_array = self.check_PNG_tail_bin(binary_array)
            #  Write data into files
            if png_binary_array:
                write_bin_to_file(png_binary_array, os.path.join(self.out_dir, file_name))
