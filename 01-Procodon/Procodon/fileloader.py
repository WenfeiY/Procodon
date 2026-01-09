"""
File loading, data processing, binary data processing
"""
import os

class Fileloader:
    
    """
    Load data from file or files, transform data into binary string
    """
    
    def __init__(
        self
    ):
        self.file_path = None
    
    @staticmethod
    def read_bin(
        in_file_path
    ) -> str:
        
        """
        Read data from file, transform into binary string
        
        Parameter:
            in_file_path: file path, str
        
        Return:
            binary_string: binary string for encoding, str, e.g. '01010111011...'
        """
        
        with open(in_file_path, 'rb') as source_file:
            # read file as binary format
            binary_data = source_file.read()
        
        # byte data to binary string
        binary_string = ''.join(format(byte, '08b') for byte in binary_data)
        
        return binary_string
    
    def read_bin_from_files(
        self,
        file_path_list: list
    ) -> dict:
        
        """
        Read data from files, transform into binary string
        
        Parameter:
            file_path_list: list containing all file path for encoding, list[str]
        
        Return:
            file_bin_dict: dictionary containing all binary string from files for encoding, str, e.g. '01010111011...'
        """

        file_bin_dict = {}
        
        for file_path in file_path_list:
            file_name = os.path.basename(file_path)
            binary_string = self.read_bin(file_path)
            file_bin_dict[file_name] = binary_string
        
        return file_bin_dict
    
    def read_bin_from_dir(
        self,
        file_dir: str
    ):
        
        """
        Read data from directory, transform into binary string
        
        Parameter:
            file_dir: file directory, str
        
        Return:
            file_bin_dict: dictionary containing all binary string from files for encoding, dict,
                {
                    filename1: '01010111011...',
                    filename2: '00101011101...',
                    ...
                }
        """
        
        file_path_list = []
        for file in os.listdir(file_dir):
            file_path_list.append(os.path.join(file_dir, file))
        
        file_bin_dict = self.read_bin_from_files(file_path_list)
        
        return file_bin_dict
    
    @staticmethod
    def text_to_bin(
        text: str,
        encoding: str = 'utf-8'
    ) -> str:
        
        """
        Read data from file, transform into binary string
        
        Parameter:
            text: text for encoding, str
            encoding: encoding format for computer, default UTF-8
        
        Return:
            binary_string: binary string for encoding, str, e.g. '01010111011...'
        """
        
        bytes_data = text.encode(encoding)
        binary_string = ''.join(format(byte, '08b') for byte in bytes_data)
        
        return binary_string

    def file_and_name_info_bin(
        self,
        file_bin_dict: dict,
        encoding: str = 'utf-8'
    ) -> str:
        
        """
        Concat all file names with their binary strings from data
        
        Parameter:
            file_bin_dict: dictionary containing all binary string from files for encoding, dict,
                {
                    filename1: '01010111011...',
                    filename2: '00101011101...',
                    ...
                }
            encoding: encoding format for computer, default UTF-8
        
        Return:
            conc_binary_string: concat binary string for encoding, str, e.g. '01010111011...'
        """
        
        conc_binary_string = ''
        
        for file_name in file_bin_dict:
            file_name_bin = self.text_to_bin(file_name, encoding)
            conc_binary_string += file_name_bin + file_bin_dict[file_name]
        
        return conc_binary_string
