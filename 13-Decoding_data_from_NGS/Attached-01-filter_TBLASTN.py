"""
Calculate threshold of bit-score and filter TBLASTN result
"""
import os
import sys
from Bio import SeqIO

def cal_cut_off(
    prot_fa_path: str
):
    """
    Calculate threshold by protein length. In this work, the reads length of NGS is set to 150 bp.

    Parameter:
        prot_fa_path: file path of protein sequence in FASTA format, str
    
    Return:
        calculated threshold, int
    """
    record = SeqIO.read(prot_fa_path, 'fasta')
    prot_len = len(record.seq)

    return int(min([prot_len, 50]) * 1.3)

def filter_BLAST_out(
    in_file_path: str, 
    out_file_path: str, 
    cut_off: int
):
    """
    Run filter program in Linux system

    Parameters:
        in_file_path: file path of input TBLASTN result
        out_file_path: output file path of filtered TBLASTN result
        cut_off: calculated threshold, int
    """
    os.system('cat ' + in_file_path + " | awk \'$12 > " + str(cut_off) + "\' > " + out_file_path)


if __name__ == "__main__":
    #  Get parameters of file paths
    prot_fa_path = sys.argv[1]
    in_file_path = sys.argv[2]
    out_file_path = sys.argv[3]
    #  Calculate threshold of bit-score
    cut_off = cal_cut_off(prot_fa_path)
    #  Filter TBLASTN
    filter_BLAST_out(in_file_path, out_file_path, cut_off)
