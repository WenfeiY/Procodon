"""
Encoding images into yeast genome with bio-index from Nakaseomyces glabratus
"""

from Bio import SeqIO
from Procodon.seqloader import Seqloader
from Procodon.codec.codec import Codecloader
from Procodon.dnaencoder import Dnaencoder
from Procodon.fileloader import Fileloader
import json
from Bio.Seq import Seq

def replace_CDS(
    record_seq: str, 
    replace_cds_seq: str, 
    gene_info: dict
):
    """
    Replace CDS with recoded sequences

    Parameters:
        record_seq: sequence of SeqRecord object in string; in this condition, chromosome sequence
        replace_cds_seq: sequence for replacing the WT in string; in this condition, recoded sequence
        gene_info: dictionary containing gene names, directions, CDS positions and introns

    Returns
    record_seq: sequences with replaced parts; in this condition, chromosome sequence with CDS replaced by recoded sequences

    """
    
    if gene_info['with_intron']:
        cds_part_start = 0
        for pos in gene_info['CDS_position']:
            cds_part = Seq(replace_cds_seq[cds_part_start: cds_part_start + pos[1] - pos[0]])
            if gene_info['direction'] == -1:
                cds_part = cds_part.reverse_complement()
            record_seq = record_seq[:pos[0]] + cds_part + record_seq[pos[1]:]
            cds_part_start = cds_part_start + pos[1] - pos[0]
    else:
        cds_part = Seq(replace_cds_seq)
        if gene_info['direction'] == -1:
            cds_part = cds_part.reverse_complement()
        record_seq = record_seq[:gene_info['CDS_position'][0]] + cds_part + record_seq[gene_info['CDS_position'][1]:]
    return record_seq

#  Five images extracted from "The arrival of a train"
file_path_list = ['./Images/frame_00240_comp_2.bmp', 
                  './Images/frame_00255_comp_2.bmp', 
                  './Images/frame_00270_comp_2.bmp', 
                  './Images/frame_00285_comp_2.bmp', 
                  './Images/frame_00300_comp_2.bmp'
]

#  Load CDSs and proteins
seqload= Seqloader(seq_dir = 'S.cerevisiae')  #Same files in directory "Codec_test", containing CDSs and proteins

#  Load codec
codecload=Codecloader()
codecload.load_codec_from_file('codec_112417.json')

#  Initialize Dnaencoder
encode = Dnaencoder(codecload.encoding_dict, 's_cerevisiae_4932.csv')

#  Load gene order from Nakaseomyces glabratus
with open('Ng.bio-index.txt', 'r') as f:
    gene_list = f.read().splitlines()

#  Load files and convert to binaries
fileload = Fileloader()
file_bin_dict = fileload.read_bin_from_files(file_path_list)

#  Add filenames into binary data
full_binary = ''
for file_name in file_bin_dict:
    full_binary += ''.join(format(byte, '08b') for byte in file_name.encode('utf-8'))
    full_binary += file_bin_dict[file_name]

#  Recode sequences
rest_bin = full_binary
recode_dict = {}
for gene in gene_list:
    recoded_seq, _, _, rest_bin = encode.binary_to_DNA(seqload.gene_seq_dict[gene]['CDS'], seqload.gene_seq_dict[gene]['Protein'], rest_bin)
    recode_dict[gene] = recoded_seq
    if rest_bin == '':
        break

#  Write recoded sequences into json file
with open('S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.json', 'w') as hd:
    json.dump(recode_dict, hd)

#  Dictionary containing yeast chromosome ID and name
chr_dict = {
    'BK006935.2': 'YA',
    'BK006936.2': 'YB',
    'BK006937.2': 'YC',
    'BK006938.2': 'YD',
    'BK006939.2': 'YE',
    'BK006940.2': 'YF',
    'BK006941.2': 'YG',
    'BK006934.2': 'YH',
    'BK006942.2': 'YI',
    'BK006943.2': 'YJ',
    'BK006944.2': 'YK',
    'BK006945.2': 'YL',
    'BK006946.2': 'YM',
    'BK006947.3': 'YN',
    'BK006948.2': 'YO',
    'BK006949.2': 'YP'  
}

#  Load WT chromosome sequences from Genbank file
gb_dict = {}
for record in SeqIO.parse('Saccharomyces.cerevisiae.gbk', 'genbank'):
    gb_dict[record.id] = record

#  Load yeast gene information containing gene names, directions, CDS positions and introns
with open('sc_gene_pos_info.json', 'r') as f:
    gene_pos_info_dict = json.loads(f.read())

#  Replace WT sequences with recoded sequences
recoded_record_list = []
for record_id in gb_dict:
    recoded_record = gb_dict[record_id]
    if record_id in chr_dict:
        gene_head = chr_dict[record_id]
        for gene in recode_dict:
            if gene[:2] == gene_head:
                recoded_record.seq = replace_CDS(recoded_record.seq, recode_dict[gene], gene_pos_info_dict[gene])
        recoded_record_list.append(recoded_record)

#  Write recoded sequences into Genbank file
SeqIO.write(recoded_record_list, 'S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.gbk', 'genbank')

#  Convert Genbank file into FASTA
records = []
for record in SeqIO.parse('S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.gbk', 'genbank'):
    records.append(record)
    SeqIO.write(records, 'S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.fasta', 'fasta')
