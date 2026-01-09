"""
Extract CDS from annotation & sequence files downloaded from https://syntheticyeast.github.io/sc2-0/data/
"""

from Bio import SeqIO
import os
import pandas as pd
import json
from Bio.SeqFeature import SeqFeature, FeatureLocation

#  Extract sequences and annotations from data downloaded from Sc2.0 website

def split_file_at_fasta(input_file, output_file1, output_file2):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    found_idx = -1
    for idx, line in enumerate(lines):
        #  Find the line "##FASTA" which separate annotations and sequences
        if line.strip() == "##FASTA":
            found_idx = idx
            break

    #  Write annotations
    with open(output_file1, 'w') as f1:
        f1.writelines(lines[:found_idx] if found_idx != -1 else lines)
    
    #  Write sequences
    if found_idx != -1:
        with open(output_file2, 'w') as f2:
            f2.writelines(lines[found_idx + 1:])
    else:
        #  Create empty file if "##FASTA" not dound
        open(output_file2, 'w').close()


#  Parse all synthetic chromosomes
for file in os.listdir('synthetic_chromosomes'):
    split_file_at_fasta('synthetic_chromosomes/' + file, 'synthetic_chromosomes/' + file.split('.')[0] + '.gff3', 'synthetic_chromosomes/' + file.split('.')[0] + '.fa')


#  Integrate annotations and sequences into genbank files
def generate_gb_from_fa_gff(in_fa, in_gff3, out_gb):
    gff_df = pd.read_csv(in_gff3, sep = '\t', header = None, comment = '#')
    gff_df.columns = ['ref', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    record = SeqIO.read(in_fa, 'fasta')
    record.annotations['molecule_type'] = 'DNA'
    record.annotations['topology'] = 'linear'
    for indexs in gff_df.index:
        if gff_df.loc[indexs, "type"] == 'CDS':
            sf_type = gff_df.loc[indexs, "type"]
            if gff_df.loc[indexs, "strand"] == '+':
                sf_strand = 1
            else:
                sf_strand = -1
            sf = SeqFeature(FeatureLocation(int(gff_df.loc[indexs, "start"]) - 1, int(gff_df.loc[indexs, "end"]), strand = sf_strand), type = sf_type)
            attr_dict = {}
            for item in gff_df.loc[indexs, "attributes"].split(';')[0:-1]:
                attr_dict[item.split("=")[0]] = item.split("=")[1]
            if "ID" in attr_dict:
                sf.qualifiers['label'] = attr_dict["ID"]
                sf.qualifiers['locus_tag'] = attr_dict["ID"]
            record.features.append(sf)
    SeqIO.write(record, out_gb, "genbank")


for file in os.listdir('synthetic_chromosomes'):
    generate_gb_from_fa_gff('synthetic_chromosomes/' + file.split('.')[0] + '.fa', 'synthetic_chromosomes/' + file.split('.')[0] + '.gff3', 'synthetic_chromosomes/genbank/' + file.split('.')[0] + '.gb')

#  Extract CDS sequences from genbank files
CDS_seq_dict = {}

for file in os.listdir('synthetic_chromosomes/genbank'):
    record = SeqIO.read('synthetic_chromosomes/genbank/' + file, 'genbank')
    gene_list = []
    for fea in record.features:
        if fea.type == 'CDS':
            gene_list.append(fea.qualifiers['locus_tag'][0])
            CDS_seq_dict[fea.qualifiers['locus_tag'][0][:-4]] = str(fea.extract(record.seq)).upper()

#  Write into json file
with open('syn_chr_CDS_seq_dict.json', 'w') as hd:
    json.dump(CDS_seq_dict, hd)

#  Finding genes with separate regions
with open('syn_chr_CDS_seq_dict.json', 'r') as f:
    CDS_seq_dict = json.loads(f.read())

for gene in CDS_seq_dict:
    if not len(CDS_seq_dict[gene]) % 3 == 0:
        print(gene)


with open('syn_chr_CDS_seq_dict.json', 'r') as f:
    CDS_seq_dict = json.loads(f.read())

for gene in CDS_seq_dict:
    if not len(CDS_seq_dict[gene]) % 3 == 0:
        print(gene)
    if '_' in gene:
        print(gene)

"""
Some CDSs were manually fixed (refer to README.md)
"""
