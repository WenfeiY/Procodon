"""
Extracting gene information (gene names, directions, CDS positions and introns) from yeast genome
"""

import json
from Bio import SeqIO

gene_pos_info_dict = {}

#  Parse all chromosomes
for record in SeqIO.parse('Saccharomyces.cerevisiae.gbk', 'genbank'):
    for fea in record.features:
        if fea.type == 'CDS':
            #  The systematic names are stored in term "locus_tag"
            gene_name = fea.qualifiers['locus_tag'][0]
            
            gene_pos_info_dict[gene_name] = {
                'CDS_position': [],
                'direction': fea.location.strand,
            }
            
            if len(fea.location.parts) == 1:  
                
                #  CDS without intron
                gene_pos_info_dict[gene_name]['with_intron'] = False
                gene_pos_info_dict[gene_name]['CDS_position'] = [fea.location.start, fea.location.end]
            
            else:
                
                #  CDS with intron
                gene_pos_info_dict[gene_name]['with_intron'] = True
                if fea.location.strand == 1:
                    loc_parts = sorted(fea.location.parts, key=lambda x: x.start)
                else:
                    loc_parts = sorted(fea.location.parts, key=lambda x: x.start, reverse = True)
                
                #  Record positions of each part of CDS
                for part in loc_parts:
                    gene_pos_info_dict[gene_name]['CDS_position'].append([part.start, part.end])

#  Write into file
with open('sc_gene_pos_info.json', 'w') as hd:
    json.dump(gene_pos_info_dict, hd)
