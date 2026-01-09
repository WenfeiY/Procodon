"""
Extract orthologous info from genbank files
"""

from Bio import SeqIO
import re

### Extract orthogs from Nakaseomyces glabratus
#----------------------------------------------------------------------------------------------------------------
#  Note: some ortholog infos in GCA_000002545.2 are incorrect. e.g. YBL091C-A was recorded as YBL091Ca
def extract_gene_names_Ng(
    text: str
):
    """
    Extract the systematic gene names from text (SeqFeature.qualifiers['note'])

    Parameter:
        text: given text, norrmally the "note" term from Genbank annotation
    
    Return:
        gene_name_list: list containing extracted systematic names
    """
    
    pattern = r'Y[A-P][LR]\d{3}[WC](?:[A-F])?'
    list1 = re.findall(pattern, text)
    pattern = r'Y[A-P][LR]\d{3}[WC]-[A-F]?'
    list2 = re.findall(pattern, text)
    return list1 + list2

def extract_ortho_info_Ng(
    genome_file_path: str, 
    ortho_gene_info_path: str
):
    """
    Extract orthologous info from genbank file of Nakaseomyces glabratus (GCA_000002545.2)

    Parameters:
        genome_file_path: genbank file path (containing all chromosomes), str
        ortho_gene_info_path: output orthologs info path, str

    Output:
        ortholog names and positions
    """
    
    ortho_gene_pos_list = []
      
    for record in SeqIO.parse(genome_file_path, 'genbank'):
        for fea in record.features:
            if fea.type == 'CDS' and 'locus_tag' in fea.qualifiers and 'note' in fea.qualifiers:
                ortho_genes = extract_gene_names_Ng(fea.qualifiers['note'][0].upper())
                gene_locus = fea.qualifiers['locus_tag'][0]
                if len(ortho_genes) == 1:
                    ortho_gene = ortho_genes[0]
                    if len(ortho_gene) == 8:
                        ortho_gene = ortho_gene[:7] + '-' + ortho_gene[-1]
                    ortho_gene_pos_list.append([record.id, gene_locus, ortho_gene, str(fea.location.start), str(fea.location.end)])

    with open(ortho_gene_info_path, 'w') as hd:
        for gene_info in ortho_gene_pos_list:
            hd.write('\t'.join(gene_info) + '\n')

#  Input and output path
genome_file_path = 'Nakaseomyces.glabratus.gbff'
ortho_gene_info_path = 'Ng.ortho_gene_info.txt'
 
#  Extract orthologs info and write into file
extract_ortho_info_Ng(genome_file_path, ortho_gene_info_path)



### Extract orthogs from Kluyveromyces lactis and Tetrapisispora phaffli
#----------------------------------------------------------------------------------------------------------------
def extract_gene_names(
    text: str
):
    """
    Extract the systematic gene names from text (SeqFeature.qualifiers['note'])

    Parameter:
        text: given text, norrmally the "note" term from Genbank annotation
    
    Return:
        gene_name_list: list containing extracted systematic names
    """
    pattern = r'Y[A-P][LR]\d{3}[WC](?:-[A-F])?'
    return re.findall(pattern, text)

def extract_ortho_info(
    genome_file_path: str, 
    ortho_gene_info_path: str
):
    """
    Extract orthologous info from genbank file

    Parameters:
        genome_file_path: genbank file path (containing all chromosomes), str
        ortho_gene_info_path: output orthologs info path, str

    Output:
        ortholog names and positions
    """
    
    ortho_gene_pos_list = []
      
    for record in SeqIO.parse(genome_file_path, 'genbank'):
        for fea in record.features:
            if fea.type == 'CDS' and 'locus_tag' in fea.qualifiers and 'note' in fea.qualifiers:
                ortho_genes = extract_gene_names(fea.qualifiers['note'][0].upper())
                gene_locus = fea.qualifiers['locus_tag'][0]
                if len(ortho_genes) == 1:
                    ortho_gene = ortho_genes[0]
                    if len(ortho_gene) == 8:
                        ortho_gene = ortho_gene[:7] + '-' + ortho_gene[-1]
                    ortho_gene_pos_list.append([record.id, gene_locus, ortho_gene, str(fea.location.start), str(fea.location.end)])

    with open(ortho_gene_info_path, 'w') as hd:
        for gene_info in ortho_gene_pos_list:
            hd.write('\t'.join(gene_info) + '\n')

#  Kluyveromyces lactis
#  Input and output path
genome_file_path = 'Kluyveromyces.lactis.gbff'
ortho_gene_info_path = 'Kl.ortho_gene_info.txt'
 
#  Extract orthologs info and write into file
extract_ortho_info(genome_file_path, ortho_gene_info_path)

#  Tetrapisispora phaffli
#  Input and output path
genome_file_path = 'Tetrapisispora.phaffli.gbff'
ortho_gene_info_path = 'Tp.ortho_gene_info.txt'
 
#  Extract orthologs info and write into file
extract_ortho_info(genome_file_path, ortho_gene_info_path)
#----------------------------------------------------------------------------------------------------------------
