"""
Load CDS and protein sequences for encoding and decoding
"""
from Bio import SeqIO
import os
import json

class Seqloader:
    def __init__(
        self,
        seq_dir: str = None,
        gene_path_dict: dict = None,
        gene_path_dict_path: str = None
    ):
        """
        Load CDS and protein sequences for protein encoding.
        
        The structure of sequence directory for automatic reading:
            /Folder(species)
                ├── CDS
                │   ├── gene1.<surfix, optional>.fa  e.g. gene1.CDS.fa; gene1.fa
                │   ├── gene2.<surfix, optional>.fa
                │   └── ...
                └── Protein
                    ├── gene1.<surfix, optional>.fa  e.g. gene1.Protein.fa; gene1.fa
                    ├── gene2.<surfix, optional>.fa
                    └── ...
        
        Alternatively, you can provide the dictionary containing CDS and protein fasta path:
            gene_path_dict:
                {
                    'gene1': {
                        'CDS': gene1_CDS_FASTA_path,
                        'Protein': gene1_Protein_FASTA_path
                        },
                    'gene2': {
                        'CDS': gene2_CDS_FASTA_path,
                        'Protein': gene2_Protein_FASTA_path
                        },
                    ...
                }
        
        Parameter:
            seq_dir: standard folder structure of CDS and protein sequences for protein encoding
            gene_path_dict: dictionary containing CDS and protein fasta path
        """
        
        if not (seq_dir or gene_path_dict, gene_path_dict_path):
            raise ValueError(
                'Please provide seq_dir or gene_path_dict or gene_path_dict_path!'
            )
        
        if not gene_path_dict:
            if gene_path_dict_path:
                with open(gene_path_dict_path, 'r') as f:
                    gene_path_dict = json.loads(f.read())
            else:
                self.seq_dir = seq_dir
                gene_path_dict = {}
                for file in os.listdir(os.path.join(self.seq_dir, 'CDS')):
                    gene_name = file.split('.')[0]
                    gene_path_dict[gene_name] = {
                        'CDS': os.path.join(self.seq_dir, 'CDS', file),
                        'Protein': os.path.join(self.seq_dir, 'Protein', gene_name + '.Protein.fa')
                    }
        
        self.gene_path_dict = gene_path_dict
        
        gene_seq_dict = {}
        for gene_name in self.gene_path_dict:
            cds_seq = str(SeqIO.read(self.gene_path_dict[gene_name]['CDS'], 'fasta').seq).upper()
            prot_seq = str(SeqIO.read(self.gene_path_dict[gene_name]['Protein'], 'fasta').seq).upper()
            gene_seq_dict[gene_name] = {
                'CDS': cds_seq,
                'Protein': prot_seq
            }
        self.gene_seq_dict = gene_seq_dict


    def get_prot_seqs(
        self,
        gene_list: list = None
    ) -> dict:
        
        """
        Extract protein sequences from gene_seq_dict by gene_list
        
        Parameter:
            gene_list: list containing gene names, e.g. ['gene1', 'gene2', 'gene3', ...]
        
        Return:
            prot_seqs: list containing extracted protein sequences, e.g. ['MHHMIC...', 'MHREEC...', ...]
        """
        
        prot_seqs = []
        if not gene_list:
            gene_list = list(self.gene_seq_dict.keys())
        for gene_name in gene_list:
            prot_seqs.append(self.gene_seq_dict[gene_name]['Protein'])
        return prot_seqs


    def get_cds_seqs(
        self,
        gene_list: list = None
    ):
        
        """
        Extract CDS sequences from gene_seq_dict by gene_list
        
        Parameter:
            gene_list: list containing gene names, e.g. ['gene1', 'gene2', 'gene3', ...]
        
        Return:
            cds_seqs: list containing extracted CDS sequences, e.g. ['ATGCGA...', 'ATGTTG...', ...]
        """
        
        cds_seqs = []
        if not gene_list:
            gene_list = list(self.gene_seq_dict.keys())
        for gene_name in gene_list:
            cds_seqs.append(self.gene_seq_dict[gene_name]['CDS'])
        return cds_seqs
    
    def get_sub_seq_dict(
        self,
        gene_list: list
    ) -> dict:
        
        """
        Extract protein and CDS sequences from gene_seq_dict by gene_list
        
        Parameter:
            gene_list: list containing gene names, e.g. ['gene1', 'gene2', 'gene3', ...]
        
        Return:
            sub_seq_dir: dictionary containing extracted protein and CDS sequences
                e.g.{
                    'gene20': {
                        'CDS': 'ATG...',
                        'Protein': 'MHHMIC'
                        },
                    'gene35': {
                        'CDS': 'ATG...',
                        'Protein': 'MHREEC'
                        },
                    ...
                    }
        """
        
        sub_seq_dir = {}
        if gene_list:
            for gene_name in gene_list:
                sub_seq_dir[gene_name] = self.gene_seq_dict[gene_name]
            return sub_seq_dir
        return sub_seq_dir
