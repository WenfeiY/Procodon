"""
Decoding binary data from contigs and computing data recovery rates
"""

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
import gzip
from Proteinstorage.decoder import Dnadecoder
from Proteinstorage.seqloader import Seqloader
from Proteinstorage.codec.codec import Codecloader
from Proteinstorage.dnaencoder import Dnaencoder
import os
import json
import sys

def reverse_complement(
    dna_seq: str
):
    """
    Generate reverse complement DNA sequence
    
    Parameter:
        dna_seq: provided DNA sequence in string
    
    Return:
        rev_seq: reverse complement DNA sequence in string
    """
    return str(Seq(dna_seq).reverse_complement())

def process_tblastn(
    xml_file: str, 
    protein_seq: str, 
    contig_file: str
):
    """
    Analyze TBLASTN result in XML format and extract AA-codon alignments
    
    Parameters:
        xml_file: path of TBLASTN XML file (supporting .gz compression)
        protein_seq: reference protein sequence in string
        contig_file: path of contigs file (FASTA)
    
    Return:
        results: list containing positions, AA (1-letter), aligned codon, alignment status (True/False)
            e.g. [[0, 'M', 'ATG', True], ...]
    """
    
    #  Load contigs
    contig_dict = {}
    for record in SeqIO.parse(contig_file, 'fasta'):
        contig_dict[record.id] = str(record.seq).upper()
    
    #  Initialize results list
    results = []
    for i, aa in enumerate(protein_seq, 1):
        results.append([i, aa, "", None])
    
    open_func = gzip.open if xml_file.endswith('.gz') else open
    
    # Open XML file
    with open_func(xml_file, 'rt') as blast_handle:
        #  blast_record_dict contains filtered hsps
        blast_record_dict = {}
        
        #  Traverse BLAST records
        for blast_record in NCBIXML.parse(blast_handle):
            #  remove alignments with low identity
            for alignment in blast_record.alignments:
                filtered_hsps = []
                for hsp in alignment.hsps:
                    if (hsp.identities / (hsp.query_end - hsp.query_start)) > 0.85:
                        filtered_hsps.append([hsp.identities / (hsp.query_end - hsp.query_start), hsp])
                filtered_hsps = sorted(filtered_hsps, key=lambda x: x[0], reverse = True)
                blast_record_dict[alignment.hit_id] = filtered_hsps
        
        #  Traverse contig records
        for contig_id in contig_dict:
            #  Extract contigs with filtered alignments
            if contig_id in blast_record_dict:
                dna_seq = contig_dict[contig_id]
                filtered_hsps = blast_record_dict[contig_id]
                
                #  Parse AA-codon alignments 
                for hsp_item in filtered_hsps:
                    hsp = hsp_item[1]
                    is_negative_strand = hsp.frame[1] < 0
                    #  Protein sequence index (0-based)
                    q_pos = hsp.query_start - 1
                    #  DNA sequence index (0-based)
                    h_pos = hsp.sbjct_start - 1
                    #  Number of traversed codons
                    codon_count = 0
                    
                    #  Traverse sequence in aligned region
                    for pos in range(len(hsp.query)):
                        #  AA in query
                        q_char = hsp.query[pos]
                        #  AA in subject
                        h_char = hsp.sbjct[pos]
                        
                        #  AA in query is not empty in this position
                        if q_char != '-':
                            aa_pos = q_pos
                            q_pos += 1
                        #  AA in query is empty in this position
                        else:
                            aa_pos = None
                        
                        #  AA in subject is not empty in this position
                        if h_char != '-':
                            
                            #  Calculate position of codon
                            #  Negative strand
                            if is_negative_strand:
                                start_idx = hsp.sbjct_end - 3 - codon_count * 3
                                codon_triplet = dna_seq[start_idx:start_idx+3]
                                codon_seq = reverse_complement(codon_triplet)
                            #  Positive strand
                            else:
                                start_idx = h_pos
                                codon_seq = dna_seq[start_idx:start_idx+3]
                                h_pos += 3
                            
                            codon_count += 1
                            
                            #  Link to AA
                            if aa_pos != None:
                                if aa_pos < len(results):
                                    if results[aa_pos][3] == None:
                                        results[aa_pos][2] = codon_seq
                                        #  Mark alignment status
                                        results[aa_pos][3] = True if results[aa_pos][1] == h_char else False
    return results


def decode_from_result(
    results: list
):
    """
    Decode binary string from result parsed from TBLASTN

    Parameter:
        results: list containing positions, AA (1-letter), aligned codon, alignment status (True/False)
            e.g. [[0, 'M', 'ATG', True], ...]

    Return
        binary_string: decoded binary string
    """
    
    binary_string = ''
    for item in results:
        if not item[1] in ['M', 'W', '*']:
            if item[3]:
                binary_string += dec.codon_decode_dict[item[2].upper()]
            else:
                #  Add 'x' to mark the incorrect or empty data
                binary_string += 'x'
    return binary_string

if __name__ == '__main__':
    
    prefix = sys.argv[1]
    cov = sys.argv[2]
    
    #  Load CDS and protein sequences
    seqload= Seqloader(seq_dir = 'S.cerevisiae')
    
    #  Load codec
    codecload=Codecloader()
    codecload.load_codec_from_file('codec_112417.json')
    
    #  Initialize Dnaencoder
    encode = Dnaencoder(codecload.encoding_dict, 'S.cerevisiae.codon_frequency.csv')
    
    #  Initialize Dnadecoder
    dec = Dnadecoder()
    
    #  Load reference of recoded sequence
    with open('S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.json', 'r') as f:
        recoded_dict = json.loads(f.read())
    
    #  Load gene order
    gene_list = list(recoded_dict.keys())
    
    #  Generate complete reference binary data
    binary = ''
    for gene in recoded_dict:
        binary += dec.decode_cds(recoded_dict[gene])

    #  Compute data recovery rate of NGS without correction
    full_binary = ''
    for gene in gene_list:
        
        #  File path containing contigs without correction
        contig_file = './' + prefix + '.contigs/' + cov + '/' + gene + '.contigs.fasta'
        #  File path containing TBLASTN results from contigs without correction
        xml_file = './' + prefix + '/TBLASTN_asm_out/' + cov + '/' + gene + '.tblastn.out'
        
        if os.path.exists(contig_file) and os.path.exists(xml_file):
            if os.path.getsize(xml_file) > 10:
                results = process_tblastn(xml_file, seqload.gene_seq_dict[gene]['Protein'], contig_file)
                full_binary += decode_from_result(results)
            else:
                full_binary += len(dec.decode_cds(seqload.gene_seq_dict[gene]['CDS'])) * 'x'
        else:
            print(gene)
            full_binary += len(dec.decode_cds(seqload.gene_seq_dict[gene]['CDS'])) * 'x'
    
    score = 0
    for i in range(len(binary)):
        if binary[i] == full_binary[i]:
            score += 1
    print(prefix + '\t' + cov + '\tWithout_correction\t' + str(score / len(binary)))

    #  Compute data recovery rate of corrected NGS
    full_binary = ''
    for gene in gene_list:
        
        #  File path containing contigs assembled by corrected reads
        contig_file = './' + prefix + '.contigs.fixed/' + cov + '/' + gene + '.contigs.fasta'
        #  File path containing TBLASTN results from contigs assembled by corrected reads
        xml_file = './' + prefix + '/TBLASTN_asm_out_fix/' + cov + '/' + gene + '.tblastn.out'

        if os.path.exists(contig_file) and os.path.exists(xml_file):
            if os.path.getsize(xml_file) > 10:
                results = process_tblastn(xml_file, seqload.gene_seq_dict[gene]['Protein'], contig_file)
                full_binary += decode_from_result(results)
            else:
                full_binary += len(dec.decode_cds(seqload.gene_seq_dict[gene]['CDS'])) * 'x'
        else:
            print(gene)
            full_binary += len(dec.decode_cds(seqload.gene_seq_dict[gene]['CDS'])) * 'x'
    
    score = 0
    for i in range(len(binary)):
        if binary[i] == full_binary[i]:
            score += 1
    print(prefix + '\t' + cov + '\Corrected\t' + str(score / len(binary)))
