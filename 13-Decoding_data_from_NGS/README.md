## Requirements (Ubuntu 18.04.6 LTS)
- **Python**: >= 3.8
- **DWGSIM**: 1.1.14
- **BLAST+**: >= 2.5.0
- **Biopython**: >= 1.85
- **seqkit**: >= 2.8.2
- **Spades**: >= 3.15.4
- **Procodon**: Python module in <01-Procodon>

## Main scripts (run directly)
 - **01-Processing_and_decoding_NGS.sh**: Build blast database and decode data from given NGS data.
<br>&emsp;&emsp;Usage: ```bash 01-Processing_and_decoding_NGS.sh <prefix> <in_R1_fastq_file> <in_R1_fastq_file>```

 - **02-Intron_removal.sh**: Remove introns from the assembled contigs and generate bare CDS.
<br>&emsp;&emsp;Usage: ```bash 04-Intron_removal.sh <prefix>```

 - **03-Copying_contigs.sh**: Copy contigs to new folders for downloading.
<br>&emsp;&emsp;Usage: ```bash 05-Copying_contigs.sh <prefix>```

 - **04-TBLASTN_contigs.sh**: Align protein sequence with contigs.
<br>&emsp;&emsp;Usage: ```bash 04-TBLASTN_contigs.sh <prefix>```

**#Note**
1. Run the next script after the previous script's task has completed.
2. Sometimes the assembly or other tasks will be terminated when too many programs are running. In this status, the attached scripts can be used for re-running.

## Attached scripts (called by main scripts)
 - **Attached-01-filter_TBLASTN.py**: Calculate threshold of bit-score and filter TBLASTN result. Called by `Attached-02-simple_gene_process.sh` and `Attached-03-intron_gene_process.sh`.
<br>&emsp;&emsp;Usage: ```python Attached-01-filter_TBLASTN.py <protein_seq_path> <TBLASTN_path> <filtered_TBLASTN_path>```

 - **Attached-02-simple_gene_process.sh**: Align protein sequence with reads by TBLASTN. Align flanking sequences with reads by BLASTN. Extract reads ID from BLAST results. Extract FASTQ subset aligned to the given gene. Assemble reads into contigs. Called by `Attached-04-simple_process.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-02-simple_gene_process.sh <prefix> <gene_name> <in_R1_fastq_file> <in_R1_fastq_file>```

 - **Attached-03-intron_gene_process.sh**: Align protein sequence with reads by TBLASTN. Align flanking sequences with reads by BLASTN. Align intron sequences with reads by BLASTN. Extract reads ID from BLAST results. Extract FASTQ subset aligned to the given gene. Assemble reads into contigs. Called by `Attached-05-intron_process.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-03-intron_gene_process.sh <prefix> <gene_name> <in_R1_fastq_file> <in_R1_fastq_file>```

 - **Attached-04-simple_process.sh**: Decode sequences (contigs) of genes without introns with given prefix and NGS data. Require `simple_gene_list.txt`. Called by `01-Processing_and_decoding_NGS.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-04-simple_process.sh <prefix> <in_R1_fastq_file> <in_R1_fastq_file>```

 - **Attached-05-intron_process.sh**: Decode sequences (contigs) of genes with introns with given prefix and NGS data. Require `intron_gene_list.txt`. Called by `01-Processing_and_decoding_NGS.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-05-intron_process.sh <prefix> <in_R1_fastq_file> <in_R1_fastq_file>```

 - **Attached-06-remove_intron_from_contigs.py**: Remove intron from contigs with given BLASTN result path, intron sequences path, contig sequences path and output intron removed contig sequences path. Called by `Attached-07-remove_intron_gene_asm.sh`.
<br>&emsp;&emsp;Usage: ```python Attached-06-remove_intron_from_contigs.py <blastn_file> <intron_file> <contig_file> <out_cds_file>```

 - **Attached-07-remove_intron_gene_asm.sh**: Remove intron from contigs with given prefix and gene name. Called by `Attached-08-remove_intron_asm.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-07-remove_intron_gene_asm.sh <prefix> <gene>```

 - **Attached-08-remove_intron_asm.sh**: Remove intron from contigs with given prefix. Called by `02-Intron_removal.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-08-remove_intron_asm.sh <prefix>```

 - **Attached-09-cp_asm.sh**: Copy contigs to a new folder. Called by `03-Copying_contigs.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-09-cp_asm.sh <prefix>```

 - **Attached-10-TBLASTN_gene_asm.sh**: Align protein sequence with contigs with given prefix and gene name. Called by `Attached-11-TBLASTN_asm.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-10-TBLASTN_gene_asm.sh <prefix> <gene>```

 - **Attached-11-TBLASTN_asm.sh**: Align protein sequence with contigs with given prefix. Called by `04-TBLASTN_contigs.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-11-TBLASTN_asm.sh <prefix>```

## Attached files
 - **IR_deletion+pSYN10.R1.fastq.gz**: R1 Fastq file of ChrIR deleted with pSYN10 (Fig.6d), which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
  - **IR_deletion+pSYN10.R2.fastq.gz**: R2 Fastq file of ChrIR deleted with pSYN10 (Fig.6d), which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **S.cerevisiae/**: directories containing all CDS and protein sequences of yeast, which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **S.cerevisiae.Protein/**: Directory containing all protein sequences extracted from *S. cerevisiae* reference genome (GCF_000146045.2), which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **S.cerevisiae.Gene_flanking/**: Directory containing all gene flanking sequences (upstream and downstream 50 bp) extracted from *S. cerevisiae* reference genome (GCF_000146045.2), which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **S.cerevisiae.Intron/**: Directory containing all intron sequences extracted from *S. cerevisiae* reference genome (GCF_000146045.2), which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **all_gene_list.txt**: All systematic names of genes used for storing images.
 - **simple_gene_list.txt**: Systematic names of genes without introns used for storing images.
 - **intron_gene_list.txt**: Systematic names of genes with introns used for storing images.
