## Requirements (Ubuntu 18.04.6 LTS)
- **Python**: >= 3.8
- **DWGSIM**: 1.1.14
- **BLAST+**: >= 2.5.0
- **Biopython**: >= 1.85
- **seqkit**: >= 2.8.2
- **Spades**: >= 3.15.4
- **Procodon**: Python module in <01-Procodon>

## Main scripts (run directly)
 - **01-Generation_and_decoding_of_mutated_NGS.sh**: Generate simulated NGS data with different mutation rates, build blast database and decode data from NGS.
<br>&emsp;&emsp;Usage: ```bash 01-Generation_and_decoding_of_mutated_NGS.sh```

 - **02-Reads_correction_by_AA.sh**: Correct reads by TBLASTN results.
<br>&emsp;&emsp;Usage: ```bash 02-Reads_correction_by_AA.sh```

 - **03-Assembly_of_corrected_reads.sh**: Assemble corrected reads into contigs.
<br>&emsp;&emsp;Usage: ```bash 03-Assembly_of_corrected_reads.sh```

 - **04-Intron_removal.sh**: Remove introns from the assembled contigs and generate bare CDS.
<br>&emsp;&emsp;Usage: ```bash 04-Intron_removal.sh```

 - **05-Copying_contigs.sh**: Copy contigs to new folders for downloading.
<br>&emsp;&emsp;Usage: ```bash 05-Copying_contigs.sh```

 - **06-TBLASTN_contigs.sh**: Align protein sequence with contigs.
<br>&emsp;&emsp;Usage: ```bash 06-TBLASTN_contigs.sh```

 - **07-Computing_recovery_rate.py**: Decode binary data from contigs and TBLASTN results, calculating data recovery rate.
<br>&emsp;&emsp;Usage: ```python 07-Computing_recovery_rate.py <prefix> <coverage>```
<br>&emsp;&emsp;&emsp;&emsp;e.g. ```python 07-Computing_recovery_rate.py Train_Ng_0.05 30x```

**#Note**
1. Run the next script after the previous script's task has completed.
2. Sometimes the assembly or other tasks will be terminated when too many programs are running. In this status, the attached scripts can be used for re-running.

## Attached scripts (called by main scripts)
 - **Attached-01-filter_TBLASTN.py**: Calculate threshold of bit-score and filter TBLASTN result. Called by `Attached-02-simple_gene_cov_process.sh` and `Attached-03-intron_gene_cov_process.sh`.
<br>&emsp;&emsp;Usage: ```python Attached-01-filter_TBLASTN.py <protein_seq_path> <TBLASTN_path> <filtered_TBLASTN_path>```

 - **Attached-02-simple_gene_cov_process.sh**: Align protein sequence with reads by TBLASTN. Align flanking sequences with reads by BLASTN. Extract reads ID from BLAST results. Extract FASTQ subset aligned to the given gene. Assemble reads into contigs. Called by `Attached-04-simple_cov_test.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-02-simple_gene_cov_process.sh <prefix> <coverage> <gene_name>```

 - **Attached-03-intron_gene_cov_process.sh**: Align protein sequence with reads by TBLASTN. Align flanking sequences with reads by BLASTN. Align intron sequences with reads by BLASTN. Extract reads ID from BLAST results. Extract FASTQ subset aligned to the given gene. Assemble reads into contigs. Called by `Attached-05-intron_cov_test.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-03-intron_gene_cov_process.sh <prefix> <coverage> <gene_name>```

 - **Attached-04-simple_cov_test.sh**: Decode sequences (contigs) of genes without introns with given prefix and coverage. Require `simple_gene_list.txt`. Called by `01-Generation_and_decoding_of_mutated_NGS.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-04-simple_cov_test.sh <prefix> <coverage>```

 - **Attached-05-intron_cov_test.sh**: Decode sequences (contigs) of genes with introns with given prefix and coverage. Require `intron_gene_list.txt`. Called by `01-Generation_and_decoding_of_mutated_NGS.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-05-intron_cov_test.sh <prefix> <coverage>```

 - **Attached-06-NGS_fix.py**: Correct NGS data in FASTQ with given protein sequence, TBLASTN result, reads in FASTQ and output FASTQ. Called by `Attached-07-NGS_fix_process.sh`.
<br>&emsp;&emsp;Usage: ```python Attached-06-NGS_fix.py <protein_seq_path> <TBLASTN_path> <FQ_path> <filtered_FQ_path>```

 - **Attached-07-NGS_fix_process.sh**: Correct NGS data with given prefix, coverage and gene name. Called by `Attached-08-NGS_cov_fix.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-07-NGS_fix_process.sh <prefix> <coverage> <gene>```

 - **Attached-08-NGS_cov_fix.sh**: Correct reads by AA with given prefix and coverage. Called by `02-Reads_correction_by_AA.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-08-NGS_cov_fix.sh <prefix> <coverage>```

 - **Attached-09-gene_fix_asm.sh**: Assemble corrected reads into contigs with given prefix, coverage and gene name. Called by `Attached-10-cov_fix_asm.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-09-gene_fix_asm.sh <prefix> <coverage> <gene>```

 - **Attached-10-cov_fix_asm.sh**: Assemble corrected reads into contigs with given prefix and coverage. Called by `03-Assembly_of_corrected_reads.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-10-cov_fix_asm.sh <prefix> <coverage>```

 - **Attached-11-remove_intron_from_contigs.py**: Remove intron from contigs with given BLASTN result path, intron sequences path, contig sequences path and output intron removed contig sequences path. Called by `Attached-12-remove_intron_asm.sh` and `Attached-14-remove_intron_corrected_asm.sh`.
<br>&emsp;&emsp;Usage: ```python Attached-11-remove_intron_from_contigs.py <blastn_file> <intron_file> <contig_file> <out_cds_file>```

 - **Attached-12-remove_intron_asm.sh**: Remove intron from contigs (without correction) with given prefix, coverage and gene name. Called by `Attached-13-remove_intron_cov.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-12-remove_intron_asm.sh <prefix> <coverage> <gene>```

 - **Attached-13-remove_intron_cov.sh**: Remove intron from contigs (without correction) with given prefix and coverage. Called by `04-Intron_removal.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-13-remove_intron_cov.sh <prefix> <coverage>```

 - **Attached-14-remove_intron_corrected_asm.sh**: Remove intron from contigs (corrected) with given prefix, coverage and gene name. Called by `Attached-15-remove_intron_corrected_cov.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-14-remove_intron_corrected_asm.sh <prefix> <coverage> <gene>```

 - **Attached-15-remove_intron_corrected_cov.sh**: Remove intron from contigs (corrected) with given prefix and coverage. Called by `04-Intron_removal.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-15-remove_intron_corrected_cov.sh <prefix> <coverage>```

 - **Attached-16-cp_asm.sh**: Copy contigs (without correction) to a new folder. Called by `05-Copying_contigs.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-16-cp_asm.sh```

 - **Attached-17-cp_corrected_asm.sh**: Copy contigs (corrected) to a new folder. Called by `05-Copying_contigs.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-17-cp_corrected_asm.sh```

 - **Attached-18-TBLASTN_asm.sh**: Align protein sequence with contigs (without correction) with given prefix, coverage and gene name. Called by `Attached-20-TBLASTN_cov.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-18-TBLASTN_asm.sh <prefix> <coverage> <gene>```

 - **Attached-19-TBLASTN_corrected_asm.sh**: Align protein sequence with contigs (corrected) with given prefix, coverage and gene name. Called by `Attached-21-TBLASTN_corrected_cov.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-19-TBLASTN_corrected_asm.sh <prefix> <coverage> <gene>```

 - **Attached-20-TBLASTN_cov.sh**: Align protein sequence with contigs (without correction) with given prefix and coverage. Called by `06-TBLASTN_contigs.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-20-TBLASTN_cov.sh <prefix> <coverage>```

 - **Attached-21-TBLASTN_corrected_cov.sh**: Align protein sequence with contigs (corrected) with given prefix and coverage. Called by `06-TBLASTN_contigs.sh`.
<br>&emsp;&emsp;Usage: ```bash Attached-21-TBLASTN_corrected_cov.sh <prefix> <coverage>```

## Attached files
 - **S.cerevisiae/**: directories containing all CDS and protein sequences of yeast, which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **S.cerevisiae.Protein/**: Directory containing all protein sequences extracted from *S. cerevisiae* reference genome (GCF_000146045.2), which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **S.cerevisiae.Gene_flanking/**: Directory containing all gene flanking sequences (upstream and downstream 50 bp) extracted from *S. cerevisiae* reference genome (GCF_000146045.2), which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **S.cerevisiae.Intron/**: Directory containing all intron sequences extracted from *S. cerevisiae* reference genome (GCF_000146045.2), which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **all_gene_list.txt**: All systematic names of genes used for storing images.
 - **simple_gene_list.txt**: Systematic names of genes without introns used for storing images.
 - **intron_gene_list.txt**: Systematic names of genes with introns used for storing images.
 - **S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.fasta**: Genome sequences with recoded CDSs.

