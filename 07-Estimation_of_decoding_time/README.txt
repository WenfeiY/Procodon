Requirements (Spydes IDE 6.0.5):
    python (3.11.11)
	numpy (2.0.1)
	pandas (2.2.3)
    musicpy (7.11)
    Procodon (Module in this work)

Directory:
	S.cerevisiae (Copy from other directories): All proteins and CDSs extracted from yeast reference genome.

Files:
	S.cerevisiae.aa_frequency_dict.json: The total number of each amino acid accross all proteins in yeast.

	S.cerevisiae.codec_0.json: The original codec of yeast.

	S.cerevisiae.codon_frequency.csv: The relative codon usage of yeast.

Scripts:
	01-Estimation_of_decoding_time_for_single_gene.py: Python script for calculating the time cost for decoding given number of genes
		Usage: python 01-Estimation_of_decoding_time_for_single_gene.py <gene_number>

	02-Computing_cracking_complexity.py: Python script for computing cracking complexity
		Usage: python 02-Computing_cracking_complexity.py
