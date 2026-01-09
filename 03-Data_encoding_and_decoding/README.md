Requirements (Spydes IDE 6.0.5):
    python (3.11.11)
	biopython (1.85)
    Procodon (Module in this work)

Directory:
	S.cerevisiae (Downloaded dataset): All proteins and CDSs extracted from yeast reference genome.

Files:
	Codec_112417.json: Codec used in this work.
	S.cerevisiae.codon_frequency.csv: Relative codon usage of yeast.

Scripts:
	01-Data_encoding_in_single_gene.py: Encode data in single gene (data and gene are included in the script).
		Usage: python 01-Data_encoding_in_single_gene.py
	
	02-Data_encoding_in_multiple_genes.py: Encode data in multiple genes (data and genes are included in the script).
		Usage: python 02-Data_encoding_in_multiple_genes.py

	03-Data_decoding_from_single_CDS.py: Decode data from single CDS (sequence is included in the script).
		Usage: python 03-Data_decoding_from_single_CDS.py
	
	04-Data_decoding_from_multiple_CDSs.py: Decode data from multiple CDSs (sequences are included in the script).
		Usage: python 04-Data_decoding_from_multiple_CDSs.py
