Requirements (Spydes IDE 6.0.5):
    python (3.11.11)
	numpy (2.0.1)
	pandas (2.2.3)
    musicpy (7.11)
    Procodon (Module in this work)

Directory:
	S.cerevisiae (Downloaded dataset): All proteins and CDSs extracted from yeast reference genome.

Files:
	Codec_112417.json: Codec used in this work.
	S.cerevisiae.codon_frequency.csv: The relative codon usage of yeast.
	The_blue_danube.YIL036WC.txt: Notes info of melodic fragment from The blue danube for encoding into YIL036WC.

Scripts:
	01-Simple_notes_encoding.py: Encode simple notes into genes (data and gene are included in the script).
		Usage: 01-Simple_notes_encoding.py

	02-Data_encoding_in_multiple_genes.py: Encode complex notes into gene (data and genes are included in the script).
		Usage: python 02-Data_encoding_in_multiple_genes.py
