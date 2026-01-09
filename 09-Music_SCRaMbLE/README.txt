Requirements (Spydes IDE 6.0.5):
    python (3.11.11)
	numpy (2.0.1)
	pandas (2.2.3)
    musicpy (7.11)
    Procodon (Module in this work)

Directories:
	S.cerevisiae (Downloaded dataset): All proteins and CDSs extracted from yeast reference genome.
	synIXR.SCRaMbLE.structure: The structures of LoxPsym units of SCRaMbLEd strains
		Strain + ".LU.txt": Each line contains a number representing LoxPsym unit (1~43). Negative numbers represent reversed LoxPsym units.
	synIXR.SCRaMbLE.gene.sort: Gene order of SCRaMbLEd strains
		Strain + ".lu_gene.txt": Each line contains a gene name. "-" represent reversed gene.

Files:
	Codec_112417.json: Codec used in this work.
	S.cerevisiae.codon_frequency.csv: The relative codon usage of yeast.
	The_blue_danube.YIL036WC.txt: Notes info of melodic fragment from The blue danube for encoding into YIL036WC.
	Canon.synIXR.150_bit_per_note.txt: Note infos of melodic fragment from Canon

Scripts:
	01-Encoding_Canon_in_synIXR_with_extended_bin.py: Read notes, extend the binary string (150-bit for 1 note) and encode complete notes into synIXR.
		Usage: 01-Encoding_Canon_in_synIXR_with_extended_bin.py

	02-Generation_of_SCRaMbLEd_music_with_extended_bin.py: Read structures of LoxPsym units and decode rearranged music from SCRaMbLEd strains.
		Usage: python 02-Generation_of_SCRaMbLEd_music_with_extended_bin.py
