Requirements (Spydes IDE 6.0.5):
    python (3.11.11)
    biopython (1.85)
    Procodon (Module in this work)

Directory:
	Images: all images used for storage simulation

Files:
	S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.json: recoded CDS sequences storing images from "The_arrival_of_a_train" in json format
	S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.gbk: recoded yeast genome storing images from "The_arrival_of_a_train" with annotations
	S.cerevisiae.The_arrival_of_a_train.Ng-bio-index.fasta: recoded yeast genome storing images from "The_arrival_of_a_train" in FASTA format

	S.cerevisiae.Girl_with_a_Pearl_Earring.Kl-bio-index.json: recoded CDS sequences storing image "Girl_with_a_Pearl_Earring" in json format
	S.cerevisiae.Girl_with_a_Pearl_Earring.Kl-bio-index.gbk: recoded yeast genome storing image "Girl_with_a_Pearl_Earring" with annotations
	S.cerevisiae.Girl_with_a_Pearl_Earring.Kl-bio-index.fasta: recoded yeast genome storing image "Girl_with_a_Pearl_Earring" in FASTA format

	S.cerevisiae.The_Scream.Tp-bio-index.json: recoded CDS sequences storing image "The_Scream" in json format
	S.cerevisiae.The_Scream.Tp-bio-index.gbk: recoded yeast genome storing image "The_Scream" with annotations
	S.cerevisiae.The_Scream.Tp-bio-index.fasta: recoded yeast genome storing image "The_Scream" in FASTA format

	Ng.bio-index.txt: gene order of bio-index from Nakaseomyces glabratus
	Kl.bio-index.txt: gene order of bio-index from Kluyveromyces lactis
	Tp.bio-index.txt: gene order of bio-index from Tetrapisispora phaffli

Scripts:
	01-Extraction_of_yeast_gene_info.py: Extract genome information
		Usage: python 01-Extraction_of_yeast_gene_info.py
		Output: sc_gene_pos_info.json containing gene information (gene names, directions, CDS positions and introns) of yeast genome

	02-Image_storage.The_arrival_of_a_train.Ng-bio-index.py: Store images from "The arrival of a train" with bio-index from N. glabratus
		Usage: python 02-Image_storage.The_arrival_of_a_train.Ng-bio-index.py
		Output: recoded yeast genome storing images from "The_arrival_of_a_train", with genes sorted as Nakaseomyces glabratus

	03-Image_storage.Girl_with_a_Pearl_Earring.Kl-bio-index.py: Store image "Girl with a Pearl Earring" with bio-index from K. lactis
		Usage: python 03-Image_storage.Girl_with_a_Pearl_Earring.Kl-bio-index.py
		Output: recoded yeast genome storing image "Girl_with_a_Pearl_Earring", with genes sorted as Kluyveromyces lactis

	04-Image_storage.The_Scream.Tp-bio-index.py: Store image "The Scream" with bio-index from T. phaffli
		Usage: python 04-Image_storage.The_Scream.Tp-bio-index.py
		Output: recoded yeast genome storing image "The_Scream", with genes sorted as Tetrapisispora phaffli
