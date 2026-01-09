Requirements for python script (Spydes IDE 6.0.5):
    python (3.11.11)
	biopython (1.85)

Requirements for R script (RStudio 2024.04.2+764):
    R (4.4.1)
	circlize (0.4.16)

Files:
	Nakaseomyces.glabratus.gbff: Genbank file containing all chromosomes of Nakaseomyces glabratus (GCA_000002545.2)
	Kluyveromyces.lactis.gbff: Genbank file containing all chromosomes of Kluyveromyces lactis (GCF_000002515.2)
	Tetrapisispora.phaffli.gbff: Genbank file containing all chromosomes of Tetrapisispora phaffli (GCF_000236905.1)

	Ng.chr_df.txt: Chromosome names, order, starts and ends of Nakaseomyces glabratus
	Kl.chr_df.txt: Chromosome names, order, starts and ends of Kluyveromyces lactis
	Tp.chr_df.txt: Chromosome names, order, starts and ends of Tetrapisispora phaffli

	Ng.sc_link.txt: Chromosomes, starts and ends of Orthologs between Nakaseomyces glabratus and Saccharomyces cerevisiae
	Kl.sc_link.txt: Chromosomes, starts and ends of Orthologs between Kluyveromyces lactis and Saccharomyces cerevisiae
	Tp.sc_link.txt: Chromosomes, starts and ends of Orthologs between Tetrapisispora phaffli and Saccharomyces cerevisiae

Scripts:
	01-Extraction_of_orthologs_info.py: Extract orthologous infos from 3 species
		Usage: python 01-Extraction_of_orthologs_info.py
	
	02-Creating_syntenic_plot.R: Create syntenic plots with given syntenic data.
		Usage: Run codes in RStudio
