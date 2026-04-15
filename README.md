# Procodon
The code and test data of "Programmable in vivo data storage in functional protein-coding sequences", Created by Wenfei Yu from Dailab@SIAT.

---

## Introduction
Current in vivo DNA storage predominantly relies on exogenous, non-functional sequences, limiting both storage capacity and long-term stability in living systems. Here we present Pro-Codon, a framework that directly writes digital information into functional protein-coding sequences by exploiting synonymous codon degeneracy. The system enables a large space of codon-binary mappings for flexible encoding, supports optional index-free data retrieval via native gene order, and provides intrinsic encryption. In Saccharomyces cerevisiae, we demonstrate stable storage of diverse data types, including text, images and music, across both essential and non-essential genes without compromising cellular fitness, even under extreme recoding conditions with sequence identities as low as ~60%. Pro-Codon further enables sequence-driven data generation, including the construction of a ~50 GB music library and the derivation of new compositions from Sc2.0 synthetic chromosomes through codec diversity and genome rearrangement. Notably, at the genome scale, we achieved functional replacement of an endogenous chromosomal arm with a neo-chromosome carrying encrypted information while preserving normal cellular performance. In the absence of a reference genome, decryption requires years of brute force effort even under highly constrained conditions, enabling robust protection of engineered strains. Together, Pro-Codon establishes a scalable paradigm for secure, function integrated data storage in living systems.

<img width="7937" height="9865" alt="Fig1 v12 1 20260325" src="https://github.com/user-attachments/assets/64df843f-bb89-45d1-aa09-5d6f6d609143" />


---
## System Requirements
- **Linux system**: Ubuntu 18.04.6 LTS
- **Python**: >= 3.8
- **Anaconda**: >=24.7.1

## Installation
Installation from source is provided here.

First, download the package via git clone
```
git clone https://github.com/WenfeiY/Procodon.git
cd Procodon
```
Create and activate anaconda environment
```
conda env create -f environment.yml
conda activate Procodon
```
We recommend to install [SPAdes](https://ablab.github.io/spades) from source and add SPAdes installation directory to the PATH variable.
```
wget https://github.com/ablab/spades/releases/download/v3.15.4/SPAdes-3.15.4-Linux.tar.gz
tar -xzf SPAdes-3.15.4-Linux.tar.gz
chmod 777 SPAdes-3.15.4-Linux/bin/*
```

## Usage
The directories with numbers include codes and test data for conducting all *in silico* progresses:
 - **01-Procodon**: <br>&emsp;&emsp;Python module named `Procodon`, provides all the basic functions of data storage via protein, including codec testing, sequence & file loading, data conversion, data encoding, info key & sequence decoding and other utilities.
 - **02-Codec_test**: <br>&emsp;&emsp;Scripts and test data for encoding binaries in all available regions of 9 species from archaea, bacteria and eukaryotes. Some datasets are too large to be uploaded, which can be downloaded from [figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **03-Data_encoding_and_decoding**: <br>&emsp;&emsp;Scripts for encoding information into given genes and decoding data from given CDSs.
 - **04-Analysis_of_orthologs**: <br>&emsp;&emsp;Python Script for extracting ortholog info of budding yeast from given Genbank file. R script and ortholog data for creating syntenic plot. The Genbank files are too large to be uploaded, which can be downloaded from [figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **05-Simulation_of_images_storage_with_bio-index**: <br>&emsp;&emsp;Scripts for extracting gene positions from budding yeast. Scripts, input images and bio-index info for encoding images into yeast genome with gene sorted as bio-index. The resulting sequences are too large to be uploaded, which can be downloaded from [figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **06-AA_corection_test**: <br>&emsp;&emsp;Scripts and data for testing reads correcting by AA, including NGS data generation, reads correction, data decoding and calculation of data recovery rates. Protein, gene flanking and intron sequences are too large to be uploaded, which can be downloaded from [figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **07-Estimation_of_decoding_time**: <br>&emsp;&emsp;Scripts and data for estimating the time cost of decoding and cracking complexity of different encryption methods.
 - **08-Music_encoding**: <br>&emsp;&emsp;Scripts and data for encoding single or complex music notes into given genes.
 - **09-Music_SCRaMbLE**: <br>&emsp;&emsp;Scripts and data for encoding "*Canon*" in synthetic chromosome and decoding rearranged music from SCRaMbLEd strains
 - **10-Generation_of_music_library**: <br>&emsp;&emsp;Scripts and data for extracting sequences and annotations of synthetic yeast genome and generating music library with all available codecs in budding yeast. Sequences & annotations of synthetic chromosomes and codecs are too large to be uploaded, which can be downloaded from [figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **11-Computing_melodic_difference**: <br>&emsp;&emsp;Scripts and data for quantifying distances between SCRaMbLEd melodies and the original one.
 - **12-Genome_encoding_space**: <br>&emsp;&emsp;Sccript for computing available encoding space from annotations in GFF3 format. The exampled GFF3 file is too large to be uploaded, which can be downloaded from [figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **13-Decoding_data_from_NGS**: <br>&emsp;&emsp;Sccript for decoding binary data from NGS data. The Fastq files are too large to be uploaded, which can be downloaded from [figshare](https://doi.org/10.6084/m9.figshare.31055530).

## Note
- The detailed requirements and descriptions of each task are in the `README.md` in relevant directory.
- To save space, some of the required input files which are the same in different tasks are not provided in every directory. Please copy them from other directories.

## License
MIT
