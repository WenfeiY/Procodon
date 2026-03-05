# Procodon
The code and test data of The Living "Hard Drive": Data Storage Programmed by Protein-coding Genes, Created by Wenfei Yu from Dailab@SIAT.

---

## Introduction
Storing information in vivo bridges digital and biological realms with encrypted, long-term data storage potential, but current methods mainly rely on exogenous, non-functional, possibly evolutionarily unstable sequences. Here, we present Pro-Codon, a in vivo data storage system embedding information into balanced protein codons to render genes both information-rich and functional. In *Saccharomyces cerevisiae*, Pro-Codon flexibly encodes multiple data streams across genes with high stability and no fitness cost. The system enables index-free retrieval via reference genome mapping and leverages codec and genomic diversity for robust encryption, resisting brute-force attacks for potentially over 10^15 years, providing a useful platform for safeguarding proprietary strains. Furthermore, by integrating synthetic genomics with music composition, we generated a 50 Gb library and “Yeast Fantasia”, where induced genomic rearrangements serve to enhance melodic diversity.  This work establishes a new paradigm for living data storage, seamlessly integrating digital information with life’s fundamental processes.

![Fig1 20260115](https://github.com/user-attachments/assets/8dfdffe3-e50a-4152-85bb-85bbd53296dd)


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
 - **Dataset-12-Genome_encoding_space**: <br>&emsp;&emsp;Sccript for computing available encoding space from annotations in GFF3 format. The exampled GFF3 file is too large to be uploaded, which can be downloaded from [figshare](https://doi.org/10.6084/m9.figshare.31055530).

## Note
- The detailed requirements and descriptions of each task are in the `README.md` in relevant directory.
- To save space, some of the required input files which are the same in different tasks are not provided in every directory. Please copy them from other directories.

## License
MIT
