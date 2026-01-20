# Procodon
The code and test data of The Living "Hard Drive": Data Storage Programmed by Protein-coding Genes, Created by Wenfei Yu from Dailab@SIAT.

---

## Introduction
Storing information in vivo bridges digital and biological realms with encrypted, long-term data storage potential, but current methods mainly rely on exogenous, non-functional, possibly evolutionarily unstable sequences. Here, we present Pro-Codon, a in vivo data storage system embedding information into balanced protein codons to render genes both information-rich and functional. In Saccharomyces cerevisiae, Pro-Codon flexibly encodes multiple data streams across genes with high stability and no fitness cost. The system enables index-free retrieval via reference genome mapping and leverages codec and genomic diversity for robust encryption, resisting brute-force attacks for potentially over 1015 years, providing a useful platform for safeguarding proprietary strains. Furthermore, by integrating synthetic genomics with music composition, we generated a 50 Gb library and “Yeast Fantasia”, where induced genomic rearrangements serve to enhance melodic diversity.  This work establishes a new paradigm for living data storage, seamlessly integrating digital information with life’s fundamental processes.

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
 - **01-Procodon**: Python module named `Procodon`, provides all the basic functions of data storage via protein, including codec testing, sequence & file loading, data conversion, data encoding, info key & sequence decoding and other utilities.
 - **02-Codec_test**: Scripts and test data for encoding binaries in all available regions of 9 species from archaea, bacteria and eukaryotes. Some datasets are too large to upload, which can be downloaded from [figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **03-Data_encoding_and_decoding**: Scripts for encoding information into given genes and decoding data from given CDSs.
