## Requirements (Ubuntu 18.04.6 LTS)
- **Python**: >= 3.8
- **Biopython**: >= 1.85
- **Procodon**: Python module in <01-Procodon>

## Species for test
1. **Archaea**: *Methanocaldococcus vulcanius*, *Methanobacterium formicicum*, *Promethearchaeum syntrophicum*
2. **Bacteria**: *Thiobacillus sedimenti*, *Escherichia coli*, *Streptomyces collinus*
3. **Eukaryote**: *Saccharomyces cerevisia*, *Physcomitrium patens*, *Homo sapiens*

## Scripts
 - **01-Generation_of_Codec_test_data.py**: 
<br>&emsp;&emsp;Usage: ```python 01-Generation_of_Codec_test_data.py```
<br>&emsp;&emsp;Output: binary test data for each species

 - **02-Generation_of_available_codecs.py**: 
<br>&emsp;&emsp;Usage: ```python 02-Generation_of_available_codecs.py```
<br>&emsp;&emsp;Output: codecs for test

 - **03-Summarizing_GC_bias.py**: 
<br>&emsp;&emsp;Usage: ```python 01-Summarizing_GC_bias.py <gc_bias_path> <out_gc_bias_csv_path>```
<br>&emsp;&emsp;Output: summarized GC bias data
<br>&emsp;&emsp;Note: this Python script is called by 04-Codec_test.py

 - **04-Codec_test.py**: 
<br>&emsp;&emsp;Usage: ```python 04-Codec_test.py <species_prefix> <species_type>```
<br>&emsp;&emsp;&emsp;&emsp;e.g. ```python 04-Codec_test.py E.coli pr```
<br>&emsp;&emsp;&emsp;&emsp;e.g. ```python 04-Codec_test.py S.cerevisiae eu```
<br>&emsp;&emsp;Output: GC bias result of each gene encoding each pattern of binary string

## Attached files
 - **Species + "/"**: directories containing all CDS and protein sequences of species for test, which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
<br>&emsp;&emsp;Note: The directories were named by abbreviations of each species. For genes with alternative splicing sites in *P. patens* and *H. sapiens*, the overlapped CDS regions with different frames were excluded for the alteration of one CDS might affect the others. Intersection of the remaining regions and incomplete codons were also removed.
- **Species + ".binary_test_data_out.json"**: files containing the binary data for codec test in this work, which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **Species + ".aa_frequency_dict.json"**: files containing the total number of each amino acid accross all proteins in each species.
 - **Species + ".codec_0.json"**: files containing the original codec of each species.
 - **Species + ".codon_frequency.csv"**: files containing the relative codon usage of each species.<br>&emsp;&emsp;Note: the codon usage of *E. coli*, *S. cerevisia* and *H. sapiens* were download from https://github.com/Edinburgh-Genome-Foundry/python_codon_tables/tree/master/python_codon_tables/codon_usage_data/tables. Other were calculated from all CDSs in each species.
 - **Species + ".codec_pass_qua_dict.json"**: files containing the first, first-quarter, middle, three-quarter and last GBS-pass codecs sorted by GBS of each species.
 - **Species + ".GBS.txt"**: files containing all potential available codecs ID and GBS of each species.
 - **Species + ".GC_bias/"**: directories containing GC bias results in this work, which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 