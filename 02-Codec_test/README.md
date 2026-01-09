# Requirements (Ubuntu 18.04.6 LTS)
- **Python**: 3.7.10
- **Biopython**: 1.85
- **Procodon**: Module in this work

---

# Species for Codec Test
## Archaea
1. `Methanocaldococcus vulcanius`
2. `Methanobacterium formicicum`
3. `Promethearchaeum syntrophicum`

## Bacteria
1. `Thiobacillus sedimenti`
2. `Escherichia coli`
3. `Streptomyces collinus`

## Eukaryote
1. `Saccharomyces cerevisia`
2. `Physcomitrium patens`  
3. `Homo sapiens`  

---

# Directories
- Species-abbreviated directories contain all proteins/CDS files from relevant species
*(Note: For genes with alternative splicing sites in P. patens and H. sapiens, the overlapped CDS regions with different frames were excluded for alteration of one CDS might affect the others. Intersection of the remaining regions and incomplete codons were also removed.)*
---

# Files
- `Species + ".binary_test_data_out.json"`: files containing the binary data for codec test in this work.
- `Species + ".aa_frequency_dict.json"`: files containing the total number of each amino acid accross all proteins in each species.
- `Species + ".codec_0.json"`: files containing the original codec of each species.
- `Species + ".codon_frequency.csv"`: files containing the relative codon usage of each species.
(Note: the codon usage of E. coli, S. cerevisia and H. sapiens were download from https://github.com/Edinburgh-Genome-Foundry/python_codon_tables/tree/master/python_codon_tables/codon_usage_data/tables. Other were calculated from all CDSs in each species.)
- `Species + ".codec_pass_qua_dict.json"`: files containing the first, first-quarter, middle, three-quarter and last GBS-pass codecs sorted by GBS of each species.
- `Species + ".GBS.txt"`: files containing all potential available codecs ID and GBS of each species.
- `Species + ".GC_bias/"`: directories containing GC bias results in this work

# Scripts
- `01-Generation_of_Codec_test_data.py`
This script is used for generating binary dataset for codec test. The binary dataset encludes a series of random binary strings for all genes (or available encoding regions) covering their full lengths with 5 proportions of "1": 0% (all- "0"), 25%, 50%, 75% and 100% (all- "1"). The pattern all- "0" and all- "1" required a single string each, while twenty randomized binary sequences were generated for the intermediate proportions (25%, 50% and 75%).
Usage: 
```bash
python 01-Generation_of_Codec_test_data.py
```

- `02-Generation_of_available_codecs.py`
