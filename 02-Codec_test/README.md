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
   *(Note: Alternative splicing regions excluded)*
3. `Homo sapiens`  
   *(Note: Alternative splicing regions excluded)*

---

# Directories
- Species-abbreviated directories contain all proteins/CDS files
- **Processing rules**:
  - Overlapped CDS regions with different frames excluded
  - Intersection regions and incomplete codons removed
  - Applies to `P. patens` and `H. sapiens`

---

# Scripts
### `01-Generation_of_Codec_test_data.py`
```bash
Usage: python 01-Generation_of_Codec_test_data.py
Output: Binary test data for all species
