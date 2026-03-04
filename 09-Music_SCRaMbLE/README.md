## Requirements (Spydes IDE 6.0.5):
- **Python**: >= 3.11.11
- **numpy**: >= 2.0.1
- **pandas**: >= 2.2.3
- **musicpy**: >= 7.11
- **Procodon**: Python module in <01-Procodon>

## Scripts
 - **01-Encoding_Canon_in_synIXR_with_extended_bin.py**: Read notes, extend the binary string (150-bit for 1 note) and encode complete notes into synIXR.
<br>&emsp;&emsp;Usage: ```python 01-Encoding_Canon_in_synIXR_with_extended_bin.py```

 - **02-Generation_of_SCRaMbLEd_music_with_extended_bin.py**: Read structures of LoxPsym units and decode rearranged music from SCRaMbLEd strains.
<br>&emsp;&emsp;Usage: ```python 02-Generation_of_SCRaMbLEd_music_with_extended_bin.py```

## Attached files
 - **S.cerevisiae/**: Directory containing all CDS and protein sequences of yeast, which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **synIXR.SCRaMbLE.structure/**: Directory containing structures of LoxPsym units of SCRaMbLEd strains.
<br>&emsp;&emsp;`Strain + ".LU.txt"`: Each line contains a number representing LoxPsym unit (1~43). Negative numbers represent reversed LoxPsym units.
 - **synIXR.SCRaMbLE.gene.sort/**: Directory containing gene order of SCRaMbLEd strains.
<br>&emsp;&emsp;`Strain + ".lu_gene.txt"`: Each line contains a gene name. "-" represent reversed gene.
 - **S.cerevisiae.codon_frequency.csv**: The relative codon usage of yeast.
 - **Canon.synIXR.150_bit_per_note.txt**: Note infos of melodic fragment from Canon.
