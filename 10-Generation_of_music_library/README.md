## Requirements (Spydes IDE 6.0.5):
- **Python**: >= 3.11.11
- **numpy**: >= 2.0.1
- **pandas**: >= 2.2.3
- **musicpy**: >= 7.11
- **Procodon**: Python module in <01-Procodon>

## Scripts
 - **01-Extraction_of_synthetic_sequences.py**: Extract synthetic CDS sequences from Genbank files and check.
<br>&emsp;&emsp;Usage: ```python 01-Extraction_of_synthetic_sequences.py```

 - **02-Generation_of_music_library.py**: Generate the music library.
<br>&emsp;&emsp;Usage: ```python 02-Generation_of_music_library.py```

## Attached files
 - **synthetic_chromosomes/**: Directory containing sequences and annotations of synthetic yeast genome, which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
<br>&emsp;&emsp;`*.gff`: annotations and sequences download from https://syntheticyeast.github.io/sc2-0/data/
<br>&emsp;&emsp;`*.gff3`: annotations extracted from `*gff`
<br>&emsp;&emsp;`*.fa`: sequences extracted from `*gff`
<br>&emsp;&emsp;`./genbank/`: genbank files generated from annotations and sequences extracted from `*.gff`
 - **S.cerevisiae.codec/**: Directory containing all potential available codecs of yeast, which can be downloaded from [Figshare](https://doi.org/10.6084/m9.figshare.31055530).
 - **syn_chr_CDS_seq_dict.json**: Sequences extracted by `01-Extraction_of_synthetic_sequences.py`, and modified as below:
<br>&emsp;&emsp;
    1. Dubious genes YCL007C, YCR018C-A, YDR278C, YDR445C, YFR052C-A, YHR180W-A, YIR023C-A, YIR030W-A, YML094C-A, YNL285W, YNL017C, YOL106W, YOL013W-A, YPR108W-A and YCR018C-A with sequences unlike CDS were removed.
    2. Uncharacterized YFL019C with sequences unlike CDS were removed.
    3. YPR169W-A was removed for being a replicate of YPR170W-B
    4. YOR072W-B and YPR159C-A with sequences unlike CDS were removed.
    5. Incomplete genes YAR053W, YDR241W, YDR366C, YDR396W, YGL182C, YJR030C, YKL118W, YLL064C and YLR256W were fixed.
    6. Separate features resulted from introns were joined.
