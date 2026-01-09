Requirements (Spydes IDE 6.0.5):
    python (3.11.11)
    numpy (2.0.1)
	pandas (2.2.3)
    musicpy (7.11)
    Procodon (Module in this work)

Directories:
    synthetic_chromosomes: sequences and annotations of synthetic yeast genome
        *.gff: annotations and sequences download from https://syntheticyeast.github.io/sc2-0/data/
        *.gff3: annotations extracted from *gff
        *.fa: sequences extracted from *gff
        ./genbank/: genbank files generated from annotations and sequences extracted from *.gff
    S.cerevisiae.codec: All potential available codecs of yeast.

Files:
    syn_chr_CDS_seq_dict.json: sequences extracted by 01-Extraction_of_synthetic_sequences.py, and modified as below:
        1. Dubious genes YCL007C, YCR018C-A, YDR278C, YDR445C, YFR052C-A, YHR180W-A, YIR023C-A, YIR030W-A, YML094C-A, YNL285W, YNL017C, YOL106W, YOL013W-A, YPR108W-A and YCR018C-A with sequences unlike CDS were removed.
        2. Uncharacterized YFL019C with sequences unlike CDS were removed.
        3. YPR169W-A was removed for being a replicate of YPR170W-B
        4. YOR072W-B and YPR159C-A with sequences unlike CDS were removed.
        5. Incomplete genes YAR053W, YDR241W, YDR366C, YDR396W, YGL182C, YJR030C, YKL118W, YLL064C and YLR256W were fixed.
        6. Separate features resulted from introns were joined.

Script:
    01-Extraction_of_synthetic_sequences.py: extracte synthetic CDS sequences from genbank files and check
        Usage: 01-Extraction_of_synthetic_sequences.py
    
    02-Generation_of_music_library.py: generate the music library
        Usage: 02-Generation_of_music_library.py
