from Procodon.seqloader import Seqloader
from Procodon.codec.testdata import Testdata

species_prefix_list = ['M.vulcanius', 'M.formicicum', 'P.syntrophicum', 'T.sedimenti', 'E.coli', 'S.collinus', 'S.cerevisiae', 'P.patens', 'H.sapiens']

for species_prefix in species_prefix_list:

    #  Load CDS and protein sequences
    seqload = Seqloader(seq_dir = species_prefix)
    
    #  Generate binary dataset
    binary_number = 20  #  For each gene and proportion (except 0 and 1), generate 20 random binary strings.
    proportion_1_list = [0, 0.25, 0.5, 0.75, 1]  #  "1" proportion
    
    datatest = Testdata(seqload.gene_seq_dict)
    datatest.prep_test_data(binary_number, proportion_1_list, species_prefix + '.binary_test_data_out.json')  #  Generate and write test dataset