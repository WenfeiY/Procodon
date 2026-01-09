import json
import sys
import os

def summary_gc_bias(
    json_file: str, 
    out_file: str
):
    '''
    Read GC bias data, Calculate average GC bias for each gene and proportion and write into Tab-Separated Values file

    Parameters:
        json_file: input GC bias json file path
        out_file: output GC bias tsv file path
    '''
    
    with open(json_file, 'r') as f:
        gc_bias_dict = json.loads(f.read())
    with open(out_file, 'w') as hd:
        hd.write('gene\t0\t0.25\t0.5\t0.75\t1\n')
        for gene in gc_bias_dict:
            bias_list = []
            for prop in gc_bias_dict[gene]:
                bias = round(sum(gc_bias_dict[gene][prop]) / len(gc_bias_dict[gene][prop]), 2)
                bias_list.append(str(bias))
            hd.write('\t'.join([gene] + bias_list) + '\n')


if __name__ == '__main__':
    json_file = sys.argv[1]
    out_file = sys.argv[2]

    summary_gc_bias(json_file, out_file)

    # Remove intermediate file to save hard drive space
    os.remove(json_file)
