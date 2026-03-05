"""
Compute rough encoding space from annotations in GFF3 format.
"""

import pandas as pd
import sys

def extract_info_from_attr(
    attributes: str, 
    key: str
):
    """
    Extract infos from attributes in GFF3 file
    
    Parameters:
        attributes: The last column separated by of a GFF3 line, str
        key: The key of attribute to be extracted, str
    
    Return:
        item_value: value of the provided key extracted from attributes, str
    """
    
    #  Split attributes
    for item in attributes.split(';'):
        if '=' in item:
            if key == item.split('=')[0]:
                return item.split('=')[1]
    return ''
 
def parse_parent_from_gff3(
    gff3_df, 
    ID_key: str = 'ID', 
    parent_key: str = 'Parent'
):
    """
    Parse GFF3 file and yield the gene-mRNA-CDS feature information in tree structure
    
    Parameters:
        gff3_df: DataFrame read from GFF3 file, pandas.DataFrame
        ID_key: The key of ID, universal identifier in GFF3 file, str
        parent_key: The key of parent, key
    
    Return:
        gene_dict: dictionary containing gene-mRNA-CDS feature information in tree structure with gene IDs as keys
    """    
    # Dictionary recording all gene info
    gene_dict = {}
    # Dictionary recording all mRNA info
    mrna_dict = {}
    # list recording all CDS info
    cds_list = []
    
    #  Traverse all lines in DataFrame, storing all gene, mRNA and CDS info
    for idx in gff3_df.index:
        if gff3_df.loc[idx, 'type'] == 'gene':
            gene_id = extract_info_from_attr(gff3_df.loc[idx, 'attributes'], ID_key)
            if gene_id:
                gene_dict[gene_id] = {'info': dict(gff3_df.loc[idx, ]), 'contain': []}
        if gff3_df.loc[idx, 'type'] == 'mRNA':
            fea_id = extract_info_from_attr(gff3_df.loc[idx, 'attributes'], ID_key)
            parent = extract_info_from_attr(gff3_df.loc[idx, 'attributes'], parent_key)
            if fea_id and parent:
                mrna_dict[fea_id] = {'info': dict(gff3_df.loc[idx, ]), 'Parent': parent, 'contain': []}
        if gff3_df.loc[idx, 'type'] == 'CDS':
            fea_id = extract_info_from_attr(gff3_df.loc[idx, 'attributes'], ID_key)
            parent = extract_info_from_attr(gff3_df.loc[idx, 'attributes'], parent_key)
            if fea_id and parent:
                cds_list.append({'info': dict(gff3_df.loc[idx, ]), 'Parent': parent})
    
    #  Link CDS info to mRNA parent
    for item in cds_list:
        if item['Parent'] in mrna_dict:
            mrna_dict[item['Parent']]['contain'].append(item)
    
    #  Link mRNA info to gene parent
    for item in mrna_dict:
        if mrna_dict[item]['Parent'] in gene_dict:
            gene_dict[mrna_dict[item]['Parent']]['contain'].append(mrna_dict[item])
    return gene_dict

#  Some annotation files include no mRNAs
def parse_parent_from_gff3_simple(
    gff3_df, 
    ID_key: str = 'ID', 
    parent_key: str = 'Parent'
):
    """
    Parse GFF3 file and yield the gene-CDS feature information in tree structure
    
    Parameters:
        gff3_df: DataFrame read from GFF3 file, pandas.DataFrame
        ID_key: The key of ID, universal identifier in GFF3 file, str
        parent_key: The key of parent, key
    
    Return:
        gene_dict: dictionary containing gene-CDS feature information in tree structure with gene IDs as keys
    """    
    # Dictionary recording all gene info
    gene_dict = {}
    # list recording all CDS info
    cds_list = []
    
    #  Traverse all lines in DataFrame, storing all gene and CDS info
    for idx in gff3_df.index:
        if gff3_df.loc[idx, 'type'] == 'gene':
            gene_id = extract_info_from_attr(gff3_df.loc[idx, 'attributes'], ID_key)
            if gene_id:
                gene_dict[gene_id] = {'info': dict(gff3_df.loc[idx, ]), 'contain': []}
        if gff3_df.loc[idx, 'type'] == 'CDS':
            fea_id = extract_info_from_attr(gff3_df.loc[idx, 'attributes'], ID_key)
            parent = extract_info_from_attr(gff3_df.loc[idx, 'attributes'], parent_key)
            if fea_id and parent:
                cds_list.append({'info': dict(gff3_df.loc[idx, ]), 'Parent': parent})
    
    #  Link CDS info to gene parent
    for item in cds_list:
        if item['Parent'] in gene_dict:
            gene_dict[item['Parent']]['contain'].append(item)
    return gene_dict

def extract_cds_region_and_phase(
    gene_info_dict: dict
):
    """
    Extract CDS regions and the "actual phase" calculated from the positions in genome

    Parameter:
        gene_info_dict: dictionary containing gene-mRNA-CDS feature information in tree structure with gene IDs as keys
    
    Returns:
    phase_dict: dictionary containing CDS regions [start, end] with their phases as keys
    actural_phase_dict: dictionary containing CDS regions [start, end] with their "actual phases" as keys
    """
    
    #  Initialize result dictionaries
    phase_dict = {0: [],
                  1: [],
                  2: []
                  }
    actural_phase_dict = {0: [],
                          1: [],
                          2: []
                          }
    
    #  Re-calculate the "actual phase" of all CDSs and record in dictionary
    for mrna in gene_info_dict['contain']:
        for cds in mrna['contain']:
            
            start = cds['info']['start'] - 1
            end = cds['info']['end']
            phase = int(cds['info']['phase'])
            if cds['info']['strand'] == '+':
                actural_phase = ((start % 3) + phase) % 3
            else:
                actural_phase = ((end % 3) - phase) % 3
            phase_dict[phase].append([start, end])
            actural_phase_dict[actural_phase].append([start, end])
    return phase_dict, actural_phase_dict

def extract_cds_region_and_phase_simple(
    gene_info_dict: dict
):
    """
    Extract CDS regions and the "actual phase" calculated from the positions in genome

    Parameter:
        gene_info_dict: dictionary containing gene-CDS feature information in tree structure with gene IDs as keys
    
    Returns:
    phase_dict: dictionary containing CDS regions [start, end] with their phases as keys
    actural_phase_dict: dictionary containing CDS regions [start, end] with their "actual phases" as keys
    """
    
    #  Initialize result dictionaries
    phase_dict = {0: [],
                  1: [],
                  2: []
                  }
    actural_phase_dict = {0: [],
                          1: [],
                          2: []
                          }
    
    #  Re-calculate the "actual phase" of all CDSs and record in dictionary
    for cds in gene_info_dict['contain']:
        start = cds['info']['start'] - 1
        end = cds['info']['end']
        phase = int(cds['info']['phase'])
        
        #  For annotations without mRNAs, NO NEED to consider the overlaps between different CDS from splicing
        actural_phase = ((start % 3) + phase) % 3
        phase_dict[phase].append([start, end])
        actural_phase_dict[actural_phase].append([start, end])
    return phase_dict, actural_phase_dict

def merge_ranges(
    ranges: list
):
    """
    Interval union

    Parameter:
        ranges: list containing intervals
    
    Return:
        merged: list after union
    """
    
    if not ranges:
        return []
    
    #  Sort the range list by positions
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    
    # Interval union
    merged = [sorted_ranges[0].copy()]
    for current in sorted_ranges[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1][1] = max(last[1], current[1])
        else:
            merged.append(current.copy())
    return merged
 
def extract_incomplete_codon_region(
    cds_region: list, 
    phase: int, 
    strand: str
):
    """
    Extract regions of incomplete codons caused by splicing from CDS region

    Parameters:
        cds_region: list containing start and end positions of CDS region
            e.g. [100, 300]
        phase: phase in integer (0, 1 or 2)
        strand: sequence strand ('+' or '-'), str
    Return:
        ex_region: regions of incomplete codons, list
    """
    
    ex_region = []
    
    if strand == '+':
        if phase == 0:
            ex_codon = (cds_region[1] - cds_region[0]) % 3
            if ex_codon == 0:
                return []
            region2 = [cds_region[1] - ex_codon, cds_region[1] + (3 - ex_codon)]
            return [region2]
        region1 = [cds_region[0] - (3 - phase), cds_region[0] + phase]
        region2 = [cds_region[1] - (cds_region[1] - (cds_region[0] + phase)) % 3, cds_region[1] + (3 - (cds_region[1] - (cds_region[0] + phase)) % 3)]
        if not region1[0] == region1[1]:
            ex_region.append(region1)
        if not region2[0] == region2[1]:
            ex_region.append(region2)
        return ex_region
    else:
        if phase == 0:
            ex_codon = (cds_region[1] - cds_region[0]) % 3
            if ex_codon == 0:
                return []
            region2 = [cds_region[0] - (3 - ex_codon), cds_region[0] + ex_codon]
            return [region2]
        region1 = [cds_region[1] - phase, cds_region[1] + (3 - phase)]
        region2 = [cds_region[0] - (3 - (cds_region[1] - cds_region[0] - phase) % 3), cds_region[0] + (cds_region[1] - cds_region[0] - phase) % 3]
        if not region1[0] == region1[1]:
            ex_region.append(region1)
        if not region2[0] == region2[1]:
            ex_region.append(region2)
        return ex_region
 
def extract_incomplete_codon_region_list(
    phase_dict: dict, 
    strand: str
):
    """
    Extract regions of incomplete codons caused by splicing from CDS region list

    Parameters
        phase_dict: dictionary containing CDS regions [start, end] with their phases as keys
        strand: sequence strand ('+' or '-'), str

    Return:
        ex_region_list: region list of incomplete codons, list
    """
    
    ex_region_list = []
    #  Parse all phases
    for phase in phase_dict:
        for cds_region in phase_dict[phase]:
            ex_region_list += extract_incomplete_codon_region(cds_region, phase, strand)
    ex_region_list = merge_ranges(ex_region_list)
    
    return ex_region_list
 
def merge_codon_region(
    phase_dict: dict
):
    """
    Do interval union of all phases

    Parameter:
        phase_dict: dictionary containing CDS regions [start, end] with their phases as keys

    Return:
        pro_phase_dict: processed dictionary with merged regions.
    """
    
    pro_phase_dict = {0: [],
                      1: [],
                      2: []
                      }
    for phase in phase_dict:
        pro_phase_dict[phase] = merge_ranges(phase_dict[phase])
    return pro_phase_dict
 
 
def subtract_range(
    current_range: list, 
    exclude_ranges: list
):
    """
    Exclude incomplete codons from CDS range

    Parameters:
        current_range: list contsining start and end positions of CDS
        exclude_ranges: list containing incomplete codon ranges
    
    Return:
        result: list containing complete codons
    """
    #  Get start and end
    start, end = current_range
    #  Interval union of incomplete codon regions
    merged_exclude = merge_ranges(exclude_ranges)
    
    result = []
    #  remove incomplete codon regions from CDS
    for ex in merged_exclude:
        ex_start, ex_end = ex
        if ex_end <= start:
            continue
        if ex_start >= end:
            break
        if ex_start > start:
            result.append([start, ex_start])
            start = ex_end
        else:
            start = max(start, ex_end)
        if start >= end:
            break
    if start < end:
        result.append([start, end])
    return result
 
def remove_diff_frame_overlap(
    pro_phase_dict: dict
):
    """
    Remove the overlapped regions of CDS regions in different "actual phases"

    Parameter:
        pro_phase_dict: processed dictionary with merged regions.
    
    Return:
        rm_overlap_phase_dict: dictionary containing processed regions with phases as keys
    """
    rm_overlap_phase_dict = {}
    for phase in pro_phase_dict:
        other_intervals = []
        for k in pro_phase_dict:
            if k != phase:
                other_intervals.extend(pro_phase_dict[k])
        merged_other = merge_ranges(other_intervals)
        new_intervals = []
        for interval in pro_phase_dict[phase]:
            subtracted = subtract_range(interval, merged_other)
            new_intervals.extend(subtracted)
        rm_overlap_phase_dict[phase] = new_intervals
    return rm_overlap_phase_dict
 
def remove_ex_region(
    rm_overlap_phase_dict: dict, 
    ex_region_list: list
):
    """
    Remove incomplete codons regions from CDS ranges

    Parameters:
        rm_overlap_phase_dict: dictionary containing processed regions (removing the overlapped regions of CDS regions in different "actual phases") with phases as keys
        ex_region_list: list containing incomplete codon regions
    
    Return:
        avail_phase_dict: dictionary with CDS regions which excluded incomplete codons
    """
    avail_phase_dict = {}
    for phase in rm_overlap_phase_dict:
        new_intervals = []
        for interval in rm_overlap_phase_dict[phase]:
            subtracted = subtract_range(interval, ex_region_list)
            new_intervals.extend(subtracted)
        avail_phase_dict[phase] = new_intervals
    return avail_phase_dict
 
def in_range(
    single_interval, 
    intervals
):
    """
    Determines whether a given single interval overlaps with a series of intervals.

    Parameters:
        single_interval: provided single inverval, list
            e.g. [100, 200]
        intervals: provided intervals, list
            e.g. [[50, 150], [250, 350], ...]
    
    Return:
        iv: overlaped interval from series of intervals
            if no overlap, return False
    """
    for iv in intervals:
        if single_interval[0] >= iv[0] and single_interval[1] <= iv[1]:
            return iv
    return False
 
def cut_incomplete_codon(
    cds_range: list, 
    containing_range: list, 
    phase: int, 
    strand: str
):
    """
    Remove incomplete codons from final CDS region

    Parameters:
        cds_range: list contsining start and end positions of CDS
        containing_range: list contsining start and end positions of incomplete codon
        phase: "actural phase" of CDS in integer (0, 1 or 2)
        strand: sequence strand ('+' or '-'), str
    
    Return:
        complete CDS region: CDS region removed incomplete codons, list
    """
    if strand == '+':
        phase_0_start = cds_range[0] + (phase + containing_range[0] - cds_range[0]) % 3
        phase_0_end = cds_range[1] - (cds_range[1] - phase_0_start) % 3
    else:
        phase_0_end = cds_range[1] - (phase + cds_range[1] - containing_range[1]) % 3
        phase_0_start = cds_range[0] + (phase_0_end - cds_range[0]) % 3
    if phase_0_end > phase_0_start:
        return [phase_0_start, phase_0_end]
    else:
        return []
 
def rm_rest_incomplete_codons(
    avail_cds_ranges: list, 
    merged_phase_dict: dict, 
    strand: str
):
    """
    Remove the rest incomplete codons in the available CDS regions

    Parameters:
        avail_cds_ranges: list contsining available CDS ranges
        merged_phase_dict: CDS ranges (interval union) with CDS phase as keys
        strand: sequence strand ('+' or '-'), str
    
    Return:
        cut_avail_cds_ranges: final available CDS ranges, list
    """
    cut_avail_cds_ranges = []
    for cds_range in avail_cds_ranges:
        for phase in merged_phase_dict:
            containing_range = in_range(cds_range, merged_phase_dict[phase])
            if containing_range:
                cut_cds_range = cut_incomplete_codon(cds_range, containing_range, phase, strand)
                if cut_cds_range:
                    cut_avail_cds_ranges.append(cut_cds_range)
    return cut_avail_cds_ranges
     
def determine_avail_CDS_ranges(
    gene_info_dict: dict, 
    mode: str
):
    """
    Summarize functions above, generate available CDS regions for encoding

    Parameter:
        gene_info_dict: dictionary containing gene-CDS feature information in tree structure with gene IDs as keys
        mode: "complex" for gene-mRNA-CDS, "simple" for gene-CDS

    Return:
        cut_avail_cds_ranges: final available CDS ranges, list
    """
    strand = gene_info_dict['info']['strand']
    if mode == 'simple':
        phase_dict, actural_phase_dict = extract_cds_region_and_phase_simple(gene_info_dict)
    else:
        phase_dict, actural_phase_dict = extract_cds_region_and_phase(gene_info_dict)
    pro_phase_dict = merge_codon_region(actural_phase_dict)
    merged_phase_dict = merge_codon_region(phase_dict)
    rm_overlap_phase_dict = remove_diff_frame_overlap(pro_phase_dict)
    ex_region_list = extract_incomplete_codon_region_list(phase_dict, strand)
    avail_phase_dict = remove_ex_region(rm_overlap_phase_dict, ex_region_list)
    avail_cds_ranges = []
    for phase in avail_phase_dict:
        avail_cds_ranges.extend(avail_phase_dict[phase])
    avail_cds_ranges = [x for x in avail_cds_ranges if x != []]
    cut_avail_cds_ranges = rm_rest_incomplete_codons(avail_cds_ranges, merged_phase_dict, strand)
    cut_avail_cds_ranges = sorted(cut_avail_cds_ranges, key=lambda x: x[0])
    if strand == '-':
        cut_avail_cds_ranges = cut_avail_cds_ranges[::-1]
    return cut_avail_cds_ranges

def extract_recode_avail_from_gff(
    gff3_df, 
    out_region_file: str, 
    mode: str
):
    """
    Summarize functions above, generate available CDS regions from GFF3 DataFrame and calculate encoding space
    
    Parameters:
        gff3_df: annotation dataframe read from GFF3 file, pandas.DataFrame
        out_region_file: output file path of available encoding regions
        mode: "complex" for gene-mRNA-CDS, "simple" for gene-CDS
    
    Output:
        Write available encoding regions
    
    Return:
        total_space: calculated rough encoding space (bits)
    """
    ref_col = gff3_df['ref']
    records = ref_col.unique()
    
    total_space = 0
    
    with open(out_region_file, 'w') as region_hd:
        for record_id in records:
            sub_df = gff3_df.loc[gff3_df['ref'] == record_id, ]
            if mode == 'simple':
                gene_dict = parse_parent_from_gff3_simple(sub_df)
            else:
                gene_dict = parse_parent_from_gff3(sub_df)
            for gene in gene_dict:
                write_regions = []
                avail_cds_ranges = determine_avail_CDS_ranges(gene_dict[gene], mode)
                if avail_cds_ranges:
                    region_hd.write(gene + ':')
                    for region in avail_cds_ranges:
                        write_regions.append(str(region[0]) + ',' + str(region[1]))
                        total_space += (region[1] - region[0]) / 3
                    region_hd.write(';'.join(write_regions) + '\n')
    return total_space

def main(
    in_gff3_file: str, 
    out_region_file: str
):
    """
    Main program, calculate available encoding space from GFF3 file

    Parameters:
        in_gff3_file: file path of GFF3, str
        out_region_file: output file path of available encoding regions
    
    Output:
        Write available encoding regions
    
    Return:
        total_space: calculated rough encoding space (bits)
    """
    gff3_df = pd.read_csv(in_gff3_file, sep = '\t', header = None, comment = '#')
    gff3_df.columns = ['ref', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    
    type_col = list(gff3_df['type'])
    if 'mRNA' in type_col:
        mode = 'complete'
    else:
        mode = 'simple'
    
    total_space = extract_recode_avail_from_gff(gff3_df, out_region_file, mode)
    
    return total_space


if __name__ == '__main__':
    #  Get file paths
    in_gff3_file = sys.argv[1]
    out_region_file = sys.argv[2]
    #  Run main program
    total_space = main(in_gff3_file, out_region_file)
    print(total_space)
