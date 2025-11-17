import argparse
import pysam
import pandas as pd

restriction_sites = {
    'MseI' : 'TTAA',
    'NlaIII' : 'CATG',
}

def parse_args():
    parser = argparse.ArgumentParser(description="Filter and check multimappers in a BAM file")
    parser.add_argument("-i", '--input', required = True, help="Input BAM file sorted by read name")
    parser.add_argument("-o", "--output", required = True, help="Output summary file")
    args = parser.parse_args()
    return args

def check_clipped(cur_read):
    """
    Check if a read has soft or hard clipping in its CIGAR string.
    Input:
     - read: pysam read object
    Output:
     - "hard" or "soft" if the read has soft or hard clipping, False otherwise
    Notes:
     - The output of cigartuples is a list of (operation, length) tuples
       where operation is an integer code:
         - 4 : Soft clipping (S)
         - 5 : Hard clipping (H)
     - Soft clipping (S) may only be at the 5' end of the read
     - Full scan of CIGAR string not needed
     """
    if cur_read.cigartuples is None:
        return False
    clip_dict = {4: 'S', 5: 'H'}
    clip_state = [ cur_read.cigartuples[0][0], cur_read.cigartuples[-1][0] ]
    clip_state = [ clip_dict[i] if i in clip_dict.keys() else 'N' for i in clip_state ]
    clip_state = ''.join(clip_state)
    return clip_state

def get_mapping_type(n, cur_read_name, next_read_name):
    """
    Determine if the current read is single-mapped or multi-mapped.
    Input:
     - n: number of times the current read name has been seen
     - cur_read_name: name of the current read
     - next_read_name: name of the next read
    Output:
     - "single_mapped" or "multi_mapped"
    """
    if n == 1 and cur_read_name != next_read_name:
        read_type = 'single_mapped'
        n = 1
    elif n > 1 or cur_read_name == next_read_name:
        read_type = 'multi_mapped'
    if cur_read_name != next_read_name:
        n = 1
    else: 
        n += 1
    return(read_type, n)

def summary_stats_dict():
    """ Create a nested dictionary to hold summary statistics."""
    summary_stats = {
        'all' : {
            'total_reads': 0,
            'mapped_reads': 0,
        },
        'single_mapped': {
            'n_reads': 0,
            'with_zero_mseI': 0,
            'with_single_mseI': 0,
            'with_multiple_mseI': 0,
            'with_zero_nlaIII': 0,
            'with_single_nlaIII': 0,
            'with_multiple_nlaIII': 0,
        },
        'multi_mapped': {
            'n_reads': 0,
            'with_zero_mseI': 0,
            'with_single_mseI': 0,
            'with_multiple_mseI': 0,
            'with_zero_nlaIII': 0,
            'with_single_nlaIII': 0,
            'with_multiple_nlaIII': 0,
        },
    }
    return(summary_stats)

def summarize_read_info(read, summary_stats):
    """ Update summary statistics for a given read type."""
    summary_stats['n_reads'] += 1
    clipped = check_clipped(read)
    if clipped:
        if clipped in summary_stats.keys():
            summary_stats[clipped] += 1
        else:
            summary_stats[clipped] = 1
    seq = read.query_sequence
    mseI_count = seq.count(restriction_sites['MseI'])
    nlaIII_count = seq.count(restriction_sites['NlaIII'])
    if mseI_count == 0:
        summary_stats['with_zero_mseI'] += 1
    elif mseI_count == 1:
        summary_stats['with_single_mseI'] += 1
    elif mseI_count > 1:
        summary_stats['with_multiple_mseI'] += 1
    if nlaIII_count == 0:
        summary_stats['with_zero_nlaIII'] += 1
    elif nlaIII_count == 1:
        summary_stats['with_single_nlaIII'] += 1
    elif nlaIII_count > 1:
        summary_stats['with_multiple_nlaIII'] += 1
    return summary_stats

def get_stats_on_bam(input_bam):
    """ Get summary statistics on a BAM file sorted by read name.
    Input:
     - input_bam: path to BAM file sorted by read name
    Output:
     - summary_stats: nested dictionary with summary statistics
    Note:
     - The BAM file must be sorted by read name.
        - The function assumes that all reads with the same name are contiguous in the BAM file.
        - The function counts total reads, mapped reads, single-mapped reads, multi-mapped reads,
          and various statistics on clipping and restriction sites for single- and multi-mapped reads.
        - The function processes the BAM file in a single pass.
    """
    first_read = True
    summary_stats = summary_stats_dict()
    with pysam.AlignmentFile(input_bam, "rb") as bam:
        for next_read in bam.fetch(until_eof = True):
            if first_read:
                first_read = False
                cur_read = next_read
                n = 1
                continue
            summary_stats['all']['total_reads'] += 1
            if next_read.is_unmapped:
                 continue
            else:
                summary_stats['all']['mapped_reads'] += 1
            # Determine if single or multi-mapped
            read_type, n = get_mapping_type(n, cur_read.query_name, next_read.query_name)
            summary_stats[read_type] = summarize_read_info(cur_read, summary_stats[read_type])
            cur_read = next_read
    return(summary_stats)

def main():
    args = parse_args()
    summary_dat = get_stats_on_bam(args.input)
    summary_df = { i : pd.DataFrame.from_dict( summary_dat[i], orient = 'index', columns = ['count'])
        for i in summary_dat.keys() }
    summary_df = pd.concat(summary_df, axis = 0)
    summary_df.index.names = ['mapping_status', 'metric']
    summary_df = summary_df.reset_index(level = 1)
    summary_df.to_csv(args.output, sep='\t', header=True, index=True)

if __name__ == "__main__":
    main()


        #'multi_mapped_reads': 0,
        #'multi_mapped_with_clipping': 0,
        #'multi_mapped_with_single_mseI': 0,
        #'multi_mapped_with_multiple_mseI': 0,
        #'multi_mapped_with_single_nlaIII': 0,
        #'multi_mapped_with_multiple_nlaIII': 0,
    #summary_df = pd.DataFrame.from_dict(summary_dat, orient='index', columns=['count'])
    #test_restruct = { i : pd.DataFrame.from_dict(test_dict[i], orient = 'index', columns = ['count'])
    #    for i in test_dict.keys()
    #}
    #test_restruct = pd.concat(test_restruct, axis = 0)
    #test_restruct.index.names = ['mapping_status', 'metric']
    #test_restruct = test_restruct.reset_index( level = 1)

#    if cur_read.cigartuples[0][0] == 4:
#        return 'soft_5prime_clipped'
#    elif cur_read.cigartuples[-1][0] == 4:
#        return 'soft_3prime_clipped'
#    elif cur_read.cigartuples[0][0] == 5:
#        return 'hard_5prime_clipped'
#    elif cur_read.cigartuples[-1][0] == 5:
#        return 'hard_3prime_clipped'
#    else:
#        return False


#def check_clipped(cur_read):
#    """
#    Check if a read has soft or hard clipping in its CIGAR string.
#    Input:
#     - read: pysam read object
#    Output:
#     - "hard" or "soft" if the read has soft or hard clipping, False otherwise
#     """
#    if cur_read.cigartuples is None:
#        return False
#    for (operation, length) in cur_read.cigartuples:
#        if operation == 4: return 'soft_clipped'
#        elif operation == 5: return 'hard_clipped'
#    return False
#
#    summary_stats = {
#        'all' : {
#            'total_reads': 0,
#            'mapped_reads': 0,
#        },
#        'single_mapped': {
#            'n_reads': 0,
#            'without_clipping': 0,
#            'hard_5prime_clipped': 0,
#            'hard_3prime_clipped': 0,
#            'soft_5prime_clipped': 0,
#            'soft_3prime_clipped': 0,
#            'with_zero_mseI': 0,
#            'with_single_mseI': 0,
#            'with_multiple_mseI': 0,
#            'with_zero_nlaIII': 0,
#            'with_single_nlaIII': 0,
#            'with_multiple_nlaIII': 0,
#        },
#        'multi_mapped': {
#            'n_reads': 0,
#            'without_clipping': 0,
#            'hard_5prime_clipped': 0,
#            'hard_3prime_clipped': 0,
#            'soft_5prime_clipped': 0,
#            'soft_3prime_clipped': 0,
#            'with_zero_mseI': 0,
#            'with_single_mseI': 0,
#            'with_multiple_mseI': 0,
#            'with_zero_nlaIII': 0,
#            'with_single_nlaIII': 0,
#            'with_multiple_nlaIII': 0,
#        }
#    }