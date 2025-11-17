# Overview 
# The purpose of this script is to filter multi-mapped bam entries and output the alignment position that is most proximal to the 5' end of the read.
# The script takes a BAM file as input and produces a filtered BAM file as output.
# The script uses the pysam library to read and write BAM files, and it processes each read to determine the best alignment based on the 5' end proximity.

# Planning Steps:
# 1. Import necessary libraries (pysam, argparse).
# 2. Define a function to process the BAM file.
# 3. Open the input BAM file using pysam.
# 4. Create an output BAM file to write the filtered results.
# 5. Iterate through each read in the input BAM file.
# 6. For each read, check if it is a multi-mapper (i.e., has multiple alignments).
#    - if only one alignment exists, write it.  
#    - multi-mappers are identified by the presence of multiple alignments for the same read name.
#    - For each multi-mapper, continue iterating through the reads until a read with a different name is encountered.
# 7. For each read with the same name, determine the alignment that is closest to the 5' end of the read.
#    - Check if the read is hard-clipped; if add the clipped bases to the alignment position.
#    - Check the strand of the read (forward or reverse).
#    - For forward strand reads, the 5' end is at the start of the alignment.
#    - For reverse strand reads, the 5' end is at the end of the
#    - The CIGAR string will be used to determine the actual alignment position.
#    - The CIGAR string is dependent on the alignment orientation (strand), so we must account for that.
# 8. Write the selected alignment to the output BAM file.
# 9. Close both the input and output BAM files

import pysam
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter BAM entries to retain the alignment closest to the 5' end of the read.")
    parser.add_argument('-i', "--input_bam", help="Input BAM file with multi-mapped reads. Must be sorted by read name.")
    parser.add_argument('-o', "--output_bam", help="Output BAM file with filtered reads")
    return parser.parse_args()

def get_and_adjust_map_pos(read):
    """
    Adjust the alignment position based on hard clipping in the CIGAR string.
    Input:
     - read: pysam AlignedSegment object
    Output:
     - adjusted map_start and map_end as integers
    """
    cigartuples = read.cigartuples
    map_start = read.query_alignment_start
    map_end = read.query_alignment_end
    is_reverse = read.is_reverse 
    read_length = read.infer_read_length()

    # Adjust positions based on hard clipping
    if cigartuples is None:
        return (map_start, map_end)
    elif cigartuples[0][0] == 5:  # if hard clipping at the start of the read, adjust map_start
        map_start += cigartuples[0][1]
        map_end += cigartuples[0][1]
    else:
        pass
    # Adjust positions for reverse strand
    if is_reverse:
        map_start_correct = read_length - map_end
        map_end_correct = read_length - map_start
        return(map_start_correct, map_end_correct)
    else:
        return(map_start, map_end)

def select_best_alignment(multi_mappers):
    """
    Select the best alignment from a list of multi-mapped reads based on proximity to the 5' end.
    
    Input:
     - multi_mappers: list of pysam AlignedSegment objects with the same read name
    Output:
     - best_read: pysam AlignedSegment object that is closest to the 5' end of the read 
    """
    best_read = None
    best_position = None
    for read in multi_mappers:
        map_pos = get_and_adjust_map_pos(read)
        if best_position is None or map_pos[0] < best_position[0]:
            best_position = map_pos
            best_read = read
        elif map_pos[0] == best_position[0]:
            # If positions are the same, we can choose to keep the longer alignment or the first one encountered
            if read.query_alignment_length > best_read.query_alignment_length:
                best_read = read
            else:   
                continue
        else:
            continue
    return best_read

def filter_multi_mappers(input_bam, output_bam):
    """Filter multi-mapped reads in a BAM file to retain the alignment closest to the 5' end of the read."""
    # Open the input BAM file
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    # Create the output BAM file
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    bam_iter = bam_in.fetch(until_eof=True)
    cur_read = None
    for next_read in bam_iter:
        if cur_read is None or cur_read.is_unmapped:
            cur_read = next_read
            continue
        multi_mappers = [cur_read]
        n = 1 # Counter for the number of reads with the same name
        while True:
            if cur_read.query_name == next_read.query_name:
                n += 1
                multi_mappers.append(next_read)
                try:
                    next_read = next(bam_iter)
                except StopIteration:
                    next_read = None
                    break
            elif next_read is None or cur_read.query_name != next_read.query_name:
                break
        # Process the previous set of multi-mappers
        # Only process if there are multiple alignments
        if len(multi_mappers) == 1:
            bam_out.write(multi_mappers[0])
        elif len(multi_mappers) > 1:
            best_read = select_best_alignment(multi_mappers)
            bam_out.write(best_read)
        else:
            pass
        # Reset for the new read name
        cur_read = next_read

    # close the BAM files
    bam_in.close()
    bam_out.close()

def main():
    args = parse_arguments()
    filter_multi_mappers(args.input_bam, args.output_bam)

if __name__ == "__main__":
    main()

    #map_start_correct = map_start
    #map_end_correct = map_end
    #return (str(map_start_correct), str(map_end_correct))


#def filter_multi_mappers(input_bam, output_bam):
#    # Open the input BAM file
#    bam_in = pysam.AlignmentFile(input_bam, "rb")
#    # Create the output BAM file
#    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)
#
#    #current_read_name = None
#    #multi_mappers = []
#
#    bam_iter = bam_in.fetch(until_eof=True)
#    cur_read = next(bam_iter)
#    for next_read in bam_iter:
#        if cur_read.is_unmapped:
#            cur_read = next_read
#            continue
#        multi_mappers = [cur_read]
#        n = 1 # Counter for the number of reads with the same name
#        while True:
#            if cur_read.query_name == next_read.query_name:
#                n += 1
#                multi_mappers.append(next_read)
#                next_read = next(bam_iter)
#            elif cur_read.query_name != next_read.query_name:
#                break
#        # Process the previous set of multi-mappers
#        # Only process if there are multiple alignments
#        if len(multi_mappers) == 1:
#            bam_out.write(multi_mappers[0])
#        elif len(multi_mappers) > 1:
#            best_read = select_best_alignment(multi_mappers)
#            bam_out.write(best_read)
#        else:
#            pass
#        # Reset for the new read name
#        cur_read = next_read
#
#    # close the BAM files
#    bam_in.close()
#    bam_out.close()









#    for read in multi_mappers:
#        if read.is_reverse:
#            # For reverse strand, 5' end is at the end of the alignment
#            alignment_end = read.reference_end
#            # Adjust for hard clipping at the end
#            if read.cigartuples:
#                if read.cigartuples[-1][0] == 5:  # Hard clip
#                    alignment_end += read.cigartuples[-1][1]
#            position = alignment_end
#        else:
#            # For forward strand, 5' end is at the start of the alignment
#            alignment_start = read.reference_start
#            # Adjust for hard clipping at the start
#            if read.cigartuples:
#                if read.cigartuples[0][0] == 5:  # Hard clip
#                    alignment_start -= read.cigartuples[0][1]
#            position = alignment_start


        #if len(multi_mappers) > 1:
        #    best_read = select_best_alignment(multi_mappers)
        #    bam_out.write(best_read)
        #elif len(multi_mappers) == 1:
        #    # If there's only one alignment, we write it as per the requirement
        #    pass
