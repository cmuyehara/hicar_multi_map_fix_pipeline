import pysam
import argparse
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter R1 and R2 BAM entries to only retain aligments found in both files.")
    parser.add_argument('-r1', "--r1_in", help="Input BAM file for R1 reads. Must be sorted by read name.")
    parser.add_argument('-r2', "--r2_in", help="Input BAM file for R2 reads. Must be sorted by read name.")
    parser.add_argument('-o1', "--r1_out", help="Output BAM file for filtered R1 reads")
    parser.add_argument('-o2', "--r2_out", help="Output BAM file for filtered R2 reads") 
    args = parser.parse_args()
    return args

def sync_bams(r1_in, r2_in, r1_out, r2_out):
    # Open input BAMs
    bam1 = pysam.AlignmentFile(r1_in, "rb")
    bam2 = pysam.AlignmentFile(r2_in, "rb")

    # Open output BAMs, copying headers
    out1 = pysam.AlignmentFile(r1_out, "wb", header=bam1.header)
    out2 = pysam.AlignmentFile(r2_out, "wb", header=bam2.header)

    # Initialize iterators
    r1_iter = bam1.fetch(until_eof=True)
    r2_iter = bam2.fetch(until_eof=True)

    try:
        r1 = next(r1_iter)
        r2 = next(r2_iter)
    except StopIteration:
        print("One of the BAMs is empty.", file=sys.stderr)
        return

    total_pairs = 0
    matched = 0

    while True:
        name1 = r1.query_name
        name2 = r2.query_name

        if name1 == name2:
            # Perfect match: write both reads
            out1.write(r1)
            out2.write(r2)
            matched += 1
            try:
                r1 = next(r1_iter)
                r2 = next(r2_iter)
            except StopIteration:
                break
        elif name1 < name2:
            # Advance R1
            try:
                r1 = next(r1_iter)
            except StopIteration:
                break
        else:
            # Advance R2
            try:
                r2 = next(r2_iter)
            except StopIteration:
                break

        total_pairs += 1
        if total_pairs % 1000000 == 0:
            print(f"Processed {total_pairs:,} read names... matched {matched:,}", file=sys.stderr)

    # Clean up
    bam1.close()
    bam2.close()
    out1.close()
    out2.close()
    print(f"Done. Wrote {matched:,} synchronized read pairs.", file=sys.stderr)
    pass

def main():
    args = parse_arguments()
    sync_bams(args.r1_in, args.r2_in, args.r1_out, args.r2_out) 

if __name__ == "__main__":
    main()
