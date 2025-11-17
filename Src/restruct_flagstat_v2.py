import pandas as pd

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Restructure flagstat output to a table")
    parser.add_argument("-i", '--input')
    parser.add_argument("-o", "--output")
    parser.add_argument("-s", "--step", default = '', help="Processing step (e.g., mapq, dedup, etc.)")
    args = parser.parse_args()
    return(args)


def restructure_flagstat(input_file, output_file, step = ''):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    data = []
    for line in lines:
        parts = line.strip().split(' + ')
        if len(parts) == 2:
            count1, rest = parts
            count2, description = rest.split(' ', 1)
            data.append((int(count1), int(count2), description.strip()))

    df = pd.DataFrame(data, columns=['Count1', 'Count2', 'Description'])
    df['Total'] = df.apply(lambda row: row['Count1'] + row['Count2'], axis=1)
    df['Count1_Percent'] = df[['Count1', 'Total']].apply(lambda x: f"{((x['Count1']/x['Total'])*100):.2f}%" if x['Total'] > 0 else "0.00%", axis=1)
    df['Count2_Percent'] = df[['Count2', 'Total']].apply(lambda x: f"{((x['Count2']/x['Total'])*100):.2f}%" if x['Total'] > 0 else "0.00%", axis=1)
    df['Step'] = step
    col_order = ['Step', 'Description', 'Count1', 'Count1_Percent', 'Count2', 'Count2_Percent', 'Total']
    df = df[col_order]
    df.to_csv(output_file, sep='\t', index=False)

def main():
    args = parse_args()
    restructure_flagstat(args.input, args.output, args.step)

if __name__ == "__main__":
    main()