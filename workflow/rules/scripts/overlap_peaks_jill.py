import pandas as pd
import bioframe as bf
import argparse
from pathlib import Path

def filter_50percent_overlap(overlaps_df):
    if overlaps_df.empty:
        return overlaps_df
    
    overlaps_df['len1'] = overlaps_df['end'] - overlaps_df['start']
    overlaps_df['len2'] = overlaps_df['end_'] - overlaps_df['start_']
    overlaps_df['overlap_len'] = overlaps_df['overlap_end'] - overlaps_df['overlap_start']
    overlaps_df['frac1'] = overlaps_df['overlap_len'] / overlaps_df['len1']
    overlaps_df['frac2'] = overlaps_df['overlap_len'] / overlaps_df['len2']
    
    mask = (overlaps_df['frac1'] >= 0.5) | (overlaps_df['frac2'] >= 0.5)
    return overlaps_df[mask]

def main(args):
    all_fragments = args.all_fragments
    subset1 = args.subset1
    subset2 = args.subset2
    output = args.output

    peak_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
                 'signalValue', 'pValue', 'qValue', 'peak']

    all_fragments = pd.read_csv(all_fragments, sep='\t', names=peak_cols, header=None)
    subset1 = pd.read_csv(subset1, sep='\t', names=peak_cols, header=None)
    subset2 = pd.read_csv(subset2, sep='\t', names=peak_cols, header=None)
    
    overlaps1 = bf.overlap(all_fragments, subset1, return_overlap=True, how='inner')
    filtered_overlaps1 = filter_50percent_overlap(overlaps1)

    if not filtered_overlaps1.empty:
        passed_peaks = all_fragments.iloc[filtered_overlaps1.index.get_level_values(0).unique()]
    else:
        raise ValueError("No peaks that overlap subset 1 by ≥50%")

    overlaps2 = bf.overlap(passed_peaks.reset_index(drop=True), subset2, return_overlap=True, how='inner')
    filtered_overlaps2 = filter_50percent_overlap(overlaps2)
    
    if not filtered_overlaps2.empty:
        final_peaks = passed_peaks.reset_index(drop=True).iloc[
            filtered_overlaps2.index.get_level_values(0).unique()
        ]
    else:
        raise ValueError("No peaks that overlap both subsets by ≥50%")
   
    final_peaks = final_peaks.drop_duplicates().sort_values(['chrom', 'start'])
    final_peaks.to_csv(output, sep='\t', header=False, index=False)

    print(f"Found {len(final_peaks)} peaks that overlap both subsets by ≥50%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("all_fragments", type=Path, help="Path to all fragments narrowPeak file")
    parser.add_argument("subset1", type=Path, help="Path to subset 1 narrowPeak file")
    parser.add_argument("subset2", type=Path, help="Path to subset 2 narrowPeak file")
    parser.add_argument("output", type=Path, help="Path to output overlap peaks file")
    args = parser.parse_args()
    main(args)
