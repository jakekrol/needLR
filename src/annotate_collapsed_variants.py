#!/usr/bin/env python3
import argparse
from cyvcf2 import VCF
import os
import pandas as pd
from jkbiolib.variant.vcf import count_alt_samples

parser = argparse.ArgumentParser(description="Annotate collapsed needLR variants")
parser.add_argument("--input_vcf", required=True, help="Path to the input VCF file")
parser.add_argument("--collapsed_vcf", required=True, help="Path to the collapsed VCF file")
parser.add_argument("--merged_vcf", required=True, help="Path to the merged VCF file")
parser.add_argument("--out", required=True, help="Output file")
args = parser.parse_args()

COLLAPSE_ID_FIELD = "CollapseId"
MATCH_ID_FIELD = "MatchId"
INVALID_POP_FREQ = -1

# map diagram
# merge SVID > CollapseId = MatchId
#      v         v            ^
# Population frequency    Collapse SVID

def main():
    print("# reading vcfs")
    vcf_input = VCF(args.input_vcf)
    vcf_collapsed = VCF(args.collapsed_vcf)
    vcf_merged = VCF(args.merged_vcf)

    print("# mapping merge svid to collapse id")
    merge_svid2collapseid = {}
    for v in vcf_merged:
        collapseid = str(v.INFO.get(COLLAPSE_ID_FIELD)) if v.INFO.get(COLLAPSE_ID_FIELD) is not None else str(-1)
        if ';' in v.ID:
            svids = v.ID.split(';')
            for svid in svids:
                merge_svid2collapseid[svid] = collapseid
        else:
            merge_svid2collapseid[v.ID] = collapseid

    print("# counting control samples with alt alleles")
    df_merged, num_control_samples = count_alt_samples(args.merged_vcf)

    print("# mapping merge svid to pop_freq")
    merge_svid2popfreq = {}
    for i, row in df_merged.iterrows():
        svid = str(row['SVID'])
        pop_freq = row['Alt_Sample_Count'] / num_control_samples
        if ';' in svid:
            svids = svid.split(';')
            for svid in svids:
                merge_svid2popfreq[svid] = pop_freq
        else:
            merge_svid2popfreq[svid] = pop_freq

    collapseid2popfreq = {}
    print("# mapping CollapseId to pop_freq")
    for svid, collapseid in merge_svid2collapseid.items():
        pop_freq = merge_svid2popfreq.get(svid, INVALID_POP_FREQ)
        if collapseid in collapseid2popfreq:
            pass
        else:
            collapseid2popfreq[collapseid] = pop_freq
    
    # matchid of collapsed vcf is a key equivalent to matchid of the merged vcf
    print("# mapping collapsed svids -> population frequency")
    collapsed_svids2popfreq = {}
    for v in vcf_collapsed:
        match_id = str(v.INFO.get(MATCH_ID_FIELD)) if v.INFO.get(MATCH_ID_FIELD) is not None else str(-1)
        pop_freq = collapseid2popfreq.get(match_id, INVALID_POP_FREQ)
        if ';' in v.ID:
            svids = v.ID.split(';')
            for svid in svids:
                collapsed_svids2popfreq[svid] = pop_freq 
        else:
            collapsed_svids2popfreq[v.ID] = pop_freq

    print("# getting population frequency for query SVIDs")
    outdata=[]
    svids_query = set(v.ID for v in vcf_input if v.ID is not None)
    for svid in svids_query:
        if svid in collapsed_svids2popfreq:
            pop_freq = collapsed_svids2popfreq[svid]
        elif svid in merge_svid2popfreq:
            pop_freq = merge_svid2popfreq[svid]
        else:
            pop_freq = INVALID_POP_FREQ
        outdata.append((svid, pop_freq))
    df = pd.DataFrame(outdata, columns=['svid', 'population_frequency'])
    df = df.sort_values(by='population_frequency', ascending=False)
    print("# writing output to {}".format(args.out))
    df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()