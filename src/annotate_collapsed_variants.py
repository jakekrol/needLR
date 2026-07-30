#!/usr/bin/env python3
import argparse
from cyvcf2 import VCF
import os
import pandas as pd

parser = argparse.ArgumentParser(description="Annotate collapsed needLR variants")
parser.add_argument("--input_vcf", required=True, help="Path to the input VCF file")
parser.add_argument("--collapsed_vcf", required=True, help="Path to the collapsed VCF file")
parser.add_argument("--table", required=True, help="Path to the needLR output table")
parser.add_argument("--out", required=True, help="Output file")
args = parser.parse_args()

COL_SVID='SVID'
COL_POP_FREQ='Pop_Freq_ALL'
COL_COLLAPSE_ID='CollapseId'
INVALID_POP_FREQ=-1


def main():
    # inputs
    df_needlr = pd.read_csv(args.table, sep="\t")
    vcf_input = VCF(args.input_vcf)
    vcf_collapsed = VCF(args.collapsed_vcf)

    # gather input SVIDs
    svids_query = set(v.ID for v in vcf_input if v.ID is not None)
    df_query = pd.DataFrame([(svid) for svid in svids_query], columns=[COL_SVID])

    # get SVID, Pop_Freq_ALL, and Collapse_Id in final output
    df_needlr=df_needlr[[COL_SVID, COL_POP_FREQ, COL_COLLAPSE_ID]]
    collapse_ids = set(df_needlr[COL_COLLAPSE_ID].unique())
    # rm '.' which is the case where the variant was not collapsed
    collapse_ids.discard('.')

    # get SVID and MatchId from collapsed VCF
    # MatchId is a key shared with Collapse_Id used to resolve variant merging
    collapsed_data = []
    for v in vcf_collapsed:
        svid = v.ID if v.ID is not None else '.'
        # we only care about query SVs, not those from their reference set which were collapsed
        if svid in svids_query:
            match_id = v.INFO.get('MatchId') if v.INFO.get('MatchId') is not None else '.'
            if match_id in collapse_ids:
                mask = df_needlr[COL_COLLAPSE_ID] == match_id
                pop_freq = df_needlr.loc[mask, COL_POP_FREQ].values[0] if mask.any() else INVALID_POP_FREQ
                collapsed_data.append((svid, pop_freq))
            else:
                collapsed_data.append((svid, match_id, INVALID_POP_FREQ))
    df_collapsed = pd.DataFrame(collapsed_data, columns=[COL_SVID, COL_POP_FREQ])

    # get SVID and Pop_Freq_ALL from needLR output table
    needlr_data=[]
    for i, row in df_needlr.iterrows():
        svid = row[COL_SVID]
        # needlr stores the SVIDs as semicolon delimited strings, so we need to split them and check each one
        if ';' in svid:
            svids = svid.split(';')
            for svid in svids:
                if svid in svids_query:
                    pop_freq = row[COL_POP_FREQ] if row[COL_POP_FREQ] is not None else INVALID_POP_FREQ
                    needlr_data.append((svid, pop_freq))
        else:
            if svid in svids_query:
                pop_freq = row[COL_POP_FREQ] if row[COL_POP_FREQ] is not None else INVALID_POP_FREQ
                needlr_data.append((svid, pop_freq))
    df_needlr_filtered = pd.DataFrame(needlr_data, columns=[COL_SVID, COL_POP_FREQ])
    # combine the two dataframes and remove duplicates
    df_concat = pd.concat([df_collapsed, df_needlr_filtered], axis=0, ignore_index=True)
    df_concat.drop_duplicates(subset=[COL_SVID], inplace=True)
    # finally force any missing query SVIDs to have an invalid Pop_Freq_ALL value
    final_svids = set(df_concat[COL_SVID].unique())
    for svid in svids_query:
        if svid not in final_svids:
            df_concat = pd.concat([df_concat, pd.DataFrame([(svid, INVALID_POP_FREQ)], columns=[COL_SVID, COL_POP_FREQ])], axis=0, ignore_index=True)
    df_concat.sort_values(by=[COL_POP_FREQ], inplace=True, ascending=False)
    df_concat.to_csv(args.out, sep="\t", index=False)
    

if __name__ == "__main__":
    main()