import pandas as pd
import sys

def main(argv):
    input = sys.argv[1]
    df = pd.read_csv(input,sep="\t")
    df_species = df.groupby(["species","genus","family","order","class","phylum","superkingdom"]).sum()
    df_genus = df.groupby(["genus","family","order","class","phylum","superkingdom"]).sum()
    samples = list(df_genus.columns.values)
    df_species = df_species.sort_values(by=samples, ascending=[0]*len(samples))
    df_genus = df_genus.sort_values(by=samples, ascending=[0]*len(samples))
    df = df.sort_values(by=samples, ascending=[0]*len(samples))
    for sample in samples:
        df_genus.insert(df_genus.columns.get_loc(sample)+1,sample+"_rel",df_genus[[sample]]*100/df_genus[[sample]].sum()) 
        df_species.insert(df_species.columns.get_loc(sample)+1,sample+"_rel",df_species[[sample]]*100/df_species[[sample]].sum()) 
    df = df.round(decimals=3)
    df.to_csv(input[:-4]+"_taxid.tsv",sep="\t")
    df_genus = df_genus.round(decimals=3)
    df_genus.to_csv(input[:-4]+"_genus.tsv",sep="\t")
    df_species = df_species.round(decimals=3)
    df_species.to_csv(input[:-4]+"_species.tsv",sep="\t")

if __name__ == "__main__":
   main(sys.argv)