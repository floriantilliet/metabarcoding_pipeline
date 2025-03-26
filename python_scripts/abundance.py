import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Calcul de l'abondance relative des phylums se terminant par 'mycota'.")
parser.add_argument("otu_table", type=str, help="Chemin du fichier OTU à analyser.")
args = parser.parse_args()

otu_table = pd.read_csv(args.otu_table, sep="\t")

otu_table["abundance"] = pd.to_numeric(otu_table["abundance"], errors="coerce")

def has_mycota(taxonomy):
    if isinstance(taxonomy, str):
        levels = taxonomy.split("|")
        for level in levels:
            if level.endswith("mycota"):
                return True
    return False

mycota_abundance = otu_table[otu_table["taxonomy"].apply(has_mycota)]["abundance"].sum()

total_abundance = otu_table["abundance"].sum()

relative_abundance = mycota_abundance / total_abundance if total_abundance > 0 else 0

print(f"Abondance totale des 'mycota' : {mycota_abundance}")
print(f"Abondance totale du jeu de données : {total_abundance}")
print(f"Abondance relative des 'mycota' : {relative_abundance:.4f} ({relative_abundance * 100:.2f}%)")
