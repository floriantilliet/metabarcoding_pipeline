import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser(description="Analyse et affichage de la distribution de l'identité des OTUs.")
parser.add_argument("otu_table", type=str, help="Chemin du fichier OTU à analyser.")
parser.add_argument("--save", type=str, default="identity_distribution.png", help="Nom du fichier de sortie (format PNG).")
args = parser.parse_args()

df = pd.read_csv(args.otu_table, sep="\t")

identity_stats = df["identity"].describe()
identity_mean = identity_stats["mean"]
identity_median = identity_stats["50%"]
identity_std = identity_stats["std"]
identity_min = identity_stats["min"]
identity_max = identity_stats["max"]
identity_q1 = identity_stats["25%"]
identity_q3 = identity_stats["75%"]

stats_text = (
    f"std : {identity_std:.2f}%\n"
    f"variance : {identity_stats['std']**2:.2f}%\n"
    f"min : {identity_min:.2f}% | max : {identity_max:.2f}%\n"
    f"Q1 : {identity_q1:.2f}% | Q3 : {identity_q3:.2f}%"
)

plt.figure(figsize=(10, 6))
plt.hist(df["identity"], bins=30, edgecolor='black', alpha=0.7)
plt.axvline(identity_mean, color='red', linestyle='dashed', linewidth=2, label=f'Moyenne: {identity_mean:.2f}')
plt.axvline(identity_median, color='blue', linestyle='dashed', linewidth=2, label=f'Médiane: {identity_median:.2f}')
plt.xlabel("Identity")
plt.ylabel("Number of OTUs")

plt.legend(loc="upper left")
plt.text(0.02, 0.85, stats_text, fontsize=12, bbox=dict(facecolor='white', alpha=0.8), transform=plt.gca().transAxes, verticalalignment='top')

plt.show()

plt.savefig(args.save, dpi=300, bbox_inches="tight")
