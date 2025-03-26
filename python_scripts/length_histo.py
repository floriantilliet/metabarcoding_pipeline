from Bio import SeqIO
import argparse
import os

import matplotlib.pyplot as plt

def plot_length_histogram(fasta_file, output_image):
    lengths = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.append(len(record.seq))
    
    plt.hist(lengths, bins=50, edgecolor='black')
    plt.title('Histogram of Sequence Lengths')
    plt.xlabel('Sequence Length')
    plt.ylabel('Frequency')
    plt.savefig(output_image)
    plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate a histogram of sequence lengths from a FASTA file.")
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("output_image", help="Output image file for the histogram")

    args = parser.parse_args()
    plot_length_histogram(args.fasta_file, args.output_image)