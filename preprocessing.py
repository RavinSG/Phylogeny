import os
import json
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Entrez.email = "test@gmail.com"

# List of proteins to analyze
proteins = ["single-stranded DNA-binding protein", "LysR family transcriptional regulator",
            "helix-turn-helix domain-containing protein", "efflux transporter outer membrane subunit"]

# proteins = ["FUSC family protein", "MFS transporter", "YraN family protein"]

file = pd.ExcelFile("Supporting_Documents/protein_tables.xlsx")
sheets = file.sheet_names
num_proteins = len(proteins)


# Extract the starting and ending location of each protein for the genome given
def extract_protein_info(protein_sheets):
    protein_database = defaultdict(dict)
    protein_profiles = []
    for sheet in protein_sheets:
        df = file.parse(sheet)
        one_hot = [0] * num_proteins
        sheet_proteins = df['Protein name'].apply(lambda x: x.split(":")[-1].strip()).values

        for i in range(num_proteins):
            if proteins[i] in sheet_proteins:
                one_hot[i] = 1
                sheet_index = np.where(sheet_proteins == proteins[i])[0][0]
                protein_row = df.iloc[sheet_index]
                protein_database[proteins[i]][sheet] = \
                    [protein_row['Start'], protein_row['Stop'], protein_row['Strand']]

        protein_profiles.append(one_hot)

    protein_profiles = np.array(protein_profiles)
    protein_count = np.sum(protein_profiles, axis=0)

    final_proteins = []
    for count, f_protein in zip(protein_count, proteins):
        if count != len(protein_sheets):
            continue
        final_proteins.append(f_protein)

    return protein_database, protein_profiles


# If not available, download the genome sequences from www.ncbi.nlm.nih.gov
def download_genome(species_id):
    print(f"Downloading {species_id}")
    handle = Entrez.efetch(db="nucleotide", id=species_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    return record.seq


# Save the downloaded genomes to a json file for future uses
def generate_gene_dict(species):
    gene_dict = {}
    for gene in species:
        gene = gene.strip()
        gene_dict[gene] = download_genome(gene)._data

    with open("genomes.json", 'w') as g_file:
        json.dump(gene_dict, g_file)

    return gene_dict


def write_to_fasta(protein_name, sequences, out_file):
    SeqIO.write(sequences[protein_name].values(), out_file, "fasta")


# If the genome dictionary is not available, download them from ncbi
if os.path.exists("genomes.json"):
    print("Genome dictionary found, loading data...")
    with open("genomes.json") as gf:
        genome_dict = json.load(gf)
else:
    print("Genome dictionary not found, downloading data...")
    genome_dict = generate_gene_dict(sheets)

p_database, _ = extract_protein_info(sheets)

# Add all the sequences relating to a protein to a seperate dictionary
extracted_sequences = defaultdict(dict)
for genome in sheets:
    genome_seq = Seq(genome_dict[genome.strip()])
    for protein in proteins:
        protein_data = p_database[protein][genome]
        if protein_data[-1] == '+':
            extracted_sequences[protein][genome] = SeqRecord(
                genome_seq[protein_data[0]:protein_data[1]], id=genome, description=protein)
        else:
            extracted_sequences[protein][genome] = SeqRecord(
                genome_seq[protein_data[0]:protein_data[1]].reverse_complement(), id=genome, description=protein)

# Save each protein in a separate file
for protein in proteins:
    write_to_fasta(protein, extracted_sequences, f"{protein}_sequences.fasta")
