import json
import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO

Entrez.email = "test@gmail.com"

proteins = ["single-stranded DNA-binding protein", "LysR family transcriptional regulator",
            "helix-turn-helix domain-containing protein", "efflux transporter outer membrane subunit"]

file = pd.ExcelFile("Supporting_Documents/protein_tables.xlsx")
sheets = file.sheet_names
protein_profiles = []

for sheet in sheets:
    df = file.parse(sheet)
    num_proteins = len(proteins)
    one_hot = [0] * num_proteins

    sheet_proteins = df['Protein name'].values

    for i in range(num_proteins):
        if proteins[i] in sheet_proteins:
            one_hot[i] = 1

    protein_profiles.append(one_hot)

protein_profiles = np.array(protein_profiles)
protein_count = np.sum(protein_profiles, axis=0)

final_proteins = []
for count, protein in zip(protein_count, proteins):
    if count != len(sheets):
        continue
    final_proteins.append(protein)


def download_genome(species_id):
    print(f"Downloading {species_id}...")
    handle = Entrez.efetch(db="nucleotide", id=species_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    return record.seq


def generate_gene_dict(species):
    gene_dict = {}
    for gene in species:
        gene = gene.strip()
        gene_dict[gene] = download_genome(gene)._data

    with open("genomes.json", 'w') as gf:
        json.dump(gene_dict, gf)


generate_gene_dict(sheets)
