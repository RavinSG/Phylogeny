import json
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = "test@gmail.com"

proteins = ["single-stranded DNA-binding protein", "LysR family transcriptional regulator",
            "helix-turn-helix domain-containing protein", "efflux transporter outer membrane subunit"]

# proteins = ["FUSC family protein", "MFS transporter", "YraN family protein"]

file = pd.ExcelFile("Supporting_Documents/protein_tables.xlsx")
sheets = file.sheet_names
num_proteins = len(proteins)


def extract_protein_info(protein_sheets):
    protein_database = defaultdict(list)
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
                protein_data = [sheet, protein_row['Start'], protein_row['Stop'], protein_row['Strand']]
                protein_database[proteins[i]].append(protein_data)

        protein_profiles.append(one_hot)

    protein_profiles = np.array(protein_profiles)
    protein_count = np.sum(protein_profiles, axis=0)

    final_proteins = []
    for count, protein in zip(protein_count, proteins):
        if count != len(protein_sheets):
            continue
        final_proteins.append(protein)

    return protein_database, protein_profiles


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


# generate_gene_dict(sheets)
a,b = extract_protein_info(sheets)
