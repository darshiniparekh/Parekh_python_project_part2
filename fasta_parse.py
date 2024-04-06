#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
from Bio import Entrez, SeqIO
import pandas as pd

def download_fasta_files(accessions, output_dir):
    Entrez.email = "your@email.com"  # Enter your email address here

    for accession in accessions:
        handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
        filename = f"{output_dir}/protein_{accession}.fasta"

        with open(filename, "w") as out_handle:
            out_handle.write(handle.read())
        print(f"Downloaded: {filename}")

def parse_fasta(filepath):
    record = SeqIO.read(filepath, "fasta")
    accession = record.id
    first_10_aa = str(record.seq[:10])
    length = len(record.seq)
    num_cysteines = str(record.seq).count("C")

    return accession, first_10_aa, length, num_cysteines

def process_fasta_files(file_paths):
    results = []

    for filepath in file_paths:
        try:
            result = parse_fasta(filepath)
            results.append(result)
        except FileNotFoundError as e:
            print(f"File not found: {e.filename}. Skipping...")

    return results

def main():
    accessions = ["AGI40145.1", "AGJ87295.1", "WVV45440.1", "WVS05366.1"]
    output_dir = "C:/Users/darsh/Parekh_python_project_part2"

    
    download_fasta_files(accessions, output_dir)

    file_paths = [f"{output_dir}/protein_{accession}.fasta" for accession in accessions]
    data = process_fasta_files(file_paths)

    
    df = pd.DataFrame(data, columns=["ID", "First_10_AA", "Length", "Number_Cs"])
    df.set_index("ID", inplace=True)

    
    df.to_csv("protein_info.csv")
    print("CSV file 'protein_info.csv' created successfully.")

if __name__ == "__main__":
    main()

