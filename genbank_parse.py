#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python

from Bio import Entrez, SeqIO
import pandas as pd

def download_genbank_files(accessions, output_dir):
    Entrez.email = "your@email.com"  

    for accession in accessions:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        filename = f"{output_dir}/sequence_{accession}.gb"

        with open(filename, "w") as out_handle:
            out_handle.write(handle.read())
        print(f"Downloaded: {filename}")

def parse_genbank(filepath):
    record = SeqIO.read(filepath, "genbank")
    accession = record.annotations["accessions"][0]
    family = record.annotations["taxonomy"][2] if len(record.annotations["taxonomy"]) >= 3 else None
    genus = record.annotations["taxonomy"][1] if len(record.annotations["taxonomy"]) >= 2 else None
    species = record.annotations["taxonomy"][0] if len(record.annotations["taxonomy"]) >= 1 else None
    num_features = len(record.features)
    source = record.annotations["source"]

    return accession, family, genus, species, num_features, source

def process_genbank_files(file_paths):
    results = []

    for filepath in file_paths:
        try:
            result = parse_genbank(filepath)
            results.append(result)
        except FileNotFoundError as e:
            print(f"File not found: {e.filename}. Skipping...")

    return results

def main():
    accessions = [
        "NZ_CALPCP010000001.1",
        "NZ_CALPCY010000130.1",
        "NZ_BHVZ01000001.1",
        "NZ_SRYA01000017.1",
        "NZ_CAJTFZ010000019.1"
    ]
    output_dir = "C:/Users/darsh/Parekh_python_project_part2"

   
    download_genbank_files(accessions, output_dir)

    file_paths = [f"{output_dir}/sequence_{accession}.gb" for accession in accessions]
    data = process_genbank_files(file_paths)

    
    df = pd.DataFrame(data, columns=["Accession", "Family", "Genus", "Species", "Num_Features", "Source"])

    
    df.to_csv("genbank_parse.csv", index=False)
    print("CSV file 'genbank_parse.csv' created successfully.")

if __name__ == "__main__":
    main()

