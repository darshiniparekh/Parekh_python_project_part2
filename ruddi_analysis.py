#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import GC
import pandas as pd


def analyze_genome(sequence):
    length_of_genome = len(sequence)
    gc_content = GC(sequence)
    atg_forward = str(sequence).count("ATG")
    atg_reverse = str(sequence.reverse_complement()).count("ATG")
    return length_of_genome, gc_content, atg_forward, atg_reverse

def main():
    
    fasta_file = "GCF_000287275.1_ASM28727v1_genomic.fna"
    record = SeqIO.read(fasta_file, "fasta")

    
    length_of_genome, gc_content, atg_forward, atg_reverse = analyze_genome(record.seq)


    df = pd.DataFrame([[length_of_genome, gc_content, atg_forward, atg_reverse]],
                      columns=["Length_of_genome", "GC_content", "ATG_forward", "ATG_reverse"])

   
    df.to_csv("ruddi.csv", index=False)
    print("CSV file 'ruddi.csv' created successfully.")

if __name__ == "__main__":
    main()

