#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 17:07:54 2019

@author: sinnafoch
"""

from Bio.ExPASy import ScanProsite
from Bio import SeqIO

results=[]
print("\ttype\tdescription\tbegin\tend")
for seq in SeqIO.parse("/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/gh5_phylo/gh5_genes_prot_test.fa", "fasta"):
    handle = ScanProsite.scan(seq=seq.seq)#, mirror="https://prosite.expasy.org/")
    record = ScanProsite.read(handle)
    record["sequence_ac"] = seq.id
    record["total_length"] = len(seq.seq)
    results.append(record)


