#!/usr/bin/env python
# This script creates codon alignments

import argparse
#import glob
import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio import codonalign
from Bio.Alphabet import IUPAC
from Bio.codonalign import default_codon_table
import os
#from Bio.Seq import Seq

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="select_complete_buscos.py", description = """This script creates a codon alignment""", epilog = """written by Philipp Resl""")
pars.add_argument('-prot', dest="prot", required=True, help="Path to amino acid alignment")
pars.add_argument('-nucl', dest="nucl", required=True, help="Path to nucleotide sequence file")
pars.add_argument('-orth', dest="orth", required=True, help="Orthogroup mapping file from Orthofinder (Orthogroups.tsv)")
pars.add_argument('-id', dest="id", required=True, help="Name mapping ID file.")
args=pars.parse_args()


# this here is a quick hack for ambiguous codons when the sequences contains ambiguous nucleotides
# this here will make a decision on how these should be handled. 
# I have decided to treat them as ambiguous as well if there is doubt about the amino acid

#extracting orthogroup from filename
filename = os.path.basename(args.prot)
orthog = filename.split(".")[0]
#print(orthog)

#read id mappings:
ids_dict = {}
for line in open(args.id, "r"):
	abbr, full = line.strip().split("\t")
	ids_dict[abbr] = full
#print(ids_dict)


#read orthogroup mapping:
# parsing of the file expects that the file is correctly formatted as produced by orthofinder. so that tha \t delimits columns
# and the the provided OG ids refer to only single copy orthologs. This should be the case if the snakemake pipeline is followed.
sequences_names = []
for line in open(args.orth, "r"):
	if line.startswith(orthog):
		sequence_names = line.strip().split("\t")[1:]

#print(sequence_names)

#read nucl sequence file and subset to the needed sequences
nucleotide_seqs = SeqIO.to_dict(SeqIO.parse(args.nucl, "fasta", alphabet=IUPAC.ambiguous_dna))#IUPACUnambiguousDNA())
#print(len(nucleotide_seqs))
seqs_list = []
for rec in nucleotide_seqs.keys():
	#print(rec.id)
	for seq_name in sequence_names:
		if rec.startswith(seq_name):
			seqs_list.append(nucleotide_seqs[rec])
renamed_seqs = []
for seq in seqs_list:
	for name in ids_dict.keys():
		if seq.id.startswith(name):
			seq.id = ids_dict[name]
			renamed_seqs.append(seq)


#sys.exit(0)


"""
default_codon_table.forward_table['GKA'] = "X"
default_codon_table.forward_table['GGR'] = "G" # aa stays the same
default_codon_table.forward_table['CAR'] = "Q" # aa stays the same
default_codon_table.forward_table['GCM'] = "A" # aa stays the same
default_codon_table.forward_table['ATY'] = "I" # aa stays the same
default_codon_table.forward_table['ACM'] = "T" # aa stays the same
default_codon_table.forward_table['GAR'] = "E" # aa stays the same
default_codon_table.forward_table['GCY'] = "A" # aa stays the same
default_codon_table.forward_table['ACW'] = "T" # aa stays the same
default_codon_table.forward_table['YTG'] = "L" # aa stays the same
default_codon_table.forward_table['AGY'] = "S" # aa stays the same
default_codon_table.forward_table['TAY'] = "Y" # aa stays the same
default_codon_table.forward_table['GGK'] = "G" # aa stays the same
default_codon_table.forward_table['CCR'] = "P" # aa stays the same
default_codon_table.forward_table['TCR'] = "S" # aa stays the same
default_codon_table.forward_table['TCY'] = "S" # aa stays the same
default_codon_table.forward_table['AGR'] = "R" # aa stays the same
default_codon_table.forward_table['CCW'] = "P" # aa stays the same
default_codon_table.forward_table['GTY'] = "V" # aa stays the same
default_codon_table.forward_table['CCM'] = "P" # aa stays the same
default_codon_table.forward_table['CAS'] = "X"
default_codon_table.forward_table['GCR'] = "A" # aa stays the same
default_codon_table.forward_table['GGY'] = "G" # aa stays the same
default_codon_table.forward_table['CTR'] = "L" # aa stays the same
default_codon_table.forward_table['CTY'] = "L" # aa stays the same
default_codon_table.forward_table['CTM'] = "L" # aa stays the same
default_codon_table.forward_table['RTT'] = "X" 
default_codon_table.forward_table['AAR'] = "K" # aa stays the same
default_codon_table.forward_table['KTC'] = "X" 
default_codon_table.forward_table['ATM'] = "I" # aa stays the same
default_codon_table.forward_table['CTW'] = "L" # aa stays the same
default_codon_table.forward_table['CGY'] = "R" # aa stays the same
default_codon_table.forward_table['ACY'] = "T" # aa stays the same
default_codon_table.forward_table['MTG'] = "X"
default_codon_table.forward_table['CRT'] = "X"
default_codon_table.forward_table['TGY'] = "C" # aa stays the same
default_codon_table.forward_table['GTR'] = "V" # aa stays the same
default_codon_table.forward_table['AAY'] = "N" # aa stays the same
default_codon_table.forward_table['CAY'] = "H" # aa stays the same
default_codon_table.forward_table['CAM'] = "X"
default_codon_table.forward_table['WCT'] = "X" 
default_codon_table.forward_table['GCK'] = "A" # aa stays the same
default_codon_table.forward_table['SCG'] = "X"
default_codon_table.forward_table['AMT'] = "X"
"""

protein = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" }

ambig = {"M":["A","C"], "R":["A", "G"], "W": ["A", "T"], "S":["C", "G"], "Y": ["C", "T"], "K":["G", "T"], "V": ["A", "C", "G"], "H":["A", "C", "T"], "D":["A", "G", "T"], "B":["C","G", "T"], "N":["G", "A", "T", "C"]}

done =  False
while not done:
	try:
		nucl_list = []
		names_nucl = []
		nucl_unordered = renamed_seqs
		prot_alignment = AlignIO.read(args.prot, "fasta", alphabet=IUPAC.protein) 
		prot_ids = [pro.id for pro in prot_alignment]
		nucl = []
		# reorder ids in transcript file to match order in protein file
		for id in prot_ids:
			for nuc in nucl_unordered:
				if id == nuc.id:
					nucl.append(nuc)
					break
		#nucl_ids = [nuc.id for nuc in nucl]
		#print(nucl_ids, file=sys.stderr)
		#print(prot_ids, file=sys.stderr)
		try:
			codon = codonalign.build(prot_alignment, nucl, max_score=100, codon_table = default_codon_table)
		except IndexError as e:
			print("ERROR: Something is wrong with the file:", args.nucl, str(e), file=sys.stderr)
			sys.exit(0) #will use success to not mess up snakemake
		done = True
	except KeyError as e:
		print("ERROR WHEN WORKING ON FILE", args.nucl, str(e), file=sys.stderr)
		sys.exit(0)	#will use success to not mess up snakemake
	except RuntimeError as e:
		print("ERROR: Something is wrong with the file:", args.nucl, str(e), file=sys.stderr)
		sys.exit(0) #will use success to not mess up snakemake
	except ValueError as e:
		print("WARNING: ", e, "(will try to fix by resolving ambiguous codon)", file=sys.stderr)
		cod = str(e).split(")")[0].split("(")[-1]
		for key in ambig.keys():
			aa_list = []
			if key in cod:
				for amb in ambig[key]:
					modified = cod.replace(key, amb)
					if modified in protein.keys():
						aa_list.append(protein[modified])
					else:
						aa_list.append("X")
				if len(set(aa_list)) == 1:
					default_codon_table.forward_table[cod] = aa_list[0]
					print("\t", cod," codon resolves to ", aa_list[0], " from ", aa_list, file=sys.stderr)
				else:
					default_codon_table.forward_table[cod] = "X"
					print("\t", cod," codon resolves to X from ", aa_list, file=sys.stderr)
	except:
		print("ERROR: Something else is wrong with the file:", args.nucl, file=sys.stderr)
		sys.exit(0) #will use success to not mess up snakemake			
		
	
aln = codon.toMultipleSeqAlignment()

for rec in aln:
	print(">"+rec.id)
	print(rec.seq)
