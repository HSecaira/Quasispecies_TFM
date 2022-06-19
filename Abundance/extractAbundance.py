"""
Script that get the abundances of all the nodes of the grandNetwork of a given region
"""

# Imports
import numpy as np;
import os, sys;

# Declare paths


seqsPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/NGS_analyses/out/Seqs_fasta_collapsed/";
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/GrandNetwork/"
dataPathOut = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/GrandNetwork/";
dicFolder = {}; 
dicFolder["r1"] = "Reg1/"; 
dicFolder["r2"] = "Reg2/"; 
dicFolder["r3"] = "Reg3/"; 

# Declare region

reg = "r2"
seqsPathIn_ = seqsPathIn + dicFolder[reg]
dataPathIn_ = dataPathIn + dicFolder[reg] + "MetaInfo/"
dataPathOut_ = dataPathOut + dicFolder[reg] + "MetaInfo/"


# Input file
fIn = os.path.join(dataPathIn_, "nodesID_" + reg + ".csv")
# Dictionary to store the abundances
abundances = {}


# Iterate over nodesID file. 
# Each nodeID corresponds to a genotype from the sequences-processed files (index 0-based)
# We want to extract the abundance from the sequenced-processed files
with open(fIn, "r") as fIn_:
	for ID, genotype in enumerate(fIn_):
		if len(genotype) != 0:
		#if ID == 0:
			print("\tSearching genotype for node", ID)
			genotype_ = genotype.strip()
			# Search the genotype in all the sequences-processed files using a command
			strOut = "grep -B 1 " + genotype_ + ' ' + seqsPathIn_ + "*.fasta" + " > " + dataPathOut_ + "temp.csv" 
			os.system(strOut)
			# Open temp.csv file and sum the abundance in case the genotypes exists in multiple files
			abundance = 0
			with open(dataPathOut_ + "/temp.csv", "r") as fInTemp:
				for line in fInTemp:
					if "->" in line:
						abundance += int(line.strip().split("-")[-1])

			# Append to dictionary
			abundances[ID] = abundance

			print("\tAbundance:", abundance)

			# Delete temp.csv
			strOut = "rm " + dataPathOut_ + "temp.csv"
			os.system(strOut)


# Save abundances
fOut = os.path.join(dataPathOut_, "abundances" + ".csv")
with open(fOut, "w") as fOut_:
	# Write into file
	fOut_.write("#NodeID,Abundance\n")
	for ID, abund in abundances.items():
		fOut_.write(str(ID) + "," + str(abund) + "\n")

















