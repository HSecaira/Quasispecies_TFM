""" 
Script to build the null model for networks
by selecting all the nucleotide sequences that translates to the same protein
"""


# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper as h
import networkx as nx
import timeit

# Variables and paths
print("\nLoading variables and paths\n")

# Dictionary containing the Reference sequences of each region (provided by Pilar)
refSeqs = {"R1": "CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCGTAAAGCGTTGAAACTTTG",
			"R2": "CGATTCATTTATCTTAAGTCGATAAATGCTTATTGCTCTCTTAGCGATATTGCGGCCTATCACGCCGATGGCGTGATAGTTGGCTTTTGGCGCGATCCATCCAGTGGTGGTGCCATACCGTTTGACTTCACTAAGTTTGATAAGACTAAATGTCCTATTCAAGCCGTGATAGTCGTTCCTCGTGCTTAGTAACTAAGGATGAAATGCATGTCTAAGACAGCATCTTCGCGTAACTCTCTCAGCGCACAATTGCGCCGAGCCGCGAACACAAGAATTGAGGTTGAAGGTAACCTCGCACTTTCCATTGCCAACGATTTACTGTTGGCCTA",
			"R3": "TTACACATTCGAGCTCGAGTCGCTTATTTTTGCTTCTCTCGCTCGTTCCGTTTGTGAGATACTGGACTTAGACTCGTCTGAGGTCACTGTTTACGGAGACGATATTATTTTACCGTCCTGTGCAGTCCCTGCCCTCCGGGAAGTTTTTAAGTATGTTGGTTTTACGACCAATACTAAAAAGACTTTTTCCGAGGGGCCGTTCAGAGAGTCGTGCGGCAAGCACTACTATTCTGGCGTAGATGTTACTCCCTTTTACATACGTCACCGTATAGTGA"}


# Translation table 11: The Bacterial, Archaeal and Plant Plastid Code from NCBI
transl_table = {"F": ["TTT", "TTC"], "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
				"I": ["ATT", "ATC", "ATA"], "M": ["ATG"], "V": ["GTT", "GTC", "GTA", "GTG"],
				"S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "P": ["CCT", "CCC", "CCA", "CCG"],
				"T": ["ACT", "ACC", "ACA", "ACG"], "A": ["GCT", "GCC", "GCA", "GCG"],
				"Y": ["TAT", "TAC"], "*": ["TAA", "TAG", "TGA"],
				"H": ["CAT", "CAC"], "Q": ["CAA", "CAG"],
				"N": ["AAT", "AAC"], "K": ["AAA", "AAG"],
				"D": ["GAT", "GAC"], "E": ["GAA", "GAG"],
				"C": ["TGT", "TGC"], "W": ["TGG"],
				"R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
				"G": ["GGT", "GGC", "GGA", "GGG"]}

# Dictionary containing the Reference protein sequences
refProts = {"A2_R1": "QQGQLYHNIDIVDGFDRRDIRLKSFTIKGERNGRPVNVSASLSAVDLFYSRLHTSNLPFATLDLDTTFSSFKHVLDSIFLLTQRVKR*",
			"A1_R2": "RFIYLKSINAYCSLSDIAAYHADGVIVGFWRDPSSGGAIPFDFTKFDKTKCPIQAVIVVPRA*",
			"Replicase_R2": "MSKTASSRNSLSAQLRRAANTRIEVEGNLALSIANDLLLA",
			"Replicase_R3": "YTFELESLIFASLARSVCEILDLDSSEVTVYGDDIILPSCAVPALREVFKYVGFTTNTKKTFSEGPFRESCGKHYYSGVDVTPFYIRHRIV"}

# Dictionary containing the positions for each region that codifies for aminoacids
refSeqsPos = {"A2_R1": [0, 263],
			"A1_R2": [0, 188],
			"Replicase_R2": [207, 326],
			"Replicase_R3": [1, 273]}

reg = "Reg3"
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/" + reg + "_wWM/"
dataPathOut = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg + "/results/"

print("\nVariables and paths loaded!\n")

##############################################################

# Star counting execution time
start = timeit.default_timer()

print("\nExtracting all the sequences from nodesID of the grandNetwork\n")
# Get all the sequence from the grandNetwork
fIn = os.path.join(dataPathIn + "nodesID_r3" + ".csv")
seqsGN = h.LoadFastaNodesGNID(fIn)
protName = "Replicase"

print("\nExtracting all the different proteins in the grandNetwork\n")
allProts = h.allDiffProts(seqsGN, refSeqsPos["Replicase_R3"], transl_table, refProts["Replicase_R3"])
# Saving into a file
fOut = os.path.join(dataPathOut + "AllProts_" + protName + "_" + reg + ".csv")
h.saveallDiffProts(fOut, allProts)
print("Number of different prots", len(allProts))

##############################################################3

print("\nFind all sequences for  neutral networks\n")
seqsNN = h.allSeqsNN(allProts, seqsGN, refSeqsPos["Replicase_R3"], transl_table)
print("\nAll sequences for neutral networks found!\n")

print("\nAll sequences checked!\n")
dataPathOut_ = dataPathOut + "NeutralNetworks" + protName + "/"
h.saveSeqsNN(seqsNN, dataPathOut_, protName, reg)
print("Number of different Prots", len(seqsNN.keys()))

#################################################################

print("\nSaving MetaInfo of Neutral Netowrks and doing a nice plot\n")
# Get the number of sequences for each file
dTimesNN = h.numSeqsPerFile(seqsNN)
# Sort the dictionary and save the Info into a file
fOut1 = os.path.join(dataPathOut + "NeutralNetworksMetaInfo_" + protName + "_" + reg + ".csv")
h.sortDictSaveInfo(fOut1, dTimesNN)

# Get the counts of sequencies across files
unique, counts = np.unique(np.array(list(dTimesNN.values())), return_counts= True)

# Plot 
fig = plt.figure()
plt.loglog(unique, counts/np.sum(counts), "bo-")
plt.xlabel("Number of Sequences")
plt.ylabel("Frequency")
plt.title("Distribution of the Number of Sequences per Neutral Network")
fig.savefig(dataPathOut + "NumberSeqsNN_distribution_" + protName + "_" + reg + ".pdf", bbox_inches = "tight")
fig.savefig(dataPathOut + "NumberSeqsNN_distribution_" + protName + "_" + reg + ".eps", bbox_inches = "tight")
#plt.show()

print("\nAll done!")

# Get execution time
stop = timeit.default_timer()
total_time = stop - start
# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)
sys.stdout.write("Total running time: %d:%d:%d. s\n" % (hours, mins, secs))

# Command to see how many sequences per neutral network file
#for f in $(find . -name "*.csv"); do cat $f | grep ">" | wc -l; done |sort -nr | uniq -c







