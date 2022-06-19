"""

	script_visualizeGrandNetworkFull.py: 

		This script loads the edges generated for the grand network and attempts a rough visualization. 
		AUthor Luiño, modified by Henry
"""

# Imports: 
import numpy as np; 
import matplotlib.pyplot as plt; 
from mpl_toolkits import mplot3d; 
import networkx as nx; 
import os, sys; 
import helper as h;


# Defining paths: 
print("\nLoading variables and paths\n")
dataPathInR1 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Reg1_wWM/GN/whereMutations.csv"
dataPathInR2 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Reg2_wWM/GN/whereMutations.csv"
dataPathInR3 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Reg3_wWM/GN/whereMutations.csv" 
picsPath = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/plotsMutations/" 

allDataR1 = np.loadtxt(dataPathInR1, skiprows=0); 
allDataR2 = np.loadtxt(dataPathInR2, skiprows=0); 
allDataR3 = np.loadtxt(dataPathInR3, skiprows=0); 


# Variables and paths
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

# Dictionary containing the neutral networks for silent mutations
network_mutations = {"F": np.array([[0,1],[1,0]], dtype=float), 
					"L": np.array([[0, 1, 1, 1, 0, 0],
									[1, 0, 1, 1, 1, 0], 
									[1, 1, 0, 1, 0, 0],
									[1, 1, 1, 0, 0, 1],
									[0, 1, 0, 0, 0, 1],
									[0, 0, 0, 1, 1, 0]], dtype=float),
					"I": np.array([[0, 1, 1],
									[1, 0, 1],
									[1, 1, 0]], dtype=float),
					"M": np.array([0]),
					"V": np.array([[0, 1, 1, 1],
									[1, 0, 1, 1],
									[1, 1, 0, 1],
									[1, 1, 1, 0]], dtype=float),
					"S": np.array([[0, 1, 1, 1, 0, 0],
									[1, 0, 1, 1, 0, 0],
									[1, 1, 0, 1, 0, 0],
									[1, 1, 1, 0, 0, 0],
									[0, 0, 0, 0, 0, 1],
									[0, 0, 0, 0, 1, 0]], dtype=float),
					"P": np.array([[0, 1, 1, 1],
									[1, 0, 1, 1],
									[1, 1, 0, 1],
									[1, 1, 1, 0]], dtype=float),
					"T": np.array([[0, 1, 1, 1],
									[1, 0, 1, 1],
									[1, 1, 0, 1],
									[1, 1, 1, 0]], dtype=float),
					"A": np.array([[0, 1, 1, 1],
									[1, 0, 1, 1],
									[1, 1, 0, 1],
									[1, 1, 1, 0]], dtype=float),
					"Y": np.array([[0, 1], [1, 0]], dtype=float),
					"*": np.array([[0, 1, 1],
									[1, 0, 0],
									[1, 0, 0]], dtype=float),
					"H": np.array([[0, 1], [1, 0]], dtype=float),
					"Q": np.array([[0, 1], [1, 0]], dtype=float),
					"N": np.array([[0, 1], [1, 0]], dtype=float),
					"K": np.array([[0, 1], [1, 0]], dtype=float),
					"D": np.array([[0, 1], [1, 0]], dtype=float),
					"E": np.array([[0, 1], [1, 0]], dtype=float),
					"C": np.array([[0, 1], [1, 0]], dtype=float), 
					"W": np.array([0]),
					"R": np.array([[0, 1, 1, 1, 0, 0],
									[1, 0, 1, 1, 1, 0], 
									[1, 1, 0, 1, 0, 0],
									[1, 1, 1, 0, 0, 1],
									[0, 1, 0, 0, 0, 1],
									[0, 0, 0, 1, 1, 0]], dtype=float),
					"G": np.array([[0, 1, 1, 1],
									[1, 0, 1, 1],
									[1, 1, 0, 1],
									[1, 1, 1, 0]], dtype=float)}

##############################################################################


reg = "Reg1"
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg + "/"
dataPathOut = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg + "/results/"

print("\nVariables and paths loaded!\n")


print("\nLoading codons info\n")

fIn = os.path.join(dataPathIn + "CodonInfo_A2_R1" + ".csv")
codonsInfo = h.loadCodonsInfo(fIn)

print("\nCodons info loaded!\n")


## Random walk
print("\nCalculating the probabilities of silent mutations\n")

# Convert the sequence into array of positions
refSeqsPos = h.convertRefSeqs(refSeqs)
R1_refSeqsPos = refSeqsPos["R1"]
lenA2 = 264
protName = "A2"

# Get the frequency of silent mutations for the protein
freqMutations = h.probSilentMutationsRegion(271, transl_table, R1_refSeqsPos, codonsInfo, network_mutations, centrality=True)

# Save mutations into a file
fOut = os.path.join(dataPathOut + "freqSMutations_" + protName + "_" + reg + ".csv")
h.saveFreqMutations(fOut, freqMutations)

# Load mutations from a file
#fIn1 = os.path.join(dataPathOut + "MutationsRandomWalk_" + reg + ".csv")
#mutations = h.loadMutations(fIn1)
###################################################



fig = plt.figure(figsize=(10,6));
width = 0.97 
h1 = plt.hist(allDataR1, 271, density=True, color="#4cc9f0", edgecolor="#adb5bd",
				linewidth=0.01); 
# Save empirical frequencies of mutation
fOutEmpFreqs = os.path.join(dataPathOut, "freqEmpMutations_" + reg + ".csv")
h.saveEmpiricalFreqMutations(fOutEmpFreqs, np.rint(h1[1]), h1[0])

plt.bar(R1_refSeqsPos + width/2, freqMutations, color = "#f72585", width = width,
		edgecolor="#adb5bd", linewidth=0.01)
plt.axvline(x = 263, color = "red", linestyle="dashed")
plt.xlabel("Position"); 
plt.ylabel("Frequency of mutations"); 
#plt.title("Region 1"); 
plt.yscale("log")
plt.ylim(1e-4, 1e-1)
#fig.savefig(picsPath + "whereMutationsNMutReg1" + ".pdf", bbox_inches = "tight"); 
#fig.savefig(picsPath + "whereMutationsNMutReg1" + ".svg", bbox_inches = "tight"); 


# Plot of empirical vs null model frequency of mutation
fig = plt.figure()
plt.plot(freqMutations[0:270], freqMutations[0:270], "ko-")
plt.scatter(freqMutations[0:270], h1[0][0:270])
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Frequency of mutations (Null model)")
plt.ylabel("Frequency of  mutations (Empirical)")


##################################################################3
##reg = "Reg2"
##dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg + "/"
##dataPathOut = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg + "/results/"
##
##print("\nVariables and paths loaded!\n")
##
##
##print("\nLoading codons info\n")
##
##fIn = os.path.join(dataPathIn + "CodonInfo_A1_R2" + ".csv")
##codonsInfo = h.loadCodonsInfo(fIn)
##fIn1 = os.path.join(dataPathIn + "CodonInfo_Replicase_R2" + ".csv")
##codonsInfo1 = h.loadCodonsInfo(fIn1)
##
##print("\nCodons info loaded!\n")
##
##
#### Random walk
##print("\nCalculating the probabilities of silent mutations\n")
##
### Convert the sequence into array of positions
##refSeqsPos = h.convertRefSeqs(refSeqs)
##R1_refSeqsPos = refSeqsPos["R2"]
##lenA2 = 189
##protName = "A1"
##
##refSeqsPos1 = h.convertRefSeqs(refSeqs)
##R1_refSeqsPos1 = refSeqsPos1["R2"]
##lenA21 = 122
##protName1 = "Replicase"
##
### Get the frequency of silent mutations for the protein
##freqMutations = h.probSilentMutationsRegion(329, transl_table, R1_refSeqsPos, codonsInfo, network_mutations)
##freqMutations1 = h.probSilentMutationsRegion(329, transl_table, R1_refSeqsPos1, codonsInfo1, network_mutations)
##
##
### Save mutations into a file
##fOut = os.path.join(dataPathOut + "freqSMutations_" + protName + "_" + reg + ".csv")
##h.saveFreqMutations(fOut, freqMutations)
##fOut1 = os.path.join(dataPathOut + "freqSMutations_" + protName1 + "_" + reg + ".csv")
##h.saveFreqMutations(fOut1, freqMutations1)
##
### Load mutations from a file
###fIn1 = os.path.join(dataPathOut + "MutationsRandomWalk_" + reg + ".csv")
###mutations = h.loadMutations(fIn1)
##
##
### Plot mutations distribution
##fig = plt.figure(figsize=(10,6))
##width = 0.97
##h2 = plt.hist(allDataR2, 329, density=True, color="#4cc9f0", edgecolor="#adb5bd",
				##linewidth=0.08); 
##
### Save empirical frequencies of mutation
##fOutEmpFreqs = os.path.join(dataPathOut, "freqEmpMutations_" + reg + ".csv")
##h.saveEmpiricalFreqMutations(fOutEmpFreqs, np.rint(h2[1]), h2[0])
##
##plt.bar(R1_refSeqsPos + width/2, freqMutations, color = "#f72585", width = width, 
		##edgecolor="#adb5bd", linewidth=0.08)
##plt.bar(R1_refSeqsPos1 + width/2, freqMutations1, color = "#f72585", width = width,
		##edgecolor="#adb5bd", linewidth=0.08)
##plt.axvline(x=188, color="red", linestyle="dashed")
##plt.axvline(x=207, color="red", linestyle="dashed")
##plt.yscale("log")
##plt.ylim(1e-4, 1e-1)
###plt.title("Region2")
##plt.xlabel("Position")
##plt.ylabel("Frequency of mutations")
##plt.savefig(picsPath + "whereMutationsNMut" + reg + ".pdf", bbox_inches = "tight")
##plt.savefig(picsPath + "whereMutationsNMut" + reg + ".svg", bbox_inches = "tight")
##
###################################################################################
##reg = "Reg3"
##dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg + "/"
##dataPathOut = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg + "/results/"
##
##print("\nVariables and paths loaded!\n")
##
##
##print("\nLoading codons info\n")
##
##fIn = os.path.join(dataPathIn + "CodonInfo_Replicase_R3" + ".csv")
##codonsInfo = h.loadCodonsInfo(fIn)
##
##print("\nCodons info loaded!\n")
##
##
#### Random walk
##print("\nCalculating the probabilities of silent mutations\n")
##
### Convert the sequence into array of positions
##refSeqsPos = h.convertRefSeqs(refSeqs)
##R1_refSeqsPos = refSeqsPos["R3"]
##lenA2 = 275
##protName = "Replicase"
##
### Get the frequency of silent mutations for the protein
##freqMutations = h.probSilentMutationsRegion(275, transl_table, R1_refSeqsPos, codonsInfo, network_mutations)
##
### Save mutations into a file
##fOut = os.path.join(dataPathOut + "freqSMutations_" + protName + "_" + reg + ".csv")
##h.saveFreqMutations(fOut, freqMutations)
##
### Load mutations from a file
###fIn1 = os.path.join(dataPathOut + "MutationsRandomWalk_" + reg + ".csv")
###mutations = h.loadMutations(fIn1)
##
##
### Plot mutations distribution
##fig = plt.figure(figsize=(10,6))
##width = 0.97
##h3 = plt.hist(allDataR3, 275, density=True, color="#4cc9f0", edgecolor="#adb5bd",
				##linewidth=0.08); 
##
### Save empirical frequencies of mutation
##fOutEmpFreqs = os.path.join(dataPathOut, "freqEmpMutations_" + reg + ".csv")
##h.saveEmpiricalFreqMutations(fOutEmpFreqs, np.rint(h3[1]), h3[0])
##
##plt.bar(R1_refSeqsPos + width/2, freqMutations, color = "#f72585", width = width,
		##edgecolor="#adb5bd", linewidth=0.08)
##plt.yscale("log")
##plt.ylim(1e-4, 1e-1)
###plt.title("Region3")
##plt.xlabel("Position")
##plt.ylabel("Frequency of Mutations")
##plt.savefig(picsPath + "whereMutationsNMut" + reg + ".pdf", bbox_inches = "tight")
##plt.savefig(picsPath + "whereMutationsNMut" + reg + ".svg", bbox_inches = "tight")


plt.show()





