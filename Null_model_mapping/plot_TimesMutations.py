""" 
Script to generate a plot of frequency of mutations 
per codon position of the three sequenced regions
"""

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper as h
import networkx as nx



# Variables and paths
print("\nLoading variables and paths\n")

picsPath = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/plotsMutations/" 

reg1 = "Reg1/";
dataPathInR1 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg1 + "/results/"

reg2 = "Reg2/";
dataPathInR2 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg2 + "/results/"

reg3 = "Reg3/";
dataPathInR3 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg3 + "/results/"

print("\nAll paths and variables loaded!")

# Load data
fInR1 = os.path.join(dataPathInR1, "freqSMutations_A2_Reg1" + ".csv")
freqMutationsR1 = h.loadFreqMutations(fInR1)

fInR2_1 = os.path.join(dataPathInR2, "freqSMutations_A1_Reg2" + ".csv")
freqMutationsR2_1 = h.loadFreqMutations(fInR2_1)
fInR2_2 = os.path.join(dataPathInR2, "freqSMutations_Replicase_Reg2" + ".csv")
freqMutationsR2_2 = h.loadFreqMutations(fInR2_2)


fInR3 = os.path.join(dataPathInR3, "freqSMutations_Replicase_Reg3" + ".csv")
freqMutationsR3 = h.loadFreqMutations(fInR3)


# Iterate over the data
mutationsInThisCodonPostition_r1 = np.zeros(3); 
for position, times in freqMutationsR1.items():
	if position <= 263: 
		iInCodon = position%3; 
		mutationsInThisCodonPostition_r1[iInCodon] += times

mutationsInThisCodonPostition_r2 = np.zeros(3); 
for position, times in freqMutationsR2_1.items():
	if position <= 188: 
		iInCodon = position%3; 
		mutationsInThisCodonPostition_r2[iInCodon] += times

for position, times in freqMutationsR2_2.items():
	if position >= 207: 
		iInCodon = position%3; 
		mutationsInThisCodonPostition_r2[iInCodon] += times

mutationsInThisCodonPostition_r3 = np.zeros(3); 
pos = 0
for position, times in freqMutationsR3.items():
	if position >= 1:
		iInCodon = pos%3; 
		mutationsInThisCodonPostition_r3[iInCodon] += times
		pos += 1

print(mutationsInThisCodonPostition_r1) 
print(mutationsInThisCodonPostition_r2)
print(mutationsInThisCodonPostition_r3)

fig = plt.figure(); 
plt.subplot(311); 
plt.plot(mutationsInThisCodonPostition_r1, 'o-k');
plt.yscale("log") 
plt.subplot(312); 
plt.plot(mutationsInThisCodonPostition_r2, 'o-k');
plt.yscale("log") 
plt.subplot(313); 
plt.plot(mutationsInThisCodonPostition_r3, 'o-k');
plt.yscale("log") 
fig.savefig(picsPath + "codonsMutModel" + ".pdf", bbox_inches = "tight"); 
fig.savefig(picsPath + "codonsMutModel" + ".svg", bbox_inches = "tight");  

plt.show()













