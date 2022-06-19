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
dataPathInR1 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Reg1_wWM/GN/whereMutations.csv"
dataPathInR2 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Reg2_wWM/GN/whereMutations.csv"
dataPathInR3 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Reg3_wWM/GN/whereMutations.csv" 
picsPath = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/plotsMutations/" 

allDataR1 = np.loadtxt(dataPathInR1, skiprows=0); 
allDataR2 = np.loadtxt(dataPathInR2, skiprows=0); 
allDataR3 = np.loadtxt(dataPathInR3, skiprows=0); 

reg1 = "Reg1/";
dataPathInR1Model = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg1 + "/results/"
reg2 = "Reg2/";
dataPathInR2Model = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg2 + "/results/"
reg3 = "Reg3/";
dataPathInR3Model = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg3 + "/results/"

##fig = plt.figure(); 
##h1 = plt.hist(allDataR1, 271, density=True, color="#212529"); 
##plt.xlabel("Position"); 
##plt.ylabel("Frequency of mutations"); 
##plt.title("Region 1"); 
###fig.savefig(picsPath + "whereMutationsReg1" + ".pdf", bbox_inches = "tight"); 
###fig.savefig(picsPath + "whereMutationsReg1" + ".svg", bbox_inches = "tight"); 
##
##
##fig = plt.figure(); 
##h2 = plt.hist(allDataR2, 329, density=True, color="#212529"); 
##plt.xlabel("Position"); 
##plt.ylabel("Frequency of mutations"); 
##plt.title("Region 2"); 
###fig.savefig(picsPath + "whereMutationsReg2" + ".pdf", bbox_inches = "tight"); 
###fig.savefig(picsPath + "whereMutationsReg2" + ".svg", bbox_inches = "tight"); 
##
##fig = plt.figure(); 
##h3 = plt.hist(allDataR3, 275, density=True, color="#212529"); 
##plt.xlabel("Position"); 
##plt.ylabel("Frequency of mutations"); 
##plt.title("Region 3"); 
###fig.savefig(picsPath + "whereMutationsReg3" + ".pdf", bbox_inches = "tight"); 
###fig.savefig(picsPath + "whereMutationsReg3" + ".svg", bbox_inches = "tight"); 


#h1 = h1[0]; 
#h2 = h2[0]; 
#h3 = h3[0]; 
#
#h1.sort(); 
#h2.sort(); 
#h3.sort(); 
#
#h1 = np.flip(h1); 
#h2 = np.flip(h2); 
#h3 = np.flip(h3); 
#
#fig = plt.figure(); 
#plt.plot(h1); 
#plt.plot(h2); 
#plt.plot(h3); 
#fig.savefig(picsPath + "cumulNMutations.pdf"); 
#
#fig = plt.figure(); 
#ax = fig.add_subplot(111); 
#plt.plot(h1); 
#plt.plot(h2); 
#plt.plot(h3); 
#ax.set_xscale("log"); 
#ax.set_yscale("log"); 
#fig.savefig(picsPath + "cumulNMutations_loglog.pdf"); 
#
#fig = plt.figure(); 
#ax = fig.add_subplot(111); 
#plt.plot(h1); 
#plt.plot(h2); 
#plt.plot(h3); 
## ax.set_xscale("log"); 
#ax.set_yscale("log"); 
#fig.savefig(picsPath + "cumulNMutations_linlog.pdf"); 


# Load data from the null model
fInR1 = os.path.join(dataPathInR1Model, "freqSMutations_A2_Reg1" + ".csv")
freqMutationsR1 = h.loadFreqMutations(fInR1)

fInR2_1 = os.path.join(dataPathInR2Model, "freqSMutations_A1_Reg2" + ".csv")
freqMutationsR2_1 = h.loadFreqMutations(fInR2_1)
fInR2_2 = os.path.join(dataPathInR2Model, "freqSMutations_Replicase_Reg2" + ".csv")
freqMutationsR2_2 = h.loadFreqMutations(fInR2_2)


fInR3 = os.path.join(dataPathInR3Model, "freqSMutations_Replicase_Reg3" + ".csv")
freqMutationsR3 = h.loadFreqMutations(fInR3)


# Iterate over the data
mutationsInThisCodonPostition_r1Model = np.zeros(3); 
for position, times in freqMutationsR1.items():
	if position <= 263: 
		iInCodon = position%3; 
		mutationsInThisCodonPostition_r1Model[iInCodon] += times
print("Total mutations as prob Model", np.sum(mutationsInThisCodonPostition_r1Model))

mutationsInThisCodonPostition_r2Model = np.zeros(3); 
for position, times in freqMutationsR2_1.items():
	if position <= 188: 
		iInCodon = position%3; 
		mutationsInThisCodonPostition_r2Model[iInCodon] += times

for position, times in freqMutationsR2_2.items():
	if position >= 207: 
		iInCodon = position%3; 
		mutationsInThisCodonPostition_r2Model[iInCodon] += times

mutationsInThisCodonPostition_r3Model = np.zeros(3); 
pos = 0
for position, times in freqMutationsR3.items():
	if position >= 1:
		iInCodon = pos%3; 
		mutationsInThisCodonPostition_r3Model[iInCodon] += times
		pos += 1




mutationsInThisCodonPostition_r1 = np.zeros(3); 
totalMuts = 0
for (iNuc, thisNuc) in enumerate(allDataR1):
	# Exclue positions that are not in the ORF 
	if thisNuc not in list(range(264, 271)):
		iInCodon = int(thisNuc)%3;
		mutationsInThisCodonPostition_r1[iInCodon] += 1; 
		totalMuts += 1
# Normalize as a probability distribution
mutationsInThisCodonPostition_r1Prob = mutationsInThisCodonPostition_r1/totalMuts
print("Total mutations as prob", np.sum(mutationsInThisCodonPostition_r1Prob),
		"total muts", totalMuts)

mutationsInThisCodonPostition_r2 = np.zeros(3); 
totalMuts = 0
for (iNuc, thisNuc) in enumerate(allDataR2):
	# Exclue positions that are not in the ORF 
	if thisNuc not in list(range(189, 207)):
		iInCodon = int(thisNuc)%3; 
		mutationsInThisCodonPostition_r2[iInCodon] += 1; 
		totalMuts += 1
# Normalize as a probability distribution
mutationsInThisCodonPostition_r2Prob = mutationsInThisCodonPostition_r2/totalMuts
print("Total mutations as prob", np.sum(mutationsInThisCodonPostition_r2Prob))

mutationsInThisCodonPostition_r3 = np.zeros(3); 
totalMuts = 0
for (iNuc, thisNuc) in enumerate(allDataR3): 
	# Exclue positions that are not in the ORF 
	if thisNuc != 0:
		iInCodon = int(thisNuc)%3; 
		mutationsInThisCodonPostition_r3[iInCodon] += 1; 
		totalMuts += 1
# Normalize as a probability distribution
mutationsInThisCodonPostition_r3Prob = mutationsInThisCodonPostition_r3/totalMuts
print("Total mutations as prob", np.sum(mutationsInThisCodonPostition_r3Prob))

fig = plt.figure(); 
plt.subplot(311); 
plt.plot(mutationsInThisCodonPostition_r1Prob, 'x-k'); 
plt.ylim(0.2, 0.4)
plt.subplot(312); 
plt.plot(mutationsInThisCodonPostition_r2Prob, 'x-k'); 
plt.ylim(0.2, 0.4)
plt.subplot(313); 
plt.plot(mutationsInThisCodonPostition_r3Prob, 'x-k');
plt.ylim(0.2, 0.41)
fig.savefig(picsPath + "codonsMutEmp" + ".pdf", bbox_inches = "tight"); 
fig.savefig(picsPath + "codonsMutEmp" + ".svg", bbox_inches = "tight");  
plt.show()

