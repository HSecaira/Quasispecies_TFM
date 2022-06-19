"""
Script that performs logarithmic binning from files of a given temperature
"""

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import matplotlib.cm as cm;
import helper_plots as h;
import helper_logBinning as hlog;

# Define paths
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/";
dataPathInGN = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/GrandNetwork/Reg1/";
dataPathInAncestor = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/Reg1_experiments/SubNetworks/";
#dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/SubnetworksConstantT/Reg1/";

# Define regions
dicFolder = {};
dicFolder["r1"] = "Reg1_experiments/";
dicFolder["r2"] = "Reg2_experiments/";
dicFolder["r3"] = "Reg3_experiments/";
regionName = "r1";

dataPathIn_ = dataPathIn + dicFolder[regionName] + "SubNetworks/"
#dataPathIn_ = dataPathIn + "SubnetworksConstantT/" + "Reg1/"
dataPathOut_ = dataPathIn_ + "Plots40C/"

# Load experimental conditions
fInExpConds = os.path.join(dataPathIn_, "filesNames40C" + ".csv")

# Load the names of the experimental conditions
expConds = h.loadExperimentalConditions(fInExpConds)

# For plots
colors = ["#ADE8F4", "#00B4D8", "#0077B6", "#023E8A"]
markers = ["^", "s", "p", "H"]
labels = [4, 10, 30, 60]

# Iterate over each expCond
i = 0
for experimentalCondition in expConds:
	# Load path names
	dataPathIn__ = dataPathIn_ + experimentalCondition + "/" 
	dataPathOut__ = dataPathIn_ + "Plots/"

##########################################################################
	# Load data
	suffix = "GrandNetwork"

	# Generate plots
	print("\nLoading data... from", experimentalCondition, "\n")

	print("\tEigenvector centrality\n")
	fIn4 = os.path.join(dataPathIn__ + "ResultsTopology/", "EigenvectorCentrality_" + suffix + ".csv")
	# Retrieve the eigenvector centrality
	EigenCentrality = h.RetrieveEigenvectorCentralityDegree(fIn4)
	# Append to lists
	EivDegrees = np.array(list(EigenCentrality.values()))
	EivCent = np.array(list(EigenCentrality.keys()))

	# LSR fitting of the central region
	bEivGN, wEivGN, ESwEivGN = hlog.calculateParametersLSR(np.log10(EivDegrees),
									np.log10(EivCent))
	# Get regression lines
	Eiv_fitGN = 10 ** (bEivGN) * EivDegrees ** wEivGN

	print("\t w:", wEivGN, "+-", ESwEivGN)

	fig = plt.figure(figsize=(9, 5))
	plt.scatter(EivDegrees, EivCent, marker = markers[i],
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2, s=85)
	
	plt.xscale("log", nonpositive='clip')
	plt.yscale("log", nonpositive='clip')
	plt.ylabel("$v_{1}(i)$")
	plt.xlabel("$k$")
	#plt.title(str(labels[i]))
	fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + str(labels[i]) + ".svg", bbox_inches = "tight")
	fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + str(labels[i]) + ".pdf", bbox_inches = "tight")

	i += 1


###############################################################

print("Grand Network")
suffix = "GrandNetwork"

fInGN4 = os.path.join(dataPathInGN + "ResultsTopology/", "EigenvectorCentrality_" + suffix + ".csv")
# Retrieve the eigenvector centrality
EigenCentralityGN = h.RetrieveEigenvectorCentralityDegree(fInGN4)
# Append to lists
EivDegreesGN = np.array(list(EigenCentralityGN.values()))
EivCentGN = np.array(list(EigenCentralityGN.keys()))

fig = plt.figure(figsize=(9, 5))
plt.scatter(EivDegreesGN, EivCentGN, marker = "o",
	alpha=0.5, edgecolor='black', facecolor = "gray", linewidth = 0.6, s=50)	
plt.xscale("log", nonpositive='clip')
plt.yscale("log", nonpositive='clip')
plt.ylabel("$v_{1}(i)$")
plt.xlabel("$k$")
#plt.title("GN")
fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + "GN" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + "GN" + ".pdf", bbox_inches = "tight")




print("Ancestor")
suffix = "GrandNetwork"
name = "EL4-reg1_S4_L001_001_all_trim_merged_filter_sort_filter_length_collapsed/"

fInAncestor5 = os.path.join(dataPathInAncestor + name + "ResultsTopology/", "EigenvectorCentrality_" + suffix + ".csv")
# Retrieve the eigenvector centrality
EigenCentralityAncestor = h.RetrieveEigenvectorCentralityDegree(fInAncestor5)
# Append to lists
EivDegreesAncestor = np.array(list(EigenCentralityAncestor.values()))
EivCentAncestor = np.array(list(EigenCentralityAncestor.keys()))

fig = plt.figure(figsize=(9, 5))
plt.scatter(EivDegreesAncestor, EivCentAncestor, marker = "o",
	alpha=0.5, edgecolor='black', facecolor = "white", linewidth = 0.6, s=50)	
plt.xscale("log", nonpositive='clip')
plt.yscale("log", nonpositive='clip')
plt.ylabel("$v_{1}(i)$")
plt.xlabel("$k$")
#plt.title("Ancestor")
fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + "Ancestor" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + "Ancestor" + ".pdf", bbox_inches = "tight")

print("DOne!")
#plt.show()
	







