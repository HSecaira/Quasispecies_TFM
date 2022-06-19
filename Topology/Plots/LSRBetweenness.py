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

colors = ["#ADE8F4", "#00B4D8", "#0077B6", "#023E8A"]
markers = ["^", "s", "p", "H"]

# Iterate over each expCond
for experimentalCondition in expConds:
	# Load path names
	dataPathIn__ = dataPathIn_ + experimentalCondition + "/" 
	dataPathOut__ = dataPathIn_ + "Plots/"

##########################################################################
	# Load data
	suffix = "GrandNetwork"

	# Generate plots
	print("\nLoading data... from", experimentalCondition, "\n")

	print("\tBetweenness centrality\n")
	fIn4 = os.path.join(dataPathIn__ + "ResultsBetweenness/", "Bewteenness_" + suffix + ".csv")
	# Retrieve the betweenness centrality
	BetCentrality = h.RetrieveBetweennessCentralityDegree(fIn4)
	# Append to lists
	BetDegrees = np.array(list(BetCentrality.values()))
	BetCent = np.array(list(BetCentrality.keys()))

	# LSR fitting of the central region
	bBetGN, wBetGN, ESwBetGN = hlog.calculateParametersLSR(np.log10(BetDegrees),
									np.log10(BetCent))
	# Get regression lines
	Bet_fitGN = 10 ** (bBetGN) * BetDegrees** wBetGN

	print("\t w:", wBetGN, "+-", ESwBetGN)

	fig = plt.figure(figsize=(9, 5))
	plt.scatter(BetDegrees, BetCent, marker = "o",
			alpha=0.7, edgecolor='black', facecolor = "blue", linewidth = 0.2, s=85)

	plt.plot(BetDegrees, Bet_fitGN, linestyle="dotted", color = "black")
	
	plt.xscale("log", nonpositive='clip')
	plt.yscale("log", nonpositive='clip')
	plt.ylabel("$B(i)$")
	plt.xlabel("$k$")
	plt.title(str(experimentalCondition))




###############################################################

print("\nGrand Network")
fInGN4 = os.path.join(dataPathInGN + "ResultsBetweenness/", "Bewteenness_" + suffix + ".csv")
# Retrieve the betweenness centrality
BetCentralityGN = h.RetrieveBetweennessCentralityDegree(fInGN4)
# Append to lists
BetDegreesGN = np.array(list(BetCentralityGN.values()))
BetCentGN = np.array(list(BetCentralityGN.keys()))

# LSR fitting of the central region
bBetGN, wBetGN, ESwBetGN = hlog.calculateParametersLSR(np.log10(BetDegreesGN),
								np.log10(BetCentGN))
# Get regression lines
Bet_fitGN = 10 ** (bBetGN) * BetDegreesGN ** wBetGN
#Bet_fitGN = 10 ** (bBetGN) * BetDegreesGN ** 1.3
print("\tGN w:", wBetGN, "+-", ESwBetGN)	

fig = plt.figure(figsize=(9, 5))

plt.scatter(BetDegreesGN, BetCentGN, marker = "o",
		alpha=0.7, edgecolor='black', facecolor = "blue", linewidth = 0.2, s=85)

plt.plot(BetDegreesGN, Bet_fitGN, linestyle="dotted", color = "black")

plt.xscale("log", nonpositive='clip')
plt.yscale("log", nonpositive='clip')
plt.ylabel("$B(i)$")
plt.xlabel("$k$")
plt.title("GN")

plt.show()
	







