# This script will generate a plot of Betweenness benchmarking

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper_plots as h;

# Define paths
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Scripts/Topology/Test/Results/Betweenness_Benchmark/";
dataPathOut = dataPathIn + "BetweennessPlots/"

##########################################################################
# Load data

suffix = "EL33-reg1"
numNodes = "544Nodes"

print("\tBetweenness centrality\n")
fIn2 = os.path.join(dataPathIn, "Bewteenness_" + numNodes + "_" + suffix + ".csv")
# Retrieve the betweenness centrality
BetCentrality = h.RetrieveBetweennessCentralityDegree(fIn2)
# Append to lists
BetDegrees = np.array(list(BetCentrality.values()))
BetCent = np.array(list(BetCentrality.keys()))

# MSE
fIn = os.path.join(dataPathIn, "Bewteenness_benchmark_" + suffix + ".csv")


##########################################################################
# Generate plots
print("\nPrepating plots...\n")

# Betweenness centrality plot
fig = plt.figure(figsize=(9, 5))
plt.scatter(BetDegrees, BetCent, c = BetDegrees, cmap = "Blues", 
			edgecolors = "black", alpha = 0.7, s=85)

plt.xscale("log", nonpositive='clip')
plt.yscale("log", nonpositive='clip')
plt.ylabel("$B(i)$")
plt.xlabel("$k$")
#plt.title("100% Nodes")
# Save
fig.savefig(dataPathOut + "BetweenessCentrality_loglin" + numNodes + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut + "BetweenessCentrality_loglin" + numNodes + ".pdf", bbox_inches = "tight")


## MSE plot
## Retrieve MSE
#fracNodes, MSE, MSESTD = h.RetrieveMSE(fIn)
#
## DO plot
#fig = plt.figure()
#plt.plot(fracNodes, MSE, marker="o", color="black", alpha = 0.7)
#plt.fill_between(fracNodes, MSE - MSESTD, MSE + MSESTD,
				#alpha = 0.1, color = "black")
#
#
#plt.xlabel("Fraction of Nodes Used")
#plt.ylabel("Mean Squared Error")
#plt.title("Betweenness Centrality Estimation")
#plt.axhline(y = 0.1, color = "red", linestyle = "dashed")
#fig.savefig(dataPathOut + "MSE_" + suffix + ".svg", bbox_inches = "tight")
#fig.savefig(dataPathOut + "MSE_" + suffix + ".pdf", bbox_inches = "tight")
#
#print("\nAll plots saved!\n")
#plt.show()



#plt.show()




	










