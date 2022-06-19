# This script will generate a plot of the betweenness centrality

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper_plots as h;


# Experimental conditions

# Define paths
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Scripts/Topology/Test/Results/Betweenness/";
dataPathOut = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Scripts/Topology/Test/Results/Betweenness/"
suffix = "7000Nodes"

################################
# Generate plots
print("\nBetweenness centrality plots\n")

# Load file
fIn = os.path.join(dataPathIn, "Bewteenness10899_" + suffix + ".csv")
BetCentralityDegree = h.RetrieveBetweennessCentralityDegree(fIn)
centrality = np.array(list(BetCentralityDegree.keys()))
degree = np.array(list(BetCentralityDegree.values()))


# Do plots
fig = plt.figure()
plt.scatter(centrality, degree, c = degree, alpha=0.7, edgecolor='black', cmap = "plasma", s = degree)
plt.colorbar()
plt.xlabel("Betweennesss Centrality")
plt.ylabel("K")
plt.title("K vs Betweennesss Centrality (lin-lin)")
fig.savefig(dataPathOut + "BetCentrality_linlin_" + suffix + ".eps", bbox_inches = "tight")

fig = plt.figure()
plt.scatter(centrality, degree, c = degree, alpha=0.7, edgecolor='black', cmap = "plasma", s = degree)
plt.colorbar()
plt.yscale("log", nonpositive='clip')
plt.xlabel("Betweennesss Centrality")
plt.ylabel("K")
plt.title("K vs Betweennesss Centrality (log-lin)")
fig.savefig(dataPathOut + "BetCentrality_loglin_" + suffix + ".eps", bbox_inches = "tight")

fig = plt.figure()
plt.scatter(centrality, degree, c = degree, alpha=0.7, edgecolor='black', cmap = "plasma", s = degree)
plt.colorbar()
plt.xscale("log", nonpositive='clip')
# If not all the points are displayed se the xlim manually
plt.xlim(1e-5, 1)
plt.xlabel("Betweennesss Centrality")
plt.ylabel("K")
plt.title("K vs Betweennesss Centrality (lin-log)")
fig.savefig(dataPathOut + "BetCentrality_linlog_" + suffix + ".eps", bbox_inches = "tight")

fig = plt.figure()
plt.scatter(centrality, degree, c = degree, alpha=0.7, edgecolor='black', cmap = "plasma", s = degree)
plt.xscale("log", nonpositive='clip')
plt.yscale("log", nonpositive='clip')
# If not all the points are displayed se the xlim manually
plt.xlim(1e-5, 1)
plt.colorbar()
plt.xlabel("Betweennesss Centrality")
plt.ylabel("K")
plt.title("K vs Betweennesss Centrality (log-log)")
fig.savefig(dataPathOut + "BetCentrality_loglog_" + suffix + ".eps", bbox_inches = "tight")

print("\nAll figures saved!")
plt.show()




















