# This script will generate a plot of the degree distribution of nodes
# Plots must be saved as .svg instead of .eps if you want to preserve transparency

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper_plots as h;
import helper_logBinning as hlog;
import brewer2mpl;
import matplotlib.cm as cm;



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


#############################################################

# Load experimental conditions
fInExpConds = os.path.join(dataPathIn_, "filesNames40C" + ".csv")

# Load the names of the experimental conditions
expConds = h.loadExperimentalConditions(fInExpConds)

# Lists to save the data loaded
degreeVertices = []
countDegrees = []
KClust = []
KDegrees = []
pCKws = []
pCKESws = []
KNN = []
KNNDegrees = []
pKNNws = []
pKNNESws = []
BetDegrees = []
BetCent = []
EivDegrees = []
EivCent = []

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
	
	print("\tDegree distribution\n")
	fIn = os.path.join(dataPathIn__ + "ResultsTopology/", "Nodes_and_Degrees_" + suffix + ".csv");
	
	# Retrieve the degrees
	degrees = h.RetrieveDegrees(fIn);
	# Get the unique counts for each degree
	degreeVertex, countDegree = np.unique(degrees, return_counts= True);
	# Append to a list
	degreeVertices.append(degreeVertex)
	countDegrees.append(countDegree/degrees.shape[0])


	print("\tKClustering distribution\n")
	fIn2 = os.path.join(dataPathIn__ + "ResultsTopology/", "KClusteringCoefficient_" + suffix + ".csv")
	# Retrieve the K clustering coefficient
	KClustCoeff = h.RetrieveKClustering(fIn2)
	# Append to list
	tempKClust = np.array(list(KClustCoeff.values()))
	tempKDegrees = np.array(list(KClustCoeff.keys()))
	KClust.append(tempKClust)
	KDegrees.append(tempKDegrees)	
	# LSR fitting of the central region
	bClust, wClust, ESwClust = hlog.calculateParametersLSR(np.log10(tempKDegrees[(tempKDegrees >= 4) & (tempKDegrees <= 55)]),
			np.log10(tempKClust[(tempKDegrees >= 4) & (tempKDegrees <= 55)]))
	# Append to list
	pCKws.append(wClust)
	pCKESws.append(ESwClust)

	print("\tKNearesNeighbors distribution\n")
	fIn3 = os.path.join(dataPathIn__ + "ResultsTopology/", "KNearestNeighborDegree_" + suffix + ".csv");
	# Retrieve the Knearest neighbor degrees
	KNNeighbors = h.RetrieveNearestNeighborDegree(fIn3)
	tempKNN = np.array(list(KNNeighbors.values()))
	tempKNNDegrees = np.array(list(KNNeighbors.keys()))
	KNN.append(tempKNN)
	KNNDegrees.append(tempKNNDegrees)	
	# LSR fitting of the central region
	bKnn, wKnn, ESwKnn = hlog.calculateParametersLSR(np.log10(tempKNNDegrees[(tempKNNDegrees >= 6) & (tempKNNDegrees <= 200)]),
			np.log10(tempKNN[(tempKNNDegrees >= 6) & (tempKNNDegrees <= 200)]))
	# Append to list
	pKNNws.append(wKnn)
	pKNNESws.append(ESwKnn)

	print("\tBetweenness centrality\n")
	fIn4 = os.path.join(dataPathIn__ + "ResultsBetweenness/", "Bewteenness_" + suffix + ".csv")
	# Retrieve the betweenness centrality
	BetCentrality = h.RetrieveBetweennessCentralityDegree(fIn4)
	# Append to lists
	BetDegrees.append(np.array(list(BetCentrality.values())))
	BetCent.append(np.array(list(BetCentrality.keys())))

	print("\tEigenvector centrality\n")
	fIn4 = os.path.join(dataPathIn__ + "ResultsTopology/", "EigenvectorCentrality_" + suffix + ".csv")
	# Retrieve the eigenvector centrality
	EigenCentrality = h.RetrieveEigenvectorCentralityDegree(fIn4)
	# Append to lists
	EivDegrees.append(np.array(list(EigenCentrality.values())))
	EivCent.append(np.array(list(EigenCentrality.keys())))
	
##################################################
# Load data from the GrandNetwork
suffix = "GrandNetwork"
print("\tLoading data from the GrandNetwork\n")
fInGN = os.path.join(dataPathInGN + "ResultsTopology/", "Nodes_and_Degrees_" + suffix + ".csv")
# Retrieve the degrees
degreesGN = h.RetrieveDegrees(fInGN);
# Get the unique counts for each degree
degreeVertexGN, countDegreeGN = np.unique(degreesGN, return_counts= True);
countDegreeGN = countDegreeGN/degreesGN.shape[0]


fInGN2 = os.path.join(dataPathInGN + "ResultsTopology/", "KClusteringCoefficient_" + suffix + ".csv")
# Retrieve the K clustering coefficient
KClustCoeffGN = h.RetrieveKClustering(fInGN2)
# Append to list
KClustGN = np.array(list(KClustCoeffGN.values()))
KDegreesGN = np.array(list(KClustCoeffGN.keys()))
# LSR fitting of the central region
bClustGN, wClustGN, ESwClustGN = hlog.calculateParametersLSR(np.log10(KDegreesGN[(KDegreesGN >= 4) & (KDegreesGN <= 55)]),
								np.log10(KClustGN[(KDegreesGN >= 4) & (KDegreesGN <= 55)]))
# Get regression lines
CK_bin_fitGN = 10 ** (bClustGN) * KDegreesGN** wClustGN


fInGN3 = os.path.join(dataPathInGN + "ResultsTopology/", "KNearestNeighborDegree_" + suffix + ".csv");
# Retrieve the Knearest neighbor degrees
KNNeighborsGN = h.RetrieveNearestNeighborDegree(fInGN3)
# Append to list
KNNGN = np.array(list(KNNeighborsGN.values()))
KNNDegreesGN = np.array(list(KNNeighborsGN.keys()))	
# LSR fitting of the central region
bKnnGN, wKnnGN, ESwKnnGN = hlog.calculateParametersLSR(np.log10(KNNDegreesGN[(KNNDegreesGN >= 6) & (KNNDegreesGN <= 200)]),
						np.log10(KNNGN[(KNNDegreesGN >= 6) & (KNNDegreesGN <= 200)]))
# Get regression lines
Knn_bin_fitGN = 10 ** (bKnnGN) * KNNDegreesGN** wKnnGN



fInGN4 = os.path.join(dataPathInGN + "ResultsBetweenness/", "Bewteenness_" + suffix + ".csv")
# Retrieve the betweenness centrality
BetCentralityGN = h.RetrieveBetweennessCentralityDegree(fInGN4)
# Append to lists
BetDegreesGN = np.array(list(BetCentralityGN.values()))
BetCentGN = np.array(list(BetCentralityGN.keys()))	


fInGN4 = os.path.join(dataPathInGN + "ResultsTopology/", "EigenvectorCentrality_" + suffix + ".csv")
# Retrieve the eigenvector centrality
EigenCentralityGN = h.RetrieveEigenvectorCentralityDegree(fInGN4)
# Append to lists
EivDegreesGN = np.array(list(EigenCentralityGN.values()))
EivCentGN = np.array(list(EigenCentralityGN.keys()))


##################################################

# Load data from the ancestor
suffix = "GrandNetwork"
print("\tLoading data from the ancestor\n")
name = "EL4-reg1_S4_L001_001_all_trim_merged_filter_sort_filter_length_collapsed/"

fInAncestor = os.path.join(dataPathInAncestor + name + "ResultsTopology/", "Nodes_and_Degrees_" + suffix + ".csv")
# Retrieve the degrees
degreesAncestor = h.RetrieveDegrees(fInAncestor);
# Get the unique counts for each degree
degreeVertexAncestor, countDegreeAncestor = np.unique(degreesAncestor, return_counts= True);
countDegreeAncestor = countDegreeAncestor/degreesAncestor.shape[0]


fInGAncestor2 = os.path.join(dataPathInAncestor + name + "ResultsTopology/", "KClusteringCoefficient_" + suffix + ".csv")
# Retrieve the K clustering coefficient
KClustCoeffAncestor = h.RetrieveKClustering(fInGAncestor2)
# Append to list
KClustAncestor = np.array(list(KClustCoeffAncestor.values()))
KDegreesAncestor = np.array(list(KClustCoeffAncestor.keys()))
# LSR fitting of the central region
bClustAncestor, wClustAncestor, ESwClustAncestor = hlog.calculateParametersLSR(np.log10(KDegreesAncestor[(KDegreesAncestor >= 4) & (KDegreesAncestor <= 55)]),
								np.log10(KClustAncestor[(KDegreesAncestor >= 4) & (KDegreesAncestor <= 55)]))
# Get regression lines
CK_bin_fitAncestor = 10 ** (bClustAncestor) * KDegreesAncestor** wClustAncestor


fInAncestor3 = os.path.join(dataPathInAncestor + name + "ResultsTopology/", "KNearestNeighborDegree_" + suffix + ".csv");
# Retrieve the Knearest neighbor degrees
KNNeighborsAncestor = h.RetrieveNearestNeighborDegree(fInAncestor3)
# Append to list
KNNAncestor = np.array(list(KNNeighborsAncestor.values()))
KNNDegreesAncestor = np.array(list(KNNeighborsAncestor.keys()))	
# LSR fitting of the central region
bKnnAncestor, wKnnAncestor, ESwKnnAncestor = hlog.calculateParametersLSR(np.log10(KNNDegreesAncestor[(KNNDegreesAncestor >= 6) & (KNNDegreesAncestor <= 200)]),
						np.log10(KNNAncestor[(KNNDegreesAncestor >= 6) & (KNNDegreesAncestor <= 200)]))
# Get regression lines
Knn_bin_fitAncestor = 10 ** (bKnnAncestor) * KNNDegreesAncestor** wKnnAncestor


fInAncestor4 = os.path.join(dataPathInAncestor + name + "ResultsBetweenness/", "Bewteenness_" + suffix + ".csv")
# Retrieve the betweenness centrality
BetCentralityAncestor = h.RetrieveBetweennessCentralityDegree(fInAncestor4)
# Append to lists
BetDegreesAncestor = np.array(list(BetCentralityAncestor.values()))
BetCentAncestor = np.array(list(BetCentralityAncestor.keys()))	


fInAncestor5 = os.path.join(dataPathInAncestor + name + "ResultsTopology/", "EigenvectorCentrality_" + suffix + ".csv")
# Retrieve the eigenvector centrality
EigenCentralityAncestor = h.RetrieveEigenvectorCentralityDegree(fInAncestor5)
# Append to lists
EivDegreesAncestor = np.array(list(EigenCentralityAncestor.values()))
EivCentAncestor = np.array(list(EigenCentralityAncestor.keys()))

##################################################


# Plots

# Set colors for each experimental conditions
#colors = cm.get_cmap("gist_rainbow")(np.linspace(0, 1, len(expConds)))
colors = ["#ADE8F4", "#00B4D8", "#0077B6", "#023E8A"]
colors_alpha = ["white", "#ADE8F4", "#00B4D8", "#0077B6", "#023E8A"]
# Set markers for each experimental condition
markers = ["^", "s", "p", "H"]
markers_alpha = ["o", "^", "s", "p", "H"]
labels = [4, 10, 30, 60]
labels_alpha = [2, 4, 10, 30, 60]


##print("\nPlot for degree distribution")
##fig = plt.figure(figsize=(9, 5))
##ax = fig.add_subplot(111)
##
### Plot for experimental conditions
##for i in range(len(expConds)):
	##label = labels[i]
	##ax.scatter(degreeVertices[i], countDegrees[i], label = str(label), marker = markers[i], 
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth=0.2, s=85)
##
### Add plot of the Grand Network
##GrandNetwork = ax.scatter(degreeVertexGN, countDegreeGN, marker = "o",
		##alpha = 0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
##
### Add plot of the ancestor
##Ancestor = ax.scatter(degreeVertexAncestor, countDegreeAncestor, marker = "o", 
		##alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
##
##ax.set_xscale("log", nonpositive='clip')
##ax.set_yscale("log", nonpositive='clip')
##ax.set_ylim(1e-6, 1e0)
##ax.set_xlabel("$k$")
##ax.set_ylabel("$p(k)$")
###ax.set_title("Degree Distribution (40 C)")
####plt.subplots_adjust(right=0.8)
### Add first legend: only labeled data
##leg1 = ax.legend(bbox_to_anchor=(0.88,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
### Add second legend
##leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(0.25,0), loc="lower right", ncol = 1, fancybox=True,
									##title="Network")
### Add the first legend back
##ax.add_artist(leg1)
### Save
##fig.savefig(dataPathOut_ + "Degree_distribution" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "Degree_distribution" + ".pdf", bbox_inches = "tight")
##


### Plot for KClustering distribution
##print("\nPlot for KCLusterin distribution")
##fig = plt.figure(figsize=(9, 5))
##ax = fig.add_subplot(111)
##for i in range(len(expConds)):
	##label = labels[i]
	##ax.scatter(KDegrees[i], KClust[i], label = str(label), marker = markers[i],
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth=0.2, s=85)
##
### Add plot of the grand network
##GrandNetwork = ax.scatter(KDegreesGN, KClustGN, marker="o",
		##alpha=0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
### Regression line
##ax.plot(KDegreesGN, CK_bin_fitGN, linestyle="dashed", color = "#6c757d")
##print("\tGrandNetwork")
##print("\tw: ", wClustGN, "+-:", ESwClustGN)
##
### Add plot of the ancestor
##Ancestor = ax.scatter(KDegreesAncestor, KClustAncestor, marker="o",
		##alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
##ax.plot(KDegreesAncestor, CK_bin_fitAncestor, linestyle="dotted", color = "#495057")
##print("\nAncestor")
##print("\tw: ", wClustAncestor, "+-:", ESwClustAncestor)
### Add alpha of the ancestor to the list of alphas
##pCKws = [wClustAncestor, *pCKws]
##pCKESws = [ESwClustAncestor, *pCKESws]
### Set limits mannually, otherwise points are not correctly plotted
##ax.set_xlim(1, 1e3)
##ax.set_ylim(1e-3, np.amax(KClust[i]) + 0.1)
##ax.set_xscale("log", nonpositive='clip')
##ax.set_yscale("log", nonpositive='clip')
##ax.set_xlabel("$k$")
##ax.set_ylabel("$C(k)$")
####plt.subplots_adjust(right=0.8)
##### Add first legend: only labeled data
####leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
##### Add second legend
####leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
##### Add the first legend back
####ax.add_artist(leg1)
### Save
##fig.savefig(dataPathOut_ + "KCustering_distribution" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "KCustering_distribution" + ".pdf", bbox_inches = "tight")
### Plot for estimated coefficients
##fig = plt.figure(figsize=(9, 5))
##for i in range(len(labels_alpha)):
	##plt.scatter(labels_alpha[i], pCKws[i], marker = markers_alpha[i], color = colors_alpha[i], alpha = 1,  s=85,
				##edgecolor="black")
##plt.plot(labels_alpha, np.array(pCKws), color = "black")
##plt.fill_between(labels_alpha, np.array(pCKws) - np.array(pCKESws), np.array(pCKws) + np.array(pCKESws),
				##alpha = 0.1, color = "black")
### Change ylims manually for a nice plot
##plt.ylim(-1.75, -1.45)
##plt.xlabel("Transfers")
##plt.ylabel("$\\hat{\\alpha}_{SN}$")
##fig.savefig(dataPathOut_ + "KCustering_distributionAlpha" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "KCustering_distributionAlpha" + ".pdf", bbox_inches = "tight")
##
##
##
##
##
##
#### Plot for KNN distribution
##print("\nPlot for KNN distribution")
##fig = plt.figure(figsize=(9, 5))
##ax = fig.add_subplot(111)
##for i in range(len(expConds)):
	##label = labels[i]
	###label = expConds[i].strip().split("_")[1]
	##ax.scatter(KNNDegrees[i], KNN[i], label = str(label), marker = markers[i],
		##alpha=0.7, edgecolor="black", facecolor = colors[i], linewidth = 0.2, s=60)
##
### Add plot of the grand network
##GrandNetwork = ax.scatter(KNNDegreesGN, KNNGN, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
### Regression line
##ax.plot(KNNDegreesGN, Knn_bin_fitGN, linestyle="dashed", color = "#6c757d")
##print("\tGrandNetwork")
##print("\tw: ", wKnnGN, "+-:", ESwKnnGN)
##
### Add plot of the ancestor
##Ancestor = ax.scatter(KNNDegreesAncestor, KNNAncestor, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
### Regression line
##ax.plot(KNNDegreesAncestor, Knn_bin_fitAncestor, linestyle="dotted", color = "#495057")
##print("\nAncestor")
##print("\tw: ", wKnnAncestor, "+-:", ESwKnnAncestor)
### Add alpha of the ancestor to the list of alphas
###labels = [2, *labels]
##pKNNws = [wKnnAncestor, *pKNNws]
##pKNNESws = [ESwKnnAncestor, *pKNNESws]
##print(pKNNws, pKNNESws)
##
##ax.set_xscale("log", nonpositive='clip')
##ax.set_yscale("log", nonpositive='clip')
##ax.set_ylim(2, 1e2 + 20)
##ax.set_xlabel("$k$")
##ax.set_ylabel("$K_{nn}(k)$")
####plt.subplots_adjust(right=0.8)
##### Add first legend: only labeled data
####leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
##### Add second legend
####leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
##### Add the first legend back
####ax.add_artist(leg1)
### Save
##fig.savefig(dataPathOut_ + "KNN_distribution" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "KNN_distribution" + ".pdf", bbox_inches = "tight")
### Plot for estimated coefficients
##fig = plt.figure(figsize=(9, 5))
##for i in range(len(labels_alpha)):
	##plt.scatter(labels_alpha[i], pKNNws[i], marker = markers_alpha[i], color = colors_alpha[i], alpha = 1,  s=85,
				##edgecolor="black")
##plt.plot(labels_alpha, np.array(pKNNws), color = "black")
##plt.fill_between(labels_alpha, np.array(pKNNws) - np.array(pKNNESws), np.array(pKNNws) + np.array(pKNNESws),
				##alpha = 0.1, color = "black")
### Change ylims manually for a nice plot
##plt.ylim(-0.85, -0.65)
##plt.xlabel("Transfers")
##plt.ylabel("$\\hat{\\alpha}_{SN}$")
##fig.savefig(dataPathOut_ + "KNN_distributionAlpha" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "KNN_distributionAlpha" + ".pdf", bbox_inches = "tight")


#### Plot for Betweenness Centrality
##print("\nPlot for Betweenness centrality")
##fig = plt.figure(figsize=(9, 5))
##ax = fig.add_subplot(111)
##for i in range(len(expConds)):
	##label = labels[i]
	###label = expConds[i].strip().split("_")[1]
	##ax.scatter(BetDegrees[i], BetCent[i], label = str(label), marker = markers[i],
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2, s=85)
##
### Add plot of grand network
##GrandNetwork = ax.scatter(BetDegreesGN, BetCentGN, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
##
### Add plot of grand network
##Ancestor = plt.scatter(BetDegreesAncestor, BetCentAncestor, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
##
##ax.set_ylabel("$B(i)$")
##ax.set_xlabel("$k$")
###ax.set_title("Betweennesss Centrality (40 C)")
####plt.subplots_adjust(right=0.8)
##### Add first legend: only labeled data
####leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
##### Add second legend
####leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
##### Add the first legend back
####ax.add_artist(leg1)
### Save
##fig.savefig(dataPathOut_ + "BetweenessCentrality_linlin" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_linlin" + ".pdf", bbox_inches = "tight")
##
##
##fig = plt.figure(figsize=(9, 5))
##ax = fig.add_subplot(111)
##for i in range(len(expConds)):
	##label = labels[i]
	##ax.scatter(BetDegrees[i], BetCent[i], label = str(label), marker = markers[i],
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2, s=85)
##
### Add plot of grand network
##GrandNetwork = ax.scatter(BetDegreesGN, BetCentGN, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
##
### Add plot of the ancestor
##Ancestor = ax.scatter(BetDegreesAncestor, BetCentAncestor, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
##
##ax.set_yscale("log", nonpositive='clip')
##ax.set_ylabel("$B(i)$")
##ax.set_xlabel("$k$")
###ax.set_title("Betweennesss Centrality (40 C)")
####plt.subplots_adjust(right=0.8)
##### Add first legend: only labeled data
####leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
##### Add second legend
####leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
##### Add the first legend back
####ax.add_artist(leg1)
### Save
##fig.savefig(dataPathOut_ + "BetweenessCentrality_loglin" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_loglin" + ".pdf", bbox_inches = "tight")
##
##
##fig = plt.figure(figsize=(9, 5))
##ax = fig.add_subplot(111)
##for i in range(len(expConds)):
	##label = labels[i]
	##ax.scatter(BetDegrees[i], BetCent[i], label = str(label), marker = markers[i],
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2, s=85)
##	
### Add plot of grand network
##GrandNetwork = ax.scatter(BetDegreesGN, BetCentGN, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=60)
##
### Add plot of the ancestor
##Ancestor = ax.scatter(BetDegreesAncestor, BetCentAncestor, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=60)
##
##ax.set_xscale("log", nonpositive='clip')
##ax.set_ylabel("$B(i)$")
##ax.set_xlabel("$k$")
###ax.set_title("Betweennesss Centrality (40 C)")
####plt.subplots_adjust(right=0.8)
##### Add first legend: only labeled data
####leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
##### Add second legend
####leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
##### Add the first legend back
####ax.add_artist(leg1)
### Save
##fig.savefig(dataPathOut_ + "BetweenessCentrality_linlog" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_linlog" + ".pdf", bbox_inches = "tight")
##
##
##fig = plt.figure(figsize=(9, 5))
##ax = fig.add_subplot(111)
##for i in range(len(expConds)):
	##label = labels[i]
	###label = expConds[i].strip().split("_")[1]
	##ax.scatter(BetDegrees[i], BetCent[i], label = str(label), marker = markers[i],
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2, s=85)
##
### Add plot of grand network
##GrandNetwork = ax.scatter(BetDegreesGN, BetCentGN, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
##
### Add plot of the ancestor
##Ancestor = ax.scatter(BetDegreesAncestor, BetCentAncestor, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
##
##ax.set_xscale("log", nonpositive='clip')
##ax.set_yscale("log", nonpositive='clip')
##ax.set_ylabel("$B(i)$")
##ax.set_xlabel("$k$")
###ax.set_title("Betweennesss Centrality (40 C)")
####plt.subplots_adjust(right=0.8)
##### Add first legend: only labeled data
####leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
##### Add second legend
####leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
##### Add the first legend back
####ax.add_artist(leg1)
### Save
##fig.savefig(dataPathOut_ + "BetweenessCentrality_loglog" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_loglog" + ".pdf", bbox_inches = "tight")


#### Plot for Eigenvector Centrality
##print("\nPlot for Eigenvector centrality")
##fig = plt.figure(figsize=(9, 5))
##ax = fig.add_subplot(111)
##for i in range(len(expConds)):
	##label = labels[i]
	###label = expConds[i].strip().split("_")[1]
	##ax.scatter(EivDegrees[i], EivCent[i], label = str(label), marker = markers[i],
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2, s=85)
##
### Add plot of grand network
##ax.scatter(EivDegreesGN, EivCentGN, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
##
### Add plot of the ancestor
##ax.scatter(EivDegreesAncestor, EivCentAncestor, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
##
##ax.set_ylabel("$v_{1}(i)$")
##ax.set_xlabel("$k$")
###ax.set_title("Eigenvector Centrality (40 C)")
####plt.subplots_adjust(right=0.8)
##### Add first legend: only labeled data
####leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
##### Add second legend
####leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
##### Add the first legend back
####ax.add_artist(leg1)
### Save
##fig.savefig(dataPathOut_ + "EigenCentrality_linlin_" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "EigenCentrality_linlin_" + ".pdf", bbox_inches = "tight")
##
##
##fig = plt.figure(figsize=(9, 5))
##ax = fig.add_subplot(111)
##for i in range(len(expConds)):
	##label = labels[i]
	###label = expConds[i].strip().split("_")[1]
	##ax.scatter(EivDegrees[i], EivCent[i], label = str(label), marker = markers[i],
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2, s=85)
##
### Add plot of grand network
##ax.scatter(EivDegreesGN, EivCentGN, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
##
### Add plot of the ancestor
##ax.scatter(EivDegreesAncestor, EivCentAncestor, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
##
##ax.set_yscale("log", nonpositive='clip')
##ax.set_ylabel("$v_{1}(i)$")
##ax.set_xlabel("$k$")
###ax.set_title("Eigenvector Centrality (40 C)")
####plt.subplots_adjust(right=0.8)
##### Add first legend: only labeled data
####leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
##### Add second legend
####leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
##### Add the first legend back
####ax.add_artist(leg1)
### Save
##fig.savefig(dataPathOut_ + "EigenCentrality_loglin_" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "EigenCentrality_loglin_" + ".pdf", bbox_inches = "tight")
##
##
##fig = plt.figure(figsize=(9, 5))
##ax = fig.add_subplot(111)
##for i in range(len(expConds)):
	##label = labels[i]
	###label = expConds[i].strip().split("_")[1]
	##ax.scatter(EivDegrees[i], EivCent[i], label = str(label), marker = markers[i],
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2, s=85)
##
### Add plot of grand network
##ax.scatter(EivDegreesGN, EivCentGN, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
##
### Add plot of the ancestor
##ax.scatter(EivDegreesAncestor, EivCentAncestor, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
##
##ax.set_xscale("log", nonpositive='clip')
##ax.set_ylabel("$v_{1}(i)$")
##ax.set_xlabel("$k$")
###ax.set_title("Eigenvector Centrality (40 C)")
####plt.subplots_adjust(right=0.8)
##### Add first legend: only labeled data
####leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
##### Add second legend
####leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
##### Add the first legend back
####ax.add_artist(leg1)
### Save
##fig.savefig(dataPathOut_ + "EigenCentrality_linlog_" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "EigenCentrality_linlog_" + ".pdf", bbox_inches = "tight")
##
##
fig = plt.figure(figsize=(9, 5))
ax = fig.add_subplot(111)
for i in range(len(expConds)):
	label = labels[i]
	plt.scatter(EivDegrees[i], EivCent[i], label = str(label), marker = markers[i],
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2, s=85)

### Add plot of grand network
##ax.scatter(EivDegreesGN, EivCentGN, marker = "o",
		##alpha=0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)

# Add plot of the ancestor
ax.scatter(EivDegreesAncestor, EivCentAncestor, marker = "o",
		alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
	
ax.set_xscale("log", nonpositive='clip')
ax.set_yscale("log", nonpositive='clip')
ax.set_ylabel("$v_{1}(i)$")
ax.set_xlabel("$k$")
#ax.set_title("Eigenvector Centrality (40 C)")
##plt.subplots_adjust(right=0.8)
### Add first legend: only labeled data
##leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
### Add second legend
##leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
### Add the first legend back
##ax.add_artist(leg1)
# Save
#fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + ".svg", bbox_inches = "tight")
#fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + ".pdf", bbox_inches = "tight")


print("\nAll plots saved!")


plt.show()