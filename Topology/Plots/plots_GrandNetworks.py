# This script will generate a plot of the degree distribution of nodes
# Plots must be saved as .svg instead of .eps if you want to preserve transparency

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper_plots as h;
import brewer2mpl;
import matplotlib.cm as cm;


# Define paths
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/GrandNetwork/";
#dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/SubnetworksConstantT/Reg1/";

# Define regions

dataPathIn_ = dataPathIn
dataPathOut_ = dataPathIn + "Plots/"


#############################################################

# Load experimental conditions
fInExpConds = os.path.join(dataPathIn_, "filesNames" + ".csv")

# Load the names of the experimental conditions
expConds = h.loadExperimentalConditions(fInExpConds)

# Lists to save the data loaded
degreeVertices = []
countDegrees = []
KClust = []
KDegrees = []
KNN = []
KNNDegrees = []
BetDegrees = []
BetCent = []
EivDegrees = []
EivCent = []
RCDegrees = []
RCCoeffs = []

# Iterate over each expCond
for experimentalCondition in expConds:
	# Load path names
	dataPathIn__ = dataPathIn_ + experimentalCondition + "/" 
	dataPathOut__ = dataPathIn_ + "Plots/"

##########################################################################

	# Load data
	suffix = "GrandNetwork"
	
	# Generate plots
	print("\n Prepating data... for", experimentalCondition, "\n")
	
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
	# Retrieve the degrees
	Kdegrees2 = h.RetrieveKDegrees(fIn2)
	# Retrieve the K clustering coefficient
	KClustCoeff = h.RetrieveKClustering(fIn2)
	# Create a dictionary
	DictKClustCoeff = dict(zip(Kdegrees2, KClustCoeff))
	# sort the dictionary, but this converts the dict into a list
	DictKClustCoeffSort = dict(sorted(DictKClustCoeff.items(), key = lambda x:x[0]))
	# Append to list
	KClust.append(np.array(list(DictKClustCoeffSort.values())))
	KDegrees.append(np.array(list(DictKClustCoeffSort.keys())))	

	print("\tKNearesNeighbors distribution\n")
	fIn3 = os.path.join(dataPathIn__ + "ResultsTopology/", "KNearestNeighborDegree_" + suffix + ".csv");
	# Retrieve the degrees
	Kdegrees = h.RetrieveKDegrees(fIn3)
	# Retrieve the Knearest neighbor degrees
	KNNeighbors = h.RetrieveNearestNeighborDegree(fIn3)
	
	# Create a dictionary
	DictKNNeighbors = dict(zip(Kdegrees, KNNeighbors))
	# sort the dictionary, but this converts the dict into a list
	DictKNNeighborsSort = dict(sorted(DictKNNeighbors.items(), key = lambda x:x[0]))
	# Append to list
	KNN.append(np.array(list(DictKNNeighborsSort.values())))
	KNNDegrees.append(np.array(list(DictKNNeighborsSort.keys())))	

	##print("\tBetweenness centrality\n")
	##betTemp = {}
	##BetCentralityDegree = {}
	##fIn4 = os.path.join(dataPathIn__ + "ResultsBetweenness/", "Bewteenness_" + suffix + ".csv")
	### Retrieve the betweenness centrality
	##BetCentrality = h.RetrieveBetweennessCentralityDegree(fIn4)
	### Append to lists
	##BetDegrees.append(np.array(list(BetCentrality.values())))
	##BetCent.append(np.array(list(BetCentrality.keys())))

	print("\tEigenvector centrality\n")
	eivTemp = {}
	EivCentralityDegree = {}
	fIn4 = os.path.join(dataPathIn__ + "ResultsTopology/", "EigenvectorCentrality_" + suffix + ".csv")
	# Retrieve the eigenvector centrality
	EigenCentrality = h.RetrieveEigenvectorCentralityDegree(fIn4)
	# Append to lists
	EivDegrees.append(np.array(list(EigenCentrality.values())))
	EivCent.append(np.array(list(EigenCentrality.keys())))

	##print("\tRichClubCoefficient\n")
	##fIn6 = os.path.join(dataPathIn__ + "ResultsTopology/", "RichClubCoeff_" + suffix + ".csv")
	### Retrieverich-club coefficient
	##rcCoeff = h.RetrieveRichClubCoeff(fIn6)
	### Append to lists
	##RCDegrees.append(np.array(list(rcCoeff.keys())))
	##RCCoeffs.append(np.array(list(rcCoeff.values())))
	
##################################################

# Plots

# Set colors for each experimental conditions
colors = cm.get_cmap("gist_rainbow")(np.linspace(0, 1, len(expConds)))
# Set markers for each experimental condition
markers = mStyles = [".",",","o","v","^","<",">","1","2","3","4","8","s","p","P","*","h","H","+","x","X","D","d","|","_",0,1,2,3,4,5,6,7,8,9,10,11]


## Plot for degree distribution
print("\nPlot for degree distribution")
fig = plt.figure(figsize=(9, 5))

for i in range(len(expConds)):
	label = expConds[i]
	
	plt.scatter(degreeVertices[i], countDegrees[i], label = str(label), 
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth=0.2)
	#plt.plot(degreeVertices[i], countDegrees[i], color = colors[i])
	plt.xscale("log", nonpositive='clip')
	plt.yscale("log", nonpositive='clip')
#plt.legend(loc = "best")
plt.subplots_adjust(right=0.75)
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
plt.xlabel("$k$")
plt.ylabel("$p(k)$")
plt.title("Degree Distribution")
fig.savefig(dataPathOut_ + "Degree_distribution" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "Degree_distribution" + ".pdf", bbox_inches = "tight")



## Plot for KClustering distribution
print("\nPlot for KCLusterin distribution")
fig = plt.figure(figsize=(9, 5))

for i in range(len(expConds)):
	#label = expConds[i]
	
	plt.scatter(KDegrees[i], KClust[i], label = str(label), 
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth=0.2)
	#plt.plot(KDegrees[i], KClust[i], color = colors[i])
	# Set limits mannually, otherwise points are not correctly plotted
	plt.xlim(1, 1.2e3)
	plt.ylim(1e-3, 1.8e-1)
	plt.xscale("log", nonpositive='clip')
	plt.yscale("log", nonpositive='clip')
#plt.legend(loc = "best")
plt.subplots_adjust(right=0.75)
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
plt.xlabel("$k$")
plt.ylabel("$C(k)$")

plt.title("Clustering Coefficient")
fig.savefig(dataPathOut_ + "KCustering_distribution" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "KCustering_distribution" + ".pdf", bbox_inches = "tight")



## Plot for KNN distribution
print("\nPlot for KNN distribution")
fig = plt.figure(figsize=(9, 5))

for i in range(len(expConds)):
	label = expConds[i]
	
	plt.scatter(KNNDegrees[i], KNN[i], label = str(label),
		alpha=0.7, edgecolor="black", facecolor = colors[i], linewidth = 0.2)
	# Set limits mannually, otherwise points are not correctly plotted
	plt.xscale("log", nonpositive='clip')
	plt.yscale("log", nonpositive='clip')
#plt.legend(loc = "best")
plt.subplots_adjust(right=0.75)
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
plt.xlabel("$k$")
plt.ylabel("$K_{nn}(k)$")

plt.title("K Nearest-Neighbors")
fig.savefig(dataPathOut_ + "KNN_distribution" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "KNN_distribution" + ".pdf", bbox_inches = "tight")

#### Plot for Betweenness Centrality
##print("\nPlot for Betweenness centrality")
##fig = plt.figure(figsize=(9, 5))
##
##for i in range(len(expConds)):
	##label = expConds[i]
	##
	##plt.scatter(BetDegrees[i], BetCent[i], label = str(label),
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2)
##plt.subplots_adjust(right=0.75)
##plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
##plt.ylabel("$B(i)$")
##plt.xlabel("$k$")
##
##plt.title("Betweennesss Centrality (lin-lin)")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_linlin" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_linlin" + ".pdf", bbox_inches = "tight")
##
##fig = plt.figure(figsize=(9, 5))
##
##for i in range(len(expConds)):
	##label = expConds[i]
	##
	##plt.scatter(BetDegrees[i], BetCent[i], label = str(label),
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2)
	##plt.yscale("log", nonpositive='clip')
##plt.subplots_adjust(right=0.75)
##plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
##plt.ylabel("$B(i)$")
##plt.xlabel("$k$")
##
##plt.title("Betweennesss Centrality (log-lin)")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_loglin" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_loglin" + ".pdf", bbox_inches = "tight")
##
##fig = plt.figure(figsize=(9, 5))
##
##for i in range(len(expConds)):
	##label = expConds[i]
	##
	##plt.scatter(BetDegrees[i], BetCent[i], label = str(label),
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2)
	##plt.xscale("log", nonpositive='clip')
##plt.subplots_adjust(right=0.75)
##plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
##plt.ylabel("$B(i)$")
##plt.xlabel("$k$")
##
##plt.title("Betweennesss Centrality (lin-log)")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_linlog" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_linlog" + ".pdf", bbox_inches = "tight")
##
##fig = plt.figure(figsize=(9, 5))
##
##for i in range(len(expConds)):
	##label = expConds[i]
	##
	##plt.scatter(BetDegrees[i], BetCent[i], label = str(label),
		##alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2)
	##plt.xscale("log", nonpositive='clip')
	##plt.yscale("log", nonpositive='clip')
##plt.subplots_adjust(right=0.75)
##plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
##plt.ylabel("$B(i)$")
##plt.xlabel("$k$")
##
##plt.title("Betweennesss Centrality (log-log)")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_loglog" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "BetweenessCentrality_loglog" + ".pdf", bbox_inches = "tight")


## Plot for Eigenvector Centrality
print("\nPlot for Eigenvector centrality")
fig = plt.figure(figsize=(9, 5))

for i in range(len(expConds)):
	label = expConds[i]
	
	plt.scatter(EivDegrees[i], EivCent[i], label = str(label),
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2)
plt.subplots_adjust(right=0.75)
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
plt.ylabel("$Eiv(i)$")
plt.xlabel("$k$")

plt.title("Eigenvector Centrality (lin-lin)")
fig.savefig(dataPathOut_ + "EigenCentrality_linlin_" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "EigenCentrality_linlin_" + ".pdf", bbox_inches = "tight")

fig = plt.figure(figsize=(9, 5))

for i in range(len(expConds)):
	label = expConds[i]
	
	plt.scatter(EivDegrees[i], EivCent[i], label = str(label),
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.1)
	plt.yscale("log", nonpositive='clip')
plt.subplots_adjust(right=0.75)
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
plt.ylabel("$Eiv(i)$")
plt.xlabel("$k$")

plt.title("Eigenvector Centrality (log-lin)")
fig.savefig(dataPathOut_ + "EigenCentrality_loglin_" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "EigenCentrality_loglin_" + ".pdf", bbox_inches = "tight")

fig = plt.figure(figsize=(9, 5))

for i in range(len(expConds)):
	label = expConds[i]
	
	plt.scatter(EivDegrees[i], EivCent[i], label = str(label),
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2)
	plt.xscale("log", nonpositive='clip')
plt.subplots_adjust(right=0.75)
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
plt.ylabel("$Eiv(i)$")
plt.xlabel("$k$")

plt.title("Eigenvector Centrality (lin-log)")
fig.savefig(dataPathOut_ + "EigenCentrality_linlog_" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "EigenCentrality_linlog_" + ".pdf", bbox_inches = "tight")

fig = plt.figure(figsize=(9, 5))

for i in range(len(expConds)):
	label = expConds[i]
	
	plt.scatter(EivDegrees[i], EivCent[i], label = str(label),
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth = 0.2)
	plt.xscale("log", nonpositive='clip')
	plt.yscale("log", nonpositive='clip')
plt.subplots_adjust(right=0.75)
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", ncol = 1)
plt.ylabel("$Eiv(i)$")
plt.xlabel("$k$")

plt.title("Eigenvector Centrality (log-log)")
fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + ".pdf", bbox_inches = "tight")

#plt.show()

print("\nAll plots saved!")


