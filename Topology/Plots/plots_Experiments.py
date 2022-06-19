# This script will generate a plot of the degree distribution of nodes

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper_plots as h;


# Experimental conditions


# Define paths
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/Reg3_experiments/SubNetworks/";
#dataPathIn = "//home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/SubnetworksConstantT/Reg1/";
#experimentalCondition = "EL01-reg1-resec_all_all_trim_merged_filter_sort_filter_length_collapsed"

# Define regions
#dicFolder = {};
#dicFolder["r1"] = "Reg1/";
#dicFolder["r2"] = "Reg2/";
#dicFolder["r3"] = "Reg3/";
#regionName = "r3";

fInExpConds = os.path.join(dataPathIn, "filesNames" + ".csv")

# Load the names of the experimental conditions
expConds = h.loadExperimentalConditions(fInExpConds)

# Iterate over each expCond
for experimentalCondition in expConds:
	# Load path names
	dataPathIn_ = dataPathIn + experimentalCondition + "/" 
	dataPathOut_ = dataPathIn_ + "Plots/"

##########################################################################
	suffix = "GrandNetwork"
	
	# Generate plots
	print("\n Prepating plots... for", experimentalCondition, "\n")
	
	## Degree distribution
	print("\tDegree distribution plot\n")
	fIn = os.path.join(dataPathIn_ + "ResultsTopology/", "Nodes_and_Degrees_" + suffix + ".csv");
	
	# Retrieve the degrees
	degrees = h.RetrieveDegrees(fIn);
	
	# Get the unique counts for each degree
	# Unique that is degrees are the keys of the dict
	# COunts of each degree is the values of the dict
	degreeVertex, countDegree = np.unique(degrees, return_counts= True);
	
	# DO plot
	fig = plt.figure()
	plt.loglog(degreeVertex, countDegree/np.sum(countDegree), "bo-")
	plt.xlabel("K")
	plt.ylabel("P(k)")
	plt.title("Degree Distribution")
	fig.savefig(dataPathOut_ + "Degree_distribution_" + suffix + ".eps", bbox_inches = "tight")
	
	
	
	## K Nearest Neighbor degree
	print("\tKNearest Neighbor plots\n")
	fIn1 = os.path.join(dataPathIn_ + "ResultsTopology/", "KNearestNeighborDegree_" + suffix + ".csv");
	# Retrieve the degrees
	Kdegrees = h.RetrieveKDegrees(fIn1)
	# Retrieve the Knearest neighbor degrees
	KNNeighbors = h.RetrieveNearestNeighborDegree(fIn1)
	
	# Create a dictionary
	DictKNNeighbors = dict(zip(Kdegrees, KNNeighbors))
	# sort the dictionary, but this converts the dict into a list
	DictKNNeighborsSort = dict(sorted(DictKNNeighbors.items(), key = lambda x:x[0]))
	
	# Do the plot
	fig = plt.figure()
	plt.loglog(np.array(list(DictKNNeighborsSort.keys())), np.array(list(DictKNNeighborsSort.values())), "bo-")
	plt.xlabel("K")
	plt.ylabel("Knn(K)")
	plt.title("K Nearest Neighbor Distribution")
	fig.savefig(dataPathOut_ + "KNearestNeighbor_" + suffix + ".eps", bbox_inches = "tight")
	
	# K Nearest Neighbor per degree
	#fIn3 = os.path.join(dataPathIn_, "KNearestNeighbor_per_Degree_" + suffix + ".csv")
	## Retrieve the Knn per degree
	#Knn = RetrieveNearestNeighborPerDegree(fIn3)
	#
	## Do the plot
	#fig = plt.figure()
	#for key, value in Knn.items():
		#if key == 813:
			#plt.scatter([key] * len(value), value, color = 'red')
		#else:
			#plt.scatter([key] * len(value), value, color = 'blue')
	#plt.xscale("log", nonpositive='clip')
	#plt.yscale("log", nonpositive='clip')
	#plt.xlabel("K")
	#plt.ylabel("K Nearest Neighbor")
	#fig.savefig(dataPathOut_ + "KNearestNeighborPerDegree_" + suffix + ".pdf", bbox_inches = "tight")
	
	
	
	## Kclustering coefficients
	print("\tKClustering plot\n")
	fIn2 = os.path.join(dataPathIn_ + "ResultsTopology/", "KClusteringCoefficient_" + suffix + ".csv")
	# Retrieve the degrees
	Kdegrees2 = h.RetrieveKDegrees(fIn2)
	# Retrieve the K clustering coefficient
	KClustCoeff = h.RetrieveKClustering(fIn2)
	
	# Create a dictionary
	DictKClustCoeff = dict(zip(Kdegrees2, KClustCoeff))
	# sort the dictionary, but this converts the dict into a list
	DictKClustCoeffSort = dict(sorted(DictKClustCoeff.items(), key = lambda x:x[0]))
	
	# Do the plot
	fig = plt.figure()
	plt.loglog(np.array(list(DictKClustCoeffSort.keys())), np.array(list(DictKClustCoeffSort.values())), "bo-")
	plt.xlabel("K")
	plt.ylabel("C(K)")
	plt.title("Clustering Coefficient")
	fig.savefig(dataPathOut_ + "KCustering_" + suffix + ".eps", bbox_inches = "tight")
	
	
	
	## Eigenvector centrality 
	print("\tEigenvector centrality plots\n")
	fIn3 = os.path.join(dataPathIn_ + "ResultsTopology/", "EigenvectorCentrality_" + suffix + ".csv")
	
	# Retrieve the eigenvector centrality
	EigenCentrality = h.RetrieveEigenvectorCentrality(fIn3)
	
	# To plot the eigenvector centrality we must calculate the fraction of nodes having a centrality x or greater
	Eigen_centrality, count_eigen_centrality = np.unique(EigenCentrality, return_counts= True)
	
	# Do the plot
	fig = plt.figure()
	plt.loglog(Eigen_centrality, count_eigen_centrality/np.sum(count_eigen_centrality), "bo-")
	plt.xlabel("Eigenvector Centrality x")
	plt.ylabel("Fraction of vertices having a centrality x or greater")
	plt.title("Eigenvector Centrality")
	fig.savefig(dataPathOut_ + "EigenCentrality_" + suffix + ".eps", bbox_inches = "tight")
	
	# dictionary of eigenvector centrality (keys) and the degrees (values)
	ev_dict = dict(zip(EigenCentrality, degrees))
	centrality = np.array(list(ev_dict.keys()))
	degree = np.array(list(ev_dict.values()))
	
	# More plots in different scales
	fig = plt.figure()
	#plt.scatter(centrality, degree, color = "blue")
	plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
	plt.colorbar()
	plt.xlabel("Eigenvector Centrality")
	plt.ylabel("K")
	plt.title("K vs Eigenvector Centrality (lin-lin)")
	fig.savefig(dataPathOut_ + "EigenCentrality_linlin_" + suffix + ".eps", bbox_inches = "tight")
	
	fig = plt.figure()
	#plt.scatter(centrality, degree, color = "blue")
	plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
	plt.colorbar()
	plt.yscale("log", nonpositive='clip')
	plt.xlabel("Eigenvector Centrality")
	plt.ylabel("K")
	plt.title("K vs Eigenvector Centrality (log-lin)")
	fig.savefig(dataPathOut_ + "EigenCentrality_loglin_" + suffix + ".eps", bbox_inches = "tight")
	
	fig = plt.figure()
	#plt.scatter(centrality, degree, color = "blue")
	plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
	plt.colorbar()
	plt.xscale("log", nonpositive='clip')
	plt.xlabel("Eigenvector Centrality")
	plt.ylabel("K")
	plt.title("K vs Eigenvector Centrality (lin-log)")
	fig.savefig(dataPathOut_ + "EigenCentrality_linlog_" + suffix + ".eps", bbox_inches = "tight")
	
	fig = plt.figure()
	#plt.scatter(centrality, degree, color = "blue")
	plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
	plt.colorbar()
	plt.xscale("log", nonpositive='clip')
	plt.yscale("log", nonpositive='clip')
	plt.xlabel("Eigenvector Centrality")
	plt.ylabel("K")
	plt.title("K vs Eigenvector Centrality (log-log)")
	fig.savefig(dataPathOut_ + "EigenCentrality_loglog_" + suffix + ".eps", bbox_inches = "tight")
	
	fIn5 = os.path.join(dataPathIn_ + "ResultsTopology/", "EigenvectorCentrality_allComponents_" + suffix + ".csv")
	
	eigenCentralityComponents = h.RetrieveEigenvectorCentralityComponents(fIn5)
	
	for component, centralities in eigenCentralityComponents.items():
		# Get the centrality and degrees of the component
		centrality = [np.array(c[0]) for c in centralities]
		degree = [np.array(c[1]) for c in centralities]
	
		# if the components has more than two nodes
		if len(centrality) != 2:
	
			print("\tPlotting eigenvector centrality for component:", component)
	
			# Distribution
			Ev_centrality_comp, count_ev_centrality_comp = np.unique(centrality, return_counts= True)
			fig = plt.figure()
			plt.loglog(Ev_centrality_comp, count_ev_centrality_comp/np.sum(count_ev_centrality_comp), "bo-")
			plt.xlabel("Eigenvector Centrality x")
			plt.ylabel("Fraction of vertices having a centrality x or greater")
			plt.title("Eigenvector Centrality")
			fig.savefig(dataPathOut_ + "EivCentrality_" + "component_" + str(component) + "_" + suffix + ".eps", bbox_inches = "tight")
		
			# lin-lin
			fig = plt.figure()
			#plt.scatter(centrality, degree, color = "blue")
			plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
			plt.colorbar()
			plt.xlabel("Eigenvector Centrality")
			plt.ylabel("K")
			plt.title("K vs Eigenvector Centrality (lin-lin)")
			fig.savefig(dataPathOut_ + "EivCentrality_linlin_" + "component_" + str(component) + "_" + suffix + ".eps", bbox_inches = "tight")
		
			# log-lin
			fig = plt.figure()
			#plt.scatter(centrality, degree, color = "blue")
			plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
			plt.colorbar()
			plt.yscale("log", nonpositive='clip')
			plt.xlabel("Eigenvector Centrality")
			plt.ylabel("K")
			plt.title("K vs Eigenvector Centrality (log-lin)")
			fig.savefig(dataPathOut_ + "EivCentrality_loglin_" + "component_" + str(component) + "_" + suffix + ".eps", bbox_inches = "tight")
	
			# lin-log
			fig = plt.figure()
			#plt.scatter(centrality, degree, color = "blue")
			plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
			plt.colorbar()
			plt.xscale("log", nonpositive='clip')
			plt.xlabel("Eigenvector Centrality")
			plt.ylabel("K")
			plt.title("K vs Eigenvector Centrality (lin-log)")
			fig.savefig(dataPathOut_ + "EivCentrality_linlog_" + "component_" + str(component) + "_" + suffix + ".eps", bbox_inches = "tight")
	
			# log-log
			fig = plt.figure()
			#plt.scatter(centrality, degree, color = "blue")
			plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
			plt.colorbar()
			plt.xscale("log", nonpositive='clip')
			plt.yscale("log", nonpositive='clip')
			plt.xlabel("Eigenvector Centrality")
			plt.ylabel("K")
			plt.title("K vs Eigenvector Centrality (log-log)")
			fig.savefig(dataPathOut_ + "EivCentrality_loglog_" + "component_" + str(component) + "_" + suffix + ".eps", bbox_inches = "tight")
	
	
	## Betweenness centrality 
	print("\tBetweenness centrality plots\n")
	fIn4 = os.path.join(dataPathIn_ + "ResultsBetweenness/", "Bewteenness_" + suffix + ".csv")
	
	# Retrieve the betweenness centrality
	BetCentrality = h.RetrieveBetweennessCentrality(fIn4)
	
	# To plot the eigenvector centrality we must calculate the fraction of nodes having a centrality x or greater
	Bet_centrality, count_bet_centrality = np.unique(BetCentrality, return_counts= True)
	
	# Do the plot
	fig = plt.figure()
	plt.loglog(Bet_centrality, count_bet_centrality/np.sum(count_bet_centrality), "bo-")
	plt.xlabel("Betweennesss Centrality x")
	plt.ylabel("Fraction of vertices having a centrality x or greater")
	plt.title("Betweennesss Centrality")
	fig.savefig(dataPathOut_ + "BetCentrality_" + suffix + ".eps", bbox_inches = "tight")
	
	# dictionary of eigenvector centrality (keys) and the degrees (values)
	bet_dict = dict(zip(BetCentrality, degrees))
	centrality = np.array(list(bet_dict.keys()))
	degree = np.array(list(bet_dict.values()))
	
	# More plots in different scales
	fig = plt.figure()
	#plt.scatter(centrality, degree, color = "blue")
	plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
	plt.colorbar()
	plt.xlabel("Betweennesss Centrality")
	plt.ylabel("K")
	plt.title("K vs Betweennesss Centrality (lin-lin)")
	fig.savefig(dataPathOut_ + "BetCentrality_linlin_" + suffix + ".eps", bbox_inches = "tight")
	
	fig = plt.figure()
	#plt.scatter(centrality, degree, color = "blue")
	plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
	plt.colorbar()
	plt.yscale("log", nonpositive='clip')
	plt.xlabel("Betweennesss Centrality")
	plt.ylabel("K")
	plt.title("K vs Betweennesss Centrality (log-lin)")
	fig.savefig(dataPathOut_ + "BetCentrality_loglin_" + suffix + ".eps", bbox_inches = "tight")
	
	fig = plt.figure()
	#plt.scatter(centrality, degree, color = "blue")
	plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
	plt.colorbar()
	plt.xscale("log", nonpositive='clip')
	plt.xlabel("Betweennesss Centrality")
	plt.ylabel("K")
	plt.title("K vs Betweennesss Centrality (lin-log)")
	fig.savefig(dataPathOut_ + "BetCentrality_linlog_" + suffix + ".eps", bbox_inches = "tight")
	
	fig = plt.figure()
	#plt.scatter(centrality, degree, color = "blue")
	plt.scatter(centrality, degree, c = degree, cmap = "plasma", s = degree, edgecolors = "black")
	plt.colorbar()
	plt.xscale("log", nonpositive='clip')
	plt.yscale("log", nonpositive='clip')
	plt.xlabel("Betweennesss Centrality")
	plt.ylabel("K")
	plt.title("K vs Betweennesss Centrality (log-log)")
	fig.savefig(dataPathOut_ + "BetCentrality_loglog_" + suffix + ".eps", bbox_inches = "tight")

	
	print("\nAll plots saved!\n")

#plt.show()












