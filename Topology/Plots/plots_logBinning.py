"""
Script that performs logarithmic binning
"""

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper_plots as h;
import helper_logBinning as hlog;
import brewer2mpl;
import matplotlib.cm as cm;
from scipy.stats import linregress;

######################################################################
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
#dataPathOut_ = dataPathIn_ + "Plots40C/"

file = "EL4-reg1_S4_L001_001_all_trim_merged_filter_sort_filter_length_collapsed/"
#file = "EL30_reg1_all_all_trim_merged_filter_sort_filter_length_collapsed/"

################################################################3
# Load degrees
suffix = "GrandNetwork"
print("\tDegree distribution\n")
#fIn = os.path.join(dataPathIn_ + file + "ResultsTopology/", "Nodes_and_Degrees_" + suffix + ".csv")
fIn = os.path.join(dataPathInGN + "ResultsTopology/", "Nodes_and_Degrees_" + suffix + ".csv")


# Retrieve the degrees
degrees = h.RetrieveDegrees(fIn)
# Get the unique counts for each degree
degreeVertex, countDegree = np.unique(degrees, return_counts= True)
# Retrieve cummulative distribution
Pk = h.RetrievePk(fIn)


###########################################################################3


# Logarithmic binning
histDegrees_norm, bins, widths = hlog.histogramDegrees(degrees)

k_ = bins[:-1]
pkBin_ = histDegrees_norm	

# Least squares regression
k_scaled = np.log10(k_[pkBin_ != 0])
pkBin_scaled = np.log10(pkBin_[pkBin_ != 0])

b, w, SEw = hlog.calculateParametersLSR(k_scaled[1:], pkBin_scaled[1:])
print("LSR estimation: ", b, w, SEw)
mSP, bSP, r_value, p_value, std_err = linregress(k_scaled[1:], pkBin_scaled[1:])
print("SCIPY estimation: ", bSP, mSP, std_err)

bCum, wCum, SEwCum = hlog.calculateParametersLSR(np.log10(degrees), np.log10(Pk))
#print(bCum, wCum, wCum - 1)

pkBin_fit = 10 ** (b) * k_** w
Pk_fit = degrees ** wCum
PkCum_fit = degrees ** (wCum - 1)

###############################

# Exponent with Newman's formula
kmin = 3
alpha, sigma = hlog.alphaNewmanClauset(degrees, kmin)
pk_alpha = degrees ** (-alpha)

# Plot
fig = plt.figure()
plt.xscale('log', nonpositive='clip')
plt.yscale('log', nonpositive='clip')
plt.scatter(degreeVertex, countDegree/degrees.shape[0], alpha = 0.5, color = "blue", label = "$p(k)$")
plt.scatter(k_, pkBin_, alpha = 0.5, color = "red", label = "Binning")
#plt.scatter(degrees, Pk, alpha = 0.5, color = "green", label = "$P(k)$")
plt.plot(k_, pkBin_fit, linestyle="dashed", color = "maroon", 
		label = "fit Binning kmin " + str(2) + " alpha: " + str(round(w, 2)) + "$\\pm: $" + str(round(SEw, 3)))
#plt.plot(degrees, PkCum_fit, linestyle="dashed", color = "darkseagreen", label = "fit Cumulative" + " alpha: " + str(round(wCum - 1, 2)))
#plt.plot(degrees, pk_alpha, linestyle="dashed", color = "black", label = "fit Newman kmin " + str(kmin) + " alpha: " + str(round(alpha, 2)) + " sigma: " + str(round(sigma, 2)))
# Set y-lim manually. Otherwise points are not plotted
#plt.ylim(1e-7, 1.8)
plt.legend()


## KClustering coefficient
# Load data
#fIn2 = os.path.join(dataPathIn_ + file + "ResultsTopology/", "KClusteringCoefficient_" + suffix + ".csv")
fIn2 = os.path.join(dataPathInGN + "ResultsTopology/", "KClusteringCoefficient_" + suffix + ".csv")

# Retrieve the K clustering coefficient
KClustCoeff = h.RetrieveKClustering(fIn2)
# Append to list
KClust = np.array(list(KClustCoeff.values()))
KDegrees= np.array(list(KClustCoeff.keys()))
PCK = h.RetrievePk(fIn2)
print("CK degrees", KDegrees)


fIn3 = os.path.join(dataPathIn_ + file + "ResultsTopology/", "ClusteringCoefficient_" + suffix + ".csv")
# Retrieve local clustering coefficient
clustCoeffDict = h.RetrieveClustCoeffKnn(fIn3)

# ALternative logarithmic binning
histCCoeff_normAlt = hlog.histogramClustCoeffKnnAlternative(KClustCoeff, bins, widths)


# exponents
b, w, ESw = hlog.calculateParametersLSR(np.log10(KDegrees[(KDegrees >= 4) & (KDegrees <= 55)]),
			np.log10(KClust[(KDegrees >= 4) & (KDegrees <= 55)]))

CK_bin_fit = 10 **(b) * KDegrees ** w
print("clustering", w, ESw)


fig = plt.figure()
plt.xscale('log', nonpositive='clip')
plt.yscale('log', nonpositive='clip')
plt.scatter(KDegrees, KClust, label = "C(k)", alpha=0.5, color = "blue")
#plt.scatter(k_, histCCoeff_norm, label="binning", alpha=0.5, color = "red")
#plt.scatter(k_, histCCoeff_normAlt, label="Binning", alpha=0.5, color = "green")
#plt.plot(k_, histCCoeff_normAlt, label="Alt binning1", color = "green")
#plt.scatter(KDegrees, PCK, label="PCK", alpha = 0.5, color = "green")
plt.plot(KDegrees, CK_bin_fit, linestyle="dashed", color = "#52b788", label = "fit Binning kmin " + str(5))
plt.ylim(1e-3, 2)
plt.legend()



## Knn
# Load data
#fIn4 = os.path.join(dataPathIn_ + file + "ResultsTopology/", "KNearestNeighborDegree_" + suffix + ".csv")
fIn4 = os.path.join(dataPathInGN + "ResultsTopology/", "KNearestNeighborDegree_" + suffix + ".csv")
# Retrieve the K clustering coefficient
KNNeighbors = h.RetrieveNearestNeighborDegree(fIn4)
# Append to list
KNN = np.array(list(KNNeighbors.values()))
KNNDegree = np.array(list(KNNeighbors.keys()))
PCKNN = h.RetrievePk(fIn4)
print("KnnDegree", KNNDegree)


fIn5 = os.path.join(dataPathIn_ + file + "ResultsTopology/", "AverageNeighborDegree_" + suffix + ".csv")
# Retrieve local clustering coefficient
KNNDict = h.RetrieveClustCoeffKnn(fIn5)

# Logarithmic binning
KNN_norm = hlog.histogramClustCoeffKnnAlternative(KNNeighbors, bins, widths)

# exponents

b, w, ESw = hlog.calculateParametersLSR(np.log10(KNNDegree[(KNNDegree >= 6) & (KNNDegree <= 200)]), 
			np.log10(KNN[(KNNDegree >= 6) & (KNNDegree <= 200)]))
KNN_bin_fit = 10 ** (b) * KNNDegree** w

print("Knn", w, ESw)


fig = plt.figure()
plt.xscale('log', nonpositive='clip')
plt.yscale('log', nonpositive='clip')
plt.scatter(KNNDegree, KNN, label = "Knn(k)", alpha=0.5, color = "blue")
#plt.scatter(k_, KNN_norm, label="binning", alpha=0.5, color = "green")
#plt.scatter(KNNDegree, PCKNN, label="PCKnn", alpha = 0.5, color = "green")
plt.plot(KNNDegree, KNN_bin_fit, linestyle="dashed", color = "#52b788", label = "fit")
plt.ylim(1, 1e2 + 20)
plt.legend()

plt.show()














