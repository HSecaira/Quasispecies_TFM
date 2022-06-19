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

######################################################################
# Some functions
def saveExponents(fIn, *args):
	"""
	Function that saves the exponents of a power law distribution 
	Inputs:
		fIn: file to save the exponents
		*args: int(s) or list of ints containing the exponents
	"""



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
dataPathOut_ = dataPathIn_ + "Plots40C/"

################################################################

# Load experimental conditions
fInExpConds = os.path.join(dataPathIn_, "filesNames40C" + ".csv")

# Load the names of the experimental conditions
expConds = h.loadExperimentalConditions(fInExpConds)

# Lists to save the data loaded
degreeVertices = []
countDegrees = []
pkBins = []
ks = []
alphas = []
sigmas = []
pkBin_fits = []
PkCum_fits = []
pk_alphas = []
pkws = []
pkESws = []
wCums = []
allDegrees = []
ESwCums = []

KClust = []
KDegrees = []
CK_bins = []
pCKws = []
pCKESws = []
CK_bin_fits = []

KNN = []
KNNDegrees = []
KNN_bins = []
pKNNws = []
pKNNESws = []
KNN_bin_fits = []

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
	allDegrees.append(degrees)
	# Retrieve cummulative distribution
	Pk = h.RetrievePk(fIn)
	# Get the unique counts for each degree
	degreeVertex, countDegree = np.unique(degrees, return_counts= True);
	# Append to a list
	degreeVertices.append(degreeVertex)
	countDegrees.append(countDegree/degrees.shape[0])

	## Logarithmic binning
	pkBin_, bins, widths = hlog.histogramDegrees(degrees)
	k_ = bins[:-1]
	# Append to lists
	pkBins.append(pkBin_)
	ks.append(k_)
	# Least squares regression
	k_scaled = np.log10(k_[pkBin_ != 0])
	pkBin_scaled = np.log10(pkBin_[pkBin_ != 0])
	# Perform regression. Starts at 1 because p(k) behaves as power law from k_[1]
	pkb, pkw, pkESw = hlog.calculateParametersLSR(k_scaled[1:], pkBin_scaled[1:])
	bCum, wCum, ESwCum = hlog.calculateParametersLSR(np.log10(degrees), np.log10(Pk))
	pkws.append(pkw)
	pkESws.append(pkESw)
	wCums.append(wCum - 1)
	ESwCums.append(ESwCum)
	# Exponent with Newman's formula
	kmin = 3
	alpha, sigma = hlog.alphaNewmanClauset(degrees, kmin)
	alphas.append(alpha)
	sigmas.append(sigma)
	# Get the regression lines and add to lists
	pkBin_fit = 10 ** (pkb) * k_** pkw
	pkBin_fits.append(pkBin_fit)
	PkCum_fit = 10 ** (bCum) * degrees ** (wCum - 1)
	PkCum_fits.append(PkCum_fit)
	pk_alpha = ((alpha - 1) * kmin ** (alpha - 1)) * degrees ** (-alpha)
	pk_alphas.append(pk_alpha)


	print("\tKClustering distribution\n")
	fIn2 = os.path.join(dataPathIn__ + "ResultsTopology/", "KClusteringCoefficient_" + suffix + ".csv")
	# Retrieve the K clustering coefficient
	KClustCoeff = h.RetrieveKClustering(fIn2)
	# Append to list
	KClust.append(np.array(list(KClustCoeff.values())))
	KDegrees.append(np.array(list(KClustCoeff.keys())))	
	# Data for logarithmic binning
	fIn3 = os.path.join(dataPathIn__ + "ResultsTopology/", "ClusteringCoefficient_" + suffix + ".csv")
	# Retrieve local clustering coefficient
	clustCoeffDict = h.RetrieveClustCoeffKnn(fIn3)

	## Logarithmic binning
	CK_bin = hlog.histogramClustCoeffKnnAlternative(KClustCoeff, bins, widths)
	CK_bins.append(CK_bin)
	# Least squares regression
	k_scaled = np.log10(k_[CK_bin != 0])

	#print("clustering", k_[CK_bin != 0])

	CK_bin_scaled = np.log10(CK_bin[CK_bin != 0])
	# Perform regression. Starts at ?? because p(k) behaves as power law from k_[1]
	pCKb, pCKw, pCKESw = hlog.calculateParametersLSR(k_scaled[1:], CK_bin_scaled[1:])
	pCKws.append(pCKw)
	pCKESws.append(pCKESw)
	# Get regression lines and add to list
	CK_bin_fit = 10 ** (pCKb) * k_** pCKw
	CK_bin_fits.append(CK_bin_fit)


	print("\tKNearesNeighbors distribution\n")
	fIn4 = os.path.join(dataPathIn__ + "ResultsTopology/", "KNearestNeighborDegree_" + suffix + ".csv");
	# Retrieve the Knearest neighbor degrees
	KNNeighbors = h.RetrieveNearestNeighborDegree(fIn4)
	KNN.append(np.array(list(KNNeighbors.values())))
	KNNDegrees.append(np.array(list(KNNeighbors.keys())))	
	# Data for logarithmic binning
	fIn5 = os.path.join(dataPathIn__ + "ResultsTopology/", "AverageNeighborDegree_" + suffix + ".csv")
	# Retrieve local clustering coefficient
	KNNDict = h.RetrieveClustCoeffKnn(fIn5)

	## Logarithmic binning
	KNN_bin = hlog.histogramClustCoeffKnnAlternative(KNNeighbors, bins, widths)
	KNN_bins.append(KNN_bin)
	# Least squares regression
	k_scaled = np.log10(k_[KNN_bin != 0])
	KNN_bin_scaled = np.log10(KNN_bin[KNN_bin != 0])
	#print("KNN", k_[KNN_bin != 0])
	# Perform regression. Starts at ?? because p(k) behaves as power law from k_[1]
	pKNNb, pKNNw, pKNNESw = hlog.calculateParametersLSR(k_scaled[1:], KNN_bin_scaled[1:])
	pKNNws.append(pKNNw)
	pKNNESws.append(pKNNESw)
	# Get regression lines and add to list
	KNN_bin_fit = 10 ** (pKNNb) * k_** pKNNw
	KNN_bin_fits.append(KNN_bin_fit)


################################################################
# Load data from the GrandNetwork
print("\nLoading data from the GrandNetwork", "\n")
suffix = "GrandNetwork"
print("\tDegree distribution\n")
fInGN = os.path.join(dataPathInGN + "ResultsTopology/", "Nodes_and_Degrees_" + suffix + ".csv");
# Retrieve the degrees
degreesGN = h.RetrieveDegrees(fInGN);
# Retrieve cummulative distribution
PkGN = h.RetrievePk(fInGN)
# Get the unique counts for each degree
degreeVertexGN, countDegreeGN = np.unique(degreesGN, return_counts= True);
countDegreeGN = countDegree/degrees.shape[0]

## Logarithmic binning
pkBin_GN, binsGN, widthsGN = hlog.histogramDegrees(degreesGN)
k_GN= binsGN[:-1]
# Least squares regression
k_scaledGN = np.log10(k_GN[pkBin_GN != 0])
pkBin_scaledGN = np.log10(pkBin_GN[pkBin_GN != 0])
# Perform regression. Starts at 1 because p(k) behaves as power law from k_[1]
pkbGN, pkwGN, pkESwGN = hlog.calculateParametersLSR(k_scaledGN[1:], pkBin_scaledGN[1:])
bCumGN, wCumGN, ESwCumGN = hlog.calculateParametersLSR(np.log10(degreesGN), np.log10(PkGN))
wCumsGN = wCumGN - 1
# Exponent with Newman's formula
kminGN = 3
alphaGN, sigmaGN = hlog.alphaNewmanClauset(degreesGN, kminGN)
# Get the regression lines and add to lists
pkBin_fitGN = 10 ** (pkbGN) * k_GN** pkwGN
PkCum_fitGN = 10 ** (bCumGN) * degreesGN ** (wCumsGN)
pk_alphaGN = ((alphaGN - 1) * kminGN ** (alphaGN - 1)) * degreesGN ** (-alphaGN)



print("\tKClustering distribution\n")
fIn2GN = os.path.join(dataPathInGN + "ResultsTopology/", "KClusteringCoefficient_" + suffix + ".csv")
# Retrieve the K clustering coefficient
KClustCoeffGN = h.RetrieveKClustering(fIn2GN)
# Append to list
KClustGN = np.array(list(KClustCoeffGN.values()))
KDegreesGN = np.array(list(KClustCoeffGN.keys()))
# Data for logarithmic binning
fIn3GN = os.path.join(dataPathInGN + "ResultsTopology/", "ClusteringCoefficient_" + suffix + ".csv")
# Retrieve local clustering coefficient
clustCoeffDictGN = h.RetrieveClustCoeffKnn(fIn3GN)
## Logarithmic binning
CK_binGN = hlog.histogramClustCoeffKnnAlternative(KClustCoeffGN, binsGN, widthsGN)
# Least squares regression
k_scaledGN = np.log10(k_GN[CK_binGN != 0])
CK_bin_scaledGN = np.log10(CK_binGN[CK_binGN != 0])
# Perform regression. Starts at ?? because p(k) behaves as power law from k_[1]
pCKbGN, pCKwGN, pCKESwGN = hlog.calculateParametersLSR(k_scaledGN[1:], CK_bin_scaledGN[1:])
# Get regression lines and add to list
CK_bin_fitGN = 10 ** (pCKbGN) * k_GN** pCKwGN


print("\tKNearesNeighbors distribution\n")
fIn4GN = os.path.join(dataPathInGN + "ResultsTopology/", "KNearestNeighborDegree_" + suffix + ".csv");
# Retrieve the Knearest neighbor degrees
KNNeighborsGN = h.RetrieveNearestNeighborDegree(fIn4GN)
KNNGN = np.array(list(KNNeighbors.values()))
KNNDegreesGN = np.array(list(KNNeighbors.keys()))
# Data for logarithmic binning
fIn5GN = os.path.join(dataPathInGN + "ResultsTopology/", "AverageNeighborDegree_" + suffix + ".csv")
# Retrieve local clustering coefficient
KNNDictGN = h.RetrieveClustCoeffKnn(fIn5GN)
## Logarithmic binning
KNN_binGN = hlog.histogramClustCoeffKnnAlternative(KNNeighborsGN, binsGN, widthsGN)
# Least squares regression
k_scaledGN = np.log10(k_[KNN_binGN != 0])
KNN_bin_scaledGN = np.log10(KNN_binGN[KNN_binGN != 0])
# Perform regression. Starts at ?? because p(k) behaves as power law from k_[1]
pKNNbGN, pKNNwGN, pKNNESwGN = hlog.calculateParametersLSR(k_scaledGN[1:], KNN_bin_scaledGN[1:])
# Get regression lines and add to list
KNN_bin_fitGN = 10 ** (pKNNbGN) * k_GN** pKNNwGN


# Load data from the ancestor
name = "EL4-reg1_S4_L001_001_all_trim_merged_filter_sort_filter_length_collapsed/"
print("\nLoading data from the Ancestor", "\n")
print("\tDegree distribution\n")
fInAncestor = os.path.join(dataPathInAncestor + name + "ResultsTopology/", "Nodes_and_Degrees_" + suffix + ".csv");
# Retrieve the degrees
degreesAncestor = h.RetrieveDegrees(fInAncestor);
# Retrieve cummulative distribution
PkAncestor = h.RetrievePk(fInAncestor)
# Get the unique counts for each degree
degreeVertexAncestor, countDegreeAncestor = np.unique(degreesAncestor, return_counts= True);
countDegreeAncestor = countDegree/degrees.shape[0]

## Logarithmic binning
pkBin_Ancestor, binsAncestor, widthsAncestor = hlog.histogramDegrees(degreesAncestor)
k_Ancestor= binsAncestor[:-1]
# Least squares regression
k_scaledAncestor = np.log10(k_Ancestor[pkBin_Ancestor != 0])
pkBin_scaledAncestor = np.log10(pkBin_Ancestor[pkBin_Ancestor != 0])
# Perform regression. Starts at 1 because p(k) behaves as power law from k_[1]
pkbAncestor, pkwAncestor, pkESwAncestor = hlog.calculateParametersLSR(k_scaledAncestor[1:], pkBin_scaledAncestor[1:])
bCumAncestor, wCumAncestor, ESwCumAncestor = hlog.calculateParametersLSR(np.log10(degreesAncestor), np.log10(PkAncestor))
wCumsAncestor = wCumAncestor - 1
# Exponent with Newman's formula
kminAncestor = 3
alphaAncestor, sigmaAncestor = hlog.alphaNewmanClauset(degreesAncestor, kminAncestor)
# Get the regression lines and add to lists
pkBin_fitAncestor = 10 ** (pkbAncestor) * k_Ancestor** pkwAncestor
PkCum_fitAncestor = 10 ** (bCumAncestor) * degreesAncestor ** (wCumsAncestor - 1)
pk_alphaAncestor = ((alphaAncestor - 1) * kminAncestor ** (alphaAncestor - 1)) * degreesAncestor ** (-alphaAncestor)



print("\tKClustering distribution\n")
fIn2Ancestor = os.path.join(dataPathInAncestor + name + "ResultsTopology/", "KClusteringCoefficient_" + suffix + ".csv")
# Retrieve the K clustering coefficient
KClustCoeffAncestor = h.RetrieveKClustering(fIn2Ancestor)
# Append to list
KClustAncestor = np.array(list(KClustCoeffAncestor.values()))
KDegreesAncestor = np.array(list(KClustCoeffAncestor.keys()))
# Data for logarithmic binning
fIn3Ancestor = os.path.join(dataPathInAncestor + name + "ResultsTopology/", "ClusteringCoefficient_" + suffix + ".csv")
# Retrieve local clustering coefficient
clustCoeffDictAncestor = h.RetrieveClustCoeffKnn(fIn3Ancestor)
## Logarithmic binning
CK_binAncestor = hlog.histogramClustCoeffKnnAlternative(KClustCoeffAncestor, binsAncestor, widthsAncestor)
# Least squares regression
k_scaledAncestor = np.log10(k_Ancestor[CK_binAncestor != 0])
CK_bin_scaledAncestor = np.log10(CK_binAncestor[CK_binAncestor != 0])
# Perform regression. Starts at ?? because p(k) behaves as power law from k_[1]
pCKbAncestor, pCKwAncestor, pCKESwAncestor = hlog.calculateParametersLSR(k_scaledAncestor[1:], CK_bin_scaledAncestor[1:])
# Get regression lines and add to list
CK_bin_fitAncestor = 10 ** (pCKbAncestor) * k_Ancestor** pCKwAncestor


print("\tKNearesNeighbors distribution\n")
fIn4Ancestor = os.path.join(dataPathInAncestor + name + "ResultsTopology/", "KNearestNeighborDegree_" + suffix + ".csv");
# Retrieve the Knearest neighbor degrees
KNNeighborsAncestor = h.RetrieveNearestNeighborDegree(fIn4Ancestor)
KNNAncestor = np.array(list(KNNeighbors.values()))
KNNDegreesAncestor = np.array(list(KNNeighbors.keys()))
# Data for logarithmic binning
fIn5Ancestor = os.path.join(dataPathInAncestor + name + "ResultsTopology/", "AverageNeighborDegree_" + suffix + ".csv")
# Retrieve local clustering coefficient
KNNDictAncestor = h.RetrieveClustCoeffKnn(fIn5Ancestor)
## Logarithmic binning
KNN_binAncestor = hlog.histogramClustCoeffKnnAlternative(KNNeighborsAncestor, binsAncestor, widthsAncestor)
# Least squares regression
k_scaledAncestor = np.log10(k_[KNN_binAncestor != 0])
KNN_bin_scaledAncestor = np.log10(KNN_binAncestor[KNN_binAncestor != 0])
# Perform regression. Starts at ?? because p(k) behaves as power law from k_[1]
pKNNbAncestor, pKNNwAncestor, pKNNESwAncestor = hlog.calculateParametersLSR(k_scaledAncestor[1:], KNN_bin_scaledAncestor[1:])
# Get regression lines and add to list
KNN_bin_fitAncestor = 10 ** (pKNNbAncestor) * k_Ancestor** pKNNwAncestor

################################################################

# Plots
# Set colors for each experimental conditions
#colors = cm.get_cmap("gist_rainbow")(np.linspace(0, 1, len(expConds)))
colors = ["#ADE8F4", "#00B4D8", "#0077B6", "#023E8A"]
# Set markers for each experimental condition
markers = ["^", "s", "p", "H"]
labels = [4, 10, 30, 60]
colors_alpha = ["white", "#ADE8F4", "#00B4D8", "#0077B6", "#023E8A"]
markers_alpha = ["o", "^", "s", "p", "H"]
labels_alpha = [2, 4, 10, 30, 60]


print("\nPlot for degree binning\n")
fig = plt.figure(figsize=(9, 5))
ax = fig.add_subplot(111)

# Plot for experimental conditions
for i in range(len(expConds)):
	label = labels[i]
	ax.scatter(ks[i], pkBins[i], marker = markers[i], 
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth=0.2, s=85)
	print("\tTransfer:", label)
	print("\tw binning:", pkws[i])
	print("\tw cumulative:", wCums[i])
	print("\talpha:", alphas[i], "sigma:", sigmas[i])

# Add plot of the Grand Network
GrandNetwork = ax.scatter(k_GN, pkBin_GN, marker = "o",
		alpha = 0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
ax.plot(k_GN, pkBin_fitGN, linestyle="dashed", color = "#6c757d", label = "GN")
print("\tGrandNetwork")
print("\tw binning: ", pkwGN, " +-:", pkESwGN)
print("\tw cumultive: ", wCumsGN)
print("\talpha: ", alphaGN, "sigma: ", sigmaGN)


# Add plot of the ancestor
Ancestor = ax.scatter(k_Ancestor, pkBin_Ancestor, marker = "o", 
		alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
ax.plot(k_Ancestor, pkBin_fitAncestor, linestyle="dotted", color = "#495057", label = "Ancestor")
print("\tAncestor")
print("\tw binning: ", pkwAncestor, " +-:", pkESwAncestor)
print("\tw cumultive: ", wCumsAncestor)
print("\talpha: ", alphaAncestor, "sigma: ", sigmaAncestor)
# Add alpha of the ancestor to the list of alphas
pkws = [pkwAncestor, *pkws]
pkESws = [pkESwAncestor, *pkESws]

ax.legend(bbox_to_anchor=(0.88,1), loc="best", ncol = 1, title="LSR Fit", fancybox=True)
ax.set_ylim(1e-7, 1e0)
ax.set_xscale("log", nonpositive='clip')
ax.set_yscale("log", nonpositive='clip')
ax.set_xlabel("$k$")
ax.set_ylabel("$p(k)$")
#ax.set_title("Logarithmic Binning (40 C)")

##plt.subplots_adjust(right=0.8)
### Add first legend: only labeled data
##leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
### Add second legend
##leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
### Add the first legend back
##ax.add_artist(leg1)
# Save
fig.savefig(dataPathOut_ + "Degree_distributionBinning" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "Degree_distributionBinning" + ".pdf", bbox_inches = "tight")
# Plot for estimated coefficients by logarithmic binnin
fig = plt.figure(figsize=(9, 5))
for i in range(len(labels_alpha)):
	plt.scatter(labels_alpha[i], pkws[i], marker = markers_alpha[i], color = colors_alpha[i], alpha = 1,  s=85,
				edgecolor="black")
plt.plot(labels_alpha, np.array(pkws), color = "#073b4c")
plt.fill_between(labels_alpha, np.array(pkws) - np.array(pkESws), np.array(pkws) + np.array(pkESws),
				alpha = 0.2, color = "#073b4c")
### Grand Network alpha
##pkwGN_ = np.array([pkwGN for i in range(len(labels))])
##pkESwGN_ = np.array([pkESwGN for i in range(len(labels))])
##plt.plot(labels, pkwGN_, color = "dimgray")
##for i in range(len(labels)):
	##plt.scatter(labels[i], pkwGN_[i], marker = "o", color = "dimgray", alpha = 1,  s=85)
##plt.fill_between(labels, pkwGN_ - pkESwGN_, pkwGN_ + pkESwGN_,
				##alpha = 0.1, color = "black")
### Ancestor alpha
##pkwAncestor_ = np.array([pkwAncestor for i in range(len(labels))])
##pkESwAncestor_ = np.array([pkESwAncestor for i in range(len(labels))])
##plt.plot(labels, pkwAncestor_, color="red")
##plt.fill_between(labels, pkwAncestor_ - pkESwAncestor_, pkwAncestor_ + pkESwAncestor_,
				##alpha = 0.1, color = "black")

# Change ylims manually for a nice plot
plt.ylim(-2.6, -2.2)
plt.xlabel("Transfers")
plt.ylabel("$\\hat{\\alpha}_{SN}$")
fig.savefig(dataPathOut_ + "Degree_distributionAlpha" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "Degree_distributionAlpha" + ".pdf", bbox_inches = "tight")

## KClustering coefficient
print("\nPlot for KCLusterin distribution\n")
fig = plt.figure(figsize=(9, 5))
ax = fig.add_subplot(111)

# Plot for experimental conditions
for i in range(len(expConds)):
	label = labels[i]
	ax.scatter(ks[i], CK_bins[i], label = str(label), marker = markers[i], 
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth=0.2, s=85)
	print("\tTransfer:", label)
	print("\tw binning:", pCKws[i], " +-:", pCKESws[i])

# Add plot of the Grand Network
GrandNetwork = ax.scatter(k_GN, CK_binGN, marker = "o",
		alpha = 0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
ax.plot(k_GN, CK_bin_fitGN, linestyle="dashed", color = "dimgray")
print("\tGrandNetwork")
print("\tw binning: ", pCKwGN, "+-:", pCKESwGN)


# Add plot of the ancestor
Ancestor = ax.scatter(k_Ancestor, CK_binAncestor, marker = "o", 
		alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
ax.plot(k_Ancestor, CK_bin_fitAncestor, linestyle="dotted", color = "dimgray")
print("\tAncestor")
print("\tw binning: ", pCKwAncestor, "+-", pCKESwAncestor)

# Set ylim manually
ax.set_ylim(1e-5, 1)
ax.set_xscale("log", nonpositive='clip')
ax.set_yscale("log", nonpositive='clip')
ax.set_xlabel("$k$")
ax.set_ylabel("$C(k)$")
#ax.set_title("Logarithmic Binning (40 C)")
##plt.subplots_adjust(right=0.8)
### Add first legend: only labeled data
##leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
### Add second legend
##leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
### Add the first legend back
##ax.add_artist(leg1)
## Save
fig.savefig(dataPathOut_ + "KCustering_distributionBinning" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "KCustering_distributionBinning" + ".pdf", bbox_inches = "tight")
# Plot for estimated coefficients by logarithmic binnin
fig = plt.figure(figsize=(9, 5))
for i in range(len(labels)):
	plt.scatter(labels[i], pCKws[i], marker = markers[i], color = colors[i], alpha = 1,  s=85)
plt.plot(labels, np.array(pCKws), color = "black")
plt.fill_between(labels, np.array(pCKws) - np.array(pCKESws), np.array(pCKws) + np.array(pCKESws),
				alpha = 0.1, color = "black")
# Change ylims manually for a nice plot
#plt.ylim(-1.75, -1.4)
plt.xlabel("Transfers")
plt.ylabel("$\\hat{\\alpha}_{SN}$")
fig.savefig(dataPathOut_ + "KCustering_distributionAlpha" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "KCustering_distributionAlpha" + ".pdf", bbox_inches = "tight")



## Plot for KNN distribution
print("\nPlot for KNN distribution\n")
fig = plt.figure(figsize=(9, 5))
ax = fig.add_subplot(111)

# Plot for experimental conditions
for i in range(len(expConds)):
	label = labels[i]
	ax.scatter(ks[i], KNN_bins[i], label = str(label), marker = markers[i], 
		alpha=0.7, edgecolor='black', facecolor = colors[i], linewidth=0.2, s=85)
	#plt.plot(k_, KNN_bin_fit, linestyle="dashed", color = "maroon", label = "fit Binning kmin " + str(14))
	print("Transfer", label)
	print("\tw binning:", pKNNws[i], "+-:", pKNNESws[i])

# Add plot of the Grand Network
GrandNetwork = ax.scatter(k_GN, KNN_binGN, marker = "o",
		alpha = 0.5, edgecolor='black', facecolor = "gray", linewidth=0.6, s=50)
ax.plot(k_GN, KNN_bin_fitGN, linestyle="dashed", color = "dimgray")
print("\tGrandNetwork")
print("\tw binning: ", pKNNwGN, "+-:", pKNNESwGN)

# Add plot of the ancestor
Ancestor = ax.scatter(k_Ancestor, KNN_binAncestor, marker = "o", 
		alpha=0.5, edgecolor='black', facecolor = "white", linewidth=0.6, s=50)
ax.plot(k_GN, KNN_bin_fitAncestor, linestyle="dotted", color = "dimgray")
print("\tAncestor")
print("\tw binning: ", pKNNwAncestor, "+-", pKNNESwAncestor)

ax.set_ylim(1e-2, 1e2 + 80)
ax.set_xscale("log", nonpositive='clip')
ax.set_yscale("log", nonpositive='clip')
ax.set_xlabel("$k$")
ax.set_ylabel("$K_{nn}(k)$")
#ax.set_title("Logarithmic Binning (40 C)")
##plt.subplots_adjust(right=0.8)
### Add first legend: only labeled data
##leg1 = ax.legend(bbox_to_anchor=(1,1), loc="upper left", ncol = 1, title="Transfers", fancybox=True)
### Add second legend
##leg2 = ax.legend([GrandNetwork,Ancestor],["Grand Network","Ancestor"], bbox_to_anchor=(1,0.6), loc="center left", ncol = 1, fancybox=True)
### Add the first legend back
##ax.add_artist(leg1)
## Save
fig.savefig(dataPathOut_ + "KNN_distributionBinning" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "KNN_distributionBinning" + ".pdf", bbox_inches = "tight")
# Plot for estimated coefficients by logarithmic binnin
fig = plt.figure(figsize=(9, 5))
for i in range(len(labels)):
	plt.scatter(labels[i], pKNNws[i], marker = markers[i], color = colors[i], alpha = 1,  s=85)
plt.plot(labels, np.array(pKNNws), color = "black")
plt.fill_between(labels, np.array(pKNNws) - np.array(pKNNESws), np.array(pKNNws) + np.array(pKNNESws),
				alpha = 0.1, color = "black")
# Change ylims manually for a nice plot
#plt.ylim(-1.65, -1.35)
plt.xlabel("Transfers")
plt.ylabel("$\\hat{\\alpha}_{SN}$")
fig.savefig(dataPathOut_ + "KNN_distributionAlpha" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "KNN_distributionAlpha" + ".pdf", bbox_inches = "tight")


plt.show()














