""" 
Script to generate a plot of frequency of mutations 
per codon position of the three sequenced regions
"""

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper as h;
import networkx as nx;
from mpl_toolkits import mplot3d;
from mpl_toolkits.mplot3d import Axes3D;

##################################################################
# Some function
def histogramArray(degrees, minDegreeInterval = 0, maxDegreeInterval = 10, base = 2.0):
	"""
	Function that performs a logarithmic binning for the degrees in a network 
	Inputs:
		degrees: array of ints containing the degrees of a network 
		minDegreeInterval: int indicating the minimum degree to create the binning intervals
			Beware! this depends on the base because minDegree = base^minDegreeInterval
		maxDegreeInterval: int indicating the maximum degree to create the binning intervals
			Beware! this depends on the base because maxDegree = base^maxDegreeInterval 
		base: int/float indicating the base for the logarithmic binning
	Outputs:
		histogram_norm: array of floats containing all the data points that falls on a given bin. The array has been
		normalized by dividing each bin by its width
		binsUnique: array of ints (since intervals are discrete) containing the intervals of all bins 
		width: array of ints/floats containing the width of each interval 
	"""

	# Create array of bins
	bins = np.logspace(minDegreeInterval, maxDegreeInterval, base = base, dtype = int)
	# Get the unique bins
	binsUnique = np.unique(bins)
	# Get the widths of each interval
	widths = (binsUnique[1:] - binsUnique[:-1])
	# Get the histogram
	histogram = np.histogram(degrees, bins = binsUnique)
	# Normalize bins by width 
	histogram_temp = histogram[0]/widths

	return histogram_temp, binsUnique, widths

def calculateParametersLSR(x, y):
	"""
	Function that calculates w* and b* for a least square regression and the standard error for b*
	Inputs:
		x: array containing the x-values in log-scale
		y: array containing the y-values in log-scale
	Outputs:
		b: float containing the intercept term of the regression line
		w: float containing the slope of the regression line
		SEw: float containing the standard error of b*
	"""

	meanTarget = np.mean(y)
	meanFeature = np.mean(x)

	centeredTarget = y - meanTarget
	centeredFeature = x - meanFeature

	w = (centeredFeature @ centeredTarget)/(centeredFeature @ centeredFeature)

	b = meanTarget - w *  meanFeature

	# Standard error
	yHat = b + w * x
	n = x.shape[0]

	SEw = np.sqrt((1/(n - 2)) * ((np.sum(np.power((y - yHat), 2)))/(np.sum(np.power((x - meanFeature), 2)))))



	return b, w, SEw

##################################################################



# Variables and paths
print("\nLoading variables and paths\n")

picsPath = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/plotsNumSeqs/" 

reg1 = "Reg1/";
dataPathInR1 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg1 + "/results/"

reg2 = "Reg2/";
dataPathInR2 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg2 + "/results/"

reg3 = "Reg3/";
dataPathInR3 = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/ORF/" + reg3 + "/results/"

print("\nAll paths and variables loaded!")

# Declare structures
numProts = []
numSeqs = []
alphas = []
SEalphas = []
numSeqsTotal = []
numProtsTotal = []
prots_fits = []
all_seqs = []
all_prots = []

# Load data
fInR1 = os.path.join(dataPathInR1, "NeutralNetworksMetaInfo_A2_Reg1" + ".csv")
seqsPerFileR1 = h.loadSeqsPerFile(fInR1)
# Get the counts of sequencies across files
numSeqs_, numProts_ = np.unique(np.array(list(seqsPerFileR1.values())), return_counts= True)
numSeqs.append(numSeqs_)
numProts.append(numProts_)
numProtsTotal.append(np.sum(numProts_))
numSeqsTotal.append(numSeqs_ @ numProts_)
print("Number of prots", np.sum(numProts_))
print("Number of seqs", numSeqs_ @ numProts_)
# Logarithmic binning
histProts_norm, bins, widths = histogramArray(np.array(list(seqsPerFileR1.values())), 
								maxDegreeInterval = 14)
seqs = bins[:-1]
all_seqs.append(seqs)
prots = histProts_norm
all_prots.append(prots)
# Least squares regression
seqs_scaled = np.log10(seqs[prots != 0])
prots_scaled = np.log10(prots[prots != 0])
b, w, SEw = calculateParametersLSR(seqs_scaled, prots_scaled)
# Fit
prots_fits.append((10 ** b) * (seqs ** w))
print("LSR estimation: ", b, w, SEw)
alphas.append(w)
SEalphas.append(SEw)


fInR2_1 = os.path.join(dataPathInR2, "NeutralNetworksMetaInfo_A1_Reg2" + ".csv")
seqsPerFileR2_1 = h.loadSeqsPerFile(fInR2_1)
# Get the counts of sequencies across files
numSeqs_, numProts_ = np.unique(np.array(list(seqsPerFileR2_1.values())), return_counts= True)
numSeqs.append(numSeqs_)
numProts.append(numProts_)
numProtsTotal.append(np.sum(numProts_))
numSeqsTotal.append(numSeqs_ @ numProts_)
print("Number of prots", np.sum(numProts_))
print("Number of seqs", numSeqs_ @ numProts_)
# Logarithmic binning
histProts_norm, bins, widths = histogramArray(np.array(list(seqsPerFileR2_1.values())), 
								maxDegreeInterval = 17)
seqs = bins[:-1]
all_seqs.append(seqs)
prots = histProts_norm
all_prots.append(prots)
# Least squares regression
seqs_scaled = np.log10(seqs[prots != 0])
prots_scaled = np.log10(prots[prots != 0])
b, w, SEw = calculateParametersLSR(seqs_scaled, prots_scaled)
# Fit
prots_fits.append((10 ** b) * (seqs ** w))
print("LSR estimation: ", b, w, SEw)
alphas.append(w)
SEalphas.append(SEw)



fInR2_2 = os.path.join(dataPathInR2, "NeutralNetworksMetaInfo_Replicase_Reg2" + ".csv")
seqsPerFileR2_2 = h.loadSeqsPerFile(fInR2_2)
# Get the counts of sequencies across files
numSeqs_, numProts_ = np.unique(np.array(list(seqsPerFileR2_2.values())), return_counts= True)
numSeqs.append(numSeqs_)
numProts.append(numProts_)
print("Number of prots", np.sum(numProts_))
print("Number of seqs", numSeqs_ @ numProts_)
numProtsTotal.append(np.sum(numProts_))
numSeqsTotal.append(numSeqs_ @ numProts_)
# Logarithmic binning
histProts_norm, bins, widths = histogramArray(np.array(list(seqsPerFileR2_2.values())), 
								maxDegreeInterval = 18)
seqs = bins[:-1]
all_seqs.append(seqs)
prots = histProts_norm
all_prots.append(prots)
# Least squares regression
seqs_scaled = np.log10(seqs[prots != 0])
prots_scaled = np.log10(prots[prots != 0])
b, w, SEw = calculateParametersLSR(seqs_scaled, prots_scaled)
# Fit
prots_fits.append((10 ** b) * (seqs ** w))
print("LSR estimation: ", b, w, SEw)
alphas.append(w)
SEalphas.append(SEw)


fInR3 = os.path.join(dataPathInR3, "NeutralNetworksMetaInfo_Replicase_Reg3" + ".csv")
seqsPerFileR3 = h.loadSeqsPerFile(fInR3)
# Get the counts of sequencies across files
numSeqs_, numProts_ = np.unique(np.array(list(seqsPerFileR3.values())), return_counts= True)
numSeqs.append(numSeqs_)
numProts.append(numProts_)
numProtsTotal.append(np.sum(numProts_))
numSeqsTotal.append(numSeqs_ @ numProts_)
print("Number of prots", np.sum(numProts_))
print("Number of seqs", numSeqs_ @ numProts_)
# Logarithmic binning
histProts_norm, bins, widths = histogramArray(np.array(list(seqsPerFileR3.values())), 
								maxDegreeInterval = 14)
seqs = bins[:-1]
all_seqs.append(seqs)
prots = histProts_norm
all_prots.append(prots)
# Least squares regression
seqs_scaled = np.log10(seqs[prots != 0])
prots_scaled = np.log10(prots[prots != 0])
b, w, SEw = calculateParametersLSR(seqs_scaled, prots_scaled)
# Fit
prots_fits.append((10 ** b) * (seqs ** w))
print("LSR estimation: ", b, w, SEw)
alphas.append(w)
SEalphas.append(SEw)


#########################################################################################
# Preparing plots
labels = ["A2 Region 1", "A1 Region 2", "Replicase Region 2", "Replicase Region 3"]
colors = ["#212529", "#495057", "#ADB5BD", "#DEE2E6"]


#fig = plt.figure()
#for i, label in enumerate(labels):
	#plt.scatter(numSeqs[i], numProts[i], color=colors[i],
				#alpha=0.7, marker="o", label=label, s=60, edgecolor="black",
				#linewidth = 0.5)
#### Plot log binning
###plt.scatter(all_seqs[0], all_prots[0], alpha=0.7, marker="o", color="red", s=60)
#### Plot fit
###plt.plot(all_seqs[0], prots_fits[0], color="red", linestyle="dashed")
#plt.legend(title="Proteins")
#plt.yscale("log")
#plt.xscale("log")
#plt.xlabel("Number of Sequences")
#plt.ylabel("Number of Proteins Codified by a Sequence")
#fig.savefig(picsPath + "NumSeqsProts" + ".pdf", bbox_inches = "tight"); 
#fig.savefig(picsPath + "NumSeqsProts" + ".svg", bbox_inches = "tight");  
#
#
#print(len(prots_fits))
#
## Log binning plot
#fig = plt.figure()
#for i, label in enumerate(labels):
	#plt.scatter(all_seqs[i], all_prots[i], color=colors[i],
				#alpha=0.7, marker="o", s=60, edgecolor="black",
				#linewidth = 0.5)
	#plt.plot(all_seqs[i], prots_fits[i], linestyle="dashed", color=colors[i])
#plt.yscale("log")
#plt.xscale("log")
#plt.ylim(1e-5, 1e5)
#plt.xlabel("Number of Sequences")
#plt.ylabel("Number of Proteins Codified by a Sequence")
#fig.savefig(picsPath + "NumSeqsProtsBinning" + ".pdf", bbox_inches = "tight"); 
#fig.savefig(picsPath + "NumSeqsProtsBinning" + ".svg", bbox_inches = "tight"); 

##fig = plt.figure()
##ax = Axes3D(fig)
##for i, label in enumerate(labels):
	##ax.scatter(numSeqsTotal[i], numProtsTotal[i], -alphas[i], color=colors[i],
		##label=label, s=65, edgecolor="black")
	### Plot error bars
	##ax.plot([numSeqsTotal[i], numSeqsTotal[i]], [numProtsTotal[i], numProtsTotal[i]], 
			##[-alphas[i]+SEalphas[i], -alphas[i]-SEalphas[i]], marker="_", color="black")
##
###ax.legend()
##ax.set_ylabel("Number of Proteins") 
##ax.set_zlabel("$\\alpha$")
##ax.set_xlabel("Number of Sequences")
###fig.savefig(picsPath + "NumSeqsProtsAlphas" + ".pdf", bbox_inches = "tight"); 
###fig.savefig(picsPath + "NumSeqsProtsAlphas" + ".svg", bbox_inches = "tight");  

fig = plt.figure()
for i, label in enumerate(labels):
	plt.scatter(numSeqsTotal[i], alphas[i], color=colors[i],
		label=label, s=65, edgecolor="black", linewidth=0.5)
	# Plot error bars
	plt.errorbar(numSeqsTotal[i], alphas[i], yerr=SEalphas[i],
			marker="_", color="grey")

plt.ylim(-2.2, -1.8)
#ax.legend()
plt.ylabel("$\\alpha$") 
plt.xlabel("Number of Sequences")
fig.savefig(picsPath + "NumSeqsProtsAlphas" + ".pdf", bbox_inches = "tight"); 
fig.savefig(picsPath + "NumSeqsProtsAlphas" + ".svg", bbox_inches = "tight");  


plt.show()





