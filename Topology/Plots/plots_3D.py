# This script will generate a 3D plot

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper_plots as h;
import matplotlib.cm as cm;
from mpl_toolkits import mplot3d;
from mpl_toolkits.mplot3d import Axes3D;
import seaborn as sns; 

# Define paths
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/";

# Define regions
dicFolder = {};
dicFolder["r1"] = "Reg1_experiments/";
dicFolder["r2"] = "Reg2_experiments/";
dicFolder["r3"] = "Reg3_experiments/";
regionName = "r1";

dataPathIn_ = dataPathIn + dicFolder[regionName] + "SubNetworks/"
#dataPathIn_ = dataPathIn + "SubnetworksConstantT/" + "Reg1/"
dataPathOut_ = dataPathIn_ + "PlotsTemperatures/"
suffix = "GrandNetwork"
print("\nAll data paths loaded!\n")

#############################################################
# Declare some variables and load data
temperature = np.array([40, 37, 33, 30])
transfer = np.array([4, 10, 30, 60])
#fileNames = np.array([["EL21", "EL20", "EL19", "EL18"],
					#["EL25", "EL24", "EL23", "EL22"],
					#["EL29", "EL28", "EL27", "EL26"],
					#["EL33", "EL32", "EL31", "EL30"]])

fileNames = np.array([["EL33", "EL32", "EL31", "EL30"],
					["EL29", "EL28", "EL27", "EL26"],
					["EL25", "EL24", "EL23", "EL22"],
					["EL21", "EL20", "EL19", "EL18"]])

fileNamesTransfers = fileNames.T
print(fileNamesTransfers)

meanDegrees = np.zeros((fileNames.shape[0], fileNames.shape[1]), dtype = float)
meanDegreesSTD = np.zeros((fileNames.shape[0], fileNames.shape[1]), dtype = float)
charPathLengths = np.zeros((fileNames.shape[0], fileNames.shape[1]), dtype = float)
charPathLengthsSTD = np.zeros((fileNames.shape[0], fileNames.shape[1]), dtype = float)
assortativities = np.zeros((fileNames.shape[0], fileNames.shape[1]), dtype = float)
globalclusterings = np.zeros((fileNames.shape[0], fileNames.shape[1]), dtype = float)
densities = np.zeros((fileNames.shape[0], fileNames.shape[1]), dtype = float)
eigenvalues = np.zeros((fileNames.shape[0], fileNames.shape[1]), dtype = float)

# Iterate over files to retrieve the mean degree and characteristic path length 
for temp in range(temperature.shape[0]):
	for trans in range(transfer.shape[0]):
		file = fileNames[temp][trans]
		fIn = os.path.join(dataPathIn_ + file + "_reg1_all_all_trim_merged_filter_sort_filter_length_collapsed/ResultsTopology/", "MetricsInformation_" + suffix + ".csv");
		# Retrieve values
		density, meanDegree, meanDegreeSTD, assortativity, globalClustering, charPathLenght, charPathLenghtSTD, eigenvalue = h.RetrieveMetric(fIn)
		# Add to matrices 
		meanDegrees[temp][trans] = meanDegree
		meanDegreesSTD[temp][trans] = meanDegreeSTD
		charPathLengths[temp][trans] = charPathLenght
		charPathLengthsSTD[temp][trans] = charPathLenghtSTD
		assortativities[temp][trans] = assortativity
		globalclusterings[temp][trans] = globalClustering
		densities[temp][trans] = density
		eigenvalues[temp][trans] = eigenvalue

print("\nAll data loaded\n")

########################################################################################

# Set colors for plots
colors = cm.get_cmap("gist_rainbow")(np.linspace(0, 1, len(temperature)))
##fig = plt.figure()
##ax = plt.axes(projection='3d')
### Iterate over conditions
##for temp in range(temperature.shape[0]):
	##ax.plot3D(temperature, transfer, meanDegrees[temp], color = colors[temp], marker="o", alpha = 0.6, 
		##linewidth=2, label = str(temperature[temp]) + "C", linestyle = "dotted")
##ax.legend(loc = "best", ncol = 1, title="Transfers", fancybox=True)
##ax.set_ylabel("Transfer") 
##ax.set_zlabel("$\\langle k \\rangle$")
##ax.set_xlabel("Temperature (C)")
##ax.set_title("Average Degree vs T vs Transfer")
##fig.savefig(dataPathOut_ + "MeanDegree3D" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "MeanDegree3D" + ".pdf", bbox_inches = "tight")

print("\nPlotting data from mean degree\n")
# 2D figure
fig = plt.figure()
# Iterate over conditions
for temp in range(temperature.shape[0]):
	plt.plot(transfer, meanDegrees[temp], color = colors[temp], marker="o", alpha = 0.6, 
		linewidth=2, label = str(temperature[temp]) + " °C", linestyle = "dotted")
	#plt.fill_between(transfer, meanDegrees[temp] - meanDegreesSTD[temp] , meanDegrees[temp] + meanDegreesSTD[temp],
					#alpha = 0.1, color = colors[temp])
plt.legend(loc = "best", ncol = 1, title="Temperature", fancybox=True)
plt.xlabel("Transfer") 
plt.ylabel("$\\langle k \\rangle$")
fig.savefig(dataPathOut_ + "MeanDegree2D" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "MeanDegree2D" + ".pdf", bbox_inches = "tight")


# HeatMap
fig, ax = plt.subplots()
im = ax.imshow(meanDegrees, cmap = "Blues")
# Show ticks
ax.set_xticks(np.arange(len(transfer)))
ax.set_yticks(np.arange(len(temperature)))
# Label ticks
ax.set_xticklabels(transfer)
ax.set_yticklabels([str(temp) + "$^\\circ$C" for temp in temperature])
# Loop over data to create text annotations
for i in range(len(temperature)):
	for j in range(len(transfer)):
		text = ax.text(j, i, str(round(meanDegrees[i][j], 1)) + " $\\pm$ " + str(round(meanDegreesSTD[i][j], 1)), 
						ha = "center", va = "center", color = "black")
# Set titles
ax.figure.colorbar(im, ax=ax, label = "$\\langle k \\rangle$")
fig.savefig(dataPathOut_ + "MeanDegreeHeatMap" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "MeanDegreeHeatMap" + ".pdf", bbox_inches = "tight")


#### Characteristic path length
##fig = plt.figure()
##ax = plt.axes(projection='3d')
### Iterate over conditions
##for temp in range(temperature.shape[0]):
	##ax.plot3D(temperature, transfer, charPathLengths[temp], color = colors[temp], marker="o", alpha = 0.6, 
		##linewidth=2, label = str(temperature[temp]) + "C", linestyle = "dotted")
##ax.legend(loc = "best", ncol = 1, title="Transfers", fancybox=True)
##ax.set_ylabel("Transfer") 
##ax.set_zlabel("$\\langle d \\rangle$")
##ax.set_xlabel("Temperature (C)")
###ax.set_title("Average Degree vs T vs Transfer")
##fig.savefig(dataPathOut_ + "CharPathLength3D" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "CharPathLength3D" + ".pdf", bbox_inches = "tight")

print("\nPlotting data from characteristic Path Length")
# 2D figure
fig = plt.figure()
# Iterate over conditions
for temp in range(temperature.shape[0]):
	plt.plot(transfer, charPathLengths[temp], color = colors[temp], marker="o", alpha = 0.6, 
		linewidth=2, label = str(temperature[temp]) + " °C", linestyle = "dotted")
	plt.fill_between(transfer, charPathLengths[temp] - charPathLengthsSTD[temp] , charPathLengths[temp] + charPathLengthsSTD[temp],
					alpha = 0.1, color = colors[temp])
plt.legend(loc = "best", ncol = 1, title="Temperature", fancybox=True)
plt.xlabel("Transfer") 
plt.ylabel("$\\langle d \\rangle$")
fig.savefig(dataPathOut_ + "CharPathLength2D" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "CharPathLength2D" + ".pdf", bbox_inches = "tight")

# HeatMap
fig, ax = plt.subplots()
im = ax.imshow(charPathLengths, cmap = "Blues")
# Show ticks
ax.set_xticks(np.arange(len(transfer)))
ax.set_yticks(np.arange(len(temperature)))
# Label ticks
ax.set_xticklabels(transfer)
ax.set_yticklabels([str(temp) + "$^\\circ$C" for temp in temperature])
# Loop over data to create text annotations
for i in range(len(temperature)):
	for j in range(len(transfer)):
		text = ax.text(j, i, str(round(charPathLengths[i][j], 1)) + " $\\pm$ " + str(round(charPathLengthsSTD[i][j], 1)), 
						ha = "center", va = "center", color = "black")
# Set titles
ax.figure.colorbar(im, ax=ax, label = "$\\langle d \\rangle$")
fig.savefig(dataPathOut_ + "CharPathLengthHeatMap" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "CharPathLengthHeatMap" + ".pdf", bbox_inches = "tight")
plt.show()



####fig = plt.figure()
####ax = plt.axes(projection='3d')
##### Iterate over conditions
####for temp in range(temperature.shape[0]):
	####ax.plot3D(temperature, transfer, assortativities[temp], color = colors[temp], marker="o", alpha = 0.6, 
		####linewidth=2, label = str(temperature[temp]) + "C", linestyle = "dotted")
####ax.legend(loc = "best", ncol = 1, title="Temperature", fancybox=True)
####ax.set_ylabel("Transfer") 
####ax.set_zlabel("$r$")
####ax.set_xlabel("Temperature (C)")
#####ax.set_title("Average Degree vs T vs Transfer")
####fig.savefig(dataPathOut_ + "Assortativities3D" + ".svg", bbox_inches = "tight")
####fig.savefig(dataPathOut_ + "Assortativities3D" + ".pdf", bbox_inches = "tight")
##
##
##print("\nPlotting data from assortativity coefficient\n")
### 2D figure
##fig = plt.figure()
### Iterate over conditions
##for temp in range(temperature.shape[0]):
	##plt.plot(transfer, assortativities[temp], color = colors[temp], marker="o", alpha = 0.6, 
		##linewidth=2, label = str(temperature[temp]) + " °C", linestyle = "dotted")
##plt.legend(loc = "best", ncol = 1, title="Temperature", fancybox=True)
##plt.xlabel("Transfer") 
##plt.ylabel("$r$")
##fig.savefig(dataPathOut_ + "Assortativities2D" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "Assortativities2D" + ".pdf", bbox_inches = "tight")
##
### HeatMap
##fig, ax = plt.subplots()
##im = ax.imshow(assortativities, cmap = "Blues")
### Show ticks
##ax.set_xticks(np.arange(len(transfer)))
##ax.set_yticks(np.arange(len(temperature)))
### Label ticks
##ax.set_xticklabels(transfer)
##ax.set_yticklabels([str(temp) + " °C" for temp in temperature])
### Loop over data to create text annotations
##for i in range(len(temperature)):
	##for j in range(len(transfer)):
		##text = ax.text(j, i, round(assortativities[i][j], 2), 
						##ha = "center", va = "center", color = "black")
### Set titles
##ax.figure.colorbar(im, ax=ax, label = "$r$")
##fig.savefig(dataPathOut_ + "AssortativitiesHeatMap" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "AssortativitiesHeatMap" + ".pdf", bbox_inches = "tight")
##
##
##print("\nPlotting data from global clustering coefficient\n")
### 2D figure
##fig = plt.figure()
### Iterate over conditions
##for temp in range(temperature.shape[0]):
	##plt.plot(transfer, globalclusterings[temp], color = colors[temp], marker="o", alpha = 0.6, 
		##linewidth=2, label = str(temperature[temp]) + " °C", linestyle = "dotted")
##plt.legend(loc = "best", ncol = 1, title="Temperature", fancybox=True)
##plt.xlabel("Transfer") 
##plt.ylabel("$C$")
##fig.savefig(dataPathOut_ + "GlobalClusterings2D" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "GlobalClusterings2D" + ".pdf", bbox_inches = "tight")
##
### HeatMap
##fig, ax = plt.subplots()
##im = ax.imshow(globalclusterings, cmap = "Blues")
### Show ticks
##ax.set_xticks(np.arange(len(transfer)))
##ax.set_yticks(np.arange(len(temperature)))
### Label ticks
##ax.set_xticklabels(transfer)
##ax.set_yticklabels([str(temp) + " °C" for temp in temperature])
### Loop over data to create text annotations
##for i in range(len(temperature)):
	##for j in range(len(transfer)):
		##text = ax.text(j, i, round(globalclusterings[i][j], 2), 
						##ha = "center", va = "center", color = "black")
### Set titles
##ax.figure.colorbar(im, ax=ax, label = "$C$")
##fig.savefig(dataPathOut_ + "GlobalClusteringsHeatMap" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "GlobalClusteringsHeatMap" + ".pdf", bbox_inches = "tight")
##
##
##print("\nPlotting data from leading eigenvalue\n")
### 2D figure
##fig = plt.figure()
### Iterate over conditions
##for temp in range(temperature.shape[0]):
	##plt.plot(transfer, eigenvalues[temp], color = colors[temp], marker="o", alpha = 0.6, 
		##linewidth=2, label = str(temperature[temp]) + " °C", linestyle = "dotted")
##plt.legend(loc = "best", ncol = 1, title="Temperature", fancybox=True)
##plt.xlabel("Transfer") 
##plt.ylabel("$\\lambda_{1}$")
##fig.savefig(dataPathOut_ + "Eigenvalues2D" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "Eigenvalues2D" + ".pdf", bbox_inches = "tight")
##
### HeatMap
##fig, ax = plt.subplots()
##im = ax.imshow(eigenvalues, cmap = "Blues")
### Show ticks
##ax.set_xticks(np.arange(len(transfer)))
##ax.set_yticks(np.arange(len(temperature)))
### Label ticks
##ax.set_xticklabels(transfer)
##ax.set_yticklabels([str(temp) + " °C" for temp in temperature])
### Loop over data to create text annotations
##for i in range(len(temperature)):
	##for j in range(len(transfer)):
		##text = ax.text(j, i, round(eigenvalues[i][j], 2), 
						##ha = "center", va = "center", color = "black")
### Set titles
##ax.figure.colorbar(im, ax=ax, label = "$\\lambda_{1}$")
##fig.savefig(dataPathOut_ + "EigenvaluesHeatMap" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "EigenvaluesHeatMap" + ".pdf", bbox_inches = "tight")


##fig = plt.figure()
##ax = plt.axes(projection='3d')
### Iterate over conditions
##for temp in range(temperature.shape[0]):
	##ax.plot3D(meanDegrees[temp], eigenvalues[temp], temperature, color = colors[temp], marker="o", alpha = 0.6, 
		##linewidth=2, label = str(temperature[temp]) + "C", linestyle = "dotted")
##ax.legend(loc = "best", ncol = 1, title="Transfers", fancybox=True)
##ax.set_ylabel("$\\lambda_{1}$") 
##ax.set_zlabel("Temperature (C)")
##ax.set_xlabel("$\\langle k \\rangle$")
###ax.set_title("Average Degree vs T vs Transfer")
##fig.savefig(dataPathOut_ + "EigenvaluesvsMeanDegree3D" + ".svg", bbox_inches = "tight")
##fig.savefig(dataPathOut_ + "EigenvaluesvsMeanDegree3D" + ".pdf", bbox_inches = "tight")
##
##plt.show()



### Mean degree and CharPathLenght plots combined
##fig, (ax1, ax2) = plt.subplots(ncols=2)
##fig.subplots_adjust(wspace=0.01)
##
### Mean degree
##im1 = ax1.imshow(meanDegrees, cmap = "Blues")
### Show ticks
##ax1.set_xticks(np.arange(len(transfer)))
##ax1.set_yticks(np.arange(len(temperature)))
### Label ticks
##ax1.set_xticklabels(transfer)
##ax1.set_yticklabels([str(temp) + " °C" for temp in temperature])
### Loop over data to create text annotations
##for i in range(len(temperature)):
	##for j in range(len(transfer)):
		##text = ax1.text(j, i, round(meanDegrees[i][j], 2), 
						##ha = "center", va = "center", color = "black")
### Set titles
###ax.figure.colorbar(im, ax=ax, label = "$\\langle k \\rangle$")	
##fig.colorbar(im1, ax=ax1, location="left", use_gridspec=False, label="$\\langle k \\rangle$",
			##fraction = 1.5)
##
### CharPathLength
##im2 = ax2.imshow(charPathLengths, cmap = "Blues")
### Show ticks
##ax2.set_xticks(np.arange(len(transfer)))
##ax2.set_yticks(np.arange(len(temperature)))
### Label ticks
##ax2.set_xticklabels(transfer)
##ax2.set_yticklabels([str(temp) + " °C" for temp in temperature])
### Loop over data to create text annotations
##for i in range(len(temperature)):
	##for j in range(len(transfer)):
		##text = ax2.text(j, i, round(charPathLengths[i][j], 2), 
						##ha = "center", va = "center", color = "black")
### Set titles
###ax.figure.colorbar(im, ax=ax, label = "$\\langle d \\rangle$")
##fig.colorbar(im2, ax=ax2, location="right", use_gridspec=False, label="$\\langle d \\rangle$",
			##fraction = 1.5)
##ax2.yaxis.tick_right()
##ax2.tick_params(rotation=0)
###fig.savefig(dataPathOut_ + "MeanDegreeCharPathHeatMap" + ".svg", bbox_inches = "tight")
###fig.savefig(dataPathOut_ + "MeanDegreeCharPathHeatMap" + ".pdf", bbox_inches = "tight")


### Eigenvalue vs mean degree
##fig = plt.figure() 
##for temp in range(temperature.shape[0]):
	##for trans in range(transfer.shape[0]):
		##mdegree = meanDegrees[temp][trans]
		##evalue = eigenvalues[temp][trans]
		##plt.scatter(mdegree, evalue, label = str(temperature[temp]) + "," + str(transfer[trans]), color=colors[temp])
		##plt.text(mdegree, evalue, str(temperature[temp]) + "," + str(transfer[trans]))
###plt.legend()
##plt.show()





#plt.show()




print("\nAll plots saved!")







