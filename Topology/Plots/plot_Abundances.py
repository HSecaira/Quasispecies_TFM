"""
Script that produces a plot of theoretical vs empirical abundances
"""

# Imports
import numpy as np;
import os, sys;
import matplotlib.pyplot as plt;
import helper_plots as h;
import brewer2mpl;
import matplotlib.cm as cm;


######################################################################
# Some functions

def MSEAbundances(theoAbundance, empAbundance):
	"""
	Function that calculates the MSE for empirical vs theoretical abundances 
	Inputs: 
		theoAbundance: array of floats containing the empirical abundances 
		empAbundance: array of floats containing the empirical abundances 
	Outputs: 
		MSE: float containing the MSE  
	"""

	n = theoAbundance.shape[0]

	MSE = (1/n) * (np.sum(np.power((theoAbundance - empAbundance), 2)))

	return MSE

def correlationAbundances(theoAbundance, empAbundance):
	"""
	Function that calculates the pearson correlation coefficient between empirical and theoretical abundances 
	Inputs: 
		theoAbundance: array of floats containing the empirical abundances 
		empAbundance: array of floats containing the empirical abundances 
	Outputs: 
		corr: float containing the correlation coefficient  
	"""

	xMean = np.mean(theoAbundance)
	yMean = np.mean(empAbundance)

	numerator = ((theoAbundance - xMean) @ (empAbundance - yMean))

	denominator = np.sqrt(np.sum(np.power((theoAbundance - xMean), 2)) * np.sum(np.power((empAbundance - yMean), 2)))


	r = numerator/denominator

	return r

######################################################################



# Define paths
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/";
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
dataPathOut_ = dataPathIn_ + "PlotsAbundances/"


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

empAbundances = []
theoAbundances = []

# Iterate over files to retrieve the mean degree and characteristic path length 
print("\nLoading data\n")
for temp in range(temperature.shape[0]):
	for trans in range(transfer.shape[0]):
		file = fileNames[temp][trans]
		print("\tRetrieveing data from file", file)

		fIn = os.path.join(dataPathIn_ + file + "_reg1_all_all_trim_merged_filter_sort_filter_length_collapsed/Abundances/", "AbundancesEmpvsTheo" + ".csv");
		# Retrieve values as dictionaries
		empAbundance, theoAbundance = h.RetrieveAbundances(fIn)

		# Add to lists
		empAbundances.append(empAbundance)
		theoAbundances.append(theoAbundance)


# Size of empAbundances and theoAbundances must be 16, as the above matrix
#print(len(empAbundances))
#print(len(theoAbundances))
print("\nData loaded!\n")




# Do plots
print("\nPreparing plots of abundances comparison\n")
MSEMatrix = np.zeros((fileNames.shape[0], fileNames.shape[1]), dtype = float)
rMatrix = np.zeros((fileNames.shape[0], fileNames.shape[1]), dtype = float)
count = 0

for temp in range(temperature.shape[0]):
	for trans in range(transfer.shape[0]):

		# Get empirical and theoretical abundances
		Eiv = np.array(list(theoAbundances[count].values()))
		Ab = np.array(list(empAbundances[count].values()))

		# Update count
		count += 1
		
		# Calculate MSE
		MSE = MSEAbundances(Eiv, Ab)
		
		# Calculate correlation coefficient
		r = correlationAbundances(Eiv, Ab)

		# Save MSE and r into a matrix
		MSEMatrix[temp][trans] = MSE 
		rMatrix[temp][trans] = r
		
		# Plot empirical abundance vs eigenvector centrality
		fig = plt.figure()
		plt.plot(Eiv, Eiv, color="#343A40")
		plt.scatter(Eiv, Eiv, color="#343A40", alpha = 0.7, edgecolor="black", linewidth=0.6)
		plt.scatter(Eiv, Ab, color="#CED4DA", alpha = 0.7, edgecolor="#495057", linewidth=0.6)
		plt.yscale("log")
		plt.xscale("log")
		plt.xlabel("$v_{1}(i)$")
		plt.ylabel("Empirical Abundance")
		# Save Figure
		fig.savefig(dataPathOut_ + "Abundances_" + str(temperature[temp]) + "C" + "_" + str(transfer[trans]) + ".svg", bbox_inches = "tight")
		fig.savefig(dataPathOut_ + "Abundances_" + str(temperature[temp]) + "C" + "_" + str(transfer[trans]) + ".pdf", bbox_inches = "tight")
		
		
print("\nAbundance plots done!\n")

# HeatMap for MSE
fig, ax = plt.subplots()
im = ax.imshow(MSEMatrix, cmap = "Blues")
# Show ticks
ax.set_xticks(np.arange(len(transfer)))
ax.set_yticks(np.arange(len(temperature)))
# Label ticks
ax.set_xticklabels(transfer)
ax.set_yticklabels([str(temp) + " C" for temp in temperature])
# Loop over data to create text annotations
for i in range(len(temperature)):
	for j in range(len(transfer)):
		text = ax.text(j, i, str(round(MSEMatrix[i][j], 6)), 
						ha = "center", va = "center", color = "black")
# Set titles
ax.figure.colorbar(im, ax=ax, label = "MSE")
ax.set_xlabel("Passages")
ax.set_ylabel("Temperature")
fig.savefig(dataPathOut_ + "MSEHeatMap" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "MSEHeatMap" + ".pdf", bbox_inches = "tight")

# HeatMap for R
fig, ax = plt.subplots()
im = ax.imshow(rMatrix, cmap = "Blues")
# Show ticks
ax.set_xticks(np.arange(len(transfer)))
ax.set_yticks(np.arange(len(temperature)))
# Label ticks
ax.set_xticklabels(transfer)
ax.set_yticklabels([str(temp) + " C" for temp in temperature])
# Loop over data to create text annotations
for i in range(len(temperature)):
	for j in range(len(transfer)):
		text = ax.text(j, i, str(round(rMatrix[i][j], 6)), 
						ha = "center", va = "center", color = "black")
# Set titles
ax.figure.colorbar(im, ax=ax, label = "$r$")
ax.set_xlabel("Passages")
ax.set_ylabel("Temperature")
fig.savefig(dataPathOut_ + "rHeatMap" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "rHeatMap" + ".pdf", bbox_inches = "tight")


# 2D plots
# Set colors for each experimental conditions
colors = ["#ADE8F4", "#00B4D8", "#0077B6", "#023E8A"]
colorsTemp = ["#212529", "#495057", "#6C757D", "#ADB5BD"]
# Set markers for each experimental condition
markers = ["^", "s", "p", "H"]
labels = [4, 10, 30, 60]


fig = plt.figure()
# Iterate over conditions
for temp in range(temperature.shape[0]):
	plt.plot(transfer, MSEMatrix[temp], color = colorsTemp[temp], 
		linestyle = "solid", label = str(temperature[temp]) + " C")
	plt.scatter(transfer, MSEMatrix[temp], facecolor=colorsTemp[temp], marker="o", alpha=0.6, s=65, 
		edgecolor='black', label = str(temperature[temp]) + " C", linestyle = "solid")
plt.legend(loc = "best", ncol = 1, title="Temperature", fancybox=True)
plt.xlabel("Passages") 
plt.ylabel("Mean Squared Error")
fig.savefig(dataPathOut_ + "MSETemp" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "MSETemp" + ".pdf", bbox_inches = "tight")


fig = plt.figure()
# Iterate over conditions
for temp in range(temperature.shape[0]):
	plt.plot(transfer, rMatrix[temp], color = colorsTemp[temp], marker="o", markersize=6.5, alpha = 0.6, 
		linewidth=2, label = str(temperature[temp]) + " C", linestyle = "solid")
plt.legend(loc = "best", ncol = 1, title="Temperature", fancybox=True)
plt.xlabel("Passage") 
plt.ylabel("Pearson Correlation Coefficient")
fig.savefig(dataPathOut_ + "rTemp" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "rTemp" + ".pdf", bbox_inches = "tight")


MSEMatrixTemp = MSEMatrix.transpose()
rMatrixTemp = rMatrix.transpose()

fig = plt.figure()
# Iterate over conditions
for trans in range(transfer.shape[0]):
	plt.plot(temperature, MSEMatrixTemp[trans], color = colors[trans], marker="o", markersize=6.5, alpha = 0.6, 
		linewidth=2, label = str(transfer[trans]), linestyle = "solid")
plt.legend(loc = "best", ncol = 1, title="Passage", fancybox=True)
plt.xlabel("Temperature") 
plt.ylabel("Mean Squared Error")
fig.savefig(dataPathOut_ + "MSETrans" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "MSETrans" + ".pdf", bbox_inches = "tight")

fig = plt.figure()
# Iterate over conditions
for trans in range(transfer.shape[0]):
	plt.plot(temperature, rMatrixTemp[trans], color = colors[trans], marker="o", markersize=6.5, alpha = 0.6, 
		linewidth=2, label = str(transfer[trans]), linestyle = "solid")
plt.legend(loc = "best", ncol = 1, title="Passage", fancybox=True)
plt.xlabel("Temperature") 
plt.ylabel("Pearson Correlation Coefficient")
fig.savefig(dataPathOut_ + "rTrans" + ".svg", bbox_inches = "tight")
fig.savefig(dataPathOut_ + "rTrans" + ".pdf", bbox_inches = "tight")

#plt.show()

print("\nAll plots done!\n")


