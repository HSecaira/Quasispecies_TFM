# Script to generate some plots for the presentation

import numpy as np
import os, sys
import matplotlib.pyplot as plt

# Declare some variables
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/NGS_analyses/plots/"
dataPathOut = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/NGS_analyses/plots/"

###############################################################################
# Some functions

def loadNumReads(fIn):
	"""
	Functionat that loads the numebr of reads per file from a file   
	Inputs:
		fIn: file to read the info  
	Outputs:
		reads: dictionary containing the file name (key) and the number of reads per file (values)
	"""
	# Initialize structures
	reads = {}

	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip header
			if line[0] != "#":
				fileName = str(line.strip().split(":")[0].split("-")[0])
				seqs = int(line.strip().split(":")[1])

				# Add to dictionary
				reads[fileName] = seqs

	return reads

###############################################################################

# Region
reg = "Reg1"

# Load data
fIn1 = os.path.join(dataPathIn, reg + "ReadsBefore" + ".csv")
R1readsBefore = loadNumReads(fIn1)
fIn2 = os.path.join(dataPathIn, reg + "ReadsAfter" + ".csv")
R1readsAfter = loadNumReads(fIn2)
fIn3 = os.path.join(dataPathIn, reg + "ReadsCollapsed" + ".csv")
R1collapsed = loadNumReads(fIn3)

# Convert to array the keys and values of the dictionaries
filesBefore = np.array(list(R1readsBefore.keys()))
numReadsBefore = np.array(list(R1readsBefore.values()))

numReadsAfter = np.array(list(R1readsAfter.values()))
numReadsCollapsed = np.array(list(R1collapsed.values()))

# Do the plots
x = np.arange(numReadsBefore.shape[0])
width = 0.3

fig, ax = plt.subplots(figsize=(12,6))
# Barplots
before = ax.bar(x - 0.30, numReadsBefore, width, label='Before Processing',
				color = "#E9ECEF", edgecolor="black")
after = ax.bar(x, numReadsAfter, width, label='After Processing',
				color = "#ADB5BD", edgecolor="black")
collapsed = ax.bar(x + 0.30, numReadsCollapsed, width, label='Collapsed',
				color = "#343A40", edgecolor="black")

# Custom plot
ax.set_ylabel('Number of Sequences')
ax.set_title('Region 1')
ax.set_ylim(1e3, 1e6)
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(filesBefore, rotation=90)
ax.legend()
fig.tight_layout()
fig.savefig(dataPathOut + "NumReads" + reg + ".pdf", bbox_inches = "tight") 
fig.savefig(dataPathOut + "NumReads" + reg + ".svg", bbox_inches = "tight")



# Region
reg = "Reg2"

# Load data
fIn1 = os.path.join(dataPathIn, reg + "ReadsBefore" + ".csv")
R1readsBefore = loadNumReads(fIn1)
fIn2 = os.path.join(dataPathIn, reg + "ReadsAfter" + ".csv")
R1readsAfter = loadNumReads(fIn2)
fIn3 = os.path.join(dataPathIn, reg + "ReadsCollapsed" + ".csv")
R1collapsed = loadNumReads(fIn3)

# Convert to array the keys and values of the dictionaries
filesBefore = np.array(list(R1readsBefore.keys()))
numReadsBefore = np.array(list(R1readsBefore.values()))

numReadsAfter = np.array(list(R1readsAfter.values()))
numReadsCollapsed = np.array(list(R1collapsed.values()))

# Do the plots
x = np.arange(numReadsBefore.shape[0])
width = 0.3

fig, ax = plt.subplots(figsize=(12,6))
# Barplots
before = ax.bar(x - 0.30, numReadsBefore, width, label='Before Processing',
				color = "#E9ECEF", edgecolor="black")
after = ax.bar(x, numReadsAfter, width, label='After Processing',
				color = "#ADB5BD", edgecolor="black")
collapsed = ax.bar(x + 0.30, numReadsCollapsed, width, label='Collapsed',
				color = "#343A40", edgecolor="black")

# Custom plot
ax.set_ylabel('Number of Sequences')
ax.set_title('Region 2')
ax.set_ylim(1e4, 1e6)
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(filesBefore, rotation=90)
#ax.legend()
fig.tight_layout()
fig.savefig(dataPathOut + "NumReads" + reg + ".pdf", bbox_inches = "tight") 
fig.savefig(dataPathOut + "NumReads" + reg + ".svg", bbox_inches = "tight")


# Region
reg = "Reg3"

# Load data
fIn1 = os.path.join(dataPathIn, reg + "ReadsBefore" + ".csv")
R1readsBefore = loadNumReads(fIn1)
fIn2 = os.path.join(dataPathIn, reg + "ReadsAfter" + ".csv")
R1readsAfter = loadNumReads(fIn2)
fIn3 = os.path.join(dataPathIn, reg + "ReadsCollapsed" + ".csv")
R1collapsed = loadNumReads(fIn3)

# Convert to array the keys and values of the dictionaries
filesBefore = np.array(list(R1readsBefore.keys()))
numReadsBefore = np.array(list(R1readsBefore.values()))

numReadsAfter = np.array(list(R1readsAfter.values()))
numReadsCollapsed = np.array(list(R1collapsed.values()))

# Do the plots
x = np.arange(numReadsBefore.shape[0])
width = 0.3

fig, ax = plt.subplots(figsize=(12,6))
# Barplots
before = ax.bar(x - 0.30, numReadsBefore, width, label='Before Processing',
				color = "#E9ECEF", edgecolor="black")
after = ax.bar(x, numReadsAfter, width, label='After Processing',
				color = "#ADB5BD", edgecolor="black")
collapsed = ax.bar(x + 0.30, numReadsCollapsed, width, label='Collapsed',
				color = "#343A40", edgecolor="black")

# Custom plot
ax.set_ylabel('Number of Sequences')
ax.set_title('Region 3')
ax.set_ylim(1e4, 1e6)
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(filesBefore, rotation=90)
#ax.legend()
fig.tight_layout()
fig.savefig(dataPathOut + "NumReads" + reg + ".pdf", bbox_inches = "tight") 
fig.savefig(dataPathOut + "NumReads" + reg + ".svg", bbox_inches = "tight")

plt.show()



##plt.bar(x, numReads, align='center', width=0.35)
##plt.bar(x, numReadsFilter, align='center', width=0.5)
##plt.yscale("log")
##plt.ylim(min(numReads) - 10000, max(numReads) + 10000)
##plt.xlabel("File name")
##plt.ylabel("Number of Reads")
##plt.legend((p1[0], p2[0]), ('BFilter', 'AFilter'))
##fig.savefig(dataPathOut + "reg1_pipeline.png")

### Reg2
##reads = LoadNumReads("./reg2.txt")
##
##readsFilter = LoadNumReads("./reg2_filters.txt")
##
### Convert to array the keys and values of the dictionary
##files = np.array(list(reads.keys()))
##numReads = np.array(list(reads.values()))
##
##filesFilter = np.array(list(readsFilter.keys()))
##numReadsFilter = np.array(list(readsFilter.values()))
##
### Do the plots
##fig = plt.figure(figsize=(17, 15))
##p1 = plt.bar(files, numReads, align='center', width=0.35)
##p2 = plt.bar(filesFilter, numReadsFilter, align='center', width=0.5)
##plt.yscale("log")
##plt.ylim(min(numReads) - 10000, max(numReads) + 10000)
##plt.xlabel("File name")
##plt.ylabel("Number of Reads")
##plt.legend((p1[0], p2[0]), ('BFilter', 'AFilter'))
##fig.savefig(dataPathOut + "reg2_pipeline.png")
##
### Reg3
##reads = LoadNumReads("./reg3.txt")
##
##readsFilter = LoadNumReads("./reg3_filters.txt")
##
### Convert to array the keys and values of the dictionary
##files = np.array(list(reads.keys()))
##numReads = np.array(list(reads.values()))
##
##filesFilter = np.array(list(readsFilter.keys()))
##numReadsFilter = np.array(list(readsFilter.values()))
##
### Do the plots
##fig = plt.figure(figsize=(17, 15))
##p1 = plt.bar(files, numReads, align='center', width=0.35)
##p2 = plt.bar(filesFilter, numReadsFilter, align='center', width=0.5)
##plt.yscale("log")
##plt.ylim(min(numReads) - 10000, max(numReads) + 10000)
##plt.xlabel("File name")
##plt.ylabel("Number of Reads")
##plt.legend((p1[0], p2[0]), ('BFilter', 'AFilter'))
##fig.savefig(dataPathOut + "reg3_pipeline.png")
##
##print("All done!")

