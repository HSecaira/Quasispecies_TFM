"""
Script that finds the empirical and theoretical abundance from files
	The empirical abundance is given by the abundance of each sequence in the sequencing files
	The theoretical abundance is given by the eigenvector centrality of the network
"""

import numpy as np;
import os, sys;
import timeit;

####################################################################
# Some functions
def loadExperimentalConditions(fIn):
	"""
	FUnction that reads all the names of experimental conditions in a file 
	Inputs: 
		fIn: file containing the names of the experimental conditions 
	Outputs:
		expConds: list of strings, containing the experimental conditions 
	"""

	# List to store the names
	expConds = []

	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			expConds.append(str(line.strip()))

	return expConds

def loadNodeIDS(fIn):
	"""
	Function that load the oldNodesIDs and newNodesIDs that were generated from mapping
	Inputs:
		fIn: file that contain the old and new node IDs
	Outputs:
		mapIDs: dictionary whose values are the oldNodesID and values are the newNodesID
	"""

	# Delcare structures
	mapIDs = {}

	# Iterate over files
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip header
			if line[0] != "#":
				oldNodeID = int(line.strip().split(",")[0])
				newNodeID = int(line.strip().split(",")[1])
	
				# Add to dictionary
				mapIDs[oldNodeID] = newNodeID

	return mapIDs

def loadGenotypesGN(fIn):
	"""
	Function that loads the genotypes of the grandNetwork as a dictionary
	Inputs:
		fIn:file containing the genotypes of the grand network
	Outputs:
		genotypesGN: dictionary whose values contains the IDs of genotypes/nodes of the grandNetwork and
					whose values contains the genotypes/sequences
	"""

	# Declare structures
	genotypesGN = {}

	# Iterate over file
	with open(fIn, "r") as fIn_:
		for ID, genotype in enumerate(fIn_):
			# Ignore empty lines
			if len(genotype) != 0:
				# Get the genotype
				genotype_ = str(genotype.strip())

				# Add to dictionary
				genotypesGN[ID] = genotype_

	return genotypesGN


def findEmpiricalAbundances(mapIDs, genotypesGN, seqsPathIn_, dataPathOut__, experimentalCondition):
	"""
	Function that fins the empirical abundance of a network from sequencing files 
	Inputs:
		mapIDs: dictionary whose values are the oldNodesID and values are the newNodesID
		genotypesGN: dictionary whose values contains the IDs of genotypes/nodes of the grandNetwork and
					whose values contains the genotypes/sequences
		seqsPathIn_: path that contains all sequencing files
		dataPathOut__: output path
		experimentalCondition: string containing the name of the experimental condition
	Outputs:
		empiricalAbundance: dictionary whose keys contains the oldNodesID and whose values contains the  
							empirical abundance
	"""

	# Declare structures
	empiricalAbundance = {}

	# Iterate over the oldIDs
	for oldID in mapIDs.keys():

		# Get the genotype
		genotype = genotypesGN[oldID]

		# Search the genotype in file using a command
		strOut = "grep -B 1 " + genotype + ' ' + seqsPathIn_ + experimentalCondition + ".fasta" + " > " + dataPathOut__ + "temp.csv" 
		os.system(strOut)

		# Open temp.csv file to extract the empirical abundance
		with open(dataPathOut__ + "/temp.csv", "r") as fInTemp:
			for line in fInTemp:
				if ">" in line:
					abundance = int(line.strip().split("-")[-1])

		# Add to dictionary
		empiricalAbundance[oldID] = abundance

		# Delete temp.csv
		strOut = "rm " + dataPathOut__ + "temp.csv"
		os.system(strOut)

	return empiricalAbundance

def loadEigenvectorCentrality(fIn):
	"""
	Function that retrieves the Eigenvector centrality of all nodes
	Inputs: 
		fIn: key-value type file
	Outputs:
		EigenCentrality: dictionary containing the nodeID (key) eigenvector centrality of each node (value)
	"""

	# List to store the eigenvector centrality 
	EigenCentrality = {}
	
	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip the header
			if line[0] != "#":
				nodeID = int(line.strip().split(",")[0])
				eiv = float(line.strip().split(",")[1])

				# Add to the dictionary
				EigenCentrality[nodeID] = eiv 

	return EigenCentrality

def findTheoreticalAbundances(mapIDs, EigenCentrality):
	"""
	Function that finds the theoretical abundance of a network from its eigenvector centrality
	Inpust:
		mapIDs: dictionary whose values are the oldNodesID and values are the newNodesID
		EigenCentrality: dictionary containing the eigenvector centralities of a network
	Outputs:
 		theoreticalAbundance: dictionary whose keys contains the newNodesID and whose values contains the  
							empirical abundance
	"""

	# Declare structures
	theoreticalAbundance = {}

	# Iterate over newIDs
	for newID in mapIDs.values():
		# If the ID exists in the eigenvector centrality. ID could not exist if the node
		# does not belong to the biggest component
		if newID in EigenCentrality.keys():
			# Get the theoretical abundance
			abundance = EigenCentrality[newID]
	
			# Add to dictionary
			theoreticalAbundance[newID] = abundance

	return theoreticalAbundance


def saveAbundances(fOut, empiricalAbundance, theoreticalAbundance, mapIDs, EigenCentrality):
	"""
	function that saves the empirical and theoretical abundances of a population 
	Inputs:
		fOut: file to save the Abundances
		empiricalAbundance: dictionary whose keys contains the oldNodesID and whose values contains the  
							empirical abundance
		theoreticalAbundance: dictionary whose keys contains the newNodesID and whose values contains the  
							empirical abundance
	Outputs:
		None
	"""

	# We must normalze the empirical abundance
	suma = 0
	for oldID, newID in mapIDs.items():
		# If the ID exists in the eigenvector centrality. ID could not exist if the node
		# does not belong to the biggest component
		if newID in EigenCentrality.keys():
			suma += empiricalAbundance[oldID]

	normEmpiricalAbundance = {}
	for oldID, newID in mapIDs.items():
		# If the ID exists in the eigenvector centrality. ID could not exist if the node
		# does not belong to the biggest component
		if newID in EigenCentrality.keys():
			normEmpiricalAbundance[oldID] = empiricalAbundance[oldID]/suma

	# Open output file
	with open(fOut, "w") as fOut_:
		# Write a header
		fOut_.write("#oldNodeID,EmpAbundance,newNodesID,TheoAbundance\n")

		# Iterate over mapIDs
		for oldID, newID in mapIDs.items():
			# If the ID exists in the eigenvector centrality. ID could not exist if the node
			# does not belong to the biggest component
			if newID in EigenCentrality.keys():
				empAb = normEmpiricalAbundance[oldID]
				theoAb = theoreticalAbundance[newID]
	
				# Write into file
				fOut_.write(str(oldID) + "," + str(empAb) + "," + str(newID) + "," + str(theoAb) + "\n")




###############################################################################

# Declare paths
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/"
GNdataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/GrandNetwork/"
seqsPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/NGS_analyses/out/Seqs_fasta_collapsed/";

dicFolder = {}; 
dicFolder["r1"] = "Reg1"; 
dicFolder["r2"] = "Reg2"; 
dicFolder["r3"] = "Reg3"; 

# Declare region
reg = "r1"
dataPathIn_ = dataPathIn + dicFolder[reg] + "_experiments/" + "SubNetworks/"
GNdataPathIn_ = GNdataPathIn + dicFolder[reg] + "/MetaInfo/"
seqsPathIn_ = seqsPathIn + dicFolder[reg] + "/"

# Star counting execution time
start = timeit.default_timer()

# Load experimental conditions
fInExpConds = os.path.join(dataPathIn_, "filesNames" + ".csv")
# Load the names of the experimental conditions
expConds = loadExperimentalConditions(fInExpConds)

# Load genotypes of the grand network as a dictionary
fInGN = os.path.join(GNdataPathIn_, "nodesID_" + reg + ".csv")
genotypesGN = loadGenotypesGN(fInGN)

# Iterate over experimental conditions
for experimentalCondition in expConds:
	print("At", experimentalCondition, "\n")
	# Redefine paths
	dataPathIn__ = dataPathIn_ + experimentalCondition + "/"
	dataPathOut__ = dataPathIn__

	# Load the mapped node IDs
	fIn = os.path.join(dataPathIn__, "mapNodesIDs" + ".csv")
	mapIDs = loadNodeIDS(fIn)

	# Get the empirical abundance from sequencing files
	empiricalAbundance = findEmpiricalAbundances(mapIDs, genotypesGN, seqsPathIn_, dataPathOut__, experimentalCondition)

	# Get the theoretical abundance
	fInEiv = os.path.join(dataPathIn__, "ResultsTopology/" + "EigenvectorCentrality_GrandNetwork" + ".csv")
	EigenCentrality = loadEigenvectorCentrality(fInEiv)
	theoreticalAbundance = findTheoreticalAbundances(mapIDs, EigenCentrality)

	# Save empirical and theoretical abundances from eigenvector centrality
	fOut = os.path.join(dataPathOut__, "Abundances/" + "AbundancesEmpvsTheo" + ".csv")
	saveAbundances(fOut, empiricalAbundance, theoreticalAbundance, mapIDs, EigenCentrality)

# Get execution time
stop = timeit.default_timer()
total_time = stop - start	
# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)	
sys.stdout.write("Total running time: %d:%d:%d. s\n" % (hours, mins, secs))



