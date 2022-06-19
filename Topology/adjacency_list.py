# This script will generate an adjacency list given an edeges file

#use grep to find all the occurrences of a vertex
#and redirect all the occurrences to a temporary file
#the adjacency list will have the following structure: a dictionary
#whose key is the node ID and the values is a list/array of all the 
#connections of the key
#to avoid repetitions in keys, we will check if the key already exists


# Imports
import numpy as np;
import os, sys;

# Some functions

def UniqueConnections(fileIn):
	"""
	Function that appends all the unique connections of a node
	Inputs:
		fileIn: CSV file containing all the connections of a node
	Outputs:
		values: list containing the unique connections
	"""
	values = []
	with open(fileIn, "r") as fIn:
		for line in fIn:
			# Append only unique connections
			if line.strip() not in values:
				values.append(line.strip())
			else:
				pass
	return values


# Define paths
dataPathIn = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/NetExamples/NetSamples/";
dataPathOut = "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/";
dicFolder = {};
dicFolder["r1"] = "Reg1/";
dicFolder["r2"] = "Reg2/";
dicFolder["r3"] = "Reg3/";

# Loading data files
regionName = "r1";
dataPathIn_ = dataPathIn; 
dataPathOut_ = dataPathOut + dicFolder[regionName];
adjacencyList = {} 
#keys = []
#values = []

# Iterate over the file that contains all the edges
fIn = os.path.join(dataPathIn, "1000Nodes" + ".csv")
fOut = os.path.join(dataPathOut_, "Adjacency_list_1000" + ".csv")
print("\nCreating the adjaceny list as a dictionary...\n")
with open(fIn, "r") as fIn_:
	for line in fIn_:
		nodes = line.strip().split(", ")
		# Append only unique node keys
		if nodes[0] not in adjacencyList.keys():
			#adjacencyList.append(nodes[0])
			# Find all the occurrences of a node, using grep (shell command) and redirect to a temporary file
			command = "grep -w " + nodes[0] + " " + fIn + " | sed 's/" + nodes[0] \
			+ "//' | tr -d ',' | tr -d ', '" + " > temp.csv"
			# Execute the command
			os.system(command)
			# Iterate over all the connections of a node
			fIn1 = os.path.join("./temp.csv")
			adjacencyList[nodes[0]] = ",".join(UniqueConnections(fIn1))
		elif nodes[1] not in adjacencyList.keys():
			# Find all the occurrences of a node, using grep (shell command) and redirect to a temporary file
			command1 = "grep -w " + nodes[1] + " " + fIn + " | sed 's/" + nodes[1] \
			+ "//' | tr -d ',' | tr -d ', '" + " > temp1.csv"
			# Execute the command
			os.system(command1)
			# Iterate over all the connections of a node
			fIn2 = os.path.join("./temp1.csv")
			adjacencyList[nodes[1]] = ",".join(UniqueConnections(fIn2))
		else:
			pass
# Delete temp.csv file
os.system("rm ./temp.csv")
os.system("rm ./temp1.csv")

# Write the dictionary to an output file
print("\nSaving the adjacency list into a file!\n")
with open(fOut, "w") as fOut_:
	# A header for the file
	fOut_.write("#NodeID:Connections\n")
	#Iterate over the dictionary
	for key in adjacencyList:
		fOut_.write(key + ":" + str(adjacencyList[key]) + "\n")

		













