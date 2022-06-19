"""
Scripts that contain functions that help in plots
"""

# Imports
import numpy as np;
import os, sys;

def RetrieveDegrees(fIn):
	"""
	Function that retireves all the degrees of all nodes
	Inputs:
		fIn: file containing the adjacency list
	Outputs:
		degrees: array with all the degrees of all nodes
	"""
	# List to store all the degrees
	degrees = []
	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# SKip the header
			if line[0] != "#":
				degrees.append(int(line.strip().split(",")[1]));
	return np.array(degrees)

def RetrieveKDegrees(fIn):
	"""
	Function that retrieves the degrees of all nodes
	Inputs:
		fIn: file containig the necessary information
	Outputs:
		Kdegrees: array with the degrees of all nodes
	"""

	# List to store the degrees
	Kdegrees = []

	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip the header
			if line[0] != "#":
				Kdegrees.append(int(line.strip().split(",")[0]))
	return np.array(Kdegrees)

def RetrievePk(fIn):
	"""
	Function that retrieves the Pk (i.e. rank/num_nodes) of an adjacency list
	INputs: 
		fIn: file containing the adjacency list
	Outputs:
		Pk: array with the Pk values for each node
	"""

	# List to store all the Pks
	Pk = []
	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# SKip the header
			if line[0] != "#":
				Pk.append(float(line.strip().split(",")[3]))
	return np.array(Pk)
	

def RetrieveNearestNeighborDegree(fIn):
	"""
	Function that retrieves the K nearest neighbor degree of all nodes
	Inputs:
		fIn: key-value type file
	Outputs:
		KNNeighbors: dictionary containing the K Nearest Neighbor for each degree
	"""
	# List to store all the neighbor degrees degrees
	KNNeighbors = {}
	temp = {}

	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# SKip the header
			if line[0] != "#":
				degree = int(line.strip().split(",")[0])
				Knn = float(line.strip().split(",")[1])

				temp[degree] = Knn

	#Sort by degrees, i.e. keys
	KNNeighbors = dict(sorted(temp.items(), key = lambda x:x[0]))

	return KNNeighbors

def RetrieveNearestNeighborPerDegree(fIn):
	"""
	Function that retrieves the K nearest neighbor degree of all degrees
	Inputs:
		fIn: key-value type file
	Outputs:
		KNNeighbors: list of lists containing the K Nearest Neighbor for each degree
	"""
	# List to store all the neighbor degrees degrees
	KNNeighbors = {}
	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# SKip the header
			if line[0] != "#":
				# Get the degree
				degree = int(line.strip().split(":")[0])
				# Get Knn as a list and convert to float
				Knn = line.strip().split(":")[1]
				Knn_floats = [float(k) for k in Knn.split(",")]
				# Add to dictionary
				KNNeighbors[degree] = Knn_floats
	return KNNeighbors



def RetrieveKClustering(fIn):
	"""
	Function that retrieves the K clustering coefficient of all nodes
	Inputs:
		fIn: key-value type file
	Outputs:
		KClustCoeff: dictionary containing the KCLustering coefficient for each degree
	"""
	# List to store all the K clustering coefficients
	KClustCoeff = {}
	temp = {}
	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# SKip the header
			if line[0] != "#":
				degree = int(line.strip().split(",")[0])
				clust = float(line.strip().split(",")[1])

				temp[degree] = clust

	#Sort by degrees, i.e. keys
	KClustCoeff = dict(sorted(temp.items(), key = lambda x:x[0]))

	return KClustCoeff

def RetrieveEigenvectorCentrality(fIn):
	"""
	Function that retrieves the Eigenvector centrality of all nodes
	Inputs: 
		fIn: key-value type file
	Outputs:
		EigenCentrality: array containing the eigenvector centrality of each node 
	"""

	# List to store the eigenvector centrality 
	EigenCentrality = []

	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip the header
			if line[0] != "#":
				EigenCentrality.append(float(line.strip().split(",")[1]))

	return np.array(EigenCentrality)

def RetrieveEigenvectorCentralityDegree(fIn):
	"""
	Function that retrieves the Eigenvector centrality of all nodes
	Inputs: 
		fIn: key-value type file
	Outputs:
		EigenCentrality: dictionary containing the eigenvector centrality of each node (key) and degree (value)
	"""

	# List to store the eigenvector centrality 
	EigenCentrality = {}
	temp = {}

	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip the header
			if line[0] != "#":
				EigenCentrality = float(line.strip().split(",")[1])
				degree = int(line.strip().split(",")[2])

				temp[EigenCentrality] = degree 

	# Order the dictionary by degree, i.e. by values
	EigenCentrality = dict(sorted(temp.items(), key = lambda x:x[1]))

	return EigenCentrality

def RetrieveEigenvectorCentrality_alternative(fIn):
	"""
	Function that retrieves the Eigenvector centrality of all nodes
	Inputs: 
		fIn: key-value type file
	Outputs:
		EigenCentrality: array containing the eigenvector centrality of each node 
	"""

	# List to store the eigenvector centrality 
	EigenCentrality = []

	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip the header
			if line[0] != "#":
				centrality = float(line.strip().split(",")[1])
				if centrality >= 1.0e-13:
					EigenCentrality.append(centrality)

	return np.array(EigenCentrality)


def RetrieveBetweennessCentrality(fIn):
	"""
	Function that retrieves the Betweenness centrality of all nodes 
	Inputs: 
		fIn: key-value type file
	Outputs:
		BetCentrality: array containing the Betweenness centrality of each node 
	"""
	# List to store the Betweenness centrality 
	BetCentrality = []

	# Iterate over the file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip the header
			if line[0] != "#":
				BetCentrality.append(float(line.strip().split(",")[1]))

	return np.array(BetCentrality)

def RetrieveBetweennessCentralityDegree(fIn):
	"""
	Function that retrieves the Betweenness centrality and degree of all nodes 
	Inputs:
		fIn: key-value type file
	Outputs: 
		BetCentralityDegree: dictionary containing the Betweenness Cetrality (key) and the degree (value)
	"""

	# Dictionary to store the betweenness centrality and degree
	BetCentralityDegree = {}
	temp = {}

	# Iterate over the file 
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip header
			if line[0] != "#":
				bet = float(line.strip().split(",")[1])
				deg = int(line.strip().split(",")[2])

				temp[bet] = deg

	# Order the dictionary by degrees, i.e. by values
	BetCentralityDegree = dict(sorted(temp.items(), key = lambda x:x[1]))

	return BetCentralityDegree

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


def RetrieveEigenvectorCentralityComponents(fIn):
	"""
	Function that retrieves the eigenvector centralities and degrees of all components in a network
	Inputs:
		fIn: file containing the eigenvector centralities of all components in the network
	Outputs:
		eigenCentralityComponents: dictionary containing the eigenvector centralities of all components
	"""

	# Dictionary to store the eigenvector centralities
	eigenCentralityComponents = {}
	oldID = -1

	# Iterate over file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip the header
			if line[0:2] != "##":
				# If line starts with only one "#" is the component ID
				if line[0] == "#":
					# Key of the dictionary is the component ID
					componentID = int(line.strip().split(",")[0][1:])

					# For the first identifier
					if oldID == -1:
						# Update oldID
						oldID = componentID
						# Initialize list 
						tempCentralityDegree = []
					else:
						if oldID != componentID:
							# Append to the list
							eigenCentralityComponents[oldID] = tempCentralityDegree
							# Update oldID
							oldID = componentID
							# Initialize list
							tempCentralityDegree = []
				else:
					nodeCentrality = float(line.strip().split(",")[1])
					degree = int(line.strip().split(",")[2])
					tempCentralityDegree.append([nodeCentrality, degree])

		# Add the last component
		eigenCentralityComponents[componentID] = tempCentralityDegree

	return eigenCentralityComponents

def RetrieveMSE(fIn):
	"""
	Function that retrieves the MSE from betweenness estimation
	Inputs: 
		fIn: file containing the MSE for all nodes used in estimation
	Outputs:
		MSEBetweenness: dictionary containing the MSE (values) for each node used in estimation (keys)
	"""
	# Declare structures
	#MSE = {}
	fracNodes = []
	MSE = []
	MSESTD = []
	# Open file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip header
			if line[0] != "#":
				# Get values
				fracNodes.append(float(line.strip().split(",")[0]))
				MSE.append(float(line.strip().split(",")[1]))
				MSESTD.append(float(line.strip().split(",")[2]))
				#MSE[numNodes] = error
	return np.array(fracNodes), np.array(MSE), np.array(MSESTD)

def RetrieveRichClubCoeff(fIn):
	"""
	Function that retrieves the normalized rich club coefficient of the biggest component of a network
	Inputs:
		fIn: file containing the MSE for all nodes used in estimation
	Outputs:
		rcCoeff: dictionary containg the normalized richClubCoeff (values) for all degrees (key)
	"""

	# Declare structures
	rcCoeff = {}
	rcCoeff_ = {}

	# Open file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip header
			if line[0] != "#":
				# Get values
				degree = int(line.strip().split(",")[0])
				rc = float(line.strip().split(",")[1])
				rcCoeff[degree] = rc 

	# Order the dictionary by centrality, i.e. by keys
	rcCoeff_ = dict(sorted(rcCoeff.items(), key = lambda x:x[0]))

	return rcCoeff_

def RetrieveMetric(fIn):
	"""
	Function that retrieves the mean degree and the characteristic path lenght of a file
	Inputs:
		fIn: file containing some topology metric information
	Outputs:
		meanDegree: float containing the mean degree 
		charPathLenght: float containing the characteristic path lenght
	"""

	# Open file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# SKip header
			if line[0] != "#" and len(line) != 0:
				# Get values
				density = float(line.strip().split(",")[0])
				meanDegree = float(line.strip().split(",")[1])
				meanDegreeSTD = float(line.strip().split(",")[2])
				assortativity = float(line.strip().split(",")[3])
				globalClustering = float(line.strip().split(",")[5])
				charPathLenght = float(line.strip().split(",")[7])
				charPathLenghtSTD = float(line.strip().split(",")[8])
				eigenvalue = float(line.strip().split(",")[9])


	return density, meanDegree, meanDegreeSTD, assortativity, globalClustering, charPathLenght, charPathLenghtSTD, eigenvalue


def RetrieveClustCoeffKnn(fIn):
	"""
	Function that retrieves the local clustering coefficient/Knn of all nodes in a network
	Inputs: 
		fIn: file containing  the local clustering coefficient
	Outputs: 
		clustCoeffDict: dictionary whose keys are degrees and values are lists containing
		the nodes with a clustering coefficient/Knn associated to a degree
	"""

	# Declare some variables
	temp = {}
	clustCoeffDict = {}

	# Open file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip header
			if line[0] != "#":
				clustCoeff = float(line.strip().split(",")[1])
				degreeClustCoeff = int(line.strip().split(",")[2])
	
				# Add to dictionary
				if degreeClustCoeff not in temp.keys():
					temp[degreeClustCoeff] = [clustCoeff]
				else:
					temp[degreeClustCoeff].append(clustCoeff)

	# Sort dictionary by degrees, i.e. by keys
	clustCoeffDict = dict(sorted(temp.items(), key = lambda x:x[0]))

	return clustCoeffDict


def RetrieveAbundances(fIn):
	"""
	Function that retrieves the empirical and theoretical abundances of the biggest component of a network 
	Inputs: 
		fIn: file containing the abundances  
	Outputs: 
		empAbundance: dictionary containing the grandNetwork node ID (key) and the relative abundance (value)
		theoAbundance: dictionary containing the subnetwork ID (key) and the relative abundance 
	"""

	# Declare structures
	empAbundance = {}
	theoAbundance = {}

	# Iterate over file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# SKip header
			if line[0] != "#":
				GNID = int(line.strip().split(",")[0])
				GNAb = float(line.strip().split(",")[1])
				SNID = int(line.strip().split(",")[2])
				SNAb = float(line.strip().split(",")[3])

				# Add to dictionaries
				empAbundance[GNID] = GNAb
				theoAbundance[SNID] = SNAb
				
	return empAbundance, theoAbundance













