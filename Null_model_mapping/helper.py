"""
	Script that contains function that help for sequence manipulation
"""

# Imports
import numpy as np;
import os, sys;

def LoadFasta(fIn):
	"""
	Function that loads a fasta file
	Inputs:
		fIn: fasta file containing the sequences
	Outputs:
		d: dictionary containing the sequences in which the keys re the IDs and the values are the sequences
	"""
	# Dictionary that will store the IDs and the seqs
	d = {}
	oldID = ''

	# Read the file
	with open (fIn, "r") as fIn_:
		for line in fIn_:
			# Skip header
			if line.strip()[0] != "#":
				# If the line starts with '>' is an identifier
				if line.strip()[0] == ">":
					newID = line.strip()[1:]
	
					# For the first identifier
					if oldID == '':
						# Update oldID
						oldID = newID
						# Initialize the seq
						seq = ''
					else:
						if oldID != newID:
							# Append to the dictionary the seq that spans over multiple lines
							d[oldID] = seq.upper()
							# Update oldID
							oldID = newID
							# Initialize the seq
							seq = ''
				else:
					seq = seq + line.strip()
		# Add the last seq
		d[newID] = seq.upper()

	return d

def SplitSeq(seq, ORF, gene):
	"""
	Function that splits a sequence into codons according to some positions of the protein coding genes
	Inputs:
		seq: string containing the sequence
		ORF: dictionary containing the positions of the protein-cding genes
		gene: string indicating the region/gene to split
	Outputs: 
		codons: dictionary containing the nucleotides separated as codons. Keys is the number of the codon and value is the codon
		codonsPos: dictionary containing the positions of codons with respect to seq. Keys is the number of the codon and value is the position of the codon

	"""

	# Dictionaries to store the codons
	codons = {}
	codonsPos = {}

	# Get the start and end codon positions for the gene
	startPos = ORF[gene][0]
	endPos = ORF[gene][1] + 1 

	# Get the sequence of the gene
	gene = seq[startPos:endPos]

	# Get the codons
	codonNumber = 0
	for i in range(0, len(gene), 3):
		codonNumber += 1
		codons[codonNumber] = gene[i:i+3]
		codonsPos[codonNumber] = [startPos + i, startPos + i+2]

	return codons, codonsPos

def ExtractSeqs(seq, refSeqsPos, region):
	"""
	Function that extract a sequence of a region from another sequence given a position
	Inputs:
		seq: string containing the sequence
		refSeqsPos: dictionary containing the positions of the regions
		region: string indicating the region for which to extract the sequence
	Outputs:
		extractSeq: string containing the sequence extracted
	"""

	# Get the start and end positions of the region
	startPos = refSeqsPos[region][0]
	endPos = refSeqsPos[region][1] + 1

	#  Extract the sequence
	extractSeq = seq[startPos:endPos]

	return extractSeq

def ExtractSequence(seq, startPos, endPos):
	"""
	Function that extract a sequence of a region given a position
	Inputs:
		seq: string containing the sequence
		startPos: int indicating the starting position
		endPos: int indicating the ending position
	Outputs:
		subSeq: string containing the sequence extracted
	"""

	# Extract the sequence
	subseq = seq[startPos:endPos + 1]

	return subseq

def ExtractSeqCodons(seq, startPos, endPos):
	"""
	Function that extracts a sequence given a position and its codons and codons position
	Inputs:
		seq: string containing the sequence
		startPos: int indicating the starting position
		endPos: int indicating the ending position
	Outputs:
		codons: dictionary containing the nucleotides separated as codons. Keys is the number of the codon and value is the codon
		codonsPos: dictionary containing the positions of codons with respect to seq. Keys is the number of the codon and value is the position of the codon

	"""

	# Extract the sequence
	subseq = seq[startPos:endPos + 1]

	# Dictionaries to store the codons
	codons = {}
	codonsPos = {}

	# Get the codons
	codonNumber = 0
	for i in range(0, len(subseq), 3):
		codons[codonNumber] = [subseq[i:i+3], [startPos + i, startPos + i + 1, startPos + i + 2]]
		codonNumber += 1

	return codons

def CheckCodonPosition(seq, startPos, codonsPos, startPosSequenced, endPosSequenced):
	"""
	Function that checks the starting position of a sequence with respect to codon positions
	Inputs:
		seq: string containing the sequence sequenced
		startPos: int indicating the starting position of the sequence according to the reference genome
		codonsPos: dictionary containing the positions of codons with respect to the reference genome. Keys is the number of the codon and value is the position of the codon
		startPosSequenced: int indicating the starting position of the seq
		endPosSequenced: int indicating the ending position of the seq
	Outputs:
		codons_, codonsPos_: dictionaries containing the codon number as key and the codon or position as values, respectively.

	"""

	# Iterate over the codons position
	for k, v in codonsPos.items():
		# Check if startPos in in the codon
		if startPos == v[0]:
			print("\tThe start position in the codons is", v[0] )
			print("\tThe region is in a codon starting position!")
			codons_ = ExtractSeqCodons(seq, startPosSequenced, endPosSequenced)
			return codons_
		elif startPos == v[0] + 1:
			print("\tThe start position in the codons is", v[0] + 1)
			print("\tThe region is one position ahead of a codon starting position!")
			print("\tExtracting the codons starting at the next two nucleotides from startPos")
			codons_ = ExtractSeqCodons(seq, startPosSequenced + 2, endPosSequenced)
			return codons_
		elif startPos == v[0] - 1:
			print("\tThe start position in the codons is", v[0] -1)
			print("\tTHe region is one position behind a codon starting position!")
			print("\tExtracting the codons starting at the next nucleotide from startPos")
			codons_ = ExtractSeqCodons(seq, startPosSequenced + 1, endPosSequenced)
			return codons_


def Translation(seq, trans_table, leaky = False, stopCodon = False):
	"""
	Function that translates a nucleotide sequence according to a translation table
	Inputs:
		seq: string containing the sequence to translate. Beware! the sequence must be in the correct ORF
			 otherwhise this does not make sense
		trans_table: dictionary containing the genetic code for trasnlation
		leaky: boolean indicating the leaky codon, only for A1 protein in QBeta Bacteriophage
	Outputs:
		prot: string containing the translated seuence without the stip codon
	"""

	# Amino acid sequence
	prot = ''
	number = 0

	# Iterate over the sequence
	for i in range(0, len(seq), 3):
		codon = seq[i:i + 3]

		# Translate using the genetic code
		for k, v in trans_table.items():
			if codon in v:
				# Add condition for leaky
				if codon == "TGA" and leaky == True:
					aminoacid_leaky = "W"
					prot += aminoacid_leaky
				else:
					prot += k

	if stopCodon:
		return prot
	else:
		return prot[0:len(prot) - 1]


def TranslationCodons(codons, trans_table, leaky = False):
	"""
	Function that translates a nucleotide sequence according to a translation table given some codons
	Inputs:
		codons: dictionary containing the codons info.
		trans_table: dictionary containing the genetic code for trasnlation
		leaky: boolean indicating the leaky codon, only for A1 protein in QBeta Bacteriophage
	Outputs:
		prot: string containing the translated seuence without the stip codon
	"""

	# Amino acid sequence
	prot = ''

	# Iterate over the sequence
	for k, v in codons.items():
		codon = v[0]

		# Check that the len of codon is appropiate
		if len(codon) == 3:
			# Translate using the genetic code
			for k_, v_ in trans_table.items():
				if codon in v_:
					# Add condition for leaky
					if codon == "TGA" and leaky == True:
						aminoacid_leaky = "W"
						prot += aminoacid_leaky
					else:
						prot += k_
		else:
			print("\tThe codon", k, " has a lenght", len(codon))

	# Check if there is an stop codon at the end of the protein
	#if prot[-1] == "*":
	#	print("\tReturning protein witout stop codon at the end")
	#	return prot[0:len(prot) - 1]
	#else: 
	#	return prot
	return prot

def proteinInfo(proteinSeq, startPos):
	"""
	Function that creates a dictionary of the aminoacids of a protein sequence 
	Inputs:
		proteinSeq: string containing the protein sequence 
		startPos: int indicating the starting position of the sequence with respect to the sequenced sequence (0-based)
	Outputs:
		protInfo: dictionary containing the protein info
			Key: int indicating the amino acid number 
			Values:
				value[0]: string containing an amino acid
				value[1]: string containing the position of the aminoacid with respect to the sequenced sequence (0-based)
	"""

	# Dictionary to store the protein Info
	protInfo = {}

	# Iterate over the protein sequence
	for i in range(len(proteinSeq)):
		protInfo[i] = [proteinSeq[i], startPos + i]

	return protInfo


def saveCodonsInfo(fOut, codonInfo, protSeq, startIndex_genome, startIndex_prot, startIndex_protaa, replicase_R2 = False):
	"""
	Function that saves the codon of a sequenced region into a file
	Inputs:
		fOut: string indicating the name of the file to save the info
		codonInfo: dictionary containing the codons information of a sequence sequenced
			key: int indicating the codon number (0-based)
			values: list
				value[0]: string containing the codon
				value[1]: list of ints containing thte positions of the codon (0-based)
		protSeq: string containing the protein sequence
		startIndex_genome: int indicating the starting position with respect to the reference Genome
		startIndex_prot: int indicating the starting position with repespect to the protein
		startIndex_protaa: int indicating the starting position with respect to the protein aminoacids
		replicase_R2: boolean indicating if replicase of R2 is the protein
	Outputs:
		None
	"""

	# Open the file
	with open(fOut, "w") as fOut_:
		# Write a header
		fOut_.write("#codonNumber,codon,codonPosSequenced,aminoacid,codonPosGenome,codonPosProtein,aaPosProtein\n")

		# Iterate over codonInfo
		counter = 0
		j = 0
		for k, v in codonInfo.items():
			positions = '-'.join(str(pos) for pos in v[1])

			posProt = startIndex_protaa + counter

			if replicase_R2:
				positionRefProt = []
				positionRefGenome = []

				for num, pos in enumerate(v[1]):
					positionRefProt.append(j + num)
					positionRefGenome.append(startIndex_genome + j + num)
				j += len(v[1])

				posRP = '-'.join(str(pos) for pos in positionRefProt)
				posRG = '-'.join(str(pos) for pos in positionRefGenome)

			else:
				positionRefGenome = [startIndex_genome + pos for pos in v[1]]
				posRG = '-'.join(str(pos) for pos in positionRefGenome)
				positionRefProt = [startIndex_prot + pos for pos in v[1]]
				posRP = '-'.join(str(pos) for pos in positionRefProt)


			# Save with \t as separator so the list is the only one that contains ","
			if counter < len(protSeq):
				fOut_.write(str(k) + "," + str(v[0]) + "," + str(positions) + "," + protSeq[counter] + "," + str(posRG) + "," + str(posRP) + "," + str(posProt) + "\n")
			else:
				fOut_.write(str(k) + "," + str(v[0]) + "," + str(positions) + "," + "NaN" + "," + str(posRG) + "," + str(posRP) + "," + str(posProt) + "\n")
			counter += 1


def loadCodonsInfo(fIn):
	"""
	Function that load the codon infor saved into a file by saveCodonsInfo()
	Inputs:
		fIn: input file containing the codons info 
	Outputs:
		codonsInfo: dictionary containing the codons info
		Keys: strings conatining the codon number
		Values: 
			value[0]: string containing the codon 
			value[1]: numpy array of ints containing the positions of codons with respect to the reference sequence 
			value[2]: string containing the aminoacid of the codon


	"""

	# Dictionary to save the info
	codonsInfo = {}

	# Open the file
	with open(fIn, "r") as fIn_:

		# Iterate over lines
		for line in fIn_:
			# Ignore header
			if line[0] != "#":
				# Ignore the \n at the end
				info = line.strip()
				# Weird tricks to get the info ;)
				codonNumber = info.split(",")[0]
				codon = info.split(",")[1]
				pos = info.split(",")[2].split("-")
				codonsPos = np.array([int(pos[0]), int(pos[1]), int(pos[2])])
				aminoacid = info.split(",")[3]

				# Add info to a dictionary
				codonsInfo[codonNumber] = [codon, codonsPos, aminoacid]

	return codonsInfo


def convertRefSeqs(refSeqs):
	"""
	Function that converts sequences into an array of positions 
	Inputs: 
		refSeqs: dictionary containing the sequences
	Outputs: 
		refSeqs: dictionary containing the array of positions
	"""

	# Dictionary to store the positions
	refSeqsPos = {}

	# Iterate over the dictionary
	for k, v in refSeqs.items():
		lenSeq = len(v)
		# Create the array
		pos = np.arange(0, lenSeq)

		# Add to dictionary
		refSeqsPos[k] = pos 

	return refSeqsPos

def checkPositionCodonORF(randomPos, codonsInfo):
	"""
	Function that checks if the random position chosen from the sequence
	corresponds to a codon position 
	Inputs: 
		randomPos: int array of size 1 containing the position selected at random 
		codonsInfo: dictionary containing the codons info
		Keys: strings conatining the codon number
		Values: 
			value[0]: string containing the codon 
			value[1]: numpy array of ints containing the positions of codons with respect to the reference sequence 
			value[2]: string containing the aminoacid of the codon
	Outputs:
		index: int array of size 1 containing the index of randomPos with respect to its codon 
		codon: string containing the codon
	"""

	# Iterate over codonsInfo
	for k, v in codonsInfo.items():
		codonsPos = v[1]
		codon = v[0]
		aminoacid = v[2]

		# Check if the length of codon is valid (important because some regions contain incomplete codons!)
		if len(codon) == 3:
			# Check if randomPos is in the current position
			if randomPos in codonsPos:
				# Get the index. This index also corresponds to the index in the current codon
				index = np.where(codonsPos == randomPos)[0]
				return index, codon, aminoacid

	print("\trandomPos not found! i.e. does not belong to the ORF!")
	# ask luino what to return
	return np.array([-1]), -1, -1

def eigenvectorCentrality(network):
	"""
	Function that calculates the eigenvector centrality of a network
	Inputs: 
		network: numpy array of ints containing the adjacency matrix of the network 
	Outputs: 
		eigenCentrality: numpy array of ints containing the eigenvector centrality, i.e. the eigenvector 
						 associated to the largest eigenvalue
	"""

	# Get eigenvalues and eigenvectors of the adjacency matrix
	evals, evecs = np.linalg.eig(network)

	# Get the largest eigenvalue
	max_ = np.argmax(evals)
	max_eigenvalue = evals[max_]

	# Get the maximum eigenvector
	max_eigenvector = evecs[:, max_]

	# Normalize the eigenvector by the sum of its values, so it is between 0 and 1 (probability)
	max_eigenvector_ = max_eigenvector/np.sum(max_eigenvector)

	return max_eigenvector_

def vecNormProb(x):
	"""
	Function that calculates the norm a vector such that it is a probability distribution 
	Inputs:
		x: array of float containing a vector  
	Outputs:
		xNorm: array of floats containing the vector  normalized 
	"""

	norm = np.sum(x)

	xNorm = x/norm 

	return xNorm

def EivCentProb(network, max_iter=100, tol=1.0e-6):
	"""
	Function that calculates the eigenvector centrality of a network
	Inputs:
		network: numpy matrix of floats containing the adjacency matrix of the network  
	Outputs:
		eivCent: array of floats containing the eigenvector ccentrality of the network     
	"""

	# Initial vector
	initVec = np.ones(network.shape[0])

	# Normalize vector
	x0 = vecNormProb(initVec)

	for i in range(max_iter):

		x1 = np.dot(network, x0)

		# Calculate eigenvalue
		maxVal = np.argmax(x1)
		eivalue = x1[maxVal]/x0[maxVal]

		# Normalize vector
		x1 = vecNormProb(x1)

		# Check convergence
		if np.sum(np.abs(x1 - x0)) < network.shape[0] * tol:  
			return x1

		# Update vector
		x0 = x1

		if i == max_iter:
			print("Convergence not reached, increase the number of iterations!")
			return x1

def probSilentMutations(index, aminoacid, codonSelected, trans_table, network_mutations, centrality):
	"""
	Function that calculates the probability of a silent mutation given the index in the codon
	Inputs:
		index: int containing the index at which a mutation can occur 
		aminoacid: string containing the amino acid
		codonSelected: string containing the codon
		trans_table: dictionary cotnaining the genetic code
		network_mutations: dictionary containing the adjacency matrices of the neutral network for an aminoacid
			Key: string indicating the aminoacid 
			Values: numpy array of integers containing the adjacency matrix of the network 
		centrality: boolean indicating if eigenvector centrality is considered for probability of silent mutations calculations, for all codons
	Outputs:
		probabilitySilentMutation: float indicating the probability of the silent mutation
	"""

	# List to store nucleotides at index that codify for the same aminoacid
	nuclSilentMutations = []
	# Set the probability
	probabilitySilentMutation = 0

	# Get all the codons of the amino acid
	allCodons = trans_table[aminoacid]
	# Get all alternative nucl that codify for the same aminoacid at index (i.e. silent mutation )
	for codon in allCodons:
		# Get only the nucleotide that is distinct to codonReal[index] (i.e. the original nucl to mutate!)
		if codon[index] not in nuclSilentMutations:
			nuclSilentMutations.append(codon[index])		
	# If no nucl are allowed to change at index
	if len(nuclSilentMutations) == 1:
		return probabilitySilentMutation
	else:
		# Get the eigenvector centrality for codon
		eigenvector = EivCentProb(network_mutations[aminoacid])

		# Check if codon has some interesting properties for silent mutations (i.e. its neutral network is not a clique)
		if aminoacid in ["L", "R"]:
			if index == 0:
				probabilitySilentMutation = (eigenvector[1] + eigenvector[4] + eigenvector[3] + eigenvector[5])*(1/3)
			elif index == 2:
				probabilitySilentMutation = ((np.sum(eigenvector[0:4]))*(1)) + ((np.sum(eigenvector[4:])*(1/3)))
		elif aminoacid == "S":
			if index == 0:
				probabilitySilentMutation = 0 
			elif index == 1:
				probabilitySilentMutation = 0
			elif index == 2:
				#print("\tCodon Selected.......", codonSelected)
				if codonSelected in ["TCT", "TCC", "TCA", "TCG"]:
					probabilitySilentMutation = ((np.sum(eigenvector[0:4]))*(1))
				elif codonSelected in ["AGT", "AGC"]:
					probabilitySilentMutation = (1/3)
		elif aminoacid == "*":
			# It seems that the power method to calculate the eigenvector centrality
			# of the stop codon never converges if the vector is normalized as a 
			# probability distribution!!
			if eigenvector == None:
					eigenvector = eigenvectorCentrality(network_mutations[aminoacid])
			if index == 1:
				probabilitySilentMutation = (eigenvector[0] + eigenvector[2])*(1/3)
			elif index == 2:
				probabilitySilentMutation = (np.sum(eigenvector[0:2]))*(1/3)
		# all other aminoacids
		else:
			if centrality:
				probabilitySilentMutation = np.sum(eigenvector) * (len(nuclSilentMutations) - 1)/3
			else:
				probabilitySilentMutation = (len(nuclSilentMutations) - 1)/3

	return probabilitySilentMutation 


def probSilentMutationsAminoacid(aminoacid, codon, trans_table, network_mutations, centrality = False):
	"""
	Function that calculates the probability of a silent mutation of an aminoacid in all its codon positions
	Inputs:
		aminoacid: string indicating the minoacid
		codon: string indicating the codon selected
		trans_table: dictionary cotnaining the genetic code
		network_mutations: dictionary containing the adjacency matrices of the neutral network for an aminoacid
			Key: string indicating the aminoacid 
			Values: numpy array of integers containing the adjacency matrix of the network
		centrality: boolean indicating if eigenvector centrality is considered for probability of silent mutations calculations, for all codons
	Outputs:
		probabilitySilentMutations: arrays of ints containing the probability of a silent mutation at each index of the codons
	"""

	# List to store the probabilities for silent mutations at all index
	probabilitySilentMutations = []

	# Iterate over all the positions of the codon for the aminoacid
	for i in range(3):
		# Calculate the probability of a silent mutation in the codon
		prob = probSilentMutations(i, aminoacid, codon, trans_table, network_mutations, centrality)
		##print("\tAt index:", i, "the probability for silent mutations is:", prob)
		probabilitySilentMutations.append(prob)

	return np.array(probabilitySilentMutations)

def probSilentMutationsRegion(lenProt, transl_table, Reg_refSeqsPos, codonsInfo, network_mutations, centrality = False):
	"""
	Function that appends the probability of silent mutations for a region
	Inputs:
		lenProt: int indicating the length of the protein in region
		trans_table: dictionary containing genetic code
		Reg_refSeqsPos: array of ints containing the positions of a region 
		codonsInfo: dictionary containing the codons info of a region
		Keys: strings conatining the codon number
		Values: 
			value[0]: string containing the codon 
			value[1]: numpy array of ints containing the positions of codons with respect to the reference sequence 
			value[2]: string containing the aminoacid of the codon
		network_mutations: dictionary containing the adjacency matrices of the neutral network for an aminoacid
			Key: string indicating the aminoacid 
			Values: numpy array of integers containing the adjacency matrix of the network
		centrality: boolean indicating if eigenvector centrality is considered for probability of silent mutations calculations, for all codons
	Outputs: 
		freqMutations: array of float containing the frequency silent mutations for a region
	"""

	# Dictionary to store the probabilities of mutations for each aminoacid
	#probSMutations = {}
	# Get the probability of a silent mutation for each aminoacid
	#for aminoacid in transl_table.keys():
		#print("\tFor aminoacid:", aminoacid)
		#probSMutations[aminoacid] = probSilentMutationsAminoacid(aminoacid, transl_table, network_mutations, centrality)

	# Dictionary to store the probabilities of mutations for each aminoacid in A2_Reg1
	#probSMutations_prot_Reg = {}
	#for k, v in probSMutations.items():
		#probSMutations_prot_Reg[k] = v*(1/lenProt)

	# Calculate the frequency of mutations
	freqMutations = np.zeros((len(Reg_refSeqsPos)), dtype = float)
	for k, v in codonsInfo.items():
		codon = v[0]
		aa = v[2]
		#print("Looking at aminoacid", aa)
		positions = v[1]
		if aa != "NaN":
			print("\tFor aminoacid:", aa, "and codon:", codon)
			# Get the probability of mutation for each aminoacid
			probSMutation = probSilentMutationsAminoacid(aa, codon, transl_table, network_mutations, centrality)
			# Calculate the probability of silent mutation in the protein
			probSMutations_prot_Reg = probSMutation * (1/lenProt)
			#probMut = probSMutations_prot_Reg[aa]
			# Iterate over the positions of each codon/aminoacid
			for i in range(positions.shape[0]):
				pos = positions[i]
				freqMutations[pos] = probSMutations_prot_Reg[i]
				#print("i", i, "position", pos, "proMut", probMut[i])
				print("position", pos, "pobMut", probSMutations_prot_Reg[i])

	return freqMutations



def saveFreqMutations(fOut, freqMutations):
	"""
	Function that saves the mutations performed by random walk into a file
	INputs:
		fOut: file in which to save the mutations
		freqMutations: array of float contining the frequency of silent mutations
	Outputs:
		None
	"""

	# Open the file
	with open(fOut, "w") as fOut_:
		# Write header
		fOut_.write("#Position,FreqSMutations\n")
		# Write into file
		for i in range(len(freqMutations)):
			fOut_.write(str(i) + "," + str(freqMutations[i]) + "\n")


def loadFreqMutations(fIn):
	"""
	Function that loads the mutations performed by random walk and saved by saveMutations()
	Inputs:
		fIn: file containing the mutations
	Outputs:
		mutations: dictionary containing the times a position has been mutated
	"""

	# Dictionary to save the mutations
	mutations = {}

	# Open file
	with open(fIn, "r") as fIn_:
		# Iterate over the lines
		for line in fIn_:
			# Skip header
			if line[0] != "#":
				position = int(line.strip().split(",")[0])
				times = float(line.strip().split(",")[1])
				# Add to dictionary
				mutations[position] = times

	return mutations


def LoadFastaNodesGNID(fIn):
	"""
	Function that loads a fasta file created for the GrandNetwork
	Inputs:
		fIn: fasta file containing the sequences
	Outputs:
		d: dictionary containing the sequences in which the keys re the IDs and the values are the sequences
	"""
	# Dictionary that will store the IDs and the seqs
	d = {}
	ID = 0

	# Read the file
	with open (fIn, "r") as fIn_:
		for line in fIn_:
			# If the line is not empty
			if len(line.strip()) != 0:
				ID += 1

				# Get the sequence and add to the dictionary
				seq = line.strip()
				d[ID] = seq.upper()

	return d


def allDiffProts(seqsGN, refSeqPos, transl_table, refProt):
	"""
	Function that finds all the different proteins in the grandNetwork
	Inputs:
		seqsGN: dictionary containing all the sequences of the grandNetwork
			Key: int indicating the ID. This ID is the same as the ID in file grandNetwork.csv
			Values: strings containing a sequence
		refSeqPos: list containing the start and end positions of the codifying region
		transl_table: dictionary containing the genetic code
		refProt: string containing the sequence of the reference protein
	Outputs:
		allProts: list containing all the different proteins found in the grandNetwork
	"""

	# List to store all the different proteins found in the grandNetwork
	allProts = []

	# Append the reference protein
	allProts.append(refProt)

	# Iterate over all sequences from the grand Network
	for k, v in seqsGN.items():
		ID = k
		startPos = refSeqPos[0]
		endPos = refSeqPos[1] + 1
		seq = v[startPos:endPos]

		# Translate each sequence
		prot = Translation(seq, transl_table, stopCodon = True)
	
		# If prot is not the same as the reference protein add to the list of allProts
		if prot not in allProts:
			allProts.append(prot)

	return allProts

def saveallDiffProts(fOut, allProts):
	"""
	Function that saves all the different protein found in a network into a fasta file
	INputs:
		fOut: file in which to save the info
		allProts: list containing all the different proteins/protein fragments
	Outputs:
		None
	"""

	# Open file
	with open(fOut, "w") as fOut_:
		for ID, prot in enumerate(allProts):
			# Write an identifier
			fOut_.write(">" + str(ID + 1) + "\n")
			# Write the prot seq
			fOut_.write(str(prot) + "\n")


def allSeqsNN(allProts, seqsGN, refSeqPos, transl_table):
	"""
	Function that finds all the sequences that translates to the same protein
	for all the proteins in the grandNetwork
	Inputs:
		allProts: list of strings containing all the proteins in the grandNetwork
		seqsGN: dictionary containing all the sequences of the grandNetwork
			Key: int indicating the ID. This ID is the same as the ID in file grandNetwork.csv
			Values: strings containing a sequence
		refSeqPos: list containing the start and end positions of the codifying region
		transl_table: dictionary containing the genetic code
	Outputs: 
		seqsNN: dictionary containing all the sequences that translates to the same protein
			keys: strings containing the protein sequence
			values: list of lists
				value[0]: ID of the sequence. This ID is the same as in the file grandNetwork.csv
				value[1]: string containing the nucleotide sequence
	"""

	# Initialize the dictionary to save all the sequences that translate to the same protein
	seqsNN = {}
	for protein in allProts:
		seqsNN[protein] = []
	
	# Iterate over all sequences from the grand Network
	for k, v in seqsGN.items():
		ID = k
		startPos = refSeqPos[0]
		endPos = refSeqPos[1] + 1
		seq = v[startPos:endPos]
	
		# Translate each sequence
		prot = Translation(seq, transl_table, stopCodon = True)
	
		# Add to the dictionary
		seqsNN[prot].append([ID, v[startPos:endPos]])

	return seqsNN

def saveSeqsNN(seqsNN, dataPathOut, protName, reg):
	"""
	Functions that saves all the neutral networks (as sequences) of the grandNetwork into multiple fasta files
	Inputs:
		seqsNN: dictionary containing all the sequences that translates to the same protein
			keys: strings containing the protein sequence
			values: list of lists
				value[0]: ID of the sequence. This ID is the same as in the file grandNetwork.csv
				value[1]: string containing the nucleotide sequence
		dataPathOut: string containing the data path output
		protName: string indicating the protein name
		ref: string indicating the region name
	Outputs:
		None 
	"""

	# Set a counter that will serve as the name for the file to save
	counter = 0

	# Iterate over all neutral networks
	for k, v in seqsNN.items():
		counter += 1

		# Open output file
		fOut = os.path.join(dataPathOut + str(counter) + "_" + protName + "_" + reg + ".csv")

		# Open output file
		with open(fOut, "w") as fOut_:
			# Write some info of the protein for the current neutral network
			fOut_.write("#" + k + "\n")

			# Iterate over all the sequences in a neutral network
			for neutralSeq in v:
				ID = neutralSeq[0]
				seq = neutralSeq[1]

				# Write all the sequences for the current neutral network
				fOut_.write(">" + str(ID) + "\n")
				fOut_.write(str(seq) + "\n")



def numSeqsPerFile(seqsNN):
	"""
	Function that calculates the number of sequence in a Neutral Network file
	Inputs:
		seqsNN: dictionary containing all the sequences that translates to the same protein
			keys: strings containing the protein sequence
			values: list of lists
				value[0]: ID of the sequence. This ID is the same as in the file grandNetwork.csv
				value[1]: string containing the nucleotide sequence
	Outputs: 
		dTimesNN:  dictionary containing the number of sequences per file. I say file because the key
				of the dicitonary if also the name of the file
	"""

	# Dictionary to save the number of sequences for each file
	dTimesNN = {}
	# Set a counter. This counter is the name of file saved by function saveSeqsNN()
	counter = 0

	# Iterate over all entries of seqsNN
	for k, v in seqsNN.items():
		counter += 1
		dTimesNN[counter] = len(v)

	return dTimesNN

def sortDictSaveInfo(fOut, dTimesNN):
	"""
	Function that sorts a dictionary by its value and saves the information into a file
	Inputs:
		Fout: name of the file to save the metaInfo
		dTimesNN: dTimesNN:  dictionary containing the number of sequences per file. I say file because the key
				of the dicitonary if also the name of the file
	Outputs: None 
	"""

	# Sort the dictionary
	dTimesNN_sorted = sorted(dTimesNN.items(), key=lambda x: x[1], reverse=True)

	# Open file
	with open(fOut, "w") as fOut_:
		# Write a header
		fOut_.write("#FileName,NumberSeqs\n")
		# Iterate over the dictionary
		for NN in dTimesNN_sorted:
			fOut_.write(str(NN[0]) + "," + str(NN[1]) + "\n")

def loadSeqsPerFile(fIn):
	"""
	Funtion that loads the number of sequences per file
	Inputs:
		fIn: file to load the data
	Output: 
		dTimesNN: dictionary containing the number of sequences per file. I say file because the key
				of the dicitonary is also the name of the file
	"""

	# Declare structures
	dTimesNN = {}

	# Iterate over file
	with open(fIn, "r") as fIn_:
		for line in fIn_:
			# Skip header
			if line[0]!= "#":
				file = int(line.strip().split(",")[0])
				numSeqs = int(line.strip().split(",")[1])

				# Add to dictionary
				dTimesNN[file] = numSeqs

	return dTimesNN


def frequencyMatrix(matrixSeqs, pseudocount = False):
	"""
	Function that calculates the frequency matrix for each nucleotide at each position.
	BEWARE: this function only works if all the sequences in the matrix has the same length!
	Inputs:
		matrixSeqs: numpy matrix of strings that contains all the sequences in a FASTA file
		pseudocount: 
	Outputs:
		freqMatrix: numpy matrix of floats that contains the frequency of each nucleotide at each position
	"""

	# Get the number of positions for the sequence
	numPos = len(matrixSeqs[0])

	# Create the atrix to store the counts of each nucleotide at each position
	matrix_counts = np.zeros((4, numPos), dtype = int)

	# Iterate over the positions of the sequences
	for pos in range(numPos):
		# Dictionary to save the counts of each nucleotide at pos
		if pseudocount:
			counts = {"A": 1, "C": 1, "G": 1, "T": 1}
		else:
			counts = {"A": 0, "C": 0, "G": 0, "T": 0}
		for seqNum in range(matrixSeqs.shape[0]):
			# Update counts
			nucl = matrixSeqs[seqNum][pos]
			counts[nucl] += 1

		# Add the counts to the pos (i.e. index) of matrix
		matrix_counts[: ,pos] = np.array(list(counts.values()))

	# Get the total of nucleotides at a position (this is the same as the number of sequences of the matrix)
	totalNucl = np.sum(matrix_counts[: ,1])

	# Calculate the frequency matrix
	matrix_freq = np.zeros((4, numPos), dtype = float)

	for j in range(matrix_freq.shape[1]):
		matrix_freq[: ,j] = matrix_counts[: ,j]/totalNucl

	return matrix_freq


def saveEmpiricalFreqMutations(fOut, positions, freqMuts):
	"""
	Function that saves the empirical frequency of mutations into a file
	Inputs:
		fOut: name of the file
		positions: array of ints containing the positions along a sequence
		freqMuts: array of floats containing the frequency of mutation of a position i
	"""

	# Open file
	with open(fOut, "w") as fOut_:

		# Write header
		fOut_.write("#Position,FreqEmpMutations\n")

		# Iterate over the data
		for i in range(freqMuts.shape[0]):
			fOut_.write(str(positions[i]) + "," + str(freqMuts[i]) + "\n")




















