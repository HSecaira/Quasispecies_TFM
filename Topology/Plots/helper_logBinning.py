"""
Script that contains functions to perform logarithmic binning
"""

# Imports
import numpy as np;
import os, sys;


def histogramDegrees(degrees, minDegreeInterval = 0, maxDegreeInterval = 10, base = 2.0):
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

	# Normalize again by the number of nodes in the network, so we have a probability distribution
	histogram_norm = histogram_temp/degrees.shape[0]

	return histogram_norm, binsUnique, widths

def histogramClustCoeffKnn(parameterDict, bins, widths):
	"""
	Function that performs a logarithmic binning over the nodes of a network by
	building a histogram and normalizing each bin by its width. Beware! This functions 
	only works for C(k) and Knn(k).
	Inputs:	
		parameterDict: dictionary whose keys are degrees and values are lists that contains 
		nodes with a parameter associated to a degree.
		bins: array of ints containing the intervals for all bins. The intervals are discrete.
		widts: array of ints: containing the width of each interval. 
	Outputs: 
		histogram_norm: array of floats containing all the data points that falls on a given bin. The array has been
		normalized by dividing each bin by its width.
	"""

	# Declare some variables
	histogram = []

	print('\ncreating histogram!')

	# Iterate over the bins
	for i in range(bins.shape[0] - 1):
		# The intervals
		minDegree = bins[i]
		# -1 since the right extreme of the interval is not taken into account and degrees are ints
		maxDegree = bins[i + 1] - 1

		print("minDegree", minDegree, "maxDegree", maxDegree)

		# Count all nodes that fall in the interval
		if minDegree == maxDegree:
			# Sum the parameter of all the point in the current interval
			suma = 0 
			for paramPoint in parameterDict[minDegree]:
				suma += paramPoint
			# Count the number of nodes at the current interval
			n = len(parameterDict[minDegree])
			histogram.append(suma/n)

		else:
			count = 0
			suma = 0
			# In case the interval comprises more than one degree
			print("interval comprises more than one degree\n")
			for j in range(minDegree, maxDegree + 1):
				#  Check if the degree exists
				if j in parameterDict.keys():
					print("\tCheck the degree in interval", i)
					# Count the number of nodes at the current interval
					count += len(parameterDict[j])
					#print("counttt", count, "degree", i)
					# Sum the parameter of all the point in the current interval
					for paramPoint in parameterDict[j]:
						suma += paramPoint

			print("\tcount outside", count)
			print("\tsum outside", suma)
			# In case no points fall in the interval
			if count == 0:
				print("\tno points fall in the interval!!!")
				histogram.append(0)
			else:
				histogram.append(suma/count)

		print("iii", i)

	# Normalize histogram
	histogram_norm = np.array(histogram)/widths

	print('shapes', len(histogram), len(widths))

	return histogram_norm

def histogramClustCoeffKnnAlternative(parameterDict, bins, widths):
	"""
	Function that performs a logarithmic binning over the nodes of a network by
	building a histogram and normalizing each bin by its width. Beware! This functions 
	only works for C(k) and Knn(k).
	Inputs:	
		parameterDict: dictionary whose keys are degrees and values are C(k) or Knn(k)
		bins: array of ints containing the intervals for all bins. The intervals are discrete.
		widts: array of ints: containing the width of each interval. 
	Outputs: 
		histogram_norm: array of floats containing all the data points that falls on a given bin. The array has been
		normalized by dividing each bin by its width.
	"""

	# Declare some variables
	#print('alternative histogram')
	histogram = []

	# Iterate over the bins
	for i in range(bins.shape[0] - 1):
		# The intervals
		minDegree = bins[i]
		# -1 since the right extreme of the interval is not taken into account and degrees are ints
		maxDegree = bins[i + 1] - 1

		#print("minDegree", minDegree, "maxDegree", maxDegree)

		# Count all nodes that fall in the interval
		if minDegree == maxDegree:
			# Sum the parameter of all the point in the current interval
			suma = parameterDict[minDegree] 
			#print("\tsumaa", suma)
			#for paramPoint in parameterDict[minDegree]:
			#	suma += paramPoint
			# Count the number of nodes at the current interval
			#n = len(parameterDict[minDegree])
			histogram.append(suma)

		else:
			count = 0
			suma = 0
			# In case the interval comprises more than one degree
			#print("interval comprises more than one degree\n")
			for j in range(minDegree, maxDegree + 1):
				#  Check if the degree exists
				if j in parameterDict.keys():
					#print("\tCheck the degree in interval", i)
					# Count the number of nodes at the current interval
					#count += len(parameterDict[j])
					#print("counttt", count, "degree", i)
					# Sum the parameter of all the point in the current interval
					suma += parameterDict[j]
					#for paramPoint in parameterDict[j]:
					#	suma += paramPoint


			# In case no points fall in the interval
			#print("\tsumaa", suma)
			if suma == 0:
				histogram.append(0)
			else:
				histogram.append(suma)

	# Normalize histogram
	histogram_norm = np.array(histogram)/widths

	return histogram_norm




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

def alphaNewmanClauset(degrees, kmin):
	"""
	Function that calculates the exponent coefficient (alpha) for the degree distribution of a network
	Inputs: 
		degrees: array of ints containing the degrees of a network
		kmin: int indicating the minimum degree for which the power law holds
	Outputs:
		alpha: float containing the exponent of the power law
			alpha = 1 + N * (sum_{i}^{N}log((k_i)/(k_min - 0.5)))^-1
		sigma: float containing the statistical error of the estimation of alpha 
			sigma = (alpha - 1)/(srqt(N))
	"""

	# Get the number of nodes with a degree greater than or equal to kmin
	N = degrees[degrees >= kmin].shape[0]
	k_min = degrees[degrees >= kmin]

	# Perform the sum
	summ = np.sum(np.log((k_min)/(kmin - 0.5)))

	# Get alpha
	alpha = 1 + N * (1/summ)
	# Get sigma
	sigma = (alpha - 1)/np.sqrt(N)

	return alpha, sigma
