// Imports: 
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib> 
#include <vector>
#include <random>
#include <algorithm>
#include <map>
#include <set>
#include <bits/stdc++.h> 
#include <boost/algorithm/string.hpp> 
#include <typeinfo>
#include <cmath>
#include <iterator>
#include <tuple>
#include <numeric>
#include <stdexcept>
#include <functional> // to divide the whole array
#include <chrono>
#include <stack> 
#include<unordered_map>
#include "gNets.h"

using namespace std; 

namespace f{


	// Functions //

	vector<string> SplitString(string line, string delim){
		// Function that splits a line given a delimiter
		// Inputs:
			//line: the line to split
			// delim: the delimiter what will be used to split the line
		// Outputs
			// lineSPlit: Vector of strings that contains the split line.
				// First element: everythong before the delimiter
				// Second element: everything after the delimiter
		// Taken from https://www.geeksforgeeks.org/boostsplit-c-library/#:~:text=This%20function%20is%20similar%20to,by%20means%20of%20the%20predicate.&text=Application%20%3A%20It%20is%20used%20to,which%20are%20separated%20by%20separators.
		vector<string> lineSplit;
		boost::split(lineSplit, line, boost::is_any_of(delim));
		return lineSplit;
	}

	vector<int> SplitStringInts(const string &line, char delim){
		// Function that splits a line given a delimiter
		// Inputs:
			//line: the line to split
			// delim: the delimiter what will be used to split the line
		// Outputs
			// lineSplit: Vector of ints that contains the split line.
				// First element: everythong before the delimiter
				// Second element: everything after the delimiter

		// Declare some variables
		vector<int> lineSplit;	
		string entry;
		istringstream entriesStream(line);

		while (std::getline(entriesStream, entry, delim)) {      
	        lineSplit.push_back(stoi(entry));     
	    }
		return lineSplit;
	}

	bool sortbyVal(const pair<string,int> &a, const pair<string, int> &b){
		// Taken from: https://www.geeksforgeeks.org/sorting-vector-of-pairs-in-c-set-1-sort-by-first-and-second/
		// Comparator function to pass to the sort() module. Sorts values in ascending order
		// Inputs: 
			// "&" uses a pointer to the values so NO inputs are required
		// Outputs
			// A boolean that evaluates an expression
		return (a.second < b.second);
	}

	bool sortbyValDesc(const pair<int,int> &a, const pair<int, int> &b){
		// Taken from: https://www.geeksforgeeks.org/sorting-vector-of-pairs-in-c-set-1-sort-by-first-and-second/
		// Comparator function to pass to the sort() module. Sorts values in descending order
		// Inputs: 
			// "&" uses a pointer to the values so NO inputs are required
		// Outputs
			// A boolean that evaluates an expression
		return (a.second > b.second);
	}

	bool sortbyValDescDouble(const pair<int,double> &a, const pair<int, double> &b){
		// Taken from: https://www.geeksforgeeks.org/sorting-vector-of-pairs-in-c-set-1-sort-by-first-and-second/
		// Comparator function to pass to the sort() module. Sorts values in descending order
		// Inputs: 
			// "&" uses a pointer to the values so NO inputs are required
		// Outputs
			// A boolean that evaluates an expression
		return (a.second > b.second);
	}


	bool maxValueDict(const pair<int, int>& a, const pair<int, int>& b){
		// Comparator function to pass to the max_element module to find the maximum value
		// of a dictionary values(). 
		// Inputs:
			// "&" uses a pointer to the values, so no inputs are required
		// Outputs:
			// A boolean that evaluates the expression
		return (a.second < b.second);
	}

	bool maxValueDictDouble(const pair<int, double>& a, const pair<int, double>& b){
		// Comparator function to pass to the max_element module to find the maximum value
		// of a dictionary values(). 
		// Inputs:
			// "&" uses a pointer to the values, so no inputs are required
		// Outputs:
			// A boolean that evaluates the expression
		return (a.second < b.second);
	}

	double vectorNorm(vector<double> const & v){
		// Function that calculates the euclidean norm of a vector
		// Inputs:
			// v: vector which uses a reference to the original vector without modifying the original vector
		// OutputsL
			// norm: norm of the vector

		// Declare some variables
		double norm, sum = 0;

		// Iterate over the elements of the vector
		for (double elem : v){
			sum += elem * elem;
		}

		norm = sqrt(sum);

		return norm;
	}

	double vectorNormProb(vector<double> const & v){
		// Function that normalizes a vector so it is a probability distribution
		// Inputs:
			// v: vector which uses a reference to the original vector without modifying the original vector
		// OutputsL
			// norm: norm of the vector

		// Declare some variables
		double norm, sum = 0;

		// Iterate over the elements of the vector
		for (double elem : v){
			sum += elem;
		}

		norm = sum;

		return norm;
	}

	vector<int> randomSampleNodes(int numNodes, vector<int> allNodesBComponent){
		/*
		Function that randomly sample a given number of nodes without repetition from the adjacency list
		Inputs:
			numNodes: int indicating the number of nodes to sample
			DictAdjList: dictionary containing the network as an adjacency list
		Outputs:
			sampledNodes: vector of ints containing the sampled nodes
		*/

		// Declare some variables
		//vector<int> nodes;
		vector<int> sampledNodes;

		// Get all the nodes from the network
		//for (auto &item : DictAdjList){
			//int node = item.first;
			//nodes.push_back(node);
		//}


		//  Do sampling withot repetition
    	size_t nelems = numNodes;
    	sample(allNodesBComponent.begin(), allNodesBComponent.end(), back_inserter(sampledNodes), nelems, mt19937{random_device{}()});

    	return sampledNodes;
	}


	map<int, vector<int>> LoadGNetDict(string fInName){
		// Function that loads a network that has been stored in the format produce by saveGNet()
		// from the GNet class from Luiño's code and save it as a dictionary which resembles and adjacency list
		// Inputs:
			// fInName: file containing the network produced by saveGNet()
		// Outputs
			// DictAdjListGNet: dictionary in which keys are the nodes (string), and values are vectors (string)

		// Declare some variables
		gN::GNet gNet; // Luiño's Network class
		map<int, vector<int>> DictAdjListGNet;

		// Load the network
		gNet.loadGNet(fInName.c_str());

		// Build the dictionary
		for (int i=0; i<gNet.getNGenotypes(); i++){
			int node = i;

			// Ignore nodes that have no connections
			if (gNet.getGConnections(i).size() != 0){
				DictAdjListGNet[node] = gNet.getGConnections(i);
			}
		}

		return DictAdjListGNet;
	}

	unordered_map<int, vector<int>> LoadGNetUnorderedDict(string fInName){
		// Function that loads a network that has been stored in the format produce by saveGNet()
		// from the GNet class from Luiño's code and save it as an unordered_map which resembles and adjacency list
		// implemented in C++ as a Hash table
		// Inputs:
			// fInName: file containing the network produced by saveGNet()
		// Outputs
			// DictAdjListGNet: dictionary in which keys are the nodes (string), and values are vectors (string)

		// Declare some variables
		gN::GNet gNet; // Luiño's Network class
		unordered_map<int, vector<int>> DictAdjListGNet;

		// Load the network
		gNet.loadGNet(fInName.c_str());

		// Build the dictionary
		for (int i=0; i<gNet.getNGenotypes(); i++){
			int node = i;

			// Ignore nodes that have no connections
			if (gNet.getGConnections(i).size() != 0){
				DictAdjListGNet[node] = gNet.getGConnections(i);
			}
		}

		return DictAdjListGNet;
	}


	map <int, vector<int>> ConvertDict(string fInName){
		// Function that creates the adjacency list as a dictionary
		// Inputs: 
			// fInName: file containing the adjacency list
		//Outputs:
			// DictAdjList: dictionary in which keys are the nodes (string), and values are vectors (string)

		// Declare the input file stream and other variables
		ifstream fIn(fInName.c_str()); 
		map <int, vector<int>> DictAdjList;
		string line;
		int node;
		vector<string> lineSplit;
		vector<int> degreeSplit;
		char delim = ',';

		// Read a the file line by line
		while (getline(fIn, line)){
			// Omit the header
			if (line[0] != '#'){
				// Split the line
				lineSplit = SplitString(line, ":");
				node = stoi(lineSplit[0]);
				// Retrieve all the connections
				degreeSplit = SplitStringInts(lineSplit[1], delim);

				// Save in a dictionary
				DictAdjList[node] = degreeSplit;
			}
		}

		// Close file
		fIn.close();

		return DictAdjList;
	}

	unordered_map <int, vector<int>> ConvertUnorderedDict(string fInName){
		// Function that creates the adjacency list as a dictionary
		// Inputs: 
			// fInName: file containing the adjacency list
		//Outputs:
			// DictAdjList: dictionary in which keys are the nodes (string), and values are vectors (string)

		// Declare the input file stream and other variables
		ifstream fIn(fInName.c_str()); 
		unordered_map <int, vector<int>> DictAdjList;
		string line;
		int node;
		vector<string> lineSplit;
		vector<int> degreeSplit;
		char delim = ',';

		// Read a the file line by line
		while (getline(fIn, line)){
			// Omit the header
			if (line[0] != '#'){
				// Split the line
				lineSplit = SplitString(line, ":");
				node = stoi(lineSplit[0]);
				// Retrieve all the connections
				degreeSplit = SplitStringInts(lineSplit[1], delim);

				// Save in a dictionary
				DictAdjList[node] = degreeSplit;
			}
		}

		// Close file
		fIn.close();

		return DictAdjList;
	}

	void saveNodesDegrees(string fOutName, map<int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// FUnction that saves the nodes and its degree, rank, and Pk, where Pk = rank/n, n = number of nodes
		// Useful to calculate the cumulative distribution of the network
		// INputs:
			// fOutName: file name to save the key:value pairs
			// DictAdjList: dictionary containing the adjacency list of the network
			// allNodesBComponent: vector of ints containing all the nodes of all components of the network

		// Open file and declare some variables: 
		ofstream fOut(fOutName.c_str());
		vector<pair<int, int>> vectorTempDegrees;
		double n = allNodesBComponent.size();

		// Get a dictionary of degrees
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];
			int degree = DictAdjList[node].size();
			vectorTempDegrees.push_back(make_pair(node, degree));
		} 

	    // Sort degrees
		sort(vectorTempDegrees.begin(), vectorTempDegrees.end(), sortbyValDesc);

		fOut << "#NodeID,Degree,Rank,Pk" << endl;

		// Iterate over the nodes of the biggest component
		for (int i = 0; i < vectorTempDegrees.size(); i++){
			int node = vectorTempDegrees[i].first;
			int degree = vectorTempDegrees[i].second;
			double Pk = (i + 1)/n;
			// Write to file
			fOut << node << "," << degree << "," << i + 1 << "," << Pk << endl;

		}
	}


	tuple<double, double> AverageDegree(map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that calculates the average degree of a network and its standard deviation
		// Inputs: 
			// DictAdjList: dictionary containing the adjacency list
			// allNodesComponents: vector of ints containing all the nodes of all components of the network
		// Outputs:
			// averageDegree: average degree of a network

		// Declare some variables
		double averageDegree, std, sumDegrees = 0, sumSTD = 0;

		// Iterate over all the nodes in the biggest component
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];
			double degreeNode = DictAdjList[node].size();
			sumDegrees += degreeNode;
		}

		double n = allNodesBComponent.size();
		//averageDegree = 1/entries * sumDegrees;
		averageDegree = sumDegrees/n;

		// Standard deviation
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];
			double degreeNode = DictAdjList[node].size();
			sumSTD += pow((degreeNode - averageDegree), 2);
		}


		std = sqrt((sumSTD)/(n));

		return {averageDegree, std};
	}

	double Assortativity(map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// FUnction that calculates the assorativity/correlation coefficient of a networks from an
		// adjacency list. The coefficient is calculated using the formulas 8.27, 8.28, and 8.29 from
		// Newman - Networks An Introduction (2010) book.
		// In the future, this function will be modified so only unique pairs of vertices are taken!!
		// Inputs:
			// DictAdjList: dictionary containting all the nodes the beighbors of each node
			// allNodesComponents: vector of ints containing all the nodes of all components of the network
		// Outputs:
			// Assortativity/correlation coefficient	

		// Declare variables
		double Se, S1, S2, S3, NodeDegree, tempDegree;
		Se = 0;
		S1 = 0;
		S2 = 0; 
		S3 = 0;

		// Iterate over the nodes of the biggest component
		for (int i = 0; i < allNodesBComponent.size(); i++){
			// Get the degree of the current node
			int node = allNodesBComponent[i];
			NodeDegree = DictAdjList[node].size();

			// Calculate S1, S2, and S3
			S1 += NodeDegree;
			S2 += pow(NodeDegree, 2);
			S3 += pow(NodeDegree, 3);

			// Iterate over the neighbors of the current node
			vector<int> neighbors = DictAdjList[node];
			for (int j = 0; j < neighbors.size(); j++){
				// Get the degree of a neighbor of the current node
				tempDegree = DictAdjList[neighbors[j]].size();

				// Get the sum for Se
				Se += NodeDegree * tempDegree;
			}
		}

		// Calculate th assortativity/correlation coefficient
		double r = ((S1*Se) - pow(S2, 2))/((S1*S3) - pow(S2, 2));
		return r;
	}

	double AverageNeighborDegreeSingleNode(map <int, vector<int>>& DictAdjList, int node){
		// Function that calculates the Average Neighbor degree of a node
		// Inputs:
			// DictAdjList: Dictionary containing the adjacency list of the network
			// node: node for which you want to calculate the average neighbor degree
			// allNodesComponents: vector of vector of ints containing all the nodes of all components of the network
			// indexBComponent: int containing the index of the biggest component of th network
		// Outputs:
			// Knn: average neighbor degree of node 

		// Declare some variables
		double Knn, degreeNeighbors = 0;

		// Get the neighbors of node
		vector<int> neighbors = DictAdjList[node]; 
		double numNeighbors = neighbors.size();

		// Iterate over the neighbors
		for (int i = 0; i < neighbors.size(); i++){
			int neighbor = neighbors[i];
			double degreeNeighbor = DictAdjList[neighbor].size();

			// Get the degree of the neighbor
			degreeNeighbors  += degreeNeighbor;
		}

		// Calculate the average neighbor degree
		Knn = 1/numNeighbors * degreeNeighbors;

		return Knn;
	}


	map<int, double> AverageNeighborDegree(map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that calculates the Average Neighbor Degree for all nodes in a network
		// Inputs:
			// DictAdjList: dictionary containing the adjacency list of the network
			// allNodesComponents: vector of ints containing all the nodes of all components of the network
		// Outputs
			// averageNeighborDegree: dictionary containing the average neighbor degree of all nodes in the network

		// Declare some variables
		map<int, double> averageNeighborDegree;

		// Iterate over all nodes in the biggest component
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];

			// Calculate the average degree node for a single node
			double Knn = AverageNeighborDegreeSingleNode(DictAdjList, node);

			// Add to dictionary
			averageNeighborDegree[node] = Knn;
		}

		return averageNeighborDegree;
	}

	void saveAverageNeighborDegree(string fOutName, map<int, vector<int>>& DictAdjList, map<int, double>& averageNeighborDegree){
		// Function that saves the average neighbor degree of all nodes into a file
		// INputs: 
			// fInName: name of the file to save
			// averageNeighborDegree: dictionary containing the average neighbor degreee of all nodes in thte network
		// Outputs: 
			// None

		// Declare some variables and files
		ofstream fOut(fOutName.c_str());
		double n = averageNeighborDegree.size();
		vector<pair<int, double>> vectorTempKnn;

		// Creater pair vector 
		for (auto &item:averageNeighborDegree){
			int node = item.first;
			double Knn = item.second;
			vectorTempKnn.push_back(make_pair(node, Knn));
		}

		// Sort vector by KClust coeff in descending order
		sort(vectorTempKnn.begin(), vectorTempKnn.end(), sortbyValDescDouble);

		// Header for file
		fOut << "#NodeID,AverageNeighborDegree,Degree,Rank,PKnn" << endl;

		for (int i = 0; i < vectorTempKnn.size(); i++){
			int node = vectorTempKnn[i].first;
			int degree = DictAdjList[node].size();
			double Knn = vectorTempKnn[i].second;
			double PKnn = (i + 1)/n;
			// Write to file
			fOut << node << "," << Knn << "," << degree << "," << i + 1 << "," << PKnn << endl;
		}
	}

	map<int, double> CountDegrees(map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that counts the times a degree node appears in the network
		// Inputs:
			// DictAdjList: dictionary containing the adjacency list of the network
			// allNodesBComponents: vector of ints containing all the nodes of all components of the network
		// Outputs:
			// degreeCounts: dictionary containing the times a degree appears in the network

		// Declare some variables
		map<int, double> degreeCounts;

		// Iterate over the nodes in the biggest component
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];
			int degree = DictAdjList[node].size();

			// If the degree is not in the keys
			if (degreeCounts.find(degree) == degreeCounts.end()){
				degreeCounts[degree] = 1;
			}
			else{
				degreeCounts[degree] += 1;
			}
		}

		return degreeCounts;
	}

	double Density(map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that calculates the density or connectance of a graph
		// The density lies between 0 and 1 (included). The the density os closer to zero, the network is sparse
		// If the density is closer to one, the network is dense (i.e. it contains the maximum possible number of edges)
		// INputs:
			// DictAdjList: Dictionary containing the adjacency list o f the network
			// allNodesComponent: vector of ints containing all the nodes of all components of the network
		// Outputs
			//density: density of the network

		// Declare some variables 
		double density, m = 0;

		// Get the size of the network (i.e. number of edges)
		map<int, double> degreeCounts = CountDegrees(DictAdjList, allNodesBComponent);
		for (auto &item : degreeCounts){
			m += item.first * item.second;
		}	

		// Get the order of the biggest component
		double n = allNodesBComponent.size();

		// Calculate the density. I do not multiply by two because m already contains all the edges of the network
		density = m/(n * (n - 1));

		return density;
	}

	map<int, double> KNearestNeighbors(map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that calculates the k nearest negihbors of a network
		// This function differs from AverageNeighborDegree() because kNearestNeighbors() calculates
		// the k nearest neighbors for every degree, which is the same as calculating the average 
		// negihbor degree but KNearestNeighbors() averages the values for every degree.
		// INputs: 
			// DictAdjList: dictionary containing the adjacency list of the network
			// allNodesBComponent: vector of ints containing all the nodes of all components of the network
		// Outputs:
			// k_nearestNeighbors: dictionary containing the k nearest neighbors for every degree in the network

		// Declare some variables
		map<int, double> knearestneighborsTemp, k_nearestNeighbors;

		// Iterate over all nodes to obtain the sum of the Average Neighbor degree for each degree
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];
			int degree = DictAdjList[node].size();

			// Calculate the average degree node for a single node
			double Knn = AverageNeighborDegreeSingleNode(DictAdjList, node);

			// Add to dictionary
			// If the degree not in keys
			if (knearestneighborsTemp.find(degree) == knearestneighborsTemp.end()){
				knearestneighborsTemp[degree] = Knn;
			}
			else{
				knearestneighborsTemp[degree] += Knn;
			}
		}

		// Get the times a degree appears in the network 
		map<int, double> degreeCounts = CountDegrees(DictAdjList, allNodesBComponent);

		// Iterate over the degrees
		for (auto &item : knearestneighborsTemp){
			int degree_ = item.first;
			// Get the times the degree appear
			double times = degreeCounts[degree_];
			double Knn_ = item.second/times;

			// Append to the dictionary
			k_nearestNeighbors[degree_] = Knn_;
		}

		return k_nearestNeighbors; 
	}



	void saveKNearestNeighbors(string fOutName, map<int, double>& k_nearestNeighbors){
		// Function that saves the K nearest neighbor degree of all nodes into a file
		// INputs: 
			// fOutName: name of the file to save
			// k_nearestNeighbors: dictionary containing the K nearest neighbor degree of all nodes in thte network
		// Outputs: 
			// None

		// Declare some variables and files
		ofstream fOut(fOutName.c_str());
		vector<pair<int, double>> vectorTempKnn;
		double k = k_nearestNeighbors.size();

		// Creater pair vector 
		for (auto &item:k_nearestNeighbors){
			int degree = item.first;
			double Knn = item.second;
			vectorTempKnn.push_back(make_pair(degree, Knn));
		}

		// Sort vector by KClust coeff in descending order
		sort(vectorTempKnn.begin(), vectorTempKnn.end(), sortbyValDescDouble);

		// Header for file
		fOut << "#Degree,KNNDegree,Rank,PKnn" << endl;

		for (int i = 0; i < vectorTempKnn.size(); i++){
			int degree = vectorTempKnn[i].first;
			double Knn = vectorTempKnn[i].second;
			double PKnn = (i + 1)/k;
			// Write to file
			fOut << degree << "," << Knn << "," << i + 1 << "," << PKnn << endl;
		}
	}

	void savekNearestNeighborsDegree(string fOutName, map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that calculates the k nearest neighbors of a network and saves them into a file
		// This function differs from AverageNeighborDegree() and KNearestNeighbors() because it calculates
		// the k nearest neighbor for every degree without taking the average of the values for every degree
		// INputs: 
			// DictAdjList: dictionary containing the adjacency list of the network
			// allNodesBComponent: vector of ints containing all the nodes of all components of the network
			// fOutName: name of the file to save

		// Declare  some variables
		ofstream fOut(fOutName.c_str());
		map<int, vector<double>> knearestneighborsTemp;

		// Header for file
		fOut << "#Degree,KNearestNeighbor" << endl;


		// Iterate over all the nodes of the network
		for (int  i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];
			int degree = DictAdjList[node].size();

			// Calculate the average degree node for a single node
			double Knn = AverageNeighborDegreeSingleNode(DictAdjList, node);

			// Add to dictionary
			// If the degree not in keys
			knearestneighborsTemp[degree].push_back(Knn);
		}

		// Write to file
		for (auto &item : knearestneighborsTemp){
			int degree_ = item.first;
			vector<double> KNN_ = knearestneighborsTemp[degree_];
			fOut << degree_ << ":";
			for (int i = 0; i < KNN_.size(); i++){
				if (i != KNN_.size() - 1){
					fOut << KNN_[i] << ",";
				}
				else{
					fOut << KNN_[i];
				}
			}
			fOut << endl;
		}
	}



	double LocalClustCoeff_Numerator(map <int, vector<int>> &DictAdjList, int node){
		// Function that calculates the numerator of the local clustering coefficient according to equation
		// 8.25/10.2 of Newman's Book: Networks and Introduction, 2010
		// INputs:
			// DIctAdjList: Dictionary containing the adjacency list
			// node: Node for which you want to calculate the clustering coefficient
		// Outputs:
			// numerator: NUmerator of the local clustering coefficient

		//Declare variables
		vector<int> neighbors;
		double numerator = 0;

		// Get the neighbors of node
		neighbors = DictAdjList[node];
		// If there is only one or zero neighbors
		if (neighbors.size() == 1 || neighbors.size() == 0){
			return 0;
		}
		else{
			// Iterate over the neighbors of node
			for (int i = 0; i < neighbors.size(); i++){
				int CurrentNode = neighbors[i];
				vector<int> neighbors_CurrentNode = DictAdjList[CurrentNode];
				// If there is only one or zero neighbors
				if (neighbors_CurrentNode.size() == 1 || neighbors_CurrentNode.size() == 0){
					numerator += 0; // A bit dumb
				}
				else{
					// Iterate over the neighbors of CurrentNode
					for (int j = 0; j < neighbors_CurrentNode.size(); j++){
						// Create a copy of neighbors of node that does not including the CurrentNode
						vector <int> temp_neighbors = vector<int>(neighbors.begin() + i + 1, neighbors.end());
						// Check if one of the temp neighbors is a neighbor of CurrentNode
						if (find(temp_neighbors.begin(), temp_neighbors.end(), neighbors_CurrentNode[j]) != temp_neighbors.end()){
							numerator += 1;
						}
					}
				}
							
			}

		}
		
		return numerator;
	}

	double LocalClustCoeff_Denominator(map <int, vector<int>>& DictAdjList, int node){
		// Function that caculates the numerator of the local clustering coefficient of a node according to equation
		// 8.25/10.2 of Newman's Book: Networks and Introduction, 2010 
		// INputs:
			// DictAdjList: Dictionary containing the adjacency list
			// node: Node for which you want to calculate the clustering coefficient
		// Outputs:
			// denominator: Denomiantor of the local clustering coefficient

		// Declare variables
		double denominator;
		int degreeNode;
		vector<int> neighbors;

		// Get the neighbors of node
		neighbors = DictAdjList[node];
		// Get the degree of node
		degreeNode = neighbors.size();

		// Calculate the denominator, this is the number of pairs of neighbors of node
		denominator = (0.5 * degreeNode) * (degreeNode - 1);

		return denominator;

	}

	double GlobalClusteringCoeff(map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that calculates the global clustering coefficient given an adjacency list as a dictionary
		// according to equation 7.41/10.3 of Newman's Book: Networks and Introduction, 2010 
		// INputs:
			// DictAdjList: Dictionary containing the adjacency list
			// allNodesBComponent: vector of ints containing all the nodes of all components of the network
		// Outputs
			// GlobalClustCoeff: Global Clustering coefficient

		// Declare variables
		double GlobalClustCoeff, numerator = 0, denominator = 0;
		int node;


		// Iterate over all nodes
		for (int i = 0; i < allNodesBComponent.size(); i++){
			node = allNodesBComponent[i];
			numerator += LocalClustCoeff_Numerator(DictAdjList, node);
			denominator += LocalClustCoeff_Denominator(DictAdjList, node);
		}

		GlobalClustCoeff = numerator/denominator;

		return GlobalClustCoeff;
	}

	double AverageClusteringCoeff(map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that calculates the average Clustering coefficient of a network given an adjacency list
		// Inputs:
			// DictAdjList: Dictionary containing the adjacency lists
			// allNodesBComponent: vector of ints containing all the nodes of all components of the network
		// Outputs:
			// AverageClusteringCoefficient: Average clustering coefficient of a network

		// Declare variables
		double LocalClustering, sumLocalClustering = 0, numerator, denominator, AverageClusteringCoefficient;
		int node;

		// Iterate over all the nodes of the Adjacency list
		for (int i = 0; i < allNodesBComponent.size(); i++){
			node = allNodesBComponent[i];
			// Calculate the local clustering coefficient
			numerator = LocalClustCoeff_Numerator(DictAdjList, node);
			denominator = LocalClustCoeff_Denominator(DictAdjList, node);
			// Beware! the denominator can return 0!!
			if (denominator != 0){
				LocalClustering = numerator/denominator;
				sumLocalClustering += LocalClustering;
			}
			else {
				sumLocalClustering += 0;
			}
		}

		double entries = DictAdjList.size();
		AverageClusteringCoefficient = 1/entries * sumLocalClustering;

		return AverageClusteringCoefficient;
	}

	map<int, double> LocalClustCoeff(map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that calculates the Local Clustering Coefficient of all nodes
		// Inputs:
			// DictAdjList: dictionary containing the adjacenchy list of the network
			// allNodesBComponent: vector of ints containing all the nodes of all components of the network
		// Outputs: 
			// clustering: dictionary containing the local clustering coefficient for all nodes

		// Declare some variables
		map<int, double> clustering; 
		double LocalClustering;

		// Iterate over all the nodes in the biggest component
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];

			// Calculate the local clustering coefficient
			double numerator = LocalClustCoeff_Numerator(DictAdjList, node);
			double denominator = LocalClustCoeff_Denominator(DictAdjList, node);

			// Beware! the denominator can return 0!!
			if (denominator != 0){
				LocalClustering = numerator/denominator;
			}
			else {
				LocalClustering = 0.0;
			}

			// Add to the dictionary
			clustering[node] = LocalClustering;
		}

		return clustering; 
	}

	void saveLocalClustCoeff(string fOutName, map<int, double>& clustering, map<int, vector<int>>& DictAdjList){
		// Function that saves the clustering coefficient of every node into a file
		// Inputs:
			// fOutName: file name in which the clustering will be saved
			// clustering: dictionary containing the clustering coefficient of all nodes
			// DictAdjList: dictionary containing the adjacency list of the network
		// Outputs: 
			// None

		// Declare some variables and files
		ofstream fOut(fOutName.c_str());
		vector<pair<int, double>> vectorTempClust;
		double n = clustering.size();

		// Creater pair vector 
		for (auto &item:clustering){
			int node = item.first;
			double clust = item.second;
			vectorTempClust.push_back(make_pair(node, clust));
		}

		// Sort vector by KClust coeff in descending order
		sort(vectorTempClust.begin(), vectorTempClust.end(), sortbyValDescDouble);

		// Header for file
		fOut << "#NodeID,LocalClusteringCoeff,Degree,Rank,PkClust" << endl;

		for (int i = 0; i < vectorTempClust.size(); i++){
			int node = vectorTempClust[i].first;
			int degree = DictAdjList[node].size();
			double clust = vectorTempClust[i].second;
			double PKClust = (i + 1)/n;
			// Write to file
			fOut << node << "," << clust << "," << degree << "," << i + 1 << "," << PKClust << endl;
		}
	}


	map<int, double> KLocalCLustCoeff(map <int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that calculates the clustering coefficient for every node
		// This function differs from LocalClustCoeff() because KLocalCLustCoeff() calculates
		// the average clustering for every degree.
		// INputs: 
			// DictAdjList: dictionary containing the adjacency list of the network
			// allNodesBComponent: vector of ints containing all the nodes of all components of the network
		// Outputs:
			// kclustering: dictionary containing the k clustering for every degree in the network

		// Declare some variables
		map<int, double> kclustTemp, kclust;
		double LocalClustering;

		// Iterate over all nodes to obtain the sum of the average local clustering for each degree
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];
			int degree = DictAdjList[node].size();

			// Calculate the local clustering coefficient
			double numerator = LocalClustCoeff_Numerator(DictAdjList, node);
			double denominator = LocalClustCoeff_Denominator(DictAdjList, node);

			// Beware! the denominator can return 0!!
			if (denominator != 0){
				LocalClustering = numerator/denominator;
			}
			else {
				LocalClustering = 0.0;
			}

			// Add to dictionary
			// If the degree not in keys
			if (kclustTemp.find(degree) == kclustTemp.end()){
				kclustTemp[degree] = LocalClustering;
			}
			else{
				kclustTemp[degree] += LocalClustering;
			}
		}

		// Get the times a degree appear in the network 
		map<int, double> degreeCounts = CountDegrees(DictAdjList, allNodesBComponent);

		// Iterate over the degrees
		for (auto item : kclustTemp){
			int degree_ = item.first;
			// Get the times the degree appear
			double times = degreeCounts[degree_];
			double LocalCLust = item.second/times;

			// Append to the dictionary
			kclust[degree_] = LocalCLust;
		}

		return kclust; 
	}


	void saveKLocalClustCoeff(string fOutName, map<int, double>& kclust){
		// Function that saves the K clustering coefficient of all nodes into a file
		// INputs: 
			// fInName: name of the file to save
			// kclust: dictionary containing the K clustering coefficient of all nodes in thte network
			// DictAdjList: dictionary containing the adjacency list of the network
		// Outputs: 
			// None

		// Declare some variables and files
		ofstream fOut(fOutName.c_str());
		vector<pair<int, double>> vectorTempKClust;
		double k = kclust.size();

		// Creater pair vector 
		for (auto &item:kclust){
			int degree = item.first;
			double clust = item.second;
			vectorTempKClust.push_back(make_pair(degree, clust));
		}

		// Sort vector by KClust coeff in descending order
		sort(vectorTempKClust.begin(), vectorTempKClust.end(), sortbyValDescDouble);

		// Header for file
		fOut << "#Degree,KLocalCLustCoeff,Rank,PCk" << endl;

		for (int i = 0; i < vectorTempKClust.size(); i++){
			int degree = vectorTempKClust[i].first;
			double clust = vectorTempKClust[i].second;
			double PCk = (i + 1)/k;
			// Write to file
			fOut << degree << "," << clust << "," << i + 1 << "," << PCk << endl;
		}
	}





	tuple<stack<int>, unordered_map<int, vector<int>>, unordered_map <int, double>, unordered_map<int, int>> BreadthFirstSearchWeights(unordered_map <int, vector<int>>& DictAdjList, int startNode){
		// Function that implements the Breadth First search algorithms given an adjacency list
		// Breadth-first search find the shortest paths from a starting vertex startNode to every other vertex
		// in the same component as startNode. This algorithm explores nodes in layers, where each layer is
		// at a distancde d from the starting vertex d. 
		// ALso, this algorithms is modified to calculate the weight of the Betweenness centrality of nodes, 
		// assuming that there are multiple shortest paths. 
		// Algortihm explanation in section 10.3.1, 10.3.2, 10.3.3, 10.3.4, 10.3.5, and 10.3.6 from
		// Newman's book: Networks an Introduction, 2010 
		// Pseudocode inspired from: https://www.youtube.com/watch?v=oDqjPvD54Ss
		// Inputs:
			// DictAdjList: Adjacency list as a dictionary
			// startNode: starting node from which find all distances
		// Outputs
			// tuple: containing the visitedNodes (vector), shortestPahts(dictionary), weights(dictionary), and distances(dictionary)

		// Declare variables and structures
		// Dictionary with keys as nodes and values as the parent node of the key.
		// The dictionary will have n entries (i.e. same as adjacency list), and will contain
		// the shortest path from startNode to all node in the component by backtracking
		unordered_map<int, vector<int>> prev;
		// Dictionary with keys as nodes and values as the distance from startNode
		unordered_map<int, int> distances;
		// Dictionary with the weights of each node. THe weights represent the number
		// of distinct shortest paths of a vertex
		unordered_map <int, double> weights; 
		// Vector to save the visited nodes
		stack<int> visitedNodes;

		// Initialize the vectors
		for (auto &item : DictAdjList){
			int node = item.first;
			// Parent nodes are unknown at the beginning
			vector<int> temp;
			prev[node] = temp;
			// Distances are -1 at the beginning
			distances[node] = -1;
			// Weights are -1 at the beginning
			weights[node] = -1;
		}
		
		// Initialize the vector of vector to visit
		queue<int> Q;
		// Append the starting node to queue vector. This is the startNode is in the queue
		Q.push(startNode);
		// Set the distance for the start node
		distances[startNode] = 0;
		// Set the weight for the start node
		weights[startNode] = 1;
		// Set the parent for start node. This is dumb but useful ;)
		prev[startNode].push_back(-1);

		// Begin search
		while (!Q.empty()){
			// 0 because we always delete the first item
			int nodeSearch = Q.front();
			// Remove nodeSearch from the queue vector
			Q.pop();
			// Append the visited node
			visitedNodes.push(nodeSearch);
			// Get the neighbors of nodeSearch
			vector<int>::iterator it;

			// Iterate over the neighbors of nodeSearch
			for (it = DictAdjList[nodeSearch].begin(); it != DictAdjList[nodeSearch].end(); it++){
				int nextNeighbor = *it;

				// Check if a neighbor has a distance higher than nodeSearch
				// This condition adds multiple shortest paths
				if (distances[nextNeighbor] == distances[nodeSearch] + 1){
					// Update the weight
					weights[nextNeighbor] += weights[nodeSearch];
					// Add that neighbor to prev
					prev[nextNeighbor].push_back(nodeSearch);
				}

				// Check if nextNeighbor has not been visited
				if (distances[nextNeighbor] == -1){
					// Add nextNeighbor to the queue
					Q.push(nextNeighbor);
					// Set the distance
					distances[nextNeighbor] = distances[nodeSearch] + 1;
					// Set the weight
					weights[nextNeighbor] = weights[nodeSearch];
					// Add nodeSearch (i.e. parent node) to prev, so the path can be reconstructed
					prev[nextNeighbor].push_back(nodeSearch);
				}
			}
		}
		return {visitedNodes, prev, weights, distances};
	}

	unordered_map<int, double> BetweennessDeependenciesStartNode(unordered_map<int, double>& betweenness, stack<int>& visitedNodes, unordered_map<int, vector<int>>& shortestPaths, unordered_map<int, double>& weights, int startNode){
		// Function that calculates all pairwise dependencies in a component, starting at startNode. 
		// Useful to calculate the Betweeness centrality of nodes
		// Multiple shortest paths are considered. 
		// ALgorithm explanation in section 10.3.6 of Newman's book: Networks an Introduction, 2010 
		// Pseudocode from Brandes, 2008: On variants of shortest-path betweenness centrality and their generic computation
		// INputs:
			// betweenness: dictionary containing the betweenness values of nodes
			// visitedNodes: vector containing the order in wchich the nodes were visited by the Breadth-First search algorithm
			// shortestPaths: dictionary containing the shortest paths off al the nodes in a component starting at startNode
			// weights: dictionary containing the weights of all nodes with respect to a starting node. This is the number of shortest paths of each node
			// startNode: string indicating the starting node. Beware! the start node must be the same as the start node of the Breadth-First search algorithm
		// Outputs:
			// betweenness: dictionary containing the betweenness values of all nodes

		// Declare some variables
		unordered_map <int, double> dependencies;

		for (auto &item : weights){
			int node = item.first;
			// Scores are zero at the beginning
			dependencies[node] = 0;
		}

		betweenness[startNode] += visitedNodes.size() - 1;

		// Main Loop: backtrack the visited nodes
		while (!visitedNodes.empty()){
			int nodeSearch = visitedNodes.top();
			// Remove nodeSearch from the visitedNodes vector
			visitedNodes.pop();

			double numerator = (1 + dependencies[nodeSearch])/weights[nodeSearch];

			// Get the parents of nodeSearch using an iterator
			vector<int>::iterator it;
			for (it = shortestPaths[nodeSearch].begin(); it != shortestPaths[nodeSearch].end(); it++){
				int parent = *it;
				if (parent != -1){
					dependencies[parent] += weights[parent] * numerator;
				}
			}

			if (nodeSearch != startNode){
				betweenness[nodeSearch] += dependencies[nodeSearch] + 1;
			}
		}

		return betweenness;	 
	}

	unordered_map<int, double> BetweennessCentrality(unordered_map<int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent, int numNodes = -1){
		// Function that calculates the betweenness centrality of all nodes in a network, according to Brandes' algorithm
		// Computed in two phases:
		// 1. Breadth-First search in which visitedNodes, weights, and shortest paths of a startNode are computed
		// 2. Calculation of dependencies of startNode
		// 3. Normalization: Scale factor as defined by the NetworkX package
		// ALgorithm explanation in section 10.3.6 of Newman's book: Networks an Introduction, 2010 
		// Pseudocode from Brandes, 2008: On variants of shortest-path betweenness centrality and their generic computation
		// Inputs: 
			// DictAdjList: dictionary containing the adjacency list of the network	
			// numNodes: int indicating the number of nodes to approximate the betweenness. The higher the better!
			// allnodesBComponent: vector of ints containing the nodes of the biggest component
		// Outputs: 
			// betweenness: dictionary containing the betweenness values of all nodes in a network

		// Declare some variables
		double entries; 
		unordered_map<int, double> betweenness;
		vector<int> sampledNodes;

		// Some info to print
		cout << "\tThe biggest component contains: " << allNodesBComponent.size() << " nodes" << endl;


		if (numNodes == -1){
			cout << "\tUsing all nodes to calculate the betweenness" << endl;

			// Get all the nodes from the Adjacency list
			for (int i = 0; i < allNodesBComponent.size(); i++){
				int node = allNodesBComponent[i];
				sampledNodes.push_back(node);
			}
		}

		else{
			cout << "\tUsing " << numNodes << " nodes to calculate the betweenness" << endl;

			// Randomly sample numNodes
			sampledNodes = randomSampleNodes(numNodes, allNodesBComponent);

		}


		// Initialize structures
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];
			// Betweenness values are zero at the beginning
			betweenness[node] = 0;
		}

		entries = allNodesBComponent.size();

		// Iterate over all the nodes in DictAdjList or the sampled nodes
		int count = 0;
		for (int i = 0; i < sampledNodes.size(); i++){
			int startNode = sampledNodes[i];

			// Use the Breadth-First Search algorithm
			auto [visitedNodes, shortestPaths, weights, distances] = BreadthFirstSearchWeights(DictAdjList, startNode);
			//cout << "\tBFS for node: " << count << " calculated!" << endl;
			
			// Calculate dependencies
			betweenness = BetweennessDeependenciesStartNode(betweenness, visitedNodes, shortestPaths, weights, startNode);
			count += 1;
			cout << "\tBetweeness for node: " << count << " calculated!" << endl;  
		}

		cout << "\tBetweeness for all nodes calculated" << endl;

		// Normalization
		// Define scale factor
		double scaleFactor = 1 / (entries * (entries - 1));	
		// Iterate over all the vectors
		for (auto &item : betweenness){
			int node = item.first;
			betweenness[node] *= scaleFactor;
		}

		cout << "\tNormalization of betweenness values done!" << endl;

		return betweenness;
	}

	void saveBetweenness(string fOutName, unordered_map<int, double>& betweenness, unordered_map<int, vector<int>>& DictAdjList){
		// Function that saves the betweenness centrality of nodes in a network into a file
		// INputs:
			// fOutName: file to save the betweenness centrality
			// betweenness: dictionary containing the betweenness centrality of all nodes
		// Outputs:
			// Node

		// Declare some variables and files
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "#NodeID,Betweenness,Degree" << endl;

		// Iterate over the dictionary
		for (auto item : betweenness){
			int node = item.first;
			double bet = item.second;
			double degree = DictAdjList[node].size();
			fOut << node << "," << bet << "," << degree << endl;
		}
	}


	vector<int> BreadthFirstSearchComponents(map<int, vector<int>>& DictAdjList, int startNode){
		// Function that finds the connected components of a network using the Breath First Search ALgorithm
		// INputs:
			// DictAdjList: dictionary containing the adjacency list of the network
			// startNode: start node for the Breadth-First Search algorithm
		// Outputs:
			// visited: boolean vector containing the visited nodes

		// Declare some variables 
		unordered_map<int, bool> visited;
		vector<int> visitedNodes;

		// Initialize structures
		for (auto &item : DictAdjList){
			int node = item.first;
			// Visited nodes are unkown at the beginning
			visited[node] = false;
		}

		// Initialize the vector of vector to visit
		queue<int> Q;
		// Append the starting node to queue vector. This is the startNode is in the queue
		Q.push(startNode);
		// Set to true the visited node
		visited[startNode] = true;

		// Begin search
		while (!Q.empty()){
			// 0 because we always delete the first item
			int nodeSearch = Q.front();
			// Remove nodeSearch from the queue vector
			Q.pop();
			// Add to visitedNodes vector
			visitedNodes.push_back(nodeSearch);
			// Get the neighbors of nodeSearch
			vector<int>::iterator it;

			// Iterate over the neighbors of nodeSearch
			for (it = DictAdjList[nodeSearch].begin(); it != DictAdjList[nodeSearch].end(); it++){
				int nextNeighbor = *it;
				// Check if nextNeighbor has not been visited
				if (visited[nextNeighbor] == false){
					// Add nextNeighbor to the queue
					Q.push(nextNeighbor);
					// Set to true the visited node
					visited[nextNeighbor] = true;
				}
			}
		}

		return visitedNodes;
	}

	vector<vector<int>> FindComponents(map<int, vector<int>>& DictAdjList){
		// Function that finds the connected components of a network given an adjacency list
		// INputs:
			// DictAdjList: dictionary containing the adjacency list of the network
		// Outputs:
			// nodesComponents: vector containing the nodes for each component

		// Declare some variables
		vector<int> visitedNodesComp;
		//vector<int> startNodesComponents;
		vector<vector<int>> allNodesComponents;

		// Iterate over all the entries of DictAdjList
		for (auto &item : DictAdjList){
			int startNode = item.first;

			// Check if the node has been visited
			if ((find(visitedNodesComp.begin(), visitedNodesComp.end(), startNode) != visitedNodesComp.end()) == false){
				// Call BFS
				vector<int> visitedNodes = BreadthFirstSearchComponents(DictAdjList, startNode);
				// Apend to allNodesComponents
				allNodesComponents.push_back(visitedNodes);
				// Update the visited nodes
				visitedNodesComp.insert(visitedNodesComp.end(), visitedNodes.begin(), visitedNodes.end());
			}
		}

		return allNodesComponents;
	}

	void saveComponents(string fOutName, vector<vector<int>>& allNodesComponents){
		// Function that saves into a file the componentID and all the nodes of a component
		// Inputs:
			// fOutName: file to save the componentID and all the nodes
			// allNodesComponents: vector of vectors containing all the nodes of each component in the network

		// Declare some variables and files
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "#ComponentID:NodesComponent" << endl;

		// Iterate over the vector
		for (int i = 0; i < allNodesComponents.size(); i++){
			int compID = i;
			// Get all the nodes in a component
			vector<int> nodesComponent = allNodesComponents[i];
			// Save compID into file
			fOut << compID << ":";
			// Iterate over the nodes of a component
			for(int j = 0; j < nodesComponent.size(); j++){
				if (j != nodesComponent.size() - 1){
					fOut << nodesComponent[j] << ",";
				}
				else{
					fOut << nodesComponent[j];
				}
			}
			fOut << endl;
		}
	}

	int getIndexBiggestComponent(vector<vector<int>>& allNodesComponents){
		/*
		Function that returns the indest of the biggest component ofa network
		Inputs:
			allNodesComponents: vector of vectors of ints containing the nodes of each component
		Outputs:
			indexBComponent: int containing the index of the biggest component
		*/

		// Declare some variables
		int indexBComponent; 
		int maxSize = 0;

		// Iterate over the vector of components
		for (int i = 0; i < allNodesComponents.size(); i++){
			// Get the size of the component
			if (allNodesComponents[i].size() > maxSize){
				maxSize = allNodesComponents[i].size();
				// Update the index
				indexBComponent = i;
			}
		}

		return indexBComponent;
	}

	unordered_map<int, int> BreadthFirstSearchDiameter(unordered_map<int, vector<int>>& DictAdjList, int startNode){
		// Function that implements the Breadth First search algorithms given an adjacency list
		//
		// Inputs:
			// DictAdjList: Adjacency list as a dictionary
			// startNode: starting node from which find all distances
		// Outputs
			// distances: dictionary containing the distances of all nodes in a component from startNode

		// Declare variables and structures
		// Dictionary with keys as nodes and values as the distance from startNode
		unordered_map<int, int> distances;

		// Initialize the vectors
		for (auto &item : DictAdjList){
			int node = item.first;
			// Distances are -1 at the beginning
			distances[node] = -1;
		}
		
		// Initialize the vector of vector to visit
		queue<int> Q;
		// Append the starting node to queue vector. This is the startNode is in the queue
		Q.push(startNode);
		// Set the distance for the start node
		distances[startNode] = 0;

		// Begin search
		while (!Q.empty()){
			// 0 because we always delete the first item
			int nodeSearch = Q.front();
			// Remove nodeSearch from the queue vector
			Q.pop();

			// Get the neighbors of nodeSearch
			vector<int>::iterator it;
			// Iterate over the neighbors of nodeSearch
			for (it = DictAdjList[nodeSearch].begin(); it != DictAdjList[nodeSearch].end(); it++){
				int nextNeighbor = *it;

				// Check if nextNeighbor has not been visited
				if (distances[nextNeighbor] == -1){
					// Add nextNeighbor to the queue
					Q.push(nextNeighbor);
					// Set the distance
					distances[nextNeighbor] = distances[nodeSearch] + 1;
				}
			}
		}
		return distances;
	}

	int DiameterSingleComponent(unordered_map<int, vector<int>>& DictAdjList, vector<int>& nodesComponent){
		// Function that calculates the diameter of a compnent using the Breadth-First Search algorithm
		// Beware! This function is only useful if there is a single component in the adjacency list
		// Inputs:
			// DictAdjList: dictionary containing the adjacency list of the network
		// Outputs:
			// diameter: integer indicating the diameter of the component

		// Declare some variables
		int diameter = 0;

		// Iterate over all the nodes in nodesComponent
		for (int i = 0; i < nodesComponent.size(); i++){
			int startNode = nodesComponent[i];

			// Execute the Breadth-First search algorithm
			unordered_map<int, int> distances = BreadthFirstSearchDiameter(DictAdjList, startNode);
			
			// Get the maximum distance of all nodes in a component starting at startNode
			unordered_map<int, int>:: iterator tempMaxDist = max_element(distances.begin(), distances.end(), maxValueDict);
			int tempMaxDist_ = tempMaxDist-> second;

			// Update the diameter
			if (tempMaxDist_ > diameter){
				diameter = tempMaxDist_;
			}

			cout << "\tCalculating diameter from node: " << i << endl;	
		}

		return diameter;
	}

	unordered_map<int, int> DiameterAllComponent(unordered_map<int, vector<int>>& DictAdjList, vector<vector<int>>& allNodesComponents){
		// Function that calculates the diameter for all components in a network
		// Inputs
			// DictAdjList: dictionary containing the adjacency list of the network
			// allNodesComponents: vector of vectors, where each element contains the nodes of a component
		// Outputs:
			// diameterComponents: dictionary containing the diameter of each component. The order of components
			// is given by the order in wchich they appear in the allNodesComponents vector

		// Declare some variables
		unordered_map<int, int> diameterComponents;

		// Iterate over the components
		for (int i = 0; i < allNodesComponents.size(); i++){

			cout << "\tAt component: " << i << endl;

			// Calculate the diameter
			int diameter = DiameterSingleComponent(DictAdjList, allNodesComponents[i]);
			int compID = i;
			// Append Diameter to dictionary
			diameterComponents[compID] = diameter;
		}

		return diameterComponents;
	}

	int DiameterBiggestComponent(unordered_map<int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		// Function that calculates the diameter for the biggest component in a network
		// Inputs
			// DictAdjList: dictionary containing the adjacency list of the network
			// allNodesBComponent: vector of ints where each element contains the nodes of a component
		// Outputs:
			// diameterBComponent: int containing the diameter of the biggest component

		// Declare some variables
		//int DiameterBComponent;

		// Get the diameter
		int DiameterBComponent = DiameterSingleComponent(DictAdjList, allNodesBComponent);

		return DiameterBComponent;
	}

	void saveDiameters(string fOutName, map<int, int>& diameterComponents){
		// Function that saves into a file the diameter of each component in a network
		// Inputs:
			// diameterComponents: dictionary containing the diameter of each component

		// Declare some variables and files
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "#ComponentID,Diameter" << endl;

		// Iterate over the vector
		for (auto &item : diameterComponents){
			int compID = item.first;
			int diameter = item.second;
			// Save into file
			fOut << compID << "," << diameter << endl;	
		}
	}

	double CharacteristicPathLenghtSingleComponent(unordered_map<int, vector<int>>& DictAdjList, vector<int>& nodesComponent){
		// Function that calculate the characteristic path lenght for a single component in the network
		// INputs: 
			// DictAdjList: dictionary containing the adjacency list of the network
			// nodesComponents: vector containing the nodes of the component
		// Outputs:
			// pathLenght: characteristic path length of the component

		// Declare some variables
		double sumDistances = 0;

		// Iterate over all the nodes in nodesComponent
		for (int i = 0; i < nodesComponent.size(); i++){
			int startNode = nodesComponent[i];

			// Execute the Breadth-First search algorithm
			unordered_map<int, int> distances = BreadthFirstSearchDiameter(DictAdjList, startNode);
			
			// Sum all the distances
			for (auto &item : distances){
				// Ignore negative distance since they do not belong to the current component
				if (item.second != -1){
					sumDistances += item.second;
				}
			}

			cout << "\tCalculating characteristic path length at node: " << i << endl;
		}

		// Get the size of the component
		double n = nodesComponent.size();
		double pathLenght = sumDistances/(n * (n - 1));

		return pathLenght;
	}

	map<int, double> CharacateristicPathLenghtAllComponent(unordered_map<int, vector<int>>& DictAdjList, vector<vector<int>>& allNodesComponents){
		// Function that calculates the characteristic path length for all components in a network
		// Inputs
			// DictAdjList: dictionary containing the adjacency list of the network
			// allNodesComponents: vector of vectors, where each element contains the nodes of a component
		// Outputs:
			// pathLenghtComponents: dictionary containing the characteristic path length of each component. The order of components
			// is given by the order in wchich they appear in the allNodesComponents vector

		// Declare some variables
		map<int, double> pathLenghtComponents;

		// Iterate over the components
		for (int i = 0; i < allNodesComponents.size(); i++){
			// Calculate the characteristic path length 

			cout << "\tAt component: " << i << endl;

			double pathLenghtComp = CharacteristicPathLenghtSingleComponent(DictAdjList, allNodesComponents[i]);
			int compID = i;
			// Append path length to dictionary
			pathLenghtComponents[compID] = pathLenghtComp;
		}

		return pathLenghtComponents;
	}

	void saveCharacteristicPathLength(string fOutName, map<int, double>& pathLenghtComponents){
		// Function that saves into a file the characteristic path lenght of each component in a network
		// Inputs:
			// pathLenghtComponents: dictionary containing the characteristic path length of each component

		// Declare some variables and files
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "#ComponentID,CharacteristicPathLenght" << endl;

		// Iterate over the vector
		for (auto &item : pathLenghtComponents){
			int compID = item.first;
			double pathLenght = item.second;
			// Save into file
			fOut << compID << "," << pathLenght << endl;	
		}
	}

	double CharacateristicPathLenghtBiggestComponent(unordered_map<int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		/*
		Function that returns the characteristic/average path lenght of the biggest component of a network
		Inputs:
			DictAdjList: dictionary containing the adjacency list of the network
			allNodesBComponent: vector of ints where each element contains the nodes of a component
		Outputs:
			pathLenghtBComponent: double containing the characteristic/average oath lenght of the biggest component
		*/

		// Get the path lenght
		double pathLenghtBComponent = CharacteristicPathLenghtSingleComponent(DictAdjList, allNodesBComponent);

		return pathLenghtBComponent;
	}

	tuple<int, double, double> PathLengthDiameterBComponent(unordered_map<int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		/*
		Function that returns the characteristic/average path lenght and the diamter of the biggest component of a network
		Inputs:
			DictAdjList: dictionary containing the adjacency list of the network
			allNodesBComponent: vector of ints where each element contains the nodes of a component
		Outputs:
			diameter: diamter of the biggest component
			pathLenghtBComponent: double containing the characteristic/average oath lenght of the biggest component
		*/

		// Declare some variables
		int diameter = 0;
		double sumDistances = 0;
		vector<double> shortestPaths;
		double sumSTD = 0;

		// Iterate over all the nodes in nodesComponent
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int startNode = allNodesBComponent[i];

			// Execute the Breadth-First search algorithm
			unordered_map<int, int> distances = BreadthFirstSearchDiameter(DictAdjList, startNode);
			
			// Get the maximum distance of all nodes in a component starting at startNode
			unordered_map<int, int>:: iterator tempMaxDist = max_element(distances.begin(), distances.end(), maxValueDict);
			int tempMaxDist_ = tempMaxDist-> second;

			// Update the diameter
			if (tempMaxDist_ > diameter){
				diameter = tempMaxDist_;
			}

			// Sum all the distances
			//double tempSumDistances = 0;
			for (auto &item : distances){
				// Ignore negative distance since they do not belong to the current component
				if (item.second != -1){
					sumDistances += item.second;
					// Append distances to vector
					shortestPaths.push_back(item.second);
				}
			}
			cout << "\tCalculating characteristic path length/diameter from node " << i << endl;

		}

		// Get the size of the component
		double n = allNodesBComponent.size();
		double pathLenght = sumDistances/(n * (n - 1));

		// Calculate the standard deviation
		for (int i = 0; i < shortestPaths.size(); i++){
			sumSTD += pow((shortestPaths[i] - pathLenght), 2);
		}
		double pathLenghtSTD = sqrt((sumSTD)/(n * (n - 1)));

		return {diameter, pathLenght, pathLenghtSTD};

	}

	tuple<map<int, double>, double> EigenvectorCentrality(map<int, vector<int>>& DictAdjList, int maxIter = 10000, double tol = 1.0e-6){
		// Function that calculates the eigenvector centrality of the network given an adjacency list
		// This function uses the iterative power method given by
			// x_{t} = A^{t} * x_{t - 1}
		// The graph must be directed/undirected and aperiodic (i.e. cannot be splitted into a bipartite, tripartite)
		// the eigenvector calculated is the vector corresponding to the leading eigenvalue. 
		// The leading eigenvalue has its magnitude bigger than any the magnitude of any other eigenvalue of the adjacency matrix
		// and the corresponding eigenvector has all its entries non-negative (i.e. >= 0)
		// Inputs: 
			// DictAdjList: dictionary containing the djacency list of the network
			// initVector: initial vector for the iterative method (add in later versions!!)
			// maxIter: maximum number of iterations
			// tol: error tolerance
		// Output:
			// x_1: leading eigenvector
			// leadingEigenvalue: leading eigenvalue

		// Declare some variables
		map<int, double> x_1;
		map<int, double> x_0;
		double leadingEigenvalue;

		// Get the number of nodes
		double n = DictAdjList.size();
		cout << "Size: " << n << endl;

		// If no initVector, create it as a vector of ones (so the vector is not orthogonal to the leafing eigenvector)
		vector<double> initVector;
		for (int a = 0; a < n; a++){
				initVector.push_back(1);
		}

		// Get the norm of initVector 
		double norm = vectorNorm(initVector);
		// Normalize the initVector to avoid problems with arbitrarily large numbers
		// Also, convert into a dictionary
		//cout << "initial vector" << endl;
		int h = 0;
		for (auto &item : DictAdjList){
			int node = item.first;
			x_0[node] = initVector[h]/norm;
			h++;
			//cout << "node: " << node << " x_0: " << x_0[node] << endl;
		}

		//cout << "check initial x_0" << endl;
		//for (auto &item:x_0){
			//cout << "node: " << item.first << ":" << item.second << endl;
//		}

		// Iterate up to maxIter
		int i; // so it can be called outside the for loop
		for (i = 0; i < maxIter; i++){
			
			for (auto &item : x_0){
				int node_x = item.first;
				// Get the neighbors of node_x
				vector<int> neighbors = DictAdjList[node_x];
				//cout << "node_x: " << node_x << " size neighbors: " << neighbors.size() << " it: " << i << endl;
				// Get the power of the neighbors of node_x
				for (int j = 0; j < neighbors.size(); j++){
					x_1[node_x] += x_0[neighbors[j]];
					//cout << "neighbors: " << neighbors[j] << endl;
					//cout << "x_0[neighbors[j]]: " << x_0[neighbors[j]] << " iteration: " << i << " neighbors[j]: " << neighbors[j] << endl;
				}

				//cout << "x_1[node_x]: " << x_1[node_x] << endl;
			}


			// Calculate the leading eigenvalue
			leadingEigenvalue = 0;
			// Get the maximum value of the vector
			map<int, double>:: iterator tempMaxDist = max_element(x_1.begin(), x_1.end(), maxValueDictDouble);
			double tempMaxDist_ = tempMaxDist-> first;
			leadingEigenvalue = x_1[tempMaxDist_]/x_0[tempMaxDist_] - 1;

			//cout << "max x_1: " << x_1[tempMaxDist_] << " max x_0: " << x_0[tempMaxDist_] << endl;
			//cout << "\tleading eigenvalue: " << leadingEigenvalue << endl;


			// Normalize the vector, so arbitrarily large numbers does not appear
			// First, get the vector
			vector<double> vectorTemp;
			vectorTemp.reserve(x_1.size());
			transform (x_1.begin(), x_1.end(), back_inserter(vectorTemp), [] (pair<int, double> const & pair){
	        	return pair.second;
	        });

	        // Get the norm 
	        norm = vectorNorm(vectorTemp);
	        //cout << "\tnorm: " << norm << endl;
	        int k = 0;
	        for (auto &item : x_1){
	        	x_1[item.first] = vectorTemp[k]/norm;
	        	k++;
	        	//cout << "normalization x_1: " << item.first << ":" << x_1[item.first] << endl;
	        }
			
			// Check for convergence
			double sum = 0;
			for (auto &tempItem : x_1){
				sum += abs(x_1[tempItem.first] - x_0[tempItem.first]);
			}
			//cout << "SUM: " << sum << endl;
			if (sum < (n * tol)){
				return {x_1, leadingEigenvalue};
			}

			for (auto &tempItem : x_1){
				//cout << "x_1 check before: " << tempItem.first << ":" << tempItem.second << endl;
			}
			// Update vectors
			x_0 = x_1;

			//cout << "checking vector" << endl;

			for (auto &tempItem : x_1){
				//cout << "x_1 check after: " << tempItem.first << ":" << tempItem.second << endl;
			}

		}

		// Check if maxIter has been reached, then convergence has not been achieved
		if ( i == maxIter){
			// use print
			cout << "\tEigenvector failed to converge! Increase number of iterations!" << endl;
			return {x_1, leadingEigenvalue};
		}

		return {x_1, leadingEigenvalue};
	}

	void saveEigenvectorCentrality(string fOutName, map<int, double>& eigenvector){
		// Function that saves the centrality of the networks into a file
		// INputs:
			// fOutName: file in which the eigenvector centrality will be saved
			// eigenvector: dictionary containing the eigenvector centrality for all nodes
		// Outputes
			// None

		// Declare some variables
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "#NodeID,EigenvectorCentrality" << endl;

		// Iterate over the vector
		for (auto &item : eigenvector){
			int node = item.first;
			double eigenValue = item.second;
			// Save into file
			fOut << node << "," << eigenValue << endl;	
		}
	}

	tuple<map<int, double>, double> EigenvectorCentralityBComponent(map<int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent, int maxIter = 10000, double tol = 1.0e-6){
		/*
		Function that calculates the EigenVector Centrality of a single component of the network
		This function uses the iterative power method given by
			x_{t} = A^{t} * x_{t - 1}
		Inputs: 
			DictAdjList: dictionary containing the djacency list of the network
			allNodesBComponent: vector containing the nodes of the biggest component
			maxIter: maximum number of iterations
			tol: error tolerance
		Output:
			x_1: leading eigenvector of the component
			leadingEigenvalue: leading eigenvalue of the component
		*/

		// Declare some variables
		map<int, double> x_1;
		map<int, double> x_0;
		double leadingEigenvalue;

		// Get the number of nodes of the component
		double n = allNodesBComponent.size();

		// If no initVector, create it as a vector of ones (so the vector is not orthogonal to the leafing eigenvector)
		vector<double> initVector;
		for (int a = 0; a < n; a++){
				initVector.push_back(1);
		}

		// Get the norm of initVector 
		double norm = vectorNorm(initVector);
		// Normalize the initVector to avoid problems with arbitrarily large numbers
		// Also, convert into a dictionary
		for (int w = 0; w < allNodesBComponent.size(); w++){
			int node = allNodesBComponent[w];
			x_0[node] = initVector[w]/norm;
			//cout << "node: " << node << " x_0: " << x_0[node] << endl;
		}

		// Iterate up to maxIter
		int i; // so it can be called outside the for loop
		for (i = 0; i < maxIter; i++){
			
			for (auto &item : x_0){
				int node_x = item.first;
				// Get the neighbors of node_x
				vector<int> neighbors = DictAdjList[node_x];
				//cout << "node_x: " << node_x << " size neighbors: " << neighbors.size() << " it: " << i << endl;
				// Get the power of the negihbors of node_x
				for (int j = 0; j < neighbors.size(); j++){
					x_1[node_x] += x_0[neighbors[j]];
				}
			}

			// Calculate the leading eigenvalue
			leadingEigenvalue = 0;
			// Get the maximum value of the vector
			map<int, double>:: iterator tempMaxDist = max_element(x_1.begin(), x_1.end(), maxValueDictDouble);
			
			double tempMaxDist_ = tempMaxDist-> first;
			leadingEigenvalue = x_1[tempMaxDist_]/x_0[tempMaxDist_] - 1;

			// Normalize the vector, so arbitrarily large numbers does not appear
			// First, get the vector
			vector<double> vectorTemp;
			vectorTemp.reserve(x_1.size());
			transform (x_1.begin(), x_1.end(), back_inserter(vectorTemp), [] (pair<int, double> const & pair){
	        	return pair.second;
	        });

	        // Get the norm 
	        norm = vectorNorm(vectorTemp);
	        //cout << "\tnorm: " << norm << endl;

	        int k = 0;
	        for (auto &item : x_1){
	        	x_1[item.first] = vectorTemp[k]/norm;
	        	k++;
	        }
			
			// Check for convergence
			double sum = 0;
			for (auto &tempItem : x_1){
				sum += abs(x_1[tempItem.first] - x_0[tempItem.first]);
			}
			if (sum < (n * tol)){
				return {x_1, leadingEigenvalue};
			}

			// Update vectors
			x_0 = x_1;
		}

		// Check if maxIter has been reached, then convergence has not been achieved
		if ( i == maxIter){
			// use print
			cout << "\tEigenvector failed to converge! Increase number of iterations!" << endl;
			return {x_1, leadingEigenvalue};
		}

		return {x_1, leadingEigenvalue};

	}


	tuple<map<int, double>, double> EigenvectorCentralityBComponentProb(map<int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent, int maxIter = 10000, double tol = 1.0e-6){
		/*
		Function that calculates the EigenVector Centrality of a single component of the network
		This function uses the iterative power method given by
			x_{t} = A^{t} * x_{t - 1}
		Inputs: 
			DictAdjList: dictionary containing the djacency list of the network
			allNodesBComponent: vector containing the nodes of the biggest component
			maxIter: maximum number of iterations
			tol: error tolerance
		Output:
			x_1: leading eigenvector of the component
			leadingEigenvalue: leading eigenvalue of the component
		*/

		// Declare some variables
		map<int, double> x_1;
		map<int, double> x_0;
		double leadingEigenvalue;

		// Get the number of nodes of the component
		double n = allNodesBComponent.size();

		// If no initVector, create it as a vector of ones (so the vector is not orthogonal to the leafing eigenvector)
		vector<double> initVector;
		for (int a = 0; a < n; a++){
				initVector.push_back(1);
		}

		// Get the norm of initVector 
		double norm = vectorNormProb(initVector);
		// Normalize the initVector to avoid problems with arbitrarily large numbers
		// Also, convert into a dictionary
		for (int w = 0; w < allNodesBComponent.size(); w++){
			int node = allNodesBComponent[w];
			x_0[node] = initVector[w]/norm;
			//cout << "node: " << node << " x_0: " << x_0[node] << endl;
		}

		// Iterate up to maxIter
		int i; // so it can be called outside the for loop
		for (i = 0; i < maxIter; i++){
			
			for (auto &item : x_0){
				int node_x = item.first;
				// Get the neighbors of node_x
				vector<int> neighbors = DictAdjList[node_x];
				//cout << "node_x: " << node_x << " size neighbors: " << neighbors.size() << " it: " << i << endl;
				// Get the power of the negihbors of node_x
				for (int j = 0; j < neighbors.size(); j++){
					x_1[node_x] += x_0[neighbors[j]];
				}
			}

			// Calculate the leading eigenvalue
			leadingEigenvalue = 0;
			// Get the maximum value of the vector
			map<int, double>:: iterator tempMaxDist = max_element(x_1.begin(), x_1.end(), maxValueDictDouble);
			
			double tempMaxDist_ = tempMaxDist-> first;
			leadingEigenvalue = x_1[tempMaxDist_]/x_0[tempMaxDist_] - 1;

			// Normalize the vector, so arbitrarily large numbers does not appear
			// First, get the vector
			vector<double> vectorTemp;
			vectorTemp.reserve(x_1.size());
			transform (x_1.begin(), x_1.end(), back_inserter(vectorTemp), [] (pair<int, double> const & pair){
	        	return pair.second;
	        });

	        // Get the norm 
	        norm = vectorNormProb(vectorTemp);
	        //cout << "\tnorm: " << norm << endl;

	        int k = 0;
	        for (auto &item : x_1){
	        		x_1[item.first] = vectorTemp[k]/norm;
	        	k++;
	        }
			
			// Check for convergence
			double sum = 0;
			for (auto &tempItem : x_1){
				sum += abs(x_1[tempItem.first] - x_0[tempItem.first]);
			}
			if (sum < (n * tol)){
				return {x_1, leadingEigenvalue};
			}

			// Update vectors
			x_0 = x_1;
		}

		// Check if maxIter has been reached, then convergence has not been achieved
		if ( i == maxIter){
			// use print
			cout << "\tEigenvector failed to converge! Increase number of iterations!!!!" << endl;
			return {x_1, leadingEigenvalue};
		}

		return {x_1, leadingEigenvalue};

	}

	void saveEigenvectorCentralityBComponent(string fOutName, map<int, double>& eigenvector, map<int, vector<int>>& DictAdjList){
		// Function that saves the centrality of the networks into a file
		// INputs:
			// fOutName: file in which the eigenvector centrality will be saved
			// eigenvector: dictionary containing the eigenvector centrality for all nodes
			//DictAdjList: dictionary containing the adjacency list of the network
		// Outputes
			// None

		// Declare some variables
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "#NodeID,EigenvectorCentrality,Degree" << endl;

		// Iterate over the vector
		for (auto &item : eigenvector){
			int node = item.first;
			double eigenVect = item.second;
			int degree = DictAdjList[node].size();
			// Save into file
			fOut << node << "," << eigenVect << "," <<	degree << endl;
		}
	}

	void saveEigenvectorCentralityAllComponents(string fOutName, map<int, vector<int>>& DictAdjList, vector<vector<int>>& allNodesComponents){
		/*
		Function that saves the eigenvector centralities for all components in a network into a file
		Inputs:
			fOutName: string indicating the name of the file i n which to save
			// DictAdjList: dictionary containing the adjacency list of the network
			// allNodesComponents: vector of vectors, where each element contains the nodes of a component
		Outputs:
			None
		*/

		// Declare some variables and files
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "##NodeID,EigenvectorCentrality,Degree" << endl;

		// Iterate over the components
		for (int i = 0; i < allNodesComponents.size(); i++){
			// Calculate the eigenvector centrality 
			auto [eigenvectorCentrality, leadingEigenvalue] = EigenvectorCentralityBComponent(DictAdjList, allNodesComponents[i]);
			int compID = i;
			cout << "\tFor component: " << compID << " the eigenvalue is: " << leadingEigenvalue << endl;

			// Write into a file
			fOut << "#" << compID << "," << leadingEigenvalue << endl;
			// Iterate over the vector
			for (auto &item : eigenvectorCentrality){
				int node = item.first;
				double eigenValue = item.second;
				// Save into file
				fOut << node << "," << eigenValue << "," << DictAdjList[node].size() << endl;	
			}
		}
	}


	void saveInfoMetrics(string fOutName, double density, double averageDegree, double stdAverageDegree, double assortativity, double GlobalClustCoeff, double AverageClustCoeff, int diameterBComponent, double pathLenghtBComponent, double pathLenghtSTD, double leadingEigenvalue){
		// Function that saves some informtion about the metrics into a file
		// Inputs:
			// fOutname: file in which the information will be saved
			// density: density of the biggest component
			// averageDegree: average degree of the biggest component
			// assortativity: assortativity of the biggest component
			// GlobalClustCoeff: Global clustering coefficient of the biggest component
			// AverageClustCoeff: Average clustering coefficient of the biggest component
			// diameterBComponent: diameter off the biggest component
			// pathLenghtBComponent: characteristic/average path length of the biggest component
			// leadingEigenvalue: leading eigenvalue of the biggest component
		// Outputs:
			// None

		// Declare some variables 
		ofstream fOut(fOutName.c_str());

		// Header for the file
		fOut << "#Density,AverageDegree,stdAverageDegree,Assortativity,TransitivityGlobalClustCoeff,AverageGlobalClustCoeff,Diameter,CharacteristicPathLenght,pathLenghtSTD,LeadingEigenvalue" << endl;

		// Save into the file
		fOut << density << "," << averageDegree << "," << stdAverageDegree << "," << assortativity << "," << GlobalClustCoeff << "," << AverageClustCoeff << "," << diameterBComponent << "," << pathLenghtBComponent << "," << pathLenghtSTD << "," << leadingEigenvalue << endl;
	}


	map <int, double> nodesDegreeGreaterK(map<int, double>& degreeCounts){
		/*
		Function that calculates the number of nodes with a degree greater than k
		Inputs:
			degreeCounts: dictionary containing the times a degree appears in a network
		Outputs:
			degreeGreaterK: dictionary containing the number of nodes that appear with a degree greater than k
		*/

		// Declare some variables
		map<int, double> nodesDegreeGK, tempCumSum;
		double totalDegrees = 0;

		// Get the "penultimo" maximum degree
		vector<int> vectorTemp;
		vectorTemp.reserve(degreeCounts.size());
		transform (degreeCounts.begin(), degreeCounts.end(), back_inserter(vectorTemp), [] (pair<int, double> const & pair){
	        return pair.first;
	       }); 
		int maxDegree = vectorTemp[vectorTemp.size() - 1];

		// Create new vector of degree counts that contain all the possible degrees in the component
		map <int, double> newdegreeCounts;
		for(int i = 0; i <= maxDegree; i++){
			int degree = i;
			// Check if degree exists in the component
			if (degreeCounts.find(degree) != degreeCounts.end()){
				newdegreeCounts[degree] = degreeCounts[degree];
			}
			else{
				newdegreeCounts[degree] = 0;
			}

		}

		// Get a cumulative sum for the times each degree appear
		for (auto &item : newdegreeCounts){
			int degree = item.first;
			double counts = item.second;
			totalDegrees += counts;
			tempCumSum[degree] += totalDegrees;
		}

		// Get the number of nodes with a degree greater thank k
		for (auto &item : newdegreeCounts){
			// Ignore nodes that are greater than the maximum degree because they are zero
			if (totalDegrees - tempCumSum[item.first] > 1){
				nodesDegreeGK[item.first] = totalDegrees - tempCumSum[item.first];
			}
		}

		// Add values for degree 0
		//nodesDegreeGK[0] = totalDegrees;

		return nodesDegreeGK;
	}


	vector<tuple<int, int>> getEdges(map<int, vector<int>>& DictAdjList, vector<int> allNodesBComponent){
		/*
		Function that retrieves all the edges of the biggest component of a network
		Inputs:
			DictAdjList: dictionary containing the adjacency list of the network
			allNodesComponents: vector of vectors, where each element contains the nodes of a component
		Outputs:
			edges: vector of tuples containing the unique edges of the component
		*/

		// Declare some variables
		vector<tuple<int, int>> edges;

		// Iterate over the nodes of the component
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];

			// Get hte neighbors of current node
			vector<int> neighbors = DictAdjList[node];
			for (int j = 0; j < neighbors.size(); j++){
				int nodeNeighbor = neighbors[j];

				// Sort nodes to check if edge already exists
				if (node < nodeNeighbor){
					edges.push_back(make_tuple(node, nodeNeighbor));
				}
				else{
					edges.push_back(make_tuple(nodeNeighbor, node));
				}
			}
		}

		// Return unique edges
		set<tuple<int, int>> s(edges.begin(), edges.end());
		edges.assign(s.begin(), s.end());

		return edges;

	}

	vector<tuple<int, int>> degreesEdges(map<int, vector<int>>& DictAdjList, vector<tuple<int, int>>& edges){
		/*
		Function that return a vecotr of tuples containing the degrees of some given nodes
		Inputs: 
			DictAdjList: dictionary containing the adjacency list of the network
			edges: vector of tuples containing the unique edges of the component
		Outputs:
			degreesEdgesComponent: vector of tuples containing the degrees of the nodes of a component
		*/

		// Declare some variables
		vector<tuple<int, int>> degreesEdgesComponent;

		// Itertate over the nodes stored in edges
		for (auto [i, j]:edges){
			int source = i;
			int sink = j;
	
			// Get the degree
			int degreeSource = DictAdjList[source].size();
			int degreeSink = DictAdjList[sink].size();
	
			// Sort the degrees in each tuple
			if (degreeSource < degreeSink){
				degreesEdgesComponent.push_back(make_tuple(degreeSource, degreeSink));
			}
			else{
				degreesEdgesComponent.push_back(make_tuple(degreeSink, degreeSource));
			}
		}

		// Sort the vector of tuples
		sort(degreesEdgesComponent.begin(), degreesEdgesComponent.end(), greater<tuple<int, int>>());

		return degreesEdgesComponent;
	}


	map<int, double> richClubCoeff(map<int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		/*
		FUnction that calculates the rich club coefficient of the biggest component of a network, given by:
			\phi(k) = (2 * E_{k})/(N_{k} * (N_{k} - 1))
			where N_k is the number of nodes with degree larger than k, and E_k is the number of edges among those nodes.
			As stated in: https://arxiv.org/abs/physics/0701290
		Inputs:
			DictAdjList: dictionary containing the adjacency list of the network
			allNodesComponents: vector of vectors, where each element contains the nodes of a component 
		Outputs:
			rcCoeff: dictionary containing the righ club coeffiecient (value) for each degree k (key)
		*/

		// Declare some variables
		map<int, double> rcCoeff;

		// Get the number of nodes with a degree greater thank k
		map<int, double> degreeCounts = CountDegrees(DictAdjList, allNodesBComponent);
		map<int, double> nodesDegreeGK = nodesDegreeGreaterK(degreeCounts);

		// Get the degrees of the edges in the network
		auto edges = getEdges(DictAdjList, allNodesBComponent);
		auto degreesEdgesComponent= degreesEdges(DictAdjList, edges);

		// Number of edges with degree k
		double edgesK = degreesEdgesComponent.size();

		// Take elements from degreesEdgesComponent
		auto elements = degreesEdgesComponent.back();
		degreesEdgesComponent.pop_back();
		int k1 = get<0>(elements);
		int k2 = get<1>(elements);

		// Calculate the rich club coefficient 
		double d = -1;
		for (auto &item : nodesDegreeGK){
			d++;
			int degree = item.first;
			double nk = item.second;

			while (k1 <= d){
				if (degreesEdgesComponent.size() == 0){
					edgesK = 0;
					break;
				}
				// Take elements from degreesEdgesComponent
				auto elements = degreesEdgesComponent.back();
				degreesEdgesComponent.pop_back();
				k1 = get<0>(elements);
				k2 = get<1>(elements);

				edgesK -= 1;
			}

			// Rich club coefficient
			rcCoeff[degree] = (2 * edgesK)/(nk * (nk - 1));
		}


		return rcCoeff;

	}

	map<int, double> richClubCoeffUnc(map<int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		/*
		Function that calculates the uncorrelated rich club coefficient as:
			\phi_{unc}(k) = (k^2)/(<k>N)
		Inputs:
			DictAdjList: dictionary containing the adjacency list of the network
			allNodesComponents: vector of vectors, where each element contains the nodes of a component 
		Outputs:
			rcUncCoeff: dictionary containing the uncorrelated righ club coeffiecient (value) for each degree k (key)
		*/

		// Declare some variables
		map<int, double> rcUncCoeff;

		// Get the number of nodes with a degree greater thank k
		map<int, double> degreeCounts = CountDegrees(DictAdjList, allNodesBComponent);
		map<int, double> nodesDegreeGK = nodesDegreeGreaterK(degreeCounts);

		// Get the average degree of the component
		auto [averageDegree, stdAverageDegree] = AverageDegree(DictAdjList, allNodesBComponent);
		// Get the number of nodes of the component
		double n = allNodesBComponent.size();

		// Calculate the uncorrelated rich club coefficient
		for (auto &item : nodesDegreeGK){
			double degree = item.first;

			rcUncCoeff[degree] = (pow(degree, 2))/(averageDegree * n);
		}

		return rcUncCoeff;

	}

	map<int, double> rhoUnc(map<int, double> &rcCoeff, map<int, double> &rcUncCoeff){
		/*
		Function that calculate the normalizaed rich club coefficient as:
			\rho_{unc} = (\phi(k))/(\phi_{unc}(k)) * \phi_{unc}(k)
		Inputs: 
			rcCoeff: dictionary containing the righ club coeffiecient (value) for each degree k (key)
			rcUncCoeff: dictionary containing the uncorrelated righ club coeffiecient (value) for each degree k (key)
		*/

		// Declare some variables
		map<int, double> rhoUncCoeff;

		// Iterate over the degrees
		for (auto &item : rcCoeff){
			int degree = item.first;
			if (degree == 0){
				rhoUncCoeff[degree] = 0;
			}
			else{
				rhoUncCoeff[degree] = ((rcCoeff[degree])/(rcUncCoeff[degree])) * rcUncCoeff[degree]; 
			}
		}

		return rhoUncCoeff;
	}

	void saverhoUncCoeff(string fOutName, map<int, double>& rhoUncCoeff){
		/*
		Function that saves the normalized uncorrelated rich club coefficient
		Inputs:
			fOutname: file in which the information will be saved
			rhoUncCoeff: dictionary containing the uncorrelated righ club coeffiecient (value) for each degree k (key)
		Outputs:
			None
		*/

		// Declare some variables
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "Degree,RCUncNorm" << endl;

		// Iterate over the degrees
		for (auto &item : rhoUncCoeff){
			fOut << item.first << "," << item.second << endl;
		}


	}

	tuple<double, double> saveOrderSizeNetwork(string fInName){
		/*
		Function that calculates the order and size of the entire network
		Inputs:
			fInname: file containing the network
		Outputs:
			orderNetwork: double containing the order of the network
			sizeNetwork: double containing the number of the network
		*/

		// Declare some variables
		gN::GNet gNet; // Luiño's Network class
		double sumDegrees = 0;

		// Load the network
		gNet.loadGNet(fInName.c_str());

		// Build the dictionary
		for (int i=0; i<gNet.getNGenotypes(); i++){
			// Ignore nodes that have no connections
			if (gNet.getGConnections(i).size() != 0){
				sumDegrees += gNet.getGConnections(i).size();
			}
		}

		double sizeNetwork = sumDegrees/2;
		double orderNetwork = gNet.getNGenotypes();

		return {orderNetwork, sizeNetwork};
	}

	tuple<double, double> saveOrderSizeBComponent(map<int, vector<int>>& DictAdjList, vector<int>& allNodesBComponent){
		/*
		Function that calculates the order (number of nodes) and size (number of edges) of the biggest component of a network
		Inputs:
			DictAdjList: dictionary containing the adjacency list of the network
			allNodesComponents: vector of vectors, where each element contains the nodes of a component 
		Outputs:
			orderBComponent: double containing the order of the biggest component
			sizeBComponent: double containing the number of the biggest component
		*/

		// Declare some variables
		double orderBComponent = 0;
		double sumDegrees = 0;

		// Get the order of the biggest component
		orderBComponent = allNodesBComponent.size();

		// Get the size of the biggest component
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];
			double degreeNode = DictAdjList[node].size();
			sumDegrees += degreeNode;
		}

		double sizeBComponent = sumDegrees/2;

		return {orderBComponent, sizeBComponent};
	}

	map<int, int> findDistancesHub(unordered_map<int, vector<int>>& DictAdjList, int startNode, vector<int>& allNodesBComponent){
		/*
		Function that finds the levels of the genotype networks. The levels are at a Hamming distance d
		from the main hub of the network (startNode). The breadth-first search algorith is useful to find the
		the neighbors at a distance d from the startNode.
		Inputs:
			DictAdjList: dictionary containing the adjacency list of the network 
			startNode: node from which the search will start
			allNodesBComponent: vector containing the nodes of the biggest component
		Ouputs:
			distances: dictionary containing all the distances of nodes from the startNode 
		*/


		// Declare variables and structures
		// Dictionary with keys as nodes and values as the distance from startNode
		map<int, int> distances;

		// Initialize the vectors
		for (int i = 0; i < allNodesBComponent.size(); i++){
			int node = allNodesBComponent[i];
			// Distances are -1 at the beginning
			distances[node] = -1;
		}

		// Initialize the vector of vector to visit
		queue<int> Q;
		// Append the starting node to queue vector. This is the startNode is in the queue
		Q.push(startNode);
		// Set the distance for the start node
		distances[startNode] = 0;

		// Begin search
		while (!Q.empty()){
			// 0 because we always delete the first item
			int nodeSearch = Q.front();
			// Remove nodeSearch from the queue vector
			Q.pop();

			// Get the neighbors of nodeSearch
			vector<int>::iterator it;
			// Iterate over the neighbors of nodeSearch
			for (it = DictAdjList[nodeSearch].begin(); it != DictAdjList[nodeSearch].end(); it++){
				int nextNeighbor = *it;

				// Check if nextNeighbor has not been visited
				if (distances[nextNeighbor] == -1){
					// Add nextNeighbor to the queue
					Q.push(nextNeighbor);
					// Set the distance
					distances[nextNeighbor] = distances[nodeSearch] + 1;
				}
			}
		}

		return distances;

	}

	void saveDistancesHub(string fOutName, map<int, int>& distances){
		/*
		Function that saves the distances from the hub
		Inputs: 
			distances: dictionary containing the distances of all nodes to the main hub
						of the network
			Outputs: None
		*/

		// Declare some variables
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "#NodeID,Distance" << endl;

		// Iterate over the degrees
		for (auto &item : distances){
			fOut << item.first << "," << item.second << endl;
		}
	}

	map<int, vector<int>> levelsHub(map<int,int>& distances){
		/*
		Function that gets all the nodes of a given level.
		Inputs:
			distances: dictionary containing the distances of all nodes
			to the main hub of the network
		Outputs:
			levels: dictionary containing the distance (key) and nodes belonging to a level
			(values). The distances defines the level. Look at the results, section 3.2.5 of Henry's TFM.
		*/

		// Declare structures
		map<int, vector<int>> levels;

		// Iterate over the distances
		for (auto &item : distances){
			int node = item.first;
			int dist = item.second;

			// Add nodes belonging to a level
			levels[dist].push_back(node);

		}

		return levels;
	}

	void saveLevelsHub(string fOutName, map<int, vector<int>>& levels){
		/*
		Functionat that saves all the nodes of the network but separated by levels
		given by the hamming distance
		Inputs:
			levels: dictionary containing the distance (key) and nodes belonging to a level
			(values).
		Outputs: None
		*/

		// Declare some variables
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "#Level,Nodes" << endl;

		// Iterate over the dictionary
		for (auto &item : levels){
			int level = item.first;
			vector<int> nodesLevel = item.second;

			fOut << level << ":";

			for (int i =0; i < nodesLevel.size(); i++){
				if (i != nodesLevel.size() - 1){
					fOut << nodesLevel[i] << ",";
				}
				else{
					fOut << nodesLevel[i] << endl;
				}
			}
		}

	}


	void averageDegreeLevels(map<int, vector<int>>& DictAdjList, map<int, vector<int>>& levels, string fOutName){
		/*
		Function that calculates the average <k> degree of each level
		Inputs:
			DictAdjList: dictionary containing the adjacency list of the network
			levels: dictionary containing the distance (key) and nodes belonging to a level
			(values).
			fOutName: name file to save the average degree of levels
		Outputs:
			averageKLevels: dictionary containing the average degree and its standard deviation for each level
		*/

		// Declare structures
		ofstream fOut(fOutName.c_str());

		// Header for file
		fOut << "#Level,AverageDegree,STD" << endl;

		// Iterate over the levels
		for (auto &item : levels){
			int level = item.first;
			vector<int> nodesLevel = item.second;

			// Calculte the average degree
			auto [averageDegree, stdAverageDegree] = AverageDegree(DictAdjList, nodesLevel);

			// Save into file
			fOut << level << "," << averageDegree << "," << stdAverageDegree << endl;
		}
	}



}