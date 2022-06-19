/*
	Script that contains the functions for topological measures
*/

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
#include "functions.h"

using namespace std; 



int main() {
	// Defining variables and loading paths:
	cout << "\nLoading paths and finding files" << endl; 
	ostringstream dataPathIn, dataPathOut, sysStr, fInName, fOutName, fOutName2, fOutName3, fOutName4, fOutName5, fOutName6, fOutName7, fOutName8, fOutName9, fOutName10, fOutName11, fOutName12, fOutName13;

	dataPathIn.str(""); 
	dataPathOut.str(""); 
	//dataPathIn << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Scripts/Topology/Test/";
	//dataPathOut << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Scripts/Topology/Test/Results/NewTopology/";
	dataPathIn << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/Reg1_experiments/SubNetworks/EL21_reg1_all_all_trim_merged_filter_sort_filter_length_collapsed/";
	dataPathOut << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/Reg1_experiments/SubNetworks/EL21_reg1_all_all_trim_merged_filter_sort_filter_length_collapsed/ResultsTopology/";
	//dataPathIn << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/GrandNetwork/Reg2/";
	//dataPathOut << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/GrandNetwork/Reg2/ResultsTopology/";
	//sysStr.str(""); 
	//sysStr << "rm -r " << dataPathOut.str(); 
	//system(sysStr.str().c_str()); 
	//sysStr.str(""); 
	//sysStr << "mkdir " << dataPathOut.str(); 
	//system(sysStr.str().c_str()); 

	fInName.str("");
	// Suffix of the file (only if needed)
	string suffix = "GrandNetwork";
	string nameAdjList = "grandNetwork";
	fInName << dataPathIn.str() << nameAdjList << ".csv";
	fOutName.str("");
	fOutName << dataPathOut.str() << "Nodes_and_Degrees_" << suffix << ".csv";
	fOutName2 << dataPathOut.str() << "Components_" << suffix << ".csv";
	fOutName3.str("");
	fOutName3 << dataPathOut.str() << "Diameters_" << suffix << ".csv";
	fOutName4.str("");
	fOutName4 << dataPathOut.str() << "AverageNeighborDegree_" << suffix << ".csv";
	fOutName5.str("");
	fOutName5 << dataPathOut.str() << "KNearestNeighborDegree_" << suffix << ".csv";
	fOutName6.str("");
	fOutName6 << dataPathOut.str() << "ClusteringCoefficient_" << suffix << ".csv";
	fOutName7.str("");
	fOutName7 << dataPathOut.str() << "KClusteringCoefficient_" << suffix << ".csv";
	fOutName8.str("");
	fOutName8 << dataPathOut.str() << "CharacteristicPathLenght_" << suffix << ".csv";
	fOutName9.str("");
	fOutName9 << dataPathOut.str() << "EigenvectorCentrality_" << suffix << ".csv";
	fOutName10.str("");
	fOutName10 << dataPathOut.str() << "MetricsInformation_" << suffix << ".csv";
	fOutName11.str("");
	fOutName11 << dataPathOut.str() << "KNearestNeighbor_per_Degree_" << suffix << ".csv";
	fOutName12.str("");
	fOutName12 << dataPathOut.str() << "EigenvectorCentrality_allComponents_" << suffix << ".csv";
	//fOutName13.str("");
	//fOutName13 << dataPathOut.str() << "RcCoeffUncNorm_" << suffix << ".csv";
	cout << "\tAll paths and files found!" << endl;

	// Start counting Execution Time
	auto startChrono = chrono::high_resolution_clock::now(); 

	// Using a dictionary to store the adjacency list
	//cout << "\nCreating a dictionary that contains the adjacency list" << endl;
	//map<int, vector<int>>  DictAdjList = f::ConvertDict(fInName.str());
	//unordered_map<int, vector<int>>  DictAdjListUnordered = f::ConvertUnorderedDict(fInName.str());
	//cout << "\tDictionary created!" << endl;

	// Using a dictionary to store the adjacency list from the GNet class
	cout << "\nCreating a dictionary that contains the adjacency list" << endl;
	map<int, vector<int>> DictAdjList = f::LoadGNetDict(fInName.str());
	unordered_map<int, vector<int>> DictAdjListUnordered = f::LoadGNetUnorderedDict(fInName.str());
	cout << "\tDictionary created!" << endl;
	cout << "\tAdjacency list size " << DictAdjList.size() << endl;

	// Components of the network
	cout << "\nFinding components in the network" << endl;
	vector<vector<int>> allNodesComponents = f::FindComponents(DictAdjList);
	f::saveComponents(fOutName2.str(), allNodesComponents);
	// Get the index of the biggest component
	int indexBComponent = f::getIndexBiggestComponent(allNodesComponents);
	// Get all the nodes of a biggest component
	vector<int> allNodesBComponent = allNodesComponents[indexBComponent];
	cout << "\tAll components found and saved into a file!" << endl;

	//// Saving degree info into a file
	//cout << "\nSaving degree info" << endl;
	//f::saveNodesDegrees(fOutName.str(), DictAdjList, allNodesBComponent);
	//cout << "\tDegrees calculated and saved into a file!" << endl;
//
//
	//// Density of the graph
	//cout << "\nCalculating the density of the graph" << endl;
	//double density = f::Density(DictAdjList, allNodesBComponent);
	//cout << "\tTHe density is " << density << endl;
//
	//// Average degree 
	//cout << "\nCalculating the average degree" << endl;
	//auto [averageDegree, stdAverageDegree] = f::AverageDegree(DictAdjList, allNodesBComponent);
	//cout << "\tThe average degree is " << averageDegree << " +-: " << stdAverageDegree << endl;
//
	//// Assortativity/correlation coefficient
	//cout << "\nCalculating the assortativity of the network" << endl;
	//double assortativity = f::Assortativity(DictAdjList, allNodesBComponent);
	//cout << "\tThe assortativity/correlation coefficient is: " << assortativity << endl;
//
	//// Assortativity: Average Neighbor Degree 
	//cout << "\nCalculating the average neighbor degree (a type of assortativity)" << endl;
	//map<int, double> averageNeighborDegree = f::AverageNeighborDegree(DictAdjList, allNodesBComponent);
	//f::saveAverageNeighborDegree(fOutName4.str(), DictAdjList, averageNeighborDegree);
	//cout << "\tAverage neighbor degree calculated and saved into a file!" << endl;
//
	//// K Nearest neighbors
	//cout << "\nCalculating the K Nearest Neighbors (a type of assortativity)" << endl;
	//map<int, double> k_nearestNeighbors = f::KNearestNeighbors(DictAdjList, allNodesBComponent);
	//f::saveKNearestNeighbors(fOutName5.str(), k_nearestNeighbors);
	//f::savekNearestNeighborsDegree(fOutName11.str(), DictAdjList, allNodesBComponent);
	//cout << "\tK Nearest Neighbors calculated and saved into a file!" << endl;
//
	////Global Clustering Coefficient - Transitivity
	//cout << "\nCalculating the global clustering coefficient (transitivity)" << endl;
	//double GlobalClustCoeff = f::GlobalClusteringCoeff(DictAdjList, allNodesBComponent);
	//cout << "\tThe Global Clustering (transitivity) coefficient of the network is " << GlobalClustCoeff << endl;
//
	//// Average Global clustering coefficient
	//cout << "\nCalculating the average global clustering coefficient" << endl;
	//double AverageClustCoeff = f::AverageClusteringCoeff(DictAdjList, allNodesBComponent);
	//cout << "\tThe Average clustering coefficient of the network is " << AverageClustCoeff << endl;
//
	//// Local clustering coefficient
	//cout << "\nCalculating the Clustering Coefficient for every node" << endl;
	//map<int, double> clustering = f::LocalClustCoeff(DictAdjList, allNodesBComponent);
	//f::saveLocalClustCoeff(fOutName6.str(), clustering, DictAdjList);
	//cout << "\tClustering coefficients calculated and saved into a file" << endl;
//
	//// K clustering coefficient
	//cout << "\nCalculating the K Local clustering coefficent" << endl;
	//map<int, double> kclust = f::KLocalCLustCoeff(DictAdjList, allNodesBComponent);
	//f::saveKLocalClustCoeff(fOutName7.str(), kclust);
	//cout << "\tK Clustering coefficients calculated and saved into a file" << endl;
//
	//// Characteristic path length and diameter of the biggest component
	//cout << "\nCalculating the characteristic path length and diameter of the biggest component" << endl;
	//auto [diameterBComponent, pathLenghtBComponent, pathLenghtSTD] = f::PathLengthDiameterBComponent(DictAdjListUnordered, allNodesBComponent);
	////int diameterBComponent = -1;
	////double pathLenghtBComponent = -1;
	//cout << "\tDiameter of the biggest component is " << diameterBComponent << endl;
	//cout << "\tCharacteristic path length of of the biggest component is " << pathLenghtBComponent << " +- " << pathLenghtSTD << endl;

	// Eigenvector centrality
	cout << "\nCalculating eigenvector centrality for the biggest component" << endl;
	auto [eigenvectorCentrality, leadingEigenvalue] = f::EigenvectorCentralityBComponentProb(DictAdjList, allNodesBComponent);
	cout << "\tThe leading eigenvalue is " << leadingEigenvalue << endl;
	f::saveEigenvectorCentralityBComponent(fOutName9.str(), eigenvectorCentrality, DictAdjList);
	cout << "\tEigenvector calculated and saved into a file!" << endl;

	//// Community detection
//
	////Save metric information
	//cout <<"\nSaving Metrics information of the network" << endl;
	//f::saveInfoMetrics(fOutName10.str(), density, averageDegree, assortativity, GlobalClustCoeff, AverageClustCoeff, diameterBComponent, pathLenghtBComponent, leadingEigenvalue);
	//cout << "\tInformation saved into a file!" << endl;

	// Measure execution time
	auto endChrono = chrono::high_resolution_clock::now();
	// Get the difference
	auto durationChrono = endChrono - startChrono;

	cout << "\nExecution time: " << chrono::duration <double, milli> (durationChrono).count()/1000 << " s" << endl; 

	return 0; 
}