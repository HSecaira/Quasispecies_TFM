/*
	Script that contains the functions for topology calculation
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
#include<unordered_map>
#include "functions.h"

using namespace std; 



int main() {
	// Defining variables and loading paths:
	cout << "\nLoading paths and findning files" << endl; 
	ostringstream dataPathIn, dataPathOut, sysStr, fInName, fOutName1;

	dataPathIn.str(""); 
	dataPathOut.str(""); 
	//dataPathIn << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Scripts/Topology/Test/";
	//dataPathOut << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Scripts/Topology/Test/Results/Betweenness/";
	dataPathIn << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/Reg1_experiments/SubNetworks/EL4-reg1_S4_L001_001_all_trim_merged_filter_sort_filter_length_collapsed/";
	dataPathOut << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/Reg1_experiments/SubNetworks/EL4-reg1_S4_L001_001_all_trim_merged_filter_sort_filter_length_collapsed/ResultsBetweenness/";
	//dataPathIn << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/GrandNetwork/Reg1/";
	//dataPathOut << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/GrandNetwork/Reg1/ResultsBetweenness/";
	//sysStr.str(""); 
	//sysStr << "rm -r " << dataPathOut.str(); 
	//system(sysStr.str().c_str()); 
	//sysStr.str(""); 
	//sysStr << "mkdir " << dataPathOut.str(); 
	//system(sysStr.str().c_str()); 
//	
	fInName.str("");
	// Suffix of the file (only if needed)
	string suffix = "GrandNetwork";
	string nameAdjList = "grandNetwork";
	fInName << dataPathIn.str() << nameAdjList << ".csv";
	fOutName1.str("");
	fOutName1 << dataPathOut.str() << "Bewteenness_" << suffix << ".csv";
	cout << "\tAll paths and files found!" << endl;

	// Start counting Execution Time
	auto startChrono = chrono::high_resolution_clock::now(); 

	// Using a dictionary to store the adjacency list
	//cout << "\nCreating a dictionary that contains the adjacency list" << endl;
	//unordered_map<int, vector<int>>  UnorderedDictAdjList = f::ConvertUnorderedDict(fInName.str());
	//map<int, vector<int>> DictAdjList = f::ConvertDict(fInName.str());
	//cout << "\tDictionary created!" << endl;

	cout << "\nCreating a dictionary that contains the adjacency list" << endl;
	map<int, vector<int>> DictAdjList = f::LoadGNetDict(fInName.str());
	unordered_map<int, vector<int>> UnorderedDictAdjList = f::LoadGNetUnorderedDict(fInName.str());
	cout << "\tDictionary created!" << endl;

	// Components of the network
	cout << "\nFinding components in the network" << endl;
	vector<vector<int>> allNodesComponents = f::FindComponents(DictAdjList);
	// Get the index of the biggest component
	int indexBComponent = f::getIndexBiggestComponent(allNodesComponents);
	// Get all the nodes of a biggest component
	vector<int> allNodesBComponent = allNodesComponents[indexBComponent];
	cout << "\tAll components found and saved into a file!" << endl;

	// Breadth-First Search and betwenness	
	cout << "\nCalculating betweenness centrality using Breadth-First search" << endl;
	// numNodes: is the number of nodes used to estimate the betweenness
	int numNodes = DictAdjList.size() * 1;
	unordered_map<int, double> betweenness = f::BetweennessCentrality(UnorderedDictAdjList, allNodesBComponent, numNodes);
	f::saveBetweenness(fOutName1.str(), betweenness, UnorderedDictAdjList);
	cout << "\tBetweenness calculated and saved into a file!" << endl;

	// Measure execution time
	auto endChrono = chrono::high_resolution_clock::now();
	// Get the difference
	auto durationChrono = endChrono - startChrono;

	cout << "\nExecution time: " << chrono::duration <double, milli> (durationChrono).count()/1000 << " s" << endl; 

	return 0; 
}