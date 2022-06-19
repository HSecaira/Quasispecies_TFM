/*
	Script for betweenness benchmark
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
	ostringstream dataPathIn, dataPathOut, sysStr, fInName, fOutName1, fOutName2;

	dataPathIn.str(""); 
	dataPathOut.str(""); 
	dataPathIn << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Scripts/Topology/Test/";
	dataPathOut << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Scripts/Topology/Test/Results/Betweenness_Benchmark/";
	//dataPathIn << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/Reg2_experiments/SubNetworks/EL13-reg2_S27_L001_001_all_trim_merged_filter_sort_filter_length_trim_collapsed/";
	//dataPathOut << "/home/henry/Downloads/Jae_Beca/NGS_all_seqs/GNR/Topology/Experiments/Reg2_experiments/SubNetworks/EL13-reg2_S27_L001_001_all_trim_merged_filter_sort_filter_length_trim_collapsed/ResultsBetweenness/";
	//sysStr.str(""); 
	//sysStr << "rm -r " << dataPathOut.str(); 
	//system(sysStr.str().c_str()); 
	//sysStr.str(""); 
	//sysStr << "mkdir " << dataPathOut.str(); 
	//system(sysStr.str().c_str()); 
	
	fInName.str("");
	// Suffix of the file (only if needed)
	string suffix = "_EL33-reg1";
	string nameAdjList = "grandNetwork_EL33-reg1";
	fInName << dataPathIn.str() << nameAdjList << ".csv";
	fOutName1.str("");
	fOutName1 << dataPathOut.str() << "Bewteenness_benchmark" << suffix << ".csv";
	fOutName2.str("");
	fOutName2 << dataPathOut.str() << "Bewteenness_allNodes" << suffix << ".csv";
	cout << "\tAll paths and files found!" << endl;

	// Start counting Execution Time
	auto startChrono = chrono::high_resolution_clock::now(); 

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
	cout << "Network contains " << DictAdjList.size() << " nodes" << endl;
	cout << "Biggest component contains " << allNodesBComponent.size() << " nodes" << endl;

	// Breadth-First Search and betwenness	
	cout << "\nCalculating betweenness centrality using Breadth-First search" << endl;
	// numNodes: is the number of nodes used to estimate the betweenness
	int numNodes = DictAdjList.size() * 1;
	unordered_map<int, double> betweenness_allNodes = f::BetweennessCentrality(UnorderedDictAdjList, allNodesBComponent, numNodes);
	f::saveBetweenness(fOutName2.str(), betweenness_allNodes, UnorderedDictAdjList);
	cout << "\tBetweenness calculated and saved into a file!\n" << endl;

	
	// Calculate betweenness using a different number of nodes for estimation
	vector<double> fracNodes{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
	map<double, double> MSE;
	map<double, double> MSEstd;

	for (int i = 0; i < fracNodes.size(); i++){
		cout << "\tUsing " << fracNodes[i] << " nodes to estimate Betweenness\n" << endl; 
		// Vector to store the calculated MSE
		vector<double> tempMSE;

		// Inner loop to calculate a mean and standard deviation
		for (int j = 0; j < 5; j++){
			// Estimate betweenness centrality
			int numNodesEst = DictAdjList.size() * fracNodes[i];
			unordered_map<int, double> betweennessEst = f::BetweennessCentrality(UnorderedDictAdjList, allNodesBComponent, numNodesEst);
			cout << "\tBetweenness calculated!" << endl;

			// Calculate the Mean Squared Error for estimation
			double error = 0;
			for (auto &item : betweenness_allNodes){
			int node = item.first;
			error += pow((betweenness_allNodes[node] - betweennessEst[node]), 2);
			}

			// Append vector to vector
			tempMSE.push_back(error);
		}
		
		// Calculate mean of MSE
		double sumMSE = 0;
		for (int k = 0; k < tempMSE.size(); k++){
			sumMSE += tempMSE[k];
		}
		double averageMSE = sumMSE/tempMSE.size();
		// Append to dictionary
		MSE[fracNodes[i]] = averageMSE;

		// Calculta standard deviation of MSE
		double sumSTD = 0;
		for (int k = 0; k < tempMSE.size(); k++){
			sumSTD += pow((tempMSE[k] - averageMSE), 2);
		}
		double STDMSE = sqrt((sumSTD)/(tempMSE.size() - 1));
		// Append to dictionary
		MSEstd[fracNodes[i]] = STDMSE;

	}

	//for (int i = 0; i < fracNodes.size(); i++){
		//cout << "\tUsing " << fracNodes[i] << " nodes to estimate Betweenness\n" << endl; 
//
		//if (fracNodes[i] != 1.0){
			//// Estimate betweenness centrality
			//int numNodesEst = DictAdjList.size() * fracNodes[i];
			//unordered_map<int, double> betweennessEst = f::BetweennessCentrality(UnorderedDictAdjList, allNodesBComponent, numNodesEst);
			//// Save into a file
			//ostringstream fOutNameEst;
			//fOutNameEst.str("");
			//fOutNameEst << dataPathOut.str() << "Bewteenness_" << numNodesEst << "Nodes" << suffix << ".csv";
			//f::saveBetweenness(fOutNameEst.str(), betweennessEst, UnorderedDictAdjList);
			//cout << "\tBetweenness calculated and saved into a file!" << endl;
//
			//// Calculate the Mean Squared Error for estimation
			//double error = 0;
			//for (auto &item : betweenness_allNodes){
				//int node = item.first;
				//error += pow((betweenness_allNodes[node] - betweennessEst[node]), 2);
			//}
//
			//MSE[fracNodes[i]] = error;
		//}
//
		//else{
			//MSE[i] = 0;
		//}
	//}



	// Measure execution time
	auto endChrono = chrono::high_resolution_clock::now();
	// Get the difference
	auto durationChrono = endChrono - startChrono;
	cout << "\nExecution time: " << chrono::duration <double, milli> (durationChrono).count()/1000 << " s" << endl; 


	// Save Estimation errors
	ofstream fOut(fOutName1.str());
	// Header for file
	fOut << "#FractionNodes,AverageMSE,STDMSE" << endl;

	// Iterate over the dictionary
	for (auto item : MSE){
		fOut << item.first << "," << item.second << "," << MSEstd[item.first] << endl;
	}


	return 0; 
}