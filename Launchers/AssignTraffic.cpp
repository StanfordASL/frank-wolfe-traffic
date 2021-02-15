#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>
#include <bits/stdc++.h> 
#include <iostream> 
#include <sys/stat.h> 
#include <sys/types.h> 

#include <boost/dynamic_bitset.hpp>
#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/nested_dissection.h>

#include "Algorithms/TrafficAssignment/Adapters/DijkstraAdapter.h"
#include "Algorithms/TrafficAssignment/ObjectiveFunctions/SystemOptimum.h"
#include "Algorithms/TrafficAssignment/ObjectiveFunctions/UserEquilibrium.h"
#include "Algorithms/TrafficAssignment/TravelCostFunctions/BprFunction.h"
#include "Algorithms/TrafficAssignment/FrankWolfeAssignment.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/CommandLine/CommandLineParser.h"

void printUsage() {
	std::cout <<
		"Usage: AssignTraffic [-f <func>] [-a <algo>] -i <file> -od <file> [-o <file>]\n"
		"This program assigns OD-pairs onto a network using the Frank-Wolfe method. It\n"
		"supports different objectives, travel cost functions and shortest-path algos.\n"
		"  -so               find the system optimum (default: user equilibrium)\n"
		"  -v                display informative messages\n"
		"  -n <num>          the number of iterations (0 means use stopping criterion)\n"
		"  -f <func>         the travel cost function\n"
		"                      possible values:\n"
		"                        bpr (default) modified_bpr custom_bpr approx_bpr unaware_bpr\n"
		"                        davidson modified_davidson inverse\n"
		"  -a <algo>         the shortest-path algorithm\n"
		"                      possible values: dijkstra (default)\n"
		"  -i <path>         the input graph folder where the edge.csv file is found\n"
		"  -od <file>        the OD-pairs to be assigned\n"
		"  -o <path>         output path\n"
		"  -help             display this help and exit\n";  
}


// Assigns all OD-flows onto the input graph.
template <typename FrankWolfeAssignmentT>
void assignTraffic(const CommandLineParser& clp) {
	const std::string infilename = clp.getValue<std::string>("i");
	const std::string odFilename = clp.getValue<std::string>("od");
	const std::string outputPath = clp.getValue<std::string>("o");
	
	Graph graph(infilename);

	if (mkdir(&outputPath[0],0777) != 0)
		std::cout << "Could not create output folder\n";
	
	const std::string patternFilename = outputPath + "/flow";
	const std::string pathFilename = outputPath + "/paths";
	const std::string weightFilename = outputPath + "/weights";
	const std::string csvFilename = outputPath + "/output";
	
	std::vector<ClusteredOriginDestination> odPairs = importClusteredODPairsFrom(odFilename);

	const int numIterations = clp.getValue<int>("n");
	if (numIterations < 0) {
		const std::string msg("negative number of iterations");
		throw std::invalid_argument(msg + " -- " + std::to_string(numIterations));
	}
	
	std::ofstream csv;
	if (!csvFilename.empty()) {
		csv.open(csvFilename + ".csv");
		if (!csv.good())
			throw std::invalid_argument("file cannot be opened -- '" + csvFilename + ".csv'");
		csv << "# Input graph: " << infilename << "\n";
		csv << "# OD-pairs: " << odFilename << "\n";
		csv << "# Objective: " << (clp.isSet("so") ? "SO" : "UE") << "\n";
		csv << "# Function: " << clp.getValue<std::string>("f", "bpr") << "\n";
		csv << "# Shortest-path algo: " << clp.getValue<std::string>("a", "dijkstra") << "\n";
		csv << std::flush;
	}
	
	std::ofstream patternFile;
	if (!patternFilename.empty()) {
		patternFile.open(patternFilename + ".csv");
		if (!patternFile.good())
			throw std::invalid_argument("file cannot be opened -- '" + patternFilename + ".csv'");
		if (!csvFilename.empty())
			patternFile << "# Main file: " << csvFilename << ".csv\n";
		patternFile << "numIteration,tail,head,freeFlowCost,actualCost,capacity,flow\n";
	}

	std::ofstream pathFile;
	if (!pathFilename.empty()) {
		pathFile.open(pathFilename + ".csv");
		if (!pathFile.good())
			throw std::invalid_argument("file cannot be opened -- '" + pathFilename + ".csv'");
		if (!csvFilename.empty())
			pathFile << "# Main file: " << csvFilename << ".csv\n";
		pathFile << "numIteration,odPair,edges\n";
	}

	std::ofstream weightFile;
	if (!weightFilename.empty()) {
		weightFile.open(weightFilename + ".csv");
		if (!weightFile.good())
			throw std::invalid_argument("file cannot be opened -- '" + weightFilename + ".csv'");
		if (!csvFilename.empty())
			weightFile << "# Main file: " << csvFilename << ".csv\n";
		weightFile << "numIteration,weight\n";
	}

	FrankWolfeAssignmentT assign(graph, odPairs, csv, patternFile, pathFile, weightFile, clp.isSet("v"));

	if (csv.is_open()) {
		csv << "# Preprocessing time: " << assign.stats.totalRunningTime << "ms\n";
		csv << "iteration,customization_time,query_time,line_search_time,total_time,";
		csv << "avg_change,max_change,obj_function_value,total_travel_cost,checksum\n";
		csv << std::flush;
	}

	assign.run(numIterations);
}

// Picks the shortest-path algorithm according to the command line options.
template <template <typename> class ObjFunctionT, typename TravelCostFunction>
void chooseShortestPathAlgo(const CommandLineParser& clp) {
	const std::string algo = clp.getValue<std::string>("a", "dijkstra");
	if (algo == "dijkstra") {
		using Assignment = FrankWolfeAssignment<ObjFunctionT, TravelCostFunction, DijkstraAdapter>;
		assignTraffic<Assignment>(clp);
	} /* else if (algo == "bidijkstra") {
		using Assignment = FrankWolfeAssignment<
			ObjFunctionT, TravelCostFunctionT, trafficassignment::BiDijkstraAdapter, Graph>;
		assignTraffic<Assignment>(clp);
	} else if (algo == "ch") {
		using Assignment = FrankWolfeAssignment<
			ObjFunctionT, TravelCostFunctionT, trafficassignment::CHAdapter, Graph>;
		assignTraffic<Assignment>(clp);
		} else if (algo == "cch") {
		using Assignment = FrankWolfeAssignment<
			ObjFunctionT, TravelCostFunctionT, trafficassignment::CCHAdapter, Graph>;
			assignTraffic<Assignment>(clp);
	} */ else {
		throw std::invalid_argument("unrecognized shortest-path algorithm -- '" + algo + "'");
	}
}

// Picks the travel cost function according to the command line options.
template <template <typename> class ObjFunctionT>
void chooseTravelCostFunction(const CommandLineParser& clp) {
	const std::string func = clp.getValue<std::string>("f", "bpr");
	if (func == "bpr")
		chooseShortestPathAlgo<ObjFunctionT, BprFunction>(clp);
	/*else if (func == "modified_bpr")
	  chooseShortestPathAlgo<ObjFunctionT, ModifiedBprFunction>(clp);
		  else if (func == "custom_bpr")
	  chooseShortestPathAlgo<ObjFunctionT, CustomBprFunction>(clp);
	  else if (func == "approx_bpr")
	  chooseShortestPathAlgo<ObjFunctionT, ApproxBprFunction>(clp);
	  else if (func == "unaware_bpr")
	  chooseShortestPathAlgo<ObjFunctionT, UnawareBprFunction>(clp);
	  else if (func == "davidson")
	  chooseShortestPathAlgo<ObjFunctionT, DavidsonFunction>(clp);
	  else if (func == "modified_davidson")
	  chooseShortestPathAlgo<ObjFunctionT, ModifiedDavidsonFunction>(clp);
	  else if (func == "inverse")
	  chooseShortestPathAlgo<ObjFunctionT, InverseFunction>(clp);*/
	else
		throw std::invalid_argument("unrecognized travel cost function -- '" + func + "'");
}

// Picks the objective function according to the command line options.
void chooseObjFunction(const CommandLineParser& clp) {
	if (clp.isSet("so"))
		chooseTravelCostFunction<SystemOptimum>(clp);
	else
		chooseTravelCostFunction<UserEquilibrium>(clp);
}

int main(int argc, char* argv[]) {
	try {
		CommandLineParser clp(argc, argv);
		if (clp.isSet("help")) {
			printUsage();
			return EXIT_SUCCESS;
		}
		chooseObjFunction(clp);
	} catch (std::invalid_argument& e) {
		std::cerr << argv[0] << ": " << e.what() << std::endl;
		std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
