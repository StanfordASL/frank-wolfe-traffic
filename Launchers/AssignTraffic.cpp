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
#include "Algorithms/TrafficAssignment/ObjectiveFunctions/CombinedEquilibrium.h"
#include "Algorithms/TrafficAssignment/TravelCostFunctions/BprFunction.h"
#include "Algorithms/TrafficAssignment/TravelCostFunctions/ModifiedBprFunction.h"
#include "Algorithms/TrafficAssignment/FrankWolfeAssignment.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/CommandLine/CommandLineParser.h"

void printUsage() {
	std::cout <<
		"Usage: AssignTraffic [-obj <objective>] [-f <func>] [-a <algo>] [-n <num>] [-ce <num>] -i <file> -od <file> [-o <path>]  \n"
		"This program assigns OD-pairs onto a network using the Frank-Wolfe method. It\n"
		"supports different objectives, travel cost functions and shortest-path algos.\n"
		"  -obj				objective function:\n"
		"						sys_opt (default), user_eq, combined_eq\n"
		"  -f <func>		travel cost function:\n"
		"						bpr (default) modified_bpr\n"
		"  -a <algo>		shortest-path algorithm:\n"
		"						dijkstra (default)\n"
		"  -n <num>         number of iterations (default = 100)\n"
		"  -ce_param <num>	combined_eq interpolation parameter in [0,1]:\n"
		"					0 for UE, 1 for SO\n"
		"  -i <path>        input graph edge CSV file\n"
		"  -od <file>       OD-pair file\n"
		"  -o <path>        output path\n"
		"  -v               display informative messages\n"
		"  -help            display this help and exit\n";  
}


// Assigns all OD-flows onto the input graph.
template <typename FrankWolfeAssignmentT>
void assignTraffic(const CommandLineParser& clp) {
	const std::string infilename = clp.getValue<std::string>("i");
	const std::string odFilename = clp.getValue<std::string>("od");
	const std::string outputPath = clp.getValue<std::string>("o");

	const double ceParameter = clp.getValue<double>("ce_param", 0.0);
	
	Graph graph(infilename, ceParameter);

	mkdir(&outputPath[0],0777); // create output folder
	
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

		const std::string objectiveFunction = clp.getValue<std::string>("obj", "sys_opt");
		
		if (objectiveFunction == "combined_eq")
			csv << "# Objective: " << objectiveFunction << "(" << ceParameter << ")\n";
		else 
			csv << "# Objective: " << objectiveFunction << "\n";
		
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
		csv << "obj_function_value,total_travel_cost\n";
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
	}
	else {
		throw std::invalid_argument("unrecognized shortest-path algorithm -- '" + algo + "'");
	}
}

// Picks the travel cost function according to the command line options.
template <template <typename> class ObjFunctionT>
void chooseTravelCostFunction(const CommandLineParser& clp) {
	const std::string func = clp.getValue<std::string>("f", "bpr");
	if (func == "bpr")
		chooseShortestPathAlgo<ObjFunctionT, BprFunction>(clp);
	else if (func == "modified_bpr")
	  chooseShortestPathAlgo<ObjFunctionT, ModifiedBprFunction>(clp);
	else
		throw std::invalid_argument("unrecognized travel cost function -- '" + func + "'");
}

// Picks the objective function according to the command line options.
void chooseObjFunction(const CommandLineParser& clp) {
	const std::string objectiveFunction = clp.getValue<std::string>("obj", "sys_opt");
	
	if (objectiveFunction == "sys_opt")
		chooseTravelCostFunction<SystemOptimum>(clp);
	else if (objectiveFunction == "user_eq")
		chooseTravelCostFunction<UserEquilibrium>(clp);
	else
		chooseTravelCostFunction<CombinedEquilibrium>(clp);
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
