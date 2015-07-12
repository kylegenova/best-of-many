#include "run.h"
#include "statistics.h"
#include "util.h"
#include "christofides.h"
#include "heldkarp.h"
#include "lrs.h"
#include "edge_split.h"
#include "edge_join.h"
#include "swapround.h"
#include "statistics.h"
#include "gamma.h"
#include "random_walk.h"
#include "col_gen.h"

bool programArguments::parseProgram(string programFilePath)
{
	ifstream prog;
	prog.open(programFilePath);
	if (!prog.is_open())
	{
		cout << "Error: Could not open " << programFilePath << endl;
		return false;
	}
	programWorkingDirectory = getDirectoryOfPath(programFilePath);
	// Parse settings here
	while (!prog.eof())
	{
		string line;
		getline(prog, line);
		if (line == "INSTANCES:")
			break;
		// line[0] can't crash b/c of short circuit eval on whitespace check
		if (stringIsWhitespace(line) || line[0] == '#')
			continue;
		vector<string> option = splitString(line, ':');
		if (option.size() != 2)
		{
			cout << "Error: Incorrectly formated program at line:" << endl;
			cout << line << endl;
			return false;
		}
		if (option[0] == "CHRISTOFIDES")
		{
			USE_CFDS = option[1] == "TRUE";
		}
		else if (option[0] == "MAX_ENTROPY")
		{
			USE_ME = option[1] == "TRUE";
		}
		else if (option[0] == "COLUMN_GENERATION")
		{
			USE_CG = option[1] == "TRUE";
		}
		else if (option[0] == "EDGE_SPLITTING")
		{
			USE_ES = option[1] == "TRUE";
		}
		else if (option[0] == "SWAPROUND")
		{
			USE_SR = option[1] == "TRUE";
		}
		else if (option[0] == "SAMPLE_COUNT")
		{
			istringstream(option[1]) >> NUM_SAMPLES;
		}
		else if (option[0] == "CG_EPSILON")
		{
			istringstream(option[1]) >> CG_EPSILON;
		}
		else if (option[0] == "CG_HISTORY_SIZE")
		{
			istringstream(option[1]) >> CG_AVERAGE_SIZE;
		}
		else if (option[0] == "SAVE_PLOTS")
		{
			SAVE_PLOTS = option[1] == "TRUE";
		}
		else if (option[0] == "SAVE_SPREADSHEET")
		{
			SAVE_SPREADSHEET = option[1] == "TRUE";
		}
		else if (option[0] == "SAVE_LOG_FILES")
		{
			SAVE_LOG_FILES = option[1] == "TRUE";
		}
	}
	//
	while (!prog.eof())
	{
		string curLine;
		getline(prog, curLine);
		if (stringIsWhitespace(curLine) || curLine[0] == '#')
			continue;
		string curF;
		int curOpt;
		istringstream(curLine) >> curF >> curOpt;
		// Check that we can open the file now; this avoid a crash later if the program is invalid
		string filePath = programWorkingDirectory + curF;
		ifstream fileCheck;
		fileCheck.open(filePath);
		if (!fileCheck.is_open())
		{
			cout << "Error: Can not find the following file:" << endl;
			cout << curF << endl;
			return false;
		}
		fileCheck.close();
		fileList.push_back(filePath);
		optimals.push_back(curOpt);
	}
	prog.close();
	if (fileList.size() == 0)
	{
		cout << "Error: Could not parse program file. Exiting." << endl;
		return false;
	}
	return true;
}

bool programArguments::parseFileList(int argc, char *argv[])
{
	// Doesn't catch all malformed input errors, but makes sure the loop below
	// won't try to access garbage data
	if (!(argc & 1))
	{
		cout << "Error: Couldn't parse arguments. Each file must be paired with an optimal value." << endl;
		return false;
	}

	for (int i = 1; i < argc; i++)
	{
		char *file = argv[i++];
		char *opt = argv[i];
		fileList.push_back(string(file));
		optimals.push_back(atoi(opt));
	}
	return true;
}

void start(string method, ofstream &f, ofstream &d)
{
	f << endl;
	f << "Method: " << method << endl;
	d << endl;
	d << "Method: " << method << endl;
	cout << "Starting " << method << "..." << endl;
}

void end(string method, ofstream &f, ofstream &d, ofstream &csv, unsigned long long cycleCount, double userTime, int best)
{
	csv << cycleCount << "," << userTime << ",";
	f << "Total Cycle time: " << cycleCount << endl;
	f << "End Method: " << method << endl;
	f << endl;
	d << "End Method: " << method << endl;
	d << endl << endl;
	cout << "Done with " << method << "." << endl;
	cout << "Result: " << best << endl;
	cout << "Total user time was " << userTime << " seconds." << endl << endl;
}

void writeSpreadsheetHeader(ofstream &csv, programArguments *prog)
{
	csv << "Instance,optimal,";
	vector<string> methods;
	if (prog->USE_CFDS)
		methods.push_back("Christofides");
	if (prog->USE_CG)
		methods.push_back("Column Generation");
	if (prog->USE_CG && prog->USE_SR)
		methods.push_back("Column Generation + swapRound");
	if (prog->USE_ME)
		methods.push_back("Max Entropy");
	if (prog->USE_ES)
		methods.push_back("Splitting Off");
	if (prog->USE_ES && prog->USE_SR)
		methods.push_back("Splitting Off + swapRound");
	for (unsigned i = 0; i < methods.size(); i++)
	{
		csv << methods[i] << ",";
		csv << "Number of Trees,";
		csv << "Average tree cost,";
		csv << "Average matching cost,";
		csv << "Average total cost,";
		csv << "Average tour cost,";
		csv << "Spread,";
		csv << "Tree Cost Standard Deviation,";
		csv << "Matching Cost Standard Deviation,";
		csv << "Total Cost Standard Deviation,";
		csv << "Tour Cost Standard Deviation,";
		csv << "Worst tour,";
		csv << "Best tour,";
		csv << "Best Error,";
		csv << "Worst Error,";
		csv << "Average Error,";
		csv << "Average matching edge count,";
		csv << "Average matching edge cost,";
		csv << "1,2,3,4,5,6,>6,";
		csv << ">6 Mean Degree,";
		csv << ">6 Degree Standard Deviation,";
		csv << ">6 Max Degree,";
		csv << "Cycle Count,";
		csv << "User Time,";
	}
	csv << endl;
}


int run(Graph *g, vector<vector<Edge>> &trees, int optimal, ofstream &f, ofstream &d, ofstream &csv, ofstream &histogram,
	string method, string *itrPlot, vector<double> betas, int hangoutTime)
{
	int best = INT_MAX;
	int worst = 0;
	int numErrors = 0;
	vector<double> tourCosts;
	vector<double> treeCosts;
	vector<double> matchingCosts;
	vector<int> numEdgesInMatchings;
	vector<double> largeBinVals;
	double distrBins[7] = { 0.0 };
	f << "Number of Trees: " << trees.size() << endl;
	vector<double> usedBetas;
	for (unsigned i = 0; i < trees.size(); i++)
	{
		d << endl;
		d << "Tree #" << i << endl;
		if (trees[i].size() == 0)
		{
			f << "Tree #" << i << " failed in generation" << endl;
			d << "FAILED" << endl;
			continue;
		}
		vector<int> tour;
		int result = 0;
		if (g->costMat.size() == 0)
		{
#ifdef DOUBLE_SUPPORT
			std::future<int> future = std::async(std::launch::async, [&i, &g, &trees, &result, &tour, &numErrors, &f]()
			{
				try {
					result = christofides(g, trees[i], tour, true);
					//cout << "Found a tour for tree #" << i << ": Length is " << result << endl;
					return result;
				}
				catch (...) {
					f << "Tree #" << i << " failed in christofides!" << endl;
					cout << "Tree #" << i << " failed in christofides!" << endl;
					numErrors++;
					return 0;
				}
			});
			std::future_status status;
			status = future.wait_for(std::chrono::seconds(hangoutTime));
			if (status == std::future_status::timeout)
			{
				cout << "Warning: BlossomV hung on tree #" << i << endl;
				f << "Warning: BlossomV hung on tree #" << i << endl;
				d << "HUNG" << endl;
				numErrors++;
				continue;
			}
#else
			try {
				result = christofides(g, trees[i], tour, true);
			}
			// This should NOT ever be happening; numErrors will be removed in an upcoming commit
			// It's here because there have been errors in the past with multithreading support leading to crashes...
			catch (...) {
				f << "Tree #" << i << " failed in christofides!" << endl;
				cout << "Tree #" << i << " failed in christofides!" << endl;
				numErrors++;
			}
#endif
		}
		else
		{
			result = christofidesMat(g, trees[i], tour, true);
		}
		//cout << "Out of async: result is: " << result << endl;
		if (result == 0)
		{
			cout << "Warning: Result 0 outside async, continuing" << endl;
			continue;
		}
		tourCosts.push_back(result);
		best = min(best, result);
		worst = max(worst, result);
		treeInfo ti = getTreeInfo(trees[i], g);
		treeCosts.push_back(ti.treeCost);
		matchingCosts.push_back(ti.matchingCost);
		numEdgesInMatchings.push_back(ti.numEdgesInMatching);
		d << "  Tree Cost: " << ti.treeCost << endl;
		d << "  Matching Cost: " << ti.matchingCost << endl;
		d << "  Total Cost: " << ti.treeCost + ti.matchingCost << endl;
		d << "  Final Tour Length: " << result << endl;
		d << "  Mean Degree: " << ti.mean << endl;
		d << "  Mode Degree: " << ti.mode << endl;
		d << "  Degree Distribution:" << endl;
		for (unsigned j = 0; j < ti.degreeDist.size(); j++)
		{
			d << "    " << ti.degreeDist[j].first << "  :  " << ti.degreeDist[j].second << endl;
			distrBins[min(ti.degreeDist[j].first - 1, 6)] += ti.degreeDist[j].second;
			if (ti.degreeDist[j].first > 6)
				largeBinVals.push_back(ti.degreeDist[j].first);
		}
		*itrPlot += to_string(i + 1) + " " + to_string(best) + "\n";
		if (betas.size() != 0)
			usedBetas.push_back(betas[i]);
	}
	if (numErrors != 0)
		f << "Warning: Tours couldn't be generated from " << numErrors << " of the trees!" << endl;
	csv << " ," << trees.size() - numErrors << ",";
	f << "Average tree cost: " << wAvg(treeCosts, usedBetas) << endl;
	csv << wAvg(treeCosts, usedBetas) << ",";
	f << "Average matching cost: " << wAvg(matchingCosts, usedBetas) << endl;
	csv << wAvg(matchingCosts, usedBetas) << ",";
	f << "Average total cost: " << wAvg(elementWiseSum(treeCosts, matchingCosts), usedBetas) << endl;
	csv << wAvg(elementWiseSum(treeCosts, matchingCosts), usedBetas) << ",";
	f << "Average tour cost: " << wAvg(tourCosts, usedBetas) << endl;
	csv << wAvg(tourCosts, usedBetas) << ",";
	f << "Spread: " << worst - best << endl;
	csv << worst - best << ",";
	f << "Tree Cost Standard Deviation: " << wStdDev(treeCosts, usedBetas) << endl;
	csv << wStdDev(treeCosts, usedBetas) << ",";
	f << "Matching Cost Standard Deviation: " << wStdDev(matchingCosts, usedBetas) << endl;
	csv << wStdDev(matchingCosts, usedBetas) << ",";
	f << "Total Cost Standard Deviation: " << wStdDev(elementWiseSum(treeCosts, matchingCosts), usedBetas) << endl;
	csv << wStdDev(elementWiseSum(treeCosts, matchingCosts), usedBetas) << ",";
	f << "Tour Cost Standard Deviation: " << stdDev(tourCosts) << endl;
	csv << wStdDev(tourCosts, usedBetas) << ",";
	f << "Worst tour: " << worst << endl;
	csv << worst << ",";
	f << "Best tour: " << best << endl;
	csv << best << ",";
	f << "Best Error: " << ((double)(best - optimal)) / ((double)optimal)  * 100.0 << "%" << endl;
	csv << ((double)(best - optimal)) / ((double)optimal)  * 100.0 << "%,";
	f << "Worst Error: " << ((double)(worst - optimal)) / ((double)optimal)  * 100.0 << "%" << endl;
	csv << ((double)(worst - optimal)) / ((double)optimal)  * 100.0 << "%,";
	f << "Average Error: " << ((double)(wAvg(tourCosts, usedBetas) - optimal)) / ((double)optimal)  * 100.0 << "%" << endl;
	csv << ((double)(wAvg(tourCosts, usedBetas) - optimal)) / ((double)optimal)  * 100.0 << "%,";
	//New CSV Contents:
	csv << wAvg(numEdgesInMatchings, usedBetas) << ","; //Avg # edges in a matching
	csv << wSum(matchingCosts, usedBetas) / wSum(numEdgesInMatchings, usedBetas) << ","; //Avg cost of an edge in a matching
	int binSum = 0;
	for (int i = 0; i < 7; i++)
		binSum += (int)distrBins[i];

	for (int i = 0; i < 7; i++)
	{
		csv << ((float)distrBins[i]) / ((float)binSum) << ",";
	}
	//
	csv << sum(largeBinVals) / max(largeBinVals.size(), 1) << ","; // >6 Mean Node Degree
	csv << stdDev(largeBinVals) << ","; // >6 Std. Dev. Node Degree
	csv << lMax(largeBinVals) << ","; // >6 Max Node Degree
	histogram << method << ",";
	for (unsigned i = 0; i < tourCosts.size(); i++)
		histogram << tourCosts[i] << ",";
	histogram << endl;
	return best;
}

int getBestTree(vector<vector<Edge>> trees, Graph *g)
{
	//int best = std::numeric_limits<int>::max();
	int best = 1000000000;
	int best_loc = -1;
	for (unsigned i = 0; i < trees.size(); i++)
	{
		if (trees[i].size() == 0)
			continue;
		vector<int> tour;
		int result = christofides(g, trees[i], tour, true);
		cout << "Result is: " << result << endl;
		if (result < best)
		{
			best = result;
			best_loc = i;
		}
	}
	cout << "Chose " << best_loc << " with cost of " << best << endl;
	return best_loc;
}

void visualizeAll(string fName, int numSamples, double colGenEps, int colGenAvg)
{
	cout << "Entered visualizer" << endl;
	Graph g(fName);
	//Held-Karp Solution:
	cout << "Solving Held-Karp" << endl;
	vector<LPEdge> hkSoln = getHKSolutions(fName, true);
	//Christofides
	cout << "Running Standard Christofides" << endl;
	vector<Edge> cfdsTree = minSpanTree(&g);

	//Init sampled trees
	vector<Edge> cgTree;
	vector<Edge> cgSrTree;
	vector<Edge> splitTree;
	vector<Edge> splitSrTree;
	vector<Edge> maxEntropyTree;

	vector<vector<Edge>> cgTrees;
	vector<vector<Edge>> cgSrTrees;
	vector<vector<Edge>> meTrees;
	vector<vector<Edge>> splitTrees;
	vector<vector<Edge>> splitSrTrees;

	//CG, CG+SR
	cout << "Starting Column Generation" << endl;
	{
		vector<vector<LPEdge>> trees;
		vector<double> betas;
		{
			int num_its = 0;
			string plot = "";
			pair<vector<vector<LPEdge>>, vector<double>> p = bestOfMany(&g, hkSoln, colGenEps, colGenAvg, &num_its, &plot);
			cout << "It took " << num_its << " iterations to solve Column Generation" << endl;
			trees = p.first;
			betas = p.second;

			for (unsigned j = 0; j < trees.size(); j++)
				cgTrees.push_back(lpToEdge(trees[j]));
			cgTree = cgTrees[getBestTree(cgTrees, &g)];
		}

		cout << "Starting swapRound for Column Generation" << endl;
		for (int i = 0; i < numSamples; i++)
		{
			cout << "Sampling Tree #" << i << " for swapRound" << endl;
			cgSrTrees.push_back(lpToEdge(swapRound(trees, betas).back()));
		}
		cgSrTree = cgSrTrees[getBestTree(cgSrTrees, &g)];
	}


	//Max Entropy
	cout << "Starting Max Entropy" << endl;
	vector<ModEdge> graph = xToZ(lpToMod(hkSoln));
	vector<double> lambda = getLambdaFast(graph);
	weightGraph(graph, lambda);
	vector<vector<Edge>> trees;
	for (int i = 0; i < numSamples; i++)
	{
		cout << "Sampling Tree #" << i << " for max entropy" << endl;
		meTrees.push_back(lpToEdge(randomWalk(graph, true)));
	}
	maxEntropyTree = meTrees[getBestTree(meTrees, &g)];

	//Splitting Off
	cout << "Starting Splitting Off..." << endl;
	int k = 0;
	{
		vector<ModEdge> gMod = lpToMod(hkSoln);
		pair<int, vector<SplitEdge>> kOut = getKLinear(gMod);
		k = kOut.first / 2;
		vector<SplitEdge> gS = kOut.second;
		historyIndex h = createLog();
		vector<SplitEdge> finalRes = edgeSplit(gS, k * 2, h);
		splitting ss = logToSplits(h);
		pair<vector<vector<JoinEdge>>, vector<JoinEdge>> ret = createTrees(finalRes, ss, k);
		vector <vector<JoinEdge>> trees = ret.first;
		vector<JoinEdge> extras = ret.second;
		joinEdges(ss, k, extras, trees);
		splitTrees = joinEForestToEForest(trees);
		splitTree = splitTrees[getBestTree(splitTrees, &g)];
	}

	cout << "Starting Splitting Off + swapRound..." << endl;
	//Splitting Off + swapRound
	{
		vector<vector<LPEdge>> treesLP = eToLP(splitTrees);
		vector<double> betas(numSamples, 1.0 / (double)k);
		for (int i = 0; i < numSamples; i++)
		{
			cout << "Sampling Tree #" << i << " for swapRound" << endl;
			splitSrTrees.push_back(lpToEdge(swapRound(treesLP, betas).back()));
		}
		splitSrTree = splitSrTrees[getBestTree(splitSrTrees, &g)];

	}

	cout << "Drawing Christofides tree: " << endl;
	g.drawTree(cfdsTree);
	cout << "Drawing Max Entropy tree: " << endl;
	g.drawTree(maxEntropyTree);
	cout << "Drawing Col. Gen. tree: " << endl;
	g.drawTree(cgTree);
	cout << "Drawing Col. Gen + SR tree: " << endl;
	g.drawTree(cgSrTree);
	cout << "Drawing Splitting Off tree: " << endl;
	g.drawTree(splitTree);
	cout << "Drawing Splitting Off + SR tree: " << endl;
	g.drawTree(splitSrTree);
	return;
}