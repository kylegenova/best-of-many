#include "christofides.h"
#include "statistics.h"
#include "edge_split.h"
#include "edge_join.h"
#include "graph.h"
#include "lrs.h"
#include "gamma.h"
#include "heldkarp.h"
#include "random_walk.h"
#include "swapround.h"
#include "col_gen.h"
#include "run.h"

/* All Windows specific definitions */
#ifdef _WIN32
#include <Windows.h>
unsigned long long startClock()
{

	HANDLE hdl = GetCurrentProcess();
	unsigned long long st;
	QueryProcessCycleTime(hdl, &st);
	return st;
}

unsigned long long endClock(unsigned long long st) {
	HANDLE hdl = GetCurrentProcess();
	unsigned long long et;
	QueryProcessCycleTime(hdl, &et);
	return et - st;
}

ULARGE_INTEGER getUserTime()
{
	FILETIME lpCreation;
	FILETIME exitTime;
	FILETIME kernelTime;
	FILETIME userTime;
	GetProcessTimes(GetCurrentProcess(), &lpCreation, &exitTime, &kernelTime, &userTime);
	ULARGE_INTEGER uT;
	uT.HighPart = userTime.dwHighDateTime;
	uT.LowPart = userTime.dwLowDateTime;
	return uT;
}

double getUserTimeDifference(ULARGE_INTEGER t1, ULARGE_INTEGER t2)
{
	return ((double)(t2.QuadPart - t1.QuadPart)) / 1.0e7;
}

double getUserTimeFromStart(ULARGE_INTEGER start)
{
	return getUserTimeDifference(start, getUserTime());
}
void cleanupDirectory(programArguments *prog)
{
	system("DEL lp lpbackup.lp *.sol");
	if (!prog->SAVE_SPREADSHEET)
		system("DEL results.csv");
	if (!prog->SAVE_PLOTS)
		system("DEL *.hst *.meip *.sosrip *.cgp *.cgsrip");
	if (!prog->SAVE_LOG_FILES)
		system("DEL log.txt details.txt ncr.txt results.txt");
}
#endif

/* Defaults */

// The number of threads to use, if OpenMP compilation is turned on. By default it is turned OFF
// in the project, so this won't affect anything. OpenMP support can be enabled in the project settings
// or with the compiler flag /openmp. The flag /openmp- disables it.
#define NUM_THREADS 4

// The program can be compiled to support either integers or integers and doubles.
// However, if double-support is enabled, blossomV may hang, even on integral instances.
// To changed the type to doubles, uncomment the definition below:
// #define DOUBLE_SUPPORT

#ifdef DOUBLE_SUPPORT
// Number of seconds to wait for a Christofides result 
// before declaring that blossomV hung on that instance
int HANGOUT_TIME_SEC = 1000;
#else 
int HANGOUT_TIME_SEC = 0; // Just a placeholder for consistent function declarations
#endif

int main(int argc, char* argv[])
{
	if (argc == 1)
	{
		cout << "Error: No input detected." << endl;
		cout << "Exiting." << endl;
		exit(1);
	}
	programArguments prog;
	if (argc == 2)
	{
		if (!prog.parseProgram(argv[1]))
		{
			cout << "Error: Couldn't parse program. Exiting." << endl;
			exit(1);
		}
	}
	else if (!prog.parseFileList(argc, argv))
	{
		cout << "Error: Couldn't parse file list. Exiting." << endl;
		exit(1);
	}
	
	// Open all the output files
	ofstream f;
	f.open("results.txt", std::ofstream::trunc);
	ofstream d;
	d.open("details.txt", std::ofstream::trunc);
	ofstream csv;
	csv.open("results.csv", std::ofstream::trunc);
	writeSpreadsheetHeader(csv, &prog);

	f << "Beginning execution on " << prog.fileList.size() << " file" << ((prog.fileList.size() > 1) ? "s" : "") << "." << endl;

	int NUM_SAMPLES = prog.NUM_SAMPLES; // Must pull out of program for OpenMP
	for (unsigned i = 0; i < prog.fileList.size(); i++)
	{
		string fName = prog.fileList[i];
		ofstream histogram;
		
		string fNoDir = getFileWithoutDirectory(fName);
		// Won't crash if < 4 characters; already checked when parsing program
		fNoDir = fNoDir.substr(0, fNoDir.size() - 4);

		histogram.open((fNoDir + ".hst"), std::ofstream::trunc);
		// This can't be done without the visualizer; that code is not currently in the repository.
		// It may be added once a good opt-out is implemented -it requires the DirectX SDK,
		// which is a 500MB download, and the visualizer itself it unlikely to be used frequently.
		if (prog.VISUALIZE)
		{
			visualizeAll(fName, prog.NUM_SAMPLES, prog.CG_EPSILON, prog.CG_AVERAGE_SIZE);
			continue;
		}
		Graph g(fName);

		
		f << endl << "Instance: " << fName << endl << endl;
		csv << fName << "," << prog.optimals[i] << ",";
		f << "Optimal (Provided): " << prog.optimals[i] << endl;

		// Held-Karp Solution:
		double lpCost = g.getTreeCost(getHKSolutions(fName, true));
		f << "LP Cost: " << lpCost << endl;
		cout << "LP Cost: " << lpCost << endl;
		
		// Christofides
		if (prog.USE_CFDS)
		{
			start("Christofides", f, d);
			unsigned long long st = startClock();
			auto userStart = getUserTime();
			vector<int> tour;
			vector<vector<Edge>> trees;
			trees.push_back(minSpanTree(&g));
			string itrPlot;
			int bestTour = run(&g, trees, prog.optimals[i], f, d, csv, histogram, "Christofides", &itrPlot, vector<double>(), HANGOUT_TIME_SEC);
			end("Christofides", f, d, csv, endClock(st), getUserTimeFromStart(userStart), bestTour);
		}
		
		// Column Generation
		if (prog.USE_CG)
		{
			vector<vector<LPEdge>> trees;
			vector<double> betas;
			vector<LPEdge> lp;
			// Column Generation
			{
				start("Column Generation", f, d);
				unsigned long long st = startClock();
				auto userStart = getUserTime();
				lp = getHKSolutions(fName, true);
				int num_its = 0;
				string *plot = new string();
				pair<vector<vector<LPEdge>>, vector<double>> p = bestOfMany(&g, lp, prog.CG_EPSILON, prog.CG_AVERAGE_SIZE, &num_its, plot);
				ofstream cgHist;
				cgHist.open(fNoDir + ".cgp", std::ofstream::trunc);
				cgHist << *plot;
				cgHist.close();
				delete plot;
				trees = p.first;
				betas = p.second;
				vector<vector<Edge>> ts;
				string itrPlot;
				for (unsigned j = 0; j < trees.size(); j++)
					ts.push_back(lpToEdge(trees[j]));
				int bestTour = run(&g, ts, prog.optimals[i], f, d, csv, histogram, "Column Generation", &itrPlot, betas, HANGOUT_TIME_SEC);
				end("Column Generation", f, d, csv, endClock(st), getUserTimeFromStart(userStart), bestTour);
			}

			//Column Generation + swapRound
			if (prog.USE_SR)
			{
				start("Column Generation + swapRound", f, d);
				unsigned long long st = startClock();
				auto userStart = getUserTime();
				cout << "Sampling..." << endl;
				vector<vector<Edge>> ts;
#pragma omp parallel num_threads(NUM_THREADS) shared(trees, betas, ts, NUM_SAMPLES)
				{
#pragma omp for 
					for (int j = 0; j < NUM_SAMPLES; j++)
					{

						try 
						{ 

							vector<Edge> t = lpToEdge(swapRound(trees, betas).back());
#pragma omp critical
							{
								ts.push_back(t);
							}
#pragma omp flush(trees)
						}
						catch (...) 
						{ 
#pragma omp critical
							{
								ts.push_back(vector<Edge>());
							}
#pragma omp flush(trees)
						}
					}
				}
				//checkOddsOfTreeSet(lpToMod(lp), ts);
				// Uncomment this for a verification of the negative correlation property of swapRound
				//checkNegativeCorrelation(lpToMod(lp), ts);
				string itrPlot;
				int bestTour = run(&g, ts, prog.optimals[i], f, d, csv, histogram, "Column Generation + swapRound", &itrPlot, vector<double>(), HANGOUT_TIME_SEC);
				end("Column Generation + swapRound", f, d, csv, endClock(st), getUserTimeFromStart(userStart), bestTour);
				ofstream cgsrImpPlot;
				cgsrImpPlot.open(fNoDir + ".cgsrip", std::ofstream::trunc);
				cgsrImpPlot << itrPlot;
				cgsrImpPlot.close();
			}
		}
		
		//Max Entropy
		if (prog.USE_ME)
		{
			start("Max Entropy", f, d);
			unsigned long long st = startClock();
			auto userStart = getUserTime();
			f << "Number of Sampled Trees: " << prog.NUM_SAMPLES << endl;
			cout << "Finding Held-Karp Solution" << endl;
			vector<LPEdge> hkSoln = getHKSolutions(fName, false);
			cout << "Done finding Held-Karp Solution" << endl;
			vector<ModEdge> graph = xToZ(lpToMod(hkSoln));
			cout << "Generating lambdas..." << endl;
			//vector<double> lambda = getLambda(graph);
			vector<double> lambda = getLambdaFast(graph);
			cout << "Done generating lambdas." << endl;
			cout << "Sampling..." << endl;
			weightGraph(graph, lambda);

			vector<vector<Edge>> trees;
#pragma omp parallel num_threads(NUM_THREADS) shared(graph, NUM_SAMPLES, trees)
			{
#pragma omp for 
				for (int j = 0; j < NUM_SAMPLES; j++)
				{
					try {
						vector<Edge> t = lpToEdge(randomWalk(graph, true));
#pragma omp critical 
						{
							trees.push_back(t);
							//cout << "Done sampling tree #" << j << endl;
						}
#pragma omp flush(trees)

					}
					catch (...) {
#pragma omp critical
						{
							trees.push_back(vector<Edge>());
						}
#pragma omp flush(trees)
					}
				}
			}
			string itrPlot;
			int bestTour = run(&g, trees, prog.optimals[i], f, d, csv, histogram, "Max Entropy", &itrPlot, vector<double>(), HANGOUT_TIME_SEC);
			end("Max Entropy", f, d, csv, endClock(st), getUserTimeFromStart(userStart), bestTour);
			ofstream meImpPlot;
			meImpPlot.open(fNoDir + ".meip", std::ofstream::trunc);
			meImpPlot << itrPlot;
			meImpPlot.close();
		}
		
		//Splitting Off
		if (prog.USE_ES)
		{
			vector<vector<Edge>> ts;
			vector<LPEdge> lp;
			int k;
			//Splitting Off
			{
				start("Splitting Off", f, d);
				unsigned long long st = startClock();
				auto userStart = getUserTime();
				lp = getHKSolutions(fName, true);
				vector<ModEdge> gMod = lpToMod(lp);
				pair<int, vector<SplitEdge>> kOut = getKLinear(gMod);
				k = kOut.first / 2;
				f << "K: " << k << endl;
				vector<SplitEdge> gS = kOut.second;
				historyIndex h = createLog();
				vector<SplitEdge> finalRes = edgeSplit(gS, k * 2, h);
				splitting ss = logToSplits(h);
				pair<vector<vector<JoinEdge>>, vector<JoinEdge>> ret = createTrees(finalRes, ss, k);
				vector <vector<JoinEdge>> trees = ret.first;
				vector<JoinEdge> extras = ret.second;
				joinEdges(ss, k, extras, trees);
				vector<vector<Edge>> cfdsTrees = joinEForestToEForest(trees);
				ts = cfdsTrees;
				checkConvComb(lp, k, pair<vector<vector<JoinEdge>>, vector<JoinEdge>>(trees, extras));
				string itrPlot;
				int bestTour = run(&g, cfdsTrees, prog.optimals[i], f, d, csv, histogram, "Splitting Off", &itrPlot, vector<double>(), HANGOUT_TIME_SEC);
				end("Splitting Off", f, d, csv, endClock(st), getUserTimeFromStart(userStart), bestTour);
			}

			//Splitting Off + swapRound
			if (prog.USE_SR)
			{
				start("Splitting Off + swapRound", f, d);
				ofstream ncf;
				ncf.open("ncr.txt", std::ofstream::app);
				unsigned long long st = startClock();
				auto userStart = getUserTime();
				cout << "Sampling..." << endl;
				vector<vector<LPEdge>> treesLP = eToLP(ts);
				vector<double> betas(ts.size(), 1.0 / (double)k);
				vector<vector<Edge>> trees;
#pragma omp parallel num_threads(NUM_THREADS) shared(NUM_SAMPLES, treesLP, betas, trees)
				{
#pragma omp for
					for (int j = 0; j < NUM_SAMPLES; j++)
					{
						try {
							vector<Edge> t = lpToEdge(swapRound(treesLP, betas).back());
#pragma omp critical
							{
								trees.push_back(t);
							}
#pragma omp flush(trees)
						}
						catch (...) {
#pragma omp critical 
							{
								trees.push_back(vector<Edge>());
							}
#pragma omp flush(trees)
						}
					}
				}
				//checkOddsOfTreeSet(lpToMod(lp), trees);
				// Uncomment this to get a verification of the negative correlation property of swapRound
				//pair<double,double> ncr = checkNegativeCorrelation(lpToMod(lp), trees);
				string itrPlot;
				int bestTour = run(&g, trees, prog.optimals[i], f, d, csv, histogram, "Splitting Off + swapRound", &itrPlot, vector<double>(),HANGOUT_TIME_SEC);
				end("Splitting Off + swapRound", f, d, csv, endClock(st), getUserTimeFromStart(userStart), bestTour);
				ofstream sosrImpPlot;
				sosrImpPlot.open(fNoDir + ".sosrip", std::ofstream::trunc);
				sosrImpPlot << itrPlot;
				sosrImpPlot.close();
			}
		}
		f << "Done with instance " << fName << endl << endl;
		csv << endl;
		histogram.close();
	}
	f.close();
	d.close();
	csv.close();
	cleanupDirectory(&prog);
	return 0;
}