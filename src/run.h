#include "graph.h"

void start(string method, ofstream &f, ofstream &d);
void end(string method, ofstream &f, ofstream &d, ofstream &csv, unsigned long long cycleCount, double userTime, int best);
void writeSpreadsheetHeader(ofstream &csv, programArguments *prog);

int run(Graph *g, vector<vector<Edge>> &trees, int optimal, ofstream &f, ofstream &d, ofstream &csv, ofstream &histogram,
	string method, string *itrPlot, vector<double> betas, int hangoutTime);
int getBestTree(vector<vector<Edge>> trees, Graph *g);

void visualizeAll(string fName, int numSamples, double colGenEps, int colGenAvg);
