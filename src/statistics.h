#ifndef STATS_H
#define STATS_H

#include "graph.h"
double stdDev(const vector<double> &l);
double stdDev(const vector<int> &l);
double wAvg(const vector<int> &l, const vector<double> &weights);
double wAvg(const vector<double> &l, const vector<double> &weights);
double wStdDev(const vector<double> &l, const vector<double> &weights);
double wStdDev(const vector<int> &l, const vector<double> &weights);

treeInfo getTreeInfo(const vector<Edge> &t, Graph *g);
void checkOddsOfTreeSet(const vector<ModEdge> &lGraph, const vector<vector<Edge>> &ts);
pair<double, double> checkNegativeCorrelation(const vector<ModEdge> &lGraph, const vector<vector<Edge>> &ts);
void checkConvComb(const vector<LPEdge> &hkSoln, int k, const pair<vector<vector<JoinEdge>>, vector<JoinEdge>> &ts);
void output(treeInfo i);

#endif