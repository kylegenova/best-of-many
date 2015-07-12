#ifndef EDGESPLIT_H
#define EDGESPLIT_H

#include "graph.h"
#include "lrs.h"
#include "util.h"



pair<unsigned int, vector<SplitEdge>> getK(const vector<ModEdge> &g);
pair<unsigned, vector<SplitEdge>> getKLinear(const vector<ModEdge> &g);
vector<SplitEdge> edgeSplit(const vector<SplitEdge> &g, int k, historyIndex &h);
vector<int> maOrdering(const vector<SplitEdge> &g, int s);
vector<int> maOrderingFast(const vector<SplitEdge> &g, int s);
historyIndex createLog();
void output(historyIndex &h);
void output(const history &h);

#endif