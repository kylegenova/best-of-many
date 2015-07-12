#ifndef EDGEJOIN_H
#define EDGEJOIN_H


#include "definitions.h"



/// Takes a history and generates a splitting from it
splitting logToSplits(historyIndex &h);


pair<vector<vector<JoinEdge>>, vector<JoinEdge>> createTrees(const vector<SplitEdge> &g, const splitting &ss, int k);
void rejoinNode(splitting &ss, int k, vector<JoinEdge> &extras, vector<vector<JoinEdge>> &trees);
void joinEdges(splitting &ss, int k, vector<JoinEdge> &extras, vector<vector<JoinEdge>> &trees);
vector<Edge> joinEtoE(const vector<JoinEdge> &t);
vector<vector<Edge>> joinEForestToEForest(const vector <vector<JoinEdge>> &ts);
#endif