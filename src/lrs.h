#ifndef LRS_H
#define LRS_H

#include "definitions.h"

vector<LPEdge> sampleLRTNew(vector<ModEdge> lGraph, bool silent);
vector<LPEdge> sampleLRTBig(vector<ModEdge> lGraph, bool silent);
vector<LPEdge> sampleLRTInverse(vector<ModEdge> lGraph, bool silent);
int getLRSTour(string fName, vector<double> gammaIn, int numNodes, bool silent = false);
int runLRS(string fName, int numSamples = 10, bool silent = false);
vector<double> getLambdaBaked(string fName, const vector<ModEdge> &g);
int runLRSBaked(string fName, int numSamples);
#endif