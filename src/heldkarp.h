//This is the interface to the interface to Concorde to retrieve held-karp solutions to the given .tsp file
#ifndef HK_H
#define HK_H

#include "definitions.h"

vector<LPEdge> getHKSolutions(string fName, bool silent = false);
vector<LPEdge> getHKSolutionsBaked(string fName, bool silent);

#endif