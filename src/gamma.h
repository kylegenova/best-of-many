#ifndef GAMMA_H
#define GAMMA_H

#include "definitions.h"
#include <math.h>

vector<double> getLambda(const vector<ModEdge> &zGraph);
vector<double> getLambdaFast(const vector<ModEdge> &zGraph);
#endif