#ifndef COL_GEN_H
#define COL_GEN_H

#include "graph.h"



/// Run Column Generation on the graph [g], using the Held-Karp solution [lp]. Terminate if the algorithm 
/// has not improved the sum of the slacks by more than [epsilon] over the last [average_over_last] iterations.
/// Set [average_over_last] to 1 to run to optimality.
/// Fills [num_its_required] with the number of iterations executed.
/// Returns a pair containing the tree and the Yt values.
pair<vector<vector<LPEdge>>, vector<double>> bestOfMany(Graph *g, vector<LPEdge> &lp, double epsilon, int average_over_last, int *num_its_required, string *plot);

#endif