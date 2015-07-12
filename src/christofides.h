#ifndef CFDS_H
#define CFDS_H

#include "graph.h"

/// [christofides g tour] takes in a Graph g and a reference to an empty vector<int> where the answer
///   will be written.
/// g: The graph on which to run the christofides algorithm
/// tour: an empty vector<int> where the result of the algorithm will be stored.
/// Returns: the size of the tour found. 0 if an error occured.
int christofides(Graph *g, vector<int> &tour, bool silent=false);

/// [minSpanTree g] takes in a Graph g
/// g: The graph on which to run the spanning tree algorithm
/// Returns: The a minimum spanning tree using a delaunay triangulation.
vector<Edge> minSpanTree(Graph *g);

/// [christofides g spanTree tour] takes in a Graph g, a vector<Edge> spanTree, and a reference to an empty 
///   vector<int> where the answer will be written.
/// g: The graph on which to run the christofides algorithm
/// spanTree: The tree to use instead of finding the default minimum spanning tree.
/// tour: an empty vector<int> where the result of the algorithm will be stored.
/// Returns: the size of the tour found. 0 If an error occured.
int christofides(Graph *g, vector<Edge> &spanTree, vector<int> &tour, bool silent=false);


/// Runs christofides on fName and writes the results to outFName.
/// silent is a verbosity flag.
int writeChristofides(string fName, string outFName, bool silent=false);

/// A wrapper for christofides
int christofides(string fName, vector<int> &tour, bool silent=false);

/// [matchingCost g spanTree numEIM] takes in a Graph g, and a Spanning Tree spanTree within
///   g, and an int ref, and returns the cost of a minimum perfect matching for the spanning tree in
///   the graph. 
/// g: The graph used for distance calculations
/// spanTree: The spanning tree which should have a minimum perfect matching.
/// numEIM: a refernce to an int where the number of edges in the matching will be saved.
/// Returns: the cost of a MPM according to the weights in g.
double matchingCost(Graph *g, vector<Edge> &spanTree, int &numEIM);


// [christofidesMat g tour] takes in a Graph g and a reference to an empty vector<int> where the answer
//   will be written.
// g: The graph on which to run the christofides algorithm
// tour: an empty vector<int> where the result of the algorithm will be stored.
// Returns: The size of the tour found. 0 if an error occured.
int christofidesMat(Graph *g, vector<int> &tour, bool silent);

// [christofidesMat g spanTree tour] takes in a Graph g, a vector<Edge> spanTree, and a reference to an empty 
//   vector<int> where the answer will be written.
// g: The graph on which to run the christofides algorithm
// spanTree: The tree to use instead of finding the default minimum spanning tree.
// tour: an empty vector<int> where the result of the algorithm will be stored.
// silent: if silent, no print statements will output
// Returns: the size of the tour found. 0 If an error occured.
int christofidesMat(Graph *g, vector<Edge> &spanTree, vector<int> &tour, bool silent);

#endif