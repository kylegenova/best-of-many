/// Kyle Genova
/// 6/16/14
/// Header for the Graph class

#ifndef GRAPH_H
#define GRAPH_H


#include "definitions.h"

class Graph
{
public:
	int numNodes = 0;
	Graph() { numNodes = 0; };

	vector<Point> nodes;
	vector<vector<double>> costMat;
	vector<Edge> getESet();
	Graph(string fName, bool silent=false);
	Graph(vector<Point> &nodeList);
	int distance(int from, int to);
	int getNumNodes();
	vector<Edge> runMST();
	vector<Edge> runDelaunayMST(vector<Edge> &in);
	//this next one probably should be in a sparse graph or Util module?
	vector<int> degrees(vector<Edge> &edges);
	vector<int> odds(vector<int> &degrees);
	vector<Edge> loadSolver(vector<int> &odds, int *numAdded, PerfectMatching *pm);
	void loadGeomSolver(vector<int> &odds, GeomPerfectMatching *gpm);
	vector<Edge> removeDups(vector<Edge> in);
	vector<Edge> createMultigraph(vector<Edge> &e1, vector<Edge> &e2);
	vector<int> createCircuit(vector<Edge> &mg, Graph *g);
	vector<int> shortcutCircuit(vector<int> &circ);
	vector<int> shortcutCircuitSelective(vector<int> &circ);
	void writeCircuit(vector<int> toWrite, string fName);
	int tourCost(vector<int> &tour);
	int edgeSetCost(const vector<Edge> &t);
	void drawTree(const vector<Edge> &t);
	double getTreeCost(const vector<LPEdge> &t);
};
#endif