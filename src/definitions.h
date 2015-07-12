#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <unordered_map>
#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <functional>
#include <fstream>
#include <chrono>
#include <string>
#include <sstream>
#include <queue>
#include <iostream>
#include <Windows.h>
#include <math.h>
#include <future>
#include "GeomPerfectMatching.h"

using std::function;
using std::vector;
using std::string;
using std::ostream;
using std::chrono::system_clock;
using std::unordered_map;
using std::mt19937;
using std::mersenne_twister_engine;
using std::uniform_real_distribution;
using std::cout;
using std::endl;
using std::ofstream;
using std::to_string;
using std::logic_error;
using std::runtime_error;
using std::invalid_argument;
using std::list;
using std::iterator;
using std::queue;
using std::ifstream;
using std::priority_queue;
using std::pair;
using std::istringstream;

/// The underlying point structure for the Graph Class. Represents a 2D point in space.
/// Switches between double and int representation with GeomPerfectMatching::REAL's header
/// definition.
struct Point {
public:
	int id;
	GeomPerfectMatching::REAL x, y;
	Point() { x = 0; y = 0; };
	Point(GeomPerfectMatching::REAL xCoord, GeomPerfectMatching::REAL yCoord) { x = xCoord; y = yCoord; };
};

/// An Unweighted edge.
struct Edge {
public:
	int u, v;
	Edge() { u = 0; v = 0; }
	Edge(int u, int v) { this->u = u; this->v = v; }
};

/// TODO: A Weighted edge with integer cost.
struct EdgeWI {
public:
	int u, v, w;
};

/// TODO: A weighted edge with double cost.
struct EdgeWD {
public:
	int u, v;
	double w;
};

struct PQEdge {
public:
	int ind1, ind2;
	int dist;
	PQEdge() { ind1 = 0; ind2 = 0; dist = 0; }
	PQEdge(int u, int v, int cost) {
		ind1 = u; ind2 = v; dist = cost;
	}
};

struct LPEdge
{
	int end0;
	int end1;
	double weight;
	LPEdge() { end0 = 0; end1 = 0; weight = 0.0; }
	LPEdge(int from, int to, double w) { end0 = from; end1 = to; weight = w; }
};

struct uEdge
{
	int u;
	int v;
	uEdge(int u, int v, bool sort = true) {
		if (u <= v || !sort)
		{
			this->u = u;
			this->v = v;
			return;
		}
		this->u = v;
		this->v = u;
	}
};
inline bool operator< (const uEdge &lhs, const uEdge &rhs) {
	return lhs.u < rhs.u || ((lhs.u == rhs.u) && lhs.v < rhs.v);
}
inline bool operator>(const uEdge &lhs, const uEdge &rhs) {
	return (rhs < lhs);
}
inline bool operator<= (const uEdge &lhs, const uEdge &rhs) {
	return !(lhs > rhs);
}
inline bool operator>= (const uEdge &lhs, const uEdge &rhs) {
	return !(lhs < rhs);
}
inline bool operator== (const uEdge &lhs, const uEdge &rhs) {
	return (lhs.u == rhs.u && lhs.v == rhs.v);
}
inline bool operator!= (const uEdge &lhs, const uEdge &rhs) {
	return !(lhs == rhs);
}

namespace std
{
	template <>
	struct hash<uEdge>
	{
		size_t operator()(const uEdge &e) const
		{
			return hash<int>()(e.u) ^ hash<int>()(e.v);
		}
	};
}

struct SplitEdge
{
	int end0;
	int end1;
	int weight;
	int orig0;
	int orig1;
	SplitEdge(int e0, int e1, int w, int was0, int was1) {
		end0 = e0; end1 = e1; weight = w; orig0 = was0; orig1 = was1;
	}
};

struct JoinEdge
{
	int u;
	int v;
	bool operator==(const JoinEdge &that)
	{
		return ((u == that.u && v == that.v) || (u == that.v && v == that.u));
	}
	JoinEdge(int end0, int end1) { u = end0; v = end1; }
};

/// This structure allows for information about the original graph to be preserved,
/// despite the graph being contracted.
struct ModEdge
{
	int end0;
	int end1;
	double weight;
	int orig0;
	int orig1;
	ModEdge(int e0, int e1, double w, int was0, int was1) {
		end0 = e0; end1 = e1; weight = w; orig0 = was0; orig1 = was1;
	}
};



struct history
{
	int s;
	int u;
	int v;
	int multiplicity;
	history *next;
	history *prev;
	history(int at, int from, int to, int weight)
	{
		s = at;
		u = from;
		v = to;
		multiplicity = weight;
		//cout << "Creating log item with s =" << s << ", u = " << u << ", v = " << v << ", weight = " << weight << endl;
		next = nullptr;
		prev = nullptr;
	}
	bool operator==(const history &that)
	{
		return (s == that.s && ((u == that.u && v == that.v) || (u == that.v && v == that.u)));
	}
	bool isFirst()
	{
		return (s == -1 && u == -1 && v == -1 && multiplicity == -1);
	}
};

struct historyIndex
{
	history *first;
	history *last;
	history *loc;
};

struct statistics
{
	int numTrials;
	double best;
	double worst;
	double mean;
	double median;
	double mode;
	double stdev;
};

struct treeInfo
{
	vector<pair<int, int>> degreeDist; //fst > Degree, snd > # occurences
	int numNodes;
	int numEdges;
	int numEdgesIncludeMatching;
	double mean;
	double mode;
	double matchingCost;
	double treeCost;
	int numEdgesInMatching;
	double avgEdgeCostInMatching;

};

struct programArguments
{
	string programWorkingDirectory;
	bool USE_CFDS; //Run Christofides
	bool USE_ME; //Run Max Entropy
	bool USE_CG; //Run Column Generation
	bool USE_ES; //Run Edge-splitting
	bool USE_SR; //Run swapRound where applicable
	bool VISUALIZE; //Don't generate data files, generate images.
	int NUM_SAMPLES; //Number of samples for Column Generation and swapRound
	double CG_EPSILON; //Early Termination Criterion for Column Generation
	int CG_AVERAGE_SIZE; //Number of previous iterations to average over
	bool SAVE_PLOTS;
	bool SAVE_SPREADSHEET;
	bool SAVE_LOG_FILES;
	vector<string> fileList;
	vector<int> optimals;
	programArguments() {
		programWorkingDirectory = "./"; // might need to be .
		USE_CFDS = true;
		USE_ME = true;
		USE_ES = true;
		USE_SR = true;
		USE_CG = true;
		VISUALIZE = false;
		NUM_SAMPLES = 100;
		CG_EPSILON = 0.1;
		CG_AVERAGE_SIZE = 100;
		fileList = vector<string>();
		optimals = vector<int>();
		SAVE_PLOTS = false;
		SAVE_SPREADSHEET = true;
		SAVE_LOG_FILES = false;
	}
	bool parseProgram(string programFilePath);
	bool parseFileList(int argc, char *argv[]);
};

typedef vector<pair<int, vector<JoinEdge>>> splitting;
typedef vector<vector<pair<int, int>>> AList;
typedef vector<vector<pair<int, double>>> AdjList;

#endif