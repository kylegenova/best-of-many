#include "statistics.h"
#include "util.h"
#include "christofides.h"

double stdDev(const vector<double> &l)
{
	double mean = sum(l) / (max(l.size(), 1));
	double sqSum = 0.0;
	for (unsigned i = 0; i < l.size(); i++)
		sqSum += pow(l[i] - mean, 2);
	return sqrt(sqSum / (max(mean, 1)));
}

double stdDev(const vector<int> &l)
{
	double mean = sum(l) / (max(l.size(), 1));
	double sqSum = 0.0;
	for (unsigned i = 0; i < l.size(); i++)
		sqSum += pow(l[i] - mean, 2);
	return sqrt(sqSum / (max(mean, 1)));
}

double wAvg(const vector<int> &l, const vector<double> &weights)
{
	if (weights.size() == 0)
		return sum(l) / l.size();
	return wSum(l, weights) / sum(weights);
}

double wAvg(const vector<double> &l, const vector<double> &weights)
{
	if (weights.size() == 0)
		return sum(l) / l.size();
	return wSum(l, weights) / sum(weights);
}

double wStdDev(const vector<double> &l, const vector<double> &weights)
{
	if (weights.size() == 0)
		return stdDev(l);
	else if (weights.size() == l.size())
	{
		double m = wAvg(l, weights);
		double acc = 0.0;
		for (unsigned int i = 0; i < l.size(); i++)
		{
			acc += weights[i] * (l[i] - m) * (l[i] - m);
		}
		return sqrt(acc / sum(weights));
	}
	else
	{
		cout << "Error in wStdDev: Weights and List do not match." << endl;
		throw logic_error("");
	}
}

double wStdDev(const vector<int> &l, const vector<double> &weights)
{
	if (weights.size() == 0)
		return stdDev(l);
	else if (weights.size() == l.size())
	{
		double m = wAvg(l, weights);
		double acc = 0.0;
		for (unsigned int i = 0; i < l.size(); i++)
		{
			acc += weights[i] * (l[i] - m) * (l[i] - m);
		}
		return sqrt(acc / sum(weights));
	}
	else
	{
		cout << "Error in wStdDev: Weights and List do not match." << endl;
		throw logic_error("");
	}
}

int numNodesInTree(const vector<Edge> &t, bool warnIfSomeMissing)
{
	int maxID = INT_MIN;
	for (unsigned i = 0; i < t.size(); i++)
	{
		if (t[i].u > maxID)
			maxID = t[i].u;
		if (t[i].v > maxID)
			maxID = t[i].v;
	}
	int initNumNodes = maxID + 1;
	vector<bool> inTree(maxID+1, false);
	for (unsigned i = 0; i < t.size(); i++)
	{
		inTree[t[i].u] = true;
		inTree[t[i].v] = true;
	}
	int numNodes = 0;
	for (unsigned i = 0; i < inTree.size(); i++)
		numNodes += (int)inTree[i];
	if (warnIfSomeMissing && initNumNodes != numNodes)
		cout << "Warning: numNodesInTree called on incomplete Tree" << endl;
	for (unsigned i = 0; i < inTree.size(); i++)
	{
		if (!inTree[i] && warnIfSomeMissing)
			cout << "Warning: Node #" << i << " is not in the tree!!" << endl;
	}
	return numNodes;
}

int isMem(int num, const vector<pair<int, int>> &degDist)
{
	for (unsigned i = 0; i < degDist.size(); i++)
	{
		if (num == degDist[i].first)
			return i;
	}
	return -1;
}

//**If called on a tree with duplicate edges it will consider them both as unique edges**
treeInfo getTreeInfo(const vector<Edge> &t, Graph *g)
{
	treeInfo i;
	i.numNodes = numNodesInTree(t, true);
	
	vector<int> adjs(i.numNodes, 0);
	for (unsigned i = 0; i < t.size(); i++)
	{
		adjs[t[i].u]++;
		adjs[t[i].v]++;
	}
	vector<pair<int, int>> dd;
	for (unsigned i = 0; i < adjs.size(); i++)
	{
		int loc = isMem(adjs[i], dd);
		if (loc == -1)
		{
			pair<int, int> p(adjs[i], 1);
			dd.push_back(p);
		}
		else
			dd[loc].second++;
	}
	i.degreeDist = dd;

	double mean = 0;
	double denom = 0;
	for (unsigned i = 0; i < dd.size(); i++)
	{
		mean += dd[i].first * dd[i].second;
		denom += dd[i].second;
	}
	mean /= denom;
	i.mean = mean;
	int mode = 0;
	int modeOccs = 0;
	for (unsigned j = 0; j < dd.size(); j++)
	{
		if (dd[j].second > modeOccs)
		{
			mode = dd[j].first;
			modeOccs = dd[j].second;
		}
	}
	i.mode = mode;
	vector<Edge> tMut = t;
	int numEIM = 0;
	i.matchingCost = matchingCost(g, tMut, numEIM);
	double treeCost = 0;
	for (unsigned i = 0; i < t.size(); i++)
		treeCost += g->distance(t[i].u, t[i].v);
	i.treeCost = treeCost;
	i.numEdges = t.size();
	i.numEdgesIncludeMatching = t.size() + numEIM;
	i.numEdgesInMatching = numEIM;
	i.avgEdgeCostInMatching = i.matchingCost / i.numEdgesInMatching;
	return i;
}


void output(treeInfo i)
{
	cout << "Number of Nodes: " << i.numNodes << endl;
	cout << "Number of Edges: " << i.numEdges << endl;
	cout << "Number of Edges including MPM: " << i.numEdgesIncludeMatching << endl;
	cout << "Degree Distribution:" << endl;
	for (unsigned j = 0; j < i.degreeDist.size(); j++)
		cout << "Degree: " << i.degreeDist[j].first << ", # Occurences: " << i.degreeDist[j].second << endl;
	cout << "Mean Degree: " << i.mean << endl;
	cout << "Mode Degree: " << i.mode << endl;
	cout << "Cost of tree: " << i.treeCost << endl;
	cout << "Cost of matching: " << i.matchingCost << endl;
	cout << "Total cost: " << i.treeCost + i.matchingCost << endl;
}

void checkOddsOfTreeSet(const vector<ModEdge> &lGraph, const vector<vector<Edge>> &ts)
{
	int numTrees = ts.size();
	vector<double> noOccs(lGraph.size(), 0.0);
	for (unsigned t = 0; t < ts.size(); t++)
	{
		for (unsigned e = 0; e < ts[t].size(); e++)
		{
			int u = ts[t][e].u;
			int v = ts[t][e].v;
			for (unsigned loc = 0; loc < lGraph.size(); loc++)
			{
				if ((lGraph[loc].end0 == u && lGraph[loc].end1 == v) || (lGraph[loc].end0 == v && lGraph[loc].end1 == u))
				{
					noOccs[loc] += 1.0;
				}
			}
		}
	}
	for (unsigned i = 0; i < noOccs.size(); i++)
		noOccs[i] /= numTrees;
	cout << "Comparison between lambda and odds: " << endl;
	for (unsigned i = 0; i < lGraph.size(); i++)
		cout << " From " << lGraph[i].end0 << " to " << lGraph[i].end1 << ": Lambda: " << lGraph[i].weight << ", Frequency: " << noOccs[i] << endl;
}

pair<double,double> checkNegativeCorrelationOld(const vector<ModEdge> &lGraph, const vector<vector<Edge>> &ts)
{
	cout << "Running Negative Correlation Check:\nPicking Subset of Edges:" << endl;
	vector<ModEdge> es;
	vector<int> added;
	while (added.size() < 5)
	{
		int next = rand() % lGraph.size();
		bool found = false;
		for (unsigned i = 0; i < added.size(); i++)
		{
			if (added[i] == next)
			{
				found = true;
				break;
			}
				
		}
		if (found)
			continue;
		added.push_back(next);
		es.push_back(lGraph[next]);
	}
	output(es);
	cout << "Finding original expectation: ";
	double p1 = 1.0;
	for (unsigned i = 0; i < es.size(); i++)
		p1 *= es[i].weight;
	cout << p1 << endl;
	cout << "Finding new expectation: ";
	double p2 = 1.0;
	for (unsigned i = 0; i < es.size(); i++)
	{
		int count = 0;
		for (unsigned j = 0; j < ts.size(); j++)
		{
			bool found = false;
			for (unsigned k = 0; k < ts[j].size(); k++)
			{
				if ((ts[j][k].u == es[i].end0 && ts[j][k].v == es[i].end1) ||
					(ts[j][k].u == es[i].end1 && ts[j][k].v == es[i].end0))
				{
					found = true;
					break;
				}
			}
			if (found)
				count++;
		}
		p2 *= ((double)count) / ((double)ts.size());
	}
	cout << p2 << endl;
	return pair<double, double>(p1, p2);
}

pair<double, double> checkNegativeCorrelation(const vector<ModEdge> &lGraph, const vector<vector<Edge>> &ts)
{
	cout << "Running Negative Correlation Check:\nPicking Subset of Edges:" << endl;
	vector<ModEdge> es;
	vector<int> added;
	unsigned int subsetSize = min(lGraph.size(), 40);
	cout << "Num Edges: " << lGraph.size() << endl;
	while (added.size() < subsetSize)
	{
		int next = rand() % lGraph.size();
		bool found = false;
		for (unsigned i = 0; i < added.size(); i++)
		{
			if (added[i] == next)
			{
				found = true;
				break;
			}

		}
		if (found)
			continue;
		added.push_back(next);
		es.push_back(lGraph[next]);
	}
	output(es);
	cout << "Finding original expectation: ";
	double p1 = 1.0;
	for (unsigned i = 0; i < es.size(); i++)
		p1 *= es[i].weight;
	cout << p1 << endl;
	cout << "Finding new expectation: ";
	double p2 = 0.0;

	double p3 = 1.0;
	for (unsigned i = 0; i < es.size(); i++)
		p3 *= 1.0 - es[i].weight;
	double p4 = 0.0;

	for (unsigned j = 0; j < ts.size(); j++)
	{
		bool foundAll = true;
		bool foundNone = true;
		bool *found = new bool[subsetSize];
		for (unsigned int i = 0; i < subsetSize; i++)
			found[i] = false;
		for (unsigned k = 0; k < ts[j].size(); k++)
		{	
			for (unsigned int i = 0; i < es.size(); i++)
			{
				if ((ts[j][k].u == es[i].end0 && ts[j][k].v == es[i].end1) ||
					(ts[j][k].u == es[i].end1 && ts[j][k].v == es[i].end0))
				{
					found[i] = true;
					foundNone = false;
					break;
				}
			}
		}
		for (unsigned int i = 0; i < subsetSize; i++)
		{
			if (!found[i])
				foundAll = false;
		}
		if (foundAll)
			p2 += 1.0 / ts.size();
		if (foundNone)
			p4 += 1.0 / ts.size();
		delete[] found;
	}
	cout << p2 << endl;
	cout << "Reverse Correlation Original: " << p3 << endl;
	cout << "Reverse Correlation New: " << p4 << endl;
	return pair<double, double>(p1, p2);
}


void checkConvComb(const vector<LPEdge> &hkSoln, int k, const pair<vector<vector<JoinEdge>>, vector<JoinEdge>> &ts)
{
	cout << "Verifying convex combination of trees...";
	vector<int> noOccs(hkSoln.size(), 0);
	for (unsigned t = 0; t < ts.first.size(); t++)
	{
		for (unsigned e = 0; e < ts.first[t].size(); e++)
		{
			int u = ts.first[t][e].u;
			int v = ts.first[t][e].v;
			for (unsigned loc = 0; loc < hkSoln.size(); loc++)
			{
				if ((hkSoln[loc].end0 == u && hkSoln[loc].end1 == v) || (hkSoln[loc].end0 == v && hkSoln[loc].end1 == u))
				{
					noOccs[loc]++;
				}
			}
		}
	}
	for (unsigned e = 0; e < ts.second.size(); e++)
	{
		int u = ts.second[e].u;
		int v = ts.second[e].v;
		for (unsigned loc = 0; loc < hkSoln.size(); loc++)
		{
			if ((hkSoln[loc].end0 == u && hkSoln[loc].end1 == v) || (hkSoln[loc].end0 == v && hkSoln[loc].end1 == u))
			{
				noOccs[loc]++;
			}
		}
	}
	for (unsigned i = 0; i < hkSoln.size(); i++)
	{
		int shouldBe = (int)((hkSoln[i].weight * k) + 0.5);
		if (noOccs[i] != shouldBe)
			cout << endl << "ERROR: should be " << shouldBe << "but is " << noOccs[i] << endl;
	}
	cout << "Done." << endl;
}


