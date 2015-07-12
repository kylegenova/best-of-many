#include "random_walk.h"


int getMaxInd(const vector<ModEdge> &t)
{
	int maxInd = 0;
	for (unsigned i = 0; i < t.size(); i++)
	{
		if (t[i].end0 > maxInd)
			maxInd = t[i].end0;
		if (t[i].end1 > maxInd)
			maxInd = t[i].end1;
	}
	return maxInd + 1;
}

AdjList getAdjs(const vector<ModEdge> &t)
{
	vector<vector<pair<int, double>>> adjs;
	int maxInd = getMaxInd(t);
	for (int i = 0; i < maxInd; i++)
		adjs.push_back(vector<pair<int, double>>());
	for (unsigned i = 0; i < t.size(); i++)
	{
		adjs[t[i].end0].push_back(pair<int, double>(t[i].end1, t[i].weight));
		adjs[t[i].end1].push_back(pair<int, double>(t[i].end0, t[i].weight));
	}
	return adjs;
}

int pickOneWeighted(vector<int> choices, vector<double> weights)
{
	double wSum = 0.0;
	for (unsigned i = 0; i < weights.size(); i++)
		wSum += weights[i];
	double frac = ((double)rand()) / ((double)RAND_MAX);
	double chosen = wSum * frac;
	for (unsigned i = 0; i < choices.size(); i++)
	{
		if (weights[i] > chosen)
			return i;
		chosen -= weights[i];
	}
	//cout << "Error: Nothing chosen!!" << endl;
	return choices.size() - 1;
}

int pickOneWeighted(const vector<pair<int, double>> &in, double uniformRand)
{
	double wSum = 0.0;
	for (unsigned i = 0; i < in.size(); i++)
		wSum += in[i].second;
	//double frac = ((double)rand()) / ((double)RAND_MAX);
	double chosen = wSum * uniformRand;
	for (unsigned i = 0; i < in.size(); i++)
	{
		if (in[i].second > chosen)
			return i;
		chosen -= in[i].second;
	}
	//cout << "Warning: Nothing chosen!!" << endl;
	return in.size() - 1;
}

vector<LPEdge> randomWalk(vector<ModEdge> lGraph, bool silent)
{
	//cout << "Entered Random Walk" << endl;
	vector<LPEdge> tree;
	int maxInd = getMaxInd(lGraph);
	tree.reserve(maxInd - 1);
	AdjList adjs = getAdjs(lGraph);
	vector<bool> added(maxInd, false);
	int numInTree = 0;
	int curLoc = 0;
	auto seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
	mt19937 mersenne = std::mt19937(seed);
	uniform_real_distribution<double> distr(0.0f, 1.0f);
	while (numInTree != maxInd-1)
	{
		std::shuffle(adjs[curLoc].begin(), adjs[curLoc].end(), mersenne);
		double randNum = distr(mersenne);
		int toCheck = pickOneWeighted(adjs[curLoc], randNum);
		int n = adjs[curLoc][toCheck].first;
		double w = adjs[curLoc][toCheck].second;
		if (!added[n])
		{
			tree.push_back(LPEdge(curLoc, n, w));
			if (!added[curLoc])
				added[curLoc] = true;
			added[n] = true;
			numInTree++;
			//cout << "Added an edge to the tree from " << numInTree - 1 << " to " << numInTree << " out of " << maxInd - 1 << endl;
		}
		curLoc = n;
	}
	return tree;
}