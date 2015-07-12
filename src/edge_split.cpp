// edgeSplit.cpp
// Implementation for the edge splitting procedure
// Takes in a graph, and returns a convex combination of trees where every edge except 1 has been split.
#include "edge_split.h"



historyIndex createLog()
{
	historyIndex h;
	h.first = new history(-1, -1, -1, -1);
	h.loc = h.first;
	h.last = h.first;
	return h;
}

void logSplit(historyIndex &h, int s, int u, int v, int weight)
{
	//cout << "logSplit called on s=" << s << ", u=" << u << ", v=" << v << ", weight=" << weight << endl;
	history *cur = h.last;
	cur->next = new history(s, u, v, weight);
	cur->next->prev = cur;
	h.last = cur->next;
}

void output(const history &h)
{
	cout << "S = " << h.s << ", U = " << h.u << ", V = " << h.v << ", weight = " << h.multiplicity << endl;
}

int historySize(historyIndex &h)
{
	history *ol = h.first;
	int size = 0;
	while (true)
	{
		if (ol->next == nullptr)
			return size;
		size++;
		ol = ol->next;
	}
}

void eraseSplit(historyIndex &h, int s, int u, int v, int weight)
{
	//cout << "Called eraseSplit on s=" << s << ", u=" << u << ", v=" << v << ", weight=" << weight << endl;
	int startSize = historySize(h);
	//cout << "Current history size is " << startSize << endl;
	h.loc = h.last;
	history toErase(s, u, v, weight);
	if (*h.last == toErase)
	{
		//cout << "Went down last erase path." << endl;
		if (h.last->multiplicity == weight)
		{
			h.last->prev->next = nullptr;
			h.last = h.last->prev;
		}
		else
			h.last->multiplicity -= weight;
		return;
	}
	//cout << "Went down full path." << endl;
	h.loc = h.last;
	//cout << "This is h" << endl;
	//cout << "last:" << endl;
	//output(*h.last);
	//cout << "first: " << endl;
	//output(*h.first);
	//cout << "current:" << endl;
	//output(*h.loc);
	while (h.loc->prev != nullptr)
	{

		h.loc = h.loc->prev;
		//cout << "This is h" << endl;
		//cout << "last:" << endl;
		//output(*h.last);
		//cout << "first: " << endl;
		//output(*h.first);
		//cout << "current:" << endl;
		//output(*h.loc);
		if (toErase == *h.loc)
		{
			if (h.loc->multiplicity == weight)
			{
				if (!(h.loc->isFirst()))
				{
					h.loc->prev->next = h.loc->next;
					h.loc->next->prev = h.loc->prev;
				}
				else
				{
					cout << "Error: Could not find split in history. Quitting." << endl;
					throw logic_error("");
				}
			}
			else
				h.loc->multiplicity -= weight;
			break;
		}
	}
	int endSize = historySize(h);
	//cout << "Ending history size is: " << endSize << endl;
}

void output(historyIndex &h)
{
	h.loc = h.first;
	while (h.loc->next != nullptr)
	{
		h.loc = h.loc->next;
		cout << "Split at " << h.loc->s << " from " << h.loc->u << " to " << h.loc->v << " with multiplicity " << h.loc->multiplicity << endl;
	}
}

vector<int> maOrdering(const vector<SplitEdge> &g, int s)
{
	int addOne = 0;
	if (cG(s, g) == 0)
		addOne = 1;
	int numNodes = getNumUsedNodes(g);
	int maxNodeInd = getMaxNodeInd(g);
	vector<bool> usedNodes = getUsedNodes(g);
	vector<int> ordering;
	vector<vector<pair<int, int>>> adjacency;
	adjacency.reserve(maxNodeInd);
	for (int i = 0; i < maxNodeInd; i++)
		adjacency.push_back(vector<pair<int, int>>());
	for (unsigned i = 0; i < g.size(); i++)
	{
		adjacency[g[i].end0].push_back(pair<int, int>(g[i].end1, g[i].weight));
		adjacency[g[i].end1].push_back(pair<int, int>(g[i].end0, g[i].weight));
	}
	//for each vertex in left, sum the weights of all the edges which end in ordering.
	//take the greatest of these, move from left to ordering; repeat while 
	//left isn't empty.
	vector<bool> inOrdering(maxNodeInd, false);
	ordering.push_back(s);
	inOrdering[s] = true;
	vector<int> dels(maxNodeInd, 0);
	while (ordering.size() != numNodes+addOne)
	{
		function<int(vector<pair<int,int>>)> f = [&inOrdering](vector<pair<int,int>> e)
		{
			int acc = 0;
			for (unsigned i = 0; i < e.size(); i++)
			{
				if (inOrdering[e[i].first])
					acc += e[i].second;
			}
			return acc;
		};
		transform(adjacency.begin(), adjacency.end(), dels.begin(), f);
		for (int i = 0; i < maxNodeInd; i++)
			dels[i] = inOrdering[i] || !usedNodes[i] ? -1 : dels[i];
		vector<int>::iterator maxLoc = max_element(dels.begin(), dels.end());
		int loc = std::distance(dels.begin(), maxLoc);
		ordering.push_back(loc);
		inOrdering[loc] = true;
	}
	//verifyMAOrdering(g, ordering);
	return ordering;
}

vector<int> maOrderingFast(const vector<SplitEdge> &g, int s)
{
	int addOne = cG(s, g) == 0 ? 1 : 0;
	int maxNodeInd = 0;
	for (unsigned i = 0; i < g.size(); i++)
		maxNodeInd = max(g[i].end0, max(g[i].end1, maxNodeInd));
	maxNodeInd++;
	vector<bool> usedNodes(maxNodeInd, false);
	for (unsigned i = 0; i < g.size(); i++)
	{
		usedNodes[g[i].end0] = true;
		usedNodes[g[i].end1] = true;
	}
	int numNodes = 0;
	for (unsigned i = 0; i < usedNodes.size(); i++)
		numNodes += usedNodes[i];
	vector<int> ordering;
	vector<bool> inOrdering;
	for (int i = 0; i < maxNodeInd; i++)
		inOrdering.push_back(false);
	ordering.push_back(s);
	inOrdering[s] = true;

	while (ordering.size() != numNodes + addOne)
	{
		vector<int> dels(maxNodeInd, 0);
		for (unsigned i = 0; i < g.size(); i++)
		{
			if (!inOrdering[g[i].end0] && inOrdering[g[i].end1])
				dels[g[i].end0] += g[i].weight;
			if (!inOrdering[g[i].end1] && inOrdering[g[i].end0])
				dels[g[i].end1] += g[i].weight;
		}
		for (int i = 0; i < maxNodeInd; i++)
		{
			if (!usedNodes[i])
				dels[i] = -1;
		}
		int loc = -1;
		int maxDel = -1;
		for (unsigned i = 0; i < dels.size(); i++)
		{
			if (dels[i] > maxDel)
			{
				maxDel = dels[i];
				loc = i;
			}
		}
		ordering.push_back(loc);
		inOrdering[loc] = true;
	}
	//verifyMAOrdering(g, ordering);
	return ordering;
}

class MACompare {
public:
	bool operator() (pair<int, int> &p1, pair<int, int> &p2) {
		return p1.first < p2.first;
	}
};

vector<int> maOrderingHeap(const vector<SplitEdge> &g, int s)
{
	vector<pair<int, int>> heap; //fst is key/incidence, snd is node ind (matched in MACompare)
	int addOne = cG(s, g) == 0 ? 1 : 0;
	int maxNodeInd = 0;
	for (unsigned i = 0; i < g.size(); i++)
		maxNodeInd = max(g[i].end0, max(g[i].end1, maxNodeInd));
	maxNodeInd++;
	vector<bool> usedNodes(maxNodeInd, false);
	for (unsigned i = 0; i < g.size(); i++)
	{
		usedNodes[g[i].end0] = true;
		usedNodes[g[i].end1] = true;
	}
	int numNodes = 0;
	for (unsigned i = 0; i < usedNodes.size(); i++)
		numNodes += usedNodes[i];
	vector<int> ordering;
	vector<bool> inOrdering(maxNodeInd, false);
	ordering.push_back(s);
	inOrdering[s] = true;

	vector<vector<pair<int, int>>> adjacency;
	adjacency.reserve(maxNodeInd);
	for (int i = 0; i < maxNodeInd; i++)
		adjacency.push_back(vector<pair<int, int>>());
	for (unsigned i = 0; i < g.size(); i++)
	{
		adjacency[g[i].end0].push_back(pair<int, int>(g[i].end1, g[i].weight));
		adjacency[g[i].end1].push_back(pair<int, int>(g[i].end0, g[i].weight));
	}

	heap.reserve(maxNodeInd);
	for (int i = 0; i < maxNodeInd; i++)
		heap.push_back(pair<int, int>(0, i));

	for (unsigned i = 0; i < g.size(); i++)
	{
		if (g[i].end0 == s)
			heap[g[i].end1].first += g[i].weight;
		else if (g[i].end1 == s)
			heap[g[i].end0].first += g[i].weight;
	}
	vector<int>realValue;
	for (unsigned i = 0; i < heap.size(); i++)
		realValue.push_back(heap[i].first);
	make_heap(heap.begin(), heap.end(), MACompare());
	while (ordering.size() != numNodes + addOne)
	{
		pair<int, int> next = heap.front();
		pop_heap(heap.begin(), heap.end());
		heap.pop_back();

		while (realValue[next.second] != next.first || inOrdering[next.second] || !usedNodes[next.second]) //i.e. this value was changed.
		{
			next = heap.front();
			pop_heap(heap.begin(), heap.end());
			heap.pop_back();
		}

		ordering.push_back(next.second);
		inOrdering[next.second] = true;

		for (unsigned i = 0; i < adjacency[next.second].size(); i++)
		{
			int n = adjacency[next.second][i].first;
			if (!inOrdering[n])
			{
				realValue[n] += adjacency[next.second][i].second;
				pair<int, int> p(realValue[n], n);
				heap.push_back(p);
				push_heap(heap.begin(), heap.end());
			}
		}
	}
	//verifyMAOrdering(g, ordering);
	return ordering;
}

vector<SplitEdge> hookUpHelper(int s, const vector<SplitEdge> &g, const vector<SplitEdge> &toHookUp, historyIndex &h)
{
	vector<SplitEdge> out = g;
	for (unsigned i = 0; i < toHookUp.size(); i++)
	{
		int u = toHookUp[i].end0;
		int v = toHookUp[i].end1;
		int subtractInd = indexOfEdge(g, u, v);
		if (subtractInd == -1)
			throw logic_error("Trying to subtract weight from nonexistant edge in hookUpHelper");
		out[subtractInd].weight -= toHookUp[i].weight;
		int usInd = indexOfEdge(g, u, s);
		int vsInd = indexOfEdge(g, v, s);
		eraseSplit(h, s, u, v, toHookUp[i].weight);
		if (usInd == -1)
		{
			SplitEdge e(u, s, toHookUp[i].weight, u, s);
			//cout << "Creating edge between " << u << " and " << s << " with weight " << toHookUp[i].weight << " in hookUpHelper" << endl;
			out.push_back(e);
		}
		else
		{
			out[usInd].weight += toHookUp[i].weight;
			//cout << "Increasing the weight from " << u << " to " << s << " from " << out[usInd].weight << " to " << out[usInd].weight + toHookUp[i].weight << endl;
		}
		if (vsInd == -1)
		{
			SplitEdge e(v, s, toHookUp[i].weight, v, s);
			out.push_back(e);
			//cout << "Creating edge between " << v << " and " << s << " with weight " << toHookUp[i].weight << " in hookUpHelper" << endl;
		}
		else
		{ 
			out[vsInd].weight += toHookUp[i].weight;
			//cout << "Increasing the weight from " << v << " to " << s << " from " << out[vsInd].weight << " to " << out[vsInd].weight + toHookUp[i].weight << endl;
		}
	}
	return out;
}

struct huReturn
{
	vector<SplitEdge> G1;
	vector<SplitEdge> BP;
	vector<vector<int>> Y;
};

huReturn hookUp(const vector<SplitEdge> &g, int s, vector<SplitEdge> &bIn, int k, vector<vector<int>> &Y, historyIndex &h)
{
	if (cG(s, g) != 0)
	{
		cout << "cG(s) problem in hookup" << endl;
		throw logic_error("");
	}
	vector<SplitEdge> H = g;
	vector<SplitEdge> G1 = g; 
	vector<SplitEdge> B = bIn;
	vector<SplitEdge> B1;
	vector<vector<int>> XS;
	int maxNodeInd = getMaxNodeInd(G1);
	for (int i = 0; i < maxNodeInd; i++)
		XS.push_back(vector<int>());
	for (int i = 0; i < maxNodeInd; i++)
		XS[i].push_back(i);
	//cout << "About to enter while loop" << endl;
	while (getNumUsedNodes(H) >= 4)
	{
		vector<int> ma = maOrderingHeap(H, s);
		int v = ma[ma.size() - 2];
		int w = ma[ma.size() - 1];
		if (v == s || w == s)
			throw logic_error("SET WAS V - S, S FOUND");
		vector<int> X1;
		H = combineVertices(H, v, w);
		H = compress(H);
		XS[v] = setUnion(XS[v], XS[w]);
		if (XS[w].size() == 0)
		{
			cout << "Error: W, " << w << " was merged twice. Quitting" << endl;
			throw logic_error("");
		}
		XS[w] = vector<int>();
		if (cG(v, H) < k)
		{
			int numToGet = (int)ceil(.5*(double(k) - double(cG(G1, XS[v]))));
			vector<SplitEdge> GX = inducedSubgraph(G1, XS[v]);
			vector<SplitEdge> delB;
			int added = 0;
			for (unsigned i = 0; i < GX.size(); i++)
			{
				SplitEdge e = SplitEdge(GX[i].end0, GX[i].end1, GX[i].weight, GX[i].orig0, GX[i].orig1);
				if (isMem(e, B))
				{
					int bW = B[indexOfEdge(B, e.end0, e.end1)].weight;
					if (bW < e.weight)
						e.weight = bW;
					if (e.weight > (numToGet - added))
					{
						e.weight = numToGet - added;
					}
					added += e.weight;
					delB.push_back(e);
				}
				if (added == numToGet)
					break;
			}
			if (added != numToGet)
			{
				cout << "Error: GX did not contain " << numToGet << " entries in B. Quitting." << endl;
				throw logic_error("");
			}
			if (!isSubset(delB, B))
			{
				cout << "ERROR: delB is not a subset of B." << endl;
				cout << "B:" << endl;
				output(B);
				cout << "delB:" << endl;
				output(delB);
				cout << "This was the GX to choose from:" << endl;
				output(GX);
				cout << "V: " << v << endl;
				cout << "W: " << w << endl;
				cout << "S: " << s << endl;
				throw logic_error("");
			}
			B = setRemove(delB, B);
			B = removeZeroWeighted(B);
			B1 = setUnion(delB, B1);
			H = removeZeroWeighted(H);
			G1 = hookUpHelper(s, G1, delB, h);
			G1 = removeZeroWeighted(G1);
			H = removeZeroWeighted(H);
			bool addedFromXSinH = false;
			numToGet *= 2;
			for (unsigned i = 0; i < H.size(); i++)
			{
				SplitEdge tester = SplitEdge(s, v, 0, 0, 0);
				if (equals(tester, H[i]))
				{
					//cout << "Increasing weight in hookUp in H between " << H[i].end0 << " and " << H[i].end1 << "from " << H[i].weight << " to " << H[i].weight + numToGet << endl;
					H[i].weight += numToGet;
					addedFromXSinH = true;
					break;
				}
			}
			if (!addedFromXSinH && numToGet != 0)
			{
				//cout << "Creating edge in hookUp in H between " << s << " and " << v << " with weight " << numToGet << endl;
				SplitEdge e(s, v, numToGet, s, v);
				H.push_back(e);
			}
			vector<vector<int>> newY;
			for (unsigned i = 0; i < Y.size(); i++)
			{
				if (!isProperSubset(Y[i], XS[v]))
					newY.push_back(Y[i]);
			}
			bool foundX1inY = false;
			for (unsigned i = 0; i < newY.size(); i++)
			{
				if (setsEqual(newY[i], XS[v]))
					foundX1inY = true;
			}
			if (!foundX1inY)
				newY.push_back(XS[v]);
			Y = newY;
		}
	}
	huReturn ret;
	ret.BP = B1;
	ret.G1 = G1;
	ret.Y = Y;
	return ret;
}

vector<SplitEdge> split(const vector<SplitEdge> &g, int u, int v, int s, int delta, historyIndex &h)
{
	vector<SplitEdge> GHat = g;
	bool addedWeight = false;
	bool subtractedU = false;
	bool subtractedV = false;
	//cout << endl;
	for (unsigned i = 0; i < GHat.size(); i++)
	{
		if (connects(GHat[i], u, v))
		{
			//cout << "Increasing the weight from " << u << " to " << v << " from " << GHat[i].weight << " to " << GHat[i].weight + delta << endl;
			GHat[i].weight += delta;
			addedWeight = true;
		}
		if (connects(GHat[i], u, s) && !subtractedU)
		{
			//cout << "Decreasing the weight from " << u << " to " << s << " from " << GHat[i].weight << " to " << GHat[i].weight - delta << endl;
			GHat[i].weight -= delta;
			subtractedU = true;
		}
		if (connects(GHat[i], s, v) && !subtractedV)
		{
			//cout << "Decreasing the weight from " << v << " to " << s << " from " << GHat[i].weight << " to " << GHat[i].weight - delta << endl;
			GHat[i].weight -= delta;
			subtractedV = true;
		}
	}
	if (!addedWeight)
	{
		SplitEdge e = SplitEdge(u, v, delta, u, v);
		//cout << "Creating edge from " << u << " to " << v << " with weight " << delta << endl;
		GHat.push_back(e);
	}
	if (!subtractedU)
	{
		cout << "Error: Couldn't find weight from " << u << " to " << s << " to subtract." << endl;
		throw logic_error("");
	}
	if (!subtractedV)
	{
		cout << "Error: Couldn't find weight from " << v << " to " << s << " to subtract." << endl;
		throw logic_error("");
	}
	//cout << endl;
	GHat = removeZeroWeighted(GHat);
	logSplit(h, s, u, v, delta);
	return GHat;
}

vector<SplitEdge> pairing(const vector<int> &Xi,const vector<int> &Xj, int delij, vector<SplitEdge> &B, const vector<SplitEdge> &g, int s, historyIndex &h)
{
	vector<SplitEdge> GHat = g;
	while (delij > 0)
	{
		int u = 0;
		int v = 0;
		vector<int> gamXi = setIntersection(Xi, neighbors(GHat, s));
		vector<int> gamXj = setIntersection(Xj, neighbors(GHat, s));
		if (gamXi.size() == 0 || gamXj.size() == 0)
		{
			cout << "Detected 0, dumping status" << endl;
			cout << "GHat:" << endl;
			output(GHat);
			cout << "B:" << endl;
			output(B);
			cout << "S: " << s << endl;
			cout << "delij: " << delij << endl;
			cout << "Xi:" << endl;
			output(Xi);
			cout << "neighbors of s in GHat:" << endl;
			output(neighbors(GHat, s));
			cout << "Xj:" << endl;
			output(Xj);
			cout << "gamXi:" << endl;
			output(gamXi);
			cout << "gamXj:" << endl;
			output(gamXj);
		}
		u = gamXi[0];
		v = gamXj[0];
		if (u == v)
		{
			cout << "U == V. Quitting" << endl;
			throw logic_error("");
		}
		int delta = min(min(cG(s, u, GHat), cG(s, v, GHat)),delij);
		GHat = split(GHat, u, v, s, delta, h);
		//Add weight delta to (u,v) in B
		bool found = false;
		for (unsigned i = 0; i < B.size(); i++)
		{
			if (connects(B[i], u, v))
			{
				B[i].weight += delta;
				found = true;
				break;
			}
		}
		if (!found)
		{
			for (unsigned i = 0; i < GHat.size(); i++)
			{
				if (connects(GHat[i], u, v))
				{
					SplitEdge e(GHat[i].end0, GHat[i].end1, delta, GHat[i].orig0, GHat[i].orig1);
					B.push_back(e);
					found = true;
					break;
				}
			}
		}
		if (!found)
		{
			SplitEdge e(u, v, delta, u, v);
			B.push_back(e);
		}
		delij -= delta;
	}
	return GHat;
}

vector<SplitEdge> eulerCSplit(const vector<SplitEdge> &g, vector<SplitEdge> &B, int s, unsigned int k, vector<vector<int>> X, historyIndex &h)
{
	vector<SplitEdge> GHat = g;
	GHat = fixEquation4(GHat, s);
	if (cG(s, g) % 2 != 0)
	{
		cout << "ERROR: odd CG(s)" << endl;
		throw logic_error("");
	}
	vector<vector<int>> X1 = X;

	while (X1.size() >= 3)
	{
		X1 = sortByCG(GHat, s, X1);
		vector<int> cgsX1;
		for (unsigned i = 0; i < X1.size(); i++)
		{
			cgsX1.push_back(cG(s, X1[i], GHat));
		}
		int del1P = 0;
		int del2P = 0;
		if ((cG(s, X1[0], GHat) - cG(s, X1[1], GHat)) >= cG(s, X1[X1.size() - 1], GHat))
			del1P = cG(s, X1[X1.size() - 1], GHat);
		else
		{
			del2P = (int)ceil(double(cG(s, X1[X1.size() - 1], GHat) - cG(s, X1[0], GHat) + cG(s, X1[1], GHat)) * .5);
			del1P = cG(s, X1[X1.size() - 1], GHat) - del2P;
		}
		GHat = removeZeroWeighted(GHat);
		GHat = pairing(X1[0], X1[X1.size() - 1], del1P, B, GHat, s, h);
		GHat = removeZeroWeighted(GHat);
		GHat = pairing(X1[1], X1[X1.size() - 1], del2P, B, GHat, s, h);
		bool erase1 = cG(s, X1[0], GHat) == 0;
		bool erase2 = cG(s, X1[1], GHat) == 0;
		bool eraseLast = cG(s, X1.back(), GHat) == 0;
		vector<vector<int>> X11;
		if (!erase1)
			X11.push_back(X1[0]);
		if (!erase2)
			X11.push_back(X1[1]);
		for (unsigned i = 2; i < X1.size() - 1; i++)
			X11.push_back(X1[i]);
		if (!eraseLast)
			X11.push_back(X1.back());
		X1 = X11;
		if (X1.size() == 1)
		{
			cout << "ERROR: X1.size() == 1" << endl;
			throw logic_error("");
		}
	}
	int del12 = 0;
	if (X1.size() != 0)
	{
		if (cG(s, X1[0], GHat) != cG(s, X1[1], GHat))
		{
			cout << "eulerCSplit sanity check failure. Make sure cG(X1) != cG(X2)" << endl;
			throw logic_error("");
		}
		del12 = cG(s, X1[0], GHat);
		GHat = removeZeroWeighted(GHat);
		GHat = pairing(X1[0], X1[1], del12, B, GHat, s, h);
	}
	return GHat;
}

vector<SplitEdge> evenSplit(const vector<SplitEdge> &g, int s, int k, historyIndex &h) 
{
	vector<SplitEdge> G = g;
	vector<vector<int>> X;
	vector<int> ns = neighbors(G, s);
	for (unsigned i = 0; i < ns.size(); i++)
	{
		vector<int> cur;
		cur.push_back(ns[i]);
		X.push_back(cur);
	}
	vector<SplitEdge>B;
	while (X.size() != 0)
	{

		vector<SplitEdge> G1 = eulerCSplit(G, B, s, k, X, h);
		G1 = removeZeroWeighted(G1);
		/*if (cG(s, G) / 2 != magnitude(B))
		{
			cout << "Error: |B| and cG(s)/2 should be the same after C Split" << endl;
			cout << "cG(s, G) /2: " << cG(s, G)/2 << endl;
			cout << "|B|: " << magnitude(B) << endl;
			throw logic_error("");
		}*/
		vector<vector<int>> Y;
		huReturn hookUpOut = hookUp(G1, s, B, k, Y, h);
		vector<SplitEdge> G2 = hookUpOut.G1;
		Y = hookUpOut.Y;
		B = hookUpOut.BP;
		G = G2;
		X = Y;
	}
	return G;
}

//currently finds GCD to 15 digits (number of guaranteed significands in a double
//reweights g as integers * k
//k is the actual 2k value
pair<unsigned int,vector<SplitEdge>> getK(const vector<ModEdge> &g)
{
	/*cout << "Generating k for these weights:" << endl;
	for (unsigned i = 0; i < g.size(); i++)
	{
		cout << g[i].weight << endl;
	}*/
	vector<SplitEdge> out;
	vector<long long int> lambda;
	for (unsigned i = 0; i < g.size(); i++)
		lambda.push_back((long long int)(ceil(g[i].weight * 1e3)));
	vector<long long int> in(lambda);
	long long int d = gcd(in);
	for (unsigned i = 0; i < g.size(); i++)
	{
		int w = (int)(lambda[i] / d);
		SplitEdge e(g[i].end0, g[i].end1, w, g[i].orig0, g[i].orig1);
		out.push_back(e);
	}
	long long int outN = (long long int)(1*2e3 / d);
	//cout << "Outputting the new weights: " << endl;
	for (unsigned i = 0; i < out.size(); i++)
	{
		//cout << out[i].weight << endl;
	}
	return pair<int, vector<SplitEdge>>((int)outN, out);
}

bool isInteger(double in)
{
	return (abs(in - round(in))) < .01;
}

pair<unsigned, vector<SplitEdge>> getKLinear(const vector<ModEdge> &g)
{
	/*cout << "Generating k for these weights:" << endl;
	for (unsigned i = 0; i < g.size(); i++)
	{
		cout << g[i].weight << endl;
	}*/
	vector<SplitEdge> out;
	int curK = 1;
	cout << "Finding K...";
	while (true)
	{
		bool done = true;
		for (unsigned i = 0; i < g.size(); i++)
		{
			if (!isInteger(g[i].weight*((double)curK)))
			{
				done = false;
				break;
			}
		}
		if (done)
			break;
		else
			curK += 1;
	}

	for (unsigned i = 0; i < g.size(); i++)
	{
		SplitEdge e(g[i].end0, g[i].end1, (int)round((g[i].weight*(double)curK)), g[i].orig0, g[i].orig1);
		out.push_back(e);
	}
	cout << "Done: " << curK << "." << endl;
	return pair<unsigned, vector<SplitEdge>>(curK * 2, out);
}

//Will need to change these headers to reflect the fact that the exact splitting history must 
//be stored
vector<SplitEdge> edgeSplit(const vector<SplitEdge> &g, int k, historyIndex &h)
{
	cout << "Beginning Edge Splitting...";
	vector<SplitEdge> cur = g;
	int maxInd = getMaxNodeInd(g);
	for (int toSplit = 0; toSplit < maxInd - 2; toSplit++)
	{
		vector<bool> used = getUsedNodes(cur);
		if (used[toSplit])
			cur = evenSplit(cur, toSplit, k, h);
	}
	cout << "Done." << endl;
	return cur;
}
