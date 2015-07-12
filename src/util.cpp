#include "util.h"



double sum(const vector<double> &l)
{
	double acc = 0.0;
	for (unsigned i = 0; i < l.size(); i++)
		acc += l[i];
	return acc;
}

//Assumes Nonnegative values
double lMax(const vector<double> &l)
{
	double max = 0.0;
	for (unsigned i = 0; i < l.size(); i++)
	{
		if (l[i] > max)
			max = l[i];
	}
	return max;
}

double sum(const vector<int> &l)
{
	double acc = 0.0;
	for (unsigned i = 0; i < l.size(); i++)
		acc += (double)l[i];
	return acc;
}

vector<double> elementWiseSum(const vector<double> &a, const vector<double> &b)
{
	if (a.size() != b.size())
		throw logic_error("elementWiseSum only defined for equal length lists");
	vector<double> out;
	for (unsigned i = 0; i < a.size(); i++)
		out.push_back(a[i] + b[i]);
	return out;
}

double wSum(const vector<int> &l, const vector<double> &weights)
{
	double acc = 0.0;
	if (weights.size() == 0)
		return sum(l);
	else if (l.size() == weights.size())
	{
		for (unsigned int i = 0; i < l.size(); i++)
			acc += l[i] * weights[i];
		return acc;
	}
	else
	{
		cout << "Error in wSum: Weights and List do not match" << endl;
		throw logic_error("");
	}
}

double wSum(const vector<double> &l, const vector<double> &weights)
{
	double acc = 0.0;
	if (weights.size() == 0)
		return sum(l);
	else if (l.size() == weights.size())
	{
		for (unsigned int i = 0; i < l.size(); i++)
			acc += l[i] * weights[i];
		return acc;
	}
	else
	{
		cout << "Error in wSum: Weights and List do not match" << endl;
		throw logic_error("");
	}
}

void output(const vector<int> &l)
{
	for_each(l.begin(), l.end(), [](int x) { cout << x << endl; });
}

int cG(int v, const vector<SplitEdge> &g)
{
	int acc = 0;
	for (unsigned i = 0; i < g.size(); i++)
	{
		if (g[i].end0 == v || g[i].end1 == v)
			acc += g[i].weight;
	}
	return acc;
}

int cG(const vector<SplitEdge> &g, const vector<int> &X)
{
	vector<bool> inX(getMaxNodeInd(g), false);
	for (unsigned i = 0; i < X.size(); i++)
		inX[X[i]] = true;
	int acc = 0;
	for (unsigned i = 0; i < g.size(); i++)
	{
		int u = g[i].end0;
		int v = g[i].end1;
		if ((inX[u] && !inX[v]) || (!inX[u] && inX[v]))
			acc += g[i].weight;
	}
	return acc;
}

int delta(const vector<SplitEdge> &g, int v)
{
	int acc = 0;
	for (unsigned i = 0; i < g.size(); i++)
	{
		if (g[i].end0 == v || g[i].end1 == v)
			acc += g[i].weight;
	}
	return acc;
}

vector<SplitEdge> inducedSubgraph(const vector<SplitEdge> &g, const vector<int> &V)
{
	vector<SplitEdge> out;
	vector<bool> inV(getMaxNodeInd(g), false);
	for (unsigned i = 0; i < V.size(); i++)
		inV[V[i]] = true;
	for (unsigned i = 0; i < g.size(); i++)
	{
		if (inV[g[i].end0] && inV[g[i].end1])
			out.push_back(g[i]);
	}
	return out;
}

bool equals(SplitEdge e1, SplitEdge e2)
{
	return (e1.end0 == e2.end0 && e1.end1 == e2.end1) ||
		(e1.end1 == e2.end0 && e1.end0 == e2.end1);
}

bool isMem(SplitEdge what, const vector<SplitEdge> &of)
{
	for (unsigned i = 0; i < of.size(); i++)
	{
		if (equals(what, of[i]))
			return true;
	}
	return false;
}

int indexOfEdge(const vector<SplitEdge> &in, int from, int to)
{
	int index = -1;
	for (unsigned i = 0; i < in.size(); i++)
	{
		if (connects(in[i], from, to))
			index = i;
	}
	return index;
}

//O(n^2)
vector<SplitEdge> setRemove(const vector<SplitEdge> &what, const vector<SplitEdge> &from)
{
	vector<SplitEdge> out;
	for (unsigned i = 0; i < from.size(); i++)
	{
		int index = indexOfEdge(what, from[i].end0, from[i].end1);
		if (index == -1)
			out.push_back(from[i]);
		else if (what[index].weight != from[i].weight)
		{
			//cout << "setRemove doing an ambiguous action" << endl;
			SplitEdge e(from[i].end0, from[i].end1, from[i].weight - what[index].weight, from[i].orig0, from[i].orig1);
			out.push_back(e);
		}
	}
	return out;
}

//O(n^2)
vector<SplitEdge> setUnion(const vector<SplitEdge> &s1, const vector<SplitEdge> &s2)
{
	vector<SplitEdge> out = s1;
	for (unsigned i = 0; i < s2.size(); i++)
	{
		int index = indexOfEdge(out, s2[i].end0, s2[i].end1);
		if (index == -1)
			out.push_back(s2[i]);
		else
		{
			out[index].weight += s2[i].weight;
		}
	}
	return out;
}

vector<SplitEdge> mergeWeights(const vector<SplitEdge> &g, const vector<SplitEdge> &h)
{
	vector<SplitEdge> out;
	vector<bool> addedFromH(h.size(), false);
	for (unsigned i = 0; i < g.size(); i++)
	{
		for (unsigned j = 0; j < h.size(); j++)
		{
			if (equals(g[i], h[j]))
			{
				SplitEdge e = SplitEdge(g[i].end0, g[i].end1, g[i].weight + h[j].weight, g[i].orig0, g[i].orig1);
				out.push_back(e);
				addedFromH[j] = true;
				break;
			}
			if (j == h.size() - 1)
			{
				SplitEdge e = SplitEdge(g[i].end0, g[i].end1, g[i].weight, g[i].orig0, g[i].orig1);
				out.push_back(e);
			}
		}
	}
	for (unsigned i = 0; i < h.size(); i++)
	{
		if (!addedFromH[i])
		{
			SplitEdge e = SplitEdge(h[i].end0, h[i].end1, h[i].weight, h[i].orig0, h[i].orig1);
			out.push_back(e);
		}
	}
	return out;
}

bool isMem(int what, const vector<int> &of)
{
	for (unsigned i = 0; i < of.size(); i++)
	{
		if (what == of[i])
			return true;
	}
	return false;
}

//O(n^2)
vector<int> setUnion(const vector<int> &s1, const vector<int> &s2)
{
	vector<int> out = s1;
	for (unsigned i = 0; i < s2.size(); i++)
	{
		if (!isMem(s2[i], out))
			out.push_back(s2[i]);
	}
	return out;
}

//O(n^2)
vector<int> setIntersection(const vector<int> &s1, const vector<int> &s2)
{
	vector<int> out;
	for (unsigned i = 0; i < s2.size(); i++)
	{
		if (isMem(s2[i], s1))
			out.push_back(s2[i]);
	}
	return out;
}

bool isSubset(const vector<int> &subset, const vector<int> &set)
{
	for (unsigned i = 0; i < subset.size(); i++)
	{
		if (!isMem(subset[i], set))
			return false;
	}
	return subset.size() <= set.size();
}

bool isProperSubset(const vector<int> &subset, const vector<int> &set)
{
	return isSubset(subset, set) && subset.size() < set.size();
}

bool setsEqual(const vector<int> &s1, const vector<int> &s2)
{
	if (s1.size() != s2.size())
		return false;
	for (unsigned i = 0; i < s1.size(); i++)
	{
		if (!isMem(s1[i], s2))
			return false;
	}
	return true;
}

bool vectorsIdentical(const vector<int> &l, const vector<int> &l2)
{
	if (l.size() != l2.size())
		return false;
	for (unsigned i = 0; i < l.size(); i++)
	{
		if (l[i] != l2[i])
			return false;
	}
	return true;
}

vector<SplitEdge> deleteEdge(const vector<SplitEdge> &g, SplitEdge e)
{
	vector<SplitEdge> out;
	for (unsigned i = 0; i < g.size(); i++)
	{
		if (!equals(g[i], e))
			out.push_back(g[i]);
	}
	return out;
}

vector<SplitEdge> combineVertices(const vector<SplitEdge> &g, int v, int w)
{
	vector<SplitEdge> H;
	for (unsigned i = 0; i < g.size(); i++)
	{
		SplitEdge cur = SplitEdge(g[i].end0, g[i].end1, g[i].weight, g[i].orig0, g[i].orig1);
		if (!((cur.end0 == v && cur.end1 == w) || (cur.end0 == w && cur.end1 == v)))
		{
			if (cur.end0 == w)
				cur.end0 = v;
			if (cur.end1 == w)
				cur.end1 = v;
			H.push_back(cur);
		}
	}
	return H;
}

//Warning: O(n) in time and space where n is the VALUE of the largest integer in what or from.
vector<int> setRemove(vector<int> what, vector<int> from)
{
	vector<int> out;
	int biggest = INT_MIN;
	for (unsigned i = 0; i < what.size(); i++)
	{
		if (what[i] > biggest)
			biggest = what[i];
	}
	vector<bool> inWhat(biggest + 1, false);
	for (unsigned i = 0; i < what.size(); i++)
	{
		inWhat[what[i]] = true;
	}
	for (unsigned i = 0; i < from.size(); i++)
	{
		if (!inWhat[from[i]])
			out.push_back(from[i]);
	}
	return out;
}


int cG(int from, int to, const vector<SplitEdge> &g)
{
	for (unsigned i = 0; i < g.size(); i++)
	{
		if (connects(g[i], from, to))
			return g[i].weight;
	}
	return 0;
}

int cG(int s, const vector<int> &X, const vector<SplitEdge> &g)
{
	int acc = 0;
	for (unsigned i = 0; i < g.size(); i++)
	{
		if (g[i].end0 == s && isMem(g[i].end1, X))
			acc += g[i].weight;
		else if (g[i].end1 == s && isMem(g[i].end0, X))
			acc += g[i].weight;
	}
	return acc;
}

vector<int> neighbors(const vector<SplitEdge> &g, int v)
{
	vector<int> n;
	for (unsigned i = 0; i < g.size(); i++)
	{
		if (g[i].end0 == v)
			n.push_back(g[i].end1);
		//changed from else if
		if (g[i].end1 == v)
			n.push_back(g[i].end0);
	}
	return n;
}

bool connects(SplitEdge e, int first, int second)
{
	SplitEdge tester(first, second, 0, 0, 0);
	return equals(e, tester);
}

void output(const vector<long long int> &l)
{
	for (unsigned i = 0; i < l.size(); i++)
		cout << l[i] << endl;
}

bool even(long long int a)
{
	return a % 2 == 0;
}

long long int gcd(long long int a, long long int b)
{
	if (a == 0)
		return b;
	if (b == 0)
		return a;
	if (even(a) && even(b))
		return 2 * gcd(a / 2, b / 2);
	if (even(a))
		return gcd(a / 2, b);
	if (even(b))
		return gcd(a, b / 2);
	if (a > b)
		return gcd((a - b) / 2, b);
	return gcd(a, (b - a) / 2);
}

long long int gcd(vector<long long int> &l)
{
	//cout << "Calculating the gcd of this list: " << endl;
	//output(l);
	if (l.size() == 2)
		return gcd(l[0], l[1]);
	long long int last = l[l.size() - 1];
	l.erase(l.end() - 1);
	return gcd(last, gcd(l));
}

void outputDif(const vector<SplitEdge> &g, const vector<SplitEdge> &h, int s)
{
	if (g.size() != h.size())
	{
		cout << "Error: g has size " << g.size() << " but h has size " << h.size() << endl;
		return;
	}
	for (unsigned i = 0; i < g.size(); i++)
	{
		if (g[i].weight != h[i].weight)
		{
			if (!equals(g[i], h[i]))
			{
				cout << "Error: graphs are not properly sorted." << endl;
				return;
			}
			cout << "The edge from node " << g[i].end0 << " to " << g[i].end1 << " used to have weight ";
			cout << g[i].weight << " but now has weight " << h[i].weight << endl;
		}
	}
	if (s != -1)
	{
		int accS = 0;
		int accOther = 0;
		for (unsigned i = 0; i < g.size(); i++)
		{
			if (h[i].end0 == s || g[i].end1 == s)
				accS += h[i].weight - g[i].weight;
			else
				accOther += h[i].weight - g[i].weight;
		}
		if (accS >= 0)
			cout << "Node #" << s << " gained total weight " << accS << endl;
		else
			cout << "Node #" << s << " lost total weight " << accS * -1 << endl;
		cout << "All other nodes ";
		if (accOther >= 0)
			cout << "gained total weight " << accOther << endl;
		else
			cout << "lost total weight " << accOther * -1 << endl;
	}
}

void output(vector<vector<int>> ll)
{
	for (unsigned i = 0; i < ll.size(); i++)
	{
		vector<int> l = ll[i];
		for (unsigned j = 0; j < l.size(); j++)
		{
			cout << l[j] << ", ";
		}
		cout << endl;
	}
}

vector<SplitEdge> removeZeroWeighted(const vector<SplitEdge> &in)
{
	vector<SplitEdge> out;
	for (unsigned i = 0; i < in.size(); i++)
	{
		if (in[i].weight != 0)
			out.push_back(in[i]);
	}
	return out;
}

vector<SplitEdge> slice(const vector<SplitEdge> &in, int from, int to)
{
	vector<SplitEdge> out;
	for (int i = from; i < to; i++)
	{
		out.push_back(in[i]);
	}
	return out;
}

int magnitude(const vector<SplitEdge> &g)
{
	int acc = 0;
	for (unsigned i = 0; i < g.size(); i++)
		acc += g[i].weight;
	return acc;
}

bool isSubset(const vector<SplitEdge> &what, const vector<SplitEdge> &of)
{
	for (unsigned i = 0; i < what.size(); i++)
	{
		int u = what[i].end0;
		int v = what[i].end1;
		bool found = false;
		for (unsigned j = 0; j < of.size(); j++)
		{
			if (connects(of[j], u, v))
			{
				found = true;
				break;
			}
		}
		if (!found)
			return false;
	}
	return what.size() <= of.size();
}

int getNumUsedNodes(const vector<SplitEdge> &in)
{
	int maxInd = 0;
	for (unsigned i = 0; i < in.size(); i++)
	{
		if (in[i].end0 > maxInd)
			maxInd = in[i].end0;
		if (in[i].end1 > maxInd)
			maxInd = in[i].end1;
	}
	vector<bool> used(maxInd+1, false);
	for (unsigned i = 0; i < in.size(); i++)
	{
		used[in[i].end0] = true;
		used[in[i].end1] = true;
	}
	int numNodes = 0;
	for (int i = 0; i < maxInd+1; i++)
	{
		if (used[i])
			numNodes++;
	}
	return numNodes;
}

//returns 1 greater since 0 is counted.
int getMaxNodeInd(const vector<SplitEdge> &in)
{
	int maxInd = 0;
	for (unsigned i = 0; i < in.size(); i++)
	{
		if (in[i].end0 > maxInd)
			maxInd = in[i].end0;
		if (in[i].end1 > maxInd)
			maxInd = in[i].end1;
	}
	return maxInd + 1;
}

vector<bool> getUsedNodes(const vector<SplitEdge> &in)
{
	int maxInd = 0;
	for (unsigned i = 0; i < in.size(); i++)
	{
		if (in[i].end0 > maxInd)
			maxInd = in[i].end0;
		if (in[i].end1 > maxInd)
			maxInd = in[i].end1;
	}
	vector<bool> used(maxInd+1, false);
	for (unsigned i = 0; i < in.size(); i++)
	{
		used[in[i].end0] = true;
		used[in[i].end1] = true;
	}
	return used;
}

void output(const vector<bool> &in)
{
	for (unsigned i = 0; i < in.size(); i++)
	{
		if (in[i])
			cout << "true for i = " << i << endl;
		else
			cout << "false for i = " << i << endl;
	}
}

vector<int> getUsedNodesSet(const vector<SplitEdge> &in)
{
	int maxInd = 0;
	for (unsigned i = 0; i < in.size(); i++)
	{
		if (in[i].end0 > maxInd)
			maxInd = in[i].end0;
		if (in[i].end1 > maxInd)
			maxInd = in[i].end1;
	}
	vector<bool> used(maxInd + 1, false);
	for (unsigned i = 0; i < in.size(); i++)
	{
		used[in[i].end0] = true;
		used[in[i].end1] = true;
	}
	vector<int> out;
	for (unsigned i = 0; i < used.size(); i++)
	{
		if (used[i])
			out.push_back(i);
	}
	return out;
}

vector<vector<int>> sortByCG(const vector<SplitEdge> &g, int s, const vector<vector<int>> &X1)
{
	vector<pair<int,int>> cgs;
	for (unsigned i = 0; i < X1.size(); i++)
	{
		int cg = cG(s, X1[i], g);
		cgs.push_back(pair<int, int>(i, cg));
	}
	vector<int> orderedInds;
	vector<bool> used(cgs.size(), false);
	for (unsigned i = 0; i < cgs.size(); i++)
	{
		int curMin = INT_MAX;
		int ind = -1;
		for (unsigned j = 0; j < cgs.size(); j++)
		{
			if (cgs[j].second < curMin && !used[j])
			{
				curMin = cgs[j].second;
				ind = cgs[j].first;
			}
		}
		used[ind] = true;
		orderedInds.push_back(ind);
	}
	vector<int> rev;
	for (int i = (int)orderedInds.size() - 1; i >= 0; i--)
		rev.push_back(orderedInds[i]);
	vector<vector<int>> out;
	for (unsigned i = 0; i < rev.size(); i++)
		out.push_back(X1[rev[i]]);
	return out;
}

void verifyEquation4(const vector<SplitEdge> &g, int s)
{
	vector<int> X = neighbors(g, s);
	if (X.size() < 2)
	{
		cout << "p < 2" << endl;
		throw logic_error("");
	}
	int acc = 0;
	int max = 0;
	for (unsigned i = 0; i < X.size(); i++)
	{
		int cur = cG(s, X[i], g);
		acc += cur;
		if (cur > max)
			max = cur;
	}
	if (!(max <= acc/2))
	{
		cout << "This X violates Equation 4" << endl;
		cout << "Max is: " << max << endl;
		cout << "Sum is: " << acc << endl;
		cout << "1/2 sum is: " << acc/2 << endl;
		output(X);
		//throw logic_error("");
	}
}

vector<SplitEdge> fixEquation4(const vector<SplitEdge> &g, int s)
{
	vector<SplitEdge> out = g;
	int max = INT_MIN;
	int xMax = 0;
	vector<int> ns = neighbors(g, s);
	for (unsigned i = 0; i < ns.size(); i++)
	{
		if (cG(s, ns[i], g) > max)
		{
			max = cG(s, ns[i], g);
			xMax = ns[i];
		}
	}
	int sum = 0;
	for (unsigned i = 0; i < ns.size(); i++)
	{
		if (ns[i] != xMax)
			sum += cG(s, ns[i], g);
	}
	int toRemove = (cG(s, xMax, g) - sum);
	toRemove = toRemove > 0 ? toRemove : 0;
	if (toRemove != 0)
		cout << "Removing " << toRemove << " edges from " << xMax << " to " << s << " to " << "verify Equation 4" << endl;
	for (unsigned i = 0; i < out.size(); i++)
	{
		if (connects(out[i], xMax, s))
			out[i].weight -= toRemove;
	}
	return out;
}

vector<SplitEdge> modToSplit(const vector<ModEdge> &in)
{
	vector<SplitEdge> out;
	for (unsigned i = 0; i < in.size(); i++)
	{
		SplitEdge e = SplitEdge(in[i].end0, in[i].end1, (int)(round(in[i].weight)), in[i].orig0, in[i].orig1);
		out.push_back(e);
	}
	return out;
}

void output(const vector<SplitEdge> &l)
{
	for (unsigned i = 0; i < l.size(); i++)
	{
		cout << "From " << l[i].end0 << " to " << l[i].end1 << " with weight ";
		cout << l[i].weight << " which was originally from " << l[i].orig0 << " to " << l[i].orig1 << endl;
	}
}

/*vector<SplitEdge> compress(const vector<SplitEdge> &g)
{
	vector<SplitEdge> out;
	for (unsigned i = 0; i < g.size(); i++)
	{
		SplitEdge tester(g[i].end0, g[i].end1, 0, 0, 0);
		if (!isMem(tester, out))
		{
			out.push_back(g[i]);
		}
		else
		{
			int j = indexOfEdge(out, g[i].end0, g[i].end1);
			out[j].weight += g[i].weight;
		}
	}
	return out;
}*/

pair<int, int> getSortedPair(int e1, int e2)
{
	if (e1 > e2)
		return pair<int, int>(e2, e1);
	return pair<int, int>(e1, e2);
}

vector<SplitEdge> compress(const vector<SplitEdge> &g)
{
	unordered_map<uEdge, int> tbl;
	for (unsigned i = 0; i < g.size(); i++)
	{
		uEdge e(g[i].end0, g[i].end1);
		if (tbl.count(e) > 0)
		{
			tbl.at(e) += g[i].weight;
		}
		else
		{
			tbl.emplace(e, g[i].weight);
		}
	}
	vector<SplitEdge> out;
	for (unordered_map<uEdge, int>::iterator e = tbl.begin(); e != tbl.end(); e++)
		out.push_back(SplitEdge(e->first.u, e->first.v, e->second, e->first.u, e->first.v));
	return out;
}

bool verifyMAOrdering(const vector<SplitEdge> &g, const vector<int> &ma)
{
	vector<int> curVs;
	for (int i = 0; i < getNumUsedNodes(g)-1; i++)
	{
		curVs.push_back(ma[i]);
		int iCG = cG(ma[i + 1], curVs, g);
		for (int j = i+1; j < getNumUsedNodes(g); j++)
		{
			int jCG = cG(ma[j], curVs, g);
			if (jCG > iCG)
			{
				cout << "NOT A VALID MA ORDERING" << endl;
				throw logic_error("");
				return false;
			}
		}
	}
	return true;
}



int locInTree(const vector<LPEdge> &t, Edge e)
{
	for (unsigned i = 0; i < t.size(); i++)
	{
		if ((t[i].end0 == e.u && t[i].end1 == e.v) || (t[i].end1 == e.u && t[i].end0 == e.v))
			return i;
	}
	return -1;
}

vector<LPEdge> treesUnion(const vector<vector<Edge>> &ts)
{
	vector<LPEdge> out;
	vector<vector<Edge>> in = ts;
	for (unsigned i = 0; i < in.back().size(); i++)
	{
		out.push_back(LPEdge(in.back()[i].u, in.back()[i].v, 1.0));
	}
	in.pop_back();
	for (unsigned i = 0; i < in.size(); i++)
	{
		for (unsigned j = 0; j < in[i].size(); j++)
		{
			int loc = locInTree(out, in[i][j]);
			if (loc == -1)
				out.push_back(LPEdge(in[i][j].u, in[i][j].v, 1));
			else
				out[loc].weight++;
		}
	}
	return out;
}

int locInTree(const vector<LPEdge> &t, LPEdge e)
{
	for (unsigned i = 0; i < t.size(); i++)
	{
		if ((t[i].end0 == e.end0 && t[i].end1 == e.end1) || (t[i].end1 == e.end0 && t[i].end0 == e.end1))
			return i;
	}
	return -1;
}

void treeError(const vector<LPEdge> &lp, const vector<LPEdge> &tu)
{
	cout << "Running Tree Error:" << endl; 
	if (lp.size() != tu.size())
	{
		cout << "treeError: lp has size " << lp.size() << " but tu has size " << tu.size() << endl;
	}
	for (unsigned i = 0; i < tu.size(); i++)
	{
		int loc = locInTree(lp, tu[i]);
		if (loc == -1)
		{
			cout << "  tu contains an edge from " << tu[i].end0 << " to " << tu[i].end1 << " which isn't in lp at all" << endl;
		}
		else if (tu[i].weight > lp[loc].weight)
		{
			cout << "  the weight of " << tu[i].end0 << " to " << tu[i].end1 << " is " << tu[i].weight << " in tu but only " << lp[loc].weight << " in lp" << endl;
		}
	}
}

SplitEdge edgeToSplitEdge(Edge e)
{
	return SplitEdge(e.u, e.v, 1, e.u, e.v);
}
vector<SplitEdge> edgeToSplitEdge(const vector<Edge> &t)
{
	vector<SplitEdge> out;
	for (unsigned i = 0; i < t.size(); i++)
		out.push_back(edgeToSplitEdge(t[i]));
	return out;
}

SplitEdge joinEdgeToSplitEdge(JoinEdge e)
{
	return SplitEdge(e.u, e.v, 1, e.u, e.v);
}

vector<SplitEdge> joinEdgeToSplitEdge(const vector<JoinEdge> &t)
{
	vector<SplitEdge> out;
	for (unsigned i = 0; i < t.size(); i++)
		out.push_back(joinEdgeToSplitEdge(t[i]));
	return out;
}

vector<int> dijkstra(const vector<SplitEdge> &t, int s)
{
	vector<vector<pair<int, int>>> adjs;
	int maxInd = 0;
	for (unsigned i = 0; i < t.size(); i++)
	{
		if (t[i].end0 > maxInd)
			maxInd = t[i].end0;
		if (t[i].end1 > maxInd)
			maxInd = t[i].end1;
	}
	maxInd++;
	for (int i = 0; i < maxInd; i++)
		adjs.push_back(vector<pair<int, int>>());
	for (unsigned i = 0; i < t.size(); i++)
	{
		adjs[t[i].end0].push_back(pair<int, int>(t[i].end1, t[i].weight));
		adjs[t[i].end1].push_back(pair<int, int>(t[i].end0, t[i].weight));
	}
	vector<bool> visited(maxInd, false);
	int cur = s;
	vector<int> distances(maxInd, INT_MAX);
	distances[cur] = 0;
	int numVisited = 0;

	vector<bool> nnFound(maxInd, false);
	for (unsigned i = 0; i < t.size(); i++)
	{
		nnFound[t[i].end0] = true;
		nnFound[t[i].end1] = true;
	}
	int numNodes = 0;
	for (unsigned i = 0; i < nnFound.size(); i++)
		numNodes += nnFound[i];
	while (numVisited != numNodes)
	{
		for (unsigned i = 0; i < adjs[cur].size(); i++)
		{
			int n = adjs[cur][i].first;
			int w = adjs[cur][i].second;
			distances[n] = min(distances[n], distances[cur] + w);
		}
		visited[cur] = true;
		numVisited++;
		int next = 0;
		int nextDist = INT_MAX;
		for (int i = 0; i < maxInd; i++)
		{
			if (distances[i] < nextDist && !visited[i])
			{
				next = i;
				nextDist = distances[i];
			}
		}
		cur = next;
	}
	return distances;
}

vector<int> dijkstra(const vector<Edge> &t, int s)
{
	return dijkstra(edgeToSplitEdge(t), s);
}

vector<int> dijkstra(const vector<JoinEdge> &t, int s)
{
	return dijkstra(joinEdgeToSplitEdge(t), s);
}

LPEdge eToLP(Edge e)
{
	return LPEdge(e.u, e.v, 1.0);
}

vector<LPEdge> eToLP(const vector<Edge> &t)
{
	vector<LPEdge> out;
	for (unsigned i = 0; i < t.size(); i++)
		out.push_back(eToLP(t[i]));
	return out;
}

vector<vector<LPEdge>> eToLP(const vector<vector<Edge>> &ts)
{
	vector<vector<LPEdge>> out;
	for (unsigned i = 0; i < ts.size(); i++)
		out.push_back(eToLP(ts[i]));
	return out;
}

Edge lpToE(LPEdge e)
{
	return Edge(e.end0, e.end1);
}

vector<Edge> lpToE(const vector<LPEdge> &t)
{
	vector<Edge> out;
	for (unsigned i = 0; i < t.size(); i++)
		out.push_back(lpToE(t[i]));
	return out;
}

vector<LPEdge> genTest()
{
	//corresponds to the following graph:
	/*
	0        x         1

	c        b

	a                  y
	b        c

	3        z         2

	*/
	vector<LPEdge> testSol;
	LPEdge a(0, 3, 1.0);
	LPEdge x(0, 1, 1.0);
	LPEdge y(1, 2, 1.0);
	LPEdge z(2, 3, 1.0);
	testSol.push_back(a);
	testSol.push_back(x);
	testSol.push_back(y);
	testSol.push_back(z);
	return testSol;
}

void weightGraph(vector<ModEdge> &g, vector<double> &ws)
{
	for (unsigned i = 0; (i < g.size() && i < ws.size()); i++)
	{
		g[i].weight = ws[i];
	}
}

//Gives numNodes currently, not originally
int getNumNodes(const vector<ModEdge> &g)
{
	int numNodes = 0;
	for (unsigned i = 0; i < g.size(); i++)
	{
		if (g[i].end0 > numNodes)
			numNodes = g[i].end0;
		if (g[i].end1 > numNodes)
			numNodes = g[i].end1;
	}
	//+1 accounts for off by one from the above method, because 0 is an edge index.
	return numNodes + 1;
}

int isMem(const vector<ModEdge> &g, ModEdge e)
{
	for (unsigned i = 0; i < g.size(); i++)
	{
		if ((g[i].end0 == e.end0 && g[i].end1 == e.end1) ||
			(g[i].end0 == e.end1 && g[i].end1 == e.end0))
		{
			return i;
		}
	}
	return -1;
}

struct link
{
	int node;
	double weight;
	int orig0;
	int orig1;
	link(int n, double w, int o0, int o1) { node = n; weight = w; orig0 = o0; orig1 = o1; }
};

void output(vector<link> l)
{
	for (unsigned i = 0; i < l.size(); i++)
	{
		cout << "Node " << l[i].node << " with weight " << l[i].weight << endl;
	}
}

vector<ModEdge> lpToMod(const vector<LPEdge> &in)
{
	vector<ModEdge> mv;
	for (unsigned i = 0; i < in.size(); i++)
	{
		ModEdge e(in[i].end0, in[i].end1, in[i].weight, in[i].end0, in[i].end1);
		mv.push_back(e);
	}
	return mv;
}

vector<Edge> lpToEdge(const vector<LPEdge> &in)
{
	vector<Edge> es;
	for (unsigned i = 0; i < in.size(); i++)
	{
		Edge e = Edge(in[i].end0, in[i].end1);
		es.push_back(e);
	}
	return es;
}

vector<ModEdge> contractEdge(const vector<ModEdge> &g, ModEdge e)
{
	//cout << "About to contract the edge from " << e.end0 << " to " << e.end1 << " with weight " << e.weight << endl;
	//cout << "The graph before contraction: " << endl;
	//output(g);
	//To contract an edge (u,v) (into u)
	int u = e.end0;
	int v = e.end1;
	double w = e.weight;
	vector<ModEdge> cg;
	vector<ModEdge> unadded;
	//1) All edges which do not have u or v at either end should be added to cg, all others (except (u,v))to unadded
	//2) Collect all nodes connected to u along with their edge weights
	//2) Collect all nodes connected to v along with their edge weights
	//3) Change these two lists into two new lists: nodes connected to either u or to v (and their edge weights),
	//   and nodes connected to both (with those edgeweights summed)
	//4) Add both lists to the graph, with the other end of the edge being u.
	//5) (Later, after other changes) go through the graph and subtract one from all nodes with id > v, crashing if there's
	// an edge with id == u.
	//6) Handle the fact that the weight that needs to be contracted should be multiplied by w at the end (after the cofactor).

	//Step 1:
	for (unsigned i = 0; i < g.size(); i++)
	{
		ModEdge ce(g[i].end0, g[i].end1, g[i].weight, g[i].orig0, g[i].orig1);
		if (g[i].end0 != u && g[i].end0 != v && g[i].end1 != u && g[i].end1 != v)
		{
			cg.push_back(ce);
		}
		else if (!(g[i].end0 == u && g[i].end1 == v))
		{
			unadded.push_back(ce);
		}
	}

	//Step 2:
	vector<ModEdge> toU;
	vector<ModEdge> toV;
	for (unsigned i = 0; i < unadded.size(); i++)
	{
		if (unadded[i].end0 == u || unadded[i].end1 == u)
		{
			toU.push_back(unadded[i]);
		}
		if (unadded[i].end0 == v || unadded[i].end1 == v)
		{
			toV.push_back(unadded[i]);
		}
	}
	//cout << "The edges connected to u: " << endl;
	//output(toU);
	//cout << "The edges connected to v: " << endl;
	//output(toV);

	//Step 3:
	//Add everything from toU to orList.
	//Go through toV. for each edge, if its node is referenced in orList, remove that node
	//from orList, add the proper element to bothList.
	vector<link> orList;
	vector<link> bothList;
	for (unsigned i = 0; i < toU.size(); i++)
	{
		int otherID = 0;
		if (toU[i].end0 != u)
			otherID = toU[i].end0;
		else
			otherID = toU[i].end1;
		link l(otherID, toU[i].weight, toU[i].end0, toU[i].end1);
		orList.push_back(l);
	}
	vector<int> toDelFromOrList;
	for (unsigned i = 0; i < toV.size(); i++)
	{
		int currentID = 0;
		if (toV[i].end0 != v)
			currentID = toV[i].end0;
		else
			currentID = toV[i].end1;
		bool already = false;
		for (unsigned j = 0; j < orList.size(); j++)
		{
			if (orList[j].node == currentID)
			{
				link bL(currentID, orList[j].weight + toV[i].weight, orList[j].orig0, orList[j].orig1);
				//if (orList[j].weight < toV[i].weight)
				if (!should((orList[j].weight / (orList[j].weight + toV[i].weight))))
				{
					bL.orig0 = toV[i].orig0;
					bL.orig1 = toV[i].orig1;
				}
				bothList.push_back(bL);
				toDelFromOrList.push_back(j);
				already = true;
				break;
			}
		}
		if (!already)
		{
			orList.push_back(link(currentID, toV[i].weight, toV[i].orig0, toV[i].orig1));
		}
	}
	vector<link> temp;
	for (unsigned i = 0; i < orList.size(); i++)
	{
		bool found = false;
		for (unsigned j = 0; j < toDelFromOrList.size(); j++)
		{
			if (toDelFromOrList[j] == i)
				found = true;
		}
		if (!found)
		{
			temp.push_back(orList[i]);
		}
	}
	orList = temp;
	//cout << "The nodes connected to u xor v: " << endl;
	//output(orList);
	//cout << "The nodes connected to u and v: " << endl;
	//output(bothList);

	//Step 4: 
	for (unsigned i = 0; i < orList.size(); i++)
	{
		ModEdge cur(u, orList[i].node, orList[i].weight, orList[i].orig0, orList[i].orig1);
		cg.push_back(cur);
	}
	for (unsigned i = 0; i < bothList.size(); i++)
	{
		ModEdge cur(u, bothList[i].node, bothList[i].weight, bothList[i].orig0, bothList[i].orig1);
		cg.push_back(cur);
	}

	/*//Step 5:
	for (unsigned i = 0; i < cg.size(); i++)
	{
	if (cg[i].end0 == v || cg[i].end1 == v)
	throw logic_error("Incorrect contraction");
	if (cg[i].end0 > v)
	cg[i].end0--;
	if (cg[i].end1 > v)
	cg[i].end1--;
	}*/

	//cout << "Contraction completed, the new graph is of size: " << cg.size() << endl;
	//cout << "Here is the contracted graph: " << endl;
	//output(cg);
	return cg;
}

vector<ModEdge> contractLastEdge(const vector<ModEdge> &g)
{
	return contractEdge(g, g.back());
}

vector<ModEdge> xToZ(const vector<ModEdge> &x)
{
	vector<ModEdge> z;
	int numNodes = getNumNodes(x);
	double scalar = ((double)numNodes - 1.0) / ((double)numNodes);
	for (unsigned i = 0; i < x.size(); i++)
	{
		ModEdge z_i = ModEdge(x[i].end0, x[i].end1, x[i].weight*scalar, x[i].orig0, x[i].orig1);
		z.push_back(z_i);
	}
	return z;
}


bool should(double pr)
{
	return (abs(1.0 - pr) <= .001) || (((float)rand()) / ((float)RAND_MAX)) < pr;
}

void output(vector<double> l)
{
	for (unsigned i = 0; i < l.size(); i++)
		cout << l[i] << endl;
}

void output(vector<LPEdge> l)
{
	for (unsigned i = 0; i < l.size(); i++)
	{
		cout << "From " << l[i].end0 << " to " << l[i].end1 << " with weight ";
		cout << l[i].weight << endl;
	}
}

void output(vector<ModEdge> l)
{
	for (unsigned i = 0; i < l.size(); i++)
	{
		cout << "From " << l[i].end0 << " to " << l[i].end1 << " with weight ";
		cout << l[i].weight << " which was originally from " << l[i].orig0 << " to " << l[i].orig1 << endl;
	}
}

//compensates for the edges that have been contracted resulting in holes in the
//vertex ordering. Does not touch any vertex which has become 0-connected from
//edge deletion (although I don't believe this should ever happen).
vector<ModEdge> adjustGraph(const vector<ModEdge> &lGraph, const vector<int> &cNodes)
{
	//First, get number of nodes referenced in the graph.
	int numRefNodes = 0;
	for (unsigned i = 0; i < lGraph.size(); i++)
	{
		if (lGraph[i].end0 > numRefNodes)
			numRefNodes = lGraph[i].end0;
		if (lGraph[i].end1 > numRefNodes)
			numRefNodes = lGraph[i].end1;
	}
	numRefNodes++;


	//Build a vector, which for each node in the graph, maps to the amount
	//by which it should be decremented to remove irrelevant nodes.
	//Nodes which are deleted could be anything, won't affect result.
	//(invariant: cNodes only contains 0-connected vertices)
	vector<int> dec;
	for (int i = 0; i < numRefNodes; i++)
	{
		int cur = 0;
		for (unsigned j = 0; j < cNodes.size(); j++)
		{
			if (cNodes[j] < i)
				cur++;
		}
		dec.push_back(cur);
	}

	//decrement each vertex by the proper amount
	vector<ModEdge> adjGraph;
	for (unsigned i = 0; i < lGraph.size(); i++)
	{
		int u = lGraph[i].end0;
		int v = lGraph[i].end1;
		ModEdge e = ModEdge(u - dec[u], v - dec[v], lGraph[i].weight, lGraph[i].orig0, lGraph[i].orig1);
		adjGraph.push_back(e);
	}
	return adjGraph;
}

//Returns a graph where any disconnected nodes from the input graph have been removed
vector<ModEdge> cleanGraph(const vector<ModEdge> &g)
{
	//get number of nodes
	//create an array mapping nodes to if they are used
	//fill the array
	//change this array to a list of unused nodes.
	//run the adjust graph method on these nodes.
	int numNodes = getNumNodes(g);
	vector<bool> used(numNodes, false);
	for (unsigned i = 0; i < g.size(); i++)
	{
		used[g[i].end0] = true;
		used[g[i].end1] = true;
	}
	vector<int> unused;
	for (int i = 0; i < numNodes; i++)
	{
		if (!used[i])
		{
			unused.push_back(i);
		}
	}
	cout << "Cleaning " << unused.size() << " nodes" << endl;
	return adjustGraph(g, unused);
}

bool edgeConnects(ModEdge e, int u, int v)
{
	return ((e.end0 == u && e.end1 == v) || (e.end1 == u && e.end0 == v));
}

int edgeLoc(const vector<ModEdge> &g, int u, int v)
{
	for (unsigned i = 0; i < g.size(); i++)
	{
		if (edgeConnects(g[i], u, v))
			return i;
	}
	return -1;
}

vector<ModEdge> contractEdgeNew(const vector<ModEdge> &g, ModEdge e)
{
	int u = e.end0;
	int v = e.end1;
	double w = e.weight;
	vector<ModEdge> cg;
	vector<ModEdge> unadded;

	//We will change node v to be node u. So: go through the entire graph, and rename all v's to u. (if you hit u,v, remove it)
	//Then, go back through the graph and see if there's any repeat edges. If so, remove both and replace with a single edge that has the sum of their weights.
	vector<ModEdge> g1;
	for (unsigned i = 0; i < g.size(); i++)
	{
		ModEdge c(g[i].end0, g[i].end1, g[i].weight, g[i].orig0, g[i].orig1);
		if (c.end0 == v)
			c.end0 = u;
		if (c.end1 == v)
			c.end1 = u;
		if (!(c.end0 == u && c.end1 == u))
			g1.push_back(c);
	}
	vector<ModEdge> g2;
	for (unsigned i = 0; i < g1.size(); i++)
	{
		int loc = edgeLoc(g2, g1[i].end0, g1[i].end1);
		if (loc == -1)
			g2.push_back(g1[i]);
		else
			g2[loc].weight += g1[i].weight;
	}
	return g2;
}

vector<string> splitString(string str, char delim)
{
	vector<string> output;
	int pos = 0;
	int len = 0;
	for (unsigned int i = 0; i < str.length(); i++)
	{
		if (str[i] == delim)
		{
			output.push_back(str.substr(pos, len));
			pos = i + 1;
			len = 0;
		}
		len++;
	}
	output.push_back(str.substr(pos, len));
	return output;
}

bool stringIsWhitespace(string str)
{
	for (unsigned int i = 0; i < str.size(); i++)
	{
		if (str[i] != ' ' && str[i] != '\t' && str[i] != '\n' && str[i] != '\0')
			return false;
	}
	return true;
}

string getFileWithoutDirectory(string fName)
{
	string fNoDir = "";
	char delim;
	if (fName.find('/') != string::npos)
		delim = '/';
	else if (fName.find('\\') != string::npos)
		delim = '\\';
	else
		return fName;

	int namePos = 0;
	for (int l = fName.size() - 1; l >= 0; l--)
	{
		if (fName[l] == delim)
			break;
		namePos = l;
	}
	return fName.substr(namePos, std::string::npos);
}

string getDirectoryOfPath(string path)
{
	return path.substr(0, path.find(getFileWithoutDirectory(path)));
}