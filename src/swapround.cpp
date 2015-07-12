#include "swapRound.h"

void print_lpv(vector<LPEdge> v)
{
	for (unsigned i = 0; i < v.size(); i++)
	{
		cout << v[i].end0 << ", " << v[i].end1 << ", w=" << v[i].weight << endl;
	}
}
double sum(const vector<double> &l, int from, int to)
{
	double acc = 0.0;
	for (int i = from; i < to; i++)
		acc += l[i];
	return acc;
}

void removeEdge(vector<LPEdge> &b, LPEdge e)
{
	for (unsigned i = 0; i < b.size(); i++)
	{
		if (b[i].end0 == e.end0 && b[i].end1 == e.end1)
		{
			b.erase(b.begin() + i);
			break;
		}
	}
}

bool edgesEqual(LPEdge e1, LPEdge e2)
{
	return ((e1.end0 == e2.end0 && e1.end1 == e2.end1) ||
		(e1.end0 == e2.end1 && e1.end1 == e2.end0));
}

bool edgesEqualStrict(LPEdge e1, LPEdge e2)
{
	return (e1.end0 == e2.end0 && e1.end1 == e2.end1);
}

bool isMem(LPEdge e, vector<LPEdge> s)
{
	for (unsigned i = 0; i < s.size(); i++)
	{
		if (edgesEqual(e, s[i]))
			return true;
	}
	return false;
}

struct setComparator {
	bool operator() (LPEdge e1, LPEdge e2) {
		/*int e1L = max(e1.end0, e1.end1);
		int e1S = min(e1.end0, e1.end1);
		int e2L = max(e2.end0, e2.end1);
		int e2S = min(e2.end0, e2.end1);
		return (e1S < e2S || (e1S == e2S && e1L < e2L));*/
		return (e1.end0 < e2.end0 || (e1.end0 == e2.end0 && e1.end1 < e2.end1));
	}
}; setComparator setComp;

// Returns s1 \ s2, or equiv. all elements of s1 not in s2.
vector<LPEdge> setDif(vector<LPEdge> s1, vector<LPEdge> s2)
{
	vector<LPEdge> dif;
	for (unsigned i = 0; i < s1.size(); i++)
	{
		if (!isMem(s1[i], s2))
			dif.push_back(s1[i]);
	}
	return dif;
}

vector<LPEdge> setDifFast(const vector<LPEdge> &s1, const vector<LPEdge> &s2)
{
	vector<LPEdge> s1c;
	vector<LPEdge> s2c;
	for (unsigned i = 0; i < s1.size(); i++)
		s1c.push_back(LPEdge(s1[i].end0, s1[i].end1, s1[i].weight));
	for (unsigned i = 0; i < s2.size(); i++)
		s2c.push_back(LPEdge(s2[i].end0, s2[i].end1, s2[i].weight));
	//for (unsigned i = 0; i < s1.size(); i++)
	//	cout << "s1[" << i << "] = " << s1[i].end0 << " " << s1[i].end1 << endl;
	std::sort(s1c.begin(), s1c.end(), setComp);
	//for (unsigned i = 0; i < s1.size(); i++)
	//	cout << "s1[" << i << "] = " << s1[i].end0 << " " << s1[i].end1 << endl;
	std::sort(s2c.begin(), s2c.end(), setComp);
	vector<LPEdge> dif;
	int loc1 = 0;
	int loc2 = 0;
	while (true)
	{
		if (loc1 == s1c.size())
			return dif;
		if (loc2 == s2c.size())
		{
			for (unsigned i = loc1; i < s1c.size(); i++)
				dif.push_back(s1c[i]);
			return dif;
		}
		if (edgesEqualStrict(s1c[loc1], s2c[loc2]))
		{
			loc1++;
			continue;
		}
		if (setComp(s1c[loc1], s2c[loc2]))
		{
			dif.push_back(s1c[loc1]);
			loc1++;
			continue;
		}
		else
		{
			loc2++;
		}
	}

}


int find(int loc, vector<int> &in)
{
	vector<int> toCompress;
	int cur = loc;
	while (true)
	{
		if (in[cur] == cur)
		{
			for (unsigned i = 0; i < toCompress.size(); i++)
				in[toCompress[i]] = cur;
			return cur;
		}
		else
		{
			toCompress.push_back(cur);
			cur = in[cur];
		}
	}
	return -1; // for compiler
}



bool hasCycle2(const vector<LPEdge> &in)
{
	int maxInd = 0;
	for (unsigned i = 0; i < in.size(); i++)
		maxInd = max(maxInd, max(in[i].end0, in[i].end1));
	maxInd++;
	vector<int> owners;
	for (int i = 0; i < maxInd; i++)
		owners.push_back(i);
	
	for (unsigned i = 0; i < in.size(); i++)
	{
		int u = in[i].end0;
		int v = in[i].end1;
		int vOwner = find(v, owners);
		int uOwner = find(u, owners);
		if (vOwner == uOwner)
			return true;
		owners[vOwner] = uOwner;
	}
	return false;
}


bool basesEqual(const vector<LPEdge> &b1, const vector<LPEdge> &b2)
{
	if (b1.size() != b2.size())
		return false;
	/*for (unsigned i = 0; i < b1.size(); i++)
	{
		if (!isMem(b1[i], b2))
			return false;
	}
	return true;*/
	vector<LPEdge> b1c;
	vector<LPEdge> b2c;
	for (unsigned i = 0; i < b1.size(); i++)
		b1c.push_back(LPEdge(b1[i].end0, b1[i].end1, b1[i].weight));
	for (unsigned i = 0; i < b2.size(); i++)
		b2c.push_back(LPEdge(b2[i].end0, b2[i].end1, b2[i].weight));
	std::sort(b1c.begin(), b1c.end(), setComp);
	std::sort(b2c.begin(), b2c.end(), setComp);
	for (unsigned i = 0; i < b1.size(); i++)
	{
		if (!edgesEqualStrict(b1c[i], b2c[i]))
			return false;
	}
	return true;
}

int findNoCycle(vector<LPEdge> b1, vector<LPEdge> b2, vector<LPEdge> dif, LPEdge i)
{
	b2.push_back(i);
	removeEdge(b1, i);

	for (unsigned k = 0; k < dif.size(); k++)
	{
		LPEdge j = dif[k];
		removeEdge(b2, j);
		b1.push_back(j);
		//cout << "Starting cyc rec" << endl;
		if (!hasCycle2(b1) && !hasCycle2(b2))
			return k;
		//cout << "Done with 2 cyc recs" << endl;
		removeEdge(b1, j);
		b2.push_back(j);
	}
	throw logic_error("Invariant broken");
}

vector<vector<int>> makeAdjs(const vector<LPEdge> &tree)
{
	vector<vector<int>> adjs;
	adjs.reserve(tree.size() + 1);
	for (unsigned int j = 0; j < tree.size() + 1; j++)
		adjs.push_back(vector<int>());
	for (unsigned j = 0; j < tree.size(); j++)
	{
		adjs[tree[j].end0].push_back(tree[j].end1);
		adjs[tree[j].end1].push_back(tree[j].end0);
	}
	return adjs;
}

//returns the unique s-t path in tree from i.end0 to i.end1
vector<LPEdge> getSTPath(const vector<LPEdge> &tree, const LPEdge &i, list<int> *adjs)
{
	//cout << "getSTPath: " << endl;
	/*vector<vector<int>> adjs;
	adjs.reserve(tree.size()+1);
	for (int j = 0; j < tree.size() + 1; j++)
		adjs.push_back(vector<int>());

	for (unsigned j = 0; j < tree.size(); j++)
	{
		adjs[tree[j].end0].push_back(tree[j].end1);
		adjs[tree[j].end1].push_back(tree[j].end0);
	}*/
	vector<int> prevs(tree.size() + 1, 0);
	//prevs.reserve(tree.size()+1);
	//for (unsigned j = 0; j < tree.size()+1; j++)
	//	prevs.push_back(0);
	//adjs.reserve(tree.size());

	int s = i.end0;
	int t = i.end1;
	queue<int> Q;
	vector<bool> discovered(tree.size() + 1, false);
	Q.push(s);
	discovered[s] = true;
	prevs[s] = s;
	while (!Q.empty())
	{
		int v = Q.front();
		Q.pop();
		for (auto it = adjs[v].begin(); it != adjs[v].end(); it++)
		{
			//cout << "Iterating on adjs[" << v << "]" << endl;
			int w = *it;
			if (!discovered[w])
			{
				Q.push(w);
				prevs[w] = v;
				discovered[w] = true;
			}
		}
	}
	int cur = t;
	vector<int> path;
	while (cur != s)
	{
		path.push_back(cur);
		cur = prevs[cur];
	}
	path.push_back(cur);
	vector<LPEdge> out;
	for (unsigned int j = 1; j < path.size(); j++)
	{

		int e0 = min(path[j - 1], path[j]);
		int e1 = max(path[j-1], path[j]);
		out.push_back(LPEdge(e0, e1, 1.0));
	}
	//cout << "Path: " << endl;
	//print_lpv(out);
	return out;

}

typedef vector<pair<int, int>> uf;

uf uf_make(unsigned int size)
{
	uf data;
	for (unsigned int i = 0; i < size; i++)
		data.push_back(pair<int, int>(i, 0));

	//for (unsigned j = 0; j < tree.size(); j++)
	//{
		//if (data.count(tree[j].end0) == 0)
		//	data.emplace(tree[j].end0, pair<int, int>(tree[j].end0, 1));
		//if (data.count(tree[j].end1) == 0)
		//	data.emplace(tree[j].end1, pair<int, int>(tree[j].end1, 1));
	//}
	/*int max_node = 0;
	for (unsigned i = 0; i < tree.size(); i++)
		max_node = max(max_node, max(tree[i].end0, tree[i].end1));
	for (int i = 0; i < max_node + 1; i++)
		data.push_back(pair<int, int>(i, 1));*/
	return data;
}

/*void uf_print(uf &u)
{
	for (auto it = u.begin(); it != u.end(); it++)
	{
		cout << "The owner of " << it->first << " is " << it->second.first << ". The size of the set is " << it->second.second << endl;
	}
}*/

int uf_find(int loc, uf &in)
{
	if (in[loc].first == loc)
		return loc;
	in[loc].first = uf_find(in[loc].first, in);
	return in[loc].first;
}

void uf_union(int a, int b, uf &in)
{
	a = uf_find(a, in);
	b = uf_find(b, in);
	if (a == b)
		return;
	if (in[a].second < in[b].second)
		in[a].first = b;
	else if (in[a].second > in[b].second)
		in[b].first = a;
	else
	{
		in[a].first = b;
		in[b].second++;
	}
}

bool uf_connected(int a, int b, uf &in)
{
	int aComp = uf_find(a, in);
	int bComp = uf_find(b, in);
	return aComp == bComp;
}

//Returns the unique cycle in the graph which is base2+i, where base2 is a tree.
/*vector<LPEdge> getCycle(const vector<LPEdge> &base2, const LPEdge &i)
{
	//cout << "get Cycle: " << endl;
	//cout << "base2: " << endl;
	//print_lpv(base2);
	vector<vector<int>> adjs = makeAdjs(base2);
	vector<LPEdge> path = getSTPath(base2, i, adjs);
	//cout << "path: " << endl;
	//print_lpv(path);
	path.push_back(LPEdge(i.end0, i.end1, i.weight));
	//cout << "i: " << i.end0 << ", " << i.end1 << endl;
	return path;
}*/
/*
bool hasCycle2(const vector<LPEdge> &in)
{
	int maxInd = 0;
	for (unsigned i = 0; i < in.size(); i++)
		maxInd = max(maxInd, max(in[i].end0, in[i].end1));
	maxInd++;
	vector<int> owners;
	for (int i = 0; i < maxInd; i++)
		owners.push_back(i);

	for (unsigned i = 0; i < in.size(); i++)
	{
		int u = in[i].end0;
		int v = in[i].end1;
		int vOwner = find(v, owners);
		int uOwner = find(u, owners);
		if (vOwner == uOwner)
			return true;
		owners[vOwner] = uOwner;
	}
	return false;
}*/
//Return the set of edges which cross the connected components of base1 - i

vector<LPEdge> getEdgesInCut(const vector<LPEdge> &base1,const LPEdge &i, const vector<LPEdge> &additionalEdges)
{
	//cout << "get edges in cut: " << endl;
	//cout << "base1: " << endl;
	//print_lpv(base1);

	//removeEdge(base1, i);
	uf uf_struct = uf_make(base1.size()+1);
	//cout << "The initial union-find structure: " << endl;
	//uf_print(uf_struct);
	for (unsigned j = 0; j < base1.size(); j++)
	{
		if (!edgesEqual(base1[j], i))
		{
			uf_union(base1[j].end0, base1[j].end1, uf_struct);
		}
	}
	vector<LPEdge> cut;
	for (unsigned j = 0; j < base1.size(); j++)
	{
		if (!uf_connected(base1[j].end0, base1[j].end1, uf_struct))
			cut.push_back(LPEdge(base1[j].end0, base1[j].end1, base1[j].weight));
	}
	//cout << "The final union-find structure: " << endl;
	//uf_print(uf_struct);
	///for (int j = 0; j < )
	//cout << "Cut without additions: " << endl << endl;
	//print_lpv(cut);

	for (unsigned j = 0; j < additionalEdges.size(); j++)
	{
		if (!uf_connected(additionalEdges[j].end0, additionalEdges[j].end1, uf_struct))
			cut.push_back(LPEdge(additionalEdges[j].end0, additionalEdges[j].end1, additionalEdges[j].weight));
	}
	return cut;
}

string str_of_lpe(const LPEdge &e)
{
	return to_string(e.end0) + "-" + to_string(e.end1);
}

vector<LPEdge> getIntersection(const vector<LPEdge> &setA,const vector<LPEdge> &setB)
{
	unordered_map<string, bool> set_hash;
	vector<LPEdge> intersection;
	for (unsigned i = 0; i < setA.size(); i++)
	{
		if (set_hash.count(str_of_lpe(setA[i])) > 0)
		{
			cout << "Error! An edge appeared twice in setA." << endl;
		}
		set_hash.emplace(str_of_lpe(setA[i]), true);
	}
	for (unsigned i = 0; i < setB.size(); i++)
	{
		if (set_hash.count(str_of_lpe(setB[i])) > 0)
		{
			intersection.push_back(LPEdge(setB[i].end0, setB[i].end1, setB[i].weight));
			//cout << "Adding the following Edge to the intersection: " << setB[i].end0 << ", " << setB[i].end1 << endl;
		}
	}
	return intersection;
}

LPEdge getEltFromDif(const vector<LPEdge> &s1, const vector<LPEdge> &s2)
{
	unordered_map<string, bool> set_hash;
	for (unsigned i = 0; i < s2.size(); i++)
		set_hash.emplace(str_of_lpe(s2[i]), true);
	for (unsigned i = 0; i < s1.size(); i++)
	{
		if (set_hash.count(str_of_lpe(s1[i])) == 0)
			return LPEdge(s1[i].end0, s1[i].end1, s1[i].weight);
	}
	return LPEdge();
}

struct edge_hash_data {
	int b;
	int a;
	int status;
	edge_hash_data()
	{
		a = 0;
		b = 0;
		status = 0;
	}
	edge_hash_data(int a, int b, int status)
	{
		this->a = a;
		this->b = b;
		this->status = status;
	}
};

/*
vector<LPEdge> mergeBasesExperiment(double beta1, vector<LPEdge> &base1, double beta2, vector<LPEdge> &base2)
{
	//Generate A hashtable which given an edge, states whether it is in a, b, or both
	const int IN_NEITHER = 0;
	const int IN_A = 1;
	const int IN_B = 2;
	const int IN_BOTH = 3;
	std::unordered_map<string, edge_hash_data> set_hash;
	int in_a_count = base1.size();
	int in_b_count = base2.size();
	int in_both_count = 0;
	for (int i = 0; i < base1.size(); i++)
		set_hash.emplace(str_of_lpe(base1[i]), edge_hash_data(i, -1, IN_A));
	for (int i = 0; i < base2.size(); i++)
	{
		if (set_hash.count(str_of_lpe(base2[i])) > 0)
		{
			in_both_count++;
			set_hash[str_of_lpe(base2[i])].b = i;
			set_hash[str_of_lpe(base2[i])].status = IN_BOTH;
		}
		else
			set_hash.emplace(str_of_lpe(base2[i]), edge_hash_data(-1, i, IN_B));
	}
	//

	

	while (in_both_count != in_a_count || in_both_count != in_b_count) //while !basesEqual
	{

		// Generate the set differences
		vector<LPEdge> BminA;
		vector<LPEdge> AminB;
		for (auto it = set_hash.begin(); it != set_hash.end(); it++)
		{
			if (it->second.status == IN_A)
				AminB.push_back(base1[it->second.a]);
			if (it->second.status == IN_B)
				BminA.push_back(base2[it->second.b]);
		}
		//

		LPEdge i = AminB[0];
		string i_hash = str_of_lpe(i);
		//Find an elt of BminA s.t. A - i + j has no cycle and 
		//B - j + i has no cycle
		for (unsigned k = 0; k < BminA.size(); k++)
		{
			LPEdge j = BminA[k];
			string j_hash = str_of_lpe(j);
			if (set_hash[j_hash].status == 
		}

	}

} */

vector<LPEdge> mergeBases(double beta1, vector<LPEdge> base1, double beta2, vector<LPEdge> base2)
{
	while (!basesEqual(base1, base2))
	{
		//LPEdge i = setDif(base1, base2)[0];
		LPEdge i = setDifFast(base1, base2)[0];
		//vector<LPEdge> jDif = setDif(base2, base1);
		vector<LPEdge> jDif = setDifFast(base2, base1);
		int jLoc = findNoCycle(base1, base2, jDif, i);
		LPEdge j = jDif[jLoc];
		double prob = beta1 / (beta1 + beta2);
		int r = rand() % 1001;
		double dr = (double)r / (1000.0);
		if (dr <= prob)
		{
			removeEdge(base2, j);
			base2.push_back(i);
		}
		else
		{
			removeEdge(base1, i);
			base1.push_back(j);
		}
	}
	return base1;
}

void makeAdjs(const vector<LPEdge> &tree, list<int> *adjs)
{
	for (unsigned int j = 0; j < tree.size() + 1; j++)
		adjs[j] = list<int>();

	for (unsigned int j = 0; j < tree.size(); j++)
	{
		adjs[tree[j].end0].push_back(tree[j].end1);
		adjs[tree[j].end1].push_back(tree[j].end0);
	}
	return;
}

vector<LPEdge> mergeBasesFast(double beta1, vector<LPEdge> base1, double beta2, vector<LPEdge> base2, mt19937 &rg)
{
	//vector<vector<int>> adjsB2 = makeAdjs(base2);
	list<int> *adjs = new list<int>[base2.size() + 1];
	makeAdjs(base2, adjs);
	uniform_real_distribution<double> distr(0.0f, 1.0f);
	/*for (int i = 0; i < base2.size() + 1; i++)
	{
		cout << "Adjacent to " << i << ": ";
		for (auto it = adjs[i].begin(); it != adjs[i].end(); it++)
			cout << *it << ", ";
		cout << endl;
	}*/
	while (!basesEqual(base1, base2))
	{
		//LPEdge i = setDifFast(base1, base2)[0];
		vector<LPEdge> AminB = setDifFast(base1, base2);
		LPEdge i = LPEdge(AminB.back().end0, AminB.back().end1, AminB.back().weight);
		//LPEdge i = getEltFromDif(base1, base2);
		vector<LPEdge> intersect = getIntersection(getSTPath(base2, i, adjs), getEdgesInCut(base1, i, base2));
		//cout << "There are " << intersect.size() << " items in the intersection." << endl;
		//if (intersect.size() != 1)
		//{
			//cout << "Error! Bug in swapRound. There are " << intersect.size() << " items in the intersection" << endl;
			//throw invalid_argument("");
		//}
		LPEdge j = LPEdge(intersect.back().end0, intersect.back().end1, intersect.back().weight);


		double prob = beta1 / (beta1 + beta2);
		//int r = rand() % 1001;
		//double dr = (double)r / (1000.0);
		double dr = distr(rg);
		//cout << "Prob: " << prob << endl;
		if (dr <= prob)
		{
			//cout << "DOING" << endl;
			removeEdge(base2, j);
			base2.push_back(i);
			bool founda = false;
			bool foundb = false;
			for (auto it = adjs[j.end0].begin(); it != adjs[j.end0].end(); it++)
			{
				if (*it == j.end1)
				{
					adjs[j.end0].erase(it);
					founda = true;
					break;
				}
			}
			for (auto it = adjs[j.end1].begin(); it != adjs[j.end1].end(); it++)
			{
				if (*it == j.end0)
				{
					adjs[j.end1].erase(it);
					foundb = true;
					break;
				}
			}
			if (!founda || !foundb)
				cout << "Warning; couldn't erase edge fully from adj list" << endl;
			adjs[i.end0].push_back(i.end1);
			adjs[i.end1].push_back(i.end0);
		}
		else
		{
			removeEdge(base1, i);
			base1.push_back(j);
		}
		//string dummy;
		//cin >> dummy;

	}
	delete[] adjs;
	return base1;
}

vector<vector<LPEdge>> swapRound(const vector<vector<LPEdge>> &trees, const vector<double> &betas)
{
	//srand((unsigned)time(NULL));
	system_clock::rep seed = system_clock::now().time_since_epoch().count();
	mt19937 mersenne_generator = mt19937((unsigned int)seed);
	vector<vector<LPEdge>> ts;
	//sort all the tries such than e.end0 < e.end1 for all edges;
	for (unsigned i = 0; i < trees.size(); i++)
	{
		vector<LPEdge> t;
		for (unsigned j = 0; j < trees[i].size(); j++)
		{
			if (trees[i][j].end0 == trees[i][j].end1)
				throw invalid_argument("");
			if (trees[i][j].end0 < trees[i][j].end1)
				t.push_back(LPEdge(trees[i][j].end0, trees[i][j].end1, trees[i][j].weight));
			else
				t.push_back(LPEdge(trees[i][j].end1, trees[i][j].end0, trees[i][j].weight));
		}
		ts.push_back(t);
	}
	//
	vector<vector<LPEdge>> C;
	C.push_back(ts[0]);
	for (int i = 0; i <= (int)ts.size() - 2; i++)
	{
		double nextBeta = sum(betas, 0, i + 1);
		//cout << "Next beta: " << nextBeta << endl;
		vector<LPEdge> CNext = mergeBasesFast(nextBeta, C.back(), betas[i + 1], ts[i + 1], mersenne_generator);
		C.push_back(CNext);
	}
	return C;
}