
#include "edge_join.h"
#include "edge_split.h"
#include "print.h"

splitting logToSplits(historyIndex &h)
{
	vector<vector<JoinEdge>> firstOut;
	//cout << "This is the first entry in the log: " << endl;
	h.loc = h.first;
	//output(*h.loc);
	//cout << "This is the second entry in the log: " << endl;
	//output(*h.loc->next);
	int maxInd = 0;
	while (h.loc->next != nullptr)
	{
		h.loc = h.loc->next;
		if (h.loc->s > maxInd)
			maxInd = h.loc->s;
	}
	for (int i = 0; i < maxInd + 1; i++)
	{
		firstOut.push_back(vector<JoinEdge>());
	}
	vector<int> order;
	h.loc = h.first;
	while (h.loc->next != nullptr)
	{
		h.loc = h.loc->next;
		int s = h.loc->s;
		int u = h.loc->u;
		int v = h.loc->v;
		for (int i = 0; i < h.loc->multiplicity; i++)
		{
			firstOut[s].push_back(JoinEdge(u, v));
		}
		order.push_back(s);
	}
	//cout << "This is the first entry in order: " << endl;
	//cout << "This is the first entry in firstOut:" << endl;
	//output(firstOut[0]);

	for (unsigned i = 0; i < order.size()-1; i++)
	{
		if (order[i] > order[i + 1])
		{
			cout << "Error: Node " << order[i] << " was split before node " << order[i + 1] << ". Quitting." << endl;
			throw logic_error("");
		}
	}
	vector<int> compressedOrder;
	compressedOrder.push_back(order[0]);
	for (unsigned i = 1; i < order.size(); i++)
	{
		if (order[i-1] != order[i])
			compressedOrder.push_back(order[i]);
	}
	if (compressedOrder.back() != order.back())
		compressedOrder.push_back(order.back());
	vector <pair<int, vector<JoinEdge>>> out;
	for (unsigned i = 0; i < compressedOrder.size(); i++)
	{
		int s = compressedOrder[i];
		pair<int, vector<JoinEdge>> cur(s, firstOut[s]);
		out.push_back(cur);
	}
	return out;
}


vector<Edge> joinEtoE(const vector<JoinEdge> &t)
{
	vector<Edge> out;
	for (unsigned i = 0; i < t.size(); i++)
	{
		Edge e(t[i].u, t[i].v);
		out.push_back(e);
	}
	return out;
}

vector<vector<Edge>> joinEForestToEForest(const vector <vector<JoinEdge>> &ts)
{
	vector<vector<Edge>> out;
	for (unsigned i = 0; i < ts.size(); i++)
	{
		vector<Edge> t = joinEtoE(ts[i]);
		out.push_back(t);
	}
	return out;
}

pair<vector<vector<JoinEdge>>,vector<JoinEdge>> createTrees(const vector<SplitEdge> &g, const splitting &ss, int k)
{
	for (unsigned i = 0; i < ss.size(); i++)
	{
		if (ss[i].second.size() != k)
		{
			cout << "Error: incorrect number of splitting operations; joining will fail. Quitting." << endl;
			cout << "K is " << k << " but there are " << ss[i].second.size() << " splittings at k." << endl;
			throw logic_error("");
		}
	}
	if (g.size() != 1 || g[0].weight != 2 * k)
	{
		cout << "Error with input graph. Quitting." << endl;
		throw logic_error("");
	}
	vector<vector<JoinEdge>> trees;
	for (int i = 0; i < k; i++)
	{
		vector<JoinEdge> tree;
		tree.push_back(JoinEdge(g[0].end0, g[0].end1));
		trees.push_back(tree);
	}
	vector<JoinEdge> extras;
	for (int i = k; i < 2 * k; i++)
	{
		extras.push_back(JoinEdge(g[0].end0, g[0].end1));
	}
	//cout << "This was the graph after the splitting:" << endl;
	//output(g);
	//cout << "This is the initial set of spanning trees:" << endl;
	//output(trees, "Tree");
	//cout << "This is the initial set of extra Edges: " << endl;
	//output(extras);
	//cout << "K is " << k << " and 2K is " << 2 * k << endl;
	return pair<vector<vector<JoinEdge>>, vector<JoinEdge>>(trees, extras);
}

bool isMem(JoinEdge what, const vector<JoinEdge> &of)
{
	for (unsigned i = 0; i < of.size(); i++)
	{
		if (what == of[i])
			return true;
	}
	return false;
}

int numEdgesSplitFromS(vector<JoinEdge> tree, vector<JoinEdge> splitsAtS)
{
	int acc = 0;
	for (unsigned i = 0; i < tree.size(); i++)
	{
		if (isMem(tree[i], splitsAtS))
			acc++;
	}
	return acc;
}

int indexOfEdge(JoinEdge e, const vector<JoinEdge> &t)
{
	for (unsigned i = 0; i < t.size(); i++)
	{
		if (e == t[i])
			return i;
	}
	return -1;
}

int indexOfEdge(JoinEdge e, const vector<JoinEdge> &t, const vector<bool> &deleg)
{
	for (unsigned i = 0; i < t.size(); i++)
	{
		if (e == t[i] && !deleg[i])
			return i;
	}
	return -1;
}

vector<JoinEdge> EdgesInSplitting(const vector<JoinEdge> &tree, const vector<JoinEdge> &splitsAtS, vector<bool> &deleg)
{
	vector<JoinEdge> acc;
	for (unsigned i = 0; i < tree.size(); i++)
	{
		if (isMem(tree[i], splitsAtS))
		{
			int ind = indexOfEdge(tree[i], splitsAtS, deleg);
			if (ind != -1)
			{
				acc.push_back(tree[i]);
				deleg[ind] = true;
			}
		}

	}
	return acc;
}

//The edge (u,v) is lifted back at s by removing (u,v) from t
//and adding (u,s) and (v,s)
void fullyLiftBack(JoinEdge e, int s, vector<JoinEdge> &t)
{

	t.erase(t.begin() + indexOfEdge(e, t));
	t.push_back(JoinEdge(e.u, s));
	t.push_back(JoinEdge(e.v, s));
}

bool isALeaf(int v, const vector<JoinEdge> &t)
{
	int acc = 0;
	for (unsigned i = 0; i < t.size(); i++)
	{
		if (t[i].u == v)
			acc++;
		if (t[i].v == v)
			acc++;
	}
	return (acc == 1);
}

int findEdgeWithTwoLeaves(const vector<JoinEdge> &tree, const vector<JoinEdge> &X)
{
	for (unsigned i = 0; i < X.size(); i++)
	{
		if (isALeaf(X[i].u, tree) && (isALeaf(X[i].v, tree)))
			return i;
	}
	return -1;
}

int locInTree(JoinEdge e, const vector<JoinEdge> &t)
{
	for (unsigned i = 0; i < t.size(); i++)
	{
		if ((t[i].u == e.u && t[i].v == e.v) || (t[i].v == e.u && t[i].u == e.v))
			return i;
	}
	return -1;
}

void rejoinNode(splitting &ss, int k, vector<JoinEdge> &extras, vector<vector<JoinEdge>> &trees)
{
	//cout << "Rejoining node " << ss.back().first << endl;
	int s = ss.back().first;
	vector<JoinEdge> splitsAtS = ss.back().second;
	ss.pop_back();
	if (splitsAtS.size() != k)
	{
		cout << "Incorrect number of splits done at node " << s << ". Quitting." << endl;
		cout << "K is " << k << " and the number of edges is " << splitsAtS.size() << endl;
		output(splitsAtS);
		throw logic_error("");
	}
	if (extras.size() != k)
	{
		cout << "Incorrect number of extra edges. Quitting." << endl;
		throw logic_error("");
	}

	vector<bool> delegatedInSplits(splitsAtS.size(), false);

	//Now consider some edge-disjoint spanning tree:
	vector<vector<JoinEdge>> crossoverEdges;

	for (unsigned i = 0; i < trees.size(); i++)
	{
		vector<JoinEdge> tree = trees[i];
		vector<JoinEdge> X = EdgesInSplitting(tree, splitsAtS, delegatedInSplits);
		crossoverEdges.push_back(X);
	}

	vector<JoinEdge> extraCrossovers = EdgesInSplitting(extras, splitsAtS, delegatedInSplits);

	int totalCrossovers = 0;
	for (int i = 0; i < k; i++)
	{
		totalCrossovers += crossoverEdges[i].size();
	}
	totalCrossovers += extraCrossovers.size();

	if (totalCrossovers != k)
	{
		cout << "Error: incorrect number of crossover edges!" << endl;
		cout << "K was " << k << endl;
		cout << "total crossovers is: " << totalCrossovers << endl;
		throw logic_error("");
	}

	//cout << "These are the sets of crossover edges: " << endl;
	//output(crossoverEdges, "Crossovers for Tree");

	//cout << "These are the crossovers in the extra edges: " << endl;
	//output(extraCrossovers);

	//store the edges to join up trees where X is empty:
	vector<JoinEdge> currentExtras;
	//

	//handle each of the extra crossovers:
	vector<bool> delegInExtras(extras.size(), false);
	for (unsigned i = 0; i < extraCrossovers.size(); i++)
	{
		int e = indexOfEdge(extraCrossovers[i], extras, delegInExtras);
		delegInExtras[e] = true;
		int u = extras[e].u;
		extras[e].u = s;
		currentExtras.push_back(JoinEdge(s, u));
	}

	vector<bool> connected(trees.size(), false);
	//do the first tree:
	for (unsigned i = 0; i < trees.size(); i++)
	{
		if (crossoverEdges[i].size() != 0)
		{
			//cout << "Size of X for tree " << i << " is " << crossoverEdges[i].size() << endl;
			connected[i] = true;
			//Handle first of X edges:
			//one of the edges MAY have both ends as leaves in its tree
			int doubleE = findEdgeWithTwoLeaves(trees[i], crossoverEdges[i]);
			//otherwise just pick one arbitrarily
			doubleE = doubleE == -1 ? 0 : doubleE;
			
			//cout << endl;
			//cout << "Removing edge from " << crossoverEdges[i][doubleE].u << " to " << crossoverEdges[i][doubleE].v << " in tree #" << i << endl;
			//cout << "Adding edge from " << crossoverEdges[i][doubleE].u << " to " << s << " in tree #" << i << endl;
			//cout << "Adding edge from " << crossoverEdges[i][doubleE].v << " to " << s << " in tree #" << i << endl;
			fullyLiftBack(crossoverEdges[i][doubleE], s, trees[i]);
			crossoverEdges[i].erase(crossoverEdges[i].begin() + doubleE);
			//
			//Handle remaining X-1 edges:
			//if (crossoverEdges[i].size() != 0)
			//	cout << "Now handling the remaining " << crossoverEdges[i].size() << " edges for tree " << i << endl;
			vector<int> distances;
			if (crossoverEdges[i].size() != 0)
				distances = dijkstra(trees[i], s);
			for (unsigned j = 0; j < crossoverEdges[i].size(); j++)
			{
				int u = crossoverEdges[i][j].u;
				int v = crossoverEdges[i][j].v;

				//at most one of u,v is a leaf in trees[0]:
				if (distances[u] > distances[v])
				{
					//cout << "Adding edge from " << s << " to " << u << " in tree #" << i << endl;
					trees[i].push_back(JoinEdge(s, u));
					//
					//cout << "Erasing edge from " << u << " to " << v << " in tree #" << i << endl;
					int loc = locInTree(JoinEdge(u, v), trees[i]);
					trees[i].erase(trees[i].begin() + loc);
					//
					currentExtras.push_back(JoinEdge(s, v));
				}
				else
				{
					//cout << "Adding edge from " << s << " to " << v << " to tree #" << i << endl;
					trees[i].push_back(JoinEdge(s, v));
					//
					//cout << "Erasing edge from " << u << " to " << v << " in tree #" << i << endl;
					int loc = locInTree(JoinEdge(u, v), trees[i]);
					trees[i].erase(trees[i].begin() + loc);
					//
					currentExtras.push_back(JoinEdge(s, u));
				}
			}
			//cout << endl;
		}
	}

	//now go back through and use the edges in currentExtras to connect s in the trees with no split edges
	for (unsigned i = 0; i < trees.size(); i++)
	{
		if (!connected[i])
		{
			int origSize = trees[i].size();
			bool found = false;
			//Not 100% sure if this check actually helps
			
			/*for (unsigned j = 0; j < currentExtras.size(); j++)
			{
				int u = currentExtras[j].u;
				int v = currentExtras[j].v;
				if (isALeaf(u, trees[i]) || isALeaf(v, trees[i]))
				{
					trees[i].push_back(JoinEdge(u, v));
					currentExtras.erase(currentExtras.begin() + j);
					found = true;
					break;
				}
			}*/
			
			if (!found)
			{
				int u = currentExtras[0].u;
				int v = currentExtras[0].v;
				//cout << "Adding edge from " << u << " to " << v << " to tree #" << i << endl;
				trees[i].push_back(JoinEdge(u, v));
				currentExtras.erase(currentExtras.begin());
			}
			if (trees[i].size() != origSize + 1)
			{
				cout << "Error: this tree was not connected to s." << endl;
				throw logic_error("");
			}
		}
	}
	if (currentExtras.size() != 0)
	{
		cout << "Error: some tree must not have been fully connected because there are more than k leftover extras." << endl;
		throw logic_error("");
	}

}

void joinEdges(splitting &ss, int k, vector<JoinEdge> &extras, vector<vector<JoinEdge>> &trees)
{
	cout << "Beginning Edge Joining...";
	int todo = ss.size();
	for (int i = 0; i < todo; i++)
	{
		//cout << "Rejoining a node (This is rejoining #" << i << "): " << endl;
		rejoinNode(ss, k, extras, trees);
		//cout << "These are the trees after rejoining #" << i << endl;
		//output(trees, "Tree");
	}
	cout << "Done." << endl;
	//cout << "This is what is left is extras:" << endl;
	//output(extras);
}