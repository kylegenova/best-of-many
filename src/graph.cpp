//Contains the graph implementations defined in graph.h
//Kyle Genova
//#include "stdafx.h"
#include "graph.h"
#include "util.h"

//Pretty print for .tsp info data
string stringAfterMatch(string source, string match)
{
	if (!source.compare(0, (match.length()), match))
		return (source.substr((source.find_first_of(':')) + 2) + "\n");
	return "";
}

Graph::Graph(vector<Point> &nodeList)
{
	numNodes = (int)nodeList.size();
	for (int i = 0; i < numNodes; i++)
		nodes.push_back(nodeList[i]);
}

double Graph::getTreeCost(const vector<LPEdge> &t)
{
	double acc = 0;
	for (unsigned i = 0; i < t.size(); i++)
		acc += t[i].weight * (this->distance(t[i].end0, t[i].end1));
	return acc;
}

int findOwner(vector<int> &owners, int i)
{
	int loc = i;
	vector<int> toCompress;
	while (owners[loc] != loc)
	{
		toCompress.push_back(loc);
		loc = owners[loc];
	}
	for (unsigned i = 0; i < toCompress.size(); i++)
		owners[toCompress[i]] = loc;
	return loc;
}
//Use Union-Find structure, oand assumes that every node in the graph is somewhere in edge set.
vector<vector<Edge>> splitIntoConnectedComponents(const vector<Edge> &in)
{
	int numNodes = 0;
	for (unsigned i = 0; i < in.size(); i++)
		numNodes = max(numNodes, max(in[i].u, in[i].v));
	numNodes++;
	cout << "NumNodes: " << numNodes << endl;
	int numEdges = in.size();
	cout << "Num Edges: " << numEdges << endl;
	vector<int> owners = vector<int>(numNodes, 0);
	cout << "Owners size: " << owners.size() << endl;
	for (int i = 0; i < numNodes; i++)
		owners[i] = i;
	for (unsigned i = 0; i < in.size(); i++)
	{
		int curOwner1 = findOwner(owners, in[i].u);
		int curOwner2 = findOwner(owners, in[i].v);
		if (curOwner1 < curOwner2)
			owners[curOwner2] = curOwner1;
		else
			owners[curOwner1] = curOwner2;
	}
	vector<vector<Edge>> roughCCs = vector<vector<Edge>>();
	for (int i = 0; i < numNodes; i++)
		roughCCs.push_back(vector<Edge>());
	for (unsigned i = 0; i < in.size(); i++)
	{
		//cout << "i is " << i << endl;
		int o1 = findOwner(owners, in[i].u);
		int o2 = findOwner(owners, in[i].v);
		if (o1 != o2)
		{
			//Something has gone horribly wrong
			cout << "Error: Bug in Union-Find structure!" << endl;
			throw runtime_error("");
		}
		roughCCs[o1].push_back(in[i]);
	}
	vector<vector<Edge>> out;
	for (unsigned i = 0; i < roughCCs.size(); i++)
	{
		if (roughCCs[i].size() != 0)
			out.push_back(roughCCs[i]);
	}
	cout << "The ccs: " << endl;
	for (unsigned i = 0; i < out.size(); i++)
	{
		cout << "Component #" << i << endl;
		for (unsigned j = 0; j < out[i].size(); j++)
		{
			cout << " " << out[i][j].u << "," << out[i][j].v << endl;
		}
	}
	return out;
}

//Rename all edge indices in the graph to count from 0 and include all naturals up to the largest index.
vector<Edge> alphaVary(const vector<Edge> &in)
{
	int maxInd = 0;
	for (unsigned i = 0; i < in.size(); i++)
		maxInd = max(maxInd, max(in[i].u, in[i].v));
	vector<bool> used = vector<bool>(maxInd+1, false);
	for (unsigned i = 0; i < in.size(); i++)
	{
		used[in[i].u] = true;
		used[in[i].v] = true;
	}
	int numNodes = 0;
	for (unsigned i = 0; i < used.size(); i++)
	{
		if (used[i])
			numNodes++;
	}
	vector<int> offsets = vector<int>(used.size(), 0);
	int numUsedBelow = 0;
	for (unsigned i = 0; i < used.size(); i++)
	{
		offsets[i] = numUsedBelow;
		if (used[i])
			numUsedBelow++;
	}
	vector<Edge> out = vector<Edge>();
	for (unsigned i = 0; i < in.size(); i++)
		out.push_back(Edge(offsets[in[i].u], offsets[in[i].v]));
	return out;
}

vector<vector<int>> getCostMat(const vector<Edge> &t, int numNodes)
{
	//SplitEdge(e.ind1, e.ind2, 1, e.ind1, e.ind2);
	vector<vector<int>> cMat;
	for (int i = 0; i < numNodes; i++)
	{
		vector<int> cur(numNodes, 0);
		cMat.push_back(cur);
	}
	//

	for (int s = 0; s < numNodes; s++)
	{
		vector<vector<pair<int, int>>> adjs;
		int maxInd = 0;
		for (unsigned i = 0; i < t.size(); i++)
		{
			if (t[i].u > maxInd)
				maxInd = t[i].u;
			if (t[i].v > maxInd)
				maxInd = t[i].v;
		}
		maxInd++;
		for (int i = 0; i < maxInd; i++)
			adjs.push_back(vector<pair<int, int>>());
		for (unsigned i = 0; i < t.size(); i++)
		{
			adjs[t[i].u].push_back(pair<int, int>(t[i].v, 1));
			adjs[t[i].v].push_back(pair<int, int>(t[i].u, 1));
		}
		vector<bool> visited(maxInd, false);
		int cur = s;
		vector<int> distances(maxInd, INT_MAX);
		distances[cur] = 0;
		int numVisited = 0;
		vector<bool> nnFound(maxInd, false);
		for (unsigned i = 0; i < t.size(); i++)
		{
			nnFound[t[i].u] = true;
			nnFound[t[i].v] = true;
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
		for (unsigned i = 0; i < distances.size(); i++)
		{
			cMat[s][i] = distances[i];
			cMat[i][s] = distances[i];
		}
	}
	return cMat;
}

Graph::Graph(string fName, bool silent)
{
	cout << "Reading file " << fName << endl;
	std::size_t isTSV = fName.find(".tsv");
	if (isTSV != std::string::npos)
	{
		cout << "FOUND TSV" << endl;
		//The file is a .tsv who edges must be generated.
		ifstream f;
		f.open(fName);
		char line[256];
		string l = "";
		f.getline(line, 1000);
		vector<Edge> es;
		while (f.rdstate() != std::ifstream::eofbit)
		{
			int e1 = 0;
			int e2 = 0;
			f >> e1 >> e2;
			cout << "(e1,e2): (" << e1 << "," << e2 << ")" << endl;
			es.push_back(Edge(e1-1, e2-1)); //TSV's count from 1
		}

		//Instead of using the graph directly, throw out all nodes not part of largest connected component
		cout << "ES Size: " << es.size() << endl;
		vector<vector<Edge>> components = splitIntoConnectedComponents(es);
		if (components.size() == 1)
			cout << "There is only 1 connected component, using entire graph" << endl;
		//Find the largest, and verify that it is unique:
		unsigned maxCompInd = 0;
		unsigned maxCompSize = 0;
		for (unsigned i = 0; i < components.size(); i++)
		{
			if (components[i].size() > maxCompSize)
			{
				maxCompInd = i;
				maxCompSize = components[i].size();
			}

		}
		cout << "The largest connected component is of size: " << maxCompSize << endl;
		for (unsigned i = 0; i < components.size(); i++)
		{
			if (components[i].size() == maxCompSize && i != maxCompInd)
			{
				cout << "Unsupported: The input graph has more than one largest connected component. Quitting." << endl;
				throw invalid_argument("");
			}
		}
		vector<Edge> lcc = components[maxCompInd];
		es = alphaVary(lcc);
		//

		int numNodes = 0;
		for (unsigned i = 0; i < es.size(); i++)
		{
			if (numNodes < es[i].u)
				numNodes = es[i].u;
			if (numNodes < es[i].v)
				numNodes = es[i].v;
		}
		numNodes++;
		vector<vector<int>> cMat = getCostMat(es, numNodes);
		for (int i = 0; i < numNodes; i++)
		{
			this->costMat.push_back(vector<double>(numNodes, 0));
			for (int j = 0; j < numNodes; j++)
			{
				this->costMat[i][j] = cMat[i][j];
				cout << costMat[i][j] << " ";
			}
			cout << endl;
		}
		this->numNodes = numNodes;
		ofstream fOut;
		string oFName = fName.substr(0, isTSV) + ".tsp";
		//string oFName = fName + ".tsp";
		fOut.open(oFName, std::ofstream::trunc);
		fOut << "NAME : " << oFName << endl;
		fOut << "TYPE : TSP\nDIMENSION : " << numNodes << endl;
		fOut << "EDGE_WEIGHT_TYPE : EXPLICIT\nEDGE_WEIGHT_FORMAT : LOWER_DIAG_ROW" << endl;
		fOut << "DISPLAY_DATA_TYPE : NO_DISPLAY" << endl;
		fOut << "EDGE_WEIGHT_SECTION" << endl;
		for (int i = 0; i < numNodes; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				fOut << costMat[i][j] << " ";
			}
			fOut << endl;
		}
		fOut << "EOF";
		fOut.close();
	}
	else
	{

		ifstream f;
		int numNodes = 0;
		//cout << "Opening " << fName << "...\n";
		f.open(fName);
		//buffer for current line
		char line[256];
		//string version of l for brevity
		string l = "";
		//if the string uses a cost matrix
		bool hasMatrix = false;
		int matType = 0;
		//print starting info and parse to node location.
		do
		{
			f.getline(line, 1000);
			l = (string)line;
			string keywords[5] = { "NAME", "COMMENT", "TYPE" };
			for (int i = 0; i < 3; i++)
			{
				//cout << stringAfterMatch(l, keywords[i]);
			}
			if (!l.compare(0, 9, "DIMENSION"))
			{
				numNodes = atoi((stringAfterMatch(l, "DIMENSION")).c_str());
				//cout << "This problem has " << numNodes << " cities\n";
			}
			if (!l.compare(0, 16, "EDGE_WEIGHT_TYPE"))
			{
				if (stringAfterMatch(l, "EDGE_WEIGHT_TYPE").compare(0, 6, "EUC_2D"))
				{
					cout << "Not EUC_2D!" << endl;
					if (!stringAfterMatch(l, "EDGE_WEIGHT_TYPE").compare(0, 8, "EXPLICIT"))
					{
						hasMatrix = true;
						f.getline(line, 1000);
						l = (string)line;
						cout << "IS: " << stringAfterMatch(l, "EDGE_WEIGHT_FORMAT") << endl;
						if (!stringAfterMatch(l, "EDGE_WEIGHT_FORMAT").compare(0, 14, "UPPER_DIAG_ROW"))
						{
							matType = 1;
						}
						if (!stringAfterMatch(l, "EDGE_WEIGHT_FORMAT").compare(0, 14, "LOWER_DIAG_ROW"))
						{
							matType = 2;
						}
						while (l.compare("EDGE_WEIGHT_SECTION"))
						{
							f.getline(line, 1000);
							l = (string)line;
						}
						break;
					}
					else
						throw invalid_argument("File not in Euclidean 2D format or Cost Matrix Format\n");
				}
				//cout << "File in Euclidean Distance Format\n";
			}
		} while (l.compare("NODE_COORD_SECTION") && l.compare("EDGE_WEIGHT_SECTION"));
		if (numNodes < 3)
			throw invalid_argument("Must have at least 3 nodes");
		this->numNodes = numNodes;
		if (!hasMatrix)
		{
			//Initialize Nodes
			for (int i = 0; i < numNodes; i++)
			{
				Point n;
				nodes.push_back(n);
			}
			//Build the Nodes
			for (int i = 0; i < numNodes; i++)
			{
				string delim = "";
				do
				{
					f.getline(line, 1000, ' ');
				} while (!delim.compare(line));
				//cout << "Node# : " << line << "\n";
				nodes[i].id = atoi(line);
				do
				{
					f.getline(line, 1000, ' ');
				} while (!delim.compare(line));
				//cout << "Dist 1: " << line << "\n";
				nodes[i].x = (GeomPerfectMatching::REAL)atof(line);
				do
				{
					f.getline(line, 1000);
				} while (!delim.compare(line));
				//cout << "Dist 2: " << line << "\n";
				nodes[i].y = (GeomPerfectMatching::REAL)atof(line);
				//cout << "Created Node with id# " << nodes[i].id << " at location " << nodes[i].x << "," << nodes[i].y << "\n";
			}
		}
		else
		{
			cout << "Entered else branch!" << endl;
			cout << "Numnodes: " << numNodes << endl;
			for (int i = 0; i < numNodes; i++)
			{
				vector<double> cr(numNodes, 0.0);
				costMat.push_back(cr);
			}
			if (matType == 2)
			{
				cout << "UPPER/LOWER_DIAG_ROW format detected" << endl;
				for (int i = 0; i < numNodes; i++)
				{
					//cout << "i is: " << i << endl;
					vector<double> found;

					for (int j = 0; j <= i; j++)
					{
						int cur = 0;
						f >> cur;
						found.push_back(cur);
					}
					for (unsigned j = 0; j < found.size(); j++)
					{
						costMat[j][i] = found[j];
						costMat[i][j] = found[j];
					}
				}
			}
			if (matType == 1)
			{
				cout << "UPPER/LOWER_DIAG_ROW format detected" << endl;
				for (int i = 0; i < numNodes; i++)
				{
					//cout << "i is: " << i << endl;
					vector<double> found;
					for (int j = numNodes; j > i; j--)
					{
						int cur = 0;
						f >> cur;
						found.push_back(cur);
					}
					for (unsigned j = 0; j < found.size(); j++)
					{
						costMat[numNodes-1-j][i] = found[found.size()-1-j];
						costMat[i][numNodes-1-j] = found[found.size()-1-j];
					}
				}
			}
		}
		/*for (unsigned i = 0; i < costMat.size(); i++)
		{
			for (unsigned j = 0; j < costMat.size(); j++)
			{
				cout << costMat[i][j] << " ";
			}
			cout << endl;
		}*/
	} 
	//if (silent)
	//	cout.rdbuf(backupBuf);
}


vector<Edge> Graph::getESet()
{
	vector<Edge> es;
	es.reserve((this->numNodes * (this->numNodes - 1)) / 2);
	for (int i = 0; i < this->numNodes; i++)
	{
		for (int j = 0; j < i; j++)
		{
			es.push_back(Edge(j, i));
		}
	}
	return es;
}

//Returns rounded Euclidean distances between nodes[from] and nodes[to].
//Rounded by .tsp file standard.
int Graph::distance(int from, int to)
{
	if (nodes.size() != 0)
	{
		double dist = 0.0f;
		dist = sqrt(pow((nodes[from].x - nodes[to].x), 2) + pow((nodes[from].y - nodes[to].y), 2));
		return (int)(dist + .5);
	}
	return (int)(costMat[from][to]+.5);
}

int Graph::getNumNodes()
{
	return numNodes;
}

//Runs the Minimum Spanning Tree algorithm on the graph.
//It returns an array of Edges, where for all Edges in the array,
//the indices point to two connected nodes in the MST.
vector<Edge> Graph::runMST()
{
	//Step 1: Start by picking a random vertex (or first).
	//Step 2: Create a matrix of size numNodes which, for each node, contains the distances to the start vertex.
	//Step 3: Create a second matrix of size numNodes which is full of the ind of the start vertex.
	//Step 4: Create a third matrix of size numNodes which marks the nodes as used or unused,
	//  marking only the start vertex as used initially.
	//Step 5: Create an empty tree (edge vector).
	//Step 6: Scan the distance matrix, and find the cheapest node to add to the tree which isn't used.
	//Step 7: Add the edge connecting the new node and the owner of the best link to the tree.
	//Step 8: Update the used matrix to true for the new node.
	//Step 9: Update the distance and owner matrices replacing the dist entry if it the cost from the new node to it is
	//  less than the current minimum cost, and replacing the owner entry with the ind of the node just added.
	//Step 10: Check if the MST is complete (size = numNodes-1). If not, repeat steps 6-10.
	
	//cout << "Calculating minimum spanning tree...\n";

	//Step 1:
	int v = 0;

	//Step 2:
	vector<int> dist;
	for (int i = 0; i < numNodes; i++)
		dist.push_back(Graph::distance(0, i));

	//Step 3:
	vector<int> owner;
	for (int i = 0; i < numNodes; i++)
		owner.push_back(0);

	//Step 4:
	vector<int> used;
	used.push_back(1);
	for (int i = 1; i < numNodes; i++)
		used.push_back(0);

	//Step 5:
	vector<Edge> tree;


	//Step 6:
STEP6:
	//get first valid index
	int cheapIndex = 0;
	for (int i = 0; i < numNodes; i++)
	{
		if (!used[i])
		{
			cheapIndex = i;
			break;
		}
	}
	//replace first valid choice with greedy min. cost choice.
	for (int i = 0; i < numNodes; i++)
	{
		if (!used[i] && dist[i] < dist[cheapIndex])
		{
			cheapIndex = i;
		}
	}

	//Step 7
	tree.push_back(Edge(cheapIndex, owner[cheapIndex]));

	//Step 8
	used[cheapIndex] = 1;

	//Step 9
	for (int i = 0; i < numNodes; i++)
	{
		if (Graph::distance(cheapIndex,i) < dist[i])
		{
			dist[i] = Graph::distance(cheapIndex, i);
			owner[i] = cheapIndex;
		}
	}

	//Step 10
	if ((int)tree.size() < numNodes - 1)
		goto STEP6;

	//cout << "Done calculating minimum spanning tree\n";
	return tree;
}

class PQEdgeCompare
{
public:
	PQEdgeCompare() {};
	bool operator() (PQEdge lhs, PQEdge rhs)
	{
		return (lhs.dist) > (rhs.dist);
	}
};

vector<Edge> Graph::runDelaunayMST(vector<Edge> &in)
{
	priority_queue<PQEdge, vector<PQEdge>, PQEdgeCompare> pq;
	//cout << "Generating minimum spanning tree from Delaunay Triangulization...\n";
	int startNode = 0;
	vector<Edge> tree;
	vector<bool> inTree(numNodes, false);
	inTree[0] = true;
	int numAdded = 1;

	//Create array
	vector<int> sizes;
	sizes.assign(numNodes, 0);
	for (unsigned i = 0; i < in.size(); i++)
	{
		sizes[(in[i].u)] += 1;
		sizes[(in[i].v)] += 1;
	}

	vector<vector<int>> adjs;
	for (int i = 0; i < numNodes; i++)
	{
		vector<int> v;
		v.assign(sizes[i], -1);
		adjs.push_back(v);
	}

	//Add edges to 2d vector. -1 means edge used.
	for (unsigned i = 0; i < in.size(); i++)
	{
		//Add edge from ind1 to ind2...
		for (int j = 0; j < sizes[in[i].u]; j++)
		{
			if (adjs[in[i].u][j] == -1)
			{
				adjs[in[i].u][j] = in[i].v;
				break;
			}
		}
		//Add edge from ind2 to ind1...
		for (int j = 0; j < sizes[in[i].v]; j++)
		{
			if (adjs[in[i].v][j] == -1)
			{
				adjs[in[i].v][j] = in[i].u;
				break;
			}
		}
	}
	//Start MST:
	//push all edges adjacent to startNode into pq
	for (unsigned i = 0; i < adjs[startNode].size(); i++)
	{
		pq.push(PQEdge(startNode, adjs[startNode][i], Graph::distance(startNode, adjs[startNode][i])));
	}
	//pop best edge
	PQEdge bestEdge = pq.top();
	pq.pop();

	tree.push_back(Edge(bestEdge.ind1, bestEdge.ind2));
	if (inTree[bestEdge.ind1])
	{
		inTree[bestEdge.ind2] = true;

		for (unsigned i = 0; i < adjs[bestEdge.ind2].size(); i++)
		{
			pq.push(PQEdge(bestEdge.ind2, adjs[bestEdge.ind2][i], Graph::distance(bestEdge.ind2, adjs[bestEdge.ind2][i])));
		}

	}
	else
	{
		inTree[bestEdge.ind1] = true;

		for (unsigned i = 0; i < adjs[bestEdge.ind1].size(); i++)
		{
			pq.push(PQEdge(bestEdge.ind1, adjs[bestEdge.ind1][i], Graph::distance(bestEdge.ind1, adjs[bestEdge.ind1][i])));
		}
	}
	numAdded++;
	//Find rest:
	while (numAdded != numNodes)
	{
		bool foundCur = false;
		PQEdge cur;
		while (!foundCur)
		{
			PQEdge candidate = pq.top();
			pq.pop();
			if ((!inTree[candidate.ind1] && inTree[candidate.ind2]) ||
				(inTree[candidate.ind1] && !inTree[candidate.ind2]))
			{
				cur = candidate;
				foundCur = true;
			}
		}
		
		tree.push_back(Edge(cur.ind1, cur.ind2));
		if (inTree[cur.ind1])
		{
			inTree[cur.ind2] = true;
			for (unsigned i = 0; i < adjs[cur.ind2].size(); i++)
			{
				pq.push(PQEdge(cur.ind2, adjs[cur.ind2][i], Graph::distance(cur.ind2, adjs[cur.ind2][i])));
			}
		}
		else
		{
			inTree[cur.ind1] = true;
			for (unsigned i = 0; i < adjs[cur.ind1].size(); i++)
			{
				pq.push(PQEdge(cur.ind1, adjs[cur.ind1][i], Graph::distance(cur.ind1, adjs[cur.ind1][i])));
			}
		}
		numAdded++;
	}
	//cout << "Done with Delaunay MST\n";
	return tree;

}

vector<int> Graph::degrees(vector<Edge> &edges)
{
	//using numNodes from graph, which already has edges, and takes in a separate edge set...
	//cout << "Calculating degrees..\n";
	vector<int> deg(numNodes, 0);
	for (unsigned i = 0; i < edges.size(); i++)
	{
		deg[edges[i].u]++;
		deg[edges[i].v]++;
	}
	int acc = 0;
	for (int i = 0; i < numNodes; i++)
		acc += deg[i];

	return deg;
}

//returns the node inds which are of odd degree
vector<int> Graph::odds(vector<int> &degrees)
{
	//cout << "Generating nodes with odd degree...\n";
	vector<int> odds;
	//cout << "Initial odds size: " << odds.size() << "\n";
	int numOdds = 0;
	for (int i = 0; i < numNodes; i++)
	{
		if (degrees[i] % 2 == 1)
		{
			//cout << degrees[i] << " is odd.\n";
			odds.push_back(i);
			numOdds++;
		}
	}
	//cout << "Done generating nodes with odd degree.\n";
	//cout << "Found " << numOdds << " vertices with odd degree\n";
	return odds;
}

//Removes duplicate and equivalent edges. That is, (i,j) is removed if (j,i) exists,
// and (i,j) is limited to 1 instance.
vector<Edge> Graph::removeDups(vector<Edge> in)
{
	vector<Edge> out;
	for (unsigned i = 0; i < in.size(); i++)
	{
		bool foundDup = false;
		for (unsigned j = i + 1; j < in.size(); j++)
		{
			if ((in[i].u == in[j].u && in[i].v == in[j].v) ||
				(in[i].u == in[j].v && in[i].v == in[j].u))
			{
				foundDup = true;
			}
		}
		if (!foundDup)
		{
			out.push_back(in[i]);
		}
	}
	return out;
}

//Loads the blossomV geometric solver with the nodes in odds
void Graph::loadGeomSolver(vector<int> &odds, GeomPerfectMatching *gpm)
{
	for (unsigned i = 0; i < odds.size(); i++)
	{
		//If typedef double real:
		GeomPerfectMatching::REAL *coord = new GeomPerfectMatching::REAL[2];
		ZeroMemory(coord, sizeof(GeomPerfectMatching::REAL)* 2);
		coord[0] = nodes[odds[i]].x;
		coord[1] = nodes[odds[i]].y;
		gpm->AddPoint(coord);
		delete[] coord;
		/*int *coord = new int[2];
		ZeroMemory(coord, sizeof(int)* 2);
		coord[0] = (int)(nodes[odds[i]].x+.5);
		coord[1] = (int)(nodes[odds[i]].y+.5);
		gpm->AddPoint(coord);
		delete[] coord;*/
	}
}

//Creates a multigraph from two edge sets
vector<Edge> Graph::createMultigraph(vector<Edge> &e1, vector<Edge> &e2)
{
	//cout << "Creating multigraph...\n";
	vector<Edge> mg;
	for (unsigned i = 0; i < e1.size(); i++)
		mg.push_back(e1[i]);
	for (unsigned i = 0; i < e2.size(); i++)
		mg.push_back(e2[i]);
	//cout << "Done creating multigraph\nSize of multigraph: " << mg.size() << "\n";
	return mg;
}

//print util for nodes
void output(vector<int> in)
{
	//for (unsigned i = 0; i < in.size(); i++)
		//cout << "   " << in[i] << "\n";
}

//print util for edges
void output(vector<Edge> in)
{
	for (unsigned i = 0; i < in.size(); i++)
	{
		//cout << "   " << in[i].ind1 << " to " << in[i].ind2 << "\n";
	}
}

//Helper for createCircuit, see its definition below for documentation
vector<int> generateSubtour(int startNode, vector<vector<int>> &adjs, vector<int> sizes, Graph *g)
{
	//cout << "Entered generatedSubtour2\n";
	vector<int> tour;
	int curNode = startNode;
	tour.push_back(curNode);
	//cout << "Added " << curNode << " to tour.\n";
	int nextNode = -1;
	do
	{
		//cout << "At start of do...while loop in generateSubtour2\n";
		nextNode = -1;
		int minCost = INT_MAX;
		//Find next node...

		for (int i = 0; i < sizes[curNode]; i++)
		{
			if (adjs[curNode][i] != -1 && adjs[curNode][i] != tour.front())
			{
				int curCost = g->distance(curNode, adjs[curNode][i]);
				if (curCost < minCost)
				{
					minCost = curCost;
					nextNode = i;
				}
			}
		}

		if (nextNode == -1)
		{
			for (int i = 0; i < sizes[curNode]; i++)
			{
				if (adjs[curNode][i] != -1)
				{
					nextNode = i;
					/*int curCost = g->distance(curNode, adjs[curNode][i]);
					if (curCost < minCost)
					{
						minCost = curCost;
						nextNode = i;
					}*/
				}
			}
		}
		//
		if (nextNode == -1)
		{
			cout << "ERROR: nextNode still -1." << endl;
			throw logic_error("");
		}
		int n = adjs[curNode][nextNode];
		adjs[curNode][nextNode] = -1;
		nextNode = n;
		for (int j = 0; j < sizes[nextNode]; j++)
		{
			if (adjs[nextNode][j] == curNode)
			{
				adjs[nextNode][j] = -1;
				break;
			}
		}
		


		
		/*for (int i = 0; i < sizes[curNode]; i++)
		{
			if (adjs[curNode][i] != -1)
			{
				nextNode = adjs[curNode][i];
				//cout << "Found that " << curNode << " is adjacent to " << nextNode << "\n";
				//remove link from curNode to nextNode and from nextNode to curNode...
				adjs[curNode][i] = -1;
				for (int j = 0; j < sizes[nextNode]; j++)
				{
					if (adjs[nextNode][j] == curNode)
					{
						adjs[nextNode][j] = -1;
						break;
					}
				}
				break;
			}
		}*/

		if (nextNode == -1)
		{
			throw logic_error("Subtour cannot be generated starting from this node\n");
		}

		//cout << "Adding " << nextNode << " to subtour.\n";
		tour.push_back(nextNode);
		curNode = nextNode;
	} while (tour.front() != tour.back());
	return tour;

}

//Helper for createCircuit, see its definition below for documentation
vector<int> mergePaths(vector<int> &mainPath, vector<int> &secondaryPath)
{
	//Find index into main path
	unsigned mergeStart;
	for (mergeStart = 0; mergeStart < mainPath.size(); mergeStart++)
	{
		if (mainPath[mergeStart] == secondaryPath[0])
			break;
	}

	//Add main path up index, exclusive
	vector<int> newPath;
	for (unsigned i = 0; i < mergeStart; i++)
		newPath.push_back(mainPath[i]);

	//Add all of secondary path after it
	for (unsigned i = 0; i < secondaryPath.size(); i++)
		newPath.push_back(secondaryPath[i]);

	//Finish by adding remainder of main path, excluding index
	for (unsigned i = mergeStart + 1; i < mainPath.size(); i++)
		newPath.push_back(mainPath[i]);

	return newPath;
}

vector<int> Graph::createCircuit(vector<Edge> &mg, Graph *g)
{
	//Step 1: Create a traversable data structure from the vector of edges (or maybe an array of vectors?)
	//Step 2: Pick the start node.
	//Step 3: Generate a subtour from the start node. This is the main path.
	//Step 4: Check if the main path is complete (i.e. length 1 greater than number of nodes in graph). 
	//  If yes, skip to step 9. If no, continue to step 5.
	//Step 5: Find a new start node by iterating over the nodes of the main path until one is found with a nonempty array of adjs.
	//  (this must exist if the path is not complete).
	//Step 6: Generate a subtour from this start node.
	//Step 7: Merge the subtour into the main path.
	//Step 8: Repeat step 4.
	//Step 9: Return the path. 

	//Step 1:
	//cout << "Creating Euler circuit...\n";
	//Create array
	vector<int> sizes;
	sizes.assign(numNodes, 0);
	for (unsigned i = 0; i < mg.size(); i++)
	{
		sizes[(mg[i].u)] += 1;
		sizes[(mg[i].v)] += 1;
	}

	vector<vector<int>> adjs;
	for (int i = 0; i < numNodes; i++)
	{
		vector<int> v;
		v.assign(sizes[i], -1);
		adjs.push_back(v);
	}

	//Add edges to 2d vector. -1 means edge used.
	for (unsigned i = 0; i < mg.size(); i++)
	{
		//Add edge from ind1 to ind2...
		for (int j = 0; j < sizes[mg[i].u]; j++)
		{
			if (adjs[mg[i].u][j] == -1)
			{
				adjs[mg[i].u][j] = mg[i].v;
				break;
			}
		}
		//Add edge from ind2 to ind1...
		for (int j = 0; j < sizes[mg[i].v]; j++)
		{
			if (adjs[mg[i].v][j] == -1)
			{
				adjs[mg[i].v][j] = mg[i].u;
				break;
			}
		}
	}

	//Step 2:
	int startNode = 0;

	//Step 3:
	//cout << "Multigraph: \n";
	//output(mg);
	//cout << "About to generate first tour\n";
	vector<int> mainPath = generateSubtour(startNode, adjs, sizes, g);
	//cout << "First mainPath: \n";
	//output(mainPath);
	/*cout << "Adjacent entries for each entry in mainPath:\n";
	for (unsigned i = 0; i < mainPath.size(); i++)
	{
		cout << "Entry for " << mainPath[i] << "\n";
		output(adjs[mainPath[i]]);
	}*/
	//Step 5:
	while (mainPath.size() != mg.size() + 1)
	{
		//cout << "Starting subtour search. Current path size: " << mainPath.size() << "\n";
		//cout << "Outputting adjs: \n";
		/*for (unsigned i = 0; i < numNodes; i++)
		{
			cout << "Entry for " << i << "\n";
			output(adjs[i], sizes[i]);
		}*/
		for (unsigned i = 0; i < mainPath.size(); i++)
		{
			int j = 0;
			for (j = 0; j < sizes[mainPath[i]]; j++)
			{
				if (adjs[(mainPath[i])][j] != -1) //then we have an unused edge from node mainpath[i] to node adjs[i][j].
				{
					startNode = mainPath[i];
					//cout << "Found " << startNode << " which is adjacent to " << adjs[(mainPath[i])][j] << "\n";
					break;
				}
			}
			if (j != sizes[mainPath[i]])
				break;
		}
		//cout << "About to call generateSubtour on " << startNode << "\n";
		vector<int> subtour = generateSubtour(startNode, adjs, sizes, g);
		//cout << "Completed generatedSubtour on " << startNode << "\n";
		vector<int> merged = mergePaths(mainPath, subtour);
		mainPath.assign(merged.begin(), merged.end());
	}
	//cout << "Done creating Euler Circuit.\n";
	//cout << "Size of mg: " << mg.size() << "\n";
	//cout << "Size of mainPath: " << mainPath.size() << "\n";
	return mainPath;

}

//Shortcuts circuit arbitrarily. For testing.
vector<int> Graph::shortcutCircuit(vector<int> &circ)
{
	//cout << "Shortcutting circuit...\n";
	vector<int> result;
	vector<int> inResult(numNodes, 0);
	
	for (unsigned i = 0; i < circ.size(); i++)
	{
		if (!inResult[circ[i]])
		{
			inResult[circ[i]] = 1;
			result.push_back(circ[i]);
		}
	}
	//make sure to append the end, because it will have incorrectly failed 
	//the isMem check above.
	result.push_back(circ.back());
	//cout << "Done shortcutting\n";
	return result;
}

//Shortcuts circuit by walking down circuit and replacing previous uses of a particular 
//vertex if the next occurance costs less.
vector<int> Graph::shortcutCircuitSelective(vector<int> &circ)
{
	vector<int> firstResult;
	vector<int> locInResult(numNodes, -1);

	for (unsigned i = 0; i < circ.size(); i++)
	{
		if (locInResult[circ[i]] == -1)
		{
			firstResult.push_back(circ[i]);
			locInResult[circ[i]] = firstResult.size() - 1;

		}
		else
		{
			int cur = locInResult[circ[i]];

			if (cur == 0)
				continue;
			if (i == circ.size() - 1)
				continue;
			if (cur == locInResult.size() - 1)
				continue;
			int backLoc = cur - 1;
			if (firstResult[backLoc] == -1)
			{
				for (int j = backLoc-1; j >= 0; j--)
				{
					if (firstResult[j] != -1)
					{
						backLoc = j;
						break;
					}
				}
			}
			int curDistBehind = Graph::distance(firstResult[cur], firstResult[backLoc]);
			int newDistBehind = Graph::distance(circ[i], firstResult.back());

			int newDistAhead = Graph::distance(circ[i], circ[i + 1]);
			//find next non -1 value
			int next = -1;
			for (unsigned j = cur + 1; j < firstResult.size(); j++)
			{
				if (firstResult[j] != -1)
				{
					next = j;
					break;
				}
			}
			if (next == -1)
				continue;
			int curDistAhead = Graph::distance(firstResult[cur], firstResult[next]);

			if (newDistBehind + newDistAhead < curDistBehind + curDistAhead)
			{
				firstResult[cur] = -1;
				firstResult.push_back(circ[i]);
				locInResult[circ[i]] = firstResult.size() - 1;
			}
		}
	}
	//clean result of -1's
	vector<int> finalResult;
	for (unsigned i = 0; i < firstResult.size(); i++)
	{
		if (firstResult[i] != -1)
			finalResult.push_back(firstResult[i]);
	}
	finalResult.push_back(firstResult[0]);
	return finalResult;

}

void Graph::writeCircuit(vector<int> toWrite, string fName)
{
	//cout << "Writing circuit to file...\n";
	ofstream f;
	f.open(fName);
	for (unsigned i = 0; i < toWrite.size(); i++)
	{
		f << toWrite[i] + 1 << "\n";
	}
	f.close();
	//cout << "Done writing circuit to file\n";
}

int Graph::tourCost(vector<int> &tour)
{
	  /*cout << "Calculating the tour cost:" << endl;*/
	int acc = 0;
	for (unsigned i = 0; i < tour.size() - 1; i++)
	{
			/*cout << "The distance from node " << tour[i] << " to node " << tour[i + 1] << " is " << Graph::distance(tour[i], tour[i + 1]) << endl;*/
		acc += Graph::distance(tour[i], tour[i + 1]);
	}
	return acc;
}

int Graph::edgeSetCost(const vector<Edge> &t)
{
	int acc = 0;
	for (unsigned i = 0; i < t.size(); i++)
		acc += this->distance(t[i].u, t[i].v);
	return acc;
}

void Graph::drawTree(const vector<Edge> &t)
{
	ofstream f;
	f.open("temp.input", std::ofstream::out | std::ofstream::trunc);
	for (unsigned i = 0; i < t.size(); i++)
	{
		Edge e = t[i];
		int x1 = (int)(this->nodes[e.u].x + .5);
		int y1 = (int)(this->nodes[e.u].y + .5);
		int x2 = (int)(this->nodes[e.v].x + .5);
		int y2 = (int)(this->nodes[e.v].y + .5);
		f << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
	}
	f.close();
	system("visualizer.exe temp.input");
}