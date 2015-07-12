//Main implementation file for the the Lambda-Random Tree Sampler
//Kyle Genova
#include "lrs.h"
#include "util.h"
#include "matrix.h"
#include "heldkarp.h"
#include "christofides.h"
#include "gamma.h"

//Description of the process:
// Step 1: Find the held-karp solution vector, x, to the TSP problem

// Step 2: Rescale this vector by (n-1)/n to get the vector z.

// Step 3: Iteratively approximate gamma values which satisfy lambda = e^gamma, for each edge

// Step 4: Calculate the lambda values for the graph

// Step 5: Iteratively generate a lambda-random tree by calculating the appropriate probability to
//         add or disregard each edge in the graph, contracting the graph for added edges and deleting 
//         edges which were not added.

//Description of Step 1:

//  Run concorde with the appropriate variables on the .tsp file to be solved.

//  Concorde writes the held-karp solution to a file. Read in this file, and create an x vector from it.


//Description of Step 2:

//  Take the x vector, and create a vector z such that z_i = ((n - 1)/n)*x_i, where n is the number of nodes in the graph.


//Description of Step 3:

//  The goal is to create a gamma vector iteratively. Start by setting gamma to the zero vector.

//  Several quantities will be needed:

//    q_e(gamma) = Sum (over T /owns e) exp(gamma(T))
//                 ---------------------------------
//                 Sum (over Ts?)       exp(gamma(T))

//    gamma(T) = Sum(over f \in T)  gamma_f

//    q_e: The probability that edge e will be included in a spanning tree T that is chosen with probability proportional
//    to exp(gamma(T))

//    epsilon: a chosen error bound

//    delta:     (  q_e(gamma)(1-(1+epsilon/2)z_e)     )
//           log (  -----------------------------      )
//               (  (1 - q_e(gamma))(1 + epsilon/2)z_e )

//    gamma' := gamma'_e = gamma_e - delta,  gamma'_f = gamma_f \forall f \in E \ {e}

//   While there exists and edge e s.t. g_e(gamma) > (1 + epsilon)z_e, compute gamma', set gamma := gamma'.


//Description of Step 4:

//  lambda_e := e^(gamma_e) for each lambda.

//Description of Step 5: 

//  At each step, the question is to find some probability p_j, which is the probability that e_j is
//  in a lambda-random tree of a graph where: 
//    
//    1) all edges which have already been added to the tree have been contracted from the original graph
//    2) all edges which were not added to the tree have been deleted from the graph.

//  Once this graph is created, create the weighted Laplacian matrix for the graph:

//    L_i,j = { -lambda_e :  e = (i,j) \in E
//			  { Sum(where e \in delta ({i})) lambda_e : i = j
//			  { 0 : o/w

//  Take the absolute value of any cofactor of this matrix, this represents:
//    Sum (over T \in Ts) (Product e \in T) lambda_e

//  NOT CERTAIN ABOUT THIS STEP
//  Then p_j is a fraction where: 
//     The numerator is the sum above for G' with edge f contracted
//     The denominator is the sum above for G'

// Alternatively, p_j is lambda_f * the effective resistance of f in G' where lambda values are edge conductances.

// Add e_j to the graph with probability p_j, updating G' accordingly.

// Do this until all edges have been considered.


//At each step, we take the current graph, and get the denominator
// as the abs(cof(wLap(graph)))
//Then, we create a new graph, cGraph, and a new coef, nCoef,
//cGraph is the graph with edge e contracted.
//nCoef is prodCoef * the weight of e.
// calculate the num as nCoef * abs(cof(wLap(cGraph)))
//Calculate the probability as num/denom
//If the edge should be added, then
//  graph <- cGraph
//  prodCoef <- nCoef
//  tree <- tree + e
//If the edge shouldn't be added, then
//  graph <- graph - e
//  prodCoef <- prodCoef
//  tree <- tree
// Repeat while there exists an edge not yet considered


vector<LPEdge> sampleLRTNew(vector<ModEdge> lGraph, bool silent)
{
	vector<LPEdge> tree;
	while (lGraph.size() != 0)
	{
		ModEdge e = lGraph.back();
		double denom = abs(cofactorSparse(wLaplacianNew(lGraph)));
		vector<ModEdge> cGraph = contractEdgeNew(lGraph, e);
		double num = e.weight * (abs(cofactorSparse(wLaplacianNew(cGraph))));
		double prob = num / denom;
		if (should(prob) || lGraph.size() == 1)
		{
			lGraph = cGraph;
			LPEdge outEdge(e.orig0, e.orig1, e.weight);
			tree.push_back(outEdge);
		}
		else
		{
			lGraph.pop_back();
		}
	}
	return tree;
}

vector<LPEdge> sampleLRTBig(vector<ModEdge> lGraph, bool silent)
{
	vector<LPEdge> tree;
	while (lGraph.size() != 0)
	{
		ModEdge e = lGraph.back();
		double denom = cofactorLog(wLaplacianNew(lGraph));
		vector<ModEdge> cGraph = contractEdgeNew(lGraph, e);
		double num = cofactorLog(wLaplacianNew(cGraph));
		double prob = e.weight * exp(num - denom);
		if (should(prob) || lGraph.size() == 1)
		{
			lGraph = cGraph;
			LPEdge outEdge(e.orig0, e.orig1, e.weight);
			tree.push_back(outEdge);
		}
		else
		{
			lGraph.pop_back();
		}
	}
	return tree;
}

vector<LPEdge> sampleLRTInverse(vector<ModEdge> lGraph, bool silent)
{
	vector<LPEdge> tree;
	while (lGraph.size() != 0)
	{
		ModEdge e = lGraph.back();
		double p = prob(lGraph, e);
		if (should(p))
		{
			lGraph = contractEdge(lGraph, e);
			tree.push_back(LPEdge(e.orig0, e.orig1, e.weight));
		}
		else
			lGraph.pop_back();
	}
	return tree;
}

int getLRSTourFast(vector<ModEdge> lGraph, vector<double> lambda, Graph g, int numNodes)
{
	vector<Edge> tree = lpToEdge(sampleLRTNew(lGraph, true));
	vector<int> tour;
	return (numNodes == tree.size() + 1) ? christofides(&g, tree, tour, true) : INT_MAX;
}

vector<double> getLambdaBaked(string fName, const vector<ModEdge> &g)
{
	vector<ModEdge> graph = g;
	{
		string lambdaF = fName.substr(0, fName.length() - 4) + ".lambda";
		ifstream f;
		f.open(lambdaF);
		vector<LPEdge> lp;
		string line;
		while (!f.eof())
		{
			int from = -1;
			int to = -1;
			double unused = -1.;
			double l = -1.;
			f >> from >> to >> unused >> l;
			cout << "From: " << from << ", To: " << to << ", Lambda: " << l << endl;
			for (unsigned i = 0; i < graph.size(); i++)
			{
				if ((graph[i].end0 == from && graph[i].end1 == to) ||
					(graph[i].end0 == to && graph[i].end1 == from))
				{
					cout << "Matched! Changing weight from " << graph[i].weight << " to " << l << endl;
					graph[i].weight = l;
				}
			}
		}
		f.close();
	}
	vector<double> lambda;
	for (unsigned i = 0; i < graph.size(); i++)
		lambda.push_back(graph[i].weight);
	return lambda;
}

int runLRSBaked(string fName, int numSamples)
{
	srand((unsigned(time(NULL))));
	cout << "numSamples: " << numSamples << endl;
	cout << "Getting gammas..." << endl;
	int best = INT_MAX;

	vector<LPEdge> hkSoln = getHKSolutionsBaked(fName, false); 
	vector<ModEdge> graph = xToZ(lpToMod(hkSoln));
	vector<double> lambda = getLambdaBaked(fName, graph);
	weightGraph(graph, lambda);
	int numNodes = getNumNodes(graph);
	Graph g = Graph(fName, false);
	cout << "Done getting gammas" << endl;

	int numErrs = 0;
	for (int i = 0; i < numSamples; i++)
	{
		int n = numNodes;
		int cur = getLRSTourFast(graph, lambda, g, n);
		best = min(cur, best);
		if (cur == INT_MAX)
			numErrs++;
	}
	cout << "Number of errors: " << numErrs << " out of " << numSamples << " samples" << endl;
	return best;
}

int runLRS(string fName, int numSamples, bool silent)
{
	srand((unsigned(time(NULL))));
	cout << "numSamples: " << numSamples << endl;
	cout << "Getting gammas..." << endl;
	int best = INT_MAX;
	vector<double> gammaIn;
	vector<LPEdge> hkSoln = getHKSolutions(fName, silent);
	vector<ModEdge> graph = xToZ(lpToMod(hkSoln));
	vector<double> lambda = getLambda(graph);
	weightGraph(graph, lambda);
	int numNodes = getNumNodes(graph);
	Graph g = Graph(fName, silent);
	cout << "Done getting gammas" << endl;

	int numErrs = 0;
#pragma omp parallel num_threads(4) shared(gammaIn, graph, lambda, g, numNodes)
	{
		#pragma omp for 
		for (int i = 0; i < numSamples; i++)
		{
			int n = numNodes;
			int cur = getLRSTourFast(graph, lambda, g, n);
#pragma omp critical
			{

				if (cur < best)
					best = cur;
				if (cur == INT_MAX)
					numErrs++;
#pragma omp flush(best)
			}
		}
	}
	cout << "Number of errors: " << numErrs << " out of " << numSamples << " samples" << endl;
	return best;
}