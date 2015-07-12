// christofides.cpp
// Function Library implementation file for Christofides algorithm
// Kyle Genova

#include "christofides.h"



void verifyTour(const vector<int> &tour, int numNodes)
{
	//cout << "Beginning tour verification:" << endl;
	vector<bool> inTour(numNodes, false);
	for (unsigned i = 0; i < tour.size(); i++)
	{
		if (inTour[tour[i]] && i != tour.size() - 1)
		{
			cout << "ERROR: Repeated node: " << tour[i] << endl;
		}
		inTour[tour[i]] = true;
	}
	for (unsigned i = 0; i < inTour.size(); i++)
	{
		if (!inTour[i])
			cout << "ERROR: Missing node: " << i << endl;
	}
	if (tour.size() != numNodes + 1)
		cout << "ERROR: Tour is incorrrect length." << endl;
	//cout << "Tour verification finished" << endl;
}
double matchingCost(Graph *g, vector<Edge> &spanTree, int &numEIM)
{
	if (g->costMat.size() == 0)
	{
		vector<int> odds = g->odds(g->degrees(spanTree));
		GeomPerfectMatching *gpm = new GeomPerfectMatching((int)odds.size(), 2);
		gpm->options.verbose = false;
		g->loadGeomSolver(odds, gpm);
		GeomPerfectMatching::GPMOptions gpmOptions;
		gpm->SolveComplete();
		vector<Edge> geoSolverOut;
		for (unsigned i = 0; i < odds.size(); i++)
			geoSolverOut.push_back(Edge(odds[i], odds[gpm->GetMatch(i)]));
		delete gpm;
		vector<Edge> compressed = g->removeDups(geoSolverOut);
		double cost = 0;
		for (unsigned i = 0; i < compressed.size(); i++)
			cost += g->distance(compressed[i].u, compressed[i].v);
		numEIM = compressed.size();
		return cost;
	}
	vector<int> odds = g->odds(g->degrees(spanTree));
	PerfectMatching *pm = new PerfectMatching(odds.size(), odds.size()*odds.size() + 1 / 2);
	pm->options.verbose = false;
	vector<Edge> eSet = g->getESet();
	vector<Edge> trackA;
	for (unsigned i = 0; i < odds.size(); i++)
	{
		for (unsigned j = 0; j < i; j++)
		{
			pm->AddEdge(i, j, g->distance(odds[i], odds[j]));
			trackA.push_back(Edge(i, j));
		}
	}
	pm->Solve();
	vector<Edge> mpm;
	for (unsigned i = 0; i < trackA.size(); i++)
	{
		if (pm->GetSolution(i))
			mpm.push_back(Edge(odds[trackA[i].u], odds[trackA[i].v]));
	}
	vector<Edge> compressed = g->removeDups(mpm);
	numEIM = compressed.size();
	double cost = 0;
	for (unsigned i = 0; i < compressed.size(); i++)
		cost += g->distance(compressed[i].u, compressed[i].v);
	return cost;
};

// [minSpanTree g] takes in a Graph g
// g: The graph on which to run the spanning tree algorithm
// Returns: The a minimum spanning tree using a delaunay triangulation.
vector<Edge> minSpanTree(Graph *g)
{
	if (g->costMat.size() == 0)
	{
		GeomPerfectMatching *delaunayWrapper = new GeomPerfectMatching(g->getNumNodes(), 2);
		delaunayWrapper->options.verbose = false;
		for (int i = 0; i < g->getNumNodes(); i++)
		{
			//IF typedef double REAL:
			GeomPerfectMatching::REAL *coord = new GeomPerfectMatching::REAL[2];
			coord[0] = (GeomPerfectMatching::REAL)g->nodes[i].x;
			coord[1] = (GeomPerfectMatching::REAL)g->nodes[i].y;
			delaunayWrapper->AddPoint(coord);
			/*int *coord = new int[2];
			coord[0] = (int)(g->nodes[i].x+.5);
			coord[1] = (int)(g->nodes[i].y+.5);
			delaunayWrapper->AddPoint(coord);*/
		}
		delaunayWrapper->InitDelaunay();
		Block<GeomPerfectMatching::Edge> * edges = delaunayWrapper->edges;
		vector<Edge> mstSparse;
		GeomPerfectMatching::Edge *firstMSTEdge = edges->ScanFirst();
		mstSparse.push_back(Edge(firstMSTEdge->head[0], firstMSTEdge->head[1]));
		while (true)
		{
			GeomPerfectMatching::Edge *cur = edges->ScanNext();
			if (cur == NULL)
				break;
			mstSparse.push_back(Edge(cur->head[0], cur->head[1]));
		}
		return g->runDelaunayMST(mstSparse);
	}
	return g->runDelaunayMST(g->getESet());
}

// [christofides g tour] takes in a Graph g and a reference to an empty vector<int> where the answer
//   will be written.
// g: The graph on which to run the christofides algorithm
// tour: an empty vector<int> where the result of the algorithm will be stored.
// Returns: The size of the tour found. 0 if an error occured.
int christofides(Graph *g, vector<int> &tour, bool silent)
{
	// Do initial delaunay triangulation
	GeomPerfectMatching *delaunayWrapper = new GeomPerfectMatching(g->getNumNodes(), 2);
	delaunayWrapper->options.verbose = !silent;

	for (int i = 0; i < g->getNumNodes(); i++)
	{
		GeomPerfectMatching::REAL *coord = new GeomPerfectMatching::REAL[2];
		coord[0] = g->nodes[i].x;
		coord[1] = g->nodes[i].y;
		delaunayWrapper->AddPoint(coord);
		delete[] coord;
	}
	delaunayWrapper->InitDelaunay();
	Block<GeomPerfectMatching::Edge>* edges = delaunayWrapper->edges;
	vector<Edge> mstSparse;
	GeomPerfectMatching::Edge *firstMSTEdge = edges->ScanFirst();
	mstSparse.push_back(Edge(firstMSTEdge->head[0], firstMSTEdge->head[1]));
	while (true)
	{
		GeomPerfectMatching::Edge *cur = edges->ScanNext();
		if (cur == NULL)
			break;
		mstSparse.push_back(Edge(cur->head[0], cur->head[1]));
	}
	// Done with initial delaunay triangulation

	vector<Edge> delMSTEdges = g->runDelaunayMST(mstSparse);
	vector<int> deg = g->degrees(delMSTEdges);
	vector<int> odds = g->odds(deg);

	// Geometric interface to blossomV; only necessary for > 10,000 nodes or so
	GeomPerfectMatching *gpm = new GeomPerfectMatching((int)odds.size(), 2);
	gpm->options.verbose = !silent;
	g->loadGeomSolver(odds, gpm);
	GeomPerfectMatching::GPMOptions gpmOptions;

	//cout << "Solving Euclidean minimum perfect matching..." << endl;

	if (g->getNumNodes() >= 2500)
	{
		//cout << "Running standard Geometric solver..." << endl;
		gpm->Solve();
	}
	else
	{
		//cout << "Running complete Geometric solver..." << endl;
		gpm->SolveComplete();
	}
	//cout << "Done solving Euclidean minimum perfect matching." << endl;
	//cout << "Loading results from BlossomV..." << endl;
	vector<Edge> geoSolverOut;
	for (unsigned i = 0; i < odds.size(); i++)
	{
		Edge e(odds[i], odds[gpm->GetMatch(i)]);
		geoSolverOut.push_back(e);
	}
	delete gpm;
	// Done unloading the solver.

	// We have to remove the duplicate edges from the vector
	vector<Edge> compressed = g->removeDups(geoSolverOut);

	// Create the multigraph of the MST edges and the MPM edges
	vector<Edge> multi = g->createMultigraph(compressed, delMSTEdges);

	// Create the Eulerian Tour
	vector<int> circ = g->createCircuit(multi, g);

	// Shortcut the duplicate vertices to get a Hamiltonian Tour
	tour = g->shortcutCircuitSelective(circ);
	int tourSize = g->tourCost(tour);

	//cout << "Done with Christofides. Final tour length: " << tourSize << endl;
	int numEIM = 0;
	double mCost = matchingCost(g, mstSparse, numEIM);
	//cout << "Matching cost for standard Christofides: " << mCost << endl;
	int tCost = g->edgeSetCost(mstSparse);
	//cout << "Tree cost for standard Christofides: " << tCost << endl;)
	//cout << "Total sum cost for standard Christofides: " << tCost + mCost << endl;)
	return tourSize;
}

// [christofidesMat g tour] takes in a Graph g and a reference to an empty vector<int> where the answer
//   will be written.
// g: The graph on which to run the christofides algorithm
// tour: an empty vector<int> where the result of the algorithm will be stored.
// Returns: The size of the tour found. 0 if an error occured.
int christofidesMat(Graph *g, vector<int> &tour, bool silent)
{
	try
	{
		vector<Edge> eSet = g->getESet();
		vector<Edge> delMSTEdges = g->runDelaunayMST(eSet);
		vector<int> deg = g->degrees(delMSTEdges);
		vector<int> odds = g->odds(deg);
		PerfectMatching *pm = new PerfectMatching(odds.size(), odds.size()*odds.size() + 1 / 2);
		pm->options.verbose = false;
		eSet = g->getESet();
		vector<Edge> trackA;
		for (unsigned i = 0; i < odds.size(); i++)
		{
			for (unsigned j = 0; j < i; j++)
			{
				pm->AddEdge(i, j, g->distance(odds[i], odds[j]));
				trackA.push_back(Edge(i, j));
			}
		}
		pm->Solve();
		vector<Edge> mpm;
		for (unsigned i = 0; i < trackA.size(); i++)
		{
			if (pm->GetSolution(i))
				mpm.push_back(Edge(odds[trackA[i].u], odds[trackA[i].v]));
		}

		vector<Edge> compressed = g->removeDups(mpm);
		vector<Edge> multi = g->createMultigraph(compressed, delMSTEdges);
		vector<int> circ = g->createCircuit(multi, g);
		vector<int> hamiltonianCirc = g->shortcutCircuitSelective(circ);
		tour = hamiltonianCirc;
		int tourSize = g->tourCost(hamiltonianCirc);
		return tourSize;
	}
	catch (invalid_argument e)
	{
		//cout << "Error: " << e.what() << "\n" << "Quitting.\n";
		return 0;
	}
}

void printDegrees(const vector<int> &vs, const vector<int> &degs)
{
	if (vs.size() != degs.size())
	{
		cout << "Major error in degree calculation" << endl;
		throw logic_error("");
	}
	for (unsigned i = 0; i < vs.size(); i++)
	{
		cout << "The degree of node " << vs[i] << " is " << degs[i] << endl;
	}
}

// [christofides g spanTree tour] takes in a Graph g, a vector<Edge> spanTree, and a reference to an empty 
//   vector<int> where the answer will be written.
// g: The graph on which to run the christofides algorithm
// spanTree: The tree to use instead of finding the default minimum spanning tree.
// tour: an empty vector<int> where the result of the algorithm will be stored.
// silent: if silent, no print statements will output
// Returns: the size of the tour found. 0 If an error occured.
int christofides(Graph *g, vector<Edge> &spanTree, vector<int> &tour, bool silent)
{
	try
	{
		vector<Edge> delMSTEdges = spanTree;
		/*cout << "Entered christofides. This is the spanning tree: " << endl;
		print(delMSTEdges);
		cout << "Calculating degrees:" << endl;
		cout << "There are " << g->numNodes << " nodes" << endl;*/
		vector<int> deg = g->degrees(delMSTEdges);
		vector<int> odds = g->odds(deg);
		/*cout << "These are the nodes with odd degree:" << endl;
		print(odds);*/
		//geom interface to blossomV, new for >10,000 nodes.
		GeomPerfectMatching *gpm = new GeomPerfectMatching((int)odds.size(), 2);
		gpm->options.verbose = !silent;
		g->loadGeomSolver(odds, gpm);
		GeomPerfectMatching::GPMOptions gpmOptions;
		//cout << "Solving Euclidean minimum perfect matching...\n";
		if (g->getNumNodes() >= 2500)
			gpm->Solve();
		else
			gpm->SolveComplete();
		//cout << "Done solving Euclidean minimum perfect matching\n";
		//end interaction with blossomV

		//load solutions
		//cout << "Unloading solver...\n";
		vector<Edge> geoSolverOut;
		for (unsigned i = 0; i < odds.size(); i++)
		{
			Edge e(odds[i], odds[gpm->GetMatch(i)]);
			geoSolverOut.push_back(e);
		}
		delete gpm;
		/*cout << "This is the minimum perfect matching:" << endl;
		  print(geoSolverOut);*/
		//cout << "Done unloading solver\n";
		//cout << "Removing duplicates...\n";
		vector<Edge> compressed = g->removeDups(geoSolverOut);
		/*cout << "This is the compressed MPM:" << endl;
		print(compressed);*/
		//cout << "Done removing duplicates\n";
		vector<Edge> multi = g->createMultigraph(compressed, delMSTEdges);
		/*cout << "This is the multigraph:" << endl;
		print(multi);*/
		vector<int> circ = g->createCircuit(multi, g);
		vector<int> hamiltonianCirc = g->shortcutCircuitSelective(circ);
		tour = hamiltonianCirc;
		int cost = g->tourCost(hamiltonianCirc);
		//cout << "Final tour length: " << cost << "\n";

		/*cout << "This is the final tour: " << endl;
		for (unsigned i = 0; i < tour.size(); i++)
		cout << "  " << tour[i] << endl;
		cout << "The tour has " << tour.size() << "stops, counting both ends" << endl;*/
		verifyTour(tour, g->numNodes);
		return cost;
	}
	catch (invalid_argument e)
	{
		//cout << "Error: " << e.what() << "\n" << "Quitting.\n";
		return 0;
	}
}

// [christofidesMat g spanTree tour] takes in a Graph g, a vector<Edge> spanTree, and a reference to an empty 
//   vector<int> where the answer will be written.
// g: The graph on which to run the christofides algorithm
// spanTree: The tree to use instead of finding the default minimum spanning tree.
// tour: an empty vector<int> where the result of the algorithm will be stored.
// silent: if silent, no print statements will output
// Returns: the size of the tour found. 0 If an error occured.
int christofidesMat(Graph *g, vector<Edge> &spanTree, vector<int> &tour, bool silent)
{
	try
	{
		//cout << "entered try, spanTree has " << spanTree.size() << " edges" << endl;
		vector<Edge> delMSTEdges = spanTree;
		//cout << "defined delMSTEdges." << endl;
		vector<int> deg = g->degrees(delMSTEdges);
		//cout << "got deg" << endl;
		vector<int> odds = g->odds(deg);
		//cout << "got odds, there are " << odds.size() << endl;
		//PerfectMatching *pm = new PerfectMatching(odds.size(), odds.size()*odds.size() + 1 / 2);
		PerfectMatching pm(odds.size(), (odds.size()*(odds.size() - 1)) / 2);
		//cout << "declared pm" << endl;
		pm.options.verbose = !silent;
		//cout << "set verbose" << endl;
		//vector<Edge> eSet = g->getESet();
		//cout << "got eSet, has " << eSet.size() << " edges" << endl;
		vector<Edge> trackA;
		trackA.reserve((g->numNodes * (g->numNodes - 1)) / 2);
		//cout << "about to enter for." << endl;
		for (unsigned i = 0; i < odds.size(); i++)
		{
			for (unsigned j = 0; j < i; j++)
			{
				pm.AddEdge(i, j, g->distance(odds[i], odds[j]));
				//cout << g->distance(odds[i], odds[j]) << endl;
				trackA.push_back(Edge(i, j));
			}
		}
		//cout << "about to solve" << endl;
		pm.Solve();
		//cout << "done solving" << endl;
		vector<Edge> mpm;
		for (unsigned i = 0; i < trackA.size(); i++)
		{
			if (pm.GetSolution(i))
				mpm.push_back(Edge(odds[trackA[i].u], odds[trackA[i].v]));
		}

		//delete pm;
		vector<Edge> compressed = g->removeDups(mpm);
		vector<Edge> multi = g->createMultigraph(compressed, delMSTEdges);
		vector<int> circ = g->createCircuit(multi, g);
		vector<int> hamiltonianCirc = g->shortcutCircuitSelective(circ);
		tour = hamiltonianCirc;
		int tourSize = g->tourCost(hamiltonianCirc);
		return tourSize;
	}
	catch (invalid_argument e)
	{
		cout << "Error: " << e.what() << "\n" << "Quitting.\n";
		return 0;
	}
}

//runs christofides on the file fName and outputs the results to outFName
int writeChristofides(string fName, string outFName, bool silent)
{
	Graph *g = new Graph(fName, silent);
	vector<int> tour;
	int r = christofides(g, tour, silent);
	g->writeCircuit(tour, outFName);
	return r;
}

//wrapper for christofides
int christofides(string fName, vector<int> &tour, bool silent)
{
	Graph *g = new Graph(fName, silent);
	return christofides(g, tour, silent);
}




