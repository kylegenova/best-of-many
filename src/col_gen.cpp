// colGen.cpp : Defines the entry point for the Column Generation standalone executeable, and contains
//   implementations of the static library functions.
//

#include "col_gen.h"
#include "windows.h"
#include "gurobi_c++.h"

//Uncomment before running a Visual Studio performance analysis
//#define PERFORMANCE_ANALYZER


#ifdef PERFORMANCE_ANALYZER
#define overrideFName "perf_analyzer.tsp"
#else
#define overrideFName ""
#endif

//Initial O(n^2) implementation, used also for testing functionality
vector<LPEdge> runMaxSpanTree(vector<LPEdge> &in, int numNodes)
{
	vector<LPEdge> tree;
	vector<bool> nodesInTree(numNodes, false);

	//Arbitrarily start with first node
	if (numNodes == 0)
	{
		vector<LPEdge> t;
		return t;
	}
	nodesInTree[0] = true;

	while (tree.size() != numNodes - 1)
	{
		int bestId = -1;
		//find most expensive valid edge to take, i.e. edge with exactly one edge in nodesInTree
		for (unsigned i = 0; i < in.size(); i++)
		{
			//check if this edge is valid
			bool fst = nodesInTree[in[i].end0];
			bool snd = nodesInTree[in[i].end1];
			bool isValid = ((fst && !snd) || (!fst && snd));
			//
			if (isValid)
			{
				if ((bestId == -1) || (in[bestId].weight < in[i].weight))
				{
					bestId = i;
				}
			}
		}
		//
		int nodeAdded = nodesInTree[in[bestId].end0] ? in[bestId].end1 : in[bestId].end0;
		nodesInTree[nodeAdded] = true;
		tree.push_back(in[bestId]);
	}
	return tree;
}

//Sum the weights of edges in a tree or graph
double edgeSum(vector<LPEdge> &in)
{
	double acc = 0.0;
	for (unsigned i = 0; i < in.size(); i++)
		acc += in[i].weight;
	return acc;
}

//Print the result of the model, if available
void output(GRBModel &model)
{
	if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
	{
		std::cout << "Found a solution" << endl;
		std::cout << "Objective: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
	}
}

//Print a histogram on the command line. Only used by the local main function
void outputHistogram(vector<int> &histogram, int bestBucket, int numEles)
{
	//output histogram of results:
	cout << "Here is where the useful trees came from: \n";
	cout << "Higher up means earlier in the search process.\n";
	cout << "The bar in red indicates the location of the best tour\n";
	HANDLE console = GetStdHandle(STD_OUTPUT_HANDLE);
	for (unsigned i = 0; i < 5; i++)
	{
		int percent = (int)((((float)histogram[i]) / ((float)numEles))*100.0);
		if (bestBucket == i)
		{
			SetConsoleTextAttribute(console, FOREGROUND_RED | FOREGROUND_INTENSITY);
		}
		cout << "    ";
		for (int j = 0; j < percent; j++)
		{
			cout << "#";
		}
		cout << endl;
		if (bestBucket == i)
		{
			SetConsoleTextAttribute(console, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
		}
	}
	//
}

//Print the graph (for sanity checks)
void output(vector<LPEdge> &e)
{
	for (unsigned i = 0; i < e.size(); i++)
		std::cout << "From " << e[i].end0 << " to " << e[i].end1 << " with weight " << e[i].weight << endl;
}

//Run Column Generation on the graph [g], using the Held-Karp solution [lp]. Terminate if the algorithm 
//has not improved the sum of the slacks by more than [epsilon] over the last [average_over_last] iterations.
//Set [average_over_last] to 1 to run to optimality.
//Fills [num_its_required] with the number of iterations executed.
//Returns a pair containing the tree and the Yt values.
pair<vector<vector<LPEdge>>, vector<double>> bestOfMany(Graph *g, vector<LPEdge> &lp, double epsilon, int average_over_last, int *num_its_required, string *plot)
{
	//plot = "";
	double last_opt = -1.0;
	int iterator_count = 0;
	bool counter_initialized = false;

	//get numNodes used in Held-Karp soln.
	int numNodes = 0;
	for (unsigned i = 0; i < lp.size(); i++)
		numNodes = max(numNodes, max(lp[i].end0, lp[i].end1));
	numNodes++;
	//
	int numEdges = lp.size();
	
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	vector<vector<LPEdge>> trees;
	vector<vector<bool>> treeContains;
	double mag = 1.0;
	vector<double> treeVars;
	GRBModel curModel = GRBModel(env);
	//Create a slack variable for each edge
	vector<GRBVar> curSlacks;
	for (int i = 0; i < numEdges; i++)
	{
		GRBVar s = curModel.addVar(0.0, 1e20 + 1, 0.0, GRB_CONTINUOUS);
		curSlacks.push_back(s);
	}
	//
	vector<GRBVar> curYTs;
	curModel.update();
	GRBLinExpr curObj = 0;
	for (int i = 0; i < numEdges; i++)
		curObj += curSlacks[i];
	curModel.setObjective(curObj, GRB_MINIMIZE);
	vector<GRBConstr> curEdgeConstraints;
	for (int i = 0; i < numEdges; i++)
	{
		GRBLinExpr expr = GRBLinExpr(curSlacks[i]);
		//curEdgeConstraints.push_back(curModel.addConstr(expr == lp[i].weight));
		curEdgeConstraints.push_back(curModel.addConstr(expr == (((double)(numNodes-1))*lp[i].weight/((double)numNodes))));
	}
	double initialObjective = -1.0;
	//Solve the LP iteratively
	std::cout << "Iteratively solving the LP..." << endl;
	while (mag > 0)
	{
		(*num_its_required)++;
		//Add a variable for the newest tree, affecting column of constraints c which
		//references those constraints which correspond to the edges of the newest tree.
		if (treeContains.size() != 0)
		{
			int lastInd = (int)treeContains.size() - 1;
			GRBColumn c = GRBColumn();
			for (int i = 0; i < numEdges; i++)
			{
				if (treeContains[lastInd][i])
				{
					c.addTerm(1.0, curEdgeConstraints[i]);
				}
			}
			GRBVar yT = curModel.addVar(0.0, 1e20 + 1, 0.0, GRB_CONTINUOUS, c);
			curYTs.push_back(yT);
		}
		curModel.update();
		curModel.optimize();
		if (treeContains.size() == 0)
			initialObjective = curModel.get(GRB_DoubleAttr_ObjVal);
		//Now get dual variables.
		vector<double> curDuals;
		for (int i = 0; i < numEdges; i++)
			curDuals.push_back(curEdgeConstraints[i].get(GRB_DoubleAttr_Pi));
		//
		//run the maximum spanning tree, and then calculate its length.
		vector<LPEdge> curMSTIn;
		for (int i = 0; i < numEdges; i++)
		{
			LPEdge cur = LPEdge(lp[i].end0, lp[i].end1, curDuals[i]);
			curMSTIn.push_back(cur);
		}
		vector<LPEdge> curTree = runMaxSpanTree(curMSTIn, numNodes);
		mag = edgeSum(curTree);
		trees.push_back(curTree);
		vector<bool> curTContains;
		for (int i = 0; i < numEdges; i++)
		{
			bool found = false;
			for (unsigned j = 0; j < curTree.size(); j++)
			{
				if ((curTree[j].end0 == lp[i].end0 && curTree[j].end1 == lp[i].end1)
					|| (curTree[j].end1 == lp[i].end0 && curTree[j].end0 == lp[i].end1))
					found = true;
			}
			curTContains.push_back(found);
		}
		treeContains.push_back(curTContains);
		
		//cout << "iteration #" << *num_its_required << endl;
		//cout << "Average over last " << average_over_last << " iterations: " << running_average(last_optimals, average_over_last, *num_its_required) << endl;
		//cout << "Current value: " << curModel.get(GRB_DoubleAttr_ObjVal) << endl;
		//cout << "Min value: " << min_value(last_optimals, average_over_last, *num_its_required) << endl;

		//Track the iterative progress for early termination
		iterator_count++;
		if (last_opt - curModel.get(GRB_DoubleAttr_ObjVal) >= epsilon || !counter_initialized)
		{
			iterator_count = 0;
			last_opt = curModel.get(GRB_DoubleAttr_ObjVal);
			counter_initialized = true;
		}
		//
		(*plot) += to_string((*num_its_required)) + " " + to_string(curModel.get(GRB_DoubleAttr_ObjVal)) + "\n";
		if (mag <= 0.0 || iterator_count == average_over_last) //cutoff Reached 
		{
			if (iterator_count == average_over_last)
				std::cout << "Objective cutoff reached at " << curModel.get(GRB_DoubleAttr_ObjVal) << endl;
			for (unsigned i = 0; i < trees.size() - 1; i++)
				treeVars.push_back(curYTs[i].get(GRB_DoubleAttr_X));
			break;
		}
	}
	std::cout << "Done iteratively solving the LP after " << *num_its_required << " iterations." << endl;
	std::cout << "Final number of trees: " << trees.size() << endl;


	vector<int> treesToCheck;
	for (unsigned i = 0; i < trees.size() - 1; i++)
	{
		if (treeVars[i] > 0.0)
			treesToCheck.push_back(i);
	}

	std::cout << "Number of positive trees: " << treesToCheck.size() << endl;

	int best = INT_MAX;
	int bestIndex = -1;
	vector<vector<LPEdge>> posTrees;
	vector<double> posYTs;
	posYTs.reserve(treesToCheck.size());
	posTrees.reserve(treesToCheck.size());
	for (unsigned i = 0; i < treesToCheck.size(); i++)
	{
		posTrees.push_back(trees[treesToCheck[i]]);
		posYTs.push_back(treeVars[treesToCheck[i]]);
	}
	return pair<vector<vector<LPEdge>>, vector<double>>(posTrees, posYTs);
}
