#include "matrix.h"
#include "graph.h"
#include "util.h"

//g is lambda-weighted
//    L_i,j = { -lambda_e :  e = (i,j) \in E
//			  { Sum(where e \in delta ({i})) lambda_e : i = j
//			  { 0 : o/w
//O(ve)
Eigen::SparseMatrix<double, Eigen::ColMajor> wLaplacianSparse(const vector<ModEdge> &lGraph, const vector<int> &cNodes)
{
	//Prep; optimize away in future
	vector<ModEdge> adjG = adjustGraph(lGraph, cNodes);
	int numNodes = getNumNodes(adjG);
	//

	vector<entry> builder;
	//Handle non-diagonal entries.
	for (unsigned i = 0; i < adjG.size(); i++)
	{
		int u = adjG[i].end0;
		int v = adjG[i].end1;
		builder.push_back(entry(u, v, adjG[i].weight*-1.0));
		builder.push_back(entry(v, u, adjG[i].weight*-1.0));
	}
	//

	//Handle diagonal entries
	//each diagonal entry is the sum of all the lambdas corresponding to the
	//edges with one end in that node.
	vector<double> diags(numNodes, 0.0);
	for (unsigned i = 0; i < adjG.size(); i++)
	{
		diags[adjG[i].end0] += adjG[i].weight;
		diags[adjG[i].end1] += adjG[i].weight;
	}

	for (int v = 0; v < numNodes; v++)
		builder.push_back(entry(v, v, diags[v]));
	//
	//cout << "The matrix: " << endl;
	//cout << lpcn << endl;
	Eigen::SparseMatrix<double, Eigen::ColMajor> lpcn(numNodes, numNodes);
	lpcn.setFromTriplets(builder.begin(), builder.end());
	return lpcn;
}

//*This is actually a minor, not a cofactor.
double cofactorSparse(const Eigen::SparseMatrix<double, Eigen::ColMajor> &m)
{
	if (m.nonZeros() == 0)
		return 0.0;
	Eigen::SparseMatrix<double, Eigen::ColMajor> b = m.block(0, 0, m.rows() - 1, m.cols() - 1);
	b.makeCompressed();
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double, Eigen::ColMajor>> cholesky(b);
	double det = cholesky.determinant();
	//cout << "Det: " << det << endl;
	//if (abs(det) > 1e250)
	//	cout << "Warning: Determinant nearing overflow!" << endl;
	return det;
}

double cofactorLog(const Eigen::SparseMatrix<double, Eigen::ColMajor> &m)
{
	if (m.nonZeros() == 0)
		return 0.0;
	Eigen::SparseMatrix<double, Eigen::ColMajor> b = m.block(0, 0, m.rows() - 1, m.cols() - 1);
	Eigen::MatrixXd d(b);
	return d.householderQr().logAbsDeterminant();
}

Eigen::SparseMatrix<double, Eigen::ColMajor> wLaplacianNew(const vector<ModEdge> &lGraph)
{
	//get max ind
	int maxInd = 0;
	for (unsigned i = 0; i < lGraph.size(); i++)
		maxInd = max(lGraph[i].end0, max(lGraph[i].end1, maxInd));
	maxInd++; //Note that maxInd is actually the number of nodes in the original graph.
	//

	//figure out which indices are used
	vector<bool> usedIndices(maxInd, false);
	for (unsigned i = 0; i < lGraph.size(); i++)
	{
		usedIndices[lGraph[i].end0] = true;
		usedIndices[lGraph[i].end1] = true;
	}
	//

	//Now, create a map from old indices to new indices...
	//newInd[oldInd] = the new number for the number oldInd.
	vector<int> newInd(maxInd, -1);
	for (int i = 0; i < maxInd; i++)
	{
		//Check if the current index is used?
		if (usedIndices[i])
		{
			int offset = 0;
			//It is, so we have to figure out how many numbers below it are unused (Optimize this in the future)
			for (int j = 0; j < i; j++)
			{
				//For each index below i, if that index is unused, it is a free slot we need to account for, so increase the offset.
				if (!usedIndices[j])
					offset++;
			}
			//Ex: i = 5, which is node 5. Nodes 2 and 3 are unused, so the offset is 2 for Node 5. The new index for Node 5 is 3; the new index for 4 will be 2.
			newInd[i] = i - offset;
		}
	}
	//

	//Get the number of actual nodes:
	int numNodes = 0;
	for (unsigned i = 0; i < usedIndices.size(); i++)
		numNodes += usedIndices[i]; //bool -> int : false -> 0, true -> 1

	//Now we have the mapping from old indices to new indices, so we will temporarily reorder all the nodes, then construct the
	// laplacian from that (then reorder again and check against the index backups).

	vector<ModEdge> g = lGraph;
	for (unsigned i = 0; i < g.size(); i++)
	{
		g[i].end0 = newInd[g[i].end0];
		g[i].end1 = newInd[g[i].end1];
	}

	//Finally time to reconstruct the laplacian:


	vector<entry> builder;
	//Handle non-diagonal entries.
	for (unsigned i = 0; i < g.size(); i++)
	{
		int u = g[i].end0;
		int v = g[i].end1;
		builder.push_back(entry(u, v, g[i].weight*-1.0));
		builder.push_back(entry(v, u, g[i].weight*-1.0));
	}
	//

	//Handle diagonal entries
	//each diagonal entry is the sum of all the lambdas corresponding to the
	//edges with one end in that node.
	vector<double> diags(numNodes, 0.0);
	for (unsigned i = 0; i < g.size(); i++)
	{
		diags[g[i].end0] += g[i].weight;
		diags[g[i].end1] += g[i].weight;
	}

	for (int v = 0; v < numNodes; v++)
		builder.push_back(entry(v, v, diags[v]));
	//
	//cout << "About to set laplacian with " << numNodes << " actual nodes" << endl;
	Eigen::SparseMatrix<double, Eigen::ColMajor> lpcn(numNodes, numNodes);
	lpcn.setFromTriplets(builder.begin(), builder.end());
	//cout << "The matrix: " << endl;
	//cout << lpcn << endl;
	return lpcn;
}

Eigen::SparseMatrix<double, Eigen::ColMajor> blocked(const Eigen::SparseMatrix<double, Eigen::ColMajor> &m)
{
	return m.block(0, 0, m.rows() - 1, m.cols() - 1);
}

Eigen::MatrixXd inv(const Eigen::SparseMatrix<double, Eigen::ColMajor> &m)
{
	return Eigen::MatrixXd(m).inverse();
}

double xTLInvX(const Eigen::MatrixXd &lInv, int u, int v)
{
	int dim = lInv.cols();
	Eigen::VectorXd x(dim);
	for (int i = 0; i < dim; i++)
		x[i] = 0.0;
	if (u < dim)
		x[u] = 1.0;
	if (v < dim)
		x[v] = -1.0;
	return (x.transpose() * lInv) * x;
}


double prob(const vector<ModEdge> &lGraph, ModEdge e)
{
	Eigen::SparseMatrix<double> laplacian = blocked(wLaplacianNew(lGraph));
	//need to find proper indexes u and v into laplacian inverse
	//get max ind
	int maxInd = 0;
	for (unsigned i = 0; i < lGraph.size(); i++)
		maxInd = max(lGraph[i].end0, max(lGraph[i].end1, maxInd));
	maxInd++; //Note that maxInd is actually the number of nodes in the original graph.
	//

	//figure out which indices are used
	vector<bool> usedIndices(maxInd, false);
	for (unsigned i = 0; i < lGraph.size(); i++)
	{
		usedIndices[lGraph[i].end0] = true;
		usedIndices[lGraph[i].end1] = true;
	}
	//

	//Now, create a map from old indices to new indices...
	//newInd[oldInd] = the new number for the number oldInd.
	vector<int> newInd(maxInd, -1);
	for (int i = 0; i < maxInd; i++)
	{
		//Check if the current index is used?
		if (usedIndices[i])
		{
			int offset = 0;
			//It is, so we have to figure out how many numbers below it are unused (Optimize this in the future)
			for (int j = 0; j < i; j++)
			{
				//For each index below i, if that index is unused, it is a free slot we need to account for, so increase the offset.
				if (!usedIndices[j])
					offset++;
			}
			//Ex: i = 5, which is node 5. Nodes 2 and 3 are unused, so the offset is 2 for Node 5. The new index for Node 5 is 3; the new index for 4 will be 2.
			newInd[i] = i - offset;
		}
	}
	//
	return e.weight * xTLInvX(inv(laplacian), newInd[e.end0], newInd[e.end1]);
}