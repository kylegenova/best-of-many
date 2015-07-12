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

//   lambda_e := e^gamma_e

//    gamma' := gamma'_e = gamma_e - delta,  gamma'_f = gamma_f \forall f \in E \ {e}

//   While there exists and edge e s.t. q_e(gamma) > (1 + epsilon)z_e, compute gamma', set gamma := gamma'.


#include "gamma.h"
#include "util.h"
#include "matrix.h"

double EPSILON = 0.01;

vector<double> gammaToLambda(const vector<double> &gam)
{
	vector<double> lambdas;
	for (unsigned i = 0; i < gam.size(); i++)
		lambdas.push_back(exp(gam[i]));
	return lambdas;
}

double q_e(int e, const vector<ModEdge> &zGraph, vector<double> gamma, double denom=0.0)
{
	vector<ModEdge> lGraph = zGraph;
	vector<double> lambda = gammaToLambda(gamma);
	weightGraph(lGraph, lambda);
	denom = denom == 0.0 ? cofactorSparse(wLaplacianNew(lGraph)) : denom;
	return lambda[e] * abs(cofactorSparse(wLaplacianNew(contractEdgeNew(lGraph, lGraph[e]))) / denom);
}

int getUpdateLoc(const vector<double> &gamma, const vector<ModEdge> &zGraph)
{
	vector<double> lambda = gammaToLambda(gamma);
	vector<ModEdge> lGraph = zGraph;
	weightGraph(lGraph, lambda);
	double denom = cofactorSparse(wLaplacianNew(lGraph));
	for (unsigned e = 0; e < zGraph.size(); e++)
	{
		if (q_e(e, zGraph, gamma, denom) > (1.0 + EPSILON)*zGraph[e].weight)
			return e;
	}
	return -1;
}

double delta(const vector<double> &gamma, const vector<ModEdge> &zGraph, int e)
{
	vector<double> lambda = gammaToLambda(gamma);
	vector<ModEdge> lGraph = zGraph;
	weightGraph(lGraph, lambda);
	double q = q_e(e, zGraph, gamma);
	double num = q * (1.0 - ((1.0 + (EPSILON / 2.0))*zGraph[e].weight));
	double denom = (1.0 - q)*(1.0 + (EPSILON / 2.0))*zGraph[e].weight;
	return log(num / denom);
}

vector<double> getLambda(const vector<ModEdge> &zGraph)
{
	vector<double> gamma(zGraph.size(), 0.0);
	while (true)
	{
		int loc = getUpdateLoc(gamma, zGraph);
		if (loc == -1)
			break;
		gamma[loc] -= delta(gamma, zGraph, loc);
	}
	return gammaToLambda(gamma);
}

vector<double> qes(const vector<ModEdge> &lGraph)
{
	Eigen::SparseMatrix<double> laplacian = blocked(wLaplacianNew(lGraph));
	//Eigen::SparseMatrix<double> laplacian = wLaplacianNew(lGraph);
	Eigen::MatrixXd invLap = inv(laplacian);
	vector<double> qes;
	for (unsigned i = 0; i < lGraph.size(); i++)
	{
		double qe = lGraph[i].weight * xTLInvX(invLap, lGraph[i].end0, lGraph[i].end1);
		qes.push_back(qe);
	}
	return qes;
}

double delta(double z, double q)
{
	double num = q * (1.0 - ((1.0 + (EPSILON / 2.0))*z));
	double denom = (1.0 - q)*(1.0 + (EPSILON / 2.0))*z;
	return log(num / denom);
}

vector<double> getLambdaFast(const vector<ModEdge> &zGraph)
{
	vector<double> gamma(zGraph.size(), 0.0);
	while (true)
	{
		bool errorExceeds = false;
		vector<ModEdge> lGraph = zGraph;
		weightGraph(lGraph, gammaToLambda(gamma));
		vector<double> qe = qes(lGraph);
		for (unsigned i = 0; i < qe.size(); i++)
		{
			if (qe[i] > (1.0 + EPSILON)*zGraph[i].weight)
			{	
				gamma[i] -= delta(zGraph[i].weight, qe[i]);
				errorExceeds = true;
			}
		}
		if (!errorExceeds)
			break;
	}
	vector<double> lambda = gammaToLambda(gamma);
	//output(lambda);
	ofstream backup;
	backup.open("lpbackup.lp", std::ofstream::trunc);
	for (unsigned i = 0; i < lambda.size(); i++)
		backup << zGraph[i].end0 << " " << zGraph[i].end1 << " " << lambda[i] << endl;
	backup.close();
	return lambda;
}