#ifndef MATRIX_H
#define MATRIX_H

#include "definitions.h"

#include <Eigen\Dense>
#include <Eigen\Sparse>

typedef Eigen::Triplet<double> entry;

Eigen::SparseMatrix<double, Eigen::ColMajor> wLaplacianSparse(const vector<ModEdge> &lGraph, const vector<int> &cNodes);
Eigen::SparseMatrix<double, Eigen::ColMajor> wLaplacianNew(const vector<ModEdge> &lGraph);
double cofactorSparse(const Eigen::SparseMatrix<double, Eigen::ColMajor> &m);
Eigen::MatrixXd inv(const Eigen::SparseMatrix<double, Eigen::ColMajor> &m);
double xTLInvX(const Eigen::MatrixXd &lInv, int u, int v);
Eigen::SparseMatrix<double, Eigen::ColMajor> blocked(const Eigen::SparseMatrix<double, Eigen::ColMajor> &m);
double prob(const vector<ModEdge> &lGraph, ModEdge e);
double cofactorLog(const Eigen::SparseMatrix<double, Eigen::ColMajor> &m);
#endif