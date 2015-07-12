#ifndef PRINT_H
#define PRINT_H

#include "definitions.h"

/// Print a splitting to std::cout
void output(const splitting &in);

/// Print a vector<JoinEdge> to std::cout
void output(const vector<JoinEdge> &in);
void output(const vector<vector<JoinEdge>> &in, string name);

void print(const vector<Edge> &t);
void print(const vector<int> &l);

#endif