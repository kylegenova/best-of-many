
#include "print.h"
#include "edge_join.h"


void output(const splitting &in)
{
	for (unsigned i = 0; i < in.size(); i++)
	{
		int s = in[i].first;
		cout << "For node " << s << ": " << endl;
		vector<JoinEdge> es = in[i].second;
		for (unsigned j = 0; j < es.size(); j++)
		{
			cout << "    From " << es[j].u << " to " << es[j].v << endl;
		}
	}
}

void output(const vector<JoinEdge> &in)
{
	for (unsigned i = 0; i < in.size(); i++)
	{
		cout << "From " << in[i].u << " to " << in[i].v << endl;
	}
}

void output(const vector<vector<JoinEdge>> &in, string name)
{
	for (unsigned i = 0; i < in.size(); i++)
	{
		cout << name << " #" << i << ": " << endl;
		for (unsigned j = 0; j < in[i].size(); j++)
		{
			cout << "    From " << in[i][j].u << " to " << in[i][j].v << endl;
		}

	}
}

void print(const vector<Edge> &t)
{
	for (unsigned i = 0; i < t.size(); i++)
	{
		cout << "From " << t[i].u << " to " << t[i].v << endl;
	}
}

void print(const vector<int> &l)
{
	for (unsigned i = 0; i < l.size(); i++)
	{
		cout << l[i] << endl;
	}
}