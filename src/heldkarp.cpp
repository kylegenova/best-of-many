//Implementation for the concorde interface to the Held-Karp solutions.
#include "heldkarp.h"

vector<LPEdge> getHKSolutions(string fName, bool silent)
{
	string outFName = fName.substr(0, fName.length() - 4) + ".txt";
	//Solve Subtour LP relaxation using concorde command-line tool, then read in the results.
	std::size_t isTSV = fName.find(".tsv");
	string command = "";
	if (isTSV != std::string::npos)
	{
		string corrected = fName.substr(0, isTSV) + ".tsp";
		command = "concorde.exe -I -s 1 -x -X lp " + corrected + " > log.txt";
	}
	else
	{
		command = "concorde.exe -I -s 1 -x -X lp " + fName + " > log.txt";
	}

	system(command.c_str());
	

	//Read results
	ifstream f;
	vector<LPEdge> lp;
	int numNodes;
	int numEdges;
	f.open("lp", ifstream::in);
	string curLine;
	getline(f, curLine);
	int split = curLine.find(" ");
	numNodes = atoi(curLine.substr(0, split).c_str());
	numEdges = atoi(curLine.substr(split + 1, string::npos).c_str());
	//std::cout << "Number of Nodes: " << numNodes << "\n";
	//std::cout << "Number of Edges: " << numEdges << "\n";
	//cout << "This is the input HK Solution in string format:" << endl;
	for (int i = 0; i < numEdges; i++)
	{
		getline(f, curLine);
		//cout << curLine << endl;
		int split1 = curLine.find(" ");
		int split2 = curLine.find_last_of(" ");
		string fst = curLine.substr(0, split1);
		string snd = curLine.substr(split1 + 1, split2);
		string tl = curLine.substr(split2 + 1, string::npos);
		int e0 = atoi(fst.c_str());
		int e1 = atoi(snd.c_str());
		double w = atof(tl.c_str());
		lp.push_back(LPEdge(e0, e1, w));
	}
	return lp;
}

vector<LPEdge> getHKSolutionsBaked(string fName, bool silent)
{
	cout << "Entered getHKSolutionsBaked" << endl;
	fName = fName.substr(0, fName.length() - 4) + ".lp";
	ifstream f;
	f.open(fName);
	vector<LPEdge> lp;
	string line;
	while (!f.eof())
	{
		int c1 = -1;
		int c2 = -1;
		double v = -1.;
		f >> c1 >> c2 >> v;
		if (v >= 0.01)
		{
			cout << "c1: " << c1 << ", c2: " << c2 << ",v: " << v << endl;
			lp.push_back(LPEdge(c1, c2, v));
		}
	}
	f.close();
	return lp;
}