#ifndef INCLUDE_BAYESNET_CPP
#define INCLUDE_BAYESNET_CPP

#include "Nde.h"
#include <cstdlib>
#include <fstream>
#include <map>
#include <cassert>

struct BayesNet
{
	long no_of_nodes;
	vector <Node*> nodes;

	map <long, string> id2str;
	map <string, long> str2id;

	BayesNet()
	{
		no_of_nodes = 0;
	}

};

void insertNode(BayesNet *b, string str, long dist_values, vector<string> vals)
{
	long id = b->no_of_nodes;
	b->no_of_nodes += 1;

	Node *node = createNode(id, str, dist_values, vals);
	b->nodes.push_back(node);
	b->str2id[str] = id;
	b->id2str[id] = str;
}

void addParents(BayesNet *b, vector<string> vars)
{
	long len = vars.size();
	if (len == 1) return;

	long child_id = b->str2id[vars[len - 1]];

	for(long i = 0; i < len-1; i++)
	{
		long par_id = b->str2id[vars[i]];

		addParent(b->nodes[child_id], b->nodes[par_id]);
		addChild(b->nodes[par_id], b->nodes[child_id]);
	}

}

void convertCPD2Factor(BayesNet *b)
{
	long n = b->no_of_nodes;

	for (long i = 0; i < n; i++)
	{
		vector <long> scp;
		map <long, long> scpMap;

		for (auto j : b->nodes[i]->parents)
		{
			scp.push_back(j->id);
			scpMap[j->id] = j->dist_values;
		}
		scp.push_back(b->nodes[i]->id);
		scpMap[b->nodes[i]->id] = b->nodes[i]->dist_values;
		b->nodes[i]->factorCPD = getFactorCPD(b->nodes[i]->CPD, scp, scpMap);
	}
}

void printToFile(BayesNet *bn, ofstream &fout)
{
	fout << "network unknown{}" << endl;
	long n = bn->no_of_nodes;

	for (long i = 0; i < n; i++)
	{
		long sz = bn->nodes[i]->attr.size();
		fout << "variable " << bn->nodes[i]->str << " " << "{" << endl;
		fout << "type discrete [ " << bn->nodes[i]->dist_values << " ] { ";

		for (long a = 0; a < sz - 1; a++)
			fout << bn->nodes[i]->attr[a] << ", ";
		fout << bn->nodes[i]->attr[sz - 1] << " };" << endl; 
		fout << "}" << endl;
	}

	for(long i = 0; i < n; i++)
	{
		fout << "probability ( " << bn->id2str[bn->nodes[i]->id]; 
		long psize = bn->nodes[i]->parents.size();

		if (psize != 0)
		{
			fout << " | ";
			for (long p = 0; p < psize - 1; p++)
			{
				fout << bn->id2str[bn->nodes[i]->parents[p]->id] << ", ";
			} 
			fout << bn->id2str[bn->nodes[i]->parents[psize - 1]->id] << " ) {" << endl;
		}
		else
		{
			fout << " ) {" << endl;
		}
		Factor *f = bn->nodes[i]->factorCPD;
		long dValues = bn->nodes[i]->dist_values;
		for (long i1 = 0; i1 < f->table.size(); i1+=dValues)
		{
			long sz = f->table[i1].first.size();

			if (sz > 1)
			{
				fout << "(";
				for (long j = 0; j < sz - 2; j++)
				{
					long c = f->table[i1].first[j];
					fout << bn->nodes[i]->parents[j]->value2str[c] << ", ";
				}
				fout << bn->nodes[i]->parents[sz - 2]->value2str[f->table[i1].first[sz-2]] << ") ";
			}
			else
				fout << "table ";

			for (long m = 0; m < dValues - 1; m++)
			{
				fout << exp(f->table[i1 + m].second) << ", ";
			}
			fout << exp(f->table[i1 + dValues - 1].second) <<";" << endl;
		}
		fout << "}" << endl;
	}
}

void printToFile2(BayesNet *bn, ofstream &fout)
{
	fout << "network unknown{}" << endl;
	long n = bn->no_of_nodes;

	for (long i = 0; i < n; i++)
	{
		long sz = bn->nodes[i]->attr.size();
		fout << "variable " << bn->nodes[i]->str << " " << "{" << endl;
		fout << "type discrete [ " << bn->nodes[i]->dist_values << " ] { ";

		for (long a = 0; a < sz - 1; a++)
			fout << bn->nodes[i]->attr[a] << ", ";
		fout << bn->nodes[i]->attr[sz - 1] << " };" << endl; 
		fout << "}" << endl;
	}

	for (long i = 0; i < n; i++)
	{
		fout << "probability ( ";
		for (long p = 0; p < bn->nodes[i]->parents.size(); p++)
		{
			fout << bn->id2str[bn->nodes[i]->parents[p]->id] << ", ";
		} 
		fout << bn->id2str[bn->nodes[i]->id] << " ) {" << endl;

		Factor *f = bn->nodes[i]->factorCPD;
		for (long j = 0; j < f->table.size(); j++)
		{
			fout << "(";
			for (long j1 = 0; j1 < f->table[j].first.size(); j1++)
			{
				long c = f->table[j].first[j1];
				if (j1 == f->table[j].first.size() - 1)
					fout << bn->nodes[i]->value2str[c] << ") ";
				else	
					fout << bn->nodes[i]->parents[j1]->value2str[c] << ", ";			
			}
			fout << exp(f->table[j].second) <<";" << endl;
		}
		fout << "}" << endl;
	}
	 
}
void printBayesNet(BayesNet *b, ofstream &fout, long num)
{
	long n = b->no_of_nodes;
	cout << n << endl;

	if(num == 1)
		printToFile2(b, fout);
	else
		printToFile(b, fout);

	for(long i = 0; i < n; i++)
	{
		Node *node = b->nodes[i];
		vector<Node*> parents = node->parents;
		cout << node->str << " " << "[";
		for(long j = 0; j < parents.size(); j++)
			cout << parents[j]->str << ((j == parents.size() - 1) ? "" : ",");
		cout << "]" << endl;
		//displayCPD(node);
		printFactor(node->factorCPD);
	}
}

#endif