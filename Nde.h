#ifndef INCLUDE_NODE_CPP
#define INCLUDE_NODE_CPP

#include <string>
#include <vector>
#include "Factor.cpp"
#include <iomanip>
#include <iostream>
#include <cassert>
#include <map>

using namespace std;


struct Node
{
	long id;
	string str;
	vector <Node*> children;
	vector <Node*> parents;

	map <vector <long>, double> CPD;
	Factor *factorCPD;
	vector <string> attr;
	map <long, string> value2str;
	map <string, long> str2value;
	long dist_values;

};

void initialize_table(long N, vector <vector <string> > &assignment, vector <string> list, long origN, 
		Node *n)
{
	if(N <= 0)
	{
		long size = list.size();
		vector <long> vec;
		for (long i = 0; i < size - 1; i++)
		{
			vec.push_back(n->parents[i]->str2value[list[i]]);
		}
		vec.push_back(n->str2value[list[size - 1]]);
		//n->CPD[list] = 0.0;
		n->CPD[vec] = 0.0;
		return;
	}
	else
	{
		for (long i = 0; i < assignment[origN - N].size(); i++)
		{
			list.push_back(assignment[origN - N][i]);
			initialize_table(N - 1, assignment, list, origN, n);
			list.pop_back();
		}
	}

}

void createTable(Node *n)
{
	vector <string> vec;
	vector <vector <string> > assignment;

	for (auto p : n->parents)
	{
		assignment.push_back(p->attr);
	}

	assignment.push_back(n->attr);
	initialize_table(assignment.size(), assignment, vec, assignment.size(), n);
}

Node *createNode(long _id, string _str, long _dist_values, vector<string> _attr)
{
	Node *n = new Node();
	n->id = _id;
	n->str = _str;
	n->dist_values = _dist_values;
	n->attr = _attr;

	for (long a = 0; a < _attr.size(); a++)
	{
		n->str2value[_attr[a]] = a;
		n->value2str[a] = _attr[a];
	}

	createTable(n);
	return n;
}

void addChild(Node *node, Node *child)
{
	node->children.push_back(child);
}

void addParent(Node *node, Node *parent)
{
	node->parents.push_back(parent);
	node->CPD.clear();
	createTable(node);
}

void displayCPD(Node *n)
{
	for (auto p : n->parents)
		cout << p->str << " ";
	cout << n->str << endl;

	for (auto m : n->CPD)
	{
		for (auto a : m.first)
			cout << a << " ";
		cout << m.second << endl;  
	}
}

#endif