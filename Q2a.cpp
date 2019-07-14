#include "parser.h"
#include <sstream>
#include <cassert>
#include <stdio.h>
#include <math.h>
using namespace std;


double LikeLihood = 0.0;
vector <double> accVector, likeVector;
long correctans = 0, totalans = 0;

bool issubset(long id, vector <long> scp)
{
	for (long i = 0; i < scp.size(); i++)
	{
		if (id == scp[i]) return true;
	}
	return false;
}

long getMaxvalGibbs(long attrId, vector <map <long, long> > &samples_gen, long attrVal, ofstream &outfile,
	BayesNet *bn, double &L)
{
	long assign;
	long maxVal = -1;
	map <long, long> val;
	for (auto sample : samples_gen)
	{
		if (val.find(sample[attrId]) == val.end()) 
			val[sample[attrId]] = 0;

		val[sample[attrId]]++;
	}
	for (auto v : val)
	{
		if (v.second > maxVal)
		{
			maxVal = v.second;
			assign = v.first;
		}
		outfile << bn->nodes[attrId]->value2str[v.first] << ":" <<(v.second/(double)samples_gen.size()) << " ";
	}
	if (val[attrVal] == 0) val[attrVal] = 1;
	LikeLihood += log((val[attrVal])/(double)(samples_gen.size()) ); // Likelihood Computation....
	L += log((val[attrVal])/(double)(samples_gen.size()) ); 
	return assign;
} 

map<long, long> Gibbspredict(vector <map <long, long> > &samples_gen, vector <long> &variables,
	map <long, long> &UnobservedMap, ofstream &outfile, BayesNet *bn)
{
	double L = 0.0;
	map <long, long> assign;

	for (auto v : variables)
	{
		outfile << bn->id2str[v] << " ";
		assign[v] = getMaxvalGibbs(v, samples_gen, UnobservedMap[v], outfile, bn, L);
		outfile << endl;
	}
	if (variables.size() != 0)
		likeVector.push_back(L);
	outfile << endl;
	return assign;
}

map <long, long> getnewSample(map <long, long> sample, vector <Factor*> &initial_facts, long turn,
	vector <long> &variables)
{
	double Z = 0;
	long assign;
	vector <Factor*> fact;

	for (auto f : initial_facts)
	{
		if (issubset(variables[turn], f->scope)) fact.push_back(f);
	}

	for (long f = 0; f < fact.size(); f++)
	{
		for (auto s : sample)
		{
			long key = s.first;
			long val = s.second;
			if (key != variables[turn] && issubset(key, fact[f]->scope))
				fact[f] = reduction(fact[f], key, val, 10);
		}
	}
	Factor *Fact_mul = compute_psinew(fact, 10);
	Fact_mul = normalize(Fact_mul);
	
	double r = ((double) rand() / (RAND_MAX));
	double CF = 0;

	for (long i = 0; i < Fact_mul->table.size(); i++)
	{
		CF += exp(Fact_mul->table[i].second);
		assign = Fact_mul->table[i].first[0];

		if (r <= CF) break;			
	}
	sample[variables[turn]] = assign;
	return sample;
}



map <long, long> GibbsSampler(vector <Factor*> &initial_facts, long burn_insamples, long gen_samples
	,vector <long> &variables, map <long, long> &UnobservedMap, ofstream &outfile, BayesNet *bn)
{
	vector < map <long, long>> samples_gen;
	map <long, long> samples, ans;
	for (long i = 0; i < variables.size(); i++)
		samples[variables[i]] = 0;

	long turn = 0;
	while(burn_insamples--)
	{
		samples = getnewSample(samples, initial_facts, turn, variables);
		turn = (turn + 1) % variables.size();
	}
	
	while(gen_samples--)
	{
		samples = getnewSample(samples, initial_facts, turn, variables);
		samples_gen.push_back(samples);
		turn = (turn + 1) % variables.size();
	}
	ans = Gibbspredict(samples_gen, variables, UnobservedMap, outfile, bn);
	return ans;
}

vector <long> getAttribute(Node *n)
{
	vector <long> s;
	for (auto p : n->parents)
	{
		s.push_back(p->id);
	}
	s.push_back(n->id);
	return s;
}

bool allEqual(map <long, long> &assignVal, Node *n, vector <long> &nums)
{
	long sz = n->parents.size();
	for(long i = 0; i < n->parents.size(); i++)
	{
		if ( assignVal[n->parents[i]->id] != nums[i]) return false;
	}
	if (assignVal[n->id] != nums[sz]) return false;

	return true;
}

void updateTable(Node *n, map <long, long> &assignment)
{

	Factor *f = n->factorCPD;

	for (long i = 0; i < f->table.size(); i++)
	{
		if (allEqual(assignment, n, f->table[i].first))
		{
			f->table[i].second += 1;
			break;
		}
	}
}
void fillTable(BayesNet *bn, vector <long> &attributes, vector <long> &values)
{
	long N = bn->no_of_nodes;
	map <long, long> assignMap;

	for (long i = 0; i < attributes.size(); i++)
	{
		assignMap[attributes[i]] = values[i];
	}
	
	for(long n = 0; n < N; n++)
		updateTable(bn->nodes[n], assignMap);

}

void normalizeTable(BayesNet *bn)
{
	long N = bn->no_of_nodes;
	for (long n = 0; n < N; n++)
	{ 
		vector <long> tmp;
		Factor *f = bn->nodes[n]->factorCPD;
		map <vector <long>, long> attrValues;
		long dist_values = bn->nodes[n]->dist_values;

		for(long i = 0; i < f->table.size(); i++)
		{
			tmp = f->table[i].first;
			tmp.pop_back();

			if (attrValues.find(tmp) == attrValues.end())
				attrValues[tmp] = 0;

			attrValues[tmp] += f->table[i].second;
		}

		for(long i = 0; i < f->table.size(); i++)
		{
			tmp = f->table[i].first;
			tmp.pop_back();

			f->table[i].second = log((f->table[i].second + 1 )/(double)(attrValues[tmp] + dist_values));
		}

	}
}

void learnBayesNet(BayesNet *bn, string trainFile)
{
	string str, s;
	ifstream f;
	ofstream fout;

	f.open(trainFile);
	fout.open("./bayesian.bif");
	vector <long> attributes, values;

	getline(f, str);
	stringstream stringS(str);

	while(stringS >> s)
		attributes.push_back(bn->str2id[s]);

	while(getline(f, str))
	{
		stringstream stringS(str);

		for (long i = 0; i < attributes.size(); i++)
		{
			stringS >> s;
			values.push_back(bn->nodes[attributes[i]]->str2value[s]);
		}
		assert(attributes.size() == values.size());
		fillTable(bn, attributes, values);
		values.clear();
	}
	normalizeTable(bn);
	printBayesNet(bn, fout, 0);

}

void predictValues(BayesNet *bn, map <long, long> &ObservedMap, map <long, long> &UnobservedMap,
	vector <long> &variables, ofstream &outfile)
{
	long N = bn->no_of_nodes;
	vector <Factor*> facts;

	for(long n = 0; n < N; n++)
	{
		Factor *fCPD = bn->nodes[n]->factorCPD;
		vector <long> table_attr = getAttribute(bn->nodes[n]);

		for (auto t : table_attr)
		{
			if(ObservedMap.find(t) == ObservedMap.end()) 
				continue;
			fCPD = reduction(fCPD, t, ObservedMap[t], 10);
		}
		facts.push_back(fCPD);		
	}

	map <long, long> assignVals = GibbsSampler(facts, 1000, 10000, variables, UnobservedMap, outfile, bn);
	double acc = 0.0;
	for (auto m : assignVals)
	{
		long key = m.first;
		long val = m.second;

		if (UnobservedMap[key] == val) 
		{
			correctans++; 
			acc++;
		}
		totalans++;
	}

	accVector.push_back(acc/(double)assignVals.size());
}

void prediction(BayesNet *bn, string testFile, string trueFile)
{
	ifstream f1, f2;
	ofstream outfile;
	f1.open(testFile);
	f2.open(trueFile);

	string str1, str2, s1, s2, s;
	outfile.open("./query.bn.out");

	vector <long> attributes, variables;
	map <long, long> ObservedMap, UnobservedMap;

	getline(f1, str1);
	getline(f2, str2);
	stringstream stringS(str1);

	while(stringS >> s)
		attributes.push_back(bn->str2id[s]);

	while(getline(f1, str1))
	{
		getline(f2, str2);
		stringstream stringS1(str1);
		stringstream stringS2(str2);

		for (long i = 0; i < attributes.size(); i++)
		{
			stringS1 >> s1;
			stringS2 >> s2;
			if (s1 == "?") 
			{
				variables.push_back(attributes[i]);
				UnobservedMap[attributes[i]] = bn->nodes[attributes[i]]->str2value[s2];
			}
			else
			{
				ObservedMap[attributes[i]] = bn->nodes[attributes[i]]->str2value[s1];
			}
		}

		predictValues(bn, ObservedMap, UnobservedMap, variables, outfile);
		ObservedMap.clear();
		UnobservedMap.clear();
		variables.clear();
	}
}

void printAccuracy()
{
	cout << "Accuracy  and LogLikelihood is :" <<endl;
	double sum = 0.0, likeSum = 0.0;
	for (auto a : accVector)
		sum += a;
	double sum1 =  (accVector.size() > 0) ? (sum/(double)accVector.size()) : 0;
	for (auto l : likeVector)
		likeSum += l;
	double likeSum1 = (likeVector.size() > 0) ? (likeSum/(double)likeVector.size()) : 0;

	printf("%lf %lfs\n\n", sum1 * 100, likeSum1);

}


int main(int argc, char const *argv[])
{
	BayesNet *bn;
	long category = atoi(argv[1]);
	string trainFile, testFile, trueFile;

	if (category == 1) 
	{
		bn = parse("./A3-data/andes.bif");
		trainFile = "./A3-data/andes.dat";
		testFile = "./A3-data/andes_test.dat";
		trueFile = "./A3-data/andes_TrueValues.dat";
	}
	else if(category == 2) 
	{
		bn = parse("./A3-data/hepar2.bif");
		trainFile = "./A3-data/hepar2.dat";
		testFile = "./A3-data/hepar_test.dat";
		trueFile = "./A3-data/hepar_TrueValues.dat";
	}
	else if (category == 3) 
	{
		bn = parse("./A3-data/insurance.bif");
		trainFile = "./A3-data/insurance.dat";
		testFile = "./A3-data/insurance_test.dat";
		trueFile = "./A3-data/insurance_TrueValues.dat";
	}
	else if(category == 4)
	{
		bn = parse("./A3-data/dummy.bif");
		trainFile = "./A3-data/dummy.dat";
	}

	convertCPD2Factor(bn);

	cout << "Learning.........." << endl;
	learnBayesNet(bn, trainFile);
	cout << "Predicting........" << endl;
	prediction(bn, testFile, trueFile);
	printAccuracy();

	return 0;
}