#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdio>
#include <set>
#include <random>
#include <chrono>
#include <sstream>
#include "base.h"
using namespace std;

struct NSFC
{
	int ns1;
	int ns2;
	int nf1;
	int nf2;
	double wcr;
	double deltaf;
	int best;
	NSFC() {
		ns1 = ns2 = nf1 = nf2 = wcr = deltaf = best = 0;
	}
};

class ccdeRepart : public base {
public:
	ccdeRepart(int, int, int, int, int, char*, string, double);
	~ccdeRepart();
	void run();
	void initialize(int);
	int recalculate(int);
	NSFC funcbody(int, int);

	int ssize;
	int fes;
	double** swarm1;
	double** swarm2;
	double** fitness;
	double* ccbest;
	double bestf;
	double** crv;
	double* crm;
	double* pf;
	vector<vector<int>> groupInfo;
	double* ns1;
	double* ns2;
	double* nf1;
	double* nf2;
	double* wcr;
	double* deltaf;
	int usedFes;
	int repInterval;

	ofstream result;
	default_random_engine gen;
	uniform_real_distribution<double> dis;
	normal_distribution<double> ndis;
	cauchy_distribution<double> cdis;
};