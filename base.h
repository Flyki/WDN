#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdio>
#include <set>
#include <cstring>
#include <cfloat>
extern "C" {
	#include "toolkit.h"
}

using namespace std;

class base {
public:
	base(char*, string, double);
	void getAllJunctions();
	void getAllRegularPipes();
	int* partitionBySourceForSingle(double*);
	int* partitionBySourceForMultiple(double*);
	int* partitionBySourceForMultiple2(double*);
	void interpret(double*);
	void interpret(int*);
	void interpret2(double*);
	double getOneSolutionPrice(double*);
	double getOneSolutionPrice(int*);
	double getOneSolutionPrice2(double*);
	pair<int, double> getViolation(double*);
	pair<int, double> getViolation(int*);
	pair<int, double> getViolation2(double*);
	double getOneFitness(double*);
	double getOneFitness(int*);
	double getOneFitness2(double*);
	void readPipeOptions(string);

	vector<int> pipes;
	vector<double> pipeLengths;
	vector<int> junctions;
	vector<int> sources;
	int dim;
//	double tocheckfeasible;
	double maxPrice;
	double PRESSURELIMIT;

	std::vector<double> pipeDiameter;
	std::vector<double> pipePrice;
	std::vector<double> pipeCoefficient;
};