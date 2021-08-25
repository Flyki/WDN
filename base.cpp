#include "base.h"

base::base(char* fileName, string pipeFile, double pre) {
	//string tempstr = fileName;
	//tempstr = tempstr.substr(0, tempstr.find(".inp"));
	char oufile[] = "100.out";
	char refile[] = "100.rep";
	int errCode = ENopen(fileName, refile, oufile);
	if (errCode > 0) exit(0);

	readPipeOptions(pipeFile);
	getAllJunctions();
	getAllRegularPipes();

	double* maxSolution = new double[dim];
	for (int i = 0; i < dim; i++) {
		maxSolution[i] = (double)pipeCoefficient.size() - 1;
	}
	this->maxPrice = getOneSolutionPrice(maxSolution);
//	this->tocheckfeasible = getOneFitness(maxSolution);
	this->PRESSURELIMIT = pre;
}

void base::getAllJunctions() {
	int temp;
	ENgetcount(EN_NODECOUNT, &temp);
	for (int i = 1; i <= temp; i++) {
		int type;
		ENgetnodetype(i, &type);
		if (type == 0) {
			float demand;
			ENgetnodevalue(i, EN_BASEDEMAND, &demand);
			if (demand != 0)
				junctions.push_back(i);
		}
		else if (type == 1) {
			sources.push_back(i);
		}
	}
}

void base::getAllRegularPipes() {
	int temp;
	ENgetcount(EN_LINKCOUNT, &temp);
	for (int i = 1; i <= temp; i++) {
		int type;
		ENgetlinktype(i, &type);
		if (type == 1) {
			pipes.push_back(i);
			float len;
			ENgetlinkvalue(i, EN_LENGTH, &len);
			pipeLengths.push_back(len);
		}
	}
	this->dim = pipes.size();
}

void base::interpret(double* oneSolution) {
	for (int i = 0; i < dim; i++) {
		ENsetlinkvalue(pipes[i], EN_DIAMETER, pipeDiameter[(int)oneSolution[i]]);
		ENsetlinkvalue(pipes[i], EN_ROUGHNESS, pipeCoefficient[(int)oneSolution[i]]);
	}
}

void base::interpret(int* oneSolution) {
	for (int i = 0; i < dim; i++) {
		ENsetlinkvalue(pipes[i], EN_DIAMETER, pipeDiameter[oneSolution[i]]);
		ENsetlinkvalue(pipes[i], EN_ROUGHNESS, pipeCoefficient[oneSolution[i]]);
	}
}

void base::interpret2(double* oneSolution) {
	for (int i = 0; i < dim; i++) {
		if (oneSolution[i] == -1) {
			ENsetlinkvalue(pipes[i], EN_DIAMETER, 0.00001);
			//ENsetlinkvalue(pipes[i], EN_ROUGHNESS, 1);
		}
		else {
			ENsetlinkvalue(pipes[i], EN_DIAMETER, pipeDiameter[(int)oneSolution[i]]);
			ENsetlinkvalue(pipes[i], EN_ROUGHNESS, pipeCoefficient[(int)oneSolution[i]]);
		}
	}
}

double base::getOneSolutionPrice(double* oneSolution) {
	double totalPrice = 0;
	for (int i = 0; i < dim; i++) {
		totalPrice += pipePrice[(int)oneSolution[i]] * pipeLengths[i];
	}
	return totalPrice;
}

double base::getOneSolutionPrice(int* oneSolution) {
	double totalPrice = 0;
	for (int i = 0; i < dim; i++) {
		totalPrice += pipePrice[oneSolution[i]] * pipeLengths[i];
	}
	return totalPrice;
}

double base::getOneSolutionPrice2(double* oneSolution) {
	double totalPrice = 0;
	for (int i = 0; i < dim; i++) {
		if (oneSolution[i] != -1)
			totalPrice += pipePrice[(int)oneSolution[i]] * pipeLengths[i];
	}
	return totalPrice;
}

pair<int, double> base::getViolation(int* oneSolution) {
	//set<int> individual;
	long t, tstep;
	int times = 0;
	double diff = 0;
	float p;
	ENopenH();
	interpret(oneSolution);
	//ENsettimeparam(EN_DURATION, 86400);
	ENinitH(0);
	do {
		ENrunH(&t);
		for (int i = 0; i < (int)junctions.size(); i++) {
			ENgetnodevalue(junctions[i], EN_PRESSURE, &p);
			if (p < PRESSURELIMIT) {
				//individual.insert(junctions[i]);
				times++;
				diff += (PRESSURELIMIT - p);
			}
		}
		ENnextH(&tstep);
	} while (tstep > 0);
	ENcloseH();
	return make_pair(times, diff);
}

pair<int, double> base::getViolation(double* oneSolution) {
	//set<int> individual;
	long t, tstep;
	int times = 0;
	double diff = 0;
	float p;
	ENopenH();
	interpret(oneSolution);
	//ENsettimeparam(EN_DURATION, 86400);
	ENinitH(0);
	do {
		ENrunH(&t);
		for (int i = 0; i < (int)junctions.size(); i++) {
			ENgetnodevalue(junctions[i], EN_PRESSURE, &p);
			if (p < PRESSURELIMIT) {
				//individual.insert(junctions[i]);
				times++;
				diff += (PRESSURELIMIT - p);
			}
		}
		ENnextH(&tstep);
	} while (tstep > 0);
	ENcloseH();
	return make_pair(times, diff);
}

pair<int, double> base::getViolation2(double* oneSolution) {
	//set<int> individual;
	long t, tstep;
	int times = 0;
	double diff = 0;
	float p;
	ENopenH();
	interpret2(oneSolution);
	//ENsettimeparam(EN_DURATION, 86400);
	ENinitH(0);
	do {
		ENrunH(&t);
		for (int i = 0; i < (int)junctions.size(); i++) {
			ENgetnodevalue(junctions[i], EN_PRESSURE, &p);
			if (p < PRESSURELIMIT) {
				//individual.insert(junctions[i]);
				times++;
				diff += (PRESSURELIMIT - p);
			}
		}
		ENnextH(&tstep);
	} while (tstep > 0);
	ENcloseH();
	return make_pair(times, diff);
}

double base::getOneFitness(double* oneSolution) {
	double onePrice = getOneSolutionPrice(oneSolution);
	pair<int, double> violation = getViolation(oneSolution);
	onePrice = onePrice / maxPrice;
	return violation.second + onePrice + violation.first;
}

double base::getOneFitness(int* oneSolution) {
	double onePrice = getOneSolutionPrice(oneSolution);
	pair<int, double> violation = getViolation(oneSolution);
	onePrice = onePrice / maxPrice;
	return violation.second + onePrice + violation.first;
}

double base::getOneFitness2(double* oneSolution) {
	double onePrice = getOneSolutionPrice2(oneSolution);
	pair<int, double> violation = getViolation2(oneSolution);
	onePrice = onePrice / maxPrice;
	return violation.second + onePrice + violation.first;
}

int* base::partitionBySourceForMultiple2(double* accordingSolution) {
	interpret(accordingSolution);
	ENsettimeparam(EN_DURATION, 604800);
	double** percentages = new double*[junctions.size()];
	int maxNodeIndex = 0;
	for (int i = 0; i < junctions.size(); i++) {
		if (maxNodeIndex < junctions[i])
			maxNodeIndex = junctions[i];
	}
	int* nodeBelong = new int[(maxNodeIndex + 1 + sources.size()) * 2];
	for (int i = 0; i < (int)junctions.size(); i++) {
		percentages[i] = new double[sources.size()];
		memset(percentages[i], 0, sizeof(double) * sources.size());
	}
	vector<vector<double>> actualDemand;
	for (int i = 0; i < junctions.size(); i++) {
		vector<double> temp;
		actualDemand.push_back(temp);
	}
	double* pipeDirection = new double[pipes.size()];
	memset(pipeDirection, 0, sizeof(double) * pipes.size());
	long t, tstep;
	ENopenH();
	ENinitH(0);
	do
	{
		ENrunH(&t);
		if (t % 3600 == 0) {
			for (int i = 0; i < (int)junctions.size(); i++) {
				float acd;
				ENgetnodevalue(junctions[i], EN_DEMAND, &acd);
				actualDemand[i].push_back(acd);
			}
			for (int i = 0; i < (int)pipes.size(); i++) {
				float fr;
				ENgetlinkvalue(pipes[i], EN_FLOW, &fr);
				pipeDirection[i] += fr;
			}
		}
		ENnextH(&tstep);
	} while (tstep > 0);
	ENcloseH();
	ENsolveH();
	for (int i = 0; i < (int)sources.size(); i++) {
		char oneId[50];
		ENgetnodeid(sources[i], oneId);
		ENsetqualtype(EN_TRACE, oneId, oneId, oneId);
		ENopenQ();
		ENinitQ(0);
		int counter = 0;
		do
		{
			ENrunQ(&t);
			if (t % 3600 == 0) {
				for (int j = 0; j < junctions.size(); j++) {
					float ra;
					ENgetnodevalue(junctions[j], EN_QUALITY, &ra);
					percentages[j][i] += (ra * actualDemand[j][counter] / 100.0);
				}
				counter++;
			}
			ENnextQ(&tstep);
		} while (tstep > 0);
		ENcloseQ();
	}
	memset(nodeBelong, 0, sizeof(int) * ((maxNodeIndex + 1 + sources.size()) * 2));
	for (int i = 0; i < junctions.size(); i++) {
		if (percentages[i][0] > percentages[i][1]) {
			nodeBelong[junctions[i] * 2] = 0;
			if (percentages[i][1] != 0) {
				nodeBelong[junctions[i] * 2 + 1] = 1;
			}
			else {
				nodeBelong[junctions[i] * 2 + 1] = -1;
			}
		}
		else {
			nodeBelong[junctions[i] * 2] = 1;
			if (percentages[i][0] != 0) {
				nodeBelong[junctions[i] * 2 + 1] = 0;
			}
			else {
				nodeBelong[junctions[i] * 2 + 1] = -1;
			}
		}
	}
	for (int i = 0; i < sources.size(); i++) {
		nodeBelong[sources[i] * 2] = i;
	}
	int* linkBelong = new int[pipes.size()];
	ofstream partitionfile("pipePartition.txt");
	for (int i = 0; i < (int)pipes.size(); i++) {
		int start, end;
		ENgetlinknodes(pipes[i], &start, &end);
		if (pipeDirection[i] < 0) {
			linkBelong[i] = nodeBelong[start * 2];
			partitionfile << linkBelong[i];
			if (nodeBelong[start * 2 + 1] != -1) {
				partitionfile << " " << nodeBelong[start * 2 + 1] << endl;
			}
			else {
				partitionfile << endl;
			}
		}
		else {
			linkBelong[i] = nodeBelong[end * 2];
			partitionfile << linkBelong[i];
			if (nodeBelong[end * 2 + 1] != -1) {
				partitionfile << " " << nodeBelong[end * 2 + 1] << endl;
			}
			else {
				partitionfile << endl;
			}
		}
	}
	for (int i = 0; i < (int)junctions.size(); i++) {
		delete[] percentages[i];
	}
	delete[] percentages;
	delete[] nodeBelong;
	delete[] pipeDirection;
	ENsettimeparam(EN_DURATION, 21600);
	partitionfile.close();
	return linkBelong;
}

int* base::partitionBySourceForMultiple(double* accordingSolution) {
	interpret(accordingSolution);
	ENsettimeparam(EN_DURATION, 604800);
	double** percentages = new double*[junctions.size()];
	int maxNodeIndex = 0;
	for (int i = 0; i < junctions.size(); i++) {
		if (maxNodeIndex < junctions[i])
			maxNodeIndex = junctions[i];
	}
	int* nodeBelong = new int[maxNodeIndex + 1 + sources.size()];
	for (int i = 0; i < (int)junctions.size(); i++) {
		percentages[i] = new double[sources.size()];
		memset(percentages[i], 0, sizeof(double) * sources.size());
	}
	vector<vector<double>> actualDemand;
	for (int i = 0; i < junctions.size(); i++) {
		vector<double> temp;
		actualDemand.push_back(temp);
	}
	double* pipeDirection = new double[pipes.size()];
	memset(pipeDirection, 0, sizeof(double) * pipes.size());
	long t, tstep;
	ENopenH();
	ENinitH(0);
	do
	{
		ENrunH(&t);
		if (t % 3600 == 0) {
			for (int i = 0; i < (int)junctions.size(); i++) {
				float acd;
				ENgetnodevalue(junctions[i], EN_DEMAND, &acd);
				actualDemand[i].push_back(acd);
			}
			for (int i = 0; i < (int)pipes.size(); i++) {
				float fr;
				ENgetlinkvalue(pipes[i], EN_FLOW, &fr);
				pipeDirection[i] += fr;
			}
		}
		ENnextH(&tstep);
	} while (tstep > 0);
	ENcloseH();
	ENsolveH();
	for (int i = 0; i < (int)sources.size(); i++) {
		char oneId[50];
		ENgetnodeid(sources[i], oneId);
		ENsetqualtype(EN_TRACE, oneId, oneId, oneId);
		ENopenQ();
		ENinitQ(0);
		int counter = 0;
		do
		{
			ENrunQ(&t);
			if (t % 3600 == 0) {
				for (int j = 0; j < junctions.size(); j++) {
					float ra;
					ENgetnodevalue(junctions[j], EN_QUALITY, &ra);
					percentages[j][i] += (ra * actualDemand[j][counter] / 100.0);
				}
				counter++;
			}
			ENnextQ(&tstep);
		} while (tstep > 0);
		ENcloseQ();
	}
	memset(nodeBelong, 0, sizeof(int) * (maxNodeIndex + 1));
	for (int i = 0; i < junctions.size(); i++) {
		for (int j = 0; j < sources.size(); j++) {
			if (percentages[i][j] > percentages[i][nodeBelong[junctions[i]]]) {
				nodeBelong[junctions[i]] = j;
			}
		}
	}
	for (int i = 0; i < (int)sources.size(); i++) {
		nodeBelong[sources[i]] = i;
	}
	int* linkBelong = new int[pipes.size()];
	for (int i = 0; i < (int)pipes.size(); i++) {
		int start, end;
		ENgetlinknodes(pipes[i], &start, &end);
		if (pipeDirection[i] < 0) {
			linkBelong[i] = nodeBelong[start];
		}
		else {
			linkBelong[i] = nodeBelong[end];
		}
	}
	for (int i = 0; i < (int)junctions.size(); i++) {
		delete[] percentages[i];
	}
	delete[] percentages;
	delete[] nodeBelong;
	delete[] pipeDirection;
	ENsettimeparam(EN_DURATION, 21600);
	return linkBelong;
}

int* base::partitionBySourceForSingle(double* accordingSolution) {
	interpret(accordingSolution);
	ENsettimeparam(EN_DURATION, 604800);
	double** percentages = new double*[junctions.size()];
	int maxNodeIndex = 0;
	for (int i = 0; i < (int)junctions.size(); i++) {
		if (maxNodeIndex < junctions[i])
			maxNodeIndex = junctions[i];
	}
	int* nodeBelong = new int[maxNodeIndex + 1 + sources.size()];
	for (int i = 0; i < (int)junctions.size(); i++) {
		percentages[i] = new double[sources.size()];
		memset(percentages[i], 0, sizeof(double) * sources.size());
	}
	double* pipeDirection = new double[pipes.size()];
	memset(pipeDirection, 0, sizeof(double) * pipes.size());
	long t, tstep;
	ENopenH();
	ENinitH(0);
	do
	{
		ENrunH(&t);
		for (int i = 0; i < (int)pipes.size(); i++) {
			float fr;
			ENgetlinkvalue(pipes[i], EN_FLOW, &fr);
			pipeDirection[i] += fr;
		}
		ENnextH(&tstep);
	} while (tstep > 0);
	ENcloseH();
	ENsolveH();
	for (int i = 0; i < (int)sources.size(); i++) {
		char oneId[50];
		ENgetnodeid(sources[i], oneId);
		ENsetqualtype(EN_TRACE, oneId, oneId, oneId);
		ENopenQ();
		ENinitQ(0);
		do
		{
			ENrunQ(&t);
			for (int j = 0; j < (int)junctions.size(); j++) {
				float ra;
				ENgetnodevalue(junctions[j], EN_QUALITY, &ra);
				percentages[j][i] += ra;
			}
			ENnextQ(&tstep);
		} while (tstep > 0);
		ENcloseQ();
	}
	memset(nodeBelong, 0, sizeof(int) * (maxNodeIndex + 1));
	for (int i = 0; i < (int)junctions.size(); i++) {
		for (int j = 0; j < (int)sources.size(); j++) {
			if (percentages[i][j] > percentages[i][nodeBelong[junctions[i]]]) {
				nodeBelong[junctions[i]] = j;
			}
		}
	}
	for (int i = 0; i < (int)sources.size(); i++) {
		nodeBelong[sources[i]] = i;
	}
	int* linkBelong = new int[pipes.size()];
	for (int i = 0; i < (int)pipes.size(); i++) {
		int start, end;
		ENgetlinknodes(pipes[i], &start, &end);
		if (pipeDirection[i] < 0) {
			linkBelong[i] = nodeBelong[start];
		}
		else {
			linkBelong[i] = nodeBelong[end];
		}
	}
	for (int i = 0; i < (int)junctions.size(); i++) {
		delete[] percentages[i];
	}
	delete[] percentages;
	delete[] nodeBelong;
	delete[] pipeDirection;
	ENsettimeparam(EN_DURATION, 0);
	return linkBelong;
}

void base::readPipeOptions(string fileName) {
	std::ifstream pipeFile(fileName.c_str());
	int option = 0;
	pipeFile >> option;
	for (int i = 0; i < option; i++) {
		double temp;
		pipeFile >> temp;
		pipeDiameter.push_back(temp);
		pipeFile >> temp;
		pipePrice.push_back(temp);
		pipeFile >> temp;
		pipeCoefficient.push_back(temp);
	}
}