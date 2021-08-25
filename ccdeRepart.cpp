#include "ccdeRepart.h"

ccdeRepart::ccdeRepart(int ssize, int fes, int repInterval, int seed, int ID, char* fileName, string pipeFile, double pre) : base(fileName, pipeFile, pre) {
	this->ssize = ssize;
	this->fes = fes;
	this->usedFes = 0;
	this->repInterval = repInterval;
	this->swarm1 = new double*[this->ssize];
	this->swarm2 = new double*[this->ssize];
	for (int i = 0; i < this->ssize; i++) {
		this->swarm1[i] = new double[dim];
		this->swarm2[i] = new double[dim];
	}
	this->ccbest = new double[dim];
	this->fitness = new double*[(int)sources.size()];
	for (int i = 0; i < (int)sources.size(); i++) {
		this->fitness[i] = new double[this->ssize];
	}
	this->crv = new double*[sources.size()];
	this->crm = new double[sources.size()];
	this->pf = new double[sources.size()];
	for (int i = 0; i < sources.size(); i++) {
		crv[i] = new double[this->ssize];
		for (int j = 0; j < this->ssize; j++) {
			crv[i][j] = 0.5;
		}
		crm[i] = 0.5;
		pf[i] = 0.5;
	}
	ns1 = new double[sources.size()];
	ns2 = new double[sources.size()];
	nf1 = new double[sources.size()];
	nf2 = new double[sources.size()];
	wcr = new double[sources.size()];
	deltaf = new double[sources.size()];
	memset(ns1, 0, sizeof(double) * sources.size());
	memset(ns2, 0, sizeof(double) * sources.size());
	memset(nf1, 0, sizeof(double) * sources.size());
	memset(nf2, 0, sizeof(double) * sources.size());
	memset(wcr, 0, sizeof(double) * sources.size());
	memset(deltaf, 0, sizeof(double) * sources.size());

	//unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine tempGen(seed);
	this->gen = tempGen;
	uniform_real_distribution<double> tempDis(0.0, 1.0);
	this->dis = tempDis;
	normal_distribution<double> tempNdis(0.5, 0.3);
	this->ndis = tempNdis;
	cauchy_distribution<double> tempCdis(0.0, 1.0);
	this->cdis = tempCdis;

	stringstream ss;
	ss << "CCDER." << ID << "." << seed << ".result.txt";
	string rfname;
	ss >> rfname;
	ss.clear();
	result.open(rfname.c_str(), ios::app);
	result << "start" << endl;
}

ccdeRepart::~ccdeRepart() {
	for (int i = 0; i < ssize; i++) {
		delete[] swarm1[i];
		delete[] swarm2[i];
	}
	delete[] swarm1;
	delete[] swarm2;
	delete[] ccbest;
	for (int i = 0; i < sources.size(); i++) {
		delete[] fitness[i];
		delete[] crv[i];
	}
	delete[] fitness;
	delete[] crv;
	delete[] crm;
	delete[] pf;
	delete[] ns1;
	delete[] ns2;
	delete[] nf1;
	delete[] nf2;
	delete[] wcr;
	delete[] deltaf;
	result << "end" << endl;
	result.close();
}

void ccdeRepart::run() {
	double* upbound = new double[this->dim];
	double* lastBest = new double[this->dim];
	int pipeType = pipeCoefficient.size() - 1;
	int former = pipeType;
	for (int i = pipeType; i >= 0; i--) {
		for (int j = 0; j < this->dim; j++) {
			upbound[j] = i;
		}
		pair<int, double> x = getViolation(upbound);
		usedFes++;
		if (x.first == 0) 
			former = i;
		else break;
	}
	for (int i = 0; i < this->dim; i++) {
		upbound[i] = former;
	}
	int* pipeBelong = partitionBySourceForSingle(upbound);
	usedFes++;
	for (int i = 0; i < (int)sources.size(); i++) {
		vector<int> temp;
		groupInfo.push_back(temp);
	}
	for (int i = 0; i < this->dim; i++) {
		groupInfo[pipeBelong[i]].push_back(i);
	}

	delete[] pipeBelong;
	initialize(former);
	memcpy(lastBest, ccbest, sizeof(double) * dim);
	usedFes++;
	int i = 0;
	while(usedFes < fes) {
		if ((i + 1) % this->repInterval == 0) {
			if (bestf <= 1) {
				int* pipeBelong = partitionBySourceForSingle(ccbest);
				usedFes++;
				for (int j = 0; j < (int)groupInfo.size(); j++) {
					groupInfo[j].clear();
				}
				for (int j = 0; j < this->dim; j++) {
					groupInfo[pipeBelong[j]].push_back(j);
				}
				delete[] pipeBelong;
			}
		}
		int rbest = 0;
		for (int j = 0; j < (int)groupInfo.size(); j++) {
			int curBest = 0;//recalculate(j);
			for (int k = 0; k < ssize; k++) {
				if (fitness[j][k] < fitness[j][curBest])
					curBest = k;
			}
			if ((i + 1) % 50 == 0) {
				pf[j] = ns1[j] * (ns2[j] + nf2[j]) / (ns2[j] * (ns1[j] + nf1[j]) + ns1[j] * (ns2[j] + nf2[j]));
				ns1[j] = ns2[j] = nf1[j] = nf2[j] = 0;
			}
			if ((i + 1) % 25 == 0) {
				crm[j] = wcr[j] / deltaf[j];
				wcr[j] = deltaf[j] = 0;
			}
			if ((i + 1) % 5 == 0) {
				normal_distribution<double> crmdis2(crm[j], 0.1);
				for (int l = 0; l < ssize; l++) {
					crv[j][l] = crmdis2(gen);
				}
			}
			NSFC re = funcbody(j, curBest);
			ns1[j] += re.ns1;
			ns2[j] += re.ns2;
			nf1[j] += re.nf1;
			nf2[j] += re.nf2;
			wcr[j] += re.wcr;
			deltaf[j] += re.deltaf;
			for (int outer = 0; outer < ssize; outer++) {
				for (int inner = 0; inner < (int)groupInfo[j].size(); inner++) {
					swarm1[outer][groupInfo[j][inner]] = swarm2[outer][groupInfo[j][inner]];
				}
			}
			for (int k = 0; k < (int)groupInfo[j].size(); k++) {
				ccbest[groupInfo[j][k]] = swarm1[re.best][groupInfo[j][k]];
			}
			for (int k = 0; k < (int)groupInfo[j].size(); k++) {
				double temp = swarm1[re.best][groupInfo[j][k]];
				swarm1[re.best][groupInfo[j][k]] = swarm1[0][groupInfo[j][k]];
				swarm1[0][groupInfo[j][k]] = temp;
			}
			double temp = fitness[j][re.best];
			fitness[j][re.best] = fitness[j][0];
			fitness[j][0] = temp;
			//rbest = re.best;
		}
		//cout << i << endl;
		/*bool allsame = true;
		for (int k = 0; k < dim; k++) {
			if (lastBest[k] != ccbest[k]) {
				allsame = false;
				break;
			}
		}
		if (allsame && bestf <= 1.0) {
			for (int k = 0; k < dim; k++) {
				if (dis(gen) < 0.08) {
					ccbest[k] = ccbest[k] - 1;
				}
			}
			double tempf = getOneFitness(ccbest);
			if (tempf < bestf) {
				bestf = tempf;
				memcpy(lastBest, ccbest, sizeof(double) * dim);
			}
			else {
				memcpy(ccbest, lastBest, sizeof(double) * dim);
			}
		}
		else {*/
		bestf = fitness[(int)groupInfo.size() - 1][0];// getOneFitness(ccbest);
			/*memcpy(lastBest, ccbest, sizeof(double) * dim);
		}*/
		result << usedFes << "," << bestf << endl;
		i++;
	}
	delete[] upbound;
	delete[] lastBest;
}

void ccdeRepart::initialize(int maxPara) {
	for (int i = 0; i < ssize; i++) {
		for (int j = 0; j < dim; j++) {
			swarm1[i][j] = swarm2[i][j] = dis(gen) * (double)pipeCoefficient.size();
		}
	}
	/*for (int j = 0; j < dim; j++) {
		swarm2[ssize - 1][j] = swarm1[ssize - 1][j] = maxPara;
	}*/
	for (int i = 0; i < dim; i++) {
		ccbest[i] = dis(gen) * (double)pipeCoefficient.size();
	}
	double* forCal = new double[dim];
	for (int i = 0; i < (int)groupInfo.size(); i++) {
		memcpy(forCal, ccbest, sizeof(double) * dim);
		int itsBest = 0;
		for (int j = 0; j < ssize; j++) {
			for (int k = 0; k < (int)groupInfo[i].size(); k++) {
				forCal[groupInfo[i][k]] = swarm1[j][groupInfo[i][k]];
			}
			fitness[i][j] = getOneFitness(forCal);
			usedFes++;
			if (fitness[i][j] < fitness[i][itsBest])
				itsBest = j;
		}
		for (int j = 0; j < (int)groupInfo[i].size(); j++) {
			ccbest[groupInfo[i][j]] = swarm1[itsBest][groupInfo[i][j]];
		}
	}
	bestf = getOneFitness(ccbest);
	usedFes++;
	delete[] forCal;
}

int ccdeRepart::recalculate(int groupNo) {
	double* temp = new double[dim];
	int bes = 0;
	memcpy(temp, ccbest, sizeof(double) * dim);
	for (int i = 0; i < ssize; i++) {
		for (int j = 0; j < (int)groupInfo[groupNo].size(); j++) {
			temp[groupInfo[groupNo][j]] = swarm1[i][groupInfo[groupNo][j]];
		}
		fitness[groupNo][i] = getOneFitness(temp);
		usedFes++;
		if (fitness[groupNo][i] < fitness[groupNo][bes]) bes = i;
	}
	delete[] temp;
	return bes;
}

NSFC ccdeRepart::funcbody(int groupNo, int curBest) {
	double* temp = new double[dim];
	memcpy(temp, ccbest, sizeof(double) * dim);
	NSFC re;
	re.best = 0;
	double* trial = new double[(int)groupInfo[groupNo].size()];
	for (int i = 0; i < ssize; i++) {
		if (dis(gen) < pf[groupNo]) {
			int a, b, c;
			do
			{
				a = dis(gen) * ssize;
			} while (a == i);
			do
			{
				b = dis(gen) * ssize;
			} while (b == a || b == i);
			do
			{
				c = dis(gen) * ssize;
			} while (c == a || c == b || c == i);
			int j = dis(gen) * (int)groupInfo[groupNo].size();
			double F = 0;
			if (dis(gen) < pf[groupNo]) F = ndis(gen);
			else F = cdis(gen);
			for (int k = 1; k <= (int)groupInfo[groupNo].size(); k++) {
				if (dis(gen) < crv[groupNo][i] || k == (int)groupInfo[groupNo].size()) {
					trial[j] = swarm1[c][groupInfo[groupNo][j]] + F * (swarm1[a][groupInfo[groupNo][j]] - swarm1[b][groupInfo[groupNo][j]]);
					if (trial[j] >= pipeCoefficient.size()) trial[j] = pipeCoefficient.size() - 1;
					if (trial[j] < 0) trial[j] = 0;
				}
				else trial[j] = swarm1[i][groupInfo[groupNo][j]];
				j = (j + 1) % groupInfo[groupNo].size();
			}
			//add some local refine techniques
			/*bool allsame = true;
			for (int k = 0; k < (int)groupInfo[groupNo].size(); k++) {
				if (trial[k] != swarm1[i][groupInfo[groupNo][k]]) {
					allsame = false;
					break;
				}
			}
			if (allsame && fitness[groupNo][i] < 1) {
				for (int k = 0; k < (int)groupInfo[groupNo].size(); k++) {
					if (dis(gen) < 0.1) {
						trial[k] = trial[k] - 1;
						if (trial[k] < 0) trial[k] = 0;
					}
				}
			}*/

			double score = DBL_MAX;
			for (int j = 0; j < (int)groupInfo[groupNo].size(); j++) {
				temp[groupInfo[groupNo][j]] = trial[j];
			}
			score = getOneFitness(temp);
			usedFes++;
			if (score <= fitness[groupNo][i]) {
				for (int j = 0; j < (int)groupInfo[groupNo].size(); j++) {
					swarm2[i][groupInfo[groupNo][j]] = trial[j];
				}
				re.ns1++;
				double de = fitness[groupNo][i] - score;
				re.deltaf += de;
				re.wcr += (de * crv[groupNo][i]);
				fitness[groupNo][i] = score;
			}
			else {
				for (int j = 0; j < (int)groupInfo[groupNo].size(); j++) {
					swarm2[i][groupInfo[groupNo][j]] = swarm1[i][groupInfo[groupNo][j]];
				}
				re.nf1++;
			}
		}
		else {
			int a, b;
			do
			{
				a = dis(gen) * ssize;
			} while (a == i);
			do
			{
				b = dis(gen) * ssize;
			} while (b == a || b == i);
			int j = dis(gen) * (int)groupInfo[groupNo].size();
			double F = 0;
			if (dis(gen) < pf[groupNo]) F = ndis(gen);
			else F = cdis(gen);
			for (int k = 1; k <= (int)groupInfo[groupNo].size(); k++) {
				if (dis(gen) < crv[groupNo][i] || k == groupInfo[groupNo].size()) {
					trial[j] = swarm1[i][groupInfo[groupNo][j]] + F * (swarm1[a][groupInfo[groupNo][j]] - swarm1[b][groupInfo[groupNo][j]]) + F * (swarm1[curBest][groupInfo[groupNo][j]] - swarm1[i][groupInfo[groupNo][j]]);
					if (trial[j] >= pipeCoefficient.size()) trial[j] = pipeCoefficient.size() - 1;
					if (trial[j] < 0) trial[j] = 0;
				}
				else {
					trial[j] = swarm1[i][groupInfo[groupNo][j]];
				}
				j = (j + 1) % groupInfo[groupNo].size();
			}
			//add some local refine techniques
			/*bool allsame = true;
			for (int k = 0; k < (int)groupInfo[groupNo].size(); k++) {
				if (trial[k] != swarm1[i][groupInfo[groupNo][k]]) {
					allsame = false;
					break;
				}
			}
			if (allsame && fitness[groupNo][i] < 1) {
				for (int k = 0; k < (int)groupInfo[groupNo].size(); k++) {
					if (dis(gen) < 0.1) {
						trial[k] = trial[k] - 1;
						if (trial[k] < 0) trial[k] = 0;
					}
				}
			}*/

			double score = DBL_MAX;
			for (int j = 0; j < (int)groupInfo[groupNo].size(); j++) {
				temp[groupInfo[groupNo][j]] = trial[j];
			}
			score = getOneFitness(temp);
			usedFes++;
			if (score <= fitness[groupNo][i]) {
				for (int j = 0; j < groupInfo[groupNo].size(); j++) {
					swarm2[i][groupInfo[groupNo][j]] = trial[j];
				}
				re.ns2++;
				double de = fitness[groupNo][i] - score;
				re.deltaf += de;
				re.wcr += (de * crv[groupNo][i]);
				fitness[groupNo][i] = score;
			}
			else {
				for (int j = 0; j < (int)groupInfo[groupNo].size(); j++) {
					swarm2[i][groupInfo[groupNo][j]] = swarm1[i][groupInfo[groupNo][j]];
				}
				re.nf2++;
			}
		}
		if (fitness[groupNo][re.best] > fitness[groupNo][i]) re.best = i;
	}
	delete[] temp;
	delete[] trial;
	return re;
}