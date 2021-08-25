#include "TSMA.h"

TSMA::TSMA(int ssize, int fes, int seed, int ID, int loc, char* infile, string pipefile, double thre) : base(infile, pipefile, thre) {
    this->ssize = ssize;
	this->fes = fes;
	this->swarm = new double* [ssize];
	this->velocity = new double* [ssize];
	this->fit = new double[ssize];
	this->gbestx = new double[dim];
	for (int i = 0; i < ssize; i++) {
		this->swarm[i] = new double[dim];
		this->velocity[i] = new double[dim];
	}
	this->gbestf = DBL_MAX;
	levelN = new int[6];
	levelN[0] = 4; levelN[1] = 6; levelN[2] = 8;
	levelN[3] = 10; levelN[4] = 20; levelN[5] = 50;
	levelP = new double[6];
	for (int i = 0; i < 6; i++) {
		levelP[i] = 1;
	}

	default_random_engine tempGen(seed);
	this->gen = tempGen;
	uniform_real_distribution<double> tempDis(0.0, 1.0);
	this->dis = tempDis;
	this->usedFes = 0;
    this->curbest = 0;
    this->flag = 0;
    this->afy = (((double)ssize / 200.0) - 1.0) * 0.05;
	this->loc = loc;

    stringstream ss;
	ss << "TSMA." << ID << "." << seed << "." << loc << ".result.txt";
	string rfname;
	ss >> rfname;
	ss.clear();
	result.open(rfname.c_str(), ios::app);
	result << "start" << endl;
}

TSMA::~TSMA() {
    for (int i = 0; i < ssize; i++) {
		delete[] swarm[i];
		delete[] velocity[i];
	}
	delete[] swarm;
	delete[] velocity;
	delete[] fit;
	delete[] gbestx;
	delete[] levelN;
	delete[] levelP;
	result << "end" << endl;
	result.close();
}

void TSMA::run() {
    initialize1();
	int counter = 0;
	list<double> ldlist;
	double lastmean = DBL_MAX;
    while (usedFes < fes)
    {
		double oldcurf = fit[curbest];
        funcbody1()

		if (oldcurf == fit[curbest]) {
			counter++;
		}
		else
			counter = 0;

		double dista = largest_difference();
		if (flag == 0) {
			ldlist.push_back(dista);
		}
		if (ldlist.size() == 30 && flag == 0) {
			double curmean = 0;
			for (list<double>::iterator itr = ldlist.begin(); itr != ldlist.end(); ++itr) {
				curmean += *itr;
			}
			curmean = curmean / 30.0;
			if (curmean < 5 && curmean >= lastmean) {
				if (loc == 1) {
					localSearchExpWidth(gbestx, gbestf);
				}
				else if (loc == 2) {
					localSearchExpDepth(gbestx, gbestf);
				}
				initialize2();
				flag = 1;
			}
			lastmean = curmean;
			ldlist.clear();
		}
		result << usedFes << ',' << gbestf << ',' << fit[curbest] << ',' << counter << ',' << dista << endl;
    }
	if (loc == 1) {
		localSearchExpWidth(swarm[curbest], fit[curbest]);
	}
	else if (loc == 2) {
		localSearchExpDepth(swarm[curbest], fit[curbest]);
	}
	result << usedFes << ',' << gbestf << endl;
}

void TSMA::initialize1() {
	curbest = 0;
	for (int i = 0; i < ssize; i++) {
		for (int j = 0; j < dim; j++) {
			swarm[i][j] = dis(gen) * (double)pipeCoefficient.size();
			velocity[i][j] = 0;
		}
		fit[i] = getOneFitness(swarm[i]);
		usedFes++;
		if (fit[i] < fit[curbest]) curbest = i;
	}
	memcpy(gbestx, swarm[curbest], sizeof(double) * dim);
	gbestf = fit[curbest];
}

void TSMA::initialize2() {
	double deviation = ((double)pipeCoefficient.size()) / 8.0;
	if (deviation < 2)
		deviation = 2;
	for (int j = 0; j < dim; j++) {
		double upper = 0;
		double lower = pipeCoefficient.size();
		for (int i = 0; i < ssize; i++) {
			if (upper < swarm[i][j]) {
				upper = swarm[i][j];
			}
			if (lower > swarm[i][j]) {
				lower = swarm[i][j];
			}
		}
		if (upper < gbestx[j] + deviation) {
			upper = gbestx[j] + deviation;
		}
		if (lower > gbestx[j] - deviation) {
			lower = gbestx[j] - deviation;
		}
		if (upper > pipeCoefficient.size())
			upper = pipeCoefficient.size();
		if (lower < 0)
			lower = 0;
		uniform_real_distribution<double> ndis(lower, upper);
		for (int i = 0; i < ssize; i++) {
			swarm[i][j] = ndis(gen);
		}
	}
	curbest = 0;
	for (int i = 0; i < ssize; i++) {
		fit[i] = getOneFitness(swarm[i]);
		usedFes++;
		if (fit[i] < fit[curbest]) curbest = i;
	}
	for (int i = 0; i < ssize; i++) {
		memset(velocity[i], 0, sizeof(double) * dim);
	}
	for (int i = 0; i < 6; i++) {
		levelP[i] = 1;
	}
	if (fit[curbest] < gbestf) {
		gbestf = fit[curbest];
		memcpy(gbestx, swarm[curbest], sizeof(double) * dim);
	}
}

void TSMA::funcbody1() {
	int* order = new int[ssize];
	for (int i = 0; i < ssize; i++) {
		order[i] = i;
	}
	for (int i = 0; i < ssize - 1; i++) {
		for (int j = 0; j < ssize - 1 - i; j++) {
			if (fit[order[j]] > fit[order[j + 1]]) {
				int temp = order[j];
				order[j] = order[j + 1];
				order[j + 1] = temp;
			}
		}
	}
	double choice[6];
	for (int i = 0; i < 6; i++)
	{
		choice[i] = std::exp(levelP[i] * 7);
	}
	for (int i = 1; i < 6; i++) {
		choice[i] += choice[i - 1];
	}
	for (int i = 0; i < 6; i++) {
		choice[i] = choice[i] / choice[5];
	}
	double thresh = dis(gen);
	int theOne = 0;
	for (int i = 0; i < 6; i++) {
		if (choice[i] >= thresh) {
			theOne = i;
			break;
		}
	}
	int levelNumber = levelN[theOne];
	int levelEach = ssize / levelNumber;

	vector<vector<double>> meanvs;
	for (int i = 0; i < levelNumber; i++) {
		vector<double> onemeav(dim, 0);
		int start = i * levelEach;
		int end = (i + 1) * levelEach;
		if (i == levelNumber - 1) {
			end = ssize;
		}
		for (int j = start; j < end; j++) {
			for (int k = 0; k < dim; k++) {
				onemeav[k] += swarm[order[j]][k];
			}
		}
		for (int k = 0; k < dim; k++) {
			onemeav[k] = onemeav[k] / (end - start);
		}
		meanvs.push_back(onemeav);
	}

	for (int i = (levelNumber - 1) * levelEach; i < ssize; i++) {
		int firstLevel, secondLevel;
		firstLevel = dis(gen) * (levelNumber - 1);
		do
		{
			secondLevel = dis(gen) * (levelNumber - 1);
		} while (secondLevel == firstLevel);
		if (firstLevel > secondLevel) {
			int temp = firstLevel;
			firstLevel = secondLevel;
			secondLevel = temp;
		}
		int firstOne = dis(gen) * levelEach + firstLevel * levelEach;
		int curOne = order[i];
		firstOne = order[firstOne];
		for (int j = 0; j < dim; j++) {
			velocity[curOne][j] = dis(gen) * velocity[curOne][j] + dis(gen) * (swarm[firstOne][j] - swarm[curOne][j]) + dis(gen) * afy * (meanvs[secondLevel][j] - swarm[curOne][j]);
			swarm[curOne][j] += velocity[curOne][j];
			if (swarm[curOne][j] >= pipeCoefficient.size()) swarm[curOne][j] = pipeCoefficient.size() - 1;
			if (swarm[curOne][j] < 0) swarm[curOne][j] = 0;
		}
	}
	for (int i = levelNumber - 2; i >= 2; i--) {
		for (int j = i * levelEach; j < (i + 1) * levelEach; j++) {
			int firstLevel, secondLevel;
			firstLevel = dis(gen) * i;
			do
			{
				secondLevel = dis(gen) * i;
			} while (firstLevel == secondLevel);
			if (firstLevel > secondLevel) {
				int temp = firstLevel;
				firstLevel = secondLevel;
				secondLevel = temp;
			}
			int firstOne = dis(gen) * levelEach + firstLevel * levelEach;
			int curOne = order[j];
			firstOne = order[firstOne];
			for (int k = 0; k < dim; k++) {
				velocity[curOne][k] = dis(gen) * velocity[curOne][k] + dis(gen) * (swarm[firstOne][k] - swarm[curOne][k]) + dis(gen) * afy * (meanvs[secondLevel][k] - swarm[curOne][k]);
				swarm[curOne][k] += velocity[curOne][k];
				if (swarm[curOne][k] >= pipeCoefficient.size()) swarm[curOne][k] = pipeCoefficient.size() - 1;
				if (swarm[curOne][k] < 0) swarm[curOne][k] = 0;
			}
		}
	}
	for (int i = levelEach; i < 2 * levelEach; i++) {
		int firstOne = dis(gen) * levelEach;
		int curOne = order[i];
		firstOne = order[firstOne];
		for (int k = 0; k < dim; k++) {
			velocity[curOne][k] = dis(gen) * velocity[curOne][k] + dis(gen) * (swarm[firstOne][k] - swarm[curOne][k]) + dis(gen) * afy * (meanvs[0][k] - swarm[curOne][k]);
			swarm[curOne][k] += velocity[curOne][k];
			if (swarm[curOne][k] >= pipeCoefficient.size()) swarm[curOne][k] = pipeCoefficient.size() - 1;
			if (swarm[curOne][k] < 0) swarm[curOne][k] = 0;
		}
	}
	int thBest = order[levelEach];
	for (int i = levelEach; i < ssize; i++) {
		fit[order[i]] = getOneFitness(swarm[order[i]]);
		usedFes++;
		if (fit[order[i]] < fit[order[thBest]]) thBest = i;
	}
    double lastf = gbestf;
	if (fit[order[thBest]] < gbestf) {
		gbestf = fit[order[thBest]];
		memcpy(gbestx, swarm[order[thBest]], sizeof(double) * dim);
	}
	levelP[theOne] = abs(gbestf - lastf) / abs(lastf);
	delete[] order;
    for (int i = 0; i < ssize; i++) {
		if (fit[curbest] > fit[i])
			curbest = i;
	}
}

void TSMA::localSearchExpWidth(double* currentSolution, double& currentFit) {
	if (currentFit <= 1) {
		double* tempb = new double[dim];
		double tempf = currentFit;
		memcpy(tempb, currentSolution, sizeof(double) * dim);
		vector<pair<double, int>> still;
		for (int i = 0; i < dim; i++) {
			if (((int)tempb[i]) > 0) {
				double pr = pipePrice[(int)tempb[i]] - pipePrice[((int)tempb[i]) - 1];
				pr = pr * pipeLengths[i];
				still.push_back(make_pair(pr, i));
			}
		}
		sort(still.begin(), still.end());
		reverse(still.begin(), still.end());
		while (!still.empty())
		{
			vector<pair<double, int>> still2;
			for (auto& e : still) {
				int p = e.second;
				tempb[p] = tempb[p] - 1;
				double tempScore = getOneFitness(tempb);
				usedFes++;
				if (tempScore <= 1) {
					tempf = tempScore;
					if (((int)tempb[p]) > 0) {
						double pr = pipePrice[(int)tempb[p]] - pipePrice[(int)tempb[p] - 1];
						pr = pr * pipeLengths[p];
						still2.push_back(make_pair(pr, p));
					}
				}
				else {
					tempb[p] = tempb[p] + 1;
				}
			}
			still = still2;
			sort(still.begin(), still.end());
			reverse(still.begin(), still.end());
		}
		currentFit = tempf;
		memcpy(currentSolution, tempb, sizeof(double) * dim);
		if (currentFit <= gbestf) {
			gbestf = currentFit;
			memcpy(gbestx, currentSolution, sizeof(double) * dim);
		}
		delete[] tempb;
	}
}

void TSMA::localSearchExpDepth(double* currentSolution, double& currentFit) {
	if (currentFit <= 1) {
		int* tempb = new int[dim];
		double tempf = currentFit;
		for (int i = 0; i < dim; i++) {
			tempb[i] = (int)currentSolution[i];
		}
		priority_queue<pair<double, int>> still;
		for (int i = 0; i < dim; i++) {
			if (tempb[i] > 0) {
				double pr = pipePrice[tempb[i]] - pipePrice[currentSolution[i] - 1];
				pr = pr * pipeLengths[i];
				still.push(make_pair(pr, i));
			}
		}
		while (!still.empty())
		{
			int p = still.top().second;
			still.pop();
			tempb[p] = tempb[p] - 1;
			double tempScore = getOneFitness(tempb);
			usedFes++;
			if (tempScore <= 1) {
				tempf = tempScore;
				if (tempb[p] > 0) {
					double pr = pipePrice[tempb[p]] - pipePrice[tempb[p] - 1];
					pr = pr * pipeLengths[p];
					still.push(make_pair(pr, p));
				}
			}
			else {
				tempb[p] = tempb[p] + 1;
			}
		}
		currentFit = tempf;
		for (int i = 0; i < dim; i++) {
			currentSolution[i] = tempb[i] + 0.5;
		}
		if (currentFit <= gbestf) {
			gbestf = currentFit;
			memcpy(gbestx, currentSolution, sizeof(double) * dim);
		}
		delete[] tempb;
	}
}

double TSMA::largest_difference() {
	double largestdif = 0;
	for (int i = 0; i < dim; i++) {
		double largestone = -1;
		double smalestone = DBL_MAX;
		for (int j = 0; j < ssize; j++) {
			if (swarm[j][i] > largestone)
				largestone = swarm[j][i];
			if (swarm[j][i] < smalestone)
				smalestone = swarm[j][i];
		}
		if (largestdif < largestone - smalestone)
			largestdif = largestone - smalestone;
	}
	return largestdif;
}
