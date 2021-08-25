#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <string>
#include <cstdio>
#include <set>
#include <random>
#include <chrono>
#include <sstream>
#include <queue>
#include <cfloat>
#include "base.h"

using namespace std;

class TSMA : public base {
public:
    TSMA(int, int, int, int, int, char* infile, string pipefile, double thre);
    ~TSMA();
    void run();
    void initialize1(); //for whole search space
    void initialize2(); //for promising regione
    void funcbody1();   //for ILLSO
    void localSearchExpWidth(double*, double&);     //1
    void localSearchExpDepth(double*, double&);     //2
    double largest_difference();

    int ssize;
    int fes;
    double** swarm;
    double** velocity;
    double* fit;
    double* gbestx;
    double gbestf;
    int* levelN;
    double* levelP;
    int usedFes;
    int curbest;
    int flag;
    double afy;
    int loc;
    ofstream result;
    default_random_engine gen;
    uniform_real_distribution<double> dis;
};
