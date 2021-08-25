#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include "base.h"
#include "TSMA.h"
#include "ccdeRepart.h"


using namespace std;


int main(int argc, char* argv[]) {
    vector<string> inpfiles = {"200_226_1_S.inp", "300_328_1_S.inp", "400_452_1_S.inp", "500_546_1_S.inp", "600_661_1_S.inp",
                               "200_226_B_S.inp", "300_328_B_S.inp", "400_452_B_S.inp", "500_546_B_S.inp", "600_661_B_S.inp",
                               "200_226_U_S.inp", "300_328_U_S.inp", "400_452_U_S.inp", "500_546_U_S.inp", "600_661_U_S.inp"};

    string pipefile = "./cases/pipeTypes.txt";
    vector<int> caseScale = {200, 300, 400, 500, 600, 200, 300, 400, 500, 600, 200, 300, 400, 500, 600};
    int caseindex = atoi(argv[1]);
    int method = atoi(argv[2]);
    int seed = atoi(argv[3]);
    string ainpfile = "./cases/" + inpfiles[caseindex];
    char* ainpfilename  = new char [ainpfile.length() + 1];
    strcpy(ainpfilename, ainpfile.c_str());
    double threshold = 16;

    if (method == 1) {
        TSMA atsma(caseScale[caseindex], caseScale[caseindex] * 4000, seed, caseindex, 1, ainpfilename, pipefile, threshold);
        atsma.run();
    }
    else if (method == 2) {
        ccdeRepart accde(100, caseScale[caseindex] * 4000, 200, seed, caseindex, ainpfilename, pipefile, threshold);
        accde.run();
    }
    

    delete[] ainpfilename;
    return 0;
}
