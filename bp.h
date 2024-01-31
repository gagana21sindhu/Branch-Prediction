#ifndef BP_H
#define BP_H

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <iomanip>
#include <vector>
#include <deque>

using namespace std;

class bp {
    public:
        unsigned int PCaddr = 0;
        unsigned int index = 0;
        int k = 0;                   // chooser table bits
        int m1 = 0;                  // gshare prediction table index bits 
        int n = 0;                   // Global BH registery bits
        int m2 = 0;                  // bimodal prediction table index bits

        int bimodalRows = 0;
        int gshareRows = 0;
        int chooserRows = 0;
        int *bimodalPredTable;
        int *gsharePredTable;
        int *chooserTable;
        unsigned int histReg = 0;
        

        unsigned int bimodalMask = 0;
        unsigned int gshareMask = 0;
        unsigned int chooserMask = 0;
        unsigned int hRegMask = 0;


        int predictions = 0;
        int mispredictions = 0;

        bp(int a, int b, int c, int d); // Constructor
        
        void initbp();
        void bimodalPrint();
        void gsharePrint();
        void hybridPrint();

        void predict(unsigned int PC, const char *tn);
};


#endif

