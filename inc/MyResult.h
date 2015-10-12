#ifndef __MyResult_h__
#define __MyResult_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "MRFitTaggerBins.h"



class	MyResult
{
private:
    TFile*             in;
    TFile*             out;
    MRFitTaggerBins    result;

protected:

			
public:
    MyResult(const char* _Name, TFile* input, TFile* output, const int _Color);
    ~MyResult();

    Bool_t  Start(const int summedTaggerChannels = 1, const int summedIMBins = 1);
};
#endif
