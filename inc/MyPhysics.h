#ifndef __MyPhysics_h__
#define __MyPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GAnalysis3Mesons.h"

class	MyPhysics  : public GTreeManager
{
private:
    GH1     EPTscalers;
    GH1     EPTscalersCor;
    TH1D    EPTscalersT;
    TH1D    EPTscalersCorT;

protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    MyPhysics();
    virtual ~MyPhysics();

    virtual Bool_t	Init(const char* configfile);

};
#endif
