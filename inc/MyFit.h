#ifndef __MyFit_h__
#define __MyFit_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GAnalysis3Mesons.h"

class	MyFit  : public GTreeManager
{
private:
    //GAnalysis3Mesons        hist_eta;
    //GAnalysis3MesonsProton  hist_eta_proton;
    GAnalysis3Mesons        hist_etap;
    GAnalysis3MesonsProton  hist_etap_proton;

    GH1     EPTscalers;
    GH1     EPTscalersCor;
    TH1D    EPTscalersT;
    TH1D    EPTscalersCorT;

    GH1             AcceptanceTrue;
    GH1             AcceptanceProtonTrue;

protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    MyFit();
    virtual ~MyFit();

    virtual Bool_t	Init(const char* configfile);

};
#endif
