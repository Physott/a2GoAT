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
    //GAnalysis3Mesons        hist_eta;
    //GAnalysis3MesonsProton  hist_eta_proton;
    //GAnalysis3Mesons        hist_etap;
    //GAnalysis3MesonsProton  hist_etap_proton;

    //GH1     EPTscalers;
    //GH1     EPTscalersCor;
    //TH1D    EPTscalersT;
    //TH1D    EPTscalersCorT;

    GHistBGSub2     TOF;
    GHistBGSub2     TOFPhoton;
    GHistBGSub2     TOFProton;
    GHistBGSub2     TOFvsTagger;
    GHistBGSub2     TOFvsTaggerPhoton;
    GHistBGSub2     TOFvsTaggerProton;
    GHistBGSub2     TOFvsTaggerPhotonTheta;

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
