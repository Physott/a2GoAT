#ifndef __MyEtap_h__
#define __MyEtap_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "PPhysics.h"
#include "PHistEvent.h"

class	MyEtap : public PPhysics
{
private:
    TH1* 	time_raw;
    PHistEvent3Meson  raw;

    TH1* 	time_cutIM;
    Double_t cutIM[3][2];
    PHistEvent3Meson  cutIMevent;

    TH1* 	time_cutMM;
    Double_t cutMM[2];
    PHistEvent3Meson  cutMMevent;
	
    Int_t 	N_eta;

protected:
    virtual Bool_t  Start();
    
    virtual void    ProcessEvent();

            Bool_t	WriteHistograms();
			
public:
    MyEtap();
    virtual ~MyEtap();

    virtual Bool_t	Init(const char* configfile);

};
#endif
