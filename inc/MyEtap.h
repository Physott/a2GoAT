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

    TH1* 	time_eta;

    PHistEvent  raw;

    Double_t cutIM[3][2];
    PHistD*	CutIM_IM_sub[3];
    PHistD*	CutIM_IM_eta;
    PHistD*	CutIM_MM_eta;
	
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
