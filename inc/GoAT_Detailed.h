#ifndef __GoAT_Detailed_h__
#define __GoAT_Detailed_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 


#include "GSort.h"
#include "GParticleReconstruction.h"
#include "GMesonReconstruction_6and7gamma.h"


class	GoAT_Detailed : public GSort, public GParticleReconstruction, public GMesonReconstruction_6and7gamma
{
private:
	Int_t	UsePeriodMacro;
	Int_t 	period;
	
	Bool_t 	UseParticleReconstruction;
    Bool_t 	UseMesonReconstruction;

	Int_t 	nEvents_written;

    GH1             AcceptanceTrue;

protected:
    virtual void 	ProcessEvent();
    virtual Bool_t	Start();


public:
    GoAT_Detailed();
    virtual ~GoAT_Detailed();

    virtual Bool_t	Init(const char* configfile);
};
#endif
