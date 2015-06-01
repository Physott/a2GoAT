#ifndef __ScaleMC_h__
#define __ScaleMC_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "GAnalysis3Mesons.h"

class	ScaleMC  : public GTreeManager
{
private:
    GHistBGSub2     CalibCB;
    GHistBGSub2     CalibTAPS;
    GHistBGSub2     CalibCBCorr;
    GHistBGSub2     CalibTAPSCorr;
    GH1             AcceptanceTrue;

protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();

public:
    ScaleMC();
    virtual ~ScaleMC();

    virtual Bool_t	Init(const char* configfile);

};
#endif
