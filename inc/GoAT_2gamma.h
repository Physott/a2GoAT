#ifndef __GoAT_2gamma_h__
#define __GoAT_2gamma_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "GH1.h"


class	GoAT_2gamma : public GTreeManager
{
private:
    Int_t 	nEvents_written;

    GH1     im;
    GH1     imCB;
    GH1     imTAPS;
    GH1     im2gamma;
    GH1     im2gammaCB;
    GH1     im2gammaTAPS;
    GH1     im2gammaProton;
    GH1     im2gammaProtonCB;
    GH1     im2gammaProtonTAPS;

protected:
    virtual void 	ProcessEvent();
    virtual Bool_t	Start();


public:
    GoAT_2gamma();
    virtual ~GoAT_2gamma();

    virtual Bool_t	Init(const char* configfile);
};
#endif
