#ifndef __P3Meson_h__
#define __P3Meson_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "PPhysics.h"
#include "PHistEvent.h"
#include "PProtonCheck.h"

class	P3Meson
{
private:
    TString name;
    Bool_t  isEtap;

    TH1D 	time_raw;
    PHistEvent3Meson  raw;

    TH1D 	time_cutIM;
    Double_t cutIM[3][2];
    PHistEvent3Meson  cutIMevent;

    TH1D 	time_cutMM;
    Double_t cutMM[2];
    PHistEvent3Meson  cutMMevent;
	
    Int_t 	nFound;

protected:
			
public:
    P3Meson(const TString& _Name, const Bool_t _IsEtap = kFALSE);
    virtual ~P3Meson();

            void	Clear();
            Int_t   GetNFound() const   {return nFound;}
            Bool_t  ProcessEvent(const GTreeMeson &meson, const GTreeTagger &tagger);
            Bool_t  Write(TDirectory& curDir);
};

#endif
