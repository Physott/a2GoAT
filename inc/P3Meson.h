#ifndef __P3Meson_h__
#define __P3Meson_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include <TH2F.h>

#include "GKinFitter.h"
#include "GKinFitterParticle.h"
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

    Double_t            cutFit3ConConfidenceLevel;
    GKinFitter          fit3Con;
    PHistEvent3MesonFit hist_fit3Con;
    Double_t            cutFit4ConConfidenceLevel;
    GKinFitter          fit4Con;
    PHistEvent3MesonFit hist_fit4Con;
	
    Int_t 	nFound;

    TFile*  GammaResFile;
    TH2F*   GammaEloss;
    TH2F*   GammaERes;
    TH2F*   GammaThetaRes;
    TH2F*   GammaPhiRes;

    Double_t    imSub[3];
    Double_t    misMass;

    void    fit(const GTreeMeson& meson, const Double_t tagger_time, const Double_t tagger_channel);

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
