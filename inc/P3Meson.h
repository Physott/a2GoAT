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

    PHistEvent3Meson  raw;

    Double_t cutIM[3][2];
    PHistEvent3Meson  cutIMevent;

    Double_t cutMM[2];
    PHistEvent3Meson  cutMMevent;

    GKinFitter          fit3Con;
    PHistEvent3MesonFit hist_fit3Con;
    Double_t            cutFit3ConConfidenceLevel;
    PHistEvent3MesonFit hist_fit3Con_cutCL;
    GKinFitter          fit4Con;
    PHistEvent3MesonFit hist_fit4Con;
    Double_t            cutFit4ConConfidenceLevel;
    PHistEvent3MesonFit hist_fit4Con_cutCL;
	
    Int_t 	nFound;

    TFile*  GammaResFile;
    TH2F*   GammaEloss;
    TH2F*   GammaERes;
    TH2F*   GammaThetaRes;
    TH2F*   GammaPhiRes;

    Double_t    im;
    Double_t    imSub[3];
    Double_t    im_fit;
    Double_t    imSub_fit[3];
    Double_t    Pull[24];
    GKinFitterParticle  pho[6];
    TLorentzVector      fittedMeson;
    TLorentzVector      fittedSubParticles[6];
    Double_t    misMass;
    Double_t    conLevel;

    Bool_t  ReconstructEvent(const GTreeMeson& meson, const GTreeTagger &tagger);
    Bool_t  ReconstructTagger(const GTreeMeson& meson, const GTreeTagger &tagger);

    Bool_t  fitInit(const GTreeMeson& meson, GKinFitter &fitter);
    Bool_t  DoFit3Con(const GTreeMeson& meson);
    Bool_t  DoFit4Con(const GTreeMeson& meson, const TLorentzVector &beamPlusTarget);

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
