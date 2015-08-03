#ifndef __MRRecEff_h__
#define __MRRecEff_h__

#include <TMath.h>
#include <TFile.h>
#include <TString.h>
#include <TH3.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "MRFitTaggerThetaBins.h"


class	MRRecEff
{
private:
    TH1D*	acquCount;
    TH1D*	acquCount6;
    TH1D*	acquCount7;
    TH1D*	acceptedCount7;

    TH2D*	acquCount2D;
    TH2D*	acquCount2D6;
    TH2D*	acquCount2D7;
    TH2D*	acceptedCount2D7;

    TH1D*   factor;
    TH1D*   factor7;
    TH1D*   recEff;
    TH1D*   recEff7;

    TH2D*   factor2D;
    TH2D*   factor2D7;
    TH2D*   recEff2D;
    TH2D*   recEff2D7;

public:
    MRRecEff(TFile* acquSignalFile, TFile* SignalFile);
    virtual ~MRRecEff();

    TH1D*   GetFactor() {return factor;}
    TH1D*   GetRecEff() {return recEff;}
    TH2D*   GetFactor2D() {return factor2D;}
    TH2D*   GetRecEff2D() {return recEff2D;}
    void    Draw(TFile *out);
    void    CalcResult();
};





#endif
