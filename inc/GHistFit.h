#ifndef __GHistFit_h__
#define __GHistFit_h__


#include <TH2.h>

#include "GHistManager.h"
#include "GH1.h"

#define GKinFitterBase_MaxSteps 10

class   GKinFitterBase;

class   GHistFit    : public    GHistLinked
{
private:
    Int_t       nPar;
    Int_t       nUnk;
    Int_t       nCon;
    GH1         invMass;
    GH1         invMassSub0;
    GH1         invMassSub1;
    GH1         invMassSub2;
    GH1         chiSq;
    GH1         confidenceLevel;
    GH1         zVertex;
    GHistBGSub2 photonsE;
    GHistBGSub2 photonsTheta;
    GHistBGSub2 photonsPhi;
    GH1         protonEnergy;
    GH1         protonTheta;
    GH1         protonPhi;
    GH1         beamEnergy;
    GH1         beamTheta;
    GH1         beamPhi;
    GHistBGSub2 con;
    GHistBGSub2 pulls;

public:
    GHistFit(const char* name, const char* title, const Int_t nPulls, Bool_t linkHistogram= kTRUE);
    ~GHistFit();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)                {}
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}

    friend class GKinFitterBase;
};


class   GHistFit2    : public    GHistLinked
{
private:
    Int_t           nPar;
    Int_t           nUnk;
    Int_t           nCon;
    Int_t           nIter;
    GHistBGSub2     invMass;
    GHistBGSub2     invMassSub0;
    GHistBGSub2     invMassSub1;
    GHistBGSub2     invMassSub2;
    GHistBGSub2     chiSq;
    GHistBGSub2     confidenceLevel;
    GHistBGSub2     zVertex;
    GHistBGSub2*    photonsEnergy[6];
    GHistBGSub2*    photonsTheta[6];
    GHistBGSub2*    photonsPhi[6];
    GHistBGSub2     protonEnergy;
    GHistBGSub2     protonTheta;
    GHistBGSub2     protonPhi;
    GHistBGSub2     beamEnergy;
    GHistBGSub2     beamTheta;
    GHistBGSub2     beamPhi;
    GHistBGSub2*    con[7];
    GHistBGSub2*    pulls[23];

public:
    GHistFit2(const char* name, const char* title, Bool_t linkHistogram= kTRUE);
    ~GHistFit2();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)                {}
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}

    friend class GKinFitterBase;
};


#endif
