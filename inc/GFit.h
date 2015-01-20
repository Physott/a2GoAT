#ifndef __GFit_h__
#define __GFit_h__


#include <TH2.h>

#include "GHistManager.h"
#include "GH1.h"
#include "GKinFitter.h"
#include "GKinFitterWithProton.h"


class   GHistFit    : public    GHistLinked
{
private:
    Int_t       nPulls;
    GH1         im;
    GH1         im_sub0;
    GH1         im_sub1;
    GH1         im_sub2;
    GH1         chiSq;
    GH1         confidenceLevel;
    GH1         vertex_X;
    GH1         vertex_Y;
    GH1         vertex_Z;
    GH1         vertex_sub0_X;
    GH1         vertex_sub0_Y;
    GH1         vertex_sub0_Z;
    GH1         vertex_sub1_X;
    GH1         vertex_sub1_Y;
    GH1         vertex_sub1_Z;
    GH1         vertex_sub2_X;
    GH1         vertex_sub2_Y;
    GH1         vertex_sub2_Z;
    GHistBGSub2 pulls;

public:
    GHistFit(const char* name, const char* title, const Int_t _NPulls, Bool_t linkHistogram= kTRUE);
    ~GHistFit();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)                {}
    virtual Int_t       Fill(GKinFitter& fitter, const Double_t taggerTime);
    virtual Int_t       Fill(GKinFitter& fitter, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};


class   GHistFit2    : public    GHistLinked
{
private:
    Int_t               count;

    Int_t               nPulls;
    GHistBGSub2         im;
    GHistBGSub2         im_sub0;
    GHistBGSub2         im_sub1;
    GHistBGSub2         im_sub2;
    GHistBGSub2         chiSq;
    GHistBGSub2         confidenceLevel;
    GHistBGSub2         vertex_X;
    GHistBGSub2         vertex_Y;
    GHistBGSub2         vertex_Z;
    GHistBGSub2         vertex_sub0_X;
    GHistBGSub2         vertex_sub0_Y;
    GHistBGSub2         vertex_sub0_Z;
    GHistBGSub2         vertex_sub1_X;
    GHistBGSub2         vertex_sub1_Y;
    GHistBGSub2         vertex_sub1_Z;
    GHistBGSub2         vertex_sub2_X;
    GHistBGSub2         vertex_sub2_Y;
    GHistBGSub2         vertex_sub2_Z;
    GHistBGSub2         pulls0;
    GHistBGSub2         pulls1;
    GHistBGSub2         pulls2;
    GHistBGSub2         pulls3;
    GHistBGSub2         pulls4;

public:
    GHistFit2(const char* name, const char* title, const Int_t _NPulls, Bool_t linkHistogram= kTRUE);
    ~GHistFit2();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)                {}
    virtual Int_t       Fill(GKinFitter& fitter, const Double_t taggerTime);
    virtual Int_t       Fill(GKinFitter& fitter, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
            void        ResetCount()                    {count = 0;}
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};


#endif
