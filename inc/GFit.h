#ifndef __GFit_h__
#define __GFit_h__


#include <TH2.h>

#include "GHistManager.h"
#include "GH1.h"
#include "GKinFitter.h"
#include "GKinFitterWithProton.h"


/*
class	GFit
{
private:

public:
    GFit()  {}
    ~GFit() {}

    virtual Double_t        ConfidenceLevel()           = 0;
    virtual TLorentzVector  GetTotalFitParticle()       = 0;
    virtual Double_t        GetChi2()                   = 0;
    virtual Double_t        GetPull(const Int_t index)  = 0;
    virtual Bool_t          IsSolved()                  = 0;
};


class	GFit3Constraints    : public GFit
{
private:
    Bool_t              solved;
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit3Constraints(const Bool_t _IsEtap);
    ~GFit3Constraints();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle()           {return fitter.GetTotalFitParticle().Get4Vector();}
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    virtual Bool_t          IsSolved()                      {return solved;}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) solved = kTRUE; else solved = kFALSE; return solved;}
};




class	GFit4Constraints    : public GFit
{
private:
    Bool_t              solved;
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4Constraints(const Bool_t _IsEtap);
    ~GFit4Constraints();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle()           {return fitter.GetTotalFitParticle().Get4Vector();}
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    virtual Bool_t          IsSolved()                      {return solved;}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) solved = kTRUE; else solved = kFALSE; return solved;}
};








class	GFit4ConstraintsBeam    : public GFit
{
private:
    Bool_t              solved;
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsBeam(const Bool_t _IsEtap);
    ~GFit4ConstraintsBeam();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    virtual Bool_t          IsSolved()                      {return solved;}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) solved = kTRUE; else solved = kFALSE; return solved;}
};












class	GFit4ConstraintsProton    : public GFit
{
private:
    Bool_t              solved;
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsProton(const Bool_t _IsEtap);
    ~GFit4ConstraintsProton();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    virtual Bool_t          IsSolved()                      {return solved;}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget,
                const TLorentzVector& proton);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) solved = kTRUE; else solved = kFALSE; return solved;}
};




class	GFit4ConstraintsProtonExact    : public GFit
{
private:
    Bool_t              solved;
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsProtonExact(const Bool_t _IsEtap);
    ~GFit4ConstraintsProtonExact();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    virtual Bool_t          IsSolved()                      {return solved;}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget,
                const TLorentzVector& proton);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) solved = kTRUE; else solved = kFALSE; return solved;}
};








class	GFit4ConstraintsBeamProton    : public GFit
{
private:
    Bool_t              solved;
    Bool_t              isEtap;
    GKinFitter          fitter;

public:
    GFit4ConstraintsBeamProton(const Bool_t _IsEtap);
    ~GFit4ConstraintsBeamProton();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    virtual Bool_t          IsSolved()                      {return solved;}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget,
                const TLorentzVector& proton);
    Bool_t  Solve()                                 {if(fitter.Solve()>0) solved = kTRUE; else solved = kFALSE; return solved;}
};



*/





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
    virtual Int_t       Fill(GKinFitterWithProton& fitter, const Double_t taggerTime);
    virtual Int_t       Fill(GKinFitterWithProton& fitter, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};


#endif
