#ifndef __GFit_h__
#define __GFit_h__


#include <TH2.h>

#include "GHistManager.h"
#include "GH1.h"
#include "GKinFitter.h"
#include "GKinFitterParticle.h"



class	GFit
{
protected:
    GKinFitter fitter;

public:
    GFit(const Int_t npart, const Int_t ncon): fitter(npart, ncon, 0)   {}
    ~GFit()                                                             {}

    virtual Double_t        ConfidenceLevel()           = 0;
            //Int_t           GetNIter()                  {return fitter.GetNIter();}
    virtual TLorentzVector  GetTotalFitParticle()       = 0;
    virtual TLorentzVector  GetSub(const int i)         = 0;
    virtual TLorentzVector  GetRecoil()                 = 0;
    virtual Double_t        GetChi2()                   = 0;
    virtual Double_t        GetPull(const Int_t index)  = 0;
    virtual Bool_t          IsSolved()                  = 0;
};


class	GFit3Constraints    : public GFit
{
private:
    Bool_t              solved;
    Bool_t              isEtap;

public:
    GFit3Constraints(const Bool_t _IsEtap);
    ~GFit3Constraints();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle()           {return fitter.GetTotalFitParticle().Get4Vector();}
    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i).Get4Vector()+fitter.GetParticle((2*i)+1).Get4Vector();}
    virtual TLorentzVector  GetRecoil()                     {return TLorentzVector(0.0, 0.0, 0.0, 938.27);}
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

public:
    GFit4Constraints(const Bool_t _IsEtap);
    ~GFit4Constraints();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle()           {return fitter.GetTotalFitParticle().Get4Vector();}

    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i).Get4Vector()+fitter.GetParticle((2*i)+1).Get4Vector();}
    virtual TLorentzVector  GetRecoil()                     {return TLorentzVector(0.0, 0.0, 0.0, 938.27);}
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

public:
    GFit4ConstraintsBeam(const Bool_t _IsEtap);
    ~GFit4ConstraintsBeam();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();

    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i).Get4Vector()+fitter.GetParticle((2*i)+1).Get4Vector();}
    virtual TLorentzVector  GetRecoil()                     {return TLorentzVector(0.0, 0.0, 0.0, 938.27);}
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












class	GFit7ConstraintsProton    : public GFit
{
private:
    Bool_t              solved;
    Bool_t              isEtap;

public:
    GFit7ConstraintsProton(const Bool_t _IsEtap);
    ~GFit7ConstraintsProton();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();

    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i).Get4Vector()+fitter.GetParticle((2*i)+1).Get4Vector();}
    virtual TLorentzVector  GetRecoil()                     {return fitter.GetParticle(6).Get4Vector();}
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





class	GFit7ConstraintsBeamProton    : public GFit
{
private:
    Bool_t              solved;
    Bool_t              isEtap;

public:
    GFit7ConstraintsBeamProton(const Bool_t _IsEtap);
    ~GFit7ConstraintsBeamProton();

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual TLorentzVector  GetTotalFitParticle();

    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i).Get4Vector()+fitter.GetParticle((2*i)+1).Get4Vector();}
    virtual TLorentzVector  GetRecoil()                     {return fitter.GetParticle(6).Get4Vector();}
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









class   GHistFit    : public    GHistLinked
{
private:
    Int_t       nPulls;
    GH1         nIter;
    GH1         im;
    GH1         sub0im;
    GH1         sub1im;
    GH1         sub2im;
    GH1         theta;
    GH1         phi;
    GH1         Pim;
    GH1         Ptheta;
    GH1         Pphi;
    GH1         chiSq;
    GH1         confidenceLevel;
    GHistBGSub2 pulls;
public:
    GHistFit(const char* name, const char* title, const Int_t _NPulls, Bool_t linkHistogram= kTRUE);
    ~GHistFit();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)                {}
    virtual Int_t       Fill(GFit& fitter, const Double_t taggerTime);
    virtual Int_t       Fill(GFit& fitter, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};


#endif
