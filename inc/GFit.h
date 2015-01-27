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
    Bool_t                  solved;
    GIterativeKinFitter     fitter;

public:
    GFit(const Int_t npart, const Int_t ncon)   : solved(kFALSE), fitter(npart, ncon)   {}
    ~GFit()                                                             {}

    virtual Double_t        ConfidenceLevel()               {return fitter.ConfidenceLevel();}
    virtual Double_t        ConstraintsConfidenceLevel()    {return fitter.ConstraintsConfidenceLevel();}
    virtual Double_t        VariablesConfidenceLevel()      {return fitter.VariablesConfidenceLevel();}
    virtual Int_t           GetIterations()                 {return fitter.GetIterations()-1;}
    virtual TLorentzVector  GetTotal()                  = 0;
    virtual TLorentzVector  GetMeson()                  = 0;
    virtual TLorentzVector  GetSub(const int i)         = 0;
    virtual TLorentzVector  GetRecoil()                 = 0;
    virtual Double_t        GetChi2()                       {return fitter.GetChi2();}
    virtual Double_t        GetVariablesChi2()              {return fitter.GetVariablesChi2();}
    virtual Double_t        GetConstraintsChi2()            {return fitter.GetConstraintsChi2();}
            Double_t        GetUnknown()                    {return fitter.GetUnknown();}
    virtual Double_t        GetPull(const Int_t index)  = 0;
            Bool_t          IsSolved()                      {return solved;}
            Bool_t          Solve()                         {if(fitter.Solve()>0) solved = kTRUE; else solved = kFALSE; return solved;}
};


class	GFit3Constraints    : public GFit
{
private:
    Bool_t              isEtap;

public:
    GFit3Constraints(const Bool_t _IsEtap);
    ~GFit3Constraints();

    virtual TLorentzVector  GetMeson()                      {return fitter.GetTotalFitParticle();}
    virtual TLorentzVector  GetTotal()                      {return GetMeson();}
    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i)+fitter.GetParticle((2*i)+1);}
    virtual TLorentzVector  GetRecoil()                     {return TLorentzVector(0.0, 0.0, 0.0, 938.27);}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5);
};




class	GFit4Constraints    : public GFit
{
private:
    Bool_t              isEtap;
    TLorentzVector      beamAndTarget;

public:
    GFit4Constraints(const Bool_t _IsEtap);
    ~GFit4Constraints();

    virtual TLorentzVector  GetMeson()                      {return fitter.GetTotalFitParticle();}
    virtual TLorentzVector  GetTotal()                      {return beamAndTarget - GetMeson();}
    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i)+fitter.GetParticle((2*i)+1);}
    virtual TLorentzVector  GetRecoil()                     {return TLorentzVector(0.0, 0.0, 0.0, 938.27);}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& _BeamAndTarget);
};








class	GFit4ConstraintsBeam    : public GFit
{
private:
    Bool_t              isEtap;

public:
    GFit4ConstraintsBeam(const Bool_t _IsEtap);
    ~GFit4ConstraintsBeam();

    virtual TLorentzVector  GetMeson();
    virtual TLorentzVector  GetTotal()                      {return -fitter.GetParticle(6) - GetMeson();}
    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i)+fitter.GetParticle((2*i)+1);}
    virtual TLorentzVector  GetRecoil()                     {return TLorentzVector(0.0, 0.0, 0.0, 938.27);}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget);
};












class	GFit7ConstraintsProton    : public GFit
{
private:
    Bool_t              isEtap;
    TLorentzVector      beamAndTarget;

public:
    GFit7ConstraintsProton(const Bool_t _IsEtap);
    ~GFit7ConstraintsProton();

    virtual TLorentzVector  GetMeson();
    virtual TLorentzVector  GetTotal()                      {return beamAndTarget - GetMeson() - fitter.GetParticle(6);}
    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i)+fitter.GetParticle((2*i)+1);}
    virtual TLorentzVector  GetRecoil()                     {return fitter.GetParticle(6);}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& _BeamAndTarget,
                const TLorentzVector& proton);
};





class	GFit7ConstraintsBeamProton    : public GFit
{
private:
    Bool_t              isEtap;

public:
    GFit7ConstraintsBeamProton(const Bool_t _IsEtap);
    ~GFit7ConstraintsBeamProton();

    virtual TLorentzVector  GetMeson();
    virtual TLorentzVector  GetTotal()                      {return -fitter.GetParticle(6) - GetMeson() - fitter.GetParticle(7);}
    virtual TLorentzVector  GetSub(const int i)             {return fitter.GetParticle(2*i)+fitter.GetParticle((2*i)+1);}
    virtual TLorentzVector  GetRecoil()                     {return fitter.GetParticle(6);}
    virtual Double_t        GetPull(const Int_t index)      {return fitter.Pull(index);}
    void    Set(const TLorentzVector& p0,
                const TLorentzVector& p1,
                const TLorentzVector& p2,
                const TLorentzVector& p3,
                const TLorentzVector& p4,
                const TLorentzVector& p5,
                const TLorentzVector& beamAndTarget,
                const TLorentzVector& proton);
};









class   GHistFit    : public    GHistLinked
{
private:
    Int_t       nPulls;
    GH1         im;
    GH1         sub0im;
    GH1         sub1im;
    GH1         sub2im;
    GH1         theta;
    GH1         phi;
    GH1         Pim;
    GH1         Ptheta;
    GH1         Pphi;
    GH1         zVertex;
    GH1         chiSq;
    GH1         VchiSq;
    GH1         CchiSq;
    GH1         confidenceLevel;
    GH1         VconfidenceLevel;
    GH1         CconfidenceLevel;
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


class   GHistIterativeFit    : public    GHistLinked
{
private:
    GHistBGSub2 im;
    GHistBGSub2 sub0im;
    GHistBGSub2 sub1im;
    GHistBGSub2 sub2im;
    GHistBGSub2 mm;
    GHistBGSub2 totE;
    GHistBGSub2 totPx;
    GHistBGSub2 totPy;
    GHistBGSub2 totPz;
    GHistBGSub2 VchiSq;
    GHistBGSub2 VconfidenceLevel;
    GHistBGSub2 CchiSq;
    GHistBGSub2 CconfidenceLevel;
    GHistBGSub2 chiSq;
    GHistBGSub2 confidenceLevel;
    GHistFit    final;

public:
    GHistIterativeFit(const char* name, const char* title, const Int_t _NPulls, const Int_t _NSteps, Bool_t linkHistogram= kTRUE);
    ~GHistIterativeFit();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)                {}
    virtual Int_t       Fill(GFit& fitter);
    virtual Int_t       FillFinal(GFit& fitter, const Double_t taggerTime)                              {final.Fill(fitter, taggerTime);}
    virtual Int_t       FillFinal(GFit& fitter, const Double_t taggerTime, const Int_t taggerChannel)   {final.Fill(fitter, taggerTime, taggerChannel);}
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};


#endif
