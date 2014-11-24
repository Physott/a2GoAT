#ifndef __GAnalysis3Mesons_h__
#define __GAnalysis3Mesons_h__


#include "GHistEvent.h"
#include "GCheckProton.h"
#include "GFit.h"
//#include "GHistFit.h"


class GTreeTagger;
class GTreeParticle;
class GTreeMeson;


class   GAnalysis3Mesons  : public GHistLinked
{
private:
    Bool_t              isEtap;

    GHistEvent3Mesons   hist_raw;
    GHistEvent3Mesons   hist_SubImCut;
    GHistEvent3Mesons   hist_MMCut;

    GFit3Constraints    fit3;
    GFit4Constraints    fit4;
    GFit4ConstraintsBeam    fit4Beam;

    GH1                 fit3_im;
    GH1                 fit3_cs;
    GH1                 fit3_cl;
    GHistBGSub2         fit3_Pull;
    GH1                 fit4_im;
    GH1                 fit4_cs;
    GH1                 fit4_cl;
    GHistBGSub2         fit4_Pull;
    GH1                 fit4Beam_im;
    GH1                 fit4Beam_cs;
    GH1                 fit4Beam_cl;
    GHistBGSub2         fit4Beam_Pull;

protected:

public:
    Double_t            cutSubIM[6];
    Double_t            cutMM[2];

    GAnalysis3Mesons(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram = kTRUE);
    ~GAnalysis3Mesons();

    virtual void    CalcResult();
            Bool_t  IsEtap()    const   {return isEtap;}
    virtual Int_t   Fill(Double_t x)    {}
    virtual void    Fill(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}

    void    SetCutSubIM(const Int_t subNumber, const Double_t min, const Double_t max);
    void    SetCutMM(const Double_t min, const Double_t max);
};


class   GAnalysis3MesonsProton  : public GHistLinked
{
private:
    Bool_t              isEtap;

    GHistEvent3Mesons   hist_raw;
    GHistEvent3Mesons   hist_SubImCut;
    GHistEvent3Mesons   hist_MMCut;

    GFit3Constraints    fit3;
    GFit4Constraints    fit4;
    GFit4ConstraintsBeam    fit4Beam;
    GFit4ConstraintsProton  fit4Proton;
    GFit4ConstraintsBeamProton    fit4BeamProton;

    GHistFit            hist_fit3;
    GHistFit            hist_fit4;
    GHistFit            hist_fit4Beam;
    GHistFit            hist_fit4Proton;
    GHistFit            hist_fit4BeamProton;

protected:

public:
    GAnalysis3MesonsProton(const char* name, const char* title, const Bool_t IsEtap, Bool_t linkHistogram = kTRUE);
    ~GAnalysis3MesonsProton();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)    {}
    virtual void        Fill(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual void        ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}

    void	SetHistMeson(const Double_t sub0_min, const Double_t sub0_max,
                         const Double_t sub1_min, const Double_t sub1_max,
                         const Double_t sub2_min, const Double_t sub2_max,
                         const Double_t mm_min, const Double_t mm_max);
    void    SetFitMeson(const Double_t fit3_CutConfidenceLevel, const Double_t fit4_CutConfidenceLevel);
    void	SetCheckProton(const Double_t maxProtonAngleDiff, const Double_t minCoplanarity,const Double_t maxCoplanarity);
    void	SetHistMesonProton(const Double_t sub0_min, const Double_t sub0_max,
                               const Double_t sub1_min, const Double_t sub1_max,
                               const Double_t sub2_min, const Double_t sub2_max,
                               const Double_t mm_min, const Double_t mm_max);
    void    SetFitMesonProton(const Double_t fit3_CutConfidenceLevel, const Double_t fit4_CutConfidenceLevel);
};

#endif
