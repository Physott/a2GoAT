#ifndef __GAnalysis3Mesons_h__
#define __GAnalysis3Mesons_h__


#include "GHistEvent.h"
#include "GCheckProton.h"
#include "GFit.h"
#include "GFitBeam.h"
#include "GFitProton.h"
#include "GFitBeamProton.h"
#include "APLCON.hpp"


class GTreeTagger;
class GTreeParticle;
class GTreeMeson;


class   GAnalysis3Mesons  : public GHistLinked
{
private:
    Bool_t              isEtap;

    GHistBGSub          hist_count;

    GHistEvent3Mesons   hist_raw;
    GHistBGSub          hist_raw_SubAll;
    GHistEvent3Mesons   hist_SubImCut;
    GHistBGSub          hist_SubImCut_SubAll;

    GHistEvent3Mesons   hist_MMCut;
    GHistBGSub          hist_SubAll;

//    GFit                fit1;
//    GHistEvent3Mesons   hist_fit1;
//    GFitVertex          fit1Vertex;
//    GHistEvent3Mesons   hist_fit1Vertex;
//    GFit                fit3;
//    GHistEvent3Mesons   hist_fit3;
//    GFitVertex          fit3Vertex;
//    GHistEvent3Mesons   hist_fit3Vertex;
    GFit                fit4;
    GHistEvent3Mesons   hist_fit4;
    GHistBGSub          hist_fit4_SubAll;
//    GFitVertex          fit4Vertex;
//    GHistEvent3Mesons   hist_fit4Vertex;

//    GFitBeam            fitBeam1;
//    GHistEvent3Mesons   hist_fitBeam1;
//    GFitBeamVertex      fitBeam1Vertex;
//    GHistEvent3Mesons   hist_fitBeam1Vertex;
//    GFitBeam            fitBeam3;
//    GHistEvent3Mesons   hist_fitBeam3;
//    GFitBeamVertex      fitBeam3Vertex;
//    GHistEvent3Mesons   hist_fitBeam3Vertex;
//    GFitBeam            fitBeam4;
//    GHistEvent3Mesons   hist_fitBeam4;
//    GFitBeamVertex      fitBeam4Vertex;
//    GHistEvent3Mesons   hist_fitBeam4Vertex;

    Double_t            cutSubIM[6];
    Double_t            cutMM[2];

    bool                success;

protected:

public:

    GAnalysis3Mesons(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram = kTRUE);
    ~GAnalysis3Mesons();

    virtual void    CalcResult();
            Bool_t  IsEtap()    const   {return isEtap;}
    virtual Int_t   Fill(Double_t x)    {return 0;}
    virtual void    Fill(const GTreeMeson& meson, GTreeParticle& photons, const GTreeTagger& tagger, const GTreeA2Geant &geantTree);
            bool    IsSuccess() const   {return success;}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {return 0;}

    void    SetCutSubIM(const Int_t subNumber, const Double_t min, const Double_t max);
    void    SetCutMM(const Double_t min, const Double_t max);
};


class   GAnalysis3MesonsProton  : public GHistLinked
{
private:
    Bool_t              isEtap;

    GHistBGSub          hist_count;

    GCheckProton        checkProton;

    GHistEvent3MesonsProton   hist_raw;
    GHistBGSub2         hist_raw_TOF;
    GHistBGSub          hist_raw_SubAll;

    GHistEvent3MesonsProton   hist_SubImCut;
    GHistBGSub2         hist_SubImCut_TOF;
    GHistBGSub          hist_SubImCut_SubAll;

    GHistEvent3MesonsProton   hist_MMCut;
    GHistBGSub2         hist_TOF;
    GHistBGSub          hist_SubAll;

//    GFit                fit4;
//    GHistEvent3Mesons   hist_fit4;
//    GFitVertex          fit4Vertex;
//    GHistEvent3Mesons   hist_fit4Vertex;

//    GFitBeam            fitBeam4;
//    GHistEvent3Mesons   hist_fitBeam4;
//    GFitBeamVertex      fitBeam4Vertex;
//    GHistEvent3Mesons   hist_fitBeam4Vertex;

    GFitProton          fitProton6;
    GHistEvent3MesonsProton   hist_fitProton6;
    GHistBGSub          hist_fitProton6_SubAll;
//    GFitProtonVertex    fitProton6Vertex;
//    GHistEvent3Mesons   hist_fitProton6Vertex;

//    GFitBeamProton          fitBeamProton6;
//    GHistEvent3Mesons       hist_fitBeamProton6;
//    GFitBeamProtonVertex    fitBeamProton6Vertex;
//    GHistEvent3Mesons       hist_fitBeamProton6Vertex;

    Double_t            cutSubIM[6];
    Double_t            cutMM[2];

    bool                success;

protected:

public:
    GAnalysis3MesonsProton(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram = kTRUE);
    ~GAnalysis3MesonsProton();

    virtual void        CalcResult();
    virtual Int_t       Fill(Double_t x)    {return 0;}
    virtual void        Fill(const GTreeMeson& meson, GTreeParticle &photons, GTreeParticle &proton, const GTreeTagger& tagger, const GTreeA2Geant &geantTree);
            bool        IsSuccess() const   {return success;}
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void        Reset(Option_t* option = "");
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {return 0;}

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

    void    SetCutSubIM(const Int_t subNumber, const Double_t min, const Double_t max);
    void    SetCutMM(const Double_t min, const Double_t max);
};

#endif
