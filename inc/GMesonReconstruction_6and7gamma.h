#ifndef __GMesonReconstruction_6and7gamma_h__
#define __GMesonReconstruction_6and7gamma_h__


#include "GTreeManager.h"
#include "GH1.h"
#include "TCutG.h"

#define DEFAULT_PI0_IM_WIDTH 20.0
#define DEFAULT_ETA_IM_WIDTH 44.0
#define DEFAULT_ETAP_IM_WIDTH 60.0


class  GMesonReconstruction_6and7gamma  : virtual public GTreeManager
{
private:
    Double_t	width_pi0;
    Double_t	width_eta;
    Double_t	width_etap;

    Double_t	meson_theta_min;
    Double_t	meson_theta_max;

    Double_t    minChiSq;
    UInt_t      minDecayIndex;
    UInt_t      minIndex;
    TLorentzVector  reconstructedEta;
    TLorentzVector  reconstructedEtap;
    Int_t           daughter_index[6];
    TLorentzVector  daughter[6];
    UInt_t          foundTaggerHitForProton;
    TLorentzVector  mesonHelp[15][3];

    static  Int_t   perm6g[15][6];

    GHistBGSub      particleTime;
    GHistBGSub      particleTimeCut;
    GHistBGSub      IMSub0Etap;
    GHistBGSub      IMSub03Pi0;
    GHistBGSub      IMSub1Etap;
    GHistBGSub      IMSub13Pi0;
    GHistBGSub      IMSub2Etap;
    GHistBGSub      IMSub23Pi0;
    GHistBGSub2     TOFAddPhoton6Hits;
    GHistBGSub2     TOFSubPhoton6Hits;
    GHistBGSub2     TOFAddPhoton7Hits;
    GHistBGSub2     TOFSubPhoton7Hits;
    GHistBGSub2     TOFAddProton;
    GHistBGSub2     TOFSubProton;

    GHistBGSub2     ChiSqDist;
    GHistBGSub2     ChiSqDist6Hits;
    GHistBGSub2     ChiSqDist7Hits;
    Double_t        ChiSqEtap;
    Double_t        ChiSq3Pi0;


    void    Reconstruct6g();
    void    Reconstruct6g(TLorentzVector** vec);
    bool Reconstruct7g();
    Bool_t  CheckProton();

    TLorentzVector  SetMass(const TLorentzVector& vec, const Double_t mass);
    TCutG*          OpenCutFile(Char_t* filename, Char_t* cutname);
protected:

            Bool_t  ProcessEventWithoutFilling();
    virtual void    ProcessEvent();
    virtual Bool_t  Start();

public:
    GMesonReconstruction_6and7gamma();
    virtual ~GMesonReconstruction_6and7gamma();

            Bool_t  Init();
            void    SetPi0Width(const Double_t width)   {width_pi0 = width;}
            void    SetEtaWidth(const Double_t width)   {width_eta = width;}
            void    SetEtapWidth(const Double_t width)  {width_etap = width;}
};





#endif
