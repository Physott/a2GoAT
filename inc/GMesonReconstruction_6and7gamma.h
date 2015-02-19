#ifndef __GMesonReconstruction_6and7gamma_h__
#define __GMesonReconstruction_6and7gamma_h__


#include "GTreeManager.h"
#include "GH1.h"

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
    Double_t    minChiSqEtap;
    Double_t    minChiSq3Pi0;
    UInt_t      minDecayIndex;
    UInt_t      minIndex;
    Int_t       bestIndex_7g;
    TLorentzVector  reconstructedEta;
    TLorentzVector  reconstructedEtap;
    Int_t           daughter_index[6];
    TLorentzVector* daughter[6];
    UInt_t          foundTaggerHitForProton;

    static  Int_t   perm6g[15][6];

    GHistBGSub2 histChiCL_6g;
    GHistBGSub2 histChiCL_7g;
    GHistBGSub  histSteps_7g;

    void    Reconstruct6g();
    void    Reconstruct6g(TLorentzVector** vec);
    void    Reconstruct7g();
    Bool_t  CheckProton();

    TLorentzVector  SetMass(const TLorentzVector& vec, const Double_t mass);

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
