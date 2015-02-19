#ifndef __GHistEvent_h__
#define __GHistEvent_h__

#include <TLorentzVector.h>

#include "GH1.h"


class	GHistEvent  : public GHistLinked
{
private:

protected:
    GH1     im;
    GH1     mm;
    GHistBGSub2 EnergyTheta;
    GHistBGSub2 phiTheta;

public:
    GHistEvent(const char* name, const char* title, Bool_t linkHistogram = kTRUE);
    virtual ~GHistEvent();

    virtual void    CalcResult()                                                                                                                                {im.CalcResult(); mm.CalcResult(); EnergyTheta.CalcResult(); phiTheta.CalcResult();}
    virtual Int_t   Fill(Double_t x)                                                                                                                            {}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t Energy, const Double_t Theta, const Double_t Phi, const Double_t taggerTime)                             {im.Fill(IM, taggerTime); mm.Fill(MM, taggerTime); EnergyTheta.Fill(Energy, Theta, taggerTime); phiTheta.Fill(Phi, Theta, taggerTime);}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t Energy, const Double_t Theta, const Double_t Phi, const Double_t taggerTime, const Int_t taggerChannel)  {im.Fill(IM, taggerTime, taggerChannel); mm.Fill(MM, taggerTime, taggerChannel); EnergyTheta.Fill(Energy, Theta, taggerTime); phiTheta.Fill(Phi, Theta, taggerTime);}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "")                                                                                                                {im.Reset(option); mm.Reset(option);}
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE)                           {im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); mm.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);}
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {}
};




class	GHistEvent3Mesons   : public GHistEvent
{
private:
    GH1     sub0_im;
    GH1     sub1_im;
    GH1     sub2_im;
    GHistBGSub2     subEtaPi0_im;
    GHistBGSub2     subPi0Pi0_im;

protected:

public:
    GHistEvent3Mesons(const char* name, const char* title, Bool_t linkHistogram = kTRUE);
    virtual ~GHistEvent3Mesons();

    virtual void    CalcResult()                                                                                                                                                                {GHistEvent::CalcResult(); sub0_im.CalcResult(); sub1_im.CalcResult(); sub2_im.CalcResult(); subEtaPi0_im.CalcResult(); subPi0Pi0_im.CalcResult();}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t Energy, const Double_t Theta, const Double_t Phi, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t taggerTime)                               {GHistEvent::Fill(IM, MM, Energy, Theta, Phi, taggerTime); sub0_im.Fill(SUB0_IM, taggerTime); sub1_im.Fill(SUB1_IM, taggerTime); sub2_im.Fill(SUB2_IM, taggerTime); subEtaPi0_im.Fill(SUB0_IM, SUB1_IM, taggerTime); subEtaPi0_im.Fill(SUB0_IM, SUB2_IM, taggerTime); subPi0Pi0_im.Fill(SUB1_IM, SUB2_IM, taggerTime);}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t Energy, const Double_t Theta, const Double_t Phi, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t taggerTime, const Int_t taggerChannel)    {GHistEvent::Fill(IM, MM, Energy, Theta, Phi, taggerTime, taggerChannel); sub0_im.Fill(SUB0_IM, taggerTime); sub1_im.Fill(SUB1_IM, taggerTime); sub2_im.Fill(SUB2_IM, taggerTime); subEtaPi0_im.Fill(SUB0_IM, SUB1_IM, taggerTime); subEtaPi0_im.Fill(SUB0_IM, SUB2_IM, taggerTime); subPi0Pi0_im.Fill(SUB1_IM, SUB2_IM, taggerTime);}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option)                                                                                                                                                     {GHistEvent::Reset(option); sub0_im.Reset(option); sub1_im.Reset(option); sub2_im.Reset(option); subEtaPi0_im.Reset(option); subPi0Pi0_im.Reset(option);}
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)                                                                    {GHistEvent::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); sub0_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); sub1_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); sub2_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); subEtaPi0_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); subPi0Pi0_im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);}
};





class	GHistEvent3MesonsProton   : public GHistEvent3Mesons
{
private:
    GH1     protonEnergy;
    GH1     protontheta;
    GHistBGSub2 protonEnergyTheta;
    GHistBGSub2 protonPhiTheta;

protected:

public:
    GHistEvent3MesonsProton(const char* name, const char* title, Bool_t linkHistogram = kTRUE);
    virtual ~GHistEvent3MesonsProton();

    virtual void    CalcResult()                                                                                                                                                                {GHistEvent3Mesons::CalcResult(); protonEnergy.CalcResult(); protontheta.CalcResult(); protonEnergyTheta.CalcResult(); protonPhiTheta.CalcResult();}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t Energy, const Double_t Theta, const Double_t Phi, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t ProtonE, const Double_t ProtonTh, const Double_t ProtonPh, const Double_t taggerTime)                               {GHistEvent3Mesons::Fill(IM, MM, Energy, Theta, Phi, SUB0_IM, SUB1_IM, SUB2_IM, taggerTime); protonEnergy.Fill(ProtonE, taggerTime); protontheta.Fill(ProtonTh, taggerTime); protonEnergyTheta.Fill(ProtonE, ProtonTh, taggerTime); protonPhiTheta.Fill(ProtonPh, ProtonTh, taggerTime);}
    virtual void    Fill(const Double_t IM, const Double_t MM, const Double_t Energy, const Double_t Theta, const Double_t Phi, const Double_t SUB0_IM, const Double_t SUB1_IM, const Double_t SUB2_IM, const Double_t ProtonE, const Double_t ProtonTh, const Double_t ProtonPh, const Double_t taggerTime, const Int_t taggerChannel)    {GHistEvent3Mesons::Fill(IM, MM, Energy, Theta, Phi, SUB0_IM, SUB1_IM, SUB2_IM, taggerTime, taggerChannel); protonEnergy.Fill(ProtonE, taggerTime, taggerChannel); protontheta.Fill(ProtonTh, taggerTime, taggerChannel); protonEnergyTheta.Fill(ProtonE, ProtonTh, taggerTime); protonPhiTheta.Fill(ProtonPh, ProtonTh, taggerTime);}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option)                                                                                                                                                     {GHistEvent3Mesons::Reset(option); protonEnergy.Reset(option); protontheta.Reset(option); protonEnergyTheta.Reset(option); protonPhiTheta.Reset(option);}
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)                                                                    {GHistEvent3Mesons::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); protonEnergy.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); protontheta.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); protonEnergyTheta.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads); protonPhiTheta.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);}
};


#endif
