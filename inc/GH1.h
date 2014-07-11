#ifndef __GH1_h__
#define __GH1_h__


#include <iostream>

#include <TROOT.h>
#include <TDirectory.h>
#include <TClonesArray.h>
#include <TH1D.h>
#include <TH1I.h>


#define GH1_MaxRandCuts 10
#define GH1_MaxBins 352

class   GTreeTagger;

class  GH1  : public TNamed
{
private:

protected:
    Int_t       nBins;

    TClonesArray   prompt;
    TClonesArray   rand;

    static  Double_t    cutPromptMin;
    static  Double_t    cutPromptMax;
    static  TArrayD     cutRandMin;
    static  TArrayD     cutRandMax;

    static  Double_t    backgroundSubstractionFactor;

    virtual Bool_t	AddPromptBin(const Int_t channel) = 0;
    virtual Bool_t	AddRandBin(const Int_t channel, const Int_t RandCut) = 0;
    virtual void	BackgroundSubtraction();
public:
    GH1(const TString& _Name, const TString& _Title, const Int_t NumberOfBins);
    GH1(const char* _Name, const char* _Title, const Int_t NumberOfBins);
    virtual ~GH1()    = 0;

    virtual void    Clear();
    virtual void    Fill(const Double_t value, const Double_t taggerTime, const Int_t taggerChannel = 0);
    virtual void    Fill(const Int_t value, const Double_t taggerTime, const Int_t taggerChannel = 0);
    static  void    InitCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t RandMin, const Double_t RandMax);
    static  void    AddRandCut(const Double_t RandMin, const Double_t RandMax);
    virtual void    Write(TDirectory& dir);

    friend class PPhysics;
};






class  GH1D   : public GH1
{
private:
    Double_t    minBin;
    Double_t    maxBin;

protected:
    virtual Bool_t	AddPromptBin(const Int_t channel);
    virtual Bool_t	AddRandBin(const Int_t channel, const Int_t RandCut);

public:
    GH1D(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max);
    GH1D(const char* name, const char* title, const Int_t nBins, const Double_t min, const Double_t max);
    virtual ~GH1D();
};






class  GH1I   : public GH1
{
private:
    Int_t    minBin;
    Int_t    maxBin;

protected:
    virtual Bool_t	AddPromptBin(const Int_t channel);
    virtual Bool_t	AddRandBin(const Int_t channel, const Int_t RandCut);

public:
    GH1I(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max);
    GH1I(const char* name, const char* title, const Int_t nBins, const Int_t min, const Int_t max);
    virtual ~GH1I();
};




#endif
