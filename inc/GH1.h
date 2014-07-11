#ifndef __GH1_h__
#define __GH1_h__


#include <TROOT.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH1I.h>

class PPhysics;

class  GH1
{
private:
    virtual void	BackgroundSubtraction();

protected:
    TH1*   prompt;
    TH1*   rand[2];
    TH1*   result;

    static  Double_t    cuts[3][2];
    static  Double_t    backgroundSubstractionFactor;

public:
    GH1();
    virtual ~GH1()    = 0;

    virtual void    Clear();
    virtual void    Fill(const Double_t value, const Double_t taggerTime);
    virtual void    Fill(const Int_t value, const Double_t taggerTime);
    static  void    SetCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t Rand0Min = 0, const Double_t Rand0Max = 0, const Double_t Rand1Min = 0, const Double_t Rand1Max = 0);
    virtual void    Write(TDirectory& dir);

    friend class PPhysics;
};






class  GH1D   : public GH1
{
private:

protected:

public:
    GH1D(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max);
    virtual ~GH1D();
};






class  GH1I   : public GH1
{
private:

protected:

public:
    GH1I(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max);
    virtual ~GH1I();
};




#endif
