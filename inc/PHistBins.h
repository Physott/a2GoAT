#ifndef __PHistBins_h__
#define __PHistBins_h__


#include "GH1.h"

#define PHistBins_maxBins   48


class PPhysics;

class  PHistBins    : public GH1
{
private:
    virtual void    BackgroundSubtraction();

protected:
    TH1*   prompt_bins[PHistBins_maxBins];
    TH1*   rand_bins[PHistBins_maxBins][2];
    TH1*   result_bins[PHistBins_maxBins];

public:
    PHistBins();
    virtual ~PHistBins()    = 0;

    virtual void    Clear();
    virtual void    Fill(const Double_t value, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void    Fill(const Int_t value, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void    Fill(const Double_t value, const Double_t taggerTime)   {Fill(value, taggerTime, 0);}
    virtual void    Fill(const Int_t value, const Double_t taggerTime)      {Fill(value, taggerTime, 0);}
    static  void    SetCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t Rand0Min = 0, const Double_t Rand0Max = 0, const Double_t Rand1Min = 0, const Double_t Rand1Max = 0);
    virtual void    Write(TDirectory& dir);

    friend class PPhysics;
};






class  PHistBinsD   : public PHistBins
{
private:

protected:

public:
    PHistBinsD(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max);
    virtual ~PHistBinsD();
};






class  PHistBinsI   : public PHistBins
{
private:

protected:

public:
    PHistBinsI(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max);
    virtual ~PHistBinsI();
};



#endif
