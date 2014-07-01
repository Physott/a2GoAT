#ifndef __PHist_h__
#define __PHist_h__


#include <TROOT.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH1I.h>

#define PHist_maxBins   48


class PPhysics;

class  PHist
{
private:
    void	BackgroundSubtraction();

protected:
    TH1*   prompt;
    TH1*   rand[2];
    TH1*   result;
    TH1*   prompt_bins[PHist_maxBins];
    TH1*   rand_bins[PHist_maxBins][2];
    TH1*   result_bins[PHist_maxBins];

    static  Double_t    cuts[3][2];

public:
    PHist();
    virtual ~PHist()    = 0;

            void    Clear();
    inline  void    Fill(const Double_t value, const Double_t taggerTime, const Int_t taggerChannel);
    inline  void    Fill(const Int_t value, const Double_t taggerTime, const Int_t taggerChannel);
    static  void    SetCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t Rand0Min = 0, const Double_t Rand0Max = 0, const Double_t Rand1Min = 0, const Double_t Rand1Max = 0);
            void    Write(TDirectory& dir);

    friend class PPhysics;
};






class  PHistD   : public PHist
{
private:

protected:

public:
    PHistD(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max);
    virtual ~PHistD();
};






class  PHistI   : public PHist
{
private:

protected:

public:
    PHistI(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max);
    virtual ~PHistI();
};






void    PHist::Fill(const Double_t value, const Double_t taggerTime, const Int_t taggerChannel)
{
    if(taggerTime>cuts[0][0] && taggerTime<cuts[0][1])
    {
        prompt->Fill(value);
        prompt_bins[taggerChannel]->Fill(value);
    }
    if(taggerTime>cuts[1][0] && taggerTime<cuts[1][1])
    {
        rand[0]->Fill(value);
        rand_bins[taggerChannel][0]->Fill(value);
    }
    if(taggerTime>cuts[2][0] && taggerTime<cuts[2][1])
    {
        rand[1]->Fill(value);
        rand_bins[taggerChannel][1]->Fill(value);
    }
}

void    PHist::Fill(const Int_t value, const Double_t taggerTime, const Int_t taggerChannel)
{
    if(taggerTime>cuts[0][0] && taggerTime<cuts[0][1])
    {
        prompt->Fill(value);
        prompt_bins[taggerChannel]->Fill(value);
    }
    if(taggerTime>cuts[1][0] && taggerTime<cuts[1][1])
    {
        rand[0]->Fill(value);
        rand_bins[taggerChannel][0]->Fill(value);
    }
    if(taggerTime>cuts[2][0] && taggerTime<cuts[2][1])
    {
        rand[1]->Fill(value);
        rand_bins[taggerChannel][1]->Fill(value);
    }
}



#endif
