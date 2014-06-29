#ifndef __PHist_h__
#define __PHist_h__


#include <TROOT.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH1I.h>

class PPhysics;

class  PHist
{
private:

protected:
    TH1*   prompt;
    TH1*   rand[2];
    TH1*   result;

    Double_t    cuts[3][2];

public:
    PHist();
    virtual ~PHist()    = 0;

    //        void    Add(const PHist *hist, const Double_t scale);
    inline  void    Fill(const Double_t taggerTime, const Double_t value);
    inline  void    Fill(const Double_t taggerTime, const Int_t value);
            void	RandomSubtraction();
            void    Reset() {prompt->Reset(); rand[0]->Reset(); rand[1]->Reset(); result->Reset();}
    inline  void    SetCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t Rand0Min = 0, const Double_t Rand0Max = 0, const Double_t Rand1Min = 0, const Double_t Rand1Max = 0);
            void    Write(TDirectory *dir);

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






void    PHist::Fill(const Double_t taggerTime, const Double_t value)
{
    if(taggerTime>cuts[0][0] && taggerTime<cuts[0][1])
        prompt->Fill(value);
    if(taggerTime>cuts[1][0] && taggerTime<cuts[1][1])
        rand[0]->Fill(value);
    if(taggerTime>cuts[2][0] && taggerTime<cuts[2][1])
        rand[1]->Fill(value);
}

void    PHist::Fill(const Double_t taggerTime, const Int_t value)
{
    if(taggerTime>cuts[0][0] && taggerTime<cuts[0][1])
        prompt->Fill(value);
    if(taggerTime>cuts[1][0] && taggerTime<cuts[1][1])
        rand[0]->Fill(value);
    if(taggerTime>cuts[2][0] && taggerTime<cuts[2][1])
        rand[1]->Fill(value);
}

void    PHist::SetCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t Rand0Min, const Double_t Rand0Max, const Double_t Rand1Min, const Double_t Rand1Max)
{
    cuts[0][0]  = PromptMin;
    cuts[0][1]  = PromptMax;
    cuts[1][0]  = Rand0Min;
    cuts[1][1]  = Rand0Max;
    cuts[2][0]  = Rand1Min;
    cuts[2][1]  = Rand1Max;
}




#endif
