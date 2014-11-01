#ifndef __GHistBGSub_h__
#define __GHistBGSub_h__


#include <TDirectory.h>

#include "GHistScaCor.h"


class   GTreeTagger;

class  GHistBGSub  : public GHistScaCor  //prompt
{
private:
    TObjArray   rand;
    GHistScaCor randSum;
    GHistScaCor prompt;

    Bool_t      writeWindows;

    static  Double_t    cutPromptMin;
    static  Double_t    cutPromptMax;
    static  std::vector<Double_t> cutRandMin;
    static  std::vector<Double_t> cutRandMax;
    static  Double_t    backgroundSubstractionFactor;


    void    CreateRandBin();
    void    ExpandRandBins(const Int_t newSize);

public:
    GHistBGSub();
    GHistBGSub(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t linkHistogram = kTRUE);
    virtual ~GHistBGSub();

    virtual Bool_t	Add(const GHistBGSub* h, Double_t c = 1);
    static  void    AddRandCut(const Double_t RandMin, const Double_t RandMax);
    virtual void    CalcResult();
    virtual Int_t   Fill(const Double_t value)                                                                      {return GHistScaCor::Fill(value);}
    virtual Int_t   Fill(const Double_t value, const Double_t taggerTime);
    virtual Int_t   Fill(const Double_t value, const GTreeTagger& tagger);
    static  Int_t   GetNRandCuts()   {cutRandMin.size();}
    static  void    InitCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t RandMin, const Double_t RandMax);
    virtual Bool_t  IsEmpty();
    static  Bool_t  IsPrompt(const Double_t value);
    static  Bool_t  IsRandom(const Double_t value);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
            void    SetWriteWindows(const Bool_t value)                                                         {writeWindows = value;}
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
};




#endif
