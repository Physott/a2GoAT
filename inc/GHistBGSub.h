#ifndef __GHistBGSub_h__
#define __GHistBGSub_h__


#include <TDirectory.h>

#include "GHistTaggerBinning.h"


class   GTreeTagger;

class  GHistBGSub  : public GHistTaggerBinning  //prompt
{
private:
    TObjArray           rand;
    GHistTaggerBinning  randSum;
    GHistTaggerBinning  prompt;

    static  Double_t    cutPromptMin;
    static  Double_t    cutPromptMax;
    static  std::vector<Double_t> cutRandMin;
    static  std::vector<Double_t> cutRandMax;
    static  Double_t    backgroundSubstractionFactor;

            void    CreateRandBin();
            void    ExpandRandBins(const Int_t newSize);
    static  void    WriteHistogram(GHistTaggerBinning *hist, const char* name, const char* title, TDirectory* dir = 0);

public:
    GHistBGSub();
    GHistBGSub(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t linkHistogram = kTRUE);
    virtual ~GHistBGSub();

    virtual Bool_t	Add(const GHistBGSub* h, Double_t c = 1);
    static  void    AddRandCut(const Double_t RandMin, const Double_t RandMax);
    virtual void    CalcResult();
    virtual Int_t   Fill(const Double_t value, const Double_t taggerTime = 0, const Int_t taggerChannel = -1, const Double_t theta = -1);
    virtual Int_t   Fill(const Double_t value, const GTreeTagger& tagger, const Double_t theta = -1, const Bool_t CreateHistogramsForTaggerBinning = kFALSE, const Bool_t CreateHistogramsForThetaBinning = kFALSE) {GHistTaggerBinning::Fill(value, tagger, theta, CreateHistogramsForTaggerBinning, CreateHistogramsForThetaBinning);}
    static  Int_t   GetNRandCuts()   {cutRandMin.size();}
    static  void    InitCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t RandMin, const Double_t RandMax);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
};




#endif
