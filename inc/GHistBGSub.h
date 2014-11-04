#ifndef __GHistBGSub_h__
#define __GHistBGSub_h__


#define GHBS_folderName         "BackgroundSubstraction"
#define GHBS_randFolderName     "RandomWindows"
#define GHBS_randNameSuffix     "_Rand"
#define GHBS_randTitleSuffix    " Rand "
#define GHBS_randSumNameSuffix  "_RandSum"
#define GHBS_randSumTitleSuffix " RandSum"
#define GHBS_promptNameSuffix   "_Prompt"
#define GHBS_promptTitleSuffix  " Prompt"




#include <TDirectory.h>

#include "GHistScaCor.h"


class   GTreeTagger;

class  GHistBGSub  : public GHistLinked
{
private:
    Bool_t      writeWindows;

protected:
    GHistScaCor*    result;
    GHistScaCor*    prompt;
    GHistScaCor*    randSum;
    TObjArray       rand;

    static  Double_t    cutPromptMin;
    static  Double_t    cutPromptMax;
    static  std::vector<Double_t> cutRandMin;
    static  std::vector<Double_t> cutRandMax;
    static  Double_t    backgroundSubstractionFactor;


    GHistBGSub(Bool_t linkHistogram);

    virtual void    CreateRandBin();
            void    ExpandRandBins(const Int_t newSize);

public:
    GHistBGSub();
    GHistBGSub(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t linkHistogram = kTRUE);
    virtual ~GHistBGSub();

    virtual Bool_t	Add(const GHistBGSub* h, Double_t c = 1);
    virtual Bool_t	Add(const GHistScaCor* _result, const GHistScaCor* _prompt, const GHistScaCor* _randSum, const TObjArray& _rand, const Double_t c = 1);
    static  void    AddRandCut(const Double_t RandMin, const Double_t RandMax);
    virtual void    CalcResult();
    virtual Int_t   Fill(const Double_t value)                                                                      {return result->Fill(value);}
    virtual Int_t   Fill(const Double_t value, const Double_t taggerTime);
    virtual Int_t   Fill(const Double_t value, const GTreeTagger& tagger);
    static  Int_t   GetNRandCuts()   {cutRandMin.size();}
    static  void    InitCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t RandMin, const Double_t RandMax);
    virtual Bool_t  IsEmpty();
    static  Bool_t  IsPrompt(const Double_t value);
    static  Bool_t  IsRandom(const Double_t value);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual void	Scale(Double_t c1 = 1, Option_t* option = "");
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
            void    SetWriteWindows(const Bool_t value)                                                         {writeWindows = value;}
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
};




#endif
