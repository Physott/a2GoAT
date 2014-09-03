#ifndef __GHistThetaBinning_h__
#define __GHistThetaBinning_h__


#include <TROOT.h>
#include <TDirectory.h>

#include "GHistScaCor.h"


class   GTreeTagger;

class  GHistThetaBinning  : public GHistScaCor
{
private:
    TObjArray   bin;

    static  Int_t       ThetaBinningCount;
    static  Double_t    ThetaBinningWidth;
    static  Double_t    ThetaBinningRangeMin;
    static  Double_t    ThetaBinningRangeMax;

    static  Int_t   GetThetaBin(const Double_t theta);
            void    CreateBin();
            void    ExpandBin(const Int_t newSize);
    static  void    WriteHistogram(GHistLinked *hist, const char* name, const char* title, TDirectory* dir = 0);

public:
    GHistThetaBinning();
    GHistThetaBinning(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t linkHistogram = kTRUE);
    virtual ~GHistThetaBinning();

    virtual Bool_t          Add(const GHistThetaBinning* h, Double_t c = 1);
    virtual void            CalcResult();
    const   GHistScaCor&    GetThetaBin(const Int_t channel)    const   {return *((GHistScaCor*)bin.At(channel));}
    const   GHistScaCor&    GetSum()                             const   {return *((GHistScaCor*)this);}
    static  void            InitTaggerBinning(const Int_t count, const Double_t min, const Double_t max);
    virtual Int_t           Fill(const Double_t value, const Double_t theta = -1);
    virtual void            PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void            Reset(Option_t* option = "");
    virtual void            Scale(Double_t c1 = 1, Option_t* option = "");
    virtual void            ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t           WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
};




#endif
