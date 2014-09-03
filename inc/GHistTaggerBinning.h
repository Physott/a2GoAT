#ifndef __GHistTaggerBinning_h__
#define __GHistTaggerBinning_h__


#include <TROOT.h>
#include <TDirectory.h>

#include "GHistThetaBinning.h"


class   GTreeTagger;

class  GHistTaggerBinning  : public GHistThetaBinning
{
private:
    TObjArray   bin;

    static  Int_t   TaggerBinningRangeMin;
    static  Int_t   TaggerBinningRangeMax;

            void    CreateBin();
            void    ExpandBin(const Int_t newSize);
    static  void    WriteHistogram(GHistLinked *hist, const char* name, const char* title, TDirectory* dir = 0);

protected:

public:
    GHistTaggerBinning();
    GHistTaggerBinning(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t linkHistogram = kTRUE);
    virtual ~GHistTaggerBinning();

    virtual Bool_t              Add(const GHistTaggerBinning* h, Double_t c = 1);
    virtual void                CalcResult();
    const   GHistThetaBinning&  GetTaggerBin(const Int_t channel)    const   {return *((GHistThetaBinning*)bin.At(channel));}
    const   GHistThetaBinning&  GetSum()                             const   {return *((GHistThetaBinning*)this);}
    static  void                InitTaggerBinning(const Int_t min, const Int_t max);
    virtual Int_t               Fill(const Double_t value, const Int_t taggerChannel = -1, const Double_t theta = -1);
    virtual Int_t               Fill(const Double_t value, const GTreeTagger& tagger, const Double_t theta = -1, const Bool_t CreateHistogramsForTaggerBinning = kFALSE, const Bool_t CreateHistogramsForThetaBinning = kFALSE);
    virtual void                PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void                Reset(Option_t* option = "");
    virtual void                Scale(Double_t c1 = 1, Option_t* option = "");
    virtual void                ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual Int_t               WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
};




#endif
