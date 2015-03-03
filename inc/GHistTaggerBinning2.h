#ifndef __GHistTaggerBinning2_h__
#define __GHistTaggerBinning2_h__


#include <TROOT.h>
#include <TDirectory.h>

#include "GHistTaggerBinning.h"
#include "GHistBGSub3.h"


class   GTreeTagger;

class  GHistTaggerBinning2  : public GHistTaggerBinning
{
private:

protected:
    GHistTaggerBinning2(Bool_t linkHistogram);

public:
    GHistTaggerBinning2();
    GHistTaggerBinning2(const char* name, const char* title,
                       const Int_t nbinsx, const Double_t xlow, const Double_t xup,
                       const Int_t nbinsy, const Double_t ylow, const Double_t yup,
                       const Int_t nTaggerBins = 352,
                       const Bool_t linkHistogram = kTRUE);
    virtual ~GHistTaggerBinning2();

    virtual void            CalcResult();
    virtual Int_t           Fill(const Double_t value);
    virtual Int_t           Fill(const Double_t x, const Double_t y)                                                          {return ((GHistBGSub2*)sum)->Fill(x, y);}
    virtual Int_t           Fill(const Double_t x, const Double_t y, const Double_t taggerTime)                               {return ((GHistBGSub2*)sum)->Fill(x, y, taggerTime);}
    virtual Int_t           Fill(const Double_t x, const Double_t y, const Double_t taggerTime, const Int_t taggerChannel)    {return ((GHistBGSub3*)array)->Fill(x, y, taggerChannel, taggerTime);}
    virtual Int_t           Fill(const Double_t x, const Double_t y, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning = kFALSE);
};




#endif
