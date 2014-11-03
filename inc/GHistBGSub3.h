#ifndef __GHistBGSub3_h__
#define __GHistBGSub3_h__


#include <TDirectory.h>

#include "GHistBGSub.h"
#include "GHistScaCor3.h"


class   GTreeTagger;

class  GHistBGSub3  : public GHistBGSub  //prompt
{
private:

protected:
    GHistBGSub3(Bool_t linkHistogram);

    virtual void    CreateRandBin();

public:
    GHistBGSub3();
    GHistBGSub3(const char* name, const char* title,
                Int_t nbinsx, Double_t xlow, Double_t xup,
                const Int_t nbinsy, const Double_t ylow, const Double_t yup,
                const Int_t nbinsz, const Double_t zlow, const Double_t zup,
                Bool_t linkHistogram = kTRUE);
    virtual ~GHistBGSub3();

    virtual Int_t   Fill(const Double_t x);
    virtual Int_t   Fill(const Double_t x, const Double_t y);
    virtual Int_t   Fill(const Double_t x, const Double_t y, const Double_t z)          {return ((GHistScaCor2*)result)->Fill(x, y);}
    virtual Int_t   Fill(const Double_t x, const Double_t y, const Double_t z, const Double_t taggerTime);
    virtual Int_t   Fill(const Double_t x, const Double_t y, const Double_t z, const GTreeTagger& tagger);
};




#endif
