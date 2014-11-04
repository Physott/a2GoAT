#ifndef __GHistBGSub2_h__
#define __GHistBGSub2_h__


#include <TDirectory.h>

#include "GHistBGSub.h"
#include "GHistScaCor2.h"


class   GTreeTagger;

class  GHistBGSub2  : public GHistBGSub  //prompt
{
private:

protected:
    GHistBGSub2(Bool_t linkHistogram);

    virtual void    CreateRandBin();

public:
    GHistBGSub2();
    GHistBGSub2(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, const Int_t nbinsy, const Double_t ylow, const Double_t yup, Bool_t linkHistogram = kTRUE);
    virtual ~GHistBGSub2();

    virtual Int_t   Fill(const Double_t x);
    virtual Int_t   Fill(const Double_t x, const Double_t y)                            {return ((GHistScaCor2*)result)->Fill(x, y);}
    virtual Int_t   Fill(const Double_t x, const Double_t y, const Double_t taggerTime);
    virtual Int_t   Fill(const Double_t x, const Double_t y, const GTreeTagger& tagger);

    virtual GHistBGSub*    ProjectionX(const char* name = "_px", Int_t firstybin = 0, Int_t lastybin = -1, Option_t* option = "");
    virtual GHistBGSub*    ProjectionY(const char* name = "_px", Int_t firstxbin = 0, Int_t lastxbin = -1, Option_t* option = "");
};




#endif
