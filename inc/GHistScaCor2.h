#ifndef __GHistScaCor2_h__
#define __GHistScaCor2_h__


#include <TH2.h>

#include "GHistScaCor.h"


class   GHistScaCor2 : public GHistScaCor
{
private:
    void    CreateSingleScalerRead();

protected:
    GHistScaCor2(const Bool_t linkHistogram);

public:
    GHistScaCor2();
    GHistScaCor2(const char* name, const char* title, const Int_t nbinsx, const Double_t xlow, const Double_t xup, const Int_t nbinsy, const Double_t ylow, const Double_t yup, const Bool_t linkHistogram = kTRUE);
    virtual ~GHistScaCor2();

    virtual Int_t	Fill(Double_t x);
    virtual Int_t	Fill(Double_t x, Double_t y)        {return ((TH2*)buffer)->Fill(x, y);}
            Int_t   GetNbinsY()                 const   {return accumulatedCorrected->GetNbinsY();}
            Int_t   GetYmin()                   const   {return accumulatedCorrected->GetYaxis()->GetXmin();}
            Int_t   GetYmax()                   const   {return accumulatedCorrected->GetYaxis()->GetXmax();}
    virtual void	SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax);

    virtual GHistScaCor*    ProjectionX(const char* name = "_px", Int_t firstybin = 0, Int_t lastybin = -1, Option_t* option = "");
    virtual GHistScaCor*    ProjectionY(const char* name = "_px", Int_t firstxbin = 0, Int_t lastxbin = -1, Option_t* option = "");
};


#endif
