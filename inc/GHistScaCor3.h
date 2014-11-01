#ifndef __GHistScaCor3_h__
#define __GHistScaCor3_h__


#include <TH3.h>

#include "GHistScaCor2.h"


class   GHistScaCor3 : public GHistScaCor2
{
private:
    void    CreateSingleScalerRead();

protected:
    GHistScaCor3(const Bool_t linkHistogram);

public:
    GHistScaCor3();
    GHistScaCor3(const char* name, const char* title,
                 const Int_t nbinsx, const Double_t xlow, const Double_t xup,
                 const Int_t nbinsy, const Double_t ylow, const Double_t yup,
                 const Int_t nbinsz, const Double_t zlow, const Double_t zup,
                 const Bool_t linkHistogram = kTRUE);
    virtual ~GHistScaCor3();

    virtual Int_t	Fill(Double_t x);
    virtual Int_t	Fill(Double_t x, Double_t y);
    virtual Int_t	Fill(Double_t x, Double_t y, Double_t z)    {return ((TH3*)buffer)->Fill(x, y, z);}
            Int_t   GetNbinsZ()                         const   {return accumulatedCorrected->GetNbinsZ();}
            Int_t   GetZmin()                           const   {return accumulatedCorrected->GetZaxis()->GetXmin();}
            Int_t   GetZmax()                           const   {return accumulatedCorrected->GetZaxis()->GetXmax();}
    virtual void	SetBins(Int_t nx, Double_t xmin, Double_t xmax,
                            Int_t ny, Double_t ymin, Double_t ymax,
                            Int_t nz, Double_t zmin, Double_t zmax);
};


#endif
