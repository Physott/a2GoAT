#ifndef __GScaCorHist_h__
#define __GScaCorHist_h__


#include <TROOT.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH1I.h>

#include "GHistManager.h"


TDirectory* GGetDirectory(TDirectory* dir, const char* name);
TDirectory* GGetDirectory(TDirectory* dir, const TString& name);



class   GScaCorHist : public GLinkedHist
{
private:
    Int_t   nCorrected;

    static  void    WriteHistogram(TH1* hist, const TString& name, const TString& title);

protected:
    TH1*    current;                //pointer to base class of corresponding type
    TH1*    accumulated;
    TH1*    accumulatedCorrected;

public:
    GScaCorHist(const Bool_t linked = true);
    virtual ~GScaCorHist() = 0;

    virtual void    Clear(Option_t* option = "");
            Int_t   GetNScalerReadCorrections() const   {return nCorrected;}
            void    ScalerReadCorrection(const Double_t CorrectionFactor, TDirectory *dir = 0);
    virtual void	Write(TDirectory *dir);
};






class   GScaCorHist1D   : public GScaCorHist, public TH1D
{
private:

protected:

public:
    GScaCorHist1D(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max, const Bool_t linked = true);
    virtual ~GScaCorHist1D();

    virtual void	Clear(Option_t* option = "")    {TH1D::Clear(option); GScaCorHist::Clear(option);}
    //virtual Int_t	Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
    virtual void	Write(TDirectory* dir)            {if(GetNScalerReadCorrections()>0) return GScaCorHist::Write(dir); TH1D::Write();}
};



class   GScaCorHist1I   : public TH1I, public GScaCorHist
{
private:

protected:

public:
    GScaCorHist1I(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max, const Bool_t linked = true);
    virtual ~GScaCorHist1I();

    virtual void	Clear(Option_t* option = "")    {TH1I::Clear(option); GScaCorHist::Clear(option);}
    //virtual Int_t	Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
    virtual void	Write(TDirectory* dir)            {if(GetNScalerReadCorrections()>0) return GScaCorHist::Write(dir); TH1I::Write();}
};






#endif
