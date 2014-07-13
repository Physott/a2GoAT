#ifndef __GScaCorHist_h__
#define __GScaCorHist_h__


#include <TROOT.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH1I.h>



TDirectory* GGetDirectory(TDirectory* dir, const char* name);
TDirectory* GGetDirectory(TDirectory* dir, const TString& name);




class   GScaCorHist
{
private:
    Int_t   nCorrected;

protected:
    TH1*    current;                //pointer to base class of corresponding type
    TH1*    accumulated;

public:
    GScaCorHist();
    virtual ~GScaCorHist() = 0;

    virtual void    Clear(Option_t* option = "");
    Int_t   GetNScalerReadCorrections() const   {return nCorrected;}
            void    ScalerReadCorrection(const Double_t CorrectionFactor);
    virtual Int_t	Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const;
    virtual Int_t	Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
};






class   GScaCorHist1D   : public TH1D, public GScaCorHist
{
private:

protected:

public:
    GScaCorHist1D(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max);
    virtual ~GScaCorHist1D();

    virtual void	Clear(Option_t* option = "")    {TH1D::Clear(option); GScaCorHist::Clear(option);}
    virtual Int_t	Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const  {if(GetNScalerReadCorrections()>0) return GScaCorHist::Write(name, option, bufsize); return TH1D::Write(name, option, bufsize);}
    virtual Int_t	Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)        {if(GetNScalerReadCorrections()>0) return GScaCorHist::Write(name, option, bufsize); return TH1D::Write(name, option, bufsize);}
};



class   GScaCorHist1I   : public TH1I, public GScaCorHist
{
private:

protected:

public:
    GScaCorHist1I(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max);
    virtual ~GScaCorHist1I();

    virtual void	Clear(Option_t* option = "")    {TH1I::Clear(option); GScaCorHist::Clear(option);}
    virtual Int_t	Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const  {if(GetNScalerReadCorrections()>0) return GScaCorHist::Write(name, option, bufsize); return TH1I::Write(name, option, bufsize);}
    virtual Int_t	Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)        {if(GetNScalerReadCorrections()>0) return GScaCorHist::Write(name, option, bufsize); return TH1I::Write(name, option, bufsize);}
};






#endif
