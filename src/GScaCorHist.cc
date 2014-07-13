#include "GScaCorHist.h"



TDirectory* GGetDirectory(TDirectory* dir, const char* name)
{
    if(!dir)
        return 0;
    TDirectory* curDir  = dir->GetDirectory(name);
    if(!curDir)
    {
        dir->cd();
        gDirectory->mkdir(name);
        curDir  = dir->GetDirectory(name);
    }
    curDir->cd();
    return curDir;
}
TDirectory* GGetDirectory(TDirectory* dir, const TString& name)
{
    GGetDirectory(dir, name.Data());
}






GScaCorHist::GScaCorHist() :
    nCorrected(0),
    current(0),
    accumulated(0)
{

}

GScaCorHist::~GScaCorHist()
{
    if(accumulated)
        delete accumulated;
}

void    GScaCorHist::Clear(Option_t* option)
{
    nCorrected  = 0;
    accumulated->Clear(option);
}

void    GScaCorHist::ScalerReadCorrection(const Double_t CorrectionFactor)
{
    current->Scale(CorrectionFactor);
    accumulated->Add(current);
    current->Reset();
    nCorrected++;
}

Int_t	GScaCorHist::Write(const char* name, Int_t option, Int_t bufsize) const
{

}

Int_t	GScaCorHist::Write(const char* name, Int_t option, Int_t bufsize)
{

}








GScaCorHist1D::GScaCorHist1D(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max) :
    GScaCorHist(),
    TH1D(TString(name).Append("_Current"), TString(title).Append(" (Current)"), nBins, min, max)
{
    current = this;
    gROOT->cd();
    accumulated = new TH1D(name, title, nBins, min, max);
}

GScaCorHist1D::~GScaCorHist1D()
{

}




GScaCorHist1I::GScaCorHist1I(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max) :
    GScaCorHist(),
    TH1I(TString(name).Append("_Current"), TString(title).Append(" (Current)"), nBins, min, max)
{
    gROOT->cd();
    current = this;
    accumulated = new TH1I(name, title, nBins, min, max);
}

GScaCorHist1I::~GScaCorHist1I()
{

}






