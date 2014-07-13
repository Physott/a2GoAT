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




void    GScaCorHist::WriteHistogram(TH1* hist, const TString& name, const TString& title)
{
    TString oldName(hist->GetName());
    TString oldTitle(hist->GetTitle());
    hist->SetNameTitle(name, title);
    hist->Write(0, TObject::kOverwrite);
    hist->SetNameTitle(oldName, oldTitle);
}



GScaCorHist::GScaCorHist() :
    nCorrected(0),
    current(0),
    accumulated(0),
    accumulatedCorrected(0)
{

}

GScaCorHist::~GScaCorHist()
{
    if(accumulated)
        delete accumulated;
    if(accumulatedCorrected)
        delete accumulatedCorrected;
}

void    GScaCorHist::Clear(Option_t* option)
{
    nCorrected  = 0;
    accumulated->Clear(option);
    accumulatedCorrected->Clear(option);
}

void    GScaCorHist::ScalerReadCorrection(const Double_t CorrectionFactor, TDirectory* dir)
{
    accumulated->Add(current);
    if(dir)
    {
        GGetDirectory(GGetDirectory(dir, "ScalerCorrection"), TString::Itoa(nCorrected, 10).Prepend("ScalerRead_"));
        WriteHistogram(current, TString(current->GetName()).Append("_Uncorrected_ScaRead").Append(TString::Itoa(nCorrected, 10)), TString(current->GetTitle()).Append("ScalerRead ").Append(TString::Itoa(nCorrected, 10)).Append(" (Uncorrected)"));
        current->Scale(CorrectionFactor);
        WriteHistogram(current, TString(current->GetName()).Append("_ScaRead").Append(TString::Itoa(nCorrected, 10)), TString(current->GetTitle()).Append("ScalerRead ").Append(TString::Itoa(nCorrected, 10)).Append(" (Corrected)"));
    }
    else
        current->Scale(CorrectionFactor);

    accumulatedCorrected->Add(current);
    current->Reset();
    nCorrected++;
}

void	GScaCorHist::Write(TDirectory* dir)
{
    if(dir)
    {
        GGetDirectory(dir, "ScalerCorrection");
        accumulated->Write(0, TObject::kOverwrite);
        dir->cd();
        WriteHistogram(accumulatedCorrected, TString(accumulatedCorrected->GetName()).Remove(strlen(accumulatedCorrected->GetName())-10), TString(accumulatedCorrected->GetTitle()).Remove(strlen(accumulatedCorrected->GetTitle())-12));
        return;
    }

    WriteHistogram(accumulatedCorrected, TString(accumulatedCorrected->GetName()).Remove(strlen(accumulatedCorrected->GetName())-10), TString(accumulatedCorrected->GetTitle()).Remove(strlen(accumulatedCorrected->GetTitle())-12));
    GGetDirectory(gDirectory, "ScalerCorrection");
    accumulated->Write(0, TObject::kOverwrite);
}










GScaCorHist1D::GScaCorHist1D(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max) :
    GScaCorHist(),
    TH1D(name, title, nBins, min, max)
{
    current = this;
    gROOT->cd();
    accumulated = new TH1D(TString(name).Append("_Uncorrected").Data(), TString(title).Append(" (Uncorrected)").Data(), nBins, min, max);
    accumulatedCorrected = new TH1D(TString(name).Append("_Corrected").Data(), TString(title).Append(" (Corrected)").Data(), nBins, min, max);
}

GScaCorHist1D::~GScaCorHist1D()
{

}
/*
Int_t	GScaCorHist1D::Write(const char* name, Int_t option, Int_t bufsize)
{
    if(GetNScalerReadCorrections()>0)
    {
        GScaCorHist::Write(gDirectory);
        return 0;
    }
    return TH1D::Write(name, option, bufsize);
}
*/



GScaCorHist1I::GScaCorHist1I(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max) :
    GScaCorHist(),
    TH1I(name, title, nBins, min, max)
{
    current = this;
    gROOT->cd();
    accumulated = new TH1I(TString(name).Append("_Uncorrected").Data(), TString(title).Append(" (Uncorrected)").Data(), nBins, min, max);
    accumulatedCorrected = new TH1I(TString(name).Append("_Corrected").Data(), TString(title).Append(" (Corrected)").Data(), nBins, min, max);
}

GScaCorHist1I::~GScaCorHist1I()
{

}
/*
Int_t	GScaCorHist1I::Write(const char* name, Int_t option, Int_t bufsize)
{
    if(GetNScalerReadCorrections()>0)
    {
        GScaCorHist::Write(gDirectory);
        return 0;
    }
    return TH1I::Write(name, option, bufsize);
}
*/





