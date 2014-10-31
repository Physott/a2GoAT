#include "GHistScaCor3.h"

#include <iostream>


GHistScaCor3::GHistScaCor3(const Bool_t linkHistogram) :
    GHistScaCor2(linkHistogram)
{
}

GHistScaCor3::GHistScaCor3() :
    GHistScaCor2(kTRUE)
{
    buffer                  = new TH3D();
    accumulated             = new TH3D();
    accumulatedCorrected    = new TH3D();
}

GHistScaCor3::GHistScaCor3(const char* name, const char* title, const Int_t nbinsx, const Double_t xlow, const Double_t xup, const Int_t nbinsy, const Double_t ylow, const Double_t yup, const Int_t nbinsz, const Double_t zlow, const Double_t zup, const Bool_t linkHistogram) :
    GHistScaCor2(linkHistogram)
{
    buffer                  = new TH3D(TString(name).Append(GHSC_bufferNameSuffix).Data(),
                                       TString(title).Append(GHSC_bufferTitleSuffix).Data(),
                                       nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
    accumulated             = new TH3D(TString(name).Append(GHSC_accumulatedNameSuffix).Data(),
                                       TString(title).Append(GHSC_accumulatedTitleSuffix).Data(),
                                       nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
    accumulatedCorrected    = new TH3D(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
}

GHistScaCor3::~GHistScaCor3()
{
}

void    GHistScaCor3::CreateSingleScalerRead()
{
    TH3D*    uncor   = new TH3D(TString(accumulatedCorrected->GetName()).Append(GHSC_singleScalerReadNameSuffix).Append(TString::Itoa(singleScalerReads.GetEntriesFast(), 10)),
                                TString(accumulatedCorrected->GetTitle()).Append(GHSC_singleScalerReadTitleSuffix).Append(TString::Itoa(singleScalerReads.GetEntriesFast(), 10)),
                                accumulatedCorrected->GetNbinsX(),
                                accumulatedCorrected->GetXaxis()->GetXmin(),
                                accumulatedCorrected->GetXaxis()->GetXmax(),
                                accumulatedCorrected->GetNbinsY(),
                                accumulatedCorrected->GetYaxis()->GetXmin(),
                                accumulatedCorrected->GetYaxis()->GetXmax(),
                                accumulatedCorrected->GetNbinsZ(),
                                accumulatedCorrected->GetZaxis()->GetXmin(),
                                accumulatedCorrected->GetZaxis()->GetXmax());
    uncor->SetDirectory(0);
    singleScalerReads.Add(uncor);

    TH3D*    cor     = new TH3D(TString(accumulatedCorrected->GetName()).Append(GHSC_singleScalerReadCorrectedNameSuffix).Append(TString::Itoa(singleScalerReadsCorrected.GetEntriesFast(), 10)),
                                TString(accumulatedCorrected->GetTitle()).Append(GHSC_singleScalerReadCorrectedTitleSuffix).Append(TString::Itoa(singleScalerReadsCorrected.GetEntriesFast(), 10)),
                                accumulatedCorrected->GetNbinsX(),
                                accumulatedCorrected->GetXaxis()->GetXmin(),
                                accumulatedCorrected->GetXaxis()->GetXmax(),
                                accumulatedCorrected->GetNbinsY(),
                                accumulatedCorrected->GetYaxis()->GetXmin(),
                                accumulatedCorrected->GetYaxis()->GetXmax(),
                                accumulatedCorrected->GetNbinsZ(),
                                accumulatedCorrected->GetZaxis()->GetXmin(),
                                accumulatedCorrected->GetZaxis()->GetXmax());
    cor->SetDirectory(0);
    singleScalerReadsCorrected.Add(cor);
}

Int_t	GHistScaCor3::Fill(Double_t x)
{
    std::cout << "WARNING: You tried to fill a 2 dim. Hist with only 1 value." << std::endl;
    return 0;
}

void	GHistScaCor3::SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax, Int_t nz, Double_t zmin, Double_t zmax)
{
    buffer->SetBins(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
    accumulated->SetBins(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
    accumulatedCorrected->SetBins(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
    for(int i=0; i<singleScalerReads.GetEntriesFast(); i++)
    {
        ((TH1D*)singleScalerReads.At(i))->SetBins(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
        ((TH1D*)singleScalerReadsCorrected.At(i))->SetBins(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax);
    }
}

