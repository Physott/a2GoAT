#include "GHistManager.h"



GHistManager*   gGHistManager = 0;


GHistManager::GHistManager()    :
    histList()
{
    gGHistManager   = this;
    histList.SetOwner();
}

GHistManager::~GHistManager()
{
    gGHistManager   = 0;
}

void    GHistManager::AddHistogramToList(GLinkedHist* hist)
{
    histList.AddAtFree(hist);
}

void    GHistManager::Clear()
{
    histList.Clear();
}

void    GHistManager::RemoveHistogramFromList(GLinkedHist* hist)
{
    histList.Remove(hist);
}

Bool_t  GHistManager::WriteLinkedHistograms(TDirectory* dir)
{
    TIter   iter(&histList);
    GLinkedHist*    hist;
    while(hist=(GLinkedHist*)iter.Next())
        hist->Write(dir);
}







GLinkedHist::GLinkedHist(const Bool_t linked)
{
    if(!linked)
        return;
    if(gGHistManager)
        gGHistManager->AddHistogramToList(this);
}

GLinkedHist::~GLinkedHist()
{
    if(gGHistManager)
        gGHistManager->RemoveHistogramFromList(this);
}







/*
GLinkedHist1D::GLinkedHist1D(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max)  :
    TH1D(name, title, nBins, min, max)
{

}

GLinkedHist1D::~GLinkedHist1D()
{

}







GLinkedHist1I::GLinkedHist1I(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max)  :
    TH1I(name, title, nBins, min, max)
{

}

GLinkedHist1I::~GLinkedHist1I()
{

}
*/
