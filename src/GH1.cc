#include "GH1.h"
#include "GTreeTagger.h"


Double_t    GH1::cutPromptMin  = -1000000;
Double_t    GH1::cutPromptMax  =  1000000;
TArrayD     GH1::cutRandMin;
TArrayD     GH1::cutRandMax;

GH1::GH1(const TString& _Name, const TString& _Title, const Int_t NumberOfBins) :
    TNamed(_Name, _Title),
    nBins(NumberOfBins),
    rand("TClonesArray", GH1_MaxRandCuts)
{

}

GH1::GH1(const char* _Name, const char* _Title, const Int_t NumberOfBins) :
    TNamed(_Name, _Title),
    nBins(NumberOfBins)
{

}

GH1::~GH1()
{
}

void    GH1::Clear()
{
    prompt.Clear("C");
    rand.Clear("C");
    //rand[0]->Reset();
    //rand[1]->Reset();
}

void    GH1::Fill(const Double_t value, const Double_t taggerTime, const Int_t taggerChannel)
{
    if(taggerTime>=cutPromptMin && taggerTime<=cutPromptMax)
    {
        if(taggerChannel>=prompt.GetSize())
            prompt.Expand(taggerChannel+1);
        if(!prompt.UncheckedAt(taggerChannel))
            AddPromptBin(taggerChannel);
        ((TH1*)prompt[taggerChannel])->Fill(value);
    }

    for(int i=0; i<cutRandMin.GetSize(); i++)
    {
        if(taggerTime>=cutRandMin[i] && taggerTime<=cutRandMax[i])
        {
            if(i>=rand.GetSize())
                rand.Expand(i+1);
            if(!rand.UncheckedAt(i))
                AddRandBin(taggerChannel, i);
            if(taggerChannel>=((TClonesArray*)rand[i])->GetSize())
                ((TClonesArray*)rand[i])->Expand(taggerChannel+1);
            if(!(((TClonesArray*)rand[i])->UncheckedAt(i)))
                AddRandBin(taggerChannel, i);
            ((TH1*)prompt[taggerChannel])->Fill(value);
        }
    }
}

void    GH1::Fill(const Int_t value, const Double_t taggerTime, const Int_t taggerChannel)
{
    if(taggerTime>cutPromptMin && taggerTime<cutPromptMax)
    {
        if(taggerChannel>=prompt.GetSize())
            prompt.Expand(taggerChannel+1);
        if(!prompt.UncheckedAt(taggerChannel))
            AddPromptBin(taggerChannel);
        ((TH1*)prompt[taggerChannel])->Fill(value);
    }
}

void    GH1::Write(TDirectory& dir)
{
    BackgroundSubtraction();

    TDirectory* curDir  = dir.GetDirectory("prompt");
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir("prompt");
        curDir  = dir.GetDirectory("prompt");
    }
    curDir->cd();
    prompt.Write();


    curDir  = dir.GetDirectory("rand");
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir("rand");
        curDir  = dir.GetDirectory("rand");
    }
    curDir->cd();
    rand.Write();

    dir.cd();
    //result->Write();
}


void	GH1::BackgroundSubtraction()
{
    /*result->Reset();
    result->Add(prompt,1);
    result->Add(rand[0],-backgroundSubstractionFactor);
    result->Add(rand[1],-backgroundSubstractionFactor);*/
}

Double_t    GH1::backgroundSubstractionFactor  = 0;

void    GH1::InitCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t RandMin, const Double_t RandMax)
{
    cutPromptMin    = PromptMin;
    cutPromptMax    = PromptMax;
    cutRandMin.Set(1, &RandMin);
    cutRandMax.Set(1, &RandMax);
    backgroundSubstractionFactor = (PromptMax - PromptMin)/(RandMax - RandMin);
}

void    GH1::AddRandCut(const Double_t RandMin, const Double_t RandMax)
{
    cutRandMin.AddAt(RandMin, cutRandMin.GetSize());
    cutRandMax.AddAt(RandMax, cutRandMax.GetSize());
    backgroundSubstractionFactor = cutRandMax[0] - cutRandMin[0];
    for(int i=1; i<cutRandMin.GetSize(); i++)
        backgroundSubstractionFactor += cutRandMax[i] - cutRandMin[i];
    backgroundSubstractionFactor    = (cutPromptMax - cutPromptMin)/backgroundSubstractionFactor;
}









GH1D::GH1D(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max) :
    GH1(name, title, nBins),
    minBin(min),
    maxBin(max)
{
    prompt.SetClass("TH1D", GH1_MaxBins);
}
GH1D::GH1D(const char* name, const char* title, const Int_t nBins, const Double_t min, const Double_t max) :
    GH1(name, title, nBins),
    minBin(min),
    maxBin(max)
{
    prompt.SetClass("TH1D", GH1_MaxBins);
}

GH1D::~GH1D()
{
}

Bool_t  GH1D::AddPromptBin(const Int_t channel)
{
    gROOT->cd();
    new(prompt[channel]) TH1D((TString(GetName()).Prepend("prompt_")).Append(TString::Itoa(channel, 10).Prepend("_Bin")).Data(), (TString(GetTitle()).Append(TString::Itoa(channel, 10).Prepend(" channel ")).Append(" (prompt)")).Data(), nBins, minBin, maxBin);
}

Bool_t  GH1D::AddRandBin(const Int_t channel, const Int_t RandCut)
{
    gROOT->cd();
    if(rand[RandCut])
        new(rand[RandCut]) TClonesArray("TH1D", GH1_MaxBins);
    new((*((TClonesArray*)rand[RandCut]))[channel]) TH1D((TString(GetName()).Prepend((TString::Itoa(RandCut,10).Prepend("rand")).Append("_"))).Append(TString::Itoa(channel, 10).Prepend("_Bin")).Data(), TString(GetTitle()).Append(TString::Itoa(channel, 10).Prepend(" channel ")).Append(TString::Itoa(RandCut,10).Prepend(" (rand")).Append(")").Data(), nBins, minBin, maxBin);
}






GH1I::GH1I(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max) :
    GH1(name, title, nBins),
    minBin(min),
    maxBin(max)
{
    prompt.SetClass("TH1I", GH1_MaxBins);
}
GH1I::GH1I(const char* name, const char* title, const Int_t nBins, const Int_t min, const Int_t max) :
    GH1(name, title, nBins),
    minBin(min),
    maxBin(max)
{
    prompt.SetClass("TH1I", GH1_MaxBins);
}

GH1I::~GH1I()
{
}

Bool_t  GH1I::AddPromptBin(const Int_t channel)
{
    gROOT->cd();
    new(prompt[channel]) TH1I((TString(GetName()).Prepend("prompt_")).Append(TString::Itoa(channel, 10).Prepend("_Bin")).Data(), (TString(GetTitle()).Append(TString::Itoa(channel, 10).Prepend(" channel ")).Append(" (prompt)")).Data(), nBins, minBin, maxBin);
}

Bool_t  GH1I::AddRandBin(const Int_t channel, const Int_t RandCut)
{
    gROOT->cd();
    if(rand[RandCut])
        new(rand[RandCut]) TClonesArray("TH1I", GH1_MaxBins);
    new((*((TClonesArray*)rand[RandCut]))[channel]) TH1I((TString(GetName()).Prepend((TString::Itoa(RandCut,10).Prepend("rand")).Append("_"))).Append(TString::Itoa(channel, 10).Prepend("_Bin")).Data(), TString(GetTitle()).Append(TString::Itoa(channel, 10).Prepend(" channel ")).Append(TString::Itoa(RandCut,10).Prepend(" (rand")).Append(")").Data(), nBins, minBin, maxBin);
}
