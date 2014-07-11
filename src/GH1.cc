#include "GH1.h"



GH1::GH1()
{

}

GH1::~GH1()
{
    if(prompt)  delete  prompt;
    if(rand[0]) delete  rand[0];
    if(rand[1]) delete  rand[1];
    if(result)  delete  result;
}

void    GH1::Clear()
{
    prompt->Reset();
    rand[0]->Reset();
    rand[1]->Reset();
    result->Reset();
}

void    GH1::Fill(const Double_t value, const Double_t taggerTime)
{
    if(taggerTime>cuts[0][0] && taggerTime<cuts[0][1])
    {
        prompt->Fill(value);
    }
    if(taggerTime>cuts[1][0] && taggerTime<cuts[1][1])
    {
        rand[0]->Fill(value);
    }
    if(taggerTime>cuts[2][0] && taggerTime<cuts[2][1])
    {
        rand[1]->Fill(value);
    }
}

void    GH1::Fill(const Int_t value, const Double_t taggerTime)
{
    if(taggerTime>cuts[0][0] && taggerTime<cuts[0][1])
    {
        prompt->Fill(value);
    }
    if(taggerTime>cuts[1][0] && taggerTime<cuts[1][1])
    {
        rand[0]->Fill(value);
    }
    if(taggerTime>cuts[2][0] && taggerTime<cuts[2][1])
    {
        rand[1]->Fill(value);
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
    prompt->Write();


    curDir  = dir.GetDirectory("rand");
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir("rand");
        curDir  = dir.GetDirectory("rand");
    }
    curDir->cd();
    rand[0]->Write();
    rand[1]->Write();

    dir.cd();
    result->Write();
}


void	GH1::BackgroundSubtraction()
{
    result->Reset();
    result->Add(prompt,1);
    result->Add(rand[0],-backgroundSubstractionFactor);
    result->Add(rand[1],-backgroundSubstractionFactor);
}

Double_t    GH1::cuts[3][2]  =
{{-1000000, 1000000},
 {1, 0},
 {1, 0}
};

Double_t    GH1::backgroundSubstractionFactor  = 0;

void    GH1::SetCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t Rand0Min, const Double_t Rand0Max, const Double_t Rand1Min, const Double_t Rand1Max)
{
    cuts[0][0]  = PromptMin;
    cuts[0][1]  = PromptMax;
    cuts[1][0]  = Rand0Min;
    cuts[1][1]  = Rand0Max;
    cuts[2][0]  = Rand1Min;
    cuts[2][1]  = Rand1Max;
    backgroundSubstractionFactor = (PromptMax - PromptMin)/((Rand0Max - Rand0Min) + (Rand1Max - Rand1Min));
}







GH1D::GH1D(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max)
{
    gROOT->cd();
    prompt  = new TH1D((TString(name).Prepend("prompt_")).Data(), (TString(title).Append(" (prompt)")).Data(), nBins, min, max);
    if(!prompt)
        printf("#ERROR GH1D::GH1D: Can not create History prompt");

    rand[0]  = new TH1D((TString(name).Prepend(TString("rand0_"))).Data(), (TString(title).Append(TString("rand0"))).Data(), nBins, min, max);
    if(!rand[0])
            printf("#ERROR GH1D::GH1D: Can not create History prompt");
    rand[1]  = new TH1D((TString(name).Prepend(TString("rand1_"))).Data(), (TString(title).Append(TString("rand1"))).Data(), nBins, min, max);
    if(!rand[1])
            printf("#ERROR GH1D::GH1D: Can not create History prompt");

    result  = new TH1D(name, title, nBins, min, max);
    if(!result)
        printf("#ERROR GH1D::GH1D: Can not create History prompt");
}

GH1D::~GH1D()
{
}








GH1I::GH1I(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max)
{
    gROOT->cd();
    prompt  = new TH1I((TString(name).Prepend("prompt_")).Data(), (TString(title).Append(" (prompt)")).Data(), nBins, min, max);
    if(!prompt)
        printf("#ERROR GH1D::GH1D: Can not create History prompt");

    rand[0]  = new TH1I((TString(name).Prepend(TString("rand0_"))).Data(), (TString(title).Append(TString("rand0"))).Data(), nBins, min, max);
    if(!rand[0])
            printf("#ERROR GH1D::GH1D: Can not create History prompt");
    rand[1]  = new TH1I((TString(name).Prepend(TString("rand1_"))).Data(), (TString(title).Append(TString("rand1"))).Data(), nBins, min, max);
    if(!rand[1])
            printf("#ERROR GH1D::GH1D: Can not create History prompt");

    result  = new TH1I(name, title, nBins, min, max);
    if(!result)
        printf("#ERROR GH1D::GH1D: Can not create History prompt");
}

GH1I::~GH1I()
{
}
