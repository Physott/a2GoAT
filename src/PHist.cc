#include "PHist.h"



PHist::PHist()
{

}

PHist::~PHist()
{
    if(prompt)  delete  prompt;
    if(rand[0]) delete  rand[0];
    if(rand[1]) delete  rand[1];
    if(result)  delete  result;
}

void    PHist::Clear()
{
    prompt->Reset();
    rand[0]->Reset();
    rand[1]->Reset();
    result->Reset();
    for(int i=0; i<PHist_maxBins; i++)
    {
        prompt_bins[i]->Reset();
        rand_bins[i][0]->Reset();
        rand_bins[i][1]->Reset();
        result_bins[i]->Reset();
    }
}

void    PHist::Write(TDirectory& dir)
{
    TDirectory* curDir  = dir.GetDirectory("prompt");
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir("prompt");
        curDir  = dir.GetDirectory("prompt");
    }
    curDir->cd();
    prompt->Write();
    TDirectory* binsDir  = curDir->GetDirectory("bins");
    if(!binsDir)
    {
        curDir->cd();
        gDirectory->mkdir("bins");
        binsDir  = curDir->GetDirectory("bins");
    }
    binsDir->cd();
    for(int i=0; i<PHist_maxBins; i++)
        prompt_bins[i]->Write();


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
    binsDir  = curDir->GetDirectory("bins");
    if(!binsDir)
    {
        curDir->cd();
        gDirectory->mkdir("bins");
        binsDir  = curDir->GetDirectory("bins");
    }
    binsDir->cd();
    for(int i=0; i<PHist_maxBins; i++)
    {
        rand_bins[i][0]->Write();
        rand_bins[i][1]->Write();
    }


    BackgroundSubtraction();
    curDir  = dir.GetDirectory("bins");
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir("bins");
        curDir  = dir.GetDirectory("bins");
    }
    curDir->cd();
    for(int i=0; i<PHist_maxBins; i++)
        result_bins[i]->Write();

    dir.cd();
    result->Write();
}


void	PHist::BackgroundSubtraction()
{
    result->Reset();
    result->Add(prompt,1);
    Double_t    ratio = (cuts[0][1] - cuts[0][0])/((cuts[1][1] - cuts[1][0]) + (cuts[2][1] - cuts[2][0]));
    result->Add(rand[0],-ratio);
    result->Add(rand[1],-ratio);
    for(int i=0; i<PHist_maxBins; i++)
    {
        result_bins[i]->Reset();
        result_bins[i]->Add(prompt_bins[i],1);
        result_bins[i]->Add(rand_bins[i][0],-ratio);
        result_bins[i]->Add(rand_bins[i][1],-ratio);
    }
}

Double_t    PHist::cuts[3][2]  =
{{-1000000, 1000000},
 {1, 0},
 {1, 0}
};

void    PHist::SetCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t Rand0Min, const Double_t Rand0Max, const Double_t Rand1Min, const Double_t Rand1Max)
{
    cuts[0][0]  = PromptMin;
    cuts[0][1]  = PromptMax;
    cuts[1][0]  = Rand0Min;
    cuts[1][1]  = Rand0Max;
    cuts[2][0]  = Rand1Min;
    cuts[2][1]  = Rand1Max;
}







PHistD::PHistD(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max)
{
    gROOT->cd();
    prompt  = new TH1D((TString(name).Prepend("prompt_")).Data(), (TString(title).Append(" (prompt)")).Data(), nBins, min, max);
    if(!prompt)
        printf("#ERROR PHistD::PHistD: Can not create History prompt");

    rand[0]  = new TH1D((TString(name).Prepend(TString("rand0_"))).Data(), (TString(title).Append(TString("rand0"))).Data(), nBins, min, max);
    if(!rand[0])
            printf("#ERROR PHistD::PHistD: Can not create History prompt");
    rand[1]  = new TH1D((TString(name).Prepend(TString("rand1_"))).Data(), (TString(title).Append(TString("rand1"))).Data(), nBins, min, max);
    if(!rand[1])
            printf("#ERROR PHistD::PHistD: Can not create History prompt");

    result  = new TH1D(name, title, nBins, min, max);
    if(!result)
        printf("#ERROR PHistD::PHistD: Can not create History prompt");


    for(int i=0; i<PHist_maxBins; i++)
    {
        prompt_bins[i]  = new TH1D(((TString(name).Prepend("prompt_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(" (prompt) for tagger channel ").Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!prompt_bins[i])
            printf("#ERROR PHistD::PHistD: Can not create History prompt");

        rand_bins[i][0]  = new TH1D(((TString(name).Prepend("rand0_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(TString(" (rand0) for tagger channel ").Append(TString::Itoa(i,10)))).Data(), nBins, min, max);
        if(!rand_bins[i][0])
                printf("#ERROR PHistD::PHistD: Can not create History prompt");
        rand_bins[i][1]  = new TH1D(((TString(name).Prepend("rand1_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(TString(" (rand1)  for tagger channel ").Append(TString::Itoa(i,10)))).Data(), nBins, min, max);
        if(!rand_bins[i][1])
                printf("#ERROR PHistD::PHistD: Can not create History prompt");

        result_bins[i]  = new TH1D((TString(name).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), TString(title).Append(TString(" for tagger channel ").Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!result_bins[i])
            printf("#ERROR PHistD::PHistD: Can not create History prompt");
    }
}

PHistD::~PHistD()
{
}








PHistI::PHistI(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max)
{
    gROOT->cd();
    prompt  = new TH1I((TString(name).Prepend("prompt_")).Data(), (TString(title).Append(" (prompt)")).Data(), nBins, min, max);
    if(!prompt)
        printf("#ERROR PHistD::PHistD: Can not create History prompt");

    rand[0]  = new TH1I((TString(name).Prepend(TString("rand0_"))).Data(), (TString(title).Append(TString("rand0"))).Data(), nBins, min, max);
    if(!rand[0])
            printf("#ERROR PHistD::PHistD: Can not create History prompt");
    rand[1]  = new TH1I((TString(name).Prepend(TString("rand1_"))).Data(), (TString(title).Append(TString("rand1"))).Data(), nBins, min, max);
    if(!rand[1])
            printf("#ERROR PHistD::PHistD: Can not create History prompt");

    result  = new TH1I(name, title, nBins, min, max);
    if(!result)
        printf("#ERROR PHistD::PHistD: Can not create History prompt");

    for(int i=0; i<PHist_maxBins; i++)
    {
        prompt_bins[i]  = new TH1I(((TString(name).Prepend("prompt_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(" (prompt) for tagger channel ").Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!prompt_bins[i])
            printf("#ERROR PHistD::PHistD: Can not create History prompt");

        rand_bins[i][0]  = new TH1I(((TString(name).Prepend("rand0_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(TString(" (rand0) for tagger channel ")).Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!rand_bins[i][0])
                printf("#ERROR PHistD::PHistD: Can not create History prompt");
        rand_bins[i][1]  = new TH1I(((TString(name).Prepend("rand1_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(TString(" (rand1)  for tagger channel ")).Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!rand_bins[i][1])
                printf("#ERROR PHistD::PHistD: Can not create History prompt");

        result_bins[i]  = new TH1I((TString(name).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), TString(title).Append(TString(" for tagger channel ").Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!result_bins[i])
            printf("#ERROR PHistD::PHistD: Can not create History prompt");
    }
}

PHistI::~PHistI()
{
}
