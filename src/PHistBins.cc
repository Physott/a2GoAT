#include "PHistBins.h"



PHistBins::PHistBins()
{

}

PHistBins::~PHistBins()
{
    for(int i=0; i<PHistBins_maxBins; i++)
    {
        if(prompt_bins[i])  delete  prompt_bins[i];
        if(rand_bins[i][0]) delete  rand_bins[i][0];
        if(rand_bins[i][1]) delete  rand_bins[i][1];
        if(result_bins[i])  delete  result_bins[i];
    }
}

void    PHistBins::Clear()
{
    GH1::Clear();
    for(int i=0; i<PHistBins_maxBins; i++)
    {
        prompt_bins[i]->Reset();
        rand_bins[i][0]->Reset();
        rand_bins[i][1]->Reset();
        result_bins[i]->Reset();
    }
}

void    PHistBins::Fill(const Double_t value, const Double_t taggerTime, const Int_t taggerChannel)
{
    if(taggerTime>cuts[0][0] && taggerTime<cuts[0][1])
    {
        prompt->Fill(value);
        prompt_bins[taggerChannel]->Fill(value);
    }
    if(taggerTime>cuts[1][0] && taggerTime<cuts[1][1])
    {
        rand[0]->Fill(value);
        rand_bins[taggerChannel][0]->Fill(value);
    }
    if(taggerTime>cuts[2][0] && taggerTime<cuts[2][1])
    {
        rand[1]->Fill(value);
        rand_bins[taggerChannel][1]->Fill(value);
    }
}

void    PHistBins::Fill(const Int_t value, const Double_t taggerTime, const Int_t taggerChannel)
{
    if(taggerTime>cuts[0][0] && taggerTime<cuts[0][1])
    {
        prompt->Fill(value);
        prompt_bins[taggerChannel]->Fill(value);
    }
    if(taggerTime>cuts[1][0] && taggerTime<cuts[1][1])
    {
        rand[0]->Fill(value);
        rand_bins[taggerChannel][0]->Fill(value);
    }
    if(taggerTime>cuts[2][0] && taggerTime<cuts[2][1])
    {
        rand[1]->Fill(value);
        rand_bins[taggerChannel][1]->Fill(value);
    }
}



void    PHistBins::Write(TDirectory& dir)
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
    for(int i=0; i<PHistBins_maxBins; i++)
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
    for(int i=0; i<PHistBins_maxBins; i++)
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
    for(int i=0; i<PHistBins_maxBins; i++)
        result_bins[i]->Write();

    dir.cd();
    result->Write();
}


void	PHistBins::BackgroundSubtraction()
{
    PHistBins::BackgroundSubtraction();
    for(int i=0; i<PHistBins_maxBins; i++)
    {
        result_bins[i]->Reset();
        result_bins[i]->Add(prompt_bins[i],1);
        result_bins[i]->Add(rand_bins[i][0],-backgroundSubstractionFactor);
        result_bins[i]->Add(rand_bins[i][1],-backgroundSubstractionFactor);
    }
}






PHistBinsD::PHistBinsD(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max)
{
    gROOT->cd();
    prompt  = new TH1D((TString(name).Prepend("prompt_")).Data(), (TString(title).Append(" (prompt)")).Data(), nBins, min, max);
    if(!prompt)
        printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");

    rand[0]  = new TH1D((TString(name).Prepend(TString("rand0_"))).Data(), (TString(title).Append(TString("rand0"))).Data(), nBins, min, max);
    if(!rand[0])
            printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");
    rand[1]  = new TH1D((TString(name).Prepend(TString("rand1_"))).Data(), (TString(title).Append(TString("rand1"))).Data(), nBins, min, max);
    if(!rand[1])
            printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");

    result  = new TH1D(name, title, nBins, min, max);
    if(!result)
        printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");


    for(int i=0; i<PHistBins_maxBins; i++)
    {
        prompt_bins[i]  = new TH1D(((TString(name).Prepend("prompt_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(" (prompt) for tagger channel ").Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!prompt_bins[i])
            printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");

        rand_bins[i][0]  = new TH1D(((TString(name).Prepend("rand0_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(TString(" (rand0) for tagger channel ").Append(TString::Itoa(i,10)))).Data(), nBins, min, max);
        if(!rand_bins[i][0])
                printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");
        rand_bins[i][1]  = new TH1D(((TString(name).Prepend("rand1_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(TString(" (rand1)  for tagger channel ").Append(TString::Itoa(i,10)))).Data(), nBins, min, max);
        if(!rand_bins[i][1])
                printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");

        result_bins[i]  = new TH1D((TString(name).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), TString(title).Append(TString(" for tagger channel ").Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!result_bins[i])
            printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");
    }
}

PHistBinsD::~PHistBinsD()
{
}








PHistBinsI::PHistBinsI(const TString& name, const TString& title, const Int_t nBins, const Int_t min, const Int_t max)
{
    gROOT->cd();
    prompt  = new TH1I((TString(name).Prepend("prompt_")).Data(), (TString(title).Append(" (prompt)")).Data(), nBins, min, max);
    if(!prompt)
        printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");

    rand[0]  = new TH1I((TString(name).Prepend(TString("rand0_"))).Data(), (TString(title).Append(TString("rand0"))).Data(), nBins, min, max);
    if(!rand[0])
            printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");
    rand[1]  = new TH1I((TString(name).Prepend(TString("rand1_"))).Data(), (TString(title).Append(TString("rand1"))).Data(), nBins, min, max);
    if(!rand[1])
            printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");

    result  = new TH1I(name, title, nBins, min, max);
    if(!result)
        printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");

    for(int i=0; i<PHistBins_maxBins; i++)
    {
        prompt_bins[i]  = new TH1I(((TString(name).Prepend("prompt_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(" (prompt) for tagger channel ").Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!prompt_bins[i])
            printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");

        rand_bins[i][0]  = new TH1I(((TString(name).Prepend("rand0_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(TString(" (rand0) for tagger channel ")).Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!rand_bins[i][0])
                printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");
        rand_bins[i][1]  = new TH1I(((TString(name).Prepend("rand1_")).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), (TString(title).Append(TString(" (rand1)  for tagger channel ")).Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!rand_bins[i][1])
                printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");

        result_bins[i]  = new TH1I((TString(name).Append(TString("_Bins").Append(TString::Itoa(i,10)))).Data(), TString(title).Append(TString(" for tagger channel ").Append(TString::Itoa(i,10))).Data(), nBins, min, max);
        if(!result_bins[i])
            printf("#ERROR PHistBinsD::PHistBinsD: Can not create History prompt");
    }
}

PHistBinsI::~PHistBinsI()
{
}
