#include "GH1.h"
#include "GTreeTagger.h"

TDirectory& GH1::GetDirectory(TDirectory& dir, const char* name)
{
    TDirectory* curDir  = dir.GetDirectory(name);
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir(name);
        curDir  = dir.GetDirectory(name);
    }
    curDir->cd();
    return *curDir;
}

Double_t    GH1::cutPromptMin  = -1000000;
Double_t    GH1::cutPromptMax  =  1000000;
std::vector<Double_t>     GH1::cutRandMin;
std::vector<Double_t>     GH1::cutRandMax;
Double_t    GH1::backgroundSubstractionFactor  = 0;

void    GH1::InitCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t RandMin, const Double_t RandMax)
{
    cutPromptMin    = PromptMin;
    cutPromptMax    = PromptMax;
    cutRandMin.assign(1, RandMin);
    cutRandMax.assign(1, RandMax);
    backgroundSubstractionFactor = (PromptMax - PromptMin)/(RandMax - RandMin);
}

void    GH1::AddRandCut(const Double_t RandMin, const Double_t RandMax)
{
    cutRandMin.push_back(RandMin);
    cutRandMax.push_back(RandMax);
    backgroundSubstractionFactor = cutRandMax[0] - cutRandMin[0];
    for(int i=1; i<cutRandMin.size(); i++)
        backgroundSubstractionFactor += cutRandMax[i] - cutRandMin[i];
    backgroundSubstractionFactor    = (cutPromptMax - cutPromptMin)/backgroundSubstractionFactor;
}






GH1::GH1(const TString& _Name, const TString& _Title, const Int_t NumberOfBins) :
    TNamed(_Name, _Title),
    nBins(NumberOfBins),
    prompt(1),
    rand(2)
{
    prompt.SetOwner();
    rand.SetOwner();
}

GH1::GH1(const char* _Name, const char* _Title, const Int_t NumberOfBins) :
    TNamed(_Name, _Title),
    nBins(NumberOfBins),
    prompt(1),
    rand(2)
{
    prompt.SetOwner();
    rand.SetOwner();
}

GH1::~GH1()
{
}

void    GH1::Clear()
{
    prompt.Clear();
    rand.Clear();
}

void  GH1::AddBin(const Int_t channel)
{
    TH1*    hist = AddHistogram((TString(GetName()).Prepend("prompt_")).Append(TString::Itoa(channel, 10).Prepend("_Bin")).Data(), (TString(GetTitle()).Append(TString::Itoa(channel, 10).Prepend(" channel ")).Append(" (prompt)")).Data());
    prompt.AddAtAndExpand(hist, channel);

    while(rand.GetEntriesFast()<GetNRandCuts())
    {
        gROOT->cd();
        TObjArray* array    = new TObjArray(channel+1);
        rand.AddAtAndExpand(array, rand.GetEntriesFast());
    }

    for(int i=0; i<cutRandMin.size(); i++)
    {
        TH1*    hist = AddHistogram((TString(GetName()).Prepend((TString::Itoa(i,10).Prepend("rand")).Append("_"))).Append(TString::Itoa(channel, 10).Prepend("_Bin")).Data(), TString(GetTitle()).Append(TString::Itoa(channel, 10).Prepend(" channel ")).Append(TString::Itoa(i,10).Prepend(" (rand")).Append(")").Data());
        ((TObjArray*)rand.At(i))->AddAtAndExpand(hist, channel);
    }
}

TH1*	GH1::BackgroundSubtractionBin(const Int_t channel, const TH1* sumRand)
{
    gROOT->cd();
    TH1*    hist    = AddHistogram(TString(GetName()).Append(TString::Itoa(channel, 10).Prepend("_Bin")).Data(), TString(GetTitle()).Append(TString::Itoa(channel, 10).Prepend(" channel ")).Data());
    hist->Add((TH1*)(prompt.At(channel)));
    hist->Add(sumRand, -backgroundSubstractionFactor);
    return hist;
}

void    GH1::Fill(const Double_t value, const Double_t taggerTime, const Int_t taggerChannel)
{
    if(taggerTime>=cutPromptMin && taggerTime<=cutPromptMax)
    {
        if(taggerChannel>=prompt.GetEntriesFast())
        {
            Int_t   start = prompt.GetEntriesFast();
            for(int i=taggerChannel; i>=start; i--)
                AddBin(i);
        }
        ((TH1*)prompt.At(taggerChannel))->Fill(value);
    }

    for(int i=0; i<cutRandMin.size(); i++)
    {
        if(taggerTime>=cutRandMin[i] && taggerTime<=cutRandMax[i])
        {
            if(i>=rand.GetEntriesFast())
            {
                TObjArray*  array = new TObjArray(1);
                rand.AddAtAndExpand(array, i);
            }
            if(taggerChannel>=((TObjArray*)rand.At(i))->GetEntriesFast())
            {
                Int_t   start = prompt.GetEntriesFast();
                for(int i=taggerChannel; i>=start; i--)
                    AddBin(i);
            }
            ((TH1*)(((TObjArray*)rand.At(i))->At(taggerChannel)))->Fill(value);
        }
    }
}

void    GH1::Fill(const Int_t value, const Double_t taggerTime, const Int_t taggerChannel)
{
    if(taggerTime>=cutPromptMin && taggerTime<=cutPromptMax)
    {
        if(taggerChannel>=prompt.GetEntriesFast())
        {
            Int_t   start = prompt.GetEntriesFast();
            for(int i=taggerChannel; i>=start; i--)
                AddBin(i);
        }
        ((TH1*)prompt.At(taggerChannel))->Fill(value);
    }

    for(int i=0; i<cutRandMin.size(); i++)
    {
        if(taggerTime>=cutRandMin[i] && taggerTime<=cutRandMax[i])
        {
            if(i>=rand.GetEntriesFast())
            {
                TObjArray*  array = new TObjArray(1);
                rand.AddAtAndExpand(array, i);
            }
            if(taggerChannel>=((TObjArray*)rand.At(i))->GetEntriesFast())
            {
                Int_t   start = prompt.GetEntriesFast();
                for(int i=taggerChannel; i>=start; i--)
                    AddBin(i);
            }
            ((TH1*)(((TObjArray*)rand.At(i))->At(taggerChannel)))->Fill(value);
        }
    }
}

void    GH1::Fill(const Double_t value, const GTreeTagger& tagger)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
        Fill(value, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
}

void    GH1::Fill(const Int_t value, const GTreeTagger& tagger)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
        Fill(value, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
}

TH1*    GH1::SumRand(const Int_t RandCut)
{
    gROOT->cd();
    TH1*    hist = AddHistogram((TString(GetName()).Prepend((TString::Itoa(RandCut,10).Prepend("rand")).Append("_"))).Data(), TString(GetTitle()).Append(TString::Itoa(RandCut,10).Prepend(" (rand")).Append(")").Data());

    if(!rand.At(RandCut))
        return 0;
    TIter   iter((TObjArray*)rand.At(RandCut));
    TH1*    help;
    while(help = (TH1*)iter.Next())
        hist->Add(help);
    return hist;
}

TH1*    GH1::SumRandAll()
{
    gROOT->cd();
    TH1*    hist    = AddHistogram((TString(GetName()).Prepend("randAll_")).Data(), TString(GetTitle()).Append(" (rand all)").Data());
    for(int i=0; i<rand.GetEntriesFast(); i++)
    {
        TH1*    help = SumRand(i);
        if(help)
            hist->Add(help);
    }
    return hist;
}

TH1*    GH1::SumRandBin(const Int_t channel)
{
    gROOT->cd();
    TH1*    hist    = AddHistogram((TString(GetName()).Prepend("randAll_")).Append(TString::Itoa(channel, 10).Prepend("_Bin")).Data(), TString(GetTitle()).Append(TString::Itoa(channel, 10).Prepend(" channel ")).Append(" (rand all)").Data());

    TIter   iter(&rand);
    TObjArray*    array;
    while(array = (TObjArray*)iter.Next())
        hist->Add((TH1*)array->At(channel));
    return hist;
}

TH1*    GH1::SumPrompt()
{
    gROOT->cd();
    TH1*    hist    = AddHistogram((TString(GetName()).Prepend("prompt_")).Data(), TString(GetTitle()).Append(" (prompt)").Data());

    TIter   iter(&prompt);
    TH1*    help;
    while(help = (TH1*)iter.Next())
        hist->Add(help);
    return hist;
}

void    GH1::Write(TDirectory& dir)
{
    if(prompt.GetEntriesFast()>1)
    {
        TH1* sum;
        TH1* res;
        for(int i=0; i<prompt.GetEntriesFast(); i++)
        {
            WriteBin(GetDirectory(dir, "bins"), i);
            sum = SumRandBin(i);
            GetDirectory(GetDirectory(dir, "bins"), "rand");
            sum->Write();
            res = BackgroundSubtractionBin(i, sum);
            GetDirectory(dir, "bins");
            res->Write();
            sum->Delete();
            res->Delete();
        }
        for(int i=0; i<rand.GetEntriesFast(); i++)
        {
            sum = SumRand(i);
            GetDirectory(dir, "rand");
            sum->Write();
            sum->Delete();
        }
        sum = SumRandAll();
        GetDirectory(dir, "rand");
        sum->Write();
        res = SumPrompt();
        GetDirectory(dir, "prompt");
        res->Write();
        gROOT->cd();
        TH1*    endResult    = AddHistogram(GetName(), GetTitle());
        endResult->Add(res);
        endResult->Add(sum, -backgroundSubstractionFactor);
        dir.cd();
        endResult->Write();

        sum->Delete();
        res->Delete();
        return;
    }
    else if(prompt.GetEntriesFast()==1)
    {
        if(!prompt.At(0))
            return;
        if(rand.GetEntriesFast()==0)
        {
            ((TH1*)prompt.At(0))->SetNameTitle(GetName(), GetTitle());
            dir.cd();
            ((TH1*)prompt.At(0))->Write();
            return;
        }
        ((TH1*)prompt.At(0))->SetNameTitle((TString(GetName()).Prepend("prompt_")).Data(), (TString(GetTitle()).Append(" (prompt)")).Data());
        GetDirectory(dir, "prompt");
        ((TH1*)prompt.At(0))->Write();

        TIter   iter(&rand);
        TObjArray*    array;
        Int_t i=-1;
        while(array = (TObjArray*)iter.Next())
        {
            i++;
            if(!array->At(0))
                continue;
            ((TH1*)(array->At(0)))->SetNameTitle((TString(GetName()).Prepend((TString::Itoa(i,10).Prepend("rand")).Append("_"))).Data(), TString(GetTitle()).Append(TString::Itoa(i,10).Prepend(" (rand")).Append(")").Data());
            GetDirectory(dir, "rand");
            ((TH1*)(array->At(0)))->Write();
        }

        if(rand.GetEntriesFast()>1)
        {
            TH1* sum = SumRandBin(0);
            sum->SetNameTitle((TString(GetName()).Prepend("randAll_")).Data(), TString(GetTitle()).Append(" (rand all)").Data());
            GetDirectory(dir, "rand");
            sum->Write();
            TH1* res = BackgroundSubtractionBin(0, sum);
            res->SetNameTitle(GetName(), GetTitle());
            dir.cd();
            res->Write();
            sum->Delete();
            res->Delete();
        }
        else if(rand.GetEntriesFast()==1)
        {
            TH1* res = BackgroundSubtractionBin(0, ((TH1*)(((TObjArray*)rand.At(0))->At(0))));
            res->SetNameTitle(GetName(), GetTitle());
            dir.cd();
            res->Write();
            res->Delete();
        }
    }
}

void    GH1::WriteBin(TDirectory& dir, const Int_t channel)
{
    GetDirectory(dir, "prompt");
    if(prompt.At(channel))
        ((TH1*)prompt.At(channel))->Write();

    TIter   iter(&rand);
    TObjArray*    array;
    Int_t i=-1;
    while(array = (TObjArray*)iter.Next())
    {
        i++;
        if(!array->At(channel))
            continue;
        GetDirectory(dir, "rand");
        ((TH1*)(array->At(channel)))->Write();
    }
}











GH1D::GH1D(const TString& name, const TString& title, const Int_t NumberOfBins, const Double_t min, const Double_t max) :
    GH1(name, title, NumberOfBins),
    minBin(min),
    maxBin(max)
{
}
GH1D::GH1D(const char* name, const char* title, const Int_t NumberOfBins, const Double_t min, const Double_t max) :
    GH1(name, title, NumberOfBins),
    minBin(min),
    maxBin(max)
{
}

GH1D::~GH1D()
{
}

TH1*    GH1D::AddHistogram(const TString& _Name, const TString& _Title)
{
    TH1*    hist = new TH1D(_Name.Data(), _Title.Data(), nBins, minBin, maxBin);
    return hist;
}









GH1I::GH1I(const TString& name, const TString& title, const Int_t NumberOfBins, const Int_t min, const Int_t max) :
    GH1(name, title, NumberOfBins),
    minBin(min),
    maxBin(max)
{
}
GH1I::GH1I(const char* name, const char* title, const Int_t NumberOfBins, const Int_t min, const Int_t max) :
    GH1(name, title, NumberOfBins),
    minBin(min),
    maxBin(max)
{
}

GH1I::~GH1I()
{
}

TH1*    GH1I::AddHistogram(const TString& _Name, const TString& _Title)
{
    gROOT->cd();
    TH1*    hist = new TH1I(_Name.Data(), _Title.Data(), nBins, minBin, maxBin);
    return hist;
}

