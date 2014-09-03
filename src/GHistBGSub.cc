#include "GHistBGSub.h"
#include "GTreeTagger.h"

using namespace std;


#define GHBS_folderName         "BackgroundSubstraction"
#define GHBS_randFolderName     "RandomWindow"
#define GHBS_subRandFolderName  "RandWin_"
#define GHBS_promptFolderName   "PromptWindow"
#define GHBS_randNameSuffix     "_Rand"
#define GHBS_randTitleSuffix    " Rand "
#define GHBS_randSumNameSuffix  "_RandSum"
#define GHBS_randSumTitleSuffix " RandSum"
#define GHBS_promptNameSuffix   "_Prompt"
#define GHBS_promptTitleSuffix  " Prompt"



Double_t    GHistBGSub::cutPromptMin  = -1000000;
Double_t    GHistBGSub::cutPromptMax  =  1000000;
std::vector<Double_t>     GHistBGSub::cutRandMin;
std::vector<Double_t>     GHistBGSub::cutRandMax;
Double_t    GHistBGSub::backgroundSubstractionFactor  = 0;

void    GHistBGSub::InitCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t RandMin, const Double_t RandMax)
{
    cutPromptMin    = PromptMin;
    cutPromptMax    = PromptMax;
    cutRandMin.assign(1, RandMin);
    cutRandMax.assign(1, RandMax);
    backgroundSubstractionFactor = (PromptMax - PromptMin)/(RandMax - RandMin);
}

void    GHistBGSub::AddRandCut(const Double_t RandMin, const Double_t RandMax)
{
    cutRandMin.push_back(RandMin);
    cutRandMax.push_back(RandMax);
    backgroundSubstractionFactor = cutRandMax[0] - cutRandMin[0];
    for(int i=1; i<cutRandMin.size(); i++)
        backgroundSubstractionFactor += cutRandMax[i] - cutRandMin[i];
    backgroundSubstractionFactor    = (cutPromptMax - cutPromptMin)/backgroundSubstractionFactor;
}





GHistBGSub::GHistBGSub() :
    GHistTaggerBinning(),
    rand(),
    randSum(),
    prompt()
{
    rand.SetOwner();
}

GHistBGSub::GHistBGSub(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t linkHistogram) :
    GHistTaggerBinning(name, title, nbinsx, xlow, xup, linkHistogram),
    rand(),
    randSum(TString(name).Append(GHBS_randSumNameSuffix), TString(title).Append(GHBS_randSumTitleSuffix), nbinsx, xlow, xup, kFALSE),
    prompt(TString(name).Append(GHBS_promptNameSuffix), TString(title).Append(GHBS_promptTitleSuffix), nbinsx, xlow, xup, kFALSE)
{
    rand.SetOwner();
}

GHistBGSub::~GHistBGSub()
{
}

Bool_t	GHistBGSub::Add(const GHistBGSub* h, Double_t c)
{
    GHistTaggerBinning::Add((GHistTaggerBinning*)h, c);
    for(int i=0; i<h->rand.GetEntriesFast(); i++)
    {
        if(i>=rand.GetEntriesFast())
            CreateRandBin();
        ((GHistTaggerBinning*)rand.At(i))->Add((GHistTaggerBinning*)h->rand.At(i), c);
    }
}

void    GHistBGSub::CalcResult()
{
    prompt.CalcResult();
    //GHistTaggerBinning::Reset();
    GHistTaggerBinning::Add(&prompt);

    if(rand.GetEntriesFast()==0)
        return;

    TIter   iter(&rand);
    GHistTaggerBinning*    hist;
    while(hist=(GHistTaggerBinning*)iter.Next())
        hist->CalcResult();

    if(rand.GetEntriesFast()>1)
    {
        //randSum.Reset();
        TIter   iter(&rand);
        GHistTaggerBinning*    hist;
        while(hist=(GHistTaggerBinning*)iter.Next())
            randSum.Add(hist);
        GHistTaggerBinning::Add(&randSum, -backgroundSubstractionFactor);
    }
    else
        GHistTaggerBinning::Add((GHistTaggerBinning*)rand.At(0), -backgroundSubstractionFactor);
}

void    GHistBGSub::CreateRandBin()
{
    GHistTaggerBinning*    hist_rand = new GHistTaggerBinning(TString(GetName()).Append(GHBS_randNameSuffix).Append(TString::Itoa(rand.GetEntriesFast(), 10)).Data(),
                                                    TString(GetTitle()).Append(GHBS_randTitleSuffix).Append(TString::Itoa(rand.GetEntriesFast(), 10)).Data(),
                                                    GetNbinsX(),
                                                    GetXmin(),
                                                    GetXmax(),
                                                    kFALSE);
    rand.AddAtFree(hist_rand);
}

void    GHistBGSub::ExpandRandBins(const Int_t newSize)
{
    while(rand.GetEntriesFast()<newSize)
        CreateRandBin();
}

Int_t   GHistBGSub::Fill(const Double_t value, const Double_t taggerTime, const Int_t taggerChannel, const Double_t theta)
{
    if(taggerTime>=cutPromptMin && taggerTime<=cutPromptMax)
        return prompt.Fill(value, taggerChannel, theta);
    for(int i=0; i<GetNRandCuts(); i++)
    {
        if(i>=rand.GetEntriesFast())
            ExpandRandBins(i+1);
        if(taggerTime>=cutRandMin[i] && taggerTime<=cutRandMax[i])
            return ((GHistTaggerBinning*)rand.At(i))->Fill(value, taggerChannel, theta);
    }
}

void    GHistBGSub::Reset(Option_t* option)
{
    GHistTaggerBinning::Reset(option);
    rand.Clear();
    randSum.Reset(option);
    prompt.Reset(option);
}


void    GHistBGSub::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    GHistTaggerBinning::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    for(int i=0; i<rand.GetEntriesFast(); i++)
        ((GHistTaggerBinning*)rand.At(i))->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    randSum.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    prompt.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void    GHistBGSub::PrepareWriteList(GHistWriteList* arr, const char *name)
{
    if(!arr)
        return;

    Int_t   res = 0;
    if(name)
    {
        if(GetNRandCuts()==0 || rand.GetEntriesFast()==0)
            return prompt.PrepareWriteList(arr, name);
        GHistTaggerBinning::PrepareWriteList(arr, name);
    }
    else
    {
        if(GetNRandCuts()==0 || rand.GetEntriesFast()==0)
            return prompt.PrepareWriteList(arr, GHistTaggerBinning::GetName());
        GHistTaggerBinning::PrepareWriteList(arr);
    }

    GHistWriteList* BackgroundSubstraction  = arr->GetDirectory(TString(GHBS_folderName));
    GHistWriteList* PromptWindow            = BackgroundSubstraction->GetDirectory(TString(GHBS_promptFolderName));

    TString nameBuffer;
    if(name)
    {
        nameBuffer  = name;
        nameBuffer.Append(GHBS_promptNameSuffix);
        prompt.PrepareWriteList(PromptWindow, nameBuffer.Data());
    }
    else
        prompt.PrepareWriteList(PromptWindow);

    GHistWriteList* RandWindow              = BackgroundSubstraction->GetDirectory(TString(GHBS_randFolderName));
    if(rand.GetEntriesFast()>1)
    {
        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHBS_randSumNameSuffix);
            randSum.PrepareWriteList(RandWindow, nameBuffer.Data());
        }
        else
            randSum.PrepareWriteList(RandWindow);

        for(int i=0; i<rand.GetEntriesFast(); i++)
        {
            GHistWriteList* subRandWindow = RandWindow->GetDirectory(TString(GHBS_subRandFolderName).Append(TString::Itoa(i, 10)));
            if(name)
            {
                nameBuffer  = name;
                nameBuffer.Append(GHBS_randNameSuffix);
                nameBuffer.Append(TString::Itoa(i, 10));
                ((GHistTaggerBinning*)rand.At(i))->PrepareWriteList(subRandWindow, nameBuffer.Data());
            }
            else
                ((GHistTaggerBinning*)rand.At(i))->PrepareWriteList(subRandWindow);
        }
    }
    else
    {
        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHBS_randNameSuffix);
            ((GHistTaggerBinning*)rand.At(0))->PrepareWriteList(RandWindow, nameBuffer.Data());
        }
        else
            ((GHistTaggerBinning*)rand.At(0))->PrepareWriteList(RandWindow);
    }
}

Int_t    GHistBGSub::WriteWithoutCalcResult(const char* name, Int_t option, Int_t bufsize)
{
    Int_t   res = 0;
    if(name)
    {
        if(GetNRandCuts()==0)
            return prompt.Write(name, option, bufsize);
        res += GHistTaggerBinning::Write(name, option, bufsize);
    }
    else
    {
        if(GetNRandCuts()==0)
            return prompt.Write(0, option, bufsize);
        res += GHistTaggerBinning::Write(0, option, bufsize);
    }

    TDirectory* parentDir   = gDirectory;
    TDirectory* dir         = GetCreateDirectory(GHBS_folderName);
    dir->cd();
    GetCreateDirectory(GHBS_promptFolderName)->cd();
    TString nameBuffer;
    if(name)
    {
        nameBuffer  = name;
        nameBuffer.Append(GHBS_promptNameSuffix);
        res += prompt.Write(nameBuffer.Data(), option, bufsize);
    }
    else
        res += prompt.Write(0, option, bufsize);

    dir->cd();
    dir = GetCreateDirectory(GHBS_randFolderName);
    dir->cd();
    if(rand.GetEntriesFast()>1)
    {
        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHBS_randSumNameSuffix);
            res += randSum.Write(nameBuffer.Data(), option, bufsize);
        }
        else
            res += randSum.Write(0, option, bufsize);

        for(int i=0; i<rand.GetEntriesFast(); i++)
        {
            dir->cd();
            GetCreateDirectory(TString(GHBS_subRandFolderName).Append(TString::Itoa(i, 10)).Data())->cd();

            if(name)
            {
                nameBuffer  = name;
                nameBuffer.Append(GHBS_randNameSuffix);
                nameBuffer.Append(TString::Itoa(i, 10));
                res += ((GHistScaCor*)rand.At(i))->Write(nameBuffer.Data(), option, bufsize);
            }
            else
                res += ((GHistScaCor*)rand.At(i))->Write(0, option, bufsize);
        }
    }
    else
    {
        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHBS_randNameSuffix);
            res += ((GHistScaCor*)rand.At(0))->Write(nameBuffer.Data(), option, bufsize);
        }
        else
            res += ((GHistScaCor*)rand.At(0))->Write(0, option, bufsize);
    }
    parentDir->cd();
    return res;
}






