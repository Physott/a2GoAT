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

Bool_t    GHistBGSub::IsPrompt(const Double_t value)
{
   if ((value >= cutPromptMin) && (value <= cutPromptMax))
       return kTRUE;
   return kFALSE;
}

Bool_t    GHistBGSub::IsRandom(const Double_t value)
{
    for(int i=0; i<cutRandMin.size(); i++)
    {
        if ((value >= cutRandMin[i]) && (value <= cutRandMax[i]))
            return kTRUE;
    }
    return kFALSE;
}

GHistBGSub::GHistBGSub() :
    GHistScaCor(),
    rand(),
    randSum(),
    prompt(),
    writeWindows(kTRUE)
{
    rand.SetOwner();
}

GHistBGSub::GHistBGSub(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t linkHistogram) :
    GHistScaCor(name, title, nbinsx, xlow, xup, linkHistogram),
    rand(),
    randSum(TString(name).Append(GHBS_randSumNameSuffix), TString(title).Append(GHBS_randSumTitleSuffix), nbinsx, xlow, xup, kFALSE),
    prompt(TString(name).Append(GHBS_promptNameSuffix), TString(title).Append(GHBS_promptTitleSuffix), nbinsx, xlow, xup, kFALSE),
    writeWindows(kTRUE)
{
    rand.SetOwner();
}

GHistBGSub::~GHistBGSub()
{
}

Bool_t	GHistBGSub::Add(const GHistBGSub* h, Double_t c)
{
    GHistScaCor::Add((GHistScaCor*)h, c);
    prompt.Add(&h->prompt, c);
    randSum.Add(&h->randSum, c);
    for(int i=0; i<h->rand.GetEntriesFast(); i++)
    {
        if(i>=rand.GetEntriesFast())
            CreateRandBin();
        ((GHistScaCor*)rand.At(i))->Add((GHistScaCor*)h->rand.At(i), c);
    }
}

void    GHistBGSub::CalcResult()
{
    GHistScaCor::Add(&prompt);

    if(rand.GetEntriesFast()==0)
        return;

    if(rand.GetEntriesFast()>1)
    {
        //randSum.Reset();
        TIter   iter(&rand);
        GHistScaCor*    hist;
        while(hist=(GHistScaCor*)iter.Next())
            randSum.Add(hist);
        GHistScaCor::Add(&randSum, -backgroundSubstractionFactor);
    }
    else
        GHistScaCor::Add((GHistScaCor*)rand.At(0), -backgroundSubstractionFactor);
}

void    GHistBGSub::CreateRandBin()
{
    GHistScaCor*    hist_rand = new GHistScaCor(TString(GetName()).Append(GHBS_randNameSuffix).Append(TString::Itoa(rand.GetEntriesFast(), 10)).Data(),
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

Bool_t  GHistBGSub::IsEmpty()
{
    if(GHistScaCor::IsEmpty()==kFALSE)
        return kFALSE;
    if(prompt.IsEmpty()==kFALSE)
        return kFALSE;
    if(randSum.IsEmpty()==kFALSE)
        return kFALSE;
    for(int i=0; i<rand.GetEntriesFast(); i++)
    {
        if(((GHistScaCor*)rand.At(i))->IsEmpty()==kFALSE)
            return kFALSE;
    }
    return kTRUE;
}

Int_t   GHistBGSub::Fill(const Double_t value, const Double_t taggerTime)
{
    if(taggerTime>=cutPromptMin && taggerTime<=cutPromptMax)
        return prompt.Fill(value);
    for(int i=0; i<GetNRandCuts(); i++)
    {
        if(i>=rand.GetEntriesFast())
            ExpandRandBins(i+1);
        if(taggerTime>=cutRandMin[i] && taggerTime<=cutRandMax[i])
            return ((GHistScaCor*)rand.At(i))->Fill(value);
    }
}

Int_t   GHistBGSub::Fill(const Double_t value, const GTreeTagger& tagger)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
        Fill(value, tagger.GetTagged_t(i));
}

void    GHistBGSub::Reset(Option_t* option)
{
    GHistScaCor::Reset(option);
    rand.Clear("C");
    randSum.Reset(option);
    prompt.Reset(option);
}


void    GHistBGSub::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    GHistScaCor::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    for(int i=0; i<rand.GetEntriesFast(); i++)
        ((GHistScaCor*)rand.At(i))->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    randSum.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    prompt.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void    GHistBGSub::PrepareWriteList(GHistWriteList* arr, const char *name)
{
    if(!arr)
        return;

    if(name)
    {
        if(GetNRandCuts()==0 || rand.GetEntriesFast()==0)
            return GHistScaCor::PrepareWriteList(arr, name);
        GHistScaCor::PrepareWriteList(arr, name);
    }
    else
    {
        if(GetNRandCuts()==0 || rand.GetEntriesFast()==0)
            return GHistScaCor::PrepareWriteList(arr);
        GHistScaCor::PrepareWriteList(arr);
    }

    if(writeWindows==kFALSE)
        return;

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
                ((GHistScaCor*)rand.At(i))->PrepareWriteList(subRandWindow, nameBuffer.Data());
            }
            else
                ((GHistScaCor*)rand.At(i))->PrepareWriteList(subRandWindow);
        }
    }
    else
    {
        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHBS_randNameSuffix);
            ((GHistScaCor*)rand.At(0))->PrepareWriteList(RandWindow, nameBuffer.Data());
        }
        else
            ((GHistScaCor*)rand.At(0))->PrepareWriteList(RandWindow);
    }
}

Int_t    GHistBGSub::WriteWithoutCalcResult(const char* name, Int_t option, Int_t bufsize)
{
    Int_t   res = 0;
    if(name)
    {
        if(GetNRandCuts()==0)
            return prompt.Write(name, option, bufsize);
        res += GHistScaCor::Write(name, option, bufsize);
    }
    else
    {
        if(GetNRandCuts()==0)
            return prompt.Write(0, option, bufsize);
        res += GHistScaCor::Write(0, option, bufsize);
    }

    if(writeWindows==kFALSE)
        return res;

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






