#include "GHistBGSub.h"
#include "GTreeTagger.h"

using namespace std;



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
    for(UInt_t i=1; i<cutRandMin.size(); i++)
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
    for(UInt_t i=0; i<cutRandMin.size(); i++)
    {
        if ((value >= cutRandMin[i]) && (value <= cutRandMax[i]))
            return kTRUE;
    }
    return kFALSE;
}

GHistBGSub::GHistBGSub() :
    GHistLinked(),
    writeWindows(kTRUE),
    result(new GHistScaCor()),
    prompt(new GHistScaCor()),
    randSum(new GHistScaCor()),
    rand()
{
    rand.SetOwner();
}

GHistBGSub::GHistBGSub(Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    writeWindows(kTRUE),
    result(0),
    prompt(0),
    randSum(0),
    rand()
{
    rand.SetOwner();
}

GHistBGSub::GHistBGSub(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    writeWindows(kTRUE),
    result(new GHistScaCor(name, title, nbinsx, xlow, xup, linkHistogram)),
    prompt(new GHistScaCor(TString(name).Append(GHBS_promptNameSuffix), TString(title).Append(GHBS_promptTitleSuffix), nbinsx, xlow, xup, kFALSE)),
    randSum(new GHistScaCor(TString(name).Append(GHBS_randSumNameSuffix), TString(title).Append(GHBS_randSumTitleSuffix), nbinsx, xlow, xup, kFALSE)),
    rand()
{
    rand.SetOwner();
}

GHistBGSub::~GHistBGSub()
{
}



Bool_t	GHistBGSub::Add(const GHistBGSub* h, Double_t c)
{
    result->Add(h->result, c);
    prompt->Add(h->prompt, c);
    randSum->Add(h->randSum, c);
    for(Int_t i=0; i<h->rand.GetEntriesFast(); i++)
    {
        if(i>=rand.GetEntriesFast())
            CreateRandBin();
        ((GHistScaCor*)rand.At(i))->Add((GHistScaCor*)h->rand.At(i), c);
    }    
    return true;
}

Bool_t	GHistBGSub::Add(const GHistScaCor* _result, const GHistScaCor* _prompt, const GHistScaCor* _randSum, const TObjArray& _rand, const Double_t c)
{
    result->Add(_result, c);
    prompt->Add(_prompt, c);
    randSum->Add(_randSum, c);
    for(Int_t i=0; i<_rand.GetEntriesFast(); i++)
    {
        if(i>=rand.GetEntriesFast())
            CreateRandBin();
        ((GHistScaCor*)rand.At(i))->Add((GHistScaCor*)_rand.At(i), c);
    }
    return true;
}

void    GHistBGSub::CalcResult()
{
    prompt->CalcResult();
    result->Add(prompt);

    if(rand.GetEntriesFast()==0)
        return;

    if(rand.GetEntriesFast()>1)
    {
        //randSum->Reset();
        TIter   iter(&rand);
        GHistScaCor*    hist;
        while((hist=(GHistScaCor*)iter.Next()))
            randSum->Add(hist);
        randSum->CalcResult();
        result->Add(randSum, -backgroundSubstractionFactor);
    }
    else
    {
        ((GHistScaCor*)rand.At(0))->CalcResult();
        result->Add((GHistScaCor*)rand.At(0), -backgroundSubstractionFactor);
    }
    result->CorrectZeroErrors();
}

void    GHistBGSub::CreateRandBin()
{
    GHistScaCor*    hist_rand = new GHistScaCor(TString(result->GetName()).Append(GHBS_randNameSuffix).Append(TString::Itoa(rand.GetEntriesFast(), 10)).Data(),
                                                    TString(result->GetTitle()).Append(GHBS_randTitleSuffix).Append(TString::Itoa(rand.GetEntriesFast(), 10)).Data(),
                                                    result->GetNbinsX(),
                                                    result->GetXmin(),
                                                    result->GetXmax(),
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
    if(result->IsEmpty()==kFALSE)
        return kFALSE;
    if(prompt->IsEmpty()==kFALSE)
        return kFALSE;
    if(randSum->IsEmpty()==kFALSE)
        return kFALSE;
    for(Int_t i=0; i<rand.GetEntriesFast(); i++)
    {
        if(((GHistScaCor*)rand.At(i))->IsEmpty()==kFALSE)
            return kFALSE;
    }
    return kTRUE;
}

Int_t   GHistBGSub::Fill(const Double_t value, const Double_t taggerTime)
{
    if(taggerTime>=cutPromptMin && taggerTime<=cutPromptMax)
        return prompt->Fill(value);
    for(Int_t i=0; i<GetNRandCuts(); i++)
    {
        if(i>=rand.GetEntriesFast())
            ExpandRandBins(i+1);
        if(taggerTime>=cutRandMin[i] && taggerTime<=cutRandMax[i])
            return ((GHistScaCor*)rand.At(i))->Fill(value);
    }
    return 0;
}

Int_t   GHistBGSub::Fill(const Double_t value, const GTreeTagger& tagger)
{
    for(Int_t i=0; i<tagger.GetNTagged(); i++)
        Fill(value, tagger.GetTaggedTime(i));
    return 0;
}

void    GHistBGSub::Reset(Option_t* option)
{
    result->Reset(option);
    rand.Clear("C");
    randSum->Reset(option);
    prompt->Reset(option);
}

void    GHistBGSub::Scale(Double_t c1, Option_t* option)
{
    result->Scale(c1, option);
    for(Int_t i=0; i<rand.GetEntriesFast(); i++)
        ((GHistScaCor*)rand.At(i))->Scale(c1, option);
    randSum->Scale(c1, option);
    prompt->Scale(c1, option);
}

void    GHistBGSub::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    result->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    for(Int_t i=0; i<rand.GetEntriesFast(); i++)
        ((GHistScaCor*)rand.At(i))->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    randSum->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    prompt->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void    GHistBGSub::PrepareWriteList(GHistWriteList* arr, const char *name)
{
    if(!arr)
        return;

    if(GetNRandCuts()==0)
        return prompt->PrepareWriteList(arr, name);
    else if(rand.GetEntriesFast()==0)
        return result->PrepareWriteList(arr, name);
    result->PrepareWriteList(arr, name);

    if(writeWindows==kFALSE)
        return;

    GHistWriteList* BackgroundSubstraction  = arr->GetDirectory(TString(GHBS_folderName));

    TString nameBuffer;
    if(name)
    {
        nameBuffer  = name;
        nameBuffer.Append(GHBS_promptNameSuffix);
        prompt->PrepareWriteList(BackgroundSubstraction, nameBuffer.Data());
    }
    else
        prompt->PrepareWriteList(BackgroundSubstraction);

    if(rand.GetEntriesFast()>1)
    {
        GHistWriteList* RandWindow  = BackgroundSubstraction->GetDirectory(TString(GHBS_randFolderName));

        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHBS_randSumNameSuffix);
            randSum->PrepareWriteList(BackgroundSubstraction, nameBuffer.Data());
        }
        else
            randSum->PrepareWriteList(BackgroundSubstraction);

        for(Int_t i=0; i<rand.GetEntriesFast(); i++)
        {
            if(name)
            {
                nameBuffer  = name;
                nameBuffer.Append(GHBS_randNameSuffix);
                nameBuffer.Append(TString::Itoa(i, 10));
                ((GHistScaCor*)rand.At(i))->PrepareWriteList(RandWindow, nameBuffer.Data());
            }
            else
                ((GHistScaCor*)rand.At(i))->PrepareWriteList(RandWindow);
        }
    }
    else if(rand.GetEntriesFast()>0)
    {
        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHBS_randNameSuffix);
            ((GHistScaCor*)rand.At(0))->PrepareWriteList(BackgroundSubstraction, nameBuffer.Data());
        }
        else
            ((GHistScaCor*)rand.At(0))->PrepareWriteList(BackgroundSubstraction);
    }
}

Int_t    GHistBGSub::WriteWithoutCalcResult(const char* name, Int_t option, Int_t bufsize)
{
    Int_t   res = 0;
    if(name)
    {
        if(GetNRandCuts()==0)
            return prompt->Write(name, option, bufsize);
        res += result->Write(name, option, bufsize);
    }
    else
    {
        if(GetNRandCuts()==0)
            return prompt->Write(0, option, bufsize);
        res += result->Write(0, option, bufsize);
    }

    if(writeWindows==kFALSE)
        return res;

    TDirectory* parentDir   = gDirectory;
    TDirectory* dir         = GetCreateDirectory(GHBS_folderName);
    dir->cd();

    TString nameBuffer;
    if(name)
    {
        nameBuffer  = name;
        nameBuffer.Append(GHBS_promptNameSuffix);
        res += prompt->Write(nameBuffer.Data(), option, bufsize);
    }
    else
        res += prompt->Write(0, option, bufsize);

    dir->cd();
    if(rand.GetEntriesFast()>1)
    {
        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHBS_randSumNameSuffix);
            res += randSum->Write(nameBuffer.Data(), option, bufsize);
        }
        else
            res += randSum->Write(0, option, bufsize);

        dir = GetCreateDirectory(GHBS_randFolderName);
        for(Int_t i=0; i<rand.GetEntriesFast(); i++)
        {
            dir->cd();
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
    else if(rand.GetEntriesFast()>0)
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






