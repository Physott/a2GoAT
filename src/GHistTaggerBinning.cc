#include "GHistTaggerBinning.h"
#include "GTreeTagger.h"



#define GHTB_folderName         "TaggerBinning"
#define GHTB_binNameSuffix      "_Bins"
#define GHTB_binTitleSuffix     " Bins "



GHistTaggerBinning::GHistTaggerBinning(Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    sum(0),
    array(0)
{
}

GHistTaggerBinning::GHistTaggerBinning() :
    GHistLinked(),
    sum(new GHistBGSub()),
    array(new GHistBGSub2())
{
}

GHistTaggerBinning::GHistTaggerBinning(const char* name, const char* title,
                                       const Int_t nbinsx, const Double_t xlow, const Double_t xup,
                                       const Int_t nTaggerBins,
                                       const Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    sum(new GHistBGSub(name, title, nbinsx, xlow, xup, kFALSE)),
    array(new GHistBGSub2(name, title, nbinsx, xlow, xup, nTaggerBins, 0, nTaggerBins, kFALSE))
{
}

GHistTaggerBinning::~GHistTaggerBinning()
{

}

Bool_t	GHistTaggerBinning::Add(const GHistTaggerBinning* h, Double_t c)
{
    sum->Add(h->sum, c);
    array->Add(h->array, c);
}

void    GHistTaggerBinning::CalcResult()
{
    GHistBGSub*    help    = ((GHistBGSub2*)array)->ProjectionX();
    sum->Add(help);
    if(help)
        delete help;
}

Int_t   GHistTaggerBinning::Fill(const Double_t value, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        if(CreateHistogramsForTaggerBinning)
            Fill(value, tagger.GetTagged_ch(i), tagger.GetTagged_t(i));
        else
            Fill(value);
    }
}

void    GHistTaggerBinning::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    sum->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    array->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void    GHistTaggerBinning::PrepareWriteList(GHistWriteList* arr, const char *name)
{
    if(!arr)
        return;

    sum->PrepareWriteList(arr, name);
    if(array->IsEmpty()==kTRUE)
        return;

    GHistWriteList* TaggerBinning    = arr->GetDirectory(TString(GHTB_folderName));
    if(name)
        array->PrepareWriteList(TaggerBinning, TString(name).Append(GHTB_binNameSuffix));
    else
        array->PrepareWriteList(TaggerBinning);
}

Int_t   GHistTaggerBinning::WriteWithoutCalcResult(const char* name, Int_t option, Int_t bufsize)
{
    Int_t res   = sum->WriteWithoutCalcResult(name, option, bufsize);
    if(array->IsEmpty()==kTRUE)
        return res;

    TDirectory* parentDir   = gDirectory;
    TDirectory* dir         = GetCreateDirectory(GHTB_folderName);
    dir->cd();

    TString nameBuffer;
    if(name)
    {
        nameBuffer  = name;
        nameBuffer.Append(GHTB_binNameSuffix);
        res += array->WriteWithoutCalcResult(nameBuffer.Data(), option, bufsize);
    }
    else
        res += array->WriteWithoutCalcResult(0, option, bufsize);

    parentDir->cd();
    return  res;
}

