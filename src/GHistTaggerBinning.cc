#include "GHistTaggerBinning.h"
#include "GTreeTagger.h"



#define GHTB_folderName         "TaggerBinning"
#define GHTB_chanelFolderName   "Channel_"
#define GHTB_binNameSuffix      "_Bin"
#define GHTB_binTitleSuffix     " Bin "



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
    //sum->Add(h->sum, c);
}

void    GHistTaggerBinning::Reset(Option_t* option)
{
    sum->Reset(option);
    array->Reset(option);
}

Int_t   GHistTaggerBinning::Fill(const Double_t value, const Int_t taggerChannel)
{
    if(taggerChannel==-1)
        return GHistScaCor::Fill(value);

    if(TaggerBinningRangeMax==-1)
    {
        if(taggerChannel>=bin.GetEntriesFast())
            ExpandBin(taggerChannel+1);
        return ((GHistScaCor*)bin.At(taggerChannel))->Fill(value);
    }
    else
    {
        if(taggerChannel<TaggerBinningRangeMin || taggerChannel>TaggerBinningRangeMax)
            return 0;
        if(taggerChannel>=(bin.GetEntriesFast()+TaggerBinningRangeMin))
            ExpandBin(taggerChannel+1-TaggerBinningRangeMin);
        return ((GHistScaCor*)bin.At(taggerChannel-TaggerBinningRangeMin))->Fill(value);
    }
}

Int_t   GHistTaggerBinning::Fill(const Double_t value, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        if(CreateHistogramsForTaggerBinning)
            Fill(value, tagger.GetTagged_ch(i));
        else
            Fill(value);
    }
}

void	GHistTaggerBinning::Scale(Double_t c1, Option_t* option)
{
    GHistScaCor::Scale(c1, option);
    for(int i=0; i<bin.GetEntriesFast(); i++)
        ((GHistScaCor*)bin.At(i))->Scale(c1, option);
}

void    GHistTaggerBinning::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    GHistScaCor::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    for(int i=0; i<bin.GetEntriesFast(); i++)
        ((GHistScaCor*)bin.At(i))->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void    GHistTaggerBinning::PrepareWriteList(GHistWriteList* arr, const char *name)
{
    if(!arr)
        return;

    GHistScaCor::PrepareWriteList(arr, name);
    if(bin.GetEntriesFast()==0)
        return;

    GHistWriteList* TaggerBinning    = arr->GetDirectory(TString(GHTB_folderName));
    for(int i=0; i<bin.GetEntriesFast(); i++)
    {
        GHistWriteList* BinDir  = TaggerBinning->GetDirectory(TString(GHTB_chanelFolderName).Append(TString::Itoa(i+TaggerBinningRangeMin, 10)));
        if(name)
            ((GHistScaCor*)bin.At(i))->PrepareWriteList(BinDir, TString(name).Append(GHTB_binNameSuffix).Append(TString::Itoa(i+TaggerBinningRangeMin, 10)));
        else
            ((GHistScaCor*)bin.At(i))->PrepareWriteList(BinDir);
    }
}

Int_t   GHistTaggerBinning::WriteWithoutCalcResult(const char* name, Int_t option, Int_t bufsize)
{
    Int_t res   = GHistScaCor::WriteWithoutCalcResult(name, option, bufsize);
    if(bin.GetEntriesFast()==0)
        return res;

    TDirectory* parentDir   = gDirectory;
    TDirectory* dir         = GetCreateDirectory(GHTB_folderName);

    TString nameBuffer;
    for(int i=0; i<bin.GetEntriesFast(); i++)
    {
        dir->cd();
        GetCreateDirectory(TString("Channel_").Append(TString::Itoa(i+TaggerBinningRangeMin, 10)).Data())->cd();

        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHTB_chanelFolderName);
            nameBuffer.Append(TString::Itoa(i+TaggerBinningRangeMin, 10));
            res += ((GHistScaCor*)bin.At(i))->WriteWithoutCalcResult(nameBuffer.Data(), option, bufsize);
        }
        else
            res += ((GHistScaCor*)bin.At(i))->WriteWithoutCalcResult(0, option, bufsize);
    }

    parentDir->cd();
    return  res;
}

