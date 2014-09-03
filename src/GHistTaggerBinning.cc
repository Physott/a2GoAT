#include "GHistTaggerBinning.h"
#include "GTreeTagger.h"



#define GHTB_folderName         "TaggerBinning"
#define GHTB_chanelFolderName   "Channel_"
#define GHTB_binNameSuffix      "_Bin"
#define GHTB_binTitleSuffix     " Bin "




Int_t   GHistTaggerBinning::TaggerBinningRangeMin = 0;
Int_t   GHistTaggerBinning::TaggerBinningRangeMax = -1;

void    GHistTaggerBinning::InitTaggerBinning(const Int_t min, const Int_t max)
{
    TaggerBinningRangeMin = min;
    TaggerBinningRangeMax = max;
}



GHistTaggerBinning::GHistTaggerBinning() :
    GHistThetaBinning(),
    bin()
{
    bin.SetOwner();
}

GHistTaggerBinning::GHistTaggerBinning(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t linkHistogram) :
    GHistThetaBinning(name, title, nbinsx, xlow, xup, linkHistogram),
    bin()
{
    bin.SetOwner();
}

GHistTaggerBinning::~GHistTaggerBinning()
{

}

void    GHistTaggerBinning::CreateBin()
{
    GHistThetaBinning*    hist = new GHistThetaBinning(TString(GetName()).Append(GHTB_binNameSuffix).Append(TString::Itoa(bin.GetEntriesFast()+TaggerBinningRangeMin, 10)).Data(),
                                           TString(GetTitle()).Append(GHTB_binTitleSuffix).Append(TString::Itoa(bin.GetEntriesFast()+TaggerBinningRangeMin, 10)).Data(),
                                           GetNbinsX(),
                                           GetXmin(),
                                           GetXmax(),
                                           kFALSE);
    bin.Add(hist);
}

void    GHistTaggerBinning::ExpandBin(const Int_t newSize)
{
    while(bin.GetEntriesFast()<newSize)
        CreateBin();
}

Bool_t	GHistTaggerBinning::Add(const GHistTaggerBinning* h, Double_t c)
{
    GHistThetaBinning::Add(h, c);
    for(int i=0; i<h->bin.GetEntriesFast(); i++)
    {
        if(i>=bin.GetEntriesFast())
            CreateBin();
        ((GHistThetaBinning*)bin.At(i))->Add((GHistThetaBinning*)h->bin.At(i), c);
    }
}

void    GHistTaggerBinning::CalcResult()
{
    if(bin.GetEntriesFast()==0)
    {
        GHistThetaBinning::CalcResult();
        return;
    }

    //GHistThetaBinning::Reset();
    for(int i=0; i<bin.GetEntriesFast(); i++)
    {
        ((GHistThetaBinning*)bin.At(i))->CalcResult();
        GHistThetaBinning::Add((GHistThetaBinning*)bin.At(i));
    }
}

Int_t   GHistTaggerBinning::Fill(const Double_t value, const Int_t taggerChannel, const Double_t theta)
{
    if(taggerChannel==-1)
        return GHistThetaBinning::Fill(value, theta);

    if(TaggerBinningRangeMax==-1)
    {
        if(taggerChannel>=bin.GetEntriesFast())
            ExpandBin(taggerChannel+1);
        return ((GHistThetaBinning*)bin.At(taggerChannel))->Fill(value, theta);
    }
    else
    {
        if(taggerChannel<TaggerBinningRangeMin || taggerChannel>TaggerBinningRangeMax)
            return 0;
        if(taggerChannel>=(bin.GetEntriesFast()+TaggerBinningRangeMin))
            ExpandBin(taggerChannel+1-TaggerBinningRangeMin);
        return ((GHistThetaBinning*)bin.At(taggerChannel-TaggerBinningRangeMin))->Fill(value, theta);
    }
}

Int_t   GHistTaggerBinning::Fill(const Double_t value, const GTreeTagger& tagger, const Double_t theta, const Bool_t CreateHistogramsForTaggerBinning, const Bool_t CreateHistogramsForThetaBinning)
{
    Int_t   res = 0;
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        if(CreateHistogramsForTaggerBinning)
        {
            if(CreateHistogramsForThetaBinning)
                res += Fill(value, tagger.GetTagged_ch(i), theta);
            else
                res += Fill(value, tagger.GetTagged_ch(i), -1);
        }
        else
        {
            if(CreateHistogramsForThetaBinning)
                res += Fill(value, -1, theta);
            else
                res += Fill(value, -1, -1);
        }
    }
}

void    GHistTaggerBinning::Reset(Option_t* option)
{
    GHistThetaBinning::Reset(option);
    bin.Clear();
}

void	GHistTaggerBinning::Scale(Double_t c1, Option_t* option)
{
    GHistThetaBinning::Scale(c1, option);
    for(int i=0; i<bin.GetEntriesFast(); i++)
        ((GHistThetaBinning*)bin.At(i))->Scale(c1, option);
}

void    GHistTaggerBinning::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    GHistThetaBinning::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    for(int i=0; i<bin.GetEntriesFast(); i++)
        ((GHistThetaBinning*)bin.At(i))->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void    GHistTaggerBinning::PrepareWriteList(GHistWriteList* arr, const char *name)
{
    if(!arr)
        return;

    GHistThetaBinning::PrepareWriteList(arr, name);
    if(bin.GetEntriesFast()==0)
        return;

    GHistWriteList* TaggerBinning    = arr->GetDirectory(TString(GHTB_folderName));
    for(int i=0; i<bin.GetEntriesFast(); i++)
    {
        GHistWriteList* BinDir  = TaggerBinning->GetDirectory(TString(GHTB_chanelFolderName).Append(TString::Itoa(i+TaggerBinningRangeMin, 10)));
        if(name)
            ((GHistThetaBinning*)bin.At(i))->PrepareWriteList(BinDir, TString(name).Append(GHTB_binNameSuffix).Append(TString::Itoa(i+TaggerBinningRangeMin, 10)));
        else
            ((GHistThetaBinning*)bin.At(i))->PrepareWriteList(BinDir);
    }
}

Int_t   GHistTaggerBinning::WriteWithoutCalcResult(const char* name, Int_t option, Int_t bufsize)
{
    Int_t res   = GHistThetaBinning::WriteWithoutCalcResult(name, option, bufsize);
    if(bin.GetEntriesFast()==0)
        return res;

    TDirectory* parentDir   = gDirectory;
    TDirectory* dir         = GetCreateDirectory(GHTB_folderName);

    TString nameBuffer;
    for(int i=0; i<bin.GetEntriesFast(); i++)
    {
        dir->cd();
        GetCreateDirectory(TString(GHTB_chanelFolderName).Append(TString::Itoa(i+TaggerBinningRangeMin, 10)).Data())->cd();

        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHTB_binTitleSuffix);
            nameBuffer.Append(TString::Itoa(i+TaggerBinningRangeMin, 10));
            res += ((GHistThetaBinning*)bin.At(i))->WriteWithoutCalcResult(nameBuffer.Data(), option, bufsize);
        }
        else
            res += ((GHistThetaBinning*)bin.At(i))->WriteWithoutCalcResult(0, option, bufsize);
    }

    parentDir->cd();
    return  res;
}

