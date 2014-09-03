#include "GHistThetaBinning.h"
#include "GTreeTagger.h"



#define GHThB_folderName         "ThetaBinning"
#define GHThB_chanelFolderName   "Theta_"
#define GHThB_binNameSuffix      "_Th"
#define GHThB_binTitleSuffix     " Theta "



Int_t       GHistThetaBinning::ThetaBinningCount    = 18;
Double_t    GHistThetaBinning::ThetaBinningWidth    = 10;
Double_t    GHistThetaBinning::ThetaBinningRangeMin = 0;
Double_t    GHistThetaBinning::ThetaBinningRangeMax = 180;

void    GHistThetaBinning::InitTaggerBinning(const Int_t count, const Double_t min, const Double_t max)
{
    ThetaBinningCount    = count;
    ThetaBinningRangeMin = min;
    ThetaBinningRangeMax = max;
    ThetaBinningWidth    = (ThetaBinningRangeMax - ThetaBinningRangeMin)/ThetaBinningCount;
}



GHistThetaBinning::GHistThetaBinning() :
    GHistScaCor(),
    bin()
{
    bin.SetOwner();
}

GHistThetaBinning::GHistThetaBinning(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Bool_t linkHistogram) :
    GHistScaCor(name, title, nbinsx, xlow, xup, linkHistogram),
    bin()
{
    bin.SetOwner();
}

GHistThetaBinning::~GHistThetaBinning()
{

}

Bool_t	GHistThetaBinning::Add(const GHistThetaBinning* h, Double_t c)
{
    GHistScaCor::Add((GHistScaCor*)h, c);
    for(int i=0; i<h->bin.GetEntriesFast(); i++)
    {
        if(i>=bin.GetEntriesFast())
            CreateBin();
        ((GHistScaCor*)bin.At(i))->Add((GHistScaCor*)h->bin.At(i), c);
    }
}

void    GHistThetaBinning::CalcResult()
{
    if(bin.GetEntriesFast()==0)
    {
        GHistScaCor::CalcResult();
        return;
    }

    //GHistScaCor::Reset();
    for(int i=0; i<bin.GetEntriesFast(); i++)
    {
        ((GHistScaCor*)bin.At(i))->CalcResult();
        GHistScaCor::Add((GHistScaCor*)bin.At(i));
    }
}

void    GHistThetaBinning::CreateBin()
{
    GHistScaCor*    hist = new GHistScaCor(TString(GetName()).Append(GHThB_binNameSuffix).Append(TString::Itoa(ThetaBinningRangeMin+(bin.GetEntriesFast()*ThetaBinningWidth), 10)).Append("_").Append(TString::Itoa(ThetaBinningRangeMin+((bin.GetEntriesFast()+1)*ThetaBinningWidth), 10)).Data(),
                                           TString(GetTitle()).Append(GHThB_binTitleSuffix).Append(TString::Itoa(ThetaBinningRangeMin+(bin.GetEntriesFast()*ThetaBinningWidth), 10)).Append("_").Append(TString::Itoa(ThetaBinningRangeMin+((bin.GetEntriesFast()+1)*ThetaBinningWidth), 10)).Data(),
                                           GetNbinsX(),
                                           GetXmin(),
                                           GetXmax(),
                                           kFALSE);
    bin.Add(hist);
}

void    GHistThetaBinning::ExpandBin(const Int_t newSize)
{
    while(bin.GetEntriesFast()<newSize)
        CreateBin();
}

Int_t   GHistThetaBinning::GetThetaBin(const Double_t theta)
{
    if(theta<ThetaBinningRangeMin || theta>=ThetaBinningRangeMax)
        return -1;
    return Int_t((theta-ThetaBinningRangeMin)/ThetaBinningWidth);
}

Int_t   GHistThetaBinning::Fill(const Double_t value, const Double_t theta)
{
    if(theta==-1)
        return GHistScaCor::Fill(value);

    Int_t   channel = GetThetaBin(theta);
    if(channel<0)
        return 0;
    if(channel>=bin.GetEntriesFast())
        ExpandBin(channel+1);
    return ((GHistScaCor*)bin.At(channel))->Fill(value);
}

void    GHistThetaBinning::Reset(Option_t* option)
{
    GHistScaCor::Reset(option);
    bin.Clear();
}

void    GHistThetaBinning::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    GHistScaCor::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    for(int i=0; i<bin.GetEntriesFast(); i++)
        ((GHistScaCor*)bin.At(i))->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void	GHistThetaBinning::Scale(Double_t c1, Option_t* option)
{
    GHistScaCor::Scale(c1, option);
    for(int i=0; i<bin.GetEntriesFast(); i++)
        ((GHistScaCor*)bin.At(i))->Scale(c1, option);
}

void    GHistThetaBinning::PrepareWriteList(GHistWriteList* arr, const char *name)
{
    if(!arr)
        return;

    GHistScaCor::PrepareWriteList(arr, name);
    if(bin.GetEntriesFast()==0)
        return;

    GHistWriteList* ThetaBinning    = arr->GetDirectory(TString(GHThB_folderName));
    for(int i=0; i<bin.GetEntriesFast(); i++)
    {
        GHistWriteList* BinDir  = ThetaBinning->GetDirectory(TString(GHThB_chanelFolderName).Append(TString::Itoa(ThetaBinningRangeMin+(i*ThetaBinningWidth), 10)).Append("_").Append(TString::Itoa(ThetaBinningRangeMin+((i+1)*ThetaBinningWidth), 10)));
        if(name)
            ((GHistScaCor*)bin.At(i))->PrepareWriteList(BinDir, TString(name).Append(GHThB_binNameSuffix).Append(TString::Itoa(ThetaBinningRangeMin+(i*ThetaBinningWidth), 10)).Append("_").Append(TString::Itoa(ThetaBinningRangeMin+((i+1)*ThetaBinningWidth), 10)));
        else
            ((GHistScaCor*)bin.At(i))->PrepareWriteList(BinDir);
    }
}

Int_t   GHistThetaBinning::WriteWithoutCalcResult(const char* name, Int_t option, Int_t bufsize)
{
    Int_t res   = GHistThetaBinning::WriteWithoutCalcResult(name, option, bufsize);
    if(bin.GetEntriesFast()==0)
        return res;

    TDirectory* parentDir   = gDirectory;
    TDirectory* dir         = GetCreateDirectory(GHThB_folderName);

    TString nameBuffer;
    for(int i=0; i<bin.GetEntriesFast(); i++)
    {
        dir->cd();
        GetCreateDirectory(TString(GHThB_chanelFolderName).Append(TString::Itoa(ThetaBinningRangeMin+(i*ThetaBinningWidth), 10)).Append("_").Append(TString::Itoa(ThetaBinningRangeMin+((i+1)*ThetaBinningWidth), 10)).Data())->cd();

        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHThB_binNameSuffix);
            nameBuffer.Append(TString::Itoa(ThetaBinningRangeMin+(i*ThetaBinningWidth), 10)).Append("_").Append(TString::Itoa(ThetaBinningRangeMin+((i+1)*ThetaBinningWidth), 10));
            res += ((GHistThetaBinning*)bin.At(i))->WriteWithoutCalcResult(nameBuffer.Data(), option, bufsize);
        }
        else
            res += ((GHistThetaBinning*)bin.At(i))->WriteWithoutCalcResult(0, option, bufsize);
    }

    parentDir->cd();
    return  res;
}
