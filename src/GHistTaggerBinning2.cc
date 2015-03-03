#include "GHistTaggerBinning2.h"
#include "GTreeTagger.h"



#define GHTB_folderName         "TaggerBinning"
#define GHTB_binNameSuffix      "_Bins"
#define GHTB_binTitleSuffix     " Bins "



GHistTaggerBinning2::GHistTaggerBinning2(Bool_t linkHistogram) :
    GHistTaggerBinning(linkHistogram)
{
}

GHistTaggerBinning2::GHistTaggerBinning2() :
    GHistTaggerBinning(kTRUE)
{
    sum     = new GHistBGSub2();
    array   = new GHistBGSub3();
}

GHistTaggerBinning2::GHistTaggerBinning2(const char* name, const char* title,
                                       const Int_t nbinsx, const Double_t xlow, const Double_t xup,
                                       const Int_t nbinsy, const Double_t ylow, const Double_t yup,
                                       const Int_t nTaggerBins,
                                       const Bool_t linkHistogram) :
    GHistTaggerBinning(linkHistogram)
{
    sum     = new GHistBGSub2(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup, kFALSE);
    array   = new GHistBGSub3(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup, nTaggerBins, 0, nTaggerBins, kFALSE);
}

GHistTaggerBinning2::~GHistTaggerBinning2()
{

}

void    GHistTaggerBinning2::CalcResult()
{
    GHistBGSub2*    help    = ((GHistBGSub3*)array)->ProjectionXY();
    sum->Add(help);
    if(help)
        delete help;

    sum->CalcResult();
    array->CalcResult();
}

Int_t	GHistTaggerBinning2::Fill(const Double_t value)
{
    std::cout << "ERROR: You tried to fill a 2 dim. Hist with only 1 value." << std::endl;
    return 0;
}

Int_t   GHistTaggerBinning2::Fill(const Double_t x, const Double_t y, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        if(CreateHistogramsForTaggerBinning)
            Fill(x, y, tagger.GetTagged_ch(i), tagger.GetTagged_t(i));
        else
            Fill(x, y);
    }
}

