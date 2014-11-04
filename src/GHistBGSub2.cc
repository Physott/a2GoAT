#include "GHistBGSub2.h"
#include "GTreeTagger.h"

using namespace std;


GHistBGSub2::GHistBGSub2(Bool_t linkHistogram) :
    GHistBGSub(linkHistogram)
{
}

GHistBGSub2::GHistBGSub2() :
    GHistBGSub(kTRUE)
{
    result      = new GHistScaCor2();
    prompt      = new GHistScaCor2();
    randSum     = new GHistScaCor2();
}

GHistBGSub2::GHistBGSub2(const char* name, const char* title,
                         const Int_t nbinsx, const Double_t xlow, const Double_t xup,
                         const Int_t nbinsy, const Double_t ylow, const Double_t yup,
                         Bool_t linkHistogram) :
    GHistBGSub(linkHistogram)
{
    result      = new GHistScaCor2(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup, linkHistogram);
    prompt      = new GHistScaCor2(TString(name).Append(GHBS_promptNameSuffix), TString(title).Append(GHBS_promptTitleSuffix),
                                   nbinsx, xlow, xup,
                                   nbinsy, ylow, yup,
                                   kFALSE);
    randSum      = new GHistScaCor2(TString(name).Append(GHBS_randSumNameSuffix), TString(title).Append(GHBS_randSumTitleSuffix),
                                    nbinsx, xlow, xup,
                                    nbinsy, ylow, yup,
                                    kFALSE);
}

GHistBGSub2::~GHistBGSub2()
{
}

void    GHistBGSub2::CreateRandBin()
{
    GHistScaCor2*    hist_rand = new GHistScaCor2(TString(result->GetName()).Append(GHBS_randNameSuffix).Append(TString::Itoa(rand.GetEntriesFast(), 10)).Data(),
                                                  TString(result->GetTitle()).Append(GHBS_randTitleSuffix).Append(TString::Itoa(rand.GetEntriesFast(), 10)).Data(),
                                                  result->GetNbinsX(),
                                                  result->GetXmin(),
                                                  result->GetXmax(),
                                                  ((GHistScaCor2*)result)->GetNbinsY(),
                                                  ((GHistScaCor2*)result)->GetYmin(),
                                                  ((GHistScaCor2*)result)->GetYmax(),
                                                  kFALSE);
    rand.AddAtFree(hist_rand);
}

Int_t   GHistBGSub2::Fill(const Double_t x)
{
    std::cout << "ERROR: You tried to fill a 2 dim. GHistBGSub2 with only 1 value." << std::endl;
    return 0;
}

Int_t   GHistBGSub2::Fill(const Double_t x, const Double_t y, const Double_t taggerTime)
{
    if(taggerTime>=cutPromptMin && taggerTime<=cutPromptMax)
        return ((GHistScaCor2*)prompt)->Fill(x, y);
    for(int i=0; i<GetNRandCuts(); i++)
    {
        if(i>=rand.GetEntriesFast())
            ExpandRandBins(i+1);
        if(taggerTime>=cutRandMin[i] && taggerTime<=cutRandMax[i])
            return ((GHistScaCor2*)rand.At(i))->Fill(x, y);
    }
}

Int_t   GHistBGSub2::Fill(const Double_t x, const Double_t y, const GTreeTagger& tagger)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
        Fill(x, y, tagger.GetTagged_t(i));
}


GHistBGSub*    GHistBGSub2::ProjectionX(const char* name, Int_t firstybin, Int_t lastybin, Option_t* option)
{
    GHistBGSub*    ret = new GHistBGSub(name, name, result->GetNbinsX(), result->GetXmin(), result->GetXmax(), kFALSE);
    ret->result->Add(((TH2D*)result)->ProjectionX());
    ret->prompt->Add(((TH2D*)prompt)->ProjectionX());
    ret->randSum->Add(((TH2D*)randSum)->ProjectionX());
    for(int i=0; i<rand.GetEntriesFast(); i++)
        ((GHistScaCor*)rand.At(i))->Add((GHistScaCor*)rand.At(i)->ProjectionX());
    return ret;
}

GHistBGSub*   GHistBGSub2::ProjectionY(const char* name, Int_t firstxbin, Int_t lastxbin, Option_t* option)
{
    GHistBGSub*    ret = new GHistBGSub(name, name, result->GetNbinsY(), result->GetYaxis()->GetXmin(), result->GetYaxis()->GetXmax(), kFALSE);
    ret->result->Add(((TH2D*)result)->ProjectionY());
    ret->prompt->Add(((TH2D*)prompt)->ProjectionY());
    ret->randSum->Add(((TH2D*)randSum)->ProjectionY());
    for(int i=0; i<rand.GetEntriesFast(); i++)
        ((GHistScaCor*)rand.At(i))->Add((GHistScaCor*)rand.At(i)->ProjectionY());
    return ret;
}
