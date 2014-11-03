#include "GHistBGSub3.h"
#include "GTreeTagger.h"

using namespace std;


GHistBGSub3::GHistBGSub3(Bool_t linkHistogram) :
    GHistBGSub(linkHistogram)
{
}

GHistBGSub3::GHistBGSub3() :
    GHistBGSub(kTRUE)
{
    result      = new GHistScaCor3();
    prompt      = new GHistScaCor3();
    randSum     = new GHistScaCor3();
}

GHistBGSub3::GHistBGSub3(const char* name, const char* title,
                         const Int_t nbinsx, const Double_t xlow, const Double_t xup,
                         const Int_t nbinsy, const Double_t ylow, const Double_t yup,
                         const Int_t nbinsz, const Double_t zlow, const Double_t zup,
                         Bool_t linkHistogram) :
    GHistBGSub(linkHistogram)
{
    result      = new GHistScaCor3(name, title,
                                   nbinsx, xlow, xup,
                                   nbinsy, ylow, yup,
                                   nbinsz, zlow, zup,
                                   linkHistogram);
    prompt      = new GHistScaCor3(TString(name).Append(GHBS_promptNameSuffix), TString(title).Append(GHBS_promptTitleSuffix),
                                   nbinsx, xlow, xup,
                                   nbinsy, ylow, yup,
                                   nbinsz, zlow, zup,
                                   kFALSE);
    randSum     = new GHistScaCor3(TString(name).Append(GHBS_randSumNameSuffix), TString(title).Append(GHBS_randSumTitleSuffix),
                                    nbinsx, xlow, xup,
                                    nbinsy, ylow, yup,
                                    nbinsz, zlow, zup,
                                    kFALSE);
}

GHistBGSub3::~GHistBGSub3()
{
}

void    GHistBGSub3::CreateRandBin()
{
    GHistScaCor3*    hist_rand = new GHistScaCor3(TString(result->GetName()).Append(GHBS_randNameSuffix).Append(TString::Itoa(rand.GetEntriesFast(), 10)).Data(),
                                                  TString(result->GetTitle()).Append(GHBS_randTitleSuffix).Append(TString::Itoa(rand.GetEntriesFast(), 10)).Data(),
                                                  result->GetNbinsX(),
                                                  result->GetXmin(),
                                                  result->GetXmax(),
                                                  ((GHistScaCor3*)result)->GetNbinsY(),
                                                  ((GHistScaCor3*)result)->GetYmin(),
                                                  ((GHistScaCor3*)result)->GetYmax(),
                                                  ((GHistScaCor3*)result)->GetNbinsZ(),
                                                  ((GHistScaCor3*)result)->GetZmin(),
                                                  ((GHistScaCor3*)result)->GetZmax(),
                                                  kFALSE);
    rand.AddAtFree(hist_rand);
}

Int_t   GHistBGSub3::Fill(const Double_t x)
{
    std::cout << "ERROR: You tried to fill a 3 dim. GHistBGSub3 with only 1 value." << std::endl;
    return 0;
}

Int_t   GHistBGSub3::Fill(const Double_t x, const Double_t y)
{
    std::cout << "ERROR: You tried to fill a 3 dim. GHistBGSub3 with only 2 value." << std::endl;
    return 0;
}

Int_t   GHistBGSub3::Fill(const Double_t x, const Double_t y, const Double_t z, const Double_t taggerTime)
{
    if(taggerTime>=cutPromptMin && taggerTime<=cutPromptMax)
        return ((GHistScaCor3*)prompt)->Fill(x, y, z);
    for(int i=0; i<GetNRandCuts(); i++)
    {
        if(i>=rand.GetEntriesFast())
            ExpandRandBins(i+1);
        if(taggerTime>=cutRandMin[i] && taggerTime<=cutRandMax[i])
            return ((GHistScaCor3*)rand.At(i))->Fill(x, y, z);
    }
}

Int_t   GHistBGSub3::Fill(const Double_t x, const Double_t y, const Double_t z, const GTreeTagger& tagger)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
        Fill(x, y, z, tagger.GetTagged_t(i));
}

