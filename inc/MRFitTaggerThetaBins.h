#ifndef __MRFitTaggerThetaBins_h__
#define __MRFitTaggerThetaBins_h__

#include <TMath.h>
#include <TFile.h>
#include <TString.h>
#include <TF1.h>
#include <TH3.h>
#include <TCanvas.h>


struct  FitValue
{
    double  value;
    double  error;
};
struct	FitValuesGauss
{
    FitValue	factor;
    FitValue	mean;
    FitValue	sigma;
    double      count;
};

class	MRFitTaggerThetaBins
{
private:
    TString         name;
    TFile*          out;
    int             color;
    TCanvas*        can;

    TH1*            thetaBins[40];
    FitValuesGauss  fitValuesBins[40];
    FitValuesGauss  fitValuesBinsHelp[40];

public:
    MRFitTaggerThetaBins(const char* _Name, TFile* output, const int _Color);
    ~MRFitTaggerThetaBins();

    void            Add(MRFitTaggerThetaBins& origin);
    TCanvas*        GetCanvas() {return can;}
    FitValuesGauss* GetFitValues()  {return fitValuesBins;}
    FitValuesGauss& GetFitValues(const int i)  {return fitValuesBins[i];}
    FitValuesGauss* GetFitValuesHelp()  {return fitValuesBinsHelp;}
    FitValuesGauss& GetFitValuesHelp(const int i)  {return fitValuesBinsHelp[i];}
    void            Draw(TCanvas* _Canvas = 0);
    void            FitGauss(const int _Color, const bool signal);
    void            FitGauss(const int _Color, MRFitTaggerThetaBins& _FitValuesSignal, MRFitTaggerThetaBins& _FitValuesBG);
    void            RebinIM(const int addedBins);
    void            Scale(const double factor);
    bool            SetFile(TFile* _File, const int minTagger, const int maxTagger);
};




#endif
