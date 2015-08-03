#ifndef __MRFitTaggerBins_h__
#define __MRFitTaggerBins_h__

#include <TMath.h>
#include <TFile.h>
#include <TString.h>
#include <TH3.h>
#include <TH2D.h>
#include <TCanvas.h>

#include "MRFitTaggerThetaBins.h"


class	MRFitTaggerBins
{
private:
    TString                 name;
    TFile*                  out;
    int                     color;
    TCanvas*                can;
    int                     nBins;
    int                     nThetaBins;
    TH1*                    bins[47];
    MRFitTaggerThetaBins*   thetaBins[47];
    TH1*                    sum;
    FitValuesGauss          fitValuesBins[47];
    FitValuesGauss          fitValuesBinsHelp[47];
    FitValuesGauss          fitValuesSum;
    FitValuesGauss          fitValuesSumHelp;

public:
    MRFitTaggerBins(const char* _Name, TFile *output, const int _Color);
    ~MRFitTaggerBins();

    static  double          TaggedEnergy[47];

    void            Add(MRFitTaggerBins& origin);
    TCanvas*        GetCanvas() {return can;}
    FitValuesGauss* GetFitValues()  {return fitValuesBins;}
    FitValuesGauss& GetFitValues(const int i)  {return fitValuesBins[i];}
    FitValuesGauss& GetFitValuesSum()  {return fitValuesSum;}
    FitValuesGauss* GetFitValuesHelp()  {return fitValuesBinsHelp;}
    FitValuesGauss& GetFitValuesHelp(const int i)  {return fitValuesBinsHelp[i];}
    FitValuesGauss& GetFitValuesSumHelp()  {return fitValuesSumHelp;}
    FitValuesGauss* GetThetaFitValues(const int i)  {return thetaBins[i]->GetFitValues();}
    int             GetNBins();
    TH1D*           GetResult();
    TH2D*           GetResult2D();
    void            Draw(TCanvas* _Canvas = 0);
    void            FitGauss(const int _Color, const bool signal);
    void            FitGauss(const int _Color, MRFitTaggerBins& _FitValuesSignal, MRFitTaggerBins& _FitValuesBG);
    void            RebinIM(const int addedBins);
    void            Scale(const double factor);
    bool            SetFile(TFile* _File, const int summedTaggerChannels = 1);
};




#endif
