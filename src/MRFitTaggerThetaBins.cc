#include "MRFitTaggerThetaBins.h"
#include "GHistManager.h"



MRFitTaggerThetaBins::MRFitTaggerThetaBins(const char* _Name, TFile* output, const int _Color) :
    name(_Name),
    out(output),
    color(_Color),
    can(0)
{
    for(int i=0; i<9; i++)
        thetaBins[i] = 0;
}

MRFitTaggerThetaBins::~MRFitTaggerThetaBins()
{
}

void    MRFitTaggerThetaBins::Add(MRFitTaggerThetaBins& origin)
{
    for(int i=0; i<9; i++)
    {
        thetaBins[i]->Add(origin.thetaBins[i]);
    }
}

bool	MRFitTaggerThetaBins::SetFile(TFile* _File, const int minTagger, const int maxTagger)
{
    TH3*		data	= (TH3*)_File->Get("TaggerBinning/IM");
    if(!data)
        return false;

    for(int i=0; i<9; i++)
    {
        if(thetaBins[i])	delete thetaBins[i];
        TString	str(name);
        str.Append("_Theta_").Append(TString().Itoa(i*20, 10));
        str.Append("_to_").Append(TString().Itoa((i+1)*20, 10));
        thetaBins[i]	= (TH1*)data->ProjectionX(str, (i*4)+1, (i+1)*4, minTagger, maxTagger);
    }

    return true;
}

void	MRFitTaggerThetaBins::RebinIM(const int addedBins)
{
    for(int i=0; i<9; i++)
    {
        thetaBins[i]->RebinX(addedBins);
        //thetaBins[i]->Scale(1/addedBins);
    }
}

void    MRFitTaggerThetaBins::FitGauss(const int _Color, const bool signal)
{
    TF1*	thetaFit[9];
    for(int i=0; i<9; i++)
    {
        TString str("fitfktGauss_");
        str.Append(thetaBins[i]->GetName());
        thetaFit[i]  = new TF1(str, "gaus(0)", 900, 1000);
        //std::cout << fit[i]->GetName() << std::endl;

        thetaFit[i]->SetLineColor(_Color);
        if(signal)
        {
            thetaFit[i]->SetParameters(thetaBins[i]->GetMaximum(), 957, 5);
            thetaFit[i]->SetParLimits(1, 950, 965);
            thetaFit[i]->SetParLimits(2, 1, 20);
        }
        else
        {
            thetaFit[i]->SetParameters(thetaBins[i]->GetMaximum(), 930, 50);
            thetaFit[i]->SetParLimits(1, 900, 1000);
            thetaFit[i]->SetParLimits(2, 25, 250);
        }
        thetaFit[i]->SetParLimits(0, thetaBins[i]->GetMaximum()*0.7, thetaBins[i]->GetMaximum()*1.5);

        thetaBins[i]->Fit(thetaFit[i], "R0");
        can->cd(i+1);
        thetaFit[i]->Draw("SAME");

        fitValuesBins[i].factor.value	= thetaFit[i]->GetParameter(0);
        fitValuesBins[i].factor.error	= thetaFit[i]->GetParError(0);
        fitValuesBins[i].mean.value     = thetaFit[i]->GetParameter(1);
        fitValuesBins[i].mean.error     = thetaFit[i]->GetParError(1);
        fitValuesBins[i].sigma.value	= thetaFit[i]->GetParameter(2);
        fitValuesBins[i].sigma.error	= thetaFit[i]->GetParError(2);
        fitValuesBins[i].count			= thetaBins[i]->Integral();
    }

    if(out)
    {
        out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
        can->Write();
    }
}

void    MRFitTaggerThetaBins::FitGauss(const int _Color, MRFitTaggerThetaBins& _FitValuesSignal, MRFitTaggerThetaBins& _FitValuesBG)
{
    TF1*	fit[9];

    for(int i=0; i<9; i++)
    {
        TString str("fitfktGauss_");
        str.Append(thetaBins[i]->GetName());
        fit[i]  = new TF1(str, "gaus(0)", 800, 1100);
        //std::cout << fit[i]->GetName() << std::endl;

        fit[i]->SetLineColor(_Color);
        fit[i]->SetParameters(_FitValuesSignal.fitValuesBins[i].factor.value, _FitValuesSignal.fitValuesBins[i].mean.value, _FitValuesSignal.fitValuesBins[i].sigma.value, _FitValuesBG.fitValuesBins[i].factor.value, _FitValuesBG.fitValuesBins[i].mean.value, _FitValuesBG.fitValuesBins[i].sigma.value);
        fit[i]->SetParLimits(0, 0, thetaBins[i]->GetMaximum()*1.5);
        fit[i]->SetParLimits(1, _FitValuesSignal.fitValuesBins[i].mean.value - _FitValuesSignal.fitValuesBins[i].mean.error, _FitValuesSignal.fitValuesBins[i].mean.value + _FitValuesSignal.fitValuesBins[i].mean.error);
        fit[i]->SetParLimits(2, _FitValuesSignal.fitValuesBins[i].sigma.value - _FitValuesSignal.fitValuesBins[i].sigma.error, _FitValuesSignal.fitValuesBins[i].sigma.value + _FitValuesSignal.fitValuesBins[i].sigma.error);
        fit[i]->SetParLimits(3, 0, thetaBins[i]->GetMaximum()*0.5);
        fit[i]->SetParLimits(4, _FitValuesBG.fitValuesBins[i].mean.value - _FitValuesBG.fitValuesBins[i].mean.error, _FitValuesBG.fitValuesBins[i].mean.value + _FitValuesBG.fitValuesBins[i].mean.error);
        fit[i]->SetParLimits(5, _FitValuesBG.fitValuesBins[i].sigma.value - _FitValuesBG.fitValuesBins[i].sigma.error, _FitValuesBG.fitValuesBins[i].sigma.value + _FitValuesBG.fitValuesBins[i].sigma.error);
        thetaBins[i]->Fit(fit[i], "R0");
        can->cd(i+1);
        fit[i]->Draw("SAME");

        fitValuesBins[i].factor.value	= fit[i]->GetParameter(0);
        fitValuesBins[i].factor.error	= fit[i]->GetParError(0);
        fitValuesBins[i].mean.value     = fit[i]->GetParameter(1);
        fitValuesBins[i].mean.error     = fit[i]->GetParError(1);
        fitValuesBins[i].sigma.value	= fit[i]->GetParameter(2);
        fitValuesBins[i].sigma.error	= fit[i]->GetParError(2);
        fitValuesBins[i].count			= fitValuesBins[i].factor.value * fitValuesBins[i].sigma.value * TMath::Sqrt(TMath::Pi());
        fitValuesBinsHelp[i].factor.value	= fit[i]->GetParameter(3);
        fitValuesBinsHelp[i].factor.error	= fit[i]->GetParError(3);
        fitValuesBinsHelp[i].mean.value     = fit[i]->GetParameter(4);
        fitValuesBinsHelp[i].mean.error     = fit[i]->GetParError(4);
        fitValuesBinsHelp[i].sigma.value	= fit[i]->GetParameter(5);
        fitValuesBinsHelp[i].sigma.error	= fit[i]->GetParError(5);
        double  help1           = fitValuesBins[i].factor.error * fitValuesBins[i].sigma.value * TMath::Sqrt(TMath::Pi());
        double  help2           = fitValuesBins[i].factor.value * fitValuesBins[i].sigma.error * TMath::Sqrt(TMath::Pi());
        fitValuesBinsHelp[i].count	= TMath::Sqrt((help1*help1) + (help2*help2));

        //thetaBins[i]->FitGauss(_Color);
    }

    if(out)
    {
        out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
        can->Write();
    }
}

void	MRFitTaggerThetaBins::Draw(TCanvas* _Canvas)
{
    if(!_Canvas)
    {
        if(can) delete can;
        can = new TCanvas(name, name, 1500, 800);
        can->Divide(TMath::Ceil(TMath::Sqrt(9)), TMath::Ceil(TMath::Sqrt(9)));

        for(int i=0; i<9; i++)
        {
            can->cd(i+1);
            thetaBins[i]->SetLineColor(color);
            thetaBins[i]->Draw();
            if(out)
            {
                out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
                thetaBins[i]->Write();
            }
        }

        if(out)
        {
            out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
            can->Write();
        }
        return;
    }

    can = _Canvas;
    for(int i=0; i<9; i++)
    {
        can->cd(i+1);
        thetaBins[i]->SetLineColor(color);
        thetaBins[i]->Draw("SAME");
        if(out)
        {
            out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
            thetaBins[i]->Write();
        }
    }

    if(out)
    {
        out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
        can->Write();
    }
}


void    MRFitTaggerThetaBins::Scale(const double factor)
{
    for(int i=0; i<9; i++)
        thetaBins[i]->Scale(factor);
}
