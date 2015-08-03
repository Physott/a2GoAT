#include "MRFitTaggerBins.h"
#include "GHistManager.h"
#include <TF1.h>



MRFitTaggerBins::MRFitTaggerBins(const char* _Name, TFile* output, const int _Color) :
    name(_Name),
    out(output),
    color(_Color),
    can(0),
    nBins(47),
    sum(0)
{
    for(int i=0; i<47; i++)
    {
        bins[i]         = 0;
        thetaBins[i]    = 0;
    }
}

MRFitTaggerBins::~MRFitTaggerBins()
{
}

void    MRFitTaggerBins::Add(MRFitTaggerBins& origin)
{
    if(origin.nBins!=nBins)
    {
        std::cout << "Can not add MRFitTaggerBins. Origin has " << origin.nBins << " bins and Destination has " << nBins <<"." << std::endl;
        return;
    }
    for(int i=0; i<nBins; i++)
    {
        bins[i]->Add(origin.bins[i]);
    }
    sum->Add(origin.sum);
}

bool	MRFitTaggerBins::SetFile(TFile* _File, const int summedTaggerChannels)
{
    TH3*		data	= (TH3*)_File->Get("TaggerBinning/IM");
    if(!data)
        return false;

    TString	str(name);
    str.Append("_TagAll");
    sum	= (TH1*)data->ProjectionX(str, 0, -1, 0, -1);

    int	overflow	= nBins%summedTaggerChannels;
    nBins   		= (nBins-overflow)/summedTaggerChannels;

    for(int i=0; i<nBins; i++)
    {
        if(bins[i])	delete bins[i];
        TString	str(name);
        str.Append("_Tag_").Append(TString().Itoa(TaggedEnergy[i*summedTaggerChannels], 10));
        str.Append("_to_").Append(TString().Itoa(TaggedEnergy[(i*summedTaggerChannels)+summedTaggerChannels], 10));
        bins[i]	= (TH1*)data->ProjectionX(str, 0, -1, (i*summedTaggerChannels)+1, (i*summedTaggerChannels)+summedTaggerChannels);

        str = name;
        str.Append("Theta");
        str.Append("_Tag_").Append(TString().Itoa(TaggedEnergy[i*summedTaggerChannels], 10));
        str.Append("_to_").Append(TString().Itoa(TaggedEnergy[(i*summedTaggerChannels)+summedTaggerChannels], 10));
        if(thetaBins[i])	delete thetaBins[i];
        thetaBins[i]    = new MRFitTaggerThetaBins(str, out, color);
        thetaBins[i]->SetFile(_File, (i*summedTaggerChannels)+1, (i*summedTaggerChannels)+summedTaggerChannels);

    }

    if(overflow)
    {
        TString	str(name);
        str.Append("_Tag_").Append(TString().Itoa(TaggedEnergy[nBins*summedTaggerChannels], 10));
        str.Append("_to_").Append(TString().Itoa(TaggedEnergy[(nBins*summedTaggerChannels)+overflow-1], 10));
        bins[nBins]	= (TH1*)data->ProjectionX(str, 0, -1, (nBins*summedTaggerChannels)+1, (nBins*summedTaggerChannels)+overflow-1);

        str = name;
        str.Append("Theta");
        str.Append("_Tag_").Append(TString().Itoa(TaggedEnergy[nBins*summedTaggerChannels], 10));
        str.Append("_to_").Append(TString().Itoa(TaggedEnergy[(nBins*summedTaggerChannels)+overflow-1], 10));
        if(thetaBins[nBins])	delete thetaBins[nBins];
        thetaBins[nBins]    = new MRFitTaggerThetaBins(str, out, color);
        thetaBins[nBins]->SetFile(_File, (nBins*summedTaggerChannels)+1, (nBins*summedTaggerChannels)+overflow-1);


        nBins++;
    }


    return true;
}

void	MRFitTaggerBins::RebinIM(const int addedBins)
{
    sum->RebinX(addedBins);
    for(int i=0; i<nBins; i++)
    {
        bins[i]->RebinX(addedBins);
        thetaBins[i]->RebinIM(addedBins);
    }
}

void    MRFitTaggerBins::FitGauss(const int _Color, const bool signal)
{
    TF1*	fit[nBins];
    TF1*	sumfit;

    TString str("fitfktGauss_");
    str.Append(sum->GetName());
    sumfit  = new TF1(str, "gaus(0)", 800, 1100);

    sumfit->SetLineColor(_Color);
    if(signal)
    {
        sumfit->SetParameters(sum->GetMaximum(), 957, 5);
        sumfit->SetParLimits(1, 950, 965);
        sumfit->SetParLimits(2, 1, 20);
    }
    else
    {
        sumfit->SetParameters(sum->GetMaximum(), 930, 50);
        sumfit->SetParLimits(1, 900, 1000);
        sumfit->SetParLimits(2, 25, 250);
    }
    sumfit->SetParLimits(0, sum->GetMaximum()*0.7, sum->GetMaximum()*1.5);
    sum->Fit(sumfit, "QR0");

    fitValuesSum.factor.value	= sumfit->GetParameter(0);
    fitValuesSum.factor.error	= sumfit->GetParError(0);
    fitValuesSum.mean.value     = sumfit->GetParameter(1);
    fitValuesSum.mean.error     = sumfit->GetParError(1);
    fitValuesSum.sigma.value	= sumfit->GetParameter(2);
    fitValuesSum.sigma.error	= sumfit->GetParError(2);
    fitValuesSum.count			= sum->Integral();

    for(int i=0; i<nBins; i++)
    {
        TString str("fitfktGauss_");
        str.Append(bins[i]->GetName());
        fit[i]  = new TF1(str, "gaus(0)", 800, 1100);
        //std::cout << fit[i]->GetName() << std::endl;

        fit[i]->SetLineColor(_Color);
        if(signal)
        {
            fit[i]->SetParameters(bins[i]->GetMaximum(), 957, 5);
            fit[i]->SetParLimits(1, 950, 965);
            fit[i]->SetParLimits(2, 1, 20);
        }
        else
        {
            fit[i]->SetParameters(bins[i]->GetMaximum(), 930, 50);
            fit[i]->SetParLimits(1, 900, 1000);
            fit[i]->SetParLimits(2, 25, 250);
        }
        fit[i]->SetParLimits(0, bins[i]->GetMaximum()*0.7, bins[i]->GetMaximum()*1.5);

        bins[i]->Fit(fit[i], "QR0");
        can->cd(i+1);
        fit[i]->Draw("SAME");

        fitValuesBins[i].factor.value	= fit[i]->GetParameter(0);
        fitValuesBins[i].factor.error	= fit[i]->GetParError(0);
        fitValuesBins[i].mean.value     = fit[i]->GetParameter(1);
        fitValuesBins[i].mean.error     = fit[i]->GetParError(1);
        fitValuesBins[i].sigma.value	= fit[i]->GetParameter(2);
        fitValuesBins[i].sigma.error	= fit[i]->GetParError(2);
        fitValuesBins[i].count			= bins[i]->Integral();

        thetaBins[i]->FitGauss(_Color, signal);
    }

    if(out)
    {
        out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
        can->Write();
        sum->Write();
        sumfit->Write();
    }
}

void    MRFitTaggerBins::FitGauss(const int _Color, MRFitTaggerBins& _FitValuesSignal, MRFitTaggerBins& _FitValuesBG)
{
    TF1*	fit[nBins];
    TF1*	sumfit;

    TString str("fitfktGaussForced_");
    str.Append(sum->GetName());
    sumfit  = new TF1(str, "gaus(0)+gaus(3)", 850, 1025);

    sumfit->SetLineColor(_Color);
    sumfit->SetParameters(_FitValuesSignal.fitValuesSum.factor.value, _FitValuesSignal.fitValuesSum.mean.value, _FitValuesSignal.fitValuesSum.sigma.value, _FitValuesBG.fitValuesSum.factor.value, _FitValuesBG.fitValuesSum.mean.value, _FitValuesBG.fitValuesSum.sigma.value);
    sumfit->SetParLimits(0, 0, sum->GetMaximum()*1.5);
    sumfit->SetParLimits(1, _FitValuesSignal.fitValuesSum.mean.value - _FitValuesSignal.fitValuesSum.mean.error, _FitValuesSignal.fitValuesSum.mean.value + _FitValuesSignal.fitValuesSum.mean.error);
    sumfit->SetParLimits(2, _FitValuesSignal.fitValuesSum.sigma.value - _FitValuesSignal.fitValuesSum.sigma.error, _FitValuesSignal.fitValuesSum.sigma.value + _FitValuesSignal.fitValuesSum.sigma.error);
    sumfit->SetParLimits(3, 0, sum->GetMaximum()*0.5);
    sumfit->SetParLimits(4, _FitValuesBG.fitValuesSum.mean.value - _FitValuesBG.fitValuesSum.mean.error, _FitValuesBG.fitValuesSum.mean.value + _FitValuesBG.fitValuesSum.mean.error);
    sumfit->SetParLimits(5, _FitValuesBG.fitValuesSum.sigma.value - _FitValuesBG.fitValuesSum.sigma.error, _FitValuesBG.fitValuesSum.sigma.value + _FitValuesBG.fitValuesSum.sigma.error);
    sum->Fit(sumfit, "QR0");

    fitValuesSum.factor.value	= sumfit->GetParameter(0);
    fitValuesSum.factor.error	= sumfit->GetParError(0);
    fitValuesSum.mean.value     = sumfit->GetParameter(1);
    fitValuesSum.mean.error     = sumfit->GetParError(1);
    fitValuesSum.sigma.value	= sumfit->GetParameter(2);
    fitValuesSum.sigma.error	= sumfit->GetParError(2);
    fitValuesSum.count			= fitValuesSum.factor.value * fitValuesSum.sigma.value * TMath::Sqrt(TMath::Pi());
    fitValuesSumHelp.factor.value	= sumfit->GetParameter(3);
    fitValuesSumHelp.factor.error	= sumfit->GetParError(3);
    fitValuesSumHelp.mean.value     = sumfit->GetParameter(4);
    fitValuesSumHelp.mean.error     = sumfit->GetParError(4);
    fitValuesSumHelp.sigma.value	= sumfit->GetParameter(5);
    fitValuesSumHelp.sigma.error	= sumfit->GetParError(5);
    double  help1           = fitValuesSum.factor.error * fitValuesSum.sigma.value * TMath::Sqrt(TMath::Pi());
    double  help2           = fitValuesSum.factor.value * fitValuesSum.sigma.error * TMath::Sqrt(TMath::Pi());
    fitValuesSumHelp.count	= TMath::Sqrt((help1*help1) + (help2*help2));

    for(int i=0; i<nBins; i++)
    {
        TString str("fitfktGauss_");
        str.Append(bins[i]->GetName());
        fit[i]  = new TF1(str, "gaus(0)", 800, 1100);
        //std::cout << fit[i]->GetName() << std::endl;

        fit[i]->SetLineColor(_Color);
        fit[i]->SetParameters(_FitValuesSignal.fitValuesBins[i].factor.value, _FitValuesSignal.fitValuesBins[i].mean.value, _FitValuesSignal.fitValuesBins[i].sigma.value, _FitValuesBG.fitValuesSum.factor.value, _FitValuesBG.fitValuesSum.mean.value, _FitValuesBG.fitValuesSum.sigma.value);
        fit[i]->SetParLimits(0, 0, bins[i]->GetMaximum()*1.5);
        fit[i]->SetParLimits(1, _FitValuesSignal.fitValuesBins[i].mean.value - _FitValuesSignal.fitValuesBins[i].mean.error, _FitValuesSignal.fitValuesBins[i].mean.value + _FitValuesSignal.fitValuesBins[i].mean.error);
        fit[i]->SetParLimits(2, _FitValuesSignal.fitValuesBins[i].sigma.value - _FitValuesSignal.fitValuesBins[i].sigma.error, _FitValuesSignal.fitValuesBins[i].sigma.value + _FitValuesSignal.fitValuesBins[i].sigma.error);
        fit[i]->SetParLimits(3, 0, bins[i]->GetMaximum()*0.5);
        fit[i]->SetParLimits(4, _FitValuesBG.fitValuesSum.mean.value - _FitValuesBG.fitValuesSum.mean.error, _FitValuesBG.fitValuesSum.mean.value + _FitValuesBG.fitValuesSum.mean.error);
        fit[i]->SetParLimits(5, _FitValuesBG.fitValuesSum.sigma.value - _FitValuesBG.fitValuesSum.sigma.error, _FitValuesBG.fitValuesSum.sigma.value + _FitValuesBG.fitValuesSum.sigma.error);
        bins[i]->Fit(fit[i], "QR0");
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

        thetaBins[i]->FitGauss(_Color, *_FitValuesSignal.thetaBins[i], *_FitValuesBG.thetaBins[i]);
    }

    if(out)
    {
        out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
        can->Write();
        sum->Write();
        sumfit->Write();
    }
}

void	MRFitTaggerBins::Draw(TCanvas *_Canvas)
{
    if(!_Canvas)
    {
        TString	str(name);
        str.Append("_Canvas");
        if(can) delete can;
        can = new TCanvas(str, name, 1500, 800);
        can->Divide(TMath::Ceil(TMath::Sqrt(nBins)), TMath::Ceil(TMath::Sqrt(nBins)));

        for(int i=0; i<nBins; i++)
        {
            can->cd(i+1);
            bins[i]->SetLineColor(color);
            bins[i]->Draw();
            thetaBins[i]->Draw();
            if(out)
            {
                out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
                bins[i]->Write();
            }
        }

        if(out)
        {
            out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
            can->Write();
            sum->Write();
        }
        return;
    }

    can = _Canvas;
    for(int i=0; i<nBins; i++)
    {
        can->cd(i+1);
        bins[i]->SetLineColor(color);
        bins[i]->Draw("SAME");
        thetaBins[i]->Draw();
        if(out)
        {
            out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
            bins[i]->Write();
        }
    }

    if(out)
    {
        out->cd(); GHistLinked::GetCreateDirectory(name)->cd();
        can->Write();
        sum->Write();
    }
}

void    MRFitTaggerBins::Scale(const double factor)
{
    for(int i=0; i<nBins; i++)
    {
        bins[i]->Scale(factor);
        thetaBins[i]->Scale(factor);
    }
}

TH1D*    MRFitTaggerBins::GetResult()
{
    TH1D*    result  = new TH1D("endresult", "endresult", 16, 0, 48);
    for(int i=0; i<16; i++)
    {
        result->SetBinContent(i+1, fitValuesBins[i].count);
        result->SetBinError(i+1, fitValuesBinsHelp[i].count);
    }
    return result;
}

TH2D*    MRFitTaggerBins::GetResult2D()
{
    TH2D*    result  = new TH2D("endresult2D", "endresult2D", 16, 0, 48, 10, -1, 1);
    for(int i=0; i<16; i++)
    {
        for(int t=0; t<10; t++)
        {
            if(thetaBins[i]->GetFitValues(t).count>0)
            {
                result->SetBinContent(i+1, t+1, thetaBins[i]->GetFitValues(t).count);
                result->SetBinError(i+1, t+1, thetaBins[i]->GetFitValuesHelp(t).count);
            }
            else
            {
                result->SetBinContent(i+1, t+1, 0);
                result->SetBinError(i+1, t+1, 1);
            }
        }
    }
    return result;
}

double	MRFitTaggerBins::TaggedEnergy[47] = {
1577.31,
1573.83,
1570.43,
1567.06,
1563.71,
1560.39,
1557.11,
1553.8,
1550.53,
1547.24,
1543.98,
1540.7,
1537.44,
1534.17,
1530.91,
1527.66,
1524.4,
1521.13,
1517.89,
1514.63,
1511.36,
1508.11,
1504.86,
1501.62,
1498.39,
1495.16,
1491.94,
1488.71,
1485.49,
1482.26,
1479.03,
1475.8,
1472.58,
1469.36,
1466.13,
1462.89,
1459.67,
1456.4,
1453.17,
1449.92,
1446.66,
1443.42,
1440.17,
1436.94,
1433.7,
1430.49,
1427.32};
