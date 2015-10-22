#ifndef __CINT__

#include "MRFitTaggerBins.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include <time.h>
#include <iostream>
#include "MRRecEff.h"


void TaggEff(Double_t* array, Double_t* darray)
{
    array[0]	= 48.0/100.0;
    array[1]	= 47.0/100.0;
    array[2]	= 47.0/100.0;
    array[3]	= 49.0/100.0;
    array[4]	= 51.0/100.0;
    array[5]	= 53.0/100.0;
    array[6]	= 54.0/100.0;
    array[7]	= 55.0/100.0;
    array[8]	= 56.0/100.0;
    array[9]	= 56.0/100.0;
    array[10]	= 58.0/100.0;
    array[11]	= 57.0/100.0;
    array[12]	= 59.0/100.0;
    array[13]	= 59.0/100.0;
    array[14]	= 60.0/100.0;
    array[15]	= 58.0/100.0;
    array[16]	= 61.0/100.0;
    array[17]	= 60.0/100.0;
    array[18]	= 60.0/100.0;
    array[19]	= 61.0/100.0;
    array[20]	= 63.0/100.0;
    array[21]	= 62.0/100.0;
    array[22]	= 64.0/100.0;
    array[23]	= 64.0/100.0;
    array[24]	= 66.0/100.0;
    array[25]	= 62.0/100.0;
    array[26]	= 65.0/100.0;
    array[27]	= 62.0/100.0;
    array[28]	= 66.0/100.0;
    array[29]	= 63.0/100.0;
    array[30]	= 65.0/100.0;
    array[31]	= 63.0/100.0;
    array[32]	= 65.0/100.0;
    array[33]	= 65.0/100.0;
    array[34]	= 65.0/100.0;
    array[35]	= 63.0/100.0;
    array[36]	= 66.0/100.0;
    array[37]	= 62.0/100.0;
    array[38]	= 68.0/100.0;
    array[39]	= 65.0/100.0;
    array[40]	= 66.0/100.0;
    array[41]	= 65.0/100.0;
    array[42]	= 65.0/100.0;
    array[43]	= 64.0/100.0;
    array[44]	= 66.0/100.0;
    array[45]	= 62.0/100.0;
    array[46]	= 66.0/100.0;
    array[47]	= 65.0/100.0;
    for(int i=0; i<48; i++)
        array[i] += 0.02;
    for(int i=0; i<48; i++)
        darray[i]	= 2.0/100.0;
}


double	MCBackgroundFactor(TFile* acquMcSignalFile, TFile* acquMcBGFile, TFile* out)
{
    TCanvas*	can = new TCanvas("canMCBackgroundFactor", "MCBackgroundFactor", 1500, 800);
    can->Divide(2, 1);

    TTree*	treeSignal	= (TTree*)acquMcSignalFile->Get("tagger");
    TTree*	treeBG		= (TTree*)acquMcBGFile->Get("tagger");

    std::cout << "calculating MCBackgroundFactor." << std::endl;
    std::cout << "acquMcSignalCount     acquMcBGCount     MCBackgroundFactor" << std::endl;

    can->cd(1);
    treeSignal->Draw("taggedChannel>>acquMcSignalCount");
    TH1D*	acquMcSignalCount = (TH1D*)gDirectory->Get("acquMcSignalCount");
    std::cout << acquMcSignalCount->GetBinContent(1) << "                  ";

    can->cd(2);
    treeBG->Draw("taggedChannel>>acquMcBGCount");
    TH1D*	acquMcBGCount = (TH1D*)gDirectory->Get("acquMcBGCount");
    std::cout << acquMcBGCount->GetBinContent(1) << "             ";

    out->cd();
    acquMcSignalCount->Write();
    acquMcBGCount->Write();
    can->Write();

    double	res	= (3*acquMcSignalCount->GetBinContent(1))/(0.082*acquMcBGCount->GetBinContent(1));
    std::cout << res << std::endl;
    return res;
}


/**
 * @brief the main routine
 * @param argc number of parameters
 * @param argv the parameters as strings
 * @return exit code
 */
int main(int argc, char *argv[])
{
    TFile*  out = TFile::Open("goatTrees/Result_CB.root", "RECREATE");
    TFile*  inAcquSimSignal = TFile::Open("goatTrees/Sim/Acqu_g4_sim_etap.root");
    TFile*  inAcquSimBG     = TFile::Open("goatTrees/Sim/Acqu_g4_sim_3pi0.root");

    double  backgroundFactor    = MCBackgroundFactor(inAcquSimSignal, inAcquSimBG, out);

    TFile*  inSimSignal = TFile::Open("goatTrees/Sim/Physics_g4_sim_etap.root");
    MRRecEff    recEff(inAcquSimSignal, inSimSignal);
    recEff.CalcResult();
    recEff.Draw(out);

    MRFitTaggerBins    resultSimSignal("simSignal", out, kRed);
    resultSimSignal.SetFile(inSimSignal, 3);
    resultSimSignal.RebinIM(5);
    resultSimSignal.Draw();
    resultSimSignal.FitGauss(kBlue, true);

    TFile*  inSimBG  = TFile::Open("goatTrees/Sim/Physics_g4_sim_3pi0.root");
    MRFitTaggerBins    resultSimBG("simBG", out, kGreen);
    resultSimBG.SetFile(inSimBG, 3);
    resultSimBG.Scale(backgroundFactor);
    resultSimBG.RebinIM(5);
    resultSimBG.Draw();
    resultSimBG.FitGauss(kBlue, false);
    //TCanvas*    can = resultSimSignal.GetCanvas()->Clone("SimBoth");

    resultSimBG.Draw(resultSimSignal.GetCanvas());

    MRFitTaggerBins    resultBoth("both", out, kBlue);
    resultBoth.SetFile(inSimSignal, 3);
    resultBoth.RebinIM(5);
    resultBoth.Add(resultSimBG);
    resultBoth.Draw();

    TFile*  in  = TFile::Open("goatTrees/Physics_CB.root");
    MRFitTaggerBins    result("data", out, kBlue);
    result.SetFile(in, 3);
    result.RebinIM(5);
    result.Draw();
    result.FitGauss(kBlack, resultSimSignal, resultSimBG);


    TCanvas*    can = new TCanvas("endresult", "endresult", 1500, 800);
    can->Divide(4, 3);
    can->cd(1);

    TH1*    etapCount   = result.GetResult();
    etapCount->Scale(1.0/5.0);
    etapCount->Draw();
    out->cd();
    etapCount->Write();
    can->cd(2);
    TH1*    etapCount2D = result.GetResult2D();
    etapCount2D->Scale(1.0/5.0);
    etapCount2D->Draw("COLZ");
    out->cd();
    etapCount2D->Write();

    can->cd(3);
    TH1*    etapCountCor = (TH1*)recEff.GetFactor();
    etapCountCor->Multiply(etapCount);
    etapCountCor->SetTitle("");
    etapCountCor->Draw();
    out->cd();
    etapCountCor->Write("etapCountCor");
    can->cd(4);
    TH2*    etapCount2DCor = (TH2*)recEff.GetFactor2D();
    etapCount2DCor->Multiply(etapCount2D);
    etapCount2DCor->Draw("COLZ");
    out->cd();
    etapCount2DCor->Write("etapCount2DCor");

    double  taggEff[48];
    double  dTaggEff[48];
    TaggEff(taggEff, dTaggEff);
    TH1*    teRaw  = new TH1D("taggEffRaw", "taggEffRaw", 48, 0, 48);
    for(int i=0; i<48; i++)
    {
        teRaw->SetBinContent(i+1, taggEff[i]);
        teRaw->SetBinError(i+1, dTaggEff[i]);
    }
    can->cd(5);
    teRaw->Draw();
    out->cd();
    teRaw->Write();
    TH1*    te  = (TH1*)teRaw->Clone("taggEff");
    te->RebinX(3);
    te->Scale(1.0/3.0);
    can->cd(6);
    te->Draw();
    out->cd();
    te->Write();

    TFile*  inScaler = TFile::Open("goatTrees/Scaler_CB.root");
    TH1*    scaler   = (TH1*)inScaler->Get("EPT_ScalerCorT");
    can->cd(9);
    scaler->Sumw2();
    scaler->Draw();
    out->cd();
    scaler->Write();

    can->cd(10);
    TH1*    scalerTeCor = (TH1*)scaler->Clone("scalerTeCor");
    //scalerTeCor->Sumw2();
    scalerTeCor->Multiply(teRaw);
    scalerTeCor->Draw();
    out->cd();
    scalerTeCor->Write();

    can->cd(11);
    TH1*    scalerTeCorRebined = (TH1*)scalerTeCor->Clone("scalerTeCorRebined");
    //scalerTeCorRebined->Sumw2();
    scalerTeCorRebined->Rebin(3);
    //scalerTeCorRebined->Scale(1.0/3.0);
    scalerTeCorRebined->Draw();
    out->cd();
    scalerTeCorRebined->Write();

    can->cd(12);
    TH1*    xSec = (TH1*)etapCountCor->Clone("xSec");
    //xSec->Sumw2();
    xSec->Divide(scalerTeCorRebined);
    xSec->Scale(1000000.0*4.2/0.08491);
    xSec->Draw();
    out->cd();
    xSec->Write();


    TCanvas*    thetacan = new TCanvas("thetaEndresult", "thetaEndresult", 1500, 800);
    thetacan->Divide(4, 4);
    for(int i=0; i<16; i++)
    {
        TH1*    hhh = (TH1*)etapCount2DCor->ProjectionY(TString("proj").Append(TString().Itoa(i, 10)), i+1, i+1)->Clone();
        double  integral    = hhh->Integral();
        hhh->Scale(xSec->GetBinContent(i+1)/integral);
        TString str("diffxSec");
        str.Append("_Tag_").Append(TString().Itoa(MRFitTaggerBins::TaggedEnergy[i*3], 10));
        str.Append("_to_").Append(TString().Itoa(MRFitTaggerBins::TaggedEnergy[(i*3)+3], 10));
        hhh->SetTitle(str);
        thetacan->cd(i+1);
        hhh->Draw();
    }

    TCanvas*    thetacanBoth = new TCanvas("thetaEndresultBoth", "thetaEndresultBoth", 1500, 800);
    thetacanBoth->Divide(4, 4);

    FILE*   otherResultsFile    = fopen("otherResults.dat", "r");
    char    str[4096];
    double  x[12][10];
    double  y[12][10];
    double  dy[12][10];
    int step = 11;
    while(!feof(otherResultsFile))
    {
        std::cout << step << std::endl;
        fgets(str, 4096, otherResultsFile); fgets(str, 4096, otherResultsFile); fgets(str, 4096, otherResultsFile);
        for(int i=0; i<10; i++)
        {
            fgets(str, 4096, otherResultsFile);
            sscanf(str, "%lf %lf %lf", &x[step][i], &y[step][i], &dy[step][i]);
            std::cout << x[step][i] << "   " << y[step][i] << "   " << dy[step][i] << std::endl;
        }
        fgets(str, 4096, otherResultsFile); fgets(str, 4096, otherResultsFile);
        step--;
    }

    for(int i=0; i<12; i++)
    {
        TH1*    hhh = (TH1*)etapCount2DCor->ProjectionY(TString("proj").Append(TString().Itoa(i, 10)), i+1, i+1)->Clone();
        double  integral    = hhh->Integral();
        hhh->Scale(xSec->GetBinContent(i+1)/integral);
        TString str("diffxSec");
        str.Append("_Tag_").Append(TString().Itoa(MRFitTaggerBins::TaggedEnergy[i*3], 10));
        str.Append("_to_").Append(TString().Itoa(MRFitTaggerBins::TaggedEnergy[(i*3)+3], 10));
        hhh->SetTitle(str);
        thetacanBoth->cd(i+1);
        hhh->Draw();

        TH1D*   hhhhh   = new TH1D(TString("otherResults").Append(TString().Itoa(i, 10)), TString("otherResults").Append(TString().Itoa(i, 10)), 10, -1, 1);
        for(int t=0; t<10; t++)
        {
            hhhhh->SetBinContent(t+1, y[i][t]);
            hhhhh->SetBinError(t+1, dy[i][t]);
        }
        thetacanBoth->cd(i+1);
        hhhhh->SetLineColor(kRed);
        hhhhh->Draw("SAME");
    }

    out->cd();
    result.GetResult()->Write();
    result.GetResult2D()->Write();
    can->Write();
    thetacan->Write();
    thetacanBoth->Write();

    return 0;
}

#endif
