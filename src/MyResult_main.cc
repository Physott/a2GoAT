#ifndef __CINT__

#include "MRFitTaggerBins.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include <time.h>
#include <TTree.h>
#include <iostream>
#include <TH2D.h>
#include <TVector3.h>
#include <TLorentzVector.h>



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


TH2* ReconstructEff(TFile* acquSignalFile, TFile* SignalFile, TFile* out)
{
    TTree*	tagger	= (TTree*)acquSignalFile->Get("tagger");
    TTree*	tracks	= (TTree*)acquSignalFile->Get("tracks");
    TTree*	pluto	= (TTree*)acquSignalFile->Get("h12");

    Int_t       nTagged;
    Int_t       taggedChannel[128];
    Double_t    taggedEnergy[128];
    Int_t       nTracks;
    Float_t     plab[20];
    Float_t     dircos[20][3];
    Float_t     elab[20];

    tagger->SetBranchAddress("nTagged", &nTagged);
    tagger->SetBranchAddress("taggedChannel", taggedChannel);
    tagger->SetBranchAddress("taggedEnergy", taggedEnergy);
    tracks->SetBranchAddress("nTracks",&nTracks);
    pluto->SetBranchAddress("plab", &plab);
    pluto->SetBranchAddress("dircos", &dircos);
    pluto->SetBranchAddress("elab", &elab);

    TH2D*	count	= new TH2D("count", "count", 48, 0, 48, 36, 0, 180);
    TH2D*	count6	= new TH2D("count6", "count6", 48, 0, 48, 36, 0, 180);
    TH2D*	count7	= new TH2D("count7", "count7", 48, 0, 48, 36, 0, 180);

    for(int i=0; i<tracks->GetEntriesFast(); i++)
    {
        tagger->GetEntry(i);
        tracks->GetEntry(i);
        pluto->GetEntry(i);

        TVector3 p(dircos[2][0], dircos[2][1], dircos[2][2]);
        p *= plab[2];
        TLorentzVector lv( p, elab[2]);
        TLorentzVector tag( 0, 0, taggedEnergy[0], taggedEnergy[0] + 938.272046);
        lv.Boost(-tag.BoostVector());
        double  etapTheta = lv.Theta() * TMath::RadToDeg();
        //std::cout << etapTheta << std::endl;

        for(int t=0; t<nTagged; t++)
        {
            count->Fill(taggedChannel[t], etapTheta);
            if(nTracks==6)
                count6->Fill(taggedChannel[t], etapTheta);
            else if(nTracks==7)
                count7->Fill(taggedChannel[t], etapTheta);
        }
    }

    TCanvas*	can = new TCanvas("canReconstructionEfficiency", "ReconstructionEfficiency", 1500, 800);
    can->Divide(3, 4);

    can->cd(1);
    count->Sumw2();
    count->Scale(1);
    count->RebinX(3);
    count->RebinY(4);
    count->SetTitle("all simulated #eta' events 6 Hits");
    count->GetXaxis()->SetTitle("tagger channel");
    count->SetStats(0);
    count->Draw("COLZ");
    can->cd(2);
    count6->Sumw2();
    count6->Scale(1);
    count6->RebinX(3);
    count6->RebinY(4);
    count6->SetTitle("all simulated #eta' events 6 Hits");
    count6->GetXaxis()->SetTitle("tagger channel");
    count6->SetStats(0);
    count6->Draw("COLZ");
    can->cd(3);
    count7->Sumw2();
    count7->Scale(1);
    count7->RebinX(3);
    count7->RebinY(4);
    count7->SetTitle("all simulated #eta' events 7 Hits");
    count7->GetXaxis()->SetTitle("tagger channel");
    count7->SetStats(0);
    count7->Draw("COLZ");

    TH3*		dataRaw	= (TH3*)SignalFile->Get("TaggerBinning/IM");
    TH2*        data	= (TH2*)dataRaw->Project3D("yze");

    can->cd(6);
    data->SetTitle("accepted simulated #eta' events 7 Hits");
    data->GetXaxis()->SetTitle("tagger channel");
    data->SetStats(0);
    data->RebinX(3);
    data->RebinY(4);
    data->Draw("COLZ");

    can->cd(4);
    TH1*        dataProjectionAcqu	= count7->ProjectionY();
    dataProjectionAcqu->SetTitle("accepted simulated #eta' events 7 Hits");
    dataProjectionAcqu->GetXaxis()->SetTitle("tagger channel");
    dataProjectionAcqu->SetStats(0);
    dataProjectionAcqu->Draw();

    can->cd(5);
    TH1*        dataProjection	= data->ProjectionY();
    dataProjection->SetTitle("accepted simulated #eta' events 7 Hits");
    dataProjection->GetXaxis()->SetTitle("tagger channel");
    dataProjection->SetStats(0);
    dataProjection->Draw();


    can->cd(7);
    TF1*        sinusFkt	= new TF1("sinusFkt", "sin(x*1.74532925199432955e-02)", 0, 180);
    sinusFkt->SetTitle("accepted simulated #eta' events 7 Hits");
    sinusFkt->GetXaxis()->SetTitle("tagger channel");
    sinusFkt->Draw();

    can->cd(8);
    TH1*        dataProjectionHelp	= (TH1*)dataProjection->Clone();
    dataProjectionHelp->SetTitle("accepted simulated #eta' events 7 Hits");
    dataProjectionHelp->GetXaxis()->SetTitle("tagger channel");
    dataProjectionHelp->Divide(sinusFkt);
    dataProjectionHelp->Draw();

    can->cd(9);
    TH2*        dataResult7	= (TH2*)count7->Clone();
    dataResult7->Divide(data);
    dataResult7->SetTitle("accepted simulated #eta' events 7 Hits");
    dataResult7->GetXaxis()->SetTitle("tagger channel");
    dataResult7->SetStats(0);
    dataResult7->Draw("COLZ");

    can->cd(10);
    TH2*        dataResult	= (TH2*)count->Clone();
    dataResult->Divide(data);
    dataResult->SetTitle("accepted simulated #eta' events 7 Hits");
    dataResult->GetXaxis()->SetTitle("tagger channel");
    dataResult->SetStats(0);
    dataResult->Draw("COLZ");


    can->cd(11);
    TH1*        dataResultProjectionY	= dataResult->ProjectionY();
    dataResultProjectionY->SetTitle("accepted simulated #eta' events 7 Hits");
    dataResultProjectionY->GetXaxis()->SetTitle("tagger channel");
    dataResultProjectionY->SetStats(0);
    dataResultProjectionY->Draw("COLZ");

    can->cd(12);
    TH1*        dataResultProjectionX	= dataResult->ProjectionX();
    dataResultProjectionX->SetTitle("accepted simulated #eta' events 7 Hits");
    dataResultProjectionX->GetXaxis()->SetTitle("tagger channel");
    dataResultProjectionX->SetStats(0);
    dataResultProjectionX->Draw("COLZ");

    out->cd();
    count->Write();
    count6->Write();
    count7->Write();
    dataResult->Write();
    can->Write();

    return dataResult;
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

    TFile*  inSimSignal  = TFile::Open("goatTrees/Sim/Physics_g4_sim_etap.root");
    TH2*    reconstructEfficiency   = ReconstructEff(inAcquSimSignal, inSimSignal, out);

    MRFitTaggerBins    resultSimSignal("simSignal", out, kRed);
    resultSimSignal.SetFile(inSimSignal, 3);
    //resultSimSignal.RebinIM(2);
    resultSimSignal.Draw();
    resultSimSignal.FitGauss(kBlue, true);


    TFile*  inSimBG  = TFile::Open("goatTrees/Sim/Physics_g4_sim_3pi0.root");
    MRFitTaggerBins    resultSimBG("simBG", out, kGreen);
    resultSimBG.SetFile(inSimBG, 3);
    resultSimBG.Scale(backgroundFactor);
    //resultSimBG.RebinIM(5);
    resultSimBG.Draw();
    resultSimBG.FitGauss(kBlue, false);
    //TCanvas*    can = resultSimSignal.GetCanvas()->Clone("SimBoth");

    resultSimBG.Draw(resultSimSignal.GetCanvas());

    MRFitTaggerBins    resultBoth("both", out, kBlue);
    resultBoth.SetFile(inSimSignal, 3);
    //resultBoth.RebinIM(5);
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
    can->cd(2);
    TH1*    etapCount2D = result.GetResult2D();
    etapCount2D->Scale(1.0/5.0);
    etapCount2D->Draw("COLZ");
    can->cd(3);
    TH1*    rrr2 = reconstructEfficiency->ProjectionX("rrr2");
    rrr2->Multiply(etapCount);
    rrr2->Draw();
    can->cd(4);
    TH2*    rrr1 = (TH2*)reconstructEfficiency->Clone("rrr1");
    rrr1->Multiply(etapCount2D);
    rrr1->Draw("COLZ");

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
    TH1*    te  = (TH1*)teRaw->Clone("taggEff");
    te->RebinX(3);
    te->Scale(1.0/3.0);
    can->cd(6);
    te->Draw();

    TFile*  inScaler = TFile::Open("goatTrees/Scaler_CB.root");
    TH1*    scaler   = (TH1*)inScaler->Get("EPT_ScalerCorT");
    can->cd(9);
    scaler->Sumw2();
    scaler->Draw();

    can->cd(10);
    TH1*    scalerTeCor = (TH1*)scaler->Clone("scalerTeCor");
    //scalerTeCor->Sumw2();
    scalerTeCor->Multiply(teRaw);
    scalerTeCor->Draw();

    can->cd(11);
    TH1*    scalerTeCorRebined = (TH1*)scalerTeCor->Clone("scalerTeCor");
    //scalerTeCorRebined->Sumw2();
    scalerTeCorRebined->Rebin(3);
    //scalerTeCorRebined->Scale(1.0/3.0);
    scalerTeCorRebined->Draw();

    can->cd(12);
    TH1*    xSec = (TH1*)rrr2->Clone("xSec");
    //xSec->Sumw2();
    xSec->Divide(scalerTeCorRebined);
    xSec->Draw();


    TCanvas*    thetacan = new TCanvas("thetaEndresult", "thetaEndresult", 1500, 800);
    thetacan->Divide(4, 4);
    for(int i=0; i<15; i++)
    {
        TH1*    hhh = (TH1*)rrr1->ProjectionY(TString("proj").Append(TString().Itoa(i, 10)), i+1, i+1)->Clone();
        double  integral    = hhh->Integral();
        hhh->Scale(xSec->GetBinContent(i+1)/integral);
        thetacan->cd(i+1);
        hhh->Draw();
    }

    out->cd();
    result.GetResult()->Write();
    result.GetResult2D()->Write();
    can->Write();
    thetacan->Write();

    return 0;
}

#endif
