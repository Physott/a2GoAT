#include "MRRecEff.h"


MRRecEff::MRRecEff(TFile *acquSignalFile, TFile *SignalFile)
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

    acquCount       = new TH1D("acquCount", "acquCount", 48, 0, 48);
    acquCount6      = new TH1D("acquCount6", "acquCount6", 48, 0, 48);
    acquCount7      = new TH1D("acquCount7", "acquCount7", 48, 0, 48);
    acquCount2D     = new TH2D("acquCount2D", "acquCount2D", 48, 0, 48, 40, -1, 1);
    acquCount2D6	= new TH2D("acquCount2D6", "acquCount2D6", 48, 0, 48, 40, -1, 1);
    acquCount2D7	= new TH2D("acquCount2D7", "acquCount2D7", 48, 0, 48, 40, -1, 1);

    acquCount->Sumw2();
    acquCount6->Sumw2();
    acquCount7->Sumw2();
    acquCount2D->Sumw2();
    acquCount2D6->Sumw2();
    acquCount2D7->Sumw2();

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
        double  etapCosTheta = TMath::Cos(lv.Theta());
        //std::cout << etapTheta << std::endl;

        for(int t=0; t<nTagged; t++)
        {
            acquCount->Fill(taggedChannel[t]);
            acquCount2D->Fill(taggedChannel[t], etapCosTheta);
            if(nTracks==6)
            {
                acquCount6->Fill(taggedChannel[t]);
                acquCount2D6->Fill(taggedChannel[t], etapCosTheta);
            }
            else if(nTracks==7)
            {
                acquCount7->Fill(taggedChannel[t]);
                acquCount2D7->Fill(taggedChannel[t], etapCosTheta);
            }
        }
    }

    TH3*    acceptedCount7Raw   = (TH3*)SignalFile->Get("TaggerBinning/IM");
            acceptedCount7      = (TH1D*)acceptedCount7Raw->Project3D("ze");
            acceptedCount2D7    = (TH2D*)acceptedCount7Raw->Project3D("yze");

    acquCount->RebinX(3);
    acquCount6->RebinX(3);
    acquCount7->RebinX(3);
    acquCount2D->RebinX(3);
    acquCount2D6->RebinX(3);
    acquCount2D7->RebinX(3);
    acceptedCount7->RebinX(3);
    acceptedCount2D7->RebinX(3);

    acquCount2D->RebinY(4);
    acquCount2D6->RebinY(4);
    acquCount2D7->RebinY(4);
    acceptedCount2D7->RebinY(4);
}

MRRecEff::~MRRecEff()
{

}

void    MRRecEff::CalcResult()
{
    factor	= (TH1D*)acquCount->Clone("factor");
    factor->Divide(acceptedCount7);
    factor7	= (TH1D*)acquCount7->Clone("factor7");
    factor7->Divide(acceptedCount7);

    recEff	= (TH1D*)acceptedCount7->Clone("recEff");
    recEff->Divide(acquCount);
    recEff7	= (TH1D*)acceptedCount7->Clone("recEff7");
    recEff7->Divide(acquCount7);

    factor2D	= (TH2D*)acquCount2D->Clone("factor2D");
    factor2D->Divide(acceptedCount2D7);
    factor2D7	= (TH2D*)acquCount2D7->Clone("factor2D7");
    factor2D7->Divide(acceptedCount2D7);

    recEff2D	= (TH2D*)acceptedCount2D7->Clone("recEff7");
    recEff2D->Divide(acquCount2D);
    recEff2D7	= (TH2D*)acceptedCount2D7->Clone("recEff2D7");
    recEff2D7->Divide(acquCount2D7);
}

void    MRRecEff::Draw(TFile *out)
{
    TCanvas*	can = new TCanvas("canReconstructionEfficiency", "ReconstructionEfficiency", 1500, 800);
    can->Divide(4, 4);

    can->cd(1);
    acquCount->SetTitle("all simulated #eta' events");
    acquCount->GetXaxis()->SetTitle("tagger channel");
    acquCount->SetStats(0);
    acquCount->Draw("COLZ");
    can->cd(2);
    acquCount6->SetTitle("all simulated #eta' events 6 Hits");
    acquCount6->GetXaxis()->SetTitle("tagger channel");
    acquCount6->SetStats(0);
    acquCount6->Draw("COLZ");
    can->cd(3);
    acquCount7->SetTitle("all simulated #eta' events 7 Hits");
    acquCount7->GetXaxis()->SetTitle("tagger channel");
    acquCount7->SetStats(0);
    acquCount7->Draw("COLZ");

    can->cd(4);
    acceptedCount7->SetTitle("accepted simulated #eta' events 7 Hits");
    acceptedCount7->GetXaxis()->SetTitle("tagger channel");
    acceptedCount7->SetStats(0);
    acceptedCount7->Draw("COLZ");

    can->cd(5);
    factor->SetTitle("factor to calc real #eta'");
    factor->GetXaxis()->SetTitle("tagger channel");
    factor->SetStats(0);
    factor->Draw("COLZ");
    can->cd(6);
    factor7->SetTitle("factor to calc real #eta' only 7 Hits");
    factor7->GetXaxis()->SetTitle("tagger channel");
    factor7->SetStats(0);
    factor7->Draw("COLZ");

    can->cd(7);
    recEff->SetTitle("rec Eff");
    recEff->GetXaxis()->SetTitle("tagger channel");
    recEff->SetStats(0);
    recEff->Draw("COLZ");
    can->cd(8);
    recEff7->SetTitle("rec Eff only 7 Hits");
    recEff7->GetXaxis()->SetTitle("tagger channel");
    recEff7->SetStats(0);
    recEff7->Draw("COLZ");






    can->cd(9);
    acquCount2D->SetTitle("all simulated #eta' events");
    acquCount2D->GetXaxis()->SetTitle("tagger channel");
    acquCount2D->GetYaxis()->SetTitle("cos(#theta)");
    acquCount2D->SetStats(0);
    acquCount2D->Draw("COLZ");
    can->cd(10);
    acquCount2D6->SetTitle("all simulated #eta' events 6 Hits");
    acquCount2D6->GetXaxis()->SetTitle("tagger channel");
    acquCount2D6->GetYaxis()->SetTitle("cos(#theta)");
    acquCount2D6->SetStats(0);
    acquCount2D6->Draw("COLZ");
    can->cd(11);
    acquCount2D7->SetTitle("all simulated #eta' events 7 Hits");
    acquCount2D7->GetXaxis()->SetTitle("tagger channel");
    acquCount2D7->GetYaxis()->SetTitle("cos(#theta)");
    acquCount2D7->SetStats(0);
    acquCount2D7->Draw("COLZ");

    can->cd(12);
    acceptedCount2D7->SetTitle("accepted simulated #eta' events 7 Hits");
    acceptedCount2D7->GetXaxis()->SetTitle("tagger channel");
    acceptedCount2D7->GetYaxis()->SetTitle("cos(#theta)");
    acceptedCount2D7->SetStats(0);
    acceptedCount2D7->Draw("COLZ");

    can->cd(13);
    factor2D->SetTitle("factor to calc real #eta'");
    factor2D->GetXaxis()->SetTitle("tagger channel");
    factor2D->GetYaxis()->SetTitle("cos(#theta)");
    factor2D->SetStats(0);
    factor2D->Draw("COLZ");
    can->cd(14);
    factor2D7->SetTitle("factor to calc real #eta' only 7 Hits");
    factor2D7->GetXaxis()->SetTitle("tagger channel");
    factor2D7->GetYaxis()->SetTitle("cos(#theta)");
    factor2D7->SetStats(0);
    factor2D7->Draw("COLZ");

    can->cd(15);
    recEff2D->SetTitle("rec Eff");
    recEff2D->GetXaxis()->SetTitle("tagger channel");
    recEff2D->GetYaxis()->SetTitle("cos(#theta)");
    recEff2D->SetStats(0);
    recEff2D->Draw("COLZ");
    can->cd(16);
    recEff2D7->SetTitle("rec Eff only 7 Hits");
    recEff2D7->GetXaxis()->SetTitle("tagger channel");
    recEff2D7->GetYaxis()->SetTitle("cos(#theta)");
    recEff2D7->SetStats(0);
    recEff2D7->Draw("COLZ");







    out->cd();
    acquCount->Write();
    acquCount6->Write();
    acquCount7->Write();
    acceptedCount7->Write();
    can->Write();
}
