#include "PProtonCheck.h"



PProtonCheck::PProtonCheck()    :
    histFoundProtons("FoundProtons", "Number of Protons found", 5, 0, 5),
    histProtonAngleDiff("ProtonAngleDiff", "Angle between meas. Proton and calc. Proton", 1800, 0, 180),
    histSmallestProtonAngleDiff("smallestProtonAngleDiff", "smallest Angle between meas. Proton and calc. Proton", 1800, 0, 180),
    histCoplanarity("Coplanarity", "Coplanarity", 3600, 0, 360)
{
    histProtonAngleDiff.SetCuts(-10, 5, -515, -15, 15, 510);
    histSmallestProtonAngleDiff.SetCuts(-10, 5, -515, -15, 15, 510);

    cutCoplanarity[0]   = 160;
    cutCoplanarity[1]   = 200;
}

PProtonCheck::~PProtonCheck()
{

}

Bool_t    PProtonCheck::ProcessEvent(const TLorentzVector &all)
{
    foundProtons = 0;
    Double_t    smallestProtonAngleDiff;
    Double_t    help;
    for(int p=0; p<protons->GetNParticles(); p++)
    {
        smallestProtonAngleDiff = TMath::RadToDeg()*protons->Particle(p).Angle((tagger->GetVector(0)+TLorentzVector(0,0,0,MASS_PROTON) - all).Vect());
        smallestProtonAngleDiffIndex[p] = 0;
        for(int i=1; i<tagger->GetNTagged(); i++)
        {

            help = TMath::RadToDeg()*protons->Particle(p).Angle((tagger->GetVector(i)+TLorentzVector(0,0,0,MASS_PROTON) - all).Vect());
            histProtonAngleDiff.Fill(tagger->GetTagged_t(i), help);
            if(help < smallestProtonAngleDiff)
            {
                smallestProtonAngleDiff         = help;
                smallestProtonAngleDiffIndex[p] = i;
            }
        }
        histSmallestProtonAngleDiff.Fill(tagger->GetTagged_t(smallestProtonAngleDiffIndex[p]), smallestProtonAngleDiff);
        if(smallestProtonAngleDiff<4)
        {
            help    =   TMath::RadToDeg()*TMath::Abs(all.Phi()-protons->Particle(p).Phi());
            histCoplanarity.Fill(help);
            if(help>cutCoplanarity[0] && help<cutCoplanarity[1])
            {
                foundProtonsIndex[foundProtons] = p;
                foundProtons++;
            }
        }
    }
    histFoundProtons.Fill(foundProtons);
}

void	PProtonCheck::Write()
{
    file_out->cd();
    TDirectory* curDir  = gDirectory->GetDirectory("WithProton");
    if(!curDir)
    {
        file_out->cd();
        gDirectory->mkdir("WithProton");
        curDir  = file_out->GetDirectory("WithProton");
    }
    curDir->cd();
    histFoundProtons.Write();
    histCoplanarity.Write();
    histProtonAngleDiff.Write(curDir);
    histSmallestProtonAngleDiff.Write(curDir);
}
