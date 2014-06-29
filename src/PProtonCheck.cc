#include "PProtonCheck.h"



PProtonCheck::PProtonCheck(const TString &_Name)    :
    name(_Name),
    histFoundProtons(TString(name).Append("_FoundProtons").Data(), TString(name).Append("  Number of Protons found").Data(),		5, 0, 5),
    histProtonAngleDiff(TString(name).Append("_ProtonAngleDiff").Data(), TString(name).Append("  Angle between meas. Proton and calc. Proton").Data(), 1800, 0, 180),
    histSmallestProtonAngleDiff(TString(name).Append("_smallestProtonAngleDiff").Data(), TString(name).Append("  smallest Angle between meas. Proton and calc. Proton").Data(), 1800, 0, 180),
    histCoplanarity(TString(name).Append("_Coplanarity").Data(), TString(name).Append("  Coplanarity").Data(), 3600, 0, 360)
{
    histProtonAngleDiff.SetCuts(-10, 5, -515, -15, 15, 510);
    histSmallestProtonAngleDiff.SetCuts(-10, 5, -515, -15, 15, 510);

    cutCoplanarity[0]   = 160;
    cutCoplanarity[1]   = 200;
}

PProtonCheck::~PProtonCheck()
{

}

Bool_t    PProtonCheck::ProcessEvent(const TLorentzVector &all, const GTreeParticle& protons, const GTreeTagger& tagger)
{
    foundProtons = 0;
    Double_t    smallestProtonAngleDiff;
    Double_t    help;
    for(int p=0; p<protons.GetNParticles(); p++)
    {
        smallestProtonAngleDiff = TMath::RadToDeg() * protons.Particle(p).Angle((tagger.GetVector(0)+TLorentzVector(0,0,0,MASS_PROTON) - all).Vect());
        smallestProtonAngleDiffIndex[p] = 0;
        for(int i=1; i<tagger.GetNTagged(); i++)
        {

            help = TMath::RadToDeg() * protons.Particle(p).Angle((tagger.GetVector(i)+TLorentzVector(0,0,0,MASS_PROTON) - all).Vect());
            histProtonAngleDiff.Fill(tagger.GetTagged_t(i), help);
            if(help < smallestProtonAngleDiff)
            {
                smallestProtonAngleDiff         = help;
                smallestProtonAngleDiffIndex[p] = i;
            }
        }
        histSmallestProtonAngleDiff.Fill(tagger.GetTagged_t(smallestProtonAngleDiffIndex[p]), smallestProtonAngleDiff);
        if(smallestProtonAngleDiff<4)
        {
            help    =   TMath::RadToDeg()*TMath::Abs(all.Phi() - protons.Particle(p).Phi());
            histCoplanarity.Fill(help);
            if(help>cutCoplanarity[0] && help<cutCoplanarity[1])
            {
                foundProtonsIndex[foundProtons] = p;
                foundProtons++;
            }
        }
    }
    histFoundProtons.Fill(foundProtons);

    if(foundProtons>0)
        return kTRUE;
    return kFALSE;
}

void	PProtonCheck::Write(TDirectory& curDir)
{
    curDir.cd();
    histFoundProtons.Write();
    histCoplanarity.Write();
    histProtonAngleDiff.Write(curDir);
    histSmallestProtonAngleDiff.Write(curDir);
}
