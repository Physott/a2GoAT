#include "GAnalysis3Mesons.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"


GAnalysis3Mesons::GAnalysis3Mesons(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(_IsEtap),
    fit4("fit4", kFALSE),
    hist_fit4(TString(name).Append("hist_fit4"), TString(title).Append("hist_fit4"), kFALSE),
    hist_fit4_SubAll(TString(name).Append("fit4_SubAll"), TString(title).Append("fit4_SubAll"), 1000, 0, 1000, kFALSE)
{
    fit4.AddConstraintMM();
    fit4.AddConstraintsIM();
}

GAnalysis3Mesons::~GAnalysis3Mesons()
{

}

void   GAnalysis3Mesons::CalcResult()
{
    fit4.CalcResult();
    hist_fit4.CalcResult();
    hist_fit4_SubAll.CalcResult();
}

void    GAnalysis3Mesons::Fill(const GTreeMeson& meson, GTreeParticle& photons, const GTreeTagger& tagger, const GTreeA2Geant& geantTree)
{
    success = false;

    Double_t    im  = meson.Particle(0).M();
    Double_t    theta  = meson.Particle(0).Theta()*TMath::RadToDeg();
    Double_t    phi  = meson.Particle(0).Phi()*TMath::RadToDeg();
    Double_t    mm;
    Double_t    sub_im_0    = (photons.Particle(0) + photons.Particle(1)).M();
    Double_t    sub_im_1    = (photons.Particle(2) + photons.Particle(3)).M();
    Double_t    sub_im_2    = (photons.Particle(4) + photons.Particle(5)).M();

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
        TLorentzVector  helpCM(meson.Particle(0));
        helpCM.Boost(-tagger.GetVectorProtonTarget(i).BoostVector());

        fit4.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
        fit4.SetBeam(tagger.GetTaggedEnergy(i));
        if(fit4.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i), geantTree))
        {
            success = true;
            hist_fit4.Fill(im, mm, theta, phi, helpCM.Theta()*TMath::RadToDeg(), sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            for(int k=0; k<6; k++)
            {
                for(int l=k+1; l<6; l++)
                    hist_fit4_SubAll.Fill((photons.Particle(k)+photons.Particle(l)).M(), tagger.GetTaggedTime(i));

                photons.AddParticle(fit4.GetFittedPhoton(k).Ek, fit4.GetFittedPhoton(k).Theta * TMath::RadToDeg(), fit4.GetFittedPhoton(k).Phi * TMath::RadToDeg(),
                                    0.0, photons.GetTime(k), photons.GetClusterSize(k), photons.GetCentralCrystal(k), photons.GetCentralVeto(k), photons.GetDetectors(k),
                                    photons.GetVetoEnergy(k), photons.GetMWPC0Energy(k), photons.GetMWPC1Energy(k), photons.GetTrackIndex(k));
            }
        }
    }
}

void    GAnalysis3Mesons::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* folder  = arr->GetDirectory("WithoutProton");

    fit4.PrepareWriteList(folder);
    hist_fit4.PrepareWriteList(folder, "fit4");
    hist_fit4_SubAll.PrepareWriteList(folder, "SubAll");
}

void    GAnalysis3Mesons::Reset(Option_t* option)
{
    fit4.Reset(option);
    hist_fit4.Reset(option);
    hist_fit4_SubAll.Reset(option);
}











GAnalysis3MesonsProton::GAnalysis3MesonsProton(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(_IsEtap),
    fitProton6("fitProton6", kFALSE),
    hist_fitProton6(TString(name).Append("hist_fitProton6"), TString(title).Append("hist_fitProton6"), kFALSE),
    hist_fitProton6_SubAll(TString(name).Append("fitProton6_SubAll"), TString(title).Append("fitProton6_SubAll"), 1000, 0, 1000, kFALSE)
{
    fitProton6.AddConstraintsTotMomentum();
    fitProton6.AddConstraintsTotEnergy();
    fitProton6.AddConstraintsIM();
}

GAnalysis3MesonsProton::~GAnalysis3MesonsProton()
{

}

void   GAnalysis3MesonsProton::CalcResult()
{
    fitProton6.CalcResult();
    hist_fitProton6.CalcResult();
    hist_fitProton6_SubAll.CalcResult();
}

void    GAnalysis3MesonsProton::Fill(const GTreeMeson& meson, GTreeParticle& photons, GTreeParticle& proton, const GTreeTagger& tagger, const GTreeA2Geant& geantTree)
{
    success = false;

    Double_t    im  = meson.Particle(0).M();
    Double_t    theta  = meson.Particle(0).Theta()*TMath::RadToDeg();
    Double_t    phi  = meson.Particle(0).Phi()*TMath::RadToDeg();
    Double_t    protonTheta  = proton.Particle(0).Theta()*TMath::RadToDeg();
    Double_t    protonPhi  = proton.Particle(0).Phi()*TMath::RadToDeg();
    Double_t    mm;
    Double_t    sub_im_0    = (photons.Particle(0) + photons.Particle(1)).M();
    Double_t    sub_im_1    = (photons.Particle(2) + photons.Particle(3)).M();
    Double_t    sub_im_2    = (photons.Particle(4) + photons.Particle(5)).M();

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
        TLorentzVector  helpCM(meson.Particle(0));
        helpCM.Boost(-tagger.GetVectorProtonTarget(i).BoostVector());
        TLorentzVector  helpProtonCM(meson.Particle(0));
        helpProtonCM.Boost(-tagger.GetVectorProtonTarget(i).BoostVector());

        fitProton6.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
        fitProton6.SetBeam(tagger.GetTaggedEnergy(i));
        fitProton6.SetProton(proton.Particle(0));

        if(fitProton6.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i), geantTree))
        {
            success = true;
            hist_fitProton6.Fill(im, mm, theta, phi, helpCM.Theta()*TMath::RadToDeg(), proton.Particle(0).E(), protonTheta, protonPhi, helpProtonCM.Theta()*TMath::RadToDeg(), sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            for(int k=0; k<6; k++)
            {
                for(int l=k+1; l<6; l++)
                    hist_fitProton6_SubAll.Fill((photons.Particle(k)+photons.Particle(l)).M(), tagger.GetTaggedTime(i));

                photons.AddParticle(fitProton6.GetFittedPhoton(k).Ek, fitProton6.GetFittedPhoton(k).Theta * TMath::RadToDeg(), fitProton6.GetFittedPhoton(k).Phi * TMath::RadToDeg(),
                                    0, photons.GetTime(k), photons.GetClusterSize(k), photons.GetCentralCrystal(k), photons.GetCentralVeto(k), photons.GetDetectors(k),
                                    photons.GetVetoEnergy(k), photons.GetMWPC0Energy(k), photons.GetMWPC1Energy(k), photons.GetTrackIndex(k));
            }
            proton.AddParticle(fitProton6.GetFittedProton().Ek, fitProton6.GetFittedProton().Theta * TMath::RadToDeg(), fitProton6.GetFittedProton().Phi * TMath::RadToDeg(),
                               MASS_PROTON, proton.GetTime(0), proton.GetClusterSize(0), proton.GetCentralCrystal(0), proton.GetCentralVeto(0), proton.GetDetectors(0),
                               proton.GetVetoEnergy(0), proton.GetMWPC0Energy(0), proton.GetMWPC1Energy(0), proton.GetTrackIndex(0));
        }
    }
}

void    GAnalysis3MesonsProton::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* folder  = arr->GetDirectory("WithProton");

    fitProton6.PrepareWriteList(folder);
    hist_fitProton6.PrepareWriteList(folder, "fitProton6");
    hist_fitProton6_SubAll.PrepareWriteList(folder, "SubAll");
}

void    GAnalysis3MesonsProton::Reset(Option_t* option)
{
    fitProton6.Reset(option);
    hist_fitProton6.Reset(option);
    hist_fitProton6_SubAll.Reset(option);
}
