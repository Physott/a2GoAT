#include "GHistPhysics.h"
#include "GTreeParticle.h"
#include "GTreeTagger.h"
#include "GTreeMeson.h"
#include "GTreeA2Geant.h"



GHistParticle::GHistParticle(const char* Name, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    name(Name),
    kinEnergy(TString(Name).Append("_kinEnergy"), "kinEnergy", 1600, 0, 1600, 48, kFALSE),
    theta(TString(Name).Append("_theta"), "theta", 180, 0, 180, 48, kFALSE),
    thetaCM(TString(Name).Append("_thetaCM"), "thetaCM", 180, 0, 180, 48, kFALSE),
    phi(TString(Name).Append("_phi"), "phi", 360, -180, 180, 48, kFALSE),
    mass(TString(Name).Append("_mass"), "mass", 1600, 0, 1600, 48, kFALSE),
    energy(TString(Name).Append("_energy"), "energy", 1000, 0, 2500, 48, kFALSE),
    px(TString(Name).Append("_px"), "px", 1000, -1000, 1000, 48, kFALSE),
    py(TString(Name).Append("_py"), "py", 1000, -1000, 1000, 48, kFALSE),
    pz(TString(Name).Append("_pz"), "pz", 1000, 0, 2000, 48, kFALSE),
    time(TString(Name).Append("_time"), "time", 1000, -100, 100, 48, kFALSE),
    clusterSize(TString(Name).Append("_clusterSize"), "clusterSize", 20, 0, 20, 48, kFALSE),
    centralCrystal(TString(Name).Append("_centralCrystal"), "centralCrystal", 800, 0, 800, 48, kFALSE),
    detectors(TString(Name).Append("_detectors"), "detectors", 100, -100, 100, 48, kFALSE),
    trackIndex(TString(Name).Append("_trackIndex"), "trackIndex", 10, 0, 10, 48, kFALSE)
{
}

void    GHistParticle::CalcResult()
{
    kinEnergy.CalcResult();
    theta.CalcResult();
    thetaCM.CalcResult();
    phi.CalcResult();
    mass.CalcResult();
    energy.CalcResult();
    px.CalcResult();
    py.CalcResult();
    pz.CalcResult();
    time.CalcResult();
    clusterSize.CalcResult();
    centralCrystal.CalcResult();
    detectors.CalcResult();
    trackIndex.CalcResult();
}

void    GHistParticle::Fill(const GTreeParticle& particle, const int index, const double beam, const double Time)
{
    try
    {
        kinEnergy.Fill(particle.GetClusterEnergy(index), Time);
        theta.Fill(particle.GetTheta(index), Time);
        TLorentzVector  lv(particle.Particle(index));
        if(beam>0)
        {
            TLorentzVector  helpCM(lv);
            TLorentzVector  helpCM3(0.0, 0.0, beam, beam + MASS_PROTON);
            helpCM.Boost(-helpCM3.BoostVector());
            thetaCM.Fill(helpCM.Theta()*TMath::RadToDeg(), Time);
        }
        phi.Fill(particle.GetPhi(index), Time);
        mass.Fill(particle.GetMass(index), Time);
        energy.Fill(lv.E(), Time);
        px.Fill(lv.Px(), Time);
        py.Fill(lv.Py(), Time);
        pz.Fill(lv.Pz(), Time);
        time.Fill(particle.GetTime(index), Time);
        clusterSize.Fill(particle.GetClusterSize(index), Time);
        centralCrystal.Fill(particle.GetCentralCrystal(index), Time);
        detectors.Fill(particle.GetDetectors(index), Time);
        trackIndex.Fill(particle.GetTrackIndex(index), Time);
    }
    catch(...)
    {
        std::cout << "Can not find particle " << index << " in tree." << std::endl;
        return;
    }
}

void    GHistParticle::Fill(const GTreeParticle& particle, const int index, const double beam, const double Time, const double channel)
{
    try
    {
        kinEnergy.Fill(particle.GetClusterEnergy(index), Time, channel);
        theta.Fill(particle.GetTheta(index), Time, channel);
        TLorentzVector  lv(particle.Particle(index));
        if(beam>0)
        {
            TLorentzVector  helpCM(lv);
            TLorentzVector  helpCM3(0.0, 0.0, beam, beam + MASS_PROTON);
            helpCM.Boost(-helpCM3.BoostVector());
            thetaCM.Fill(helpCM.Theta()*TMath::RadToDeg(), Time, channel);
        }
        phi.Fill(particle.GetPhi(index), Time, channel);
        mass.Fill(particle.GetMass(index), Time, channel);
        energy.Fill(lv.E(), Time, channel);
        px.Fill(lv.Px(), Time, channel);
        py.Fill(lv.Py(), Time, channel);
        pz.Fill(lv.Pz(), Time, channel);
        time.Fill(particle.GetTime(index), Time, channel);
        clusterSize.Fill(particle.GetClusterSize(index), Time, channel);
        centralCrystal.Fill(particle.GetCentralCrystal(index), Time, channel);
        detectors.Fill(particle.GetDetectors(index), Time, channel);
        trackIndex.Fill(particle.GetTrackIndex(index), Time, channel);
    }
    catch(...)
    {
        std::cout << "Can not find particle " << index << " in tree." << std::endl;
        return;
    }
}

void    GHistParticle::Fill(const TLorentzVector& particle, const double beam, const double Time, const double channel)
{
    double  m   = particle.M();
    kinEnergy.Fill(particle.E()-m, Time, channel);
    theta.Fill(particle.Theta()*TMath::RadToDeg(), Time, channel);
    if(beam>0)
    {
        TLorentzVector  helpCM(particle);
        TLorentzVector  helpCM3(0.0, 0.0, beam, beam + MASS_PROTON);
        helpCM.Boost(-helpCM3.BoostVector());
        thetaCM.Fill(helpCM.Theta()*TMath::RadToDeg(), Time, channel);
    }
    phi.Fill(particle.Phi()*TMath::RadToDeg(), Time, channel);
    mass.Fill(m, Time, channel);
    energy.Fill(particle.E(), Time, channel);
    px.Fill(particle.Px(), Time, channel);
    py.Fill(particle.Py(), Time, channel);
    pz.Fill(particle.Pz(), Time, channel);
}

void    GHistParticle::PrepareWriteList(GHistWriteList* arr, const char* Name)
{
    if(!arr)
        return;

    GHistWriteList* folder  = arr;
    if(Name)
        folder  = arr->GetDirectory(Name);

    kinEnergy.PrepareWriteList(folder, "kinEnergy");
    theta.PrepareWriteList(folder, "theta");
    thetaCM.PrepareWriteList(folder, "thetaCM");
    phi.PrepareWriteList(folder, "phi");
    mass.PrepareWriteList(folder, "mass");
    energy.PrepareWriteList(folder, "energy");
    px.PrepareWriteList(folder, "px");
    py.PrepareWriteList(folder, "py");
    pz.PrepareWriteList(folder, "pz");
    time.PrepareWriteList(folder, "time");
    clusterSize.PrepareWriteList(folder, "clusterSize");
    centralCrystal.PrepareWriteList(folder, "centralCrystal");
    detectors.PrepareWriteList(folder, "detectors");
    trackIndex.PrepareWriteList(folder, "trackIndex");
}

void    GHistParticle::Reset(Option_t* option)
{
    kinEnergy.Reset(option);
    theta.Reset(option);
    thetaCM.Reset(option);
    phi.Reset(option);
    mass.Reset(option);
    energy.Reset(option);
    px.Reset(option);
    py.Reset(option);
    pz.Reset(option);
    time.Reset(option);
    clusterSize.Reset(option);
    centralCrystal.Reset(option);
    detectors.Reset(option);
    trackIndex.Reset(option);
}

void    GHistParticle::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    kinEnergy.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    theta.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    thetaCM.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    phi.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    mass.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    energy.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    px.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    py.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pz.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    time.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    clusterSize.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    centralCrystal.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    detectors.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    trackIndex.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}



















GHistPhysics::GHistPhysics(const char* Name, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    name(Name),
    proton(TString(Name).Append("_proton"), kFALSE),
    etap(TString(Name).Append("_etap"), kFALSE),
    etaPhotons(TString(Name).Append("_etaPhotons"), kFALSE),
    pi0Photons(TString(Name).Append("_pi0Photons"), kFALSE),
    allPhotons(TString(Name).Append("_allPhotons"), kFALSE)
{
}

void    GHistPhysics::CalcResult()
{
    proton.CalcResult();
    etap.CalcResult();
    etaPhotons.CalcResult();
    pi0Photons.CalcResult();
    allPhotons.CalcResult();
}

void    GHistPhysics::Fill(const GTreeMeson& meson, const GTreeParticle& photons, const GTreeParticle& protons, const GTreeTagger& tagger)
{
    try
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            etap.Fill(meson, 0, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            if(protons.GetNParticles()>0)
                proton.Fill(protons, 0, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            etaPhotons.Fill(photons, 0, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            etaPhotons.Fill(photons, 1, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(photons, 2, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(photons, 3, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(photons, 4, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(photons, 5, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 0, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 1, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 2, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 3, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 4, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 5, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
        }
    }
    catch(...)
    {
        std::cout << "Can not find all particles in GHistPhysics!" << std::endl;
        return;
    }
}

void    GHistPhysics::FillFitted(const GTreeParticle& photons, const GTreeParticle &protons, const GTreeTagger &tagger)
{
    try
    {
        TLorentzVector  lv(photons.Particle(6));
        for(int i=7; i<12; i++)
            lv  +=  photons.Particle(i);

        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            etap.Fill(lv, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            if(protons.GetNParticles()>0)
                proton.Fill(protons, 1, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            etaPhotons.Fill(photons, 6, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            etaPhotons.Fill(photons, 7, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(photons, 8, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(photons, 9, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(photons, 10, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(photons, 11, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 6, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 7, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 8, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 9, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 10, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(photons, 11, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
        }
    }
    catch(...)
    {
        //std::cout << "no true tree!" << std::endl;
        return;
    }
}

void    GHistPhysics::FillTrue(const GTreeA2Geant& geant, const GTreeTagger &tagger)
{
    try
    {
        TLorentzVector  ph[6];
        TLorentzVector  lv(0.0, 0.0, 0.0, 0.0);
        for(int i=0; i<6; i++)
        {
            ph[i]   = geant.GetTrueVector(i+6)*1000.0;
            lv     += ph[i];
        }
        TLorentzVector  pr(geant.GetTrueVector(1)*1000.0);

        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            etap.Fill(lv, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            proton.Fill(pr, tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            etaPhotons.Fill(ph[0], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            etaPhotons.Fill(ph[1], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(ph[2], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(ph[3], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(ph[4], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            pi0Photons.Fill(ph[5], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(ph[0], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(ph[1], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(ph[2], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(ph[3], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(ph[4], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            allPhotons.Fill(ph[5], tagger.GetTaggedEnergy(i), tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
        }
    }
    catch(...)
    {
        //std::cout << "no true tree!" << std::endl;
        return;
    }
}

void    GHistPhysics::PrepareWriteList(GHistWriteList* arr, const char* Name)
{
    if(!arr)
        return;


    GHistWriteList* folder  = arr;
    if(Name)
        folder  = arr->GetDirectory(Name);

    proton.PrepareWriteList(folder, "proton");
    etap.PrepareWriteList(folder, "etap");
    etaPhotons.PrepareWriteList(folder, "etaPhotons");
    pi0Photons.PrepareWriteList(folder, "pi0Photons");
    allPhotons.PrepareWriteList(folder, "allPhotons");
}

void    GHistPhysics::Reset(Option_t* option)
{
    proton.Reset(option);
    etap.Reset(option);
    etaPhotons.Reset(option);
    pi0Photons.Reset(option);
    allPhotons.Reset(option);
}

void    GHistPhysics::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    proton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    etap.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    etaPhotons.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pi0Photons.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    allPhotons.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}



















GHistPhysicsFitted::GHistPhysicsFitted(const char* Name, Bool_t linkHistogram)  :
    GHistPhysics(TString(Name).Append("_unfitted"), linkHistogram),
    name(Name),
    fitted(TString(Name).Append("_fitted"), kFALSE),
    trueValues(TString(Name).Append("_true"), kFALSE)
{

}

void    GHistPhysicsFitted::CalcResult()
{
    GHistPhysics::CalcResult();
    fitted.CalcResult();
    trueValues.CalcResult();
}

void    GHistPhysicsFitted::Fill(const GTreeMeson& meson, const GTreeParticle& photons, const GTreeParticle& protons, const GTreeA2Geant& geant, const GTreeTagger& tagger)
{
    GHistPhysics::Fill(meson, photons, protons, tagger);
    fitted.FillFitted(photons, protons, tagger);
    trueValues.FillTrue(geant, tagger);
}

void    GHistPhysicsFitted::PrepareWriteList(GHistWriteList* arr, const char* Name)
{
    if(!arr)
        return;

    GHistWriteList* folder  = arr->GetDirectory(name);

    GHistPhysics::PrepareWriteList(folder, "raw");
    fitted.PrepareWriteList(folder, "fit");
    trueValues.PrepareWriteList(folder, "true");
}

void    GHistPhysicsFitted::Reset(Option_t* option)
{
    GHistPhysics::Reset(option);
    fitted.Reset(option);
    trueValues.Reset(option);
}

void    GHistPhysicsFitted::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    GHistPhysics::ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    fitted.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    trueValues.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}
