#include "GHistPhysics.h"
#include "GTreeParticle.h"



GHistParticle::GHistParticle(const char* Name, const char* title, Bool_t linkHistogram) :
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

void    GHistParticle::PrepareWriteList(GHistWriteList* arr, const char* Name)
{
    if(!arr)
        return;

    GHistWriteList* folder  = arr->GetDirectory(name);

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
