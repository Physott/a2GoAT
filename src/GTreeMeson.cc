#include "GTreeMeson.h"
#include "GTreeManager.h"


GTreeMeson::GTreeMeson(GTreeManager *Manager, const TString& _Name)    :
    GTreeParticle(Manager, _Name)
{
    for(Int_t i=0; i<GTreeTrack_MAX; i++)
    {
        nSubParticles[i]  = 0;
        nSubRootinos[i]   = 0;
        nSubPhotons[i]    = 0;
        nSubChargedPions[i]  = 0;
        trackIndex[i]     = 0;
    }
}

GTreeMeson::~GTreeMeson()
{
}

void    GTreeMeson::Clear()
{
    GTreeParticle::Clear();
}

void    GTreeMeson::SetBranchAdresses()
{
    GTreeParticle::SetBranchAdresses();
    inputTree->SetBranchAddress("nSubParticles", nSubParticles);
    inputTree->SetBranchAddress("nSubRootinos", nSubRootinos);
    inputTree->SetBranchAddress("nSubPhotons", nSubPhotons);
    inputTree->SetBranchAddress("nSubChargedPions", nSubChargedPions);
}

void    GTreeMeson::SetBranches()
{
    GTreeParticle::SetBranches();
    outputTree->Branch("nSubParticles", nSubParticles, "nSubParticles[nParticles]/I");
    outputTree->Branch("nSubRootinos", nSubRootinos,"nSubRootinos[nParticles]/I");
    outputTree->Branch("nSubPhotons", nSubPhotons,"nSubParticles[nParticles]/I");
    outputTree->Branch("nSubChargedPions", nSubChargedPions,"nSubParticles[nParticles]/I");
}

void    GTreeMeson::AddParticle(const Int_t _NSubRootinos, const Int_t _NSubPhotons, const Int_t _NSubChargedPions, Int_t* subParticles_index, TLorentzVector* subParticles_list)
{
    nSubParticles[nParticles]  = _NSubRootinos + _NSubPhotons + _NSubChargedPions;
    nSubRootinos[nParticles]   = _NSubRootinos;
    nSubPhotons[nParticles]    = _NSubPhotons;
    nSubChargedPions[nParticles]  = _NSubChargedPions;

    TLorentzVector sum;

    for(Int_t i=0; i<_NSubRootinos; i++)
    {
        sum += subParticles_list[i];
        detectors[nParticles]    = detectors[nParticles] | manager->rootinos->GetDetectors(subParticles_index[i]);
        time[nParticles]         += manager->rootinos->GetTime(subParticles_index[i]);
        clusterSize[nParticles]  += manager->rootinos->GetClusterSize(subParticles_index[i]);
        vetoEnergy[nParticles]   += manager->rootinos->GetVetoEnergy(subParticles_index[i]);
        MWPC0Energy[nParticles]  += manager->rootinos->GetMWPC0Energy(subParticles_index[i]);
        MWPC1Energy[nParticles]  += manager->rootinos->GetMWPC1Energy(subParticles_index[i]);
        trackIndex[nParticles]   += 1<<(manager->rootinos->GetTrackIndex(subParticles_index[i]));
    }
    for(Int_t i=_NSubRootinos; i<_NSubRootinos+_NSubPhotons; i++)
    {
        sum += subParticles_list[i];
        detectors[nParticles]    = detectors[nParticles] | manager->photons->GetDetectors(subParticles_index[i]);
        time[nParticles]         += manager->photons->GetTime(subParticles_index[i]);
        clusterSize[nParticles]  += manager->photons->GetClusterSize(subParticles_index[i]);
        vetoEnergy[nParticles]   += manager->photons->GetVetoEnergy(subParticles_index[i]);
        MWPC0Energy[nParticles]  += manager->photons->GetMWPC0Energy(subParticles_index[i]);
        MWPC1Energy[nParticles]  += manager->photons->GetMWPC1Energy(subParticles_index[i]);
        trackIndex[nParticles]   += 1<<(manager->photons->GetTrackIndex(subParticles_index[i]));
    }
    for(Int_t i=_NSubRootinos+_NSubPhotons; i<_NSubRootinos+_NSubPhotons+_NSubChargedPions; i++)
    {
        sum += subParticles_list[i];
        detectors[nParticles]    = detectors[nParticles] | manager->chargedPions->GetDetectors(subParticles_index[i]);
        time[nParticles]         += manager->chargedPions->GetTime(subParticles_index[i]);
        clusterSize[nParticles]  += manager->chargedPions->GetClusterSize(subParticles_index[i]);
        vetoEnergy[nParticles]   += manager->chargedPions->GetVetoEnergy(subParticles_index[i]);
        MWPC0Energy[nParticles]  += manager->chargedPions->GetMWPC0Energy(subParticles_index[i]);
        MWPC1Energy[nParticles]  += manager->chargedPions->GetMWPC1Energy(subParticles_index[i]);
        trackIndex[nParticles]   += 1<<(manager->chargedPions->GetTrackIndex(subParticles_index[i]));
    }

    clusterEnergy[nParticles] = sum.E()-sum.M();
    theta[nParticles] = sum.Theta()*TMath::RadToDeg();
    phi[nParticles] = sum.Phi()*TMath::RadToDeg();
    mass[nParticles] = sum.M();
    time[nParticles]         /= nSubParticles[nParticles];
    nParticles++;
    manager->countReconstructed++;
}


void    GTreeMeson::AddParticle(const Int_t subParticle_index0, const TLorentzVector subParticle0, const Int_t pdg0, const Int_t subParticle_index1, const TLorentzVector subParticle1, const Int_t pdg1)
{
    nSubParticles[nParticles]  = 2;
    nSubRootinos[nParticles]   = 0;
    nSubPhotons[nParticles]    = 0;
    nSubChargedPions[nParticles]  = 0;

    TLorentzVector sum = subParticle0 + subParticle1;

    if(pdg0 == manager->pdgDB->GetParticle("gamma")->PdgCode())
    {
        detectors[nParticles]    = manager->photons->GetDetectors(subParticle_index0);
        time[nParticles]         = manager->photons->GetTime(subParticle_index0);
        clusterSize[nParticles]  = manager->photons->GetClusterSize(subParticle_index0);
        vetoEnergy[nParticles]   = manager->photons->GetVetoEnergy(subParticle_index0);
        MWPC0Energy[nParticles]  = manager->photons->GetMWPC0Energy(subParticle_index0);
        MWPC1Energy[nParticles]  = manager->photons->GetMWPC1Energy(subParticle_index0);
        trackIndex[nParticles]   = 1<<(manager->photons->GetTrackIndex(subParticle_index0));
        nSubPhotons[nParticles]++;
    }
    else if(pdg0 == manager->pdgDB->GetParticle("pi+")->PdgCode())
    {
        detectors[nParticles]    = manager->chargedPions->GetDetectors(subParticle_index0);
        time[nParticles]         = manager->chargedPions->GetTime(subParticle_index0);
        clusterSize[nParticles]  = manager->chargedPions->GetClusterSize(subParticle_index0);
        vetoEnergy[nParticles]   = manager->chargedPions->GetVetoEnergy(subParticle_index0);
        MWPC0Energy[nParticles]  = manager->chargedPions->GetMWPC0Energy(subParticle_index0);
        MWPC1Energy[nParticles]  = manager->chargedPions->GetMWPC1Energy(subParticle_index0);
        trackIndex[nParticles]   = 1<<(manager->chargedPions->GetTrackIndex(subParticle_index0));
        nSubChargedPions[nParticles]++;
    }
    else
    {
        detectors[nParticles]    = manager->rootinos->GetDetectors(subParticle_index0);
        time[nParticles]         = manager->rootinos->GetTime(subParticle_index0);
        clusterSize[nParticles]  = manager->rootinos->GetClusterSize(subParticle_index0);
        vetoEnergy[nParticles]   = manager->rootinos->GetVetoEnergy(subParticle_index0);
        MWPC0Energy[nParticles]  = manager->rootinos->GetMWPC0Energy(subParticle_index0);
        MWPC1Energy[nParticles]  = manager->rootinos->GetMWPC1Energy(subParticle_index0);
        trackIndex[nParticles]   = 1<<(manager->rootinos->GetTrackIndex(subParticle_index0));
        nSubRootinos[nParticles]++;
    }
    
    if(pdg1 == manager->pdgDB->GetParticle("gamma")->PdgCode())
    {
        detectors[nParticles]    = detectors[nParticles] | manager->photons->GetDetectors(subParticle_index1);
        time[nParticles]         += manager->photons->GetTime(subParticle_index1);
        clusterSize[nParticles]  += manager->photons->GetClusterSize(subParticle_index1);
        vetoEnergy[nParticles]   += manager->photons->GetVetoEnergy(subParticle_index1);
        MWPC0Energy[nParticles]  += manager->photons->GetMWPC0Energy(subParticle_index1);
        MWPC1Energy[nParticles]  += manager->photons->GetMWPC1Energy(subParticle_index1);
        trackIndex[nParticles]   += 1<<(manager->photons->GetTrackIndex(subParticle_index1));
        nSubPhotons[nParticles]++;
    }
    else if(pdg1 == manager->pdgDB->GetParticle("pi+")->PdgCode())
    {
        detectors[nParticles]    = detectors[nParticles] | manager->chargedPions->GetDetectors(subParticle_index1);
        time[nParticles]         += manager->chargedPions->GetTime(subParticle_index1);
        clusterSize[nParticles]  += manager->chargedPions->GetClusterSize(subParticle_index1);
        vetoEnergy[nParticles]   += manager->chargedPions->GetVetoEnergy(subParticle_index1);
        MWPC0Energy[nParticles]  += manager->chargedPions->GetMWPC0Energy(subParticle_index1);
        MWPC1Energy[nParticles]  += manager->chargedPions->GetMWPC1Energy(subParticle_index1);
        trackIndex[nParticles]   += 1<<(manager->chargedPions->GetTrackIndex(subParticle_index1));
        nSubChargedPions[nParticles]++;
    }
    else
    {
        detectors[nParticles]    = detectors[nParticles] | manager->rootinos->GetDetectors(subParticle_index1);
        time[nParticles]         += manager->rootinos->GetTime(subParticle_index1);
        clusterSize[nParticles]  += manager->rootinos->GetClusterSize(subParticle_index1);
        vetoEnergy[nParticles]   += manager->rootinos->GetVetoEnergy(subParticle_index1);
        MWPC0Energy[nParticles]  += manager->rootinos->GetMWPC0Energy(subParticle_index1);
        MWPC1Energy[nParticles]  += manager->rootinos->GetMWPC1Energy(subParticle_index1);
        trackIndex[nParticles]   += 1<<(manager->rootinos->GetTrackIndex(subParticle_index1));
        nSubRootinos[nParticles]++;
    }

    clusterEnergy[nParticles] = sum.E()-sum.M();
    theta[nParticles] = sum.Theta()*TMath::RadToDeg();
    phi[nParticles] = sum.Phi()*TMath::RadToDeg();
    mass[nParticles] = sum.M();
    time[nParticles]         /= 2;
    nParticles++;
    manager->countReconstructed++;
}


TLorentzVector GTreeMeson::SubParticles(const Int_t meson, const Int_t particleIndex)
{
    if(meson>=nParticles || meson<0)    return TLorentzVector();
    if(particleIndex>=nSubParticles[meson] || particleIndex<0)    return TLorentzVector();
    manager->GetTracks()->GetVector(trackIndex[particleIndex]).Print();
    return manager->GetTracks()->GetVector(trackIndex[particleIndex]);
}

const   TLorentzVector GTreeMeson::SubParticles(const Int_t meson, const Int_t particleIndex) const
{
    if(meson>=nParticles || meson<0)    return TLorentzVector();
    if(particleIndex>=nSubParticles[meson] || particleIndex<0)    return TLorentzVector();
    return manager->GetTracks()->GetVector(trackIndex[particleIndex]);
}

TLorentzVector GTreeMeson::SubRootinos(const Int_t meson, const Int_t rootinoIndex)
{
    if(meson>=nParticles || meson<0)    return TLorentzVector();
    if(rootinoIndex>=nSubRootinos[meson] || rootinoIndex<0)    return TLorentzVector();
    return manager->GetTracks()->GetVector(trackIndex[rootinoIndex]);
}

const   TLorentzVector GTreeMeson::SubRootinos(const Int_t meson, const Int_t rootinoIndex) const
{
    if(meson>=nParticles || meson<0)    return TLorentzVector();
    if(rootinoIndex>=nSubRootinos[meson] || rootinoIndex<0)    return TLorentzVector();
    return manager->GetTracks()->GetVector(trackIndex[rootinoIndex]);
}

TLorentzVector GTreeMeson::SubPhotons(const Int_t meson, const Int_t photonIndex)
{
    if(meson>=nParticles || meson<0)    return TLorentzVector();
    if(photonIndex>=nSubPhotons[meson] || photonIndex<0)    return TLorentzVector();
    return manager->GetTracks()->GetVector(trackIndex[nSubRootinos[meson]+photonIndex]);
}

const   TLorentzVector GTreeMeson::SubPhotons(const Int_t meson, const Int_t photonIndex) const
{
    if(meson>=nParticles || meson<0)    return TLorentzVector();
    if(photonIndex>=nSubPhotons[meson] || photonIndex<0)    return TLorentzVector();
    return manager->GetTracks()->GetVector(trackIndex[nSubRootinos[meson]+photonIndex]);
}

TLorentzVector GTreeMeson::SubChargedPi(const Int_t meson, const Int_t chargedPiIndex)
{
    if(meson>=nParticles || meson<0)    return TLorentzVector();
    if(chargedPiIndex>=nSubChargedPions[meson] || chargedPiIndex<0)    return TLorentzVector();
    return manager->GetTracks()->GetVector(trackIndex[nSubRootinos[meson]+nSubPhotons[meson]+chargedPiIndex]);
}

const   TLorentzVector GTreeMeson::SubChargedPi(const Int_t meson, const Int_t chargedPiIndex) const
{
    if(meson>=nParticles || meson<0)    return TLorentzVector();
    if(chargedPiIndex>=nSubChargedPions[meson] || chargedPiIndex<0)    return TLorentzVector();
    return manager->GetTracks()->GetVector(trackIndex[nSubRootinos[meson]+nSubPhotons[meson]+chargedPiIndex]);
}


void    GTreeMeson::Print() const
{
    GTree::Print();
    cout << "nParticles: " << nParticles << endl;
    for(Int_t i=0; i<nParticles; i++)
    {
        GTreeParticle::PrintParticle(i);
        cout << "nSubParticles: " << nSubParticles[i] << "   nSubRootinos: " << nSubRootinos[i] << "   nSubPhotons: " << nSubPhotons[i] << "   nSubChargedPions: " << nSubChargedPions[i] << endl;
    }
}

Int_t               nSubParticles[GTreeTrack_MAX];
Int_t               nSubRootinos[GTreeTrack_MAX];
Int_t               nSubPhotons[GTreeTrack_MAX];
Int_t               nSubChargedPions[GTreeTrack_MAX];
