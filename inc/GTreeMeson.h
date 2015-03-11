#ifndef __GTreeMeson_h__
#define __GTreeMeson_h__


#include "Rtypes.h"
#include "GTreeParticle.h"


#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046


using namespace std;

class  GTreeMeson   : public GTreeParticle
{
private:
    Int_t             nSubParticles[GTreeTrack_MAX];
    Int_t             nSubRootinos[GTreeTrack_MAX];
    Int_t             nSubPhotons[GTreeTrack_MAX];
    Int_t             nSubChargedPions[GTreeTrack_MAX];

protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeMeson(GTreeManager *Manager, const TString& _Name);
    virtual ~GTreeMeson();

            void            AddParticle(Int_t subParticle_index0, const TLorentzVector subParticle0, const Int_t pdg0, Int_t subParticle_index1, const TLorentzVector subParticle1, const Int_t pdg1);
            void            AddParticle(const Int_t _NSubRootinos, const Int_t _NSubPhotons, const Int_t _NSubChargedPions, Int_t* subParticles_index, TLorentzVector* subParticles_list);
    virtual void            Clear();
            Int_t          GetNMesons()                            const	{return GetNParticles();}
            Int_t          GetNSubParticles(const Int_t meson)     const	{return nSubParticles[meson];}
            Int_t          GetNSubRootinos(const Int_t meson)      const	{return nSubRootinos[meson];}
            Int_t          GetNSubPhotons(const Int_t meson)       const	{return nSubPhotons[meson];}
            Int_t          GetNSubChargedPi(const Int_t meson)     const	{return nSubChargedPions[meson];}
    inline  const   std::vector<Int_t> GetTrackIndexList(const Int_t meson)    const;
            Bool_t         HasTrack(const Int_t meson, const Int_t track) const {return (trackIndex[meson] & 1<<track);}
            TLorentzVector Meson(const Int_t meson)             {return Particle(meson);}
    const   TLorentzVector Meson(const Int_t meson) const       {return Particle(meson);}
            TLorentzVector SubParticles(const Int_t meson, const Int_t particleIndex);
    const   TLorentzVector SubParticles(const Int_t meson, const Int_t particleIndex) const;
            TLorentzVector SubPhotons(const Int_t meson, const Int_t photonIndex);
    const   TLorentzVector SubPhotons(const Int_t meson, const Int_t photonIndex) const;
            TLorentzVector SubRootinos(const Int_t meson, const Int_t rootinoIndex);
    const   TLorentzVector SubRootinos(const Int_t meson, const Int_t rootinoIndex) const;
            TLorentzVector SubChargedPi(const Int_t meson, const Int_t chargedPiIndex);
    const   TLorentzVector SubChargedPi(const Int_t meson, const Int_t chargedPiIndex) const;
    virtual void            Print() const;
};

const std::vector<Int_t> GTreeMeson::GetTrackIndexList(const Int_t meson) const
{
    std::vector<Int_t> trackIndexList;
    Int_t i=0;

    for(Int_t j=0; j<nSubParticles[meson]; j++)
    {
        while((trackIndex[meson] & 1<<i) == 0) i++;
        trackIndexList.push_back(i++);
    }

    return trackIndexList;
}

#endif
