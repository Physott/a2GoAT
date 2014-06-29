#ifndef __PProtonCheck_h__
#define __PProtonCheck_h__

#include "GTreeManager.h"
#include "PHist.h"

class	PProtonCheck : virtual public GTreeManager
{
private:
    UInt_t      foundProtons;
    UInt_t      foundProtonsIndex[GTreeParticle_MaxEntries];
    UInt_t      smallestProtonAngleDiffIndex[GTreeParticle_MaxEntries];

    TH1I        histFoundProtons;
    PHistD      histProtonAngleDiff;
    PHistD      histSmallestProtonAngleDiff;
    TH1D        histCoplanarity;

    Double_t    cutCoplanarity[2];

protected:


    inline  void            Clear();
            TLorentzVector& GetProton(const int i)          {return protons->Particle(foundProtonsIndex[i]);}
    const   TLorentzVector& GetProton(const int i)  const   {return protons->Particle(foundProtonsIndex[i]);}
            Bool_t          ProcessEvent(const TLorentzVector& all);
    inline  void            RandomSubtraction();
            void            Write();

public:
    PProtonCheck();
    virtual ~PProtonCheck();

};


void    PProtonCheck::Clear()
{
    histFoundProtons.Reset();
    histProtonAngleDiff.Reset();
    histSmallestProtonAngleDiff.Reset();
    histCoplanarity.Reset();
}

void    PProtonCheck::RandomSubtraction()
{
    histProtonAngleDiff.RandomSubtraction();
    histSmallestProtonAngleDiff.RandomSubtraction();
}

#endif
