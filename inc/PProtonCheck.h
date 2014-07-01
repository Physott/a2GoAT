#ifndef __PProtonCheck_h__
#define __PProtonCheck_h__

#include "GTreeManager.h"
#include "PHist.h"

class	PProtonCheck
{
private:
    TString name;

    UInt_t      foundProtons;
    UInt_t      foundProtonsIndex[GTreeParticle_MaxEntries];
    UInt_t      smallestProtonAngleDiffIndex[GTreeParticle_MaxEntries];

    TH1I        histFoundProtons;
    PHistD      histProtonAngleDiff;
    PHistD      histSmallestProtonAngleDiff;
    TH1D        histCoplanarity;

    Double_t    cutCoplanarity[2];

protected:

//            TLorentzVector& GetProton(const int i)          {return protons->Particle(foundProtonsIndex[i]);}
//    const   TLorentzVector& GetProton(const int i)  const   {return protons->Particle(foundProtonsIndex[i]);}

public:
    PProtonCheck(const TString& _Name);
    virtual ~PProtonCheck();

    inline  void            Clear();
            Bool_t          ProcessEvent(const TLorentzVector& all, const GTreeParticle& protons, const GTreeTagger& tagger);
    inline  void            RandomSubtraction();
            void            Write(TDirectory& curDir);
};


void    PProtonCheck::Clear()
{
    histFoundProtons.Clear();
    histProtonAngleDiff.Clear();
    histSmallestProtonAngleDiff.Clear();
    histCoplanarity.Clear();
}

void    PProtonCheck::RandomSubtraction()
{
    histProtonAngleDiff.RandomSubtraction();
    histSmallestProtonAngleDiff.RandomSubtraction();
}

#endif
