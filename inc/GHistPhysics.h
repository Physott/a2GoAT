#ifndef __GHistPhysics_h__
#define __GHistPhysics_h__


#include "GHistManager.h"
#include "GH1.h"

#define MASS_PROTON 938.272046



class   GTreeParticle;

class	GHistParticle  : public GHistLinked
{
private:
    TString name;
    GH1     kinEnergy;
    GH1     theta;
    GH1     thetaCM;
    GH1     phi;
    GH1     mass;
    GH1     energy;
    GH1     px;
    GH1     py;
    GH1     pz;
    GH1     time;
    GH1     clusterSize;
    GH1     centralCrystal;
    GH1     detectors;
    GH1     trackIndex;

protected:

public:
    GHistParticle(const char* Name, const char* title, Bool_t linkHistogram = kTRUE);
    virtual ~GHistParticle()                                                                                                    {}

    virtual void    CalcResult();
    virtual Int_t   Fill(Double_t x)                                                                                            {return 0;}
    virtual void    Fill(const GTreeParticle& particle, const int index, const double beam, const double Time);
    virtual void    Fill(const GTreeParticle& particle, const int index, const double beam, const double Time, const double channel);
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* Name = 0);
    virtual void    Reset(Option_t* option);
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads);
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)                           {return 0;}
};



#endif
