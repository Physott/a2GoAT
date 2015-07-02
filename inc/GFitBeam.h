#ifndef __GFitBeam_h__
#define __GFitBeam_h__


#include <TLorentzVector.h>

#include "GFit.h"



class	GFitBeam    : public GFit
{
private:
    FitParticle         aplconBeam;

    GH1                 beamTheta;
    GHistBGSub          beamPhi;

    virtual int GetNDOF()   {return 21;}

public:
    GFitBeam(const char* _Name, const Bool_t linkHistogram);
    ~GFitBeam()             {}

    virtual void    AddConstraintsIM();
    virtual void    AddConstraintMM();
    virtual void    CalcResult();
    virtual Int_t   Fill(Double_t x)    {return 0;}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {return 0;}
    virtual bool    Solve(const double time, const int channel, const GTreeA2Geant &geantTree);
    virtual void    SetBeam(const double _Beam)    {aplconBeam.SetFromBeam(_Beam);}
};




class	GFitBeamVertex    : public GFitBeam
{
private:
    GH1     vertex;

    void    ConvertTheta(std::vector<double>& p, std::vector<double>& v)  {p[1] = std::atan2( 25.4*sin(p[1]), 25.4*cos(p[1]) - v[0]);}

    virtual int GetNDOF()   {return 21;}

public:
    GFitBeamVertex(const char* _Name, const Bool_t linkHistogram);
    ~GFitBeamVertex()             {}

    virtual void    AddConstraintsIM();
    virtual void    AddConstraintMM();
    virtual void    CalcResult();
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual bool    Solve(const double time, const int channel, const GTreeA2Geant& geantTree);
};



#endif
