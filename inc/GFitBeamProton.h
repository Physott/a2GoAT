#ifndef __GFitBeamProton_h__
#define __GFitBeamProton_h__


#include "GFitBeam.h"



class	GFitBeamProton    : public GFitBeam
{
private:
    double              aplconProtonTheta;
    double              aplconProtonPhi;
    double              aplconProtonThetaSigma;
    double              aplconProtonPhiSigma;

    GH1                 protonEnergy;
    GH1                 protonTheta;
    GHistBGSub          protonPhi;

    virtual int GetNDOF()   {return 24;}
    std::vector<double*>    GetLink()       {return {&aplconProtonTheta, &aplconProtonPhi};}
    std::vector<double*>    GetLinkSigma()  {return {&aplconProtonThetaSigma, &aplconProtonPhiSigma};}


public:
    GFitBeamProton(const char* _Name, const Bool_t linkHistogram);
    virtual ~GFitBeamProton()               {}

    virtual void    AddConstraintsTotMomentum();
    virtual void    AddConstraintsTotEnergy();
    virtual void    CalcResult();
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual void    SetProton(const TLorentzVector proton){
        aplconProtonTheta = proton.Theta();
        aplconProtonPhi = proton.Phi();
        aplconProtonThetaSigma = 2.5*TMath::DegToRad();
        if(aplconProtonTheta>20*TMath::DegToRad() && aplconProtonTheta<160*TMath::DegToRad()) {
            aplconProtonPhiSigma = aplconProtonThetaSigma/sin(aplconProtonTheta);
        }
        else {
            aplconProtonPhiSigma = 1*TMath::DegToRad();
        }
    }
    virtual bool    Solve(const double time, const int channel, const GTreeA2Geant &geantTree);
};





class	GFitBeamProtonVertex    : public GFitBeamProton
{
private:
    GH1     vertex;

    void    ConvertTheta(std::vector<double>& p, std::vector<double>& v)  {p[1] = std::atan2( 25.4*sin(p[1]), 25.4*cos(p[1]) - v[0]);}

    virtual int GetNDOF()   {return 24;}

public:
    GFitBeamProtonVertex(const char* _Name, const Bool_t linkHistogram);
    ~GFitBeamProtonVertex()             {}

    virtual void    AddConstraintsIM();
    virtual void    AddConstraintMM();
    virtual void    AddConstraintsTotMomentum();
    virtual void    AddConstraintsTotEnergy();
    virtual void    CalcResult();
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual bool    Solve(const double time, const int channel, const GTreeA2Geant &geantTree);
};


#endif
