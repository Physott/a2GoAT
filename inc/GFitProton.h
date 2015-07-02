#ifndef __GFitProton_h__
#define __GFitProton_h__


#include "GFit.h"



class	GFitProton    : public GFit
{
private:
    double              aplconProtonTheta;
    double              aplconProtonPhi;
    double              aplconProtonThetaSigma;
    double              aplconProtonPhiSigma;

    GHistBGSub          protonEnergy;
    GHistBGSub          protonTheta;
    GHistBGSub          protonPhi;

    virtual int GetNDOF()   {return 21;}
    std::vector<double*>    GetLink()       {return {&aplconProtonTheta, &aplconProtonPhi};}
    std::vector<double*>    GetLinkSigma()  {return {&aplconProtonThetaSigma, &aplconProtonPhiSigma};}

public:
    GFitProton(const char* _Name, const Bool_t linkHistogram);
    virtual ~GFitProton()               {}

    virtual void    AddConstraintsTotMomentum();
    virtual void    AddConstraintsTotEnergy();
    virtual void    CalcResult();
    inline  const GFit::FitParticle& GetFittedProton() const;
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

const   GFit::FitParticle&    GFitProton::GetFittedProton() const
{
    const APLCON::Result_Variable_t& pe = result.Variables.at("PE");
    static FitParticle p;
    p.Ek            = pe.Value.After;
    p.Ek_Sigma      = pe.Sigma.After;
    p.Theta         = aplconProtonTheta;
    p.Theta_Sigma   = aplconProtonThetaSigma;
    p.Phi         = aplconProtonPhi;
    p.Phi_Sigma   = aplconProtonPhiSigma;
    return p;
}

class	GFitProtonVertex    : public GFitProton
{
private:
    GH1     vertex;

    void    ConvertTheta(std::vector<double>& p, std::vector<double>& v)  {p[1] = std::atan2( 25.4*sin(p[1]), 25.4*cos(p[1]) - v[0]);}

    virtual int GetNDOF()   {return 21;}

public:
    GFitProtonVertex(const char* _Name, const Bool_t linkHistogram);
    ~GFitProtonVertex()             {}

    virtual void    AddConstraintsIM();
    virtual void    AddConstraintMM();
    virtual void    AddConstraintsTotMomentum();
    virtual void    AddConstraintsTotEnergy();
    virtual void    CalcResult();
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual bool    Solve(const double time, const int channel, const GTreeA2Geant& geantTree);
};


#endif
