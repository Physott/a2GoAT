#ifndef __GFit_h__
#define __GFit_h__


#include <TLorentzVector.h>

#include "GH1.h"
#include "APLCON.hpp"

class   GTreeA2Geant;


class	GFit    : public GHistLinked
{
protected:
    // lightweight structure for linking to fitter
    struct FitParticle {
        void SetFromBeam(const Double_t beam) {
            Ek = beam;
            Theta = 0;
            Phi = 0;
            Ek_Sigma = 1;
            Theta_Sigma = 1*TMath::DegToRad();
            Phi_Sigma = 10;
        }
        void SetFromVector(const TLorentzVector& p_) {
            Ek = p_.E()-p_.M();
            Theta = p_.Theta();
            Phi = p_.Phi();
            Ek_Sigma = 0.01*Ek*pow(Ek,-0.36);
            Theta_Sigma = 2.5*TMath::DegToRad();
            if(Theta>20*TMath::DegToRad() && Theta<160*TMath::DegToRad()) {
                Phi_Sigma = Theta_Sigma/sin(Theta);
            }
            else {
                Phi_Sigma = 1*TMath::DegToRad();
            }
        }
        void SetToProtonResolution()
        {
            Ek_Sigma = 0.5*Ek*pow(Ek,-0.36);
        }
        static TLorentzVector Make(const std::vector<double>& EkThetaPhi,
                                   const Double_t m){
            const double E = EkThetaPhi[0] + m;
            const Double_t p = sqrt( E*E - m*m );
            TVector3 pv(1,0,0);
            pv.SetMagThetaPhi(p, EkThetaPhi[1], EkThetaPhi[2]);
            TLorentzVector l(pv, E);
            return l;
        }
        static TLorentzVector Make(const Double_t Ek,
                                   const std::vector<double>& ThetaPhi,
                                   const Double_t m){
            const double E = Ek + m;
            const Double_t p = sqrt( E*E - m*m );
            TVector3 pv(1,0,0);
            pv.SetMagThetaPhi(p, ThetaPhi[0], ThetaPhi[1]);
            TLorentzVector l(pv, E);
            return l;
        }
        static TLorentzVector Make(const Double_t Ek,
                                   const Double_t Theta,
                                   const Double_t Phi,
                                   const Double_t m){
            const double E = Ek + m;
            const Double_t p = sqrt( E*E - m*m );
            TVector3 pv(1,0,0);
            pv.SetMagThetaPhi(p, Theta, Phi);
            TLorentzVector l(pv, E);
            return l;
        }
        static TLorentzVector Make(const FitParticle& p,
                                   const Double_t m) {
            return Make(std::vector<double>{p.Ek, p.Theta, p.Phi}, m);
        }
        std::vector<double*> Link() {
            return {std::addressof(Ek),
                        std::addressof(Theta),
                        std::addressof(Phi)};
        }
        std::vector<double*> LinkSigma() {
            return {std::addressof(Ek_Sigma),
                        std::addressof(Theta_Sigma),
                        std::addressof(Phi_Sigma)};
        }


        double Ek;
        double Ek_Sigma;
        double Theta;
        double Theta_Sigma;
        double Phi;
        double Phi_Sigma;

    };

    APLCON                      fitter;
    Double_t                    beam;
    APLCON::Result_t            result;

    GHistBGSub2             pulls;

private:
    std::string                 name;
    std::vector<FitParticle>    aplconPhotons;

    bool                    success;

    GHistBGSub              steps;
    GH1                     im;
    GHistBGSub              sub0Im;
    GHistBGSub              sub1Im;
    GHistBGSub              sub2Im;
    GHistBGSub              theta;
    GH1                     thetaCM;
    GH1                     thetaCMTrue;
    GHistBGSub              CMtest;
    GHistBGSub              phi;
    GH1                     chiSq;
    GH1                     confidenceLevel;

    virtual int GetNDOF()   {return 18;}

public:
    GFit(const char* _Name, const Bool_t linkHistogram);
    ~GFit()             {}

    virtual void    AddConstraintsIM();
    virtual void    AddConstraintMM();
    virtual void    CalcResult();
    virtual Int_t   Fill(Double_t x)    {return 0;}
    inline  const FitParticle& GetFittedPhoton(const int i) const;
            bool    IsSuccess() const   {return success;}
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)   {return 0;}
    virtual void    Set(const TLorentzVector& p0,
                        const TLorentzVector& p1,
                        const TLorentzVector& p2,
                        const TLorentzVector& p3,
                        const TLorentzVector& p4,
                        const TLorentzVector& p5)
    {
        aplconPhotons[0].SetFromVector(p0);
        aplconPhotons[1].SetFromVector(p1);
        aplconPhotons[2].SetFromVector(p2);
        aplconPhotons[3].SetFromVector(p3);
        aplconPhotons[4].SetFromVector(p4);
        aplconPhotons[5].SetFromVector(p5);
    }
    virtual void    SetBeam(const Double_t _Beam)    {beam=_Beam;}
    virtual bool    Solve(const double time, const int channel, const GTreeA2Geant& geantTree);
};

const GFit::FitParticle& GFit::GetFittedPhoton(const int i) const
{
    try
    {
        return  aplconPhotons[i];
    }
    catch(...)
    {
        std::cout << "There is no " << i << "th photon in GFit!" << std::endl;
        return aplconPhotons[0];
    }
}



class	GFitVertex    : public GFit
{
private:
    GH1     vertex;

    void    ConvertTheta(std::vector<double>& p, std::vector<double>& v)  {p[1] = std::atan2( 25.4*sin(p[1]), 25.4*cos(p[1]) - v[0]);}

    virtual int GetNDOF()   {return 18;}

public:
    GFitVertex(const char* _Name, const Bool_t linkHistogram);
    ~GFitVertex()             {}

    virtual void    AddConstraintsIM();
    virtual void    AddConstraintMM();
    virtual void    CalcResult();
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual bool    Solve(const double time, const int channel, const GTreeA2Geant &geantTree);
};



#endif
