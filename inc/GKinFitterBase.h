#ifndef _GKinFitterBase_h
#define _GKinFitterBase_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "GHistFit.h"


#define     GKinFitterBase_ParametersPerParticle    3



class   GKinFitterBaseNewTheta
{
public:
    GKinFitterBaseNewTheta()    {}
    ~GKinFitterBaseNewTheta()   {}

            Double_t    Eval(const Double_t theta, const Double_t zVertex)              {return TMath::ATan2(TMath::Sin(theta) ,(4*zVertex) + TMath::Cos(theta));}
    inline  Double_t    DerivateTheta(const Double_t theta, const Double_t zVertex);
    inline  Double_t    DerivateZVertex(const Double_t theta, const Double_t zVertex);
};
Double_t    GKinFitterBaseNewTheta::DerivateTheta(const Double_t theta, const Double_t zVertex)
{
    Double_t    cosTh       = TMath::Cos(theta);
    Double_t    denominator = 16*zVertex*zVertex;
    denominator += 4*zVertex*cosTh;
    denominator += 1;
    return ((4*cosTh*zVertex) + 1)/denominator;
}
Double_t    GKinFitterBaseNewTheta::DerivateZVertex(const Double_t theta, const Double_t zVertex)
{
    Double_t    sinTh       = TMath::Sin(theta);
    Double_t    denominator = (4*zVertex) + TMath::Cos(theta);
    denominator *= denominator;
    denominator += (sinTh*sinTh);
    return -4*sinTh/denominator;
}

class   GKinFitterBaseGamma
{
private:
    inline TLorentzVector  RawEval(const Double_t energy, const Double_t theta, const Double_t phi);
    inline TLorentzVector  RawDerivateEnergy(const Double_t energy, const Double_t theta, const Double_t phi);
    inline TLorentzVector  RawDerivateTheta(const Double_t energy, const Double_t theta, const Double_t phi);
    inline TLorentzVector  RawDerivatePhi(const Double_t energy, const Double_t theta, const Double_t phi);

public:
    GKinFitterBaseGamma()    {}
    ~GKinFitterBaseGamma()   {}

    inline TLorentzVector  GammaEval(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi);
    inline TLorentzVector  GammaDerivateEnergy(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi);
    inline TLorentzVector  GammaDerivateTheta(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi);
    inline TLorentzVector  GammaDerivateZVertex(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi);
    inline TLorentzVector  GammaDerivatePhi(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi);
/*           TLorentzVector  BeamAndTarget(const TMatrixD& par, const TMatrixD& unk, const Double_t targetMass)  {return Beam(par, unk)+TLorentzVector(0.0, 0.0, 0.0, targetMass);}
           TLorentzVector  BeamAndTargetDerivateEnergy(const TMatrixD& par, const TMatrixD& unk)               {return BeamDerivateEnergy(par, unk);}
           TLorentzVector  BeamAndTargetDerivateTheta(const TMatrixD& par, const TMatrixD& unk)                {return BeamDerivateTheta(par, unk);}
           TLorentzVector  BeamAndTargetDerivatePhi(const TMatrixD& par, const TMatrixD& unk)                  {return BeamDerivatePhi(par, unk);}
           TLorentzVector  BeamAndTargetDerivateZVertex(const TMatrixD& par, const TMatrixD& unk)              {return BeamDerivateZVertex(par, unk);}
    inline TLorentzVector  Particle(const int particleIndex, const TMatrixD& par, const Double_t unk);
    inline TLorentzVector  ParticleDerivateEnergy(const int particleIndex, const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  ParticleDerivateTheta(const int particleIndex, const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  ParticleDerivatePhi(const int particleIndex, const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  ParticleDerivateZVertex(const int particleIndex, const TMatrixD& par, const TMatrixD& unk);*/
};
TLorentzVector  GKinFitterBaseGamma::RawEval(const Double_t energy, const Double_t theta, const Double_t phi)
{
    Double_t sinTheta   = TMath::Sin(theta);
    Double_t cosTheta   = TMath::Cos(theta);
    Double_t sinPhi     = TMath::Sin(phi);
    Double_t cosPhi     = TMath::Cos(phi);
    return TLorentzVector(energy * sinTheta * cosPhi,
                          energy * sinTheta * sinPhi,
                          energy * cosTheta,
                          energy);
}
TLorentzVector  GKinFitterBaseGamma::RawDerivateEnergy(const Double_t energy, const Double_t theta, const Double_t phi)
{
    Double_t sinTheta   = TMath::Sin(theta);
    Double_t cosTheta   = TMath::Cos(theta);
    Double_t sinPhi     = TMath::Sin(phi);
    Double_t cosPhi     = TMath::Cos(phi);
    return TLorentzVector(sinTheta * cosPhi,
                          sinTheta * sinPhi,
                          cosTheta,
                          1);
}
TLorentzVector  GKinFitterBaseGamma::RawDerivateTheta(const Double_t energy, const Double_t theta, const Double_t phi)
{
    Double_t sinTheta   = TMath::Sin(theta);
    Double_t cosTheta   = TMath::Cos(theta);
    Double_t sinPhi     = TMath::Sin(phi);
    Double_t cosPhi     = TMath::Cos(phi);
    return TLorentzVector(energy * cosTheta * cosPhi,
                          energy * cosTheta * sinPhi,
                          -energy * sinTheta,
                          0);
}
TLorentzVector  GKinFitterBaseGamma::RawDerivatePhi(const Double_t energy, const Double_t theta, const Double_t phi)
{
    Double_t sinTheta   = TMath::Sin(theta);
    Double_t sinPhi     = TMath::Sin(phi);
    Double_t cosPhi     = TMath::Cos(phi);
    return TLorentzVector(-energy * sinTheta * sinPhi,
                          energy * sinTheta * cosPhi,
                          0,
                          0);
}
TLorentzVector  GKinFitterBaseGamma::GammaEval(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi)
{
    return RawEval(energy, GKinFitterBaseNewTheta().Eval(theta, zVertex), phi);
}
TLorentzVector  GKinFitterBaseGamma::GammaDerivateEnergy(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi)
{
    return RawDerivateEnergy(energy, GKinFitterBaseNewTheta().Eval(theta, zVertex), phi);
}
TLorentzVector  GKinFitterBaseGamma::GammaDerivateTheta(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi)
{
    return RawDerivateTheta(energy, GKinFitterBaseNewTheta().Eval(theta, zVertex), phi) * GKinFitterBaseNewTheta().DerivateTheta(theta, zVertex);
}
TLorentzVector  GKinFitterBaseGamma::GammaDerivateZVertex(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi)
{
    return RawDerivateTheta(energy, GKinFitterBaseNewTheta().Eval(theta, zVertex), phi) * GKinFitterBaseNewTheta().DerivateZVertex(theta, zVertex);
}
TLorentzVector  GKinFitterBaseGamma::GammaDerivatePhi(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi)
{
    return RawDerivatePhi(energy, GKinFitterBaseNewTheta().Eval(theta, zVertex), phi);
}
/*TLorentzVector  GKinFitterBaseGamma::Beam(const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownEval(par[0][0], par[1][0], unk[0][0], par[2][0]);
}
TLorentzVector  GKinFitterBaseGamma::BeamDerivateEnergy(const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateEnergy(par[0][0], par[1][0], unk[0][0], par[2][0]);
}
TLorentzVector  GKinFitterBaseGamma::BeamDerivateTheta(const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateTheta(par[0][0], par[1][0], unk[0][0], par[2][0]);
}
TLorentzVector  GKinFitterBaseGamma::BeamDerivatePhi(const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivatePhi(par[0][0], par[1][0], unk[0][0], par[2][0]);
}
TLorentzVector  GKinFitterBaseGamma::BeamDerivateZVertex(const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateZVertex(par[0][0], par[1][0], unk[0][0], par[2][0]);
}
TLorentzVector  GKinFitterBaseGamma::Particle(const int particleIndex, const TMatrixD& par, const Double_t unk)
{
    return WithUnknownEval(par[ (particleIndex+1)*GKinFitterBase_ParametersPerParticle   ][0],
                           par[((particleIndex+1)*GKinFitterBase_ParametersPerParticle)+1][0],
                           unk,
                           par[((particleIndex+1)*GKinFitterBase_ParametersPerParticle)+2][0]);
}
TLorentzVector  GKinFitterBaseGamma::ParticleDerivateEnergy(const int particleIndex, const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateEnergy(par[ (particleIndex+1)*GKinFitterBase_ParametersPerParticle   ][0],
                                     par[((particleIndex+1)*GKinFitterBase_ParametersPerParticle)+1][0],
                                     unk[0][0],
                                     par[((particleIndex+1)*GKinFitterBase_ParametersPerParticle)+2][0]);
}
TLorentzVector  GKinFitterBaseGamma::ParticleDerivateTheta(const int particleIndex, const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateTheta(par[ (particleIndex+1)*GKinFitterBase_ParametersPerParticle   ][0],
                                    par[((particleIndex+1)*GKinFitterBase_ParametersPerParticle)+1][0],
                                    unk[0][0],
                                    par[((particleIndex+1)*GKinFitterBase_ParametersPerParticle)+2][0]);
}
TLorentzVector  GKinFitterBaseGamma::ParticleDerivatePhi(const int particleIndex, const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivatePhi(par[ (particleIndex+1)*GKinFitterBase_ParametersPerParticle   ][0],
                                  par[((particleIndex+1)*GKinFitterBase_ParametersPerParticle)+1][0],
                                  unk[0][0],
                                  par[((particleIndex+1)*GKinFitterBase_ParametersPerParticle)+2][0]);
}
TLorentzVector  GKinFitterBaseGamma::ParticleDerivateZVertex(const int particleIndex, const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateZVertex(par[ (particleIndex+1)*GKinFitterBase_ParametersPerParticle   ][0],
                                      par[((particleIndex+1)*GKinFitterBase_ParametersPerParticle)+1][0],
                                      unk[0][0],
                                      par[((particleIndex+1)*GKinFitterBase_ParametersPerParticle)+2][0]);
}*/
/*class   GKinFitterBaseInvMass
{
public:
    GKinFitterBaseInvMass()    {}
    ~GKinFitterBaseInvMass()   {}

    inline Double_t Eval(const TMatrixD& par, const Double_t unk, const int* partList, const int nPartList, const Double_t mass);
    inline void     DerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk, const int* partList, const int nPartList);
    inline Double_t DerivateUnk(const TMatrixD& par, const TMatrixD& unk, const int* partList, const int nPartList);
};
Double_t    GKinFitterBaseInvMass::Eval(const TMatrixD& par, const Double_t unk, const int* partList, const int nPartList, const Double_t mass)
{
    TLorentzVector    tot(0.0, 0.0, 0.0, 0.0);
    for(int i=0; i<nPartList; i++)
        tot += GKinFitterBaseGamma().GammaEval(partList[i], par, unk);
    return tot.M2()-(mass*mass);
}
void        GKinFitterBaseInvMass::DerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk, const int* partList, const int nPartList)
{
    TLorentzVector    tot(0.0, 0.0, 0.0, 0.0);
    for(int i=0; i<nPartList; i++)
        tot += GKinFitterBaseGamma().Particle(partList[i], par, unk);

    for(int i=0; i<nPartList; i++)
    {
        ret[0][ (partList[i]+1)*GKinFitterBase_ParametersPerParticle   ]    = 2*tot*GKinFitterBaseGamma().ParticleDerivateEnergy(partList[i], par, unk);
        ret[0][((partList[i]+1)*GKinFitterBase_ParametersPerParticle)+1]    = 2*tot*GKinFitterBaseGamma().ParticleDerivateTheta(partList[i], par, unk);
        ret[0][((partList[i]+1)*GKinFitterBase_ParametersPerParticle)+2]    = 2*tot*GKinFitterBaseGamma().ParticleDerivatePhi(partList[i], par, unk);
    }
}
Double_t    GKinFitterBaseInvMass::DerivateUnk(const TMatrixD& par, const TMatrixD& unk, const int* partList, const int nPartList)
{
    TLorentzVector    tot(0.0, 0.0, 0.0, 0.0);
    for(int i=0; i<nPartList; i++)
        tot += GKinFitterBaseGamma().Particle(partList[i], par, unk);
    Double_t    ret = 0.0;
    for(int i=0; i<nPartList; i++)
        ret += 2*tot*GKinFitterBaseGamma().ParticleDerivateZVertex(partList[i], par, unk);
    return  ret;
}*/


struct  GKinFitterBaseIMConstraint
{
    int         nPartList;
    int         PartList[10];
    Double_t    mass;
};

class   GKinFitterBase
{
private:
    Int_t               nPart; //Number of particles
    Int_t               nCon; //Number of constraints
    Int_t               nPar; //Number of parameters
    Double_t            targetMass;
    Int_t               countPart; //Number of particles
    Int_t               countCon; //Number of constraints
    Int_t               countIter; // Number of times Solve has been called
    TMatrixD            par0;      //vector of measured parameters
    TMatrixD            par1;       //vector of fitted parameters
    TMatrixD            par2;       //vector of fitted parameters
    Double_t            unk1;      //vector of unknowns calculated from constrain and par0
    Double_t            unk2;       //vector of fitted unknowns
    TMatrixD            lambda;  	//Vector of lagrangian multipliers
    TMatrixD            V0;        //Covariance matrix of measured parameters
    TMatrixD            V;         //Covariance matrix of fitted parameters
    Double_t            chiSq;

    //Constraints
    int                         nImConstraint;
    GKinFitterBaseIMConstraint  imConstraint[10];
    Bool_t                      isMmConstraint;
    Double_t                    mmConstraint;

    //Matrices for Calculation
    TMatrixD            G;
    TMatrixD            GPar;
    TMatrixD            GUnk;
    TMatrixD            SInv;
    Double_t            U;

    TLorentzVector  GetBeam(const TMatrixD& par, const Double_t unk)                const   {return GKinFitterBaseGamma().GammaEval(par[0][0], par[1][0], unk, par[2][0]);}
    TLorentzVector  GetBeamDerivateEnergy(const TMatrixD& par, const Double_t unk)  const   {return GKinFitterBaseGamma().GammaDerivateEnergy(par[0][0], par[1][0], unk, par[2][0]);}
    TLorentzVector  GetBeamDerivateTheta(const TMatrixD& par, const Double_t unk)   const   {return GKinFitterBaseGamma().GammaDerivateTheta(par[0][0], par[1][0], unk, par[2][0]);}
    TLorentzVector  GetBeamDerivatePhi(const TMatrixD& par, const Double_t unk)     const   {return GKinFitterBaseGamma().GammaDerivatePhi(par[0][0], par[1][0], unk, par[2][0]);}
    TLorentzVector  GetBeamDerivateZVertex(const TMatrixD& par, const Double_t unk) const   {return GKinFitterBaseGamma().GammaDerivateZVertex(par[0][0], par[1][0], unk, par[2][0]);}
    TLorentzVector  GetPhoton(const TMatrixD& par, const Double_t unk, const int i)                 const   {return GKinFitterBaseGamma().GammaEval(par[(i+1)*GKinFitterBase_ParametersPerParticle][0], par[((i+1)*GKinFitterBase_ParametersPerParticle)+1][0], unk, par[((i+1)*GKinFitterBase_ParametersPerParticle)+2][0]);}
    TLorentzVector  GetPhotonDerivateEnergy(const TMatrixD& par, const Double_t unk, const int i)   const   {return GKinFitterBaseGamma().GammaDerivateEnergy(par[(i+1)*GKinFitterBase_ParametersPerParticle][0], par[((i+1)*GKinFitterBase_ParametersPerParticle)+1][0], unk, par[((i+1)*GKinFitterBase_ParametersPerParticle)+2][0]);}
    TLorentzVector  GetPhotonDerivateTheta(const TMatrixD& par, const Double_t unk, const int i)    const   {return GKinFitterBaseGamma().GammaDerivateTheta(par[(i+1)*GKinFitterBase_ParametersPerParticle][0], par[((i+1)*GKinFitterBase_ParametersPerParticle)+1][0], unk, par[((i+1)*GKinFitterBase_ParametersPerParticle)+2][0]);}
    TLorentzVector  GetPhotonDerivatePhi(const TMatrixD& par, const Double_t unk, const int i)      const   {return GKinFitterBaseGamma().GammaDerivatePhi(par[(i+1)*GKinFitterBase_ParametersPerParticle][0], par[((i+1)*GKinFitterBase_ParametersPerParticle)+1][0], unk, par[((i+1)*GKinFitterBase_ParametersPerParticle)+2][0]);}
    TLorentzVector  GetPhotonDerivateZVertex(const TMatrixD& par, const Double_t unk, const int i)  const   {return GKinFitterBaseGamma().GammaDerivateZVertex(par[(i+1)*GKinFitterBase_ParametersPerParticle][0], par[((i+1)*GKinFitterBase_ParametersPerParticle)+1][0], unk, par[((i+1)*GKinFitterBase_ParametersPerParticle)+2][0]);}

    Double_t        GetIMConstraint(const TMatrixD& par, const Double_t unk,const int i)                        const;
    void            GetIMConstraintDerPar(TMatrixD& ret, const TMatrixD& par, const Double_t unk,const int i);
    Double_t        GetIMConstraintDerUnk(const TMatrixD& par, const Double_t unk,const int i)                  const;
    Double_t        GetMMConstraint(const TMatrixD& par, const Double_t unk)                                    const;
    void            GetMMConstraintDerPar(TMatrixD& ret, const TMatrixD& par, const Double_t unk);
    Double_t        GetMMConstraintDerUnk(const TMatrixD& par, const Double_t unk)                              const;

    void    GetInitialG(TMatrixD& ret);
    void    GetInitialGDerPar(TMatrixD& ret);
    void    GetInitialGDerUnk(TMatrixD& ret);
    void    GetG(TMatrixD& ret, const TMatrixD &par, const Double_t unk);
    void    GetGDerPar(TMatrixD& ret, const TMatrixD& par, const Double_t unk);
    void    GetGDerUnk(TMatrixD& ret, const TMatrixD& par, const Double_t unk);

    void    GetR(TMatrixD& ret, const TMatrixD &par, const Double_t unk);


    Bool_t  CalcChiSq(const TMatrixD& par, const TMatrixD& unk);
    Bool_t  CalcV(const TMatrixD& par, const TMatrixD& unk);
    void    FillHists(GHistFit2& hists);
    Bool_t  SolveStep(GHistFit2& hists);
    Bool_t  SolveStart(GHistFit2& hists);

protected:

public:
    GKinFitterBase(const Int_t nParticles, const Int_t nConstraints);
    ~GKinFitterBase();

    TLorentzVector  GetBeam()               const   {return GetBeam(par2, unk2);}
    TLorentzVector  GetPhoton(const int i)  const   {return GetPhoton(par2, unk2, i);}

    void            AddInvMassConstraint(const int* partList, const int nPartList, const Double_t mass);
    void            AddMisMassConstraint(const Double_t mass);
    void            AddBeam(const Double_t beamEnergy, const Double_t _targetMass, const Double_t beamEnergyError, const Double_t beamSpotRadius);
    void            AddGamma(const Double_t energy, const Double_t theta, const Double_t phi, const Double_t energyError, const Double_t thetaError, const Double_t phiError);
    Double_t        ConfidenceLevel()           {return TMath::Prob(chiSq, nCon-1);}
    Double_t        Pull(const Int_t i)         {return (par0[i][0] - par2[i][0])/TMath::Sqrt(V0[i][i] - V[i][i]);}
    void            Print(const char* option = "");
    void            Reset()     {countPart=0; countCon=0; countIter=0; nImConstraint=0; isMmConstraint=kFALSE;}
    Bool_t          Solve(GHistFit2& hists);
    //Bool_t          ReSolve();*/
};
/*TLorentzVector  GKinFitterBase::Recoil(const TMatrixD& par, const TMatrixD& unk)
{
   return GKinFitterBaseGamma().WithUnknownEval(unk[1][0], par[(nPart+1)*GKinFitterBase_ParametersPerParticle][0], unk[0][0], par[((nPart+1)*GKinFitterBase_ParametersPerParticle)+1][0]);
}
TLorentzVector  GKinFitterBase::RecoilDerivateEnergy(const TMatrixD& par, const TMatrixD& unk)
{
    return GKinFitterBaseGamma().WithUnknownDerivateEnergy(unk[1][0], par[(nPart+1)*GKinFitterBase_ParametersPerParticle][0], unk[0][0], par[((nPart+1)*GKinFitterBase_ParametersPerParticle)+1][0]);
}
TLorentzVector  GKinFitterBase::RecoilDerivateTheta(const TMatrixD& par, const TMatrixD& unk)
{
    return GKinFitterBaseGamma().WithUnknownDerivateTheta(unk[1][0], par[(nPart+1)*GKinFitterBase_ParametersPerParticle][0], unk[0][0], par[((nPart+1)*GKinFitterBase_ParametersPerParticle)+1][0]);
}
TLorentzVector  GKinFitterBase::RecoilDerivatePhi(const TMatrixD& par, const TMatrixD& unk)
{
    return GKinFitterBaseGamma().WithUnknownDerivatePhi(unk[1][0], par[(nPart+1)*GKinFitterBase_ParametersPerParticle][0], unk[0][0], par[((nPart+1)*GKinFitterBase_ParametersPerParticle)+1][0]);
}
TLorentzVector  GKinFitterBase::RecoilDerivateZVertex(const TMatrixD& par, const TMatrixD& unk)
{
    return GKinFitterBaseGamma().WithUnknownDerivateZVertex(unk[1][0], par[(nPart+1)*GKinFitterBase_ParametersPerParticle][0], unk[0][0], par[((nPart+1)*GKinFitterBase_ParametersPerParticle)+1][0]);
}
Double_t        GKinFitterBase::MM(const TMatrixD& par, const TMatrixD& unk, const Double_t mass)
{
    TLorentzVector    tot(GKinFitterBaseGamma().BeamAndTarget(par, unk, targetMass));
    for(int i=0; i<nPart; i++)
        tot -= GKinFitterBaseGamma().Particle(i, par, unk);
    return tot.M2()-(mass*mass);
}
void            GKinFitterBase::MMDerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    TLorentzVector    tot(GKinFitterBaseGamma().BeamAndTarget(par, unk, targetMass));
    for(int i=0; i<nPart; i++)
        tot -= GKinFitterBaseGamma().Particle(i, par, unk);


    ret[0][0]    = 2*tot*GKinFitterBaseGamma().BeamAndTargetDerivateEnergy(par, unk);
    ret[0][1]    = 2*tot*GKinFitterBaseGamma().BeamAndTargetDerivateTheta(par, unk);
    ret[0][2]    = 2*tot*GKinFitterBaseGamma().BeamAndTargetDerivatePhi(par, unk);

    for(int i=0; i<nPart; i++)
    {
        ret[0][ (i+1)*GKinFitterBase_ParametersPerParticle   ]    = 2*tot*GKinFitterBaseGamma().ParticleDerivateEnergy(i, par, unk);
        ret[0][((i+1)*GKinFitterBase_ParametersPerParticle)+1]    = 2*tot*GKinFitterBaseGamma().ParticleDerivateTheta(i, par, unk);
        ret[0][((i+1)*GKinFitterBase_ParametersPerParticle)+2]    = 2*tot*GKinFitterBaseGamma().ParticleDerivatePhi(i, par, unk);
    }
}
Double_t        GKinFitterBase::MMDerivateUnk(const TMatrixD& par, const TMatrixD& unk)
{
    TLorentzVector    tot(GKinFitterBaseGamma().BeamAndTarget(par, unk, targetMass));
    for(int i=0; i<nPart; i++)
        tot -= GKinFitterBaseGamma().Particle(i, par, unk);


    Double_t    ret    = 2*tot*GKinFitterBaseGamma().BeamAndTargetDerivateZVertex(par, unk);

    for(int i=0; i<nPart; i++)
        ret   += 2*tot*GKinFitterBaseGamma().ParticleDerivateZVertex(i, par, unk);

    return ret;
}
TLorentzVector  GKinFitterBase::I4Vec(const TMatrixD& par, const TMatrixD& unk, const TLorentzVector &vec)
{
    TLorentzVector    tot(GKinFitterBaseGamma().BeamAndTarget(par, unk, targetMass));
    for(int i=0; i<nPart; i++)
        tot -= GKinFitterBaseGamma().Particle(i, par, unk);
    tot -= Recoil(par, unk);
    tot -= vec;
    return tot;
}
void            GKinFitterBase::I4VecDerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    TLorentzVector  help(GKinFitterBaseGamma().BeamAndTargetDerivateEnergy(par, unk));
    ret[0][0]   = help.Px();
    ret[1][0]   = help.Py();
    ret[2][0]   = help.Pz();
    ret[3][0]   = help.E();

    help    = GKinFitterBaseGamma().BeamAndTargetDerivateTheta(par, unk);
    ret[0][1]   = help.Px();
    ret[1][1]   = help.Py();
    ret[2][1]   = help.Pz();
    ret[3][1]   = help.E();

    help    = GKinFitterBaseGamma().BeamAndTargetDerivatePhi(par, unk);
    ret[0][2]   = help.Px();
    ret[1][2]   = help.Py();
    ret[2][2]   = help.Pz();
    ret[3][2]   = help.E();

    for(int i=0; i<nPart; i++)
    {
        help = GKinFitterBaseGamma().ParticleDerivateEnergy(i, par, unk);
        ret[0][ (i+1)*GKinFitterBase_ParametersPerParticle   ]   = help.Px();
        ret[1][ (i+1)*GKinFitterBase_ParametersPerParticle   ]   = help.Py();
        ret[2][ (i+1)*GKinFitterBase_ParametersPerParticle   ]   = help.Pz();
        ret[3][ (i+1)*GKinFitterBase_ParametersPerParticle   ]   = help.E();

        help = GKinFitterBaseGamma().ParticleDerivateTheta(i, par, unk);
        ret[0][((i+1)*GKinFitterBase_ParametersPerParticle)+1]   = help.Px();
        ret[1][((i+1)*GKinFitterBase_ParametersPerParticle)+1]   = help.Py();
        ret[2][((i+1)*GKinFitterBase_ParametersPerParticle)+1]   = help.Pz();
        ret[3][((i+1)*GKinFitterBase_ParametersPerParticle)+1]   = help.E();

        help = GKinFitterBaseGamma().ParticleDerivatePhi(i, par, unk);
        ret[0][((i+1)*GKinFitterBase_ParametersPerParticle)+2]   = help.Px();
        ret[1][((i+1)*GKinFitterBase_ParametersPerParticle)+2]   = help.Py();
        ret[2][((i+1)*GKinFitterBase_ParametersPerParticle)+2]   = help.Pz();
        ret[3][((i+1)*GKinFitterBase_ParametersPerParticle)+2]   = help.E();
    }

   help = RecoilDerivateTheta(par, unk);
   ret[0][ (nPart+1)*GKinFitterBase_ParametersPerParticle   ]   = help.Px();
   ret[1][ (nPart+1)*GKinFitterBase_ParametersPerParticle   ]   = help.Py();
   ret[2][ (nPart+1)*GKinFitterBase_ParametersPerParticle   ]   = help.Pz();
   ret[3][ (nPart+1)*GKinFitterBase_ParametersPerParticle   ]   = help.E();

  help = RecoilDerivatePhi(par, unk);
  ret[0][((nPart+1)*GKinFitterBase_ParametersPerParticle)+1]   = help.Px();
  ret[1][((nPart+1)*GKinFitterBase_ParametersPerParticle)+1]   = help.Py();
  ret[2][((nPart+1)*GKinFitterBase_ParametersPerParticle)+1]   = help.Pz();
  ret[3][((nPart+1)*GKinFitterBase_ParametersPerParticle)+1]   = help.E();
}
void            GKinFitterBase::I4VecDerivateUnk(TMatrixD &ret, const TMatrixD& par, const TMatrixD& unk)
{
    TLorentzVector  help(GKinFitterBaseGamma().BeamAndTargetDerivateZVertex(par, unk));
    ret[0][0]   = help.Px();
    ret[1][0]   = help.Py();
    ret[2][0]   = help.Pz();
    ret[3][0]   = help.E();

    for(int i=0; i<nPart; i++)
    {
        help = GKinFitterBaseGamma().ParticleDerivateZVertex(i, par, unk);
        ret[0][0]   += help.Px();
        ret[1][0]   += help.Py();
        ret[2][0]   += help.Pz();
        ret[3][0]   += help.E();
    }

    help = RecoilDerivateZVertex(par, unk);
    ret[0][0]   += help.Px();
    ret[1][0]   += help.Py();
    ret[2][0]   += help.Pz();
    ret[3][0]   += help.E();

    help    = RecoilDerivateEnergy(par, unk);
    ret[0][1]   = help.Px();
    ret[1][1]   = help.Py();
    ret[2][1]   = help.Pz();
    ret[3][1]   = help.E();

}
void            GKinFitterBase::g(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    int    indices[2];
    indices[0]  = 0;
    indices[1]  = 1;
    ret[0][0]   = GKinFitterBaseInvMass().Eval(par, unk, indices, 2, 547);
    indices[0]  = 2;
    indices[1]  = 3;
    ret[1][0]   = GKinFitterBaseInvMass().Eval(par, unk, indices, 2, 135);
    indices[0]  = 4;
    indices[1]  = 5;
    ret[2][0]   = GKinFitterBaseInvMass().Eval(par, unk, indices, 2, 135);

    TLorentzVector  help(I4Vec(par, unk, TLorentzVector(0.0, 0.0, 0.0, 0.0)));
    ret[3][0]   = help.Px();
    ret[4][0]   = help.Py();
    ret[5][0]   = help.Pz();
    ret[6][0]   = help.E();
}
void            GKinFitterBase::gDerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    int    indices[2];
    TMatrixD    help(4, nPar);
    indices[0]  = 0;
    indices[1]  = 1;
    for(int i=0; i<nPar; i++)
        help[0][i]   = 0;
    GKinFitterBaseInvMass().DerivatePar(help, par, unk, indices, 2);
    for(int i=0; i<nPar; i++)
        ret[0][i]   = help[0][i];
    indices[0]  = 2;
    indices[1]  = 3;
    for(int i=0; i<nPar; i++)
        help[0][i]   = 0;
    GKinFitterBaseInvMass().DerivatePar(help, par, unk, indices, 2);
    for(int i=0; i<nPar; i++)
        ret[1][i]   = help[0][i];
    indices[0]  = 4;
    indices[1]  = 5;
    for(int i=0; i<nPar; i++)
        help[0][i]   = 0;
    GKinFitterBaseInvMass().DerivatePar(help, par, unk, indices, 2);
    for(int i=0; i<nPar; i++)
        ret[2][i]   = help[0][i];

    I4VecDerivatePar(help, par, unk);
    for(int i=0; i<nPart; i++)
    {
        ret[3][i]   = help[0][i];
        ret[4][i]   = help[1][i];
        ret[5][i]   = help[2][i];
        ret[6][i]   = help[3][i];
    }
}
void            GKinFitterBase::gDerivateUnk(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    int    indices[2];
    TMatrixD    help(4, 1);
    indices[0]  = 0;
    indices[1]  = 1;
    ret[0][0]   = GKinFitterBaseInvMass().DerivateUnk(par, unk, indices, 2);
    indices[0]  = 2;
    indices[1]  = 3;
    ret[1][0]   = GKinFitterBaseInvMass().DerivateUnk(par, unk, indices, 2);
    indices[0]  = 4;
    indices[1]  = 5;
    ret[2][0]   = GKinFitterBaseInvMass().DerivateUnk(par, unk, indices, 2);

    ret[0][1]   = 0;
    ret[1][1]   = 0;
    ret[2][1]   = 0;
    for(int i=0; i<4; i++)
    {
        help[i][0]  = 0;
    }
    I4VecDerivateUnk(help, par, unk);
    ret[3][0]   = help[0][0];
    ret[4][0]   = help[1][0];
    ret[5][0]   = help[2][0];
    ret[6][0]   = help[3][0];
}
void            GKinFitterBase::r(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    ret = G;
    TMatrixD    gPar(nCon, nPar);
    gDerivatePar(gPar, par, unk);
    TMatrixD    diff(par0);
    diff    -= par;
    ret     += gPar*diff;
}*/

#endif

