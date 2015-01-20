#ifndef _GKinFitter_h
#define _GKinFitter_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "GKinFitterConverter.h"


#define     GKinFitter_ParametersPerParticle    3


/*
class GKinFitter4VectorX: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

        return x[3]*helpX/r;
    }
    unsigned int NDim() const                               {return 7;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitter4VectorX();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));
        switch(ipar)
        {
        case 0:
            return x[3]*(helpY + helpZ)/(r*r*r);
        case 1:
            return -x[3]*helpX*helpY/(r*r*r);
        case 2:
            return -x[3]*helpX*helpZ/(r*r*r);
        case 3:
            return helpX/r;
        case 4:
            return x[3]*(helpY + helpZ)/(r*r*r);
        case 5:
            return x[3]*helpX*helpY/(r*r*r);
        case 6:
            return x[3]*helpX*helpZ/(r*r*r);
        default:
            return DoEval(x);
        }
    }
};
class GKinFitter4VectorY: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

        return x[3]*helpY/r;
    }
    unsigned int NDim() const                               {return 7;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitter4VectorY();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));
        switch(ipar)
        {
        case 0:
            return -x[3]*helpY*helpX/(r*r*r);
        case 1:
            return x[3]*(helpX + helpZ)/(r*r*r);
        case 2:
            return -x[3]*helpY*helpZ/(r*r*r);
        case 3:
            return helpY/r;
        case 4:
            return x[3]*helpY*helpX/(r*r*r);
        case 5:
            return x[3]*(helpX + helpZ)/(r*r*r);
        case 6:
            return x[3]*helpY*helpZ/(r*r*r);
        default:
            return DoEval(x);
        }
    }
};
class GKinFitter4VectorZ: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

        return x[3]*helpZ/r;
    }
    unsigned int NDim() const                               {return 7;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitter4VectorZ();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));
        switch(ipar)
        {
        case 0:
            return -x[3]*helpZ*helpX/(r*r*r);
        case 1:
            return -x[3]*helpZ*helpY/(r*r*r);
        case 2:
            return x[3]*(helpX + helpY)/(r*r*r);
        case 3:
            return helpZ/r;
        case 4:
            return x[3]*helpZ*helpX/(r*r*r);
        case 5:
            return x[3]*helpZ*helpY/(r*r*r);
        case 6:
            return x[3]*(helpX*helpY)/(r*r*r);
        default:
            return DoEval(x);
        }
    }
};*/
/*
class GKinFitterConstrainInvMass: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        Double_t    helpX = x[0] + x[4];
        Double_t    helpY = x[1] + x[5];
        Double_t    helpZ = x[2] + x[6];
        Double_t    helpE = x[3] + x[7];

        return (helpE*helpE)-(helpX*helpX)-(helpY*helpY)-(helpZ*helpZ)-(x[8]*x[8]);
    }
    unsigned int NDim() const                               {return 9;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitterConstrainInvMass();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        if(ipar==8)
            return -2*x[8];
        switch(ipar)
        {
        case 0:
        case 4:
            return -2*(x[0] + x[4]);
        case 1:
        case 5:
            return -2*(x[1] + x[5]);
        case 2:
        case 6:
            return -2*(x[2] + x[6]);
        case 3:
        case 7:
            return 2*(x[3] + x[7]);
        }
    }
};
class GKinFitterConstrainMisMass: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        Double_t    helpX = x[0] - x[4] - x[8]  - x[12] - x[16] - x[20] - x[24];
        Double_t    helpY = x[1] - x[5] - x[9]  - x[13] - x[17] - x[21] - x[25];
        Double_t    helpZ = x[2] - x[6] - x[10] - x[14] - x[18] - x[22] - x[26];
        Double_t    helpE = x[3] - x[7] - x[11] - x[15] - x[19] - x[23] - x[27];

        return (helpE*helpE)-(helpX*helpX)-(helpY*helpY)-(helpZ*helpZ)-(x[28]*x[28]);
    }
    unsigned int NDim() const                               {return 29;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitterConstrainMisMass();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        if(ipar==28)
            return 2*x[28];
        switch(ipar)
        {
        case 0:
        case 4:
        case 8:
        case 12:
        case 16:
        case 20:
        case 24:
            return -2*(x[0] - x[4] - x[8]  - x[12] - x[16] - x[20] - x[24]);
        case 1:
        case 5:
        case 9:
        case 13:
        case 17:
        case 21:
        case 25:
            return -2*(x[1] - x[5] - x[9]  - x[13] - x[17] - x[21] - x[25]);
        case 2:
        case 6:
        case 10:
        case 14:
        case 18:
        case 22:
        case 26:
            return -2*(x[2] - x[6] - x[10] - x[14] - x[18] - x[22] - x[26]);
        case 3:
        case 7:
        case 11:
        case 15:
        case 19:
        case 23:
        case 27:
            return 2*(x[3] - x[7] - x[11] - x[15] - x[19] - x[23] - x[27]);
        }
    }
};
class GKinFitterConstrainBeam: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        return ((x[3]-938.272046)*(x[3]-938.272046))-(x[0]*x[0])-(x[1]*x[1])-(x[2]*x[2]);
    }
    unsigned int NDim() const                               {return 4;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitterConstrainBeam();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        switch(ipar)
        {
        case 0:
            return -2*x[0];
        case 1:
            return -2*x[1];
        case 2:
            return -2*x[2];
        case 3:
            return 2*(x[3]-938.272046);
        default:
            DoEval(x);
        }
    }
};*/

class   GKinFitterNewTheta
{
public:
    GKinFitterNewTheta()    {}
    ~GKinFitterNewTheta()   {}

            Double_t    Eval(const Double_t theta, const Double_t zVertex)              {return TMath::ATan2(TMath::Sin(theta) ,(4*zVertex) + TMath::Cos(theta));}
    inline  Double_t    DerivateTheta(const Double_t theta, const Double_t zVertex);
    inline  Double_t    DerivateZVertex(const Double_t theta, const Double_t zVertex);
};
Double_t    GKinFitterNewTheta::DerivateTheta(const Double_t theta, const Double_t zVertex)
{
    Double_t    cosTh       = TMath::Cos(theta);
    Double_t    denominator = 16*zVertex*zVertex;
    denominator += 4*zVertex*cosTh;
    denominator += 1;
    return ((4*cosTh*zVertex) + 1)/denominator;
}
Double_t    GKinFitterNewTheta::DerivateZVertex(const Double_t theta, const Double_t zVertex)
{
    Double_t    sinTh       = TMath::Sin(theta);
    Double_t    denominator = (4*zVertex) + TMath::Cos(theta);
    denominator *= denominator;
    denominator += (sinTh*sinTh);
    return -4*sinTh/denominator;
}

class   GKinFitterGamma
{
private:
    inline TLorentzVector  RawEval(const Double_t energy, const Double_t theta, const Double_t phi);
    inline TLorentzVector  RawDerivateEnergy(const Double_t energy, const Double_t theta, const Double_t phi);
    inline TLorentzVector  RawDerivateTheta(const Double_t energy, const Double_t theta, const Double_t phi);
    inline TLorentzVector  RawDerivatePhi(const Double_t energy, const Double_t theta, const Double_t phi);

public:
    GKinFitterGamma()    {}
    ~GKinFitterGamma()   {}

    inline TLorentzVector  WithUnknownEval(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi);
    inline TLorentzVector  WithUnknownDerivateEnergy(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi);
    inline TLorentzVector  WithUnknownDerivateTheta(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi);
    inline TLorentzVector  WithUnknownDerivateZVertex(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi);
    inline TLorentzVector  WithUnknownDerivatePhi(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi);
    inline TLorentzVector  Beam(const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  BeamDerivateEnergy(const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  BeamDerivateTheta(const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  BeamDerivatePhi(const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  BeamDerivateZVertex(const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  Particle(const int particleIndex, const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  ParticleDerivateEnergy(const int particleIndex, const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  ParticleDerivateTheta(const int particleIndex, const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  ParticleDerivatePhi(const int particleIndex, const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector  ParticleDerivateZVertex(const int particleIndex, const TMatrixD& par, const TMatrixD& unk);
};
TLorentzVector  GKinFitterGamma::RawEval(const Double_t energy, const Double_t theta, const Double_t phi)
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
TLorentzVector  GKinFitterGamma::RawDerivateEnergy(const Double_t energy, const Double_t theta, const Double_t phi)
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
TLorentzVector  GKinFitterGamma::RawDerivateTheta(const Double_t energy, const Double_t theta, const Double_t phi)
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
TLorentzVector  GKinFitterGamma::RawDerivatePhi(const Double_t energy, const Double_t theta, const Double_t phi)
{
    Double_t sinTheta   = TMath::Sin(theta);
    Double_t sinPhi     = TMath::Sin(phi);
    Double_t cosPhi     = TMath::Cos(phi);
    return TLorentzVector(-energy * sinTheta * sinPhi,
                          energy * sinTheta * cosPhi,
                          0,
                          0);
}
TLorentzVector  GKinFitterGamma::WithUnknownEval(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi)
{
    return RawEval(energy, GKinFitterNewTheta().Eval(theta, zVertex), phi);
}
TLorentzVector  GKinFitterGamma::WithUnknownDerivateEnergy(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi)
{
    return RawDerivateEnergy(energy, GKinFitterNewTheta().Eval(theta, zVertex), phi);
}
TLorentzVector  GKinFitterGamma::WithUnknownDerivateTheta(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi)
{
    return RawDerivateTheta(energy, GKinFitterNewTheta().Eval(theta, zVertex), phi) * GKinFitterNewTheta().DerivateTheta(theta, zVertex);
}
TLorentzVector  GKinFitterGamma::WithUnknownDerivateZVertex(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi)
{
    return RawDerivateTheta(energy, GKinFitterNewTheta().Eval(theta, zVertex), phi) * GKinFitterNewTheta().DerivateZVertex(theta, zVertex);
}
TLorentzVector  GKinFitterGamma::WithUnknownDerivatePhi(const Double_t energy, const Double_t theta, const Double_t zVertex, const Double_t phi)
{
    return RawDerivatePhi(energy, GKinFitterNewTheta().Eval(theta, zVertex), phi);
}
TLorentzVector  GKinFitterGamma::Beam(const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownEval(par[0][0], par[1][0], unk[0][0], par[2][0]);
}
TLorentzVector  GKinFitterGamma::BeamDerivateEnergy(const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateEnergy(par[0][0], par[1][0], unk[0][0], par[2][0]);
}
TLorentzVector  GKinFitterGamma::BeamDerivateTheta(const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateTheta(par[0][0], par[1][0], unk[0][0], par[2][0]);
}
TLorentzVector  GKinFitterGamma::BeamDerivatePhi(const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivatePhi(par[0][0], par[1][0], unk[0][0], par[2][0]);
}
TLorentzVector  GKinFitterGamma::BeamDerivateZVertex(const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateZVertex(par[0][0], par[1][0], unk[0][0], par[2][0]);
}
TLorentzVector  GKinFitterGamma::Particle(const int particleIndex, const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownEval(par[ (particleIndex+1)*GKinFitter_ParametersPerParticle   ][0],
                           par[((particleIndex+1)*GKinFitter_ParametersPerParticle)+1][0],
                           unk[0][0],
                           par[((particleIndex+1)*GKinFitter_ParametersPerParticle)+2][0]);
}
TLorentzVector  GKinFitterGamma::ParticleDerivateEnergy(const int particleIndex, const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateEnergy(par[ (particleIndex+1)*GKinFitter_ParametersPerParticle   ][0],
                                     par[((particleIndex+1)*GKinFitter_ParametersPerParticle)+1][0],
                                     unk[0][0],
                                     par[((particleIndex+1)*GKinFitter_ParametersPerParticle)+2][0]);
}
TLorentzVector  GKinFitterGamma::ParticleDerivateTheta(const int particleIndex, const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateTheta(par[ (particleIndex+1)*GKinFitter_ParametersPerParticle   ][0],
                                    par[((particleIndex+1)*GKinFitter_ParametersPerParticle)+1][0],
                                    unk[0][0],
                                    par[((particleIndex+1)*GKinFitter_ParametersPerParticle)+2][0]);
}
TLorentzVector  GKinFitterGamma::ParticleDerivatePhi(const int particleIndex, const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivatePhi(par[ (particleIndex+1)*GKinFitter_ParametersPerParticle   ][0],
                                  par[((particleIndex+1)*GKinFitter_ParametersPerParticle)+1][0],
                                  unk[0][0],
                                  par[((particleIndex+1)*GKinFitter_ParametersPerParticle)+2][0]);
}
TLorentzVector  GKinFitterGamma::ParticleDerivateZVertex(const int particleIndex, const TMatrixD& par, const TMatrixD& unk)
{
    return WithUnknownDerivateZVertex(par[ (particleIndex+1)*GKinFitter_ParametersPerParticle   ][0],
                                      par[((particleIndex+1)*GKinFitter_ParametersPerParticle)+1][0],
                                      unk[0][0],
                                      par[((particleIndex+1)*GKinFitter_ParametersPerParticle)+2][0]);
}
class   GKinFitterInvMass
{
public:
    GKinFitterInvMass()    {}
    ~GKinFitterInvMass()   {}

    inline Double_t Eval(const TMatrixD& par, const TMatrixD& unk, const int* partList, const int nPartList, const Double_t mass);
    inline void     DerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk, const int* partList, const int nPartList);
    inline Double_t DerivateUnk(const TMatrixD& par, const TMatrixD& unk, const int* partList, const int nPartList);
};
Double_t    GKinFitterInvMass::Eval(const TMatrixD& par, const TMatrixD& unk, const int* partList, const int nPartList, const Double_t mass)
{
    TLorentzVector    tot(0.0, 0.0, 0.0, 0.0);
    for(int i=0; i<nPartList; i++)
        tot += GKinFitterGamma().Particle(partList[i], par, unk);
    return tot.M2()-(mass*mass);
}
void        GKinFitterInvMass::DerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk, const int* partList, const int nPartList)
{
    TLorentzVector    tot(0.0, 0.0, 0.0, 0.0);
    for(int i=0; i<nPartList; i++)
        tot += GKinFitterGamma().Particle(partList[i], par, unk);

    for(int i=0; i<nPartList; i++)
    {
        ret[0][ (partList[i]+1)*GKinFitter_ParametersPerParticle   ]    = 2*tot*GKinFitterGamma().ParticleDerivateEnergy(partList[i], par, unk);
        ret[0][((partList[i]+1)*GKinFitter_ParametersPerParticle)+1]    = 2*tot*GKinFitterGamma().ParticleDerivateTheta(partList[i], par, unk);
        ret[0][((partList[i]+1)*GKinFitter_ParametersPerParticle)+2]    = 2*tot*GKinFitterGamma().ParticleDerivatePhi(partList[i], par, unk);
    }
}
Double_t    GKinFitterInvMass::DerivateUnk(const TMatrixD& par, const TMatrixD& unk, const int* partList, const int nPartList)
{
    TLorentzVector    tot(0.0, 0.0, 0.0, 0.0);
    for(int i=0; i<nPartList; i++)
        tot += GKinFitterGamma().Particle(partList[i], par, unk);
    Double_t    ret = 0.0;
    for(int i=0; i<nPartList; i++)
        ret += 2*tot*GKinFitterGamma().ParticleDerivateZVertex(partList[i], par, unk);
    return  ret;
}




class   GKinFitter
{
public:
    enum    GKinFitterFitType   //value is type of unknown
    {
        flagNoRecoil        =1,
        flagRecoilAngles    =2,
        flagUnknownRecoil   =4
    };
private:
    Int_t               nPart; //Number of particles
    Int_t               nUnk; //Number of unknown parameters
    Int_t               nCon; //Number of constraints
    GKinFitterFitType   fitType;
    Int_t               nPar; //Number of parameters
    Double_t            targetMass;
    Int_t               countPart; //Number of particles
    Int_t               countCon; //Number of constraints
    Int_t               countIter; // Number of times Solve has been called
    TMatrixD            par0;      //vector of measured parameters
    TMatrixD            par;       //vector of fitted parameters
    TMatrixD            unk0;      //vector of unknowns calculated from constrain and par0
    TMatrixD            unk;       //vector of fitted unknowns
    TMatrixD            lambda;  	//Vector of lagrangian multipliers
    TMatrixD            V0;        //Covariance matrix of measured parameters
    TMatrixD            V;         //Covariance matrix of fitted parameters
    Double_t            chiSq;

    Int_t   GetNParameters();

    inline TLorentzVector   Recoil(const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector   RecoilDerivateEnergy(const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector   RecoilDerivateTheta(const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector   RecoilDerivatePhi(const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector   RecoilDerivateZVertex(const TMatrixD& par, const TMatrixD& unk);
    inline Double_t         MM(const TMatrixD& par, const TMatrixD& unk, const Double_t mass);
    inline void             MMDerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk);
    inline Double_t         MMDerivateUnk(const TMatrixD& par, const TMatrixD& unk);
    inline TLorentzVector   I4Vec(const TMatrixD& par, const TMatrixD& unk, const TLorentzVector& vec);
    inline void             I4VecDerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk);
    inline void             I4VecDerivateUnk(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk);
    inline void             g(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk);
    inline void             gDerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk);
    inline void             gDerivateUnk(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk);
    inline void             r(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk, const TMatrixD& G);

    Bool_t  SolveStep(TMatrixD &par, TMatrixD &unk, const TMatrixD& oldPar, const TMatrixD& oldUnk);
    Bool_t  CalcChiSq(const TMatrixD& par, const TMatrixD& unk, const TMatrixD &g);

protected:

public:
    GKinFitter(const Int_t nParticles, const Int_t nConstraints, const GKinFitterFitType type);
    ~GKinFitter();

    void            AddBeam(const Double_t beamEnergy, const Double_t _targetMass, const Double_t beamEnergyError, const Double_t beamSpotRadius);
    void            AddGamma(const Double_t energy, const Double_t theta, const Double_t phi, const Double_t energyError, const Double_t thetaError, const Double_t phiError);
    void            AddRecoilAngles(const Double_t theta, const Double_t phi, const Double_t thetaError, const Double_t phiError);
    void            Print(const char* option = "");
    void            Reset()     {countPart=0; countCon=0; countIter=0;}
    Bool_t          Solve();
    //Bool_t          ReSolve();*/
};
TLorentzVector  GKinFitter::Recoil(const TMatrixD& par, const TMatrixD& unk)
{
    switch(fitType)
    {
    case flagNoRecoil:
        return TLorentzVector(0.0, 0.0, 0.0, 0.0);
    case flagRecoilAngles:
        return GKinFitterGamma().WithUnknownEval(unk[1][0], par[(nPart+1)*GKinFitter_ParametersPerParticle][0], unk[0][0], par[((nPart+1)*GKinFitter_ParametersPerParticle)+1][0]);
    case flagUnknownRecoil:
        return GKinFitterGamma().WithUnknownEval(unk[1][0], unk[2][0], unk[0][0], unk[3][0]);
    }
}
TLorentzVector  GKinFitter::RecoilDerivateEnergy(const TMatrixD& par, const TMatrixD& unk)
{
    switch(fitType)
    {
    case flagNoRecoil:
        return TLorentzVector(0.0, 0.0, 0.0, 0.0);
    case flagRecoilAngles:
        return GKinFitterGamma().WithUnknownDerivateEnergy(unk[1][0], par[(nPart+1)*GKinFitter_ParametersPerParticle][0], unk[0][0], par[((nPart+1)*GKinFitter_ParametersPerParticle)+1][0]);
    case flagUnknownRecoil:
        return GKinFitterGamma().WithUnknownDerivateEnergy(unk[1][0], unk[2][0], unk[0][0], unk[3][0]);
    }
}
TLorentzVector  GKinFitter::RecoilDerivateTheta(const TMatrixD& par, const TMatrixD& unk)
{
    switch(fitType)
    {
    case flagNoRecoil:
        return TLorentzVector(0.0, 0.0, 0.0, 0.0);
    case flagRecoilAngles:
        return GKinFitterGamma().WithUnknownDerivateTheta(unk[1][0], par[(nPart+1)*GKinFitter_ParametersPerParticle][0], unk[0][0], par[((nPart+1)*GKinFitter_ParametersPerParticle)+1][0]);
    case flagUnknownRecoil:
        return GKinFitterGamma().WithUnknownDerivateTheta(unk[1][0], unk[2][0], unk[0][0], unk[3][0]);
    }
}
TLorentzVector  GKinFitter::RecoilDerivatePhi(const TMatrixD& par, const TMatrixD& unk)
{
    switch(fitType)
    {
    case flagNoRecoil:
        return TLorentzVector(0.0, 0.0, 0.0, 0.0);
    case flagRecoilAngles:
        return GKinFitterGamma().WithUnknownDerivatePhi(unk[1][0], par[(nPart+1)*GKinFitter_ParametersPerParticle][0], unk[0][0], par[((nPart+1)*GKinFitter_ParametersPerParticle)+1][0]);
    case flagUnknownRecoil:
        return GKinFitterGamma().WithUnknownDerivatePhi(unk[1][0], unk[2][0], unk[0][0], unk[3][0]);
    }
}
TLorentzVector  GKinFitter::RecoilDerivateZVertex(const TMatrixD& par, const TMatrixD& unk)
{
    switch(fitType)
    {
    case flagNoRecoil:
        return TLorentzVector(0.0, 0.0, 0.0, 0.0);
    case flagRecoilAngles:
        return GKinFitterGamma().WithUnknownDerivateZVertex(unk[1][0], par[(nPart+1)*GKinFitter_ParametersPerParticle][0], unk[0][0], par[((nPart+1)*GKinFitter_ParametersPerParticle)+1][0]);
    case flagUnknownRecoil:
        return GKinFitterGamma().WithUnknownDerivateZVertex(unk[1][0], unk[2][0], unk[0][0], unk[3][0]);
    }
}
Double_t        GKinFitter::MM(const TMatrixD& par, const TMatrixD& unk, const Double_t mass)
{
    TLorentzVector    tot(GKinFitterGamma().Beam(par, unk));
    for(int i=0; i<nPart; i++)
        tot -= GKinFitterGamma().Particle(i, par, unk);
    return tot.M2()-(mass*mass);
}
void            GKinFitter::MMDerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    TLorentzVector    tot(GKinFitterGamma().Beam(par, unk));
    for(int i=0; i<nPart; i++)
        tot -= GKinFitterGamma().Particle(i, par, unk);


    ret[0][0]    = 2*tot*GKinFitterGamma().BeamDerivateEnergy(par, unk);
    ret[0][1]    = 2*tot*GKinFitterGamma().BeamDerivateTheta(par, unk);
    ret[0][2]    = 2*tot*GKinFitterGamma().BeamDerivatePhi(par, unk);

    for(int i=0; i<nPart; i++)
    {
        ret[0][ (i+1)*GKinFitter_ParametersPerParticle   ]    = 2*tot*GKinFitterGamma().ParticleDerivateEnergy(i, par, unk);
        ret[0][((i+1)*GKinFitter_ParametersPerParticle)+1]    = 2*tot*GKinFitterGamma().ParticleDerivateTheta(i, par, unk);
        ret[0][((i+1)*GKinFitter_ParametersPerParticle)+2]    = 2*tot*GKinFitterGamma().ParticleDerivatePhi(i, par, unk);
    }
}
Double_t        GKinFitter::MMDerivateUnk(const TMatrixD& par, const TMatrixD& unk)
{
    TLorentzVector    tot(GKinFitterGamma().Beam(par, unk));
    for(int i=0; i<nPart; i++)
        tot -= GKinFitterGamma().Particle(i, par, unk);


    Double_t    ret    = 2*tot*GKinFitterGamma().BeamDerivateZVertex(par, unk);

    for(int i=0; i<nPart; i++)
        ret   += 2*tot*GKinFitterGamma().ParticleDerivateZVertex(i, par, unk);

    return ret;
}
TLorentzVector  GKinFitter::I4Vec(const TMatrixD& par, const TMatrixD& unk, const TLorentzVector &vec)
{
    TLorentzVector    tot(GKinFitterGamma().Beam(par, unk));
    for(int i=0; i<nPart; i++)
        tot -= GKinFitterGamma().Particle(i, par, unk);
    tot -= Recoil(par, unk);
    tot -= vec;
    return tot;
}
void            GKinFitter::I4VecDerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    TLorentzVector  help(GKinFitterGamma().BeamDerivateEnergy(par, unk));
    ret[0][0]   = help.Px();
    ret[1][0]   = help.Py();
    ret[2][0]   = help.Pz();
    ret[3][0]   = help.E();

    help    = GKinFitterGamma().BeamDerivateTheta(par, unk);
    ret[0][1]   = help.Px();
    ret[1][1]   = help.Py();
    ret[2][1]   = help.Pz();
    ret[3][1]   = help.E();

    help    = GKinFitterGamma().BeamDerivatePhi(par, unk);
    ret[0][2]   = help.Px();
    ret[1][2]   = help.Py();
    ret[2][2]   = help.Pz();
    ret[3][2]   = help.E();

    for(int i=0; i<nPart; i++)
    {
        help = GKinFitterGamma().ParticleDerivateEnergy(i, par, unk);
        ret[0][ (i+1)*GKinFitter_ParametersPerParticle   ]   = help.Px();
        ret[1][ (i+1)*GKinFitter_ParametersPerParticle   ]   = help.Py();
        ret[2][ (i+1)*GKinFitter_ParametersPerParticle   ]   = help.Pz();
        ret[3][ (i+1)*GKinFitter_ParametersPerParticle   ]   = help.E();

        help = GKinFitterGamma().ParticleDerivateTheta(i, par, unk);
        ret[0][((i+1)*GKinFitter_ParametersPerParticle)+1]   = help.Px();
        ret[1][((i+1)*GKinFitter_ParametersPerParticle)+1]   = help.Py();
        ret[2][((i+1)*GKinFitter_ParametersPerParticle)+1]   = help.Pz();
        ret[3][((i+1)*GKinFitter_ParametersPerParticle)+1]   = help.E();

        help = GKinFitterGamma().ParticleDerivatePhi(i, par, unk);
        ret[0][((i+1)*GKinFitter_ParametersPerParticle)+2]   = help.Px();
        ret[1][((i+1)*GKinFitter_ParametersPerParticle)+2]   = help.Py();
        ret[2][((i+1)*GKinFitter_ParametersPerParticle)+2]   = help.Pz();
        ret[3][((i+1)*GKinFitter_ParametersPerParticle)+2]   = help.E();
    }

    if(fitType==flagRecoilAngles)
    {
        help = RecoilDerivateTheta(par, unk);
        ret[0][ (nPart+1)*GKinFitter_ParametersPerParticle   ]   = help.Px();
        ret[1][ (nPart+1)*GKinFitter_ParametersPerParticle   ]   = help.Py();
        ret[2][ (nPart+1)*GKinFitter_ParametersPerParticle   ]   = help.Pz();
        ret[3][ (nPart+1)*GKinFitter_ParametersPerParticle   ]   = help.E();

        help = RecoilDerivatePhi(par, unk);
        ret[0][((nPart+1)*GKinFitter_ParametersPerParticle)+1]   = help.Px();
        ret[1][((nPart+1)*GKinFitter_ParametersPerParticle)+1]   = help.Py();
        ret[2][((nPart+1)*GKinFitter_ParametersPerParticle)+1]   = help.Pz();
        ret[3][((nPart+1)*GKinFitter_ParametersPerParticle)+1]   = help.E();
    }
}
void            GKinFitter::I4VecDerivateUnk(TMatrixD &ret, const TMatrixD& par, const TMatrixD& unk)
{
    TLorentzVector  help(GKinFitterGamma().BeamDerivateZVertex(par, unk));
    ret[0][0]   = help.Px();
    ret[1][0]   = help.Py();
    ret[2][0]   = help.Pz();
    ret[3][0]   = help.E();

    for(int i=0; i<nPart; i++)
    {
        help = GKinFitterGamma().ParticleDerivateZVertex(i, par, unk);
        ret[0][0]   += help.Px();
        ret[1][0]   += help.Py();
        ret[2][0]   += help.Pz();
        ret[3][0]   += help.E();
    }

    if(fitType==flagRecoilAngles)
    {
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
    else if(fitType==flagUnknownRecoil)
    {
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

        help    = RecoilDerivateTheta(par, unk);
        ret[0][2]   = help.Px();
        ret[1][2]   = help.Py();
        ret[2][2]   = help.Pz();
        ret[3][2]   = help.E();

        help    = RecoilDerivatePhi(par, unk);
        ret[0][3]   = help.Px();
        ret[1][3]   = help.Py();
        ret[2][3]   = help.Pz();
        ret[3][3]   = help.E();
    }
}
void            GKinFitter::g(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    int    indices[2];
    indices[0]  = 0;
    indices[1]  = 1;
    ret[0][0]   = GKinFitterInvMass().Eval(par, unk, indices, 2, 547);
    indices[0]  = 2;
    indices[1]  = 3;
    ret[1][0]   = GKinFitterInvMass().Eval(par, unk, indices, 2, 135);
    indices[0]  = 4;
    indices[1]  = 5;
    ret[2][0]   = GKinFitterInvMass().Eval(par, unk, indices, 2, 135);

    switch(fitType)
    {
    case flagNoRecoil:
        ret[3][0]   = MM(par, unk, 938);
        break;
    case flagUnknownRecoil:
    case flagRecoilAngles:
        {
        TLorentzVector  help(I4Vec(par, unk, TLorentzVector(0.0, 0.0, 0.0, 0.0)));
        ret[3][0]   = help.Px();
        ret[4][0]   = help.Py();
        ret[5][0]   = help.Pz();
        ret[6][0]   = help.E();
        }
        break;
    }
}
void            GKinFitter::gDerivatePar(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    int    indices[2];
    TMatrixD    help(4, nPart);
    indices[0]  = 0;
    indices[1]  = 1;
    GKinFitterInvMass().DerivatePar(help, par, unk, indices, 2);
    for(int i=0; i<nPar; i++)
        ret[0][i]   = help[0][i];
    indices[0]  = 2;
    indices[1]  = 3;
    GKinFitterInvMass().DerivatePar(help, par, unk, indices, 2);
    for(int i=0; i<nPar; i++)
        ret[1][i]   = help[0][i];
    indices[0]  = 4;
    indices[1]  = 5;
    GKinFitterInvMass().DerivatePar(help, par, unk, indices, 2);
    for(int i=0; i<nPar; i++)
        ret[2][i]   = help[0][i];

    switch(fitType)
    {
    case flagNoRecoil:
        MMDerivatePar(*((TMatrixD*)&ret.GetSub(3, 3, 0, nPart-1, help)), par, unk);
        for(int i=0; i<nPar; i++)
            ret[3][i]   = help[0][i];
        break;
    case flagUnknownRecoil:
    case flagRecoilAngles:
        MMDerivatePar(*((TMatrixD*)&ret.GetSub(3, 6, 0, nPart-1, help)), par, unk);
        for(int i=0; i<nPart; i++)
        {
            ret[3][i]   = help[0][i];
            ret[4][i]   = help[1][i];
            ret[5][i]   = help[2][i];
            ret[6][i]   = help[3][i];
        }
        break;
    }
}
void            GKinFitter::gDerivateUnk(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk)
{
    int    indices[2];
    TMatrixD    help;
    indices[0]  = 0;
    indices[1]  = 1;
    ret[0][0]   = GKinFitterInvMass().DerivateUnk(par, unk, indices, 2);
    indices[0]  = 2;
    indices[1]  = 3;
    ret[1][0]   = GKinFitterInvMass().DerivateUnk(par, unk, indices, 2);
    indices[0]  = 4;
    indices[1]  = 5;
    ret[2][0]   = GKinFitterInvMass().DerivateUnk(par, unk, indices, 2);

    switch(fitType)
    {
    case flagNoRecoil:
        ret[3][0]   = MMDerivateUnk(par, unk);
        break;
    case flagUnknownRecoil:
        ret[0][1]   = 0;
        ret[1][1]   = 0;
        ret[2][1]   = 0;
        ret[0][2]   = 0;
        ret[1][2]   = 0;
        ret[2][2]   = 0;
        ret[0][3]   = 0;
        ret[1][3]   = 0;
        ret[2][3]   = 0;
        I4VecDerivateUnk(help, par, unk);
        for(int i=0; i<nUnk; i++)
        {
            ret[3][i]   = help[0][i];
            ret[4][i]   = help[1][i];
            ret[5][i]   = help[2][i];
            ret[6][i]   = help[3][i];
        }
        break;
    case flagRecoilAngles:
        ret[0][1]   = 0;
        ret[1][1]   = 0;
        ret[2][1]   = 0;
        I4VecDerivateUnk(help, par, unk);
        for(int i=0; i<nUnk; i++)
        {
            ret[3][i]   = help[0][i];
            ret[4][i]   = help[1][i];
            ret[5][i]   = help[2][i];
            ret[6][i]   = help[3][i];
        }
        break;
    }
}
void            GKinFitter::r(TMatrixD& ret, const TMatrixD& par, const TMatrixD& unk, const TMatrixD& G)
{
    ret = G;
    TMatrixD    gPar(nCon, nPar);
    gDerivatePar(gPar, par, unk);
    TMatrixD    diff(par0);
    diff    -= par;
    ret     += gPar*diff;
}

#endif

