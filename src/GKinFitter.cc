//////////////////////////////////////////////////////////////////
//GKinFitter
//////////////////////////////////////////////////////////////////


#include "GKinFitter.h"
#include "TDecompLU.h"
#include <iostream>


#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046






GKinFitter::GKinFitter(const Int_t nParticles, const Int_t nConstraints, const GKinFitterFitType type)    :
    nPart(nParticles),
    nUnk(Int_t(type)),
    nCon(nConstraints),
    fitType(type),
    nPar(GetNParameters()),
    countPart(nParticles),
    countCon(nConstraints),
    countIter(0),
    par0(nPar,1),
    par(nPar,1),
    unk0(nUnk,1),
    unk(nUnk,1),
    lambda(nCon,1),
    V0(nPar, nPar),
    V(nPar, nPar),
    //fmD(nCon,nPar),
    //fmd(nCon,1),
    //fmlamda(nCon,1),
    //fmV_D(nCon,nCon),
    chiSq(0)//,
    //fPtot(0, 0, 0, 0),
    //solved(kFALSE)
{
}

GKinFitter::~GKinFitter()
{

}

Int_t   GKinFitter::GetNParameters()
{
    switch(fitType)
    {
    case flagNoRecoil:
    case flagUnknownRecoil:
        return ((nPart+1)*GKinFitter_ParametersPerParticle);    //+1 for the beam
        return ((nPart+1)*GKinFitter_ParametersPerParticle);    //+1 for the beam
    case flagRecoilAngles:
        return ((nPart+1)*GKinFitter_ParametersPerParticle) + 2;    //+1 for the beam +2 for the recoil angles
    }
}


//public members


void            GKinFitter::AddBeam(const Double_t beamEnergy, const Double_t _targetMass, const Double_t beamEnergyError, const Double_t beamSpotRadius)
{
    targetMass  = _targetMass;
    par0[0][0]  = beamEnergy;
    par0[1][0]  = 0;
    par0[2][0]  = 0;
    V0[0][0]    = beamEnergyError*beamEnergyError;
    V0[1][1]    = beamSpotRadius*beamSpotRadius/100;
    V0[2][2]    = 2*TMath::Pi();
}
void            GKinFitter::AddGamma(const Double_t energy, const Double_t theta, const Double_t phi, const Double_t energyError, const Double_t thetaError, const Double_t phiError)
{
    if(countPart==nPart)
    {
        std::cout << "Can not add Particle. Already " << countPart << " there." << std::endl;
        return;
    }
    par0[ (countPart+1)*GKinFitter_ParametersPerParticle   ][0]  = energy;
    par0[((countPart+1)*GKinFitter_ParametersPerParticle)+1][0]  = theta;
    par0[((countPart+1)*GKinFitter_ParametersPerParticle)+2][0]  = phi;
    V0[ (countPart+1)*GKinFitter_ParametersPerParticle   ][ (countPart+1)*GKinFitter_ParametersPerParticle   ]    = energyError*energyError;
    V0[((countPart+1)*GKinFitter_ParametersPerParticle)+1][((countPart+1)*GKinFitter_ParametersPerParticle)+1]    = thetaError*thetaError;
    V0[((countPart+1)*GKinFitter_ParametersPerParticle)+2][((countPart+1)*GKinFitter_ParametersPerParticle)+2]    = phiError*phiError;
    countPart++;
}
void            GKinFitter::AddRecoilAngles(const Double_t theta, const Double_t phi, const Double_t thetaError, const Double_t phiError)
{
    par0[ (nPart+1)*GKinFitter_ParametersPerParticle   ][0]  = theta;
    par0[((nPart+1)*GKinFitter_ParametersPerParticle)+1][0]  = phi;
    V0[ (nPart+1)*GKinFitter_ParametersPerParticle   ][ (nPart+1)*GKinFitter_ParametersPerParticle   ]    = thetaError*thetaError;
    V0[((nPart+1)*GKinFitter_ParametersPerParticle)+1][((nPart+1)*GKinFitter_ParametersPerParticle)+1]    = phiError*phiError;
}

void            GKinFitter::Print(const char* option)
{
    std::cout << "GKinFitter: nPart=" << nPart << ", nUnk=" << nUnk << ", nCon=" << nCon;
    if(fitType==flagNoRecoil)
        std::cout << "   FitType: NoRecoil" << std::endl;
    else if(fitType==flagUnknownRecoil)
        std::cout << "   FitType: UnknownRecoil" << std::endl;
    else if(fitType==flagRecoilAngles)
        std::cout << "   FitType: RecoilAngles" << std::endl;

    TString str(option);
    str.ToLower();
    if(strcmp(str.Data(), "input")==0)
    {
        if(countPart<nPart)
            std::cout << "Input not completed. Only " << countPart << " Particles added yet." << std::endl;
        std::cout << "Input Parameters GKinFitter::par0     beam, particles, recoil angles(optional)" << std::endl;
        par0.Print();
        std::cout << "Input Covariance Matrix GKinFitter::V0     beam, particles, recoil angles(optional)" << std::endl;
        V0.Print();
    }
}

Bool_t          GKinFitter::SolveStep(TMatrixD& par, TMatrixD& unk, const TMatrixD& oldPar, const TMatrixD& oldUnk)
{
    TMatrixD    gPar(nCon, nPar);
    gDerivatePar(gPar, oldPar, oldUnk);
    TMatrixD    gParT(gPar);
    gParT.T();

    TMatrixD    S(gPar*V0*gParT);
    Double_t    determinant    = S.Determinant();
    if(determinant==0)
        return kFALSE;
    TMatrixD    SInv(S);
    SInv.Invert();


    TMatrixD    gUnk(nCon, nPar);
    gDerivateUnk(gUnk, oldPar, oldUnk);
    TMatrixD    gUnkT(gUnk);
    gUnkT.T();
    TMatrixD    T(gUnkT*SInv*gUnk);
    determinant    = T.Determinant();
    if(determinant==0)
        return kFALSE;
    TMatrixD    TInv(T);
    TInv.Invert();

    TMatrixD    G(nCon, 1);
    TMatrixD    R(nCon, 1);
    g(G, oldPar, oldUnk);
    r(R, oldPar, oldUnk, G);
    unk  = oldUnk;
    unk -= TInv*gUnkT*SInv*R;

    TMatrixD    diff(unk);
    diff    -= oldUnk;
    TMatrixD    help(R);
    help    += gUnk*diff;
    lambda   = SInv*help;

    par  = oldPar;
    par -= V0*gParT*lambda;

    if(CalcChiSq(par, unk, G)==kFALSE)
        return kFALSE;
    return kTRUE;
}

Bool_t          GKinFitter::Solve()
{
    //Calc Unk0
    TLorentzVector    tot(-GKinFitterGamma().Beam(par, unk));
    unk0[0][0]   = 0;
    switch(fitType)
    {
    case flagNoRecoil:
        break;
    case flagUnknownRecoil:
        for(int i=0; i<nPart; i++)
            tot += GKinFitterGamma().Particle(i, par, unk);
        unk0[0][1]   = tot.E();
        unk0[0][2]   = tot.Theta();
        unk0[0][3]   = tot.Phi();
        break;
    case flagRecoilAngles:
        for(int i=0; i<nPart; i++)
            tot += GKinFitterGamma().Particle(i, par, unk);
        unk[0][1]   = tot.E();
        break;
    }

    SolveStep(par, unk, par0, unk0);
    Double_t    oldChiSq    = chiSq;
    TMatrixD    newPar(par);
    TMatrixD    newUnk(unk);
    SolveStep(newPar, newUnk, par, unk);
    countIter   = 1;
    while(chiSq<oldChiSq)
    {
        oldChiSq    = chiSq;
        par         = newPar;
        unk         = newUnk;
        countIter++;
        SolveStep(newPar, newUnk, par, unk);
    }
}

Bool_t    GKinFitter::CalcChiSq(const TMatrixD& par, const TMatrixD& unk, const TMatrixD& g)
{
    TMatrixD    lambdaT(lambda);
    lambdaT.T();
    TMatrixD    res(lambdaT*g);
    res *= 2;

    TMatrixD    diff(par0);
    diff    -=  par;
    TMatrixD    diffT(diff);
    diffT.T();

    TMatrixD    V0Inv(V0);
    Double_t    determinant    = V0Inv.Determinant();
    if(determinant==0)
        return kFALSE;
    V0Inv.Invert();
    TMatrixD    help(diffT*V0Inv*diff);
    res += help;

    chiSq   = res[0][0];
    return kTRUE;
}

/*Double_t    GKinFitter::GetReactionPx(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    Double_t h[7] = {u[0][0],u[1][0],u[2][0], x[0][0],0,0,-10500};
    if(derivate_i<0)
        return GKinFitter4VectorX().DoEval(h);
    return GKinFitter4VectorX().DoDerivative(h, derivate_i);
}
Double_t    GKinFitter::GetReactionPy(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    Double_t h[7] = {u[0][0],u[1][0],u[2][0], x[0][0],0,0,-10500};
    if(derivate_i<0)
        return GKinFitter4VectorY().DoEval(h);
    return GKinFitter4VectorY().DoDerivative(h, derivate_i);
}
Double_t    GKinFitter::GetReactionPz(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    Double_t h[7] = {u[0][0],u[1][0],u[2][0], x[0][0],0,0,-10500};
    if(derivate_i<0)
        return GKinFitter4VectorZ().DoEval(h);
    return GKinFitter4VectorZ().DoDerivative(h, derivate_i);
}
Double_t    GKinFitter::GetPhotonPx(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    Double_t h[7] = {x[(4*i)+1][0],x[(4*i)+2][0],x[(4*i)+3][0], x[(4*i)+4][0],u[0][0],u[1][0],u[2][0]};
    if(derivate_i<0)
        return GKinFitter4VectorX().DoEval(h);
    return GKinFitter4VectorX().DoDerivative(h, derivate_i);
}
Double_t    GKinFitter::GetPhotonPy(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    Double_t h[7] = {x[(4*i)+1][0],x[(4*i)+2][0],x[(4*i)+3][0], x[(4*i)+4][0],u[0][0],u[1][0],u[2][0]};
    if(derivate_i<0)
        return GKinFitter4VectorY().DoEval(h);
    return GKinFitter4VectorY().DoDerivative(h, derivate_i);
}
Double_t    GKinFitter::GetPhotonPz(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    Double_t h[7] = {x[(4*i)+1][0],x[(4*i)+2][0],x[(4*i)+3][0], x[(4*i)+4][0],u[0][0],u[1][0],u[2][0]};
    if(derivate_i<0)
        return GKinFitter4VectorZ().DoEval(h);
    return GKinFitter4VectorZ().DoDerivative(h, derivate_i);
}
Double_t    GKinFitter::GetConstraint(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    switch(i)
    {
    case 0:
        {
            Double_t    h[9] = {GetPhotonPx(0, x, u), GetPhotonPy(0, x, u), GetPhotonPz(0, x, u), GetPhotonE(0, x),
                                GetPhotonPx(1, x, u), GetPhotonPy(1, x, u), GetPhotonPz(1, x, u), GetPhotonE(1, x), MASS_ETA};
            if(derivate_i<0)
                return GKinFitterConstrainInvMass().DoEval(h);
            else
                return GKinFitterConstrainInvMass().DoDerivative(h, derivate_i);
        }
    case 1:
        {
            Double_t    h[9] = {GetPhotonPx(2, x, u), GetPhotonPy(2, x, u), GetPhotonPz(2, x, u), GetPhotonE(2, x),
                                GetPhotonPx(3, x, u), GetPhotonPy(3, x, u), GetPhotonPz(3, x, u), GetPhotonE(3, x), MASS_PI0};
            if(derivate_i<0)
                return GKinFitterConstrainInvMass().DoEval(h);
            else
                return GKinFitterConstrainInvMass().DoDerivative(h, derivate_i);
        }
    case 2:
        {
            Double_t    h[9] = {GetPhotonPx(4, x, u), GetPhotonPy(4, x, u), GetPhotonPz(4, x, u), GetPhotonE(4, x),
                                GetPhotonPx(5, x, u), GetPhotonPy(5, x, u), GetPhotonPz(5, x, u), GetPhotonE(5, x), MASS_PI0};
            if(derivate_i<0)
                return GKinFitterConstrainInvMass().DoEval(h);
            else
                return GKinFitterConstrainInvMass().DoDerivative(h, derivate_i);
        }
    case 3:
        {
            Double_t    h[29] = {GetReactionPx(x, u),  GetReactionPy(x, u),  GetReactionPz(x, u),  GetReactionE(x),
                                 GetPhotonPx(0, x, u), GetPhotonPy(0, x, u), GetPhotonPz(0, x, u), GetPhotonE(0, x),
                                 GetPhotonPx(1, x, u), GetPhotonPy(1, x, u), GetPhotonPz(1, x, u), GetPhotonE(1, x),
                                 GetPhotonPx(2, x, u), GetPhotonPy(2, x, u), GetPhotonPz(2, x, u), GetPhotonE(2, x),
                                 GetPhotonPx(3, x, u), GetPhotonPy(3, x, u), GetPhotonPz(3, x, u), GetPhotonE(3, x),
                                 GetPhotonPx(4, x, u), GetPhotonPy(4, x, u), GetPhotonPz(4, x, u), GetPhotonE(4, x),
                                 GetPhotonPx(5, x, u), GetPhotonPy(5, x, u), GetPhotonPz(5, x, u), GetPhotonE(5, x)};
            if(derivate_i<0)
                return GKinFitterConstrainInvMass().DoEval(h);
            else
                return GKinFitterConstrainInvMass().DoDerivative(h, derivate_i);
        }
    case 4:
        {
            Double_t    h[9] = {GetReactionPx(x, u),  GetReactionPy(x, u),  GetReactionPz(x, u),  GetReactionE(x), MASS_PROTON};
            if(derivate_i<0)
                return GKinFitterConstrainBeam().DoEval(h);
            else
                return GKinFitterConstrainBeam().DoDerivative(h, derivate_i);
        }
    }
}

Double_t    GKinFitter::GetReactionPx_DeriPar(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    if(derivate_i==0)
        return GetReactionPx(x, u, 3);
    return 0;
}
Double_t    GKinFitter::GetReactionPy_DeriPar(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    if(derivate_i==0)
        return GetReactionPy(x, u, 3);
    return 0;
}
Double_t    GKinFitter::GetReactionPz_DeriPar(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    if(derivate_i==0)
        return GetReactionPz(x, u, 3);
    return 0;
}
Double_t    GKinFitter::GetPhotonPx_DeriPar(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    if(derivate_i<(4*i)+1)
        return 0;
    if(derivate_i>(4*i)+4)
        return 0;
    return GetPhotonPx(i, x, u, derivate_i-((4*i)+1));
}
Double_t    GKinFitter::GetPhotonPy_DeriPar(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    if(derivate_i<(4*i)+1)
        return 0;
    if(derivate_i>(4*i)+4)
        return 0;
    return GetPhotonPy(i, x, u, derivate_i-((4*i)+1));
}
Double_t    GKinFitter::GetPhotonPz_DeriPar(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    if(derivate_i<(4*i)+1)
        return 0;
    if(derivate_i>(4*i)+4)
        return 0;
    return GetPhotonPz(i, x, u, derivate_i-((4*i)+1));
}
Double_t    GKinFitter::GetConstraint_DeriPar(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    switch(i)
    {
    case 0:
        {
            Double_t    h[9] = {GetPhotonPx(0, x, u), GetPhotonPy(0, x, u), GetPhotonPz(0, x, u), GetPhotonE(0, x),
                                GetPhotonPx(1, x, u), GetPhotonPy(1, x, u), GetPhotonPz(1, x, u), GetPhotonE(1, x), MASS_ETA};
            Double_t    dh[8] = {GetPhotonPx_DeriPar(0, x, u, derivate_i), GetPhotonPy_DeriPar(0, x, u, derivate_i), GetPhotonPz_DeriPar(0, x, u, derivate_i), GetPhotonE_DeriPar(0, x, derivate_i),
                                 GetPhotonPx_DeriPar(1, x, u, derivate_i), GetPhotonPy_DeriPar(1, x, u, derivate_i), GetPhotonPz_DeriPar(1, x, u, derivate_i), GetPhotonE_DeriPar(1, x, derivate_i)};
            Double_t    ret = 0;
            for(int i=0; i<8; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            return ret;
        }
    case 1:
        {
            Double_t    h[9] = {GetPhotonPx(2, x, u), GetPhotonPy(2, x, u), GetPhotonPz(2, x, u), GetPhotonE(2, x),
                                GetPhotonPx(3, x, u), GetPhotonPy(3, x, u), GetPhotonPz(3, x, u), GetPhotonE(3, x), MASS_PI0};
            Double_t    dh[8] = {GetPhotonPx_DeriPar(2, x, u, derivate_i), GetPhotonPy_DeriPar(2, x, u, derivate_i), GetPhotonPz_DeriPar(2, x, u, derivate_i), GetPhotonE_DeriPar(2, x, derivate_i),
                                 GetPhotonPx_DeriPar(3, x, u, derivate_i), GetPhotonPy_DeriPar(3, x, u, derivate_i), GetPhotonPz_DeriPar(3, x, u, derivate_i), GetPhotonE_DeriPar(3, x, derivate_i)};
            Double_t    ret = 0;
            for(int i=0; i<8; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            return ret;
        }
    case 2:
        {
            Double_t    h[9] = {GetPhotonPx(4, x, u), GetPhotonPy(4, x, u), GetPhotonPz(4, x, u), GetPhotonE(4, x),
                                GetPhotonPx(5, x, u), GetPhotonPy(5, x, u), GetPhotonPz(5, x, u), GetPhotonE(5, x), MASS_PI0};
            Double_t    dh[8] = {GetPhotonPx_DeriPar(4, x, u, derivate_i), GetPhotonPy_DeriPar(4, x, u, derivate_i), GetPhotonPz_DeriPar(4, x, u, derivate_i), GetPhotonE_DeriPar(4, x, derivate_i),
                                 GetPhotonPx_DeriPar(5, x, u, derivate_i), GetPhotonPy_DeriPar(5, x, u, derivate_i), GetPhotonPz_DeriPar(5, x, u, derivate_i), GetPhotonE_DeriPar(5, x, derivate_i)};
            Double_t    ret = 0;
            for(int i=0; i<8; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            return ret;
        }
    case 3:
        {
            Double_t    h[29] = {GetReactionPx(x, u),  GetReactionPy(x, u),  GetReactionPz(x, u),  GetReactionE(x),
                                 GetPhotonPx(0, x, u), GetPhotonPy(0, x, u), GetPhotonPz(0, x, u), GetPhotonE(0, x),
                                 GetPhotonPx(1, x, u), GetPhotonPy(1, x, u), GetPhotonPz(1, x, u), GetPhotonE(1, x),
                                 GetPhotonPx(2, x, u), GetPhotonPy(2, x, u), GetPhotonPz(2, x, u), GetPhotonE(2, x),
                                 GetPhotonPx(3, x, u), GetPhotonPy(3, x, u), GetPhotonPz(3, x, u), GetPhotonE(3, x),
                                 GetPhotonPx(4, x, u), GetPhotonPy(4, x, u), GetPhotonPz(4, x, u), GetPhotonE(4, x),
                                 GetPhotonPx(5, x, u), GetPhotonPy(5, x, u), GetPhotonPz(5, x, u), GetPhotonE(5, x), MASS_PROTON};
            Double_t    dh[28] = {GetReactionPx_DeriPar(x, u, derivate_i), GetReactionPy_DeriPar(x, u, derivate_i), GetReactionPz_DeriPar(x, u, derivate_i), GetReactionE_DeriPar(x, derivate_i),
                                  GetPhotonPx_DeriPar(0, x, u, derivate_i), GetPhotonPy_DeriPar(0, x, u, derivate_i), GetPhotonPz_DeriPar(0, x, u, derivate_i), GetPhotonE_DeriPar(0, x, derivate_i),
                                  GetPhotonPx_DeriPar(1, x, u, derivate_i), GetPhotonPy_DeriPar(1, x, u, derivate_i), GetPhotonPz_DeriPar(1, x, u, derivate_i), GetPhotonE_DeriPar(1, x, derivate_i),
                                  GetPhotonPx_DeriPar(2, x, u, derivate_i), GetPhotonPy_DeriPar(2, x, u, derivate_i), GetPhotonPz_DeriPar(2, x, u, derivate_i), GetPhotonE_DeriPar(2, x, derivate_i),
                                  GetPhotonPx_DeriPar(3, x, u, derivate_i), GetPhotonPy_DeriPar(3, x, u, derivate_i), GetPhotonPz_DeriPar(3, x, u, derivate_i), GetPhotonE_DeriPar(3, x, derivate_i),
                                  GetPhotonPx_DeriPar(4, x, u, derivate_i), GetPhotonPy_DeriPar(4, x, u, derivate_i), GetPhotonPz_DeriPar(4, x, u, derivate_i), GetPhotonE_DeriPar(4, x, derivate_i),
                                  GetPhotonPx_DeriPar(5, x, u, derivate_i), GetPhotonPy_DeriPar(5, x, u, derivate_i), GetPhotonPz_DeriPar(5, x, u, derivate_i), GetPhotonE_DeriPar(5, x, derivate_i)};
            Double_t    ret = 0;
            for(int i=0; i<28; i++)
                ret += GKinFitterConstrainMisMass().DoDerivative(h, i) * dh[i];
            return ret;
        }
    case 4:
        {
            Double_t    h[9] = {GetReactionPx(x, u),  GetReactionPy(x, u),  GetReactionPz(x, u),  GetReactionE(x), MASS_PROTON};
            Double_t    dh[8] = {GetReactionPx_DeriPar(x, u, derivate_i), GetReactionPy_DeriPar(x, u, derivate_i), GetReactionPz_DeriPar(x, u, derivate_i), GetReactionE_DeriPar(x, derivate_i)};
            Double_t    ret = 0;
            for(int i=0; i<8; i++)
                ret += GKinFitterConstrainBeam().DoDerivative(h, i) * dh[i];
            return ret;
        }
    }
}

Double_t    GKinFitter::GetPhotonPx_DeriUnk(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    if(derivate_i>2)
        return 0;
    return GetPhotonPx(i, x, u, derivate_i+4);
}
Double_t    GKinFitter::GetPhotonPy_DeriUnk(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    if(derivate_i>2)
        return 0;
    return GetPhotonPy(i, x, u, derivate_i+4);
}
Double_t    GKinFitter::GetPhotonPz_DeriUnk(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    if(derivate_i>2)
        return 0;
    return GetPhotonPz(i, x, u, derivate_i+4);
}
Double_t    GKinFitter::GetConstraint_DeriUnk(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)
{
    switch(i)
    {
    case 0:
        {
            Double_t    h[9] = {GetPhotonPx(0, x, u), GetPhotonPy(0, x, u), GetPhotonPz(0, x, u), GetPhotonE(0, x),
                                GetPhotonPx(1, x, u), GetPhotonPy(1, x, u), GetPhotonPz(1, x, u), GetPhotonE(1, x), MASS_ETA};
            Double_t    dh[8] = {GetPhotonPx_DeriUnk(0, x, u, derivate_i), GetPhotonPy_DeriUnk(0, x, u, derivate_i), GetPhotonPz_DeriUnk(0, x, u, derivate_i), 0,
                                 GetPhotonPx_DeriUnk(1, x, u, derivate_i), GetPhotonPy_DeriUnk(1, x, u, derivate_i), GetPhotonPz_DeriUnk(1, x, u, derivate_i), 0};
            Double_t    ret = 0;
            for(int i=0; i<3; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            for(int i=4; i<7; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            return ret;
        }
    case 1:
        {
            Double_t    h[9] = {GetPhotonPx(2, x, u), GetPhotonPy(2, x, u), GetPhotonPz(2, x, u), GetPhotonE(2, x),
                                GetPhotonPx(3, x, u), GetPhotonPy(3, x, u), GetPhotonPz(3, x, u), GetPhotonE(3, x), MASS_PI0};
            Double_t    dh[8] = {GetPhotonPx_DeriUnk(2, x, u, derivate_i), GetPhotonPy_DeriUnk(2, x, u, derivate_i), GetPhotonPz_DeriUnk(2, x, u, derivate_i), 0,
                                 GetPhotonPx_DeriUnk(3, x, u, derivate_i), GetPhotonPy_DeriUnk(3, x, u, derivate_i), GetPhotonPz_DeriUnk(3, x, u, derivate_i), 0};
            Double_t    ret = 0;
            for(int i=0; i<3; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            for(int i=4; i<7; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            return ret;
        }
    case 2:
        {
            Double_t    h[9] = {GetPhotonPx(4, x, u), GetPhotonPy(4, x, u), GetPhotonPz(4, x, u), GetPhotonE(4, x),
                                GetPhotonPx(5, x, u), GetPhotonPy(5, x, u), GetPhotonPz(5, x, u), GetPhotonE(5, x), MASS_PI0};
            Double_t    dh[8] = {GetPhotonPx_DeriUnk(4, x, u, derivate_i), GetPhotonPy_DeriUnk(4, x, u, derivate_i), GetPhotonPz_DeriUnk(4, x, u, derivate_i), 0,
                                 GetPhotonPx_DeriUnk(5, x, u, derivate_i), GetPhotonPy_DeriUnk(5, x, u, derivate_i), GetPhotonPz_DeriUnk(5, x, u, derivate_i), 0};
            Double_t    ret = 0;
            for(int i=0; i<3; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            for(int i=4; i<7; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            return ret;
        }
    case 3:
        {
            Double_t    h[29] = {GetReactionPx(x, u),  GetReactionPy(x, u),  GetReactionPz(x, u),  GetReactionE(x),
                                 GetPhotonPx(0, x, u), GetPhotonPy(0, x, u), GetPhotonPz(0, x, u), GetPhotonE(0, x),
                                 GetPhotonPx(1, x, u), GetPhotonPy(1, x, u), GetPhotonPz(1, x, u), GetPhotonE(1, x),
                                 GetPhotonPx(2, x, u), GetPhotonPy(2, x, u), GetPhotonPz(2, x, u), GetPhotonE(2, x),
                                 GetPhotonPx(3, x, u), GetPhotonPy(3, x, u), GetPhotonPz(3, x, u), GetPhotonE(3, x),
                                 GetPhotonPx(4, x, u), GetPhotonPy(4, x, u), GetPhotonPz(4, x, u), GetPhotonE(4, x),
                                 GetPhotonPx(5, x, u), GetPhotonPy(5, x, u), GetPhotonPz(5, x, u), GetPhotonE(5, x), MASS_PROTON};
            Double_t    dh[28] = {GetReactionPx_DeriUnk(x, u, derivate_i), GetReactionPy_DeriUnk(x, u, derivate_i), GetReactionPz_DeriUnk(x, u, derivate_i), 0,
                                  GetPhotonPx_DeriUnk(0, x, u, derivate_i), GetPhotonPy_DeriUnk(0, x, u, derivate_i), GetPhotonPz_DeriUnk(0, x, u, derivate_i), 0,
                                  GetPhotonPx_DeriUnk(1, x, u, derivate_i), GetPhotonPy_DeriUnk(1, x, u, derivate_i), GetPhotonPz_DeriUnk(1, x, u, derivate_i), 0,
                                  GetPhotonPx_DeriUnk(2, x, u, derivate_i), GetPhotonPy_DeriUnk(2, x, u, derivate_i), GetPhotonPz_DeriUnk(2, x, u, derivate_i), 0,
                                  GetPhotonPx_DeriUnk(3, x, u, derivate_i), GetPhotonPy_DeriUnk(3, x, u, derivate_i), GetPhotonPz_DeriUnk(3, x, u, derivate_i), 0,
                                  GetPhotonPx_DeriUnk(4, x, u, derivate_i), GetPhotonPy_DeriUnk(4, x, u, derivate_i), GetPhotonPz_DeriUnk(4, x, u, derivate_i), 0,
                                  GetPhotonPx_DeriUnk(5, x, u, derivate_i), GetPhotonPy_DeriUnk(5, x, u, derivate_i), GetPhotonPz_DeriUnk(5, x, u, derivate_i), 0};
            Double_t    ret = 0;
            for(int i=0; i<3; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            for(int i=4; i<7; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            for(int i=8; i<11; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            for(int i=12; i<15; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            for(int i=16; i<19; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            for(int i=20; i<23; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            for(int i=24; i<27; i++)
                ret += GKinFitterConstrainInvMass().DoDerivative(h, i) * dh[i];
            return ret;
        }
    /*case 4:
        {
            Double_t    h[5] = {GetReactionPx(x, u),  GetReactionPy(x, u),  GetReactionPz(x, u),  GetReactionE(x), MASS_PROTON};
            Double_t    dh[3] = {GetReactionPx_DeriUnk(x, u, derivate_i), GetReactionPy_DeriUnk(x, u, derivate_i), GetReactionPz_DeriUnk(x, u, derivate_i)};
            Double_t    ret = 0;
            for(int i=0; i<3; i++)
            {
                ret += GKinFitterConstrainBeam().DoDerivative(h, i) * dh[i];
                std::cout << "gfkf " << ret << std::endl;
            }
            return ret;
        }    //}
//}



void        GKinFitter::Get_g(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret)
{
    for(int i=0; i<nCon; i++)
        ret[i][0] = GetConstraint(i, x, u);
}
void        GKinFitter::Get_g_DeriPar(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret)
{
    for(int i=0; i<nCon; i++)
    {
        for(int p=0; p<nPar; p++)
            ret[i][p] = GetConstraint_DeriPar(i, x, u, p);
    }
}
void        GKinFitter::Get_g_DeriUnk(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret)
{
    for(int i=0; i<nCon; i++)
    {
        for(int p=0; p<nUnk; p++)
            ret[i][p] = GetConstraint_DeriUnk(i, x, u, p);
    }
}

void        GKinFitter::Get_r(const TMatrixD& x, const TMatrixD& u, const TMatrixD& g_DeriPar, const TMatrixD& g, TMatrixD& ret)
{
    TMatrixD    dpar(par0);
                dpar -= x;
    ret = g + (g_DeriPar * dpar);
}

void        GKinFitter::Get_S(const TMatrixD& g_Derivated_Par, TMatrixD& ret)
{
    TMatrixD    g_Derivated_Par_T(g_Derivated_Par);
    g_Derivated_Par_T.T();

    ret = g_Derivated_Par * V0 * g_Derivated_Par_T;
}

Double_t    GKinFitter::Calc_Chi2(const TMatrixD& x, const TMatrixD& g)
{
    TMatrixD    dpar(par0-x);
    TMatrixD    dpar_T(dpar);
    dpar_T.T();
    TMatrixD    lambda_T(lambda);
    lambda_T.T();
    TMatrixD    help0(dpar_T * V0 * dpar);
    TMatrixD    help1(lambda_T * g);
    return  help0[0][0] + (2 * help1[0][0]);
}


/*
TVector3        GKinFitter::GetMainVertex_Derivated_Unk(const Int_t index_m, const TMatrixD &u)
{
    switch(index_m)
    {
    case 0:
        return TVector3(1, 0, 0);
    case 1:
        return TVector3(0, 1, 0);
    case 2:
        return TVector3(0, 0, 1);
    default:
        return TVector3(0, 0, 0);
    }
}

TVector3        GKinFitter::GetVertex_Derivated_Unk(const Int_t index, const Int_t index_m, const TMatrixD &u)
{
    switch(index)
    {
    case 0:
        switch(index_m)
        {
        case 3:
            return TVector3(1, 0, 0);
        case 4:
            return TVector3(0, 1, 0);
        case 5:
            return TVector3(0, 0, 1);
        default:
            return TVector3(0, 0, 0);
        }
    case 1:
        switch(index_m)
        {
        case 6:
            return TVector3(1, 0, 0);
        case 7:
            return TVector3(0, 1, 0);
        case 8:
            return TVector3(0, 0, 1);
        default:
            return TVector3(0, 0, 0);
        }
    case 2:
        switch(index_m)
        {
        case 9:
            return TVector3(1, 0, 0);
        case 10:
            return TVector3(0, 1, 0);
        case 11:
            return TVector3(0, 0, 1);
        default:
            return TVector3(0, 0, 0);
        }
    default:
        return TVector3(0, 0, 0);
    }
}

TLorentzVector  GKinFitter::GetBeam(const TMatrixD x, const TMatrixD u)
{
    TVector3    beamPos(u[12][0], u[13][0], -GKinFitter_CBRadius);
    beamPos     = GetMainVertex(u)-beamPos;
    beamPos     = x[0][0] * beamPos * (1/beamPos.Mag());
    return TLorentzVector(beamPos, x[0][0]);
}

TLorentzVector  GKinFitter::GetBeam_Derivated_Par(const Int_t index_n, const TMatrixD x, const TMatrixD u)
{
    if(index_n==0 || index_n==1 || index_n==2)
    {
        TVector3    beamPos(u[12][0], u[13][0], -GKinFitter_CBRadius);
        beamPos     = GetMainVertex(u)-beamPos;
        beamPos     = beamPos * (1/beamPos.Mag());
        return TLorentzVector(beamPos, 1);
    }
    else
        return TLorentzVector(0, 0, 0, 0);
}

TLorentzVector  GKinFitter::GetBeam_Derivated_Unk(const Int_t index_m, const TMatrixD x, const TMatrixD u)
{
    if(index_n==12)
    {
        TVector3    beamPos(u[12][0], u[13][0], -GKinFitter_CBRadius);
        beamPos     = GetMainVertex(u)-beamPos;
        beamPos     = beamPos * (1/beamPos.Mag());
        return TLorentzVector(beamPos, 1);
    }
    else
        return TLorentzVector(0, 0, 0, 0);
}

TLorentzVector  GKinFitter::GetParticle(const Int_t index, const TMatrixD x, const TMatrixD u)
{
    TVector3    beamPos(x[(index*4)+1][0], x[(index*4)+2][0], x[(index*4)+3][0]);
    switch(index)
    {
    case 0:
    case 1:
        beamPos     = GetVertex(0, u)-beamPos;
        break;
    case 2:
    case 3:
        beamPos     = GetVertex(1, u)-beamPos;
        break;
    case 4:
    case 5:
        beamPos     = GetVertex(2, u)-beamPos;
        break;
    }
    beamPos     = x[(index*4)+4][0] * beamPos * (1/beamPos.Mag());
    return TLorentzVector(beamPos, x[(index*4)+4][0]);
}

Double_t    GKinFitter::Get_g(const Int_t index_k, const TMatrixD x, const TMatrixD u)
{
    switch(index_k)
    {
    case 0:                     //inv. Mass eta
        return (GetParticle(0, x, u)+GetParticle(1, x, u)).M2()-(MASS_ETA*MASS_ETA);
    case 1:                     //inv. Mass pi0 a
        return (GetParticle(2, x, u)+GetParticle(3, x, u)).M2()-(MASS_PI0*MASS_PI0);
    case 2:                     //inv. Mass pi0 b
        return (GetParticle(4, x, u)+GetParticle(5, x, u)).M2()-(MASS_PI0*MASS_PI0);
    case 3:                     //mis. Mass proton
        return (GetBeam(x, u)-GetParticle(0, x, u)-GetParticle(1, x, u)-GetParticle(2, x, u)-GetParticle(3, x, u)-GetParticle(4, x, u)-GetParticle(5, x, u)).M2()-(MASS_PROTON*MASS_PROTON);
    }
}

TMatrixD        GKinFitter::Get_g_Derivated_Par(const Int_t index_k, const TMatrixD x, const TMatrixD u)
{
    TMatrixD    ret(nPar, 1);

    switch(index_k)
    {
    case 0:                     //inv. Mass eta
        {
            TLorentzVector vec0(GetParticle(0, x, u));
            TLorentzVector vec1(GetParticle(1, x, u));
            ret[0][0]   = 0;
            ret[1][0]   = -2*(vec0.Px()+vec1.Px());
            ret[2][0]   = -2*(vec0.Py()+vec1.Py());
            ret[3][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[4][0]   = 2*(vec0.E()+vec1.E());
            ret[5][0]   = -2*(vec0.Px()+vec1.Px());
            ret[6][0]   = -2*(vec0.Py()+vec1.Py());
            ret[7][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[8][0]   = 2*(vec0.E()+vec1.E());
            for(int i=9; i<nPar; i++)
                ret[i][0]   = 0;
        }
        return ret;
    case 1:                     //inv. Mass pi0 a
        {
            TLorentzVector vec0(GetParticle(2, x, u));
            TLorentzVector vec1(GetParticle(3, x, u));
            for(int i=0; i<9; i++)
                ret[i][0]   = 0;
            ret[9][0]   = -2*(vec0.Px()+vec1.Px());
            ret[10][0]   = -2*(vec0.Py()+vec1.Py());
            ret[11][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[12][0]   = 2*(vec0.E()+vec1.E());
            ret[13][0]   = -2*(vec0.Px()+vec1.Px());
            ret[14][0]   = -2*(vec0.Py()+vec1.Py());
            ret[15][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[16][0]   = 2*(vec0.E()+vec1.E());
            for(int i=17; i<nPar; i++)
                ret[i][0]   = 0;
        }
        return ret;
    case 2:                     //inv. Mass pi0 b
        {
            TLorentzVector vec0(GetParticle(4, x, u));
            TLorentzVector vec1(GetParticle(5, x, u));
            for(int i=0; i<17; i++)
                ret[i][0]   = 0;
            ret[17][0]   = -2*(vec0.Px()+vec1.Px());
            ret[18][0]   = -2*(vec0.Py()+vec1.Py());
            ret[19][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[20][0]   = 2*(vec0.E()+vec1.E());
            ret[21][0]   = -2*(vec0.Px()+vec1.Px());
            ret[22][0]   = -2*(vec0.Py()+vec1.Py());
            ret[23][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[24][0]   = 2*(vec0.E()+vec1.E());
        }
        return ret;
    case 3:                     //mis. Mass proton
    {
        TLorentzVector beam(GetBeam(x, u));
        TLorentzVector vec0(GetParticle(0, x, u));
        TLorentzVector vec1(GetParticle(1, x, u));
        TLorentzVector vec2(GetParticle(2, x, u));
        TLorentzVector vec3(GetParticle(3, x, u));
        TLorentzVector vec4(GetParticle(4, x, u));
        TLorentzVector vec5(GetParticle(5, x, u));
        ret[0][0]   = 2*(beam.E()-vec0.E()-vec1.E()-vec2.E()-vec3.E()-vec4.E()-vec5.E())-2*(beam.Pz()-vec0.Pz()-vec1.Pz()-vec2.Pz()-vec3.Pz()-vec4.Pz()-vec5.Pz());
        ret[1][0]   = -2*(beam.Px()-vec0.Px()-vec1.Px()-vec2.Px()-vec3.Px()-vec4.Px()-vec5.Px());
        ret[2][0]   = -2*(beam.Py()-vec0.Py()-vec1.Py()-vec2.Py()-vec3.Py()-vec4.Py()-vec5.Py());
        ret[3][0]   = -2*(beam.Pz()-vec0.Pz()-vec1.Pz()-vec2.Pz()-vec3.Pz()-vec4.Pz()-vec5.Pz());
        ret[4][0]   = 2*(vec0.E()+vec1.E());
        ret[5][0]   = -2*(vec0.Px()+vec1.Px());
        ret[6][0]   = -2*(vec0.Py()+vec1.Py());
        ret[7][0]   = -2*(vec0.Pz()+vec1.Pz());
        ret[8][0]   = 2*(vec0.E()+vec1.E());
    }
    return ret;
    }
}

TMatrixD        GKinFitter::Get_g_Derivated_Unk(const Int_t index_k, const TMatrixD x, const TMatrixD u)
{

}*/


/*void GKinFitter::Constraints()
{
    /*
  //Add up the particle 4 vectors
    TLorentzVector eta(fmAlpha0[16][0], fmAlpha0[17][0], fmAlpha0[18][0], fmAlpha0[19][0]);
    eta     += TLorentzVector(fmAlpha0[20][0], fmAlpha0[21][0], fmAlpha0[22][0], fmAlpha0[23][0]);
    TLorentzVector pi0a(fmAlpha0[24][0], fmAlpha0[25][0], fmAlpha0[26][0], fmAlpha0[27][0]);
    pi0a    += TLorentzVector(fmAlpha0[28][0], fmAlpha0[29][0], fmAlpha0[30][0], fmAlpha0[31][0]);
    TLorentzVector pi0b(fmAlpha0[32][0], fmAlpha0[33][0], fmAlpha0[34][0], fmAlpha0[35][0]);
    pi0b    += TLorentzVector(fmAlpha0[36][0], fmAlpha0[37][0], fmAlpha0[38][0], fmAlpha0[39][0]);
    TLorentzVector all(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]);
    all     -= eta;
    all     -= pi0a;
    all     -= pi0b;*/
    /*all.Print();
    (all-TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0])).Print();
    TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]).Print();
    TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]).Print();
    std::cout << all.M() << "   " << (all-TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0])).M() << "   " << TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]).M() << "   " << std::endl;
    std::cout << eta.M() << "   " << pi0a.M() << "   " << pi0b.M() << "   " << std::endl;
    std::cout << TLorentzVector(fmAlpha0[16][0], fmAlpha0[17][0], fmAlpha0[18][0], fmAlpha0[19][0]).M() << "   " << TLorentzVector(fmAlpha0[20][0], fmAlpha0[21][0], fmAlpha0[22][0], fmAlpha0[23][0]).M() << "   " << TLorentzVector(fmAlpha0[24][0], fmAlpha0[25][0], fmAlpha0[26][0], fmAlpha0[27][0]).M() << "   " << std::endl;
*//*
  //d matrix (evaluate constraint eqn.)
    fmd[0][0]   =   eta.M2() - (MASS_ETA*MASS_ETA);
    fmd[1][0]   =   pi0a.M2() - (MASS_PI0*MASS_PI0);
    fmd[2][0]   =   pi0b.M2() - (MASS_PI0*MASS_PI0);
    fmd[3][0]   =   all.M2() - (MASS_PROTON*MASS_PROTON);

  //D matrix (derivitives of constraint eqn)
    for(int i=0; i<4; i++)
    {
        for(int j=0; j<40; j++)
            fmD[i][j] = 0;
    }

    //[Cons Number][Var Number]
    for(int i=0; i<2; i++)
    {
        fmD[0][4*i+16] =-2*eta.X();
        fmD[0][4*i+17] =-2*eta.Y();
        fmD[0][4*i+18] =-2*eta.Z();
        fmD[0][4*i+19] = 2*eta.T();

        fmD[1][4*i+24] =-2*pi0a.X();
        fmD[1][4*i+25] =-2*pi0a.Y();
        fmD[1][4*i+26] =-2*pi0a.Z();
        fmD[1][4*i+27] = 2*pi0a.T();

        fmD[2][4*i+32] =-2*pi0b.X();
        fmD[2][4*i+33] =-2*pi0b.Y();
        fmD[2][4*i+34] =-2*pi0b.Z();
        fmD[2][4*i+35] = 2*pi0b.T();
    }

    fmD[3][12] =-2*all.X();
    fmD[3][13] =-2*all.Y();
    fmD[3][14] =-2*all.Z();
    fmD[3][15] = 2*all.T();
    for(int i=0; i<6; i++)
    {
        fmD[3][4*i+16] =-2*all.X();
        fmD[3][4*i+17] =-2*all.Y();
        fmD[3][4*i+18] =-2*all.Z();
        fmD[3][4*i+19] = 2*all.T();
    }*/
//}

/*TLorentzVector  GKinFitter::GetEta()
{
    TLorentzVector ret(fmAlpha[16][0], fmAlpha[17][0], fmAlpha[18][0], fmAlpha[19][0]);
        ret +=  TLorentzVector(fmAlpha[20][0], fmAlpha[21][0], fmAlpha[22][0], fmAlpha[23][0]);
    return ret;
}*/

/*TLorentzVector  GKinFitter::GetEtap()
{
    //TLorentzVector(fmAlpha[12][0], fmAlpha[13][0], fmAlpha[14][0], fmAlpha[15][0]).Print();
    //TLorentzVector(fmAlpha0[21][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]).Print();

    TLorentzVector ret(0.0, 0.0, 0.0, 0.0);
    for(int i=0; i<6; i++)
    {
        ret +=  TLorentzVector(fmAlpha[4*i+16][0], fmAlpha[4*i+17][0], fmAlpha[4*i+18][0], fmAlpha[4*i+19][0]);
        //TLorentzVector(fmAlpha[4*i+16][0], fmAlpha[4*i+17][0], fmAlpha[4*i+18][0], fmAlpha[4*i+19][0]).Print();
        //TLorentzVector(fmAlpha0[4*i+16][0], fmAlpha0[4*i+17][0], fmAlpha0[4*i+18][0], fmAlpha0[4*i+19][0]).Print();
    }
    //ret.Print();
    //std::cout << std::endl;
    return ret;
}

TLorentzVector  GKinFitter::GetPi0a()
{
    TLorentzVector ret(fmAlpha[24][0], fmAlpha[25][0], fmAlpha[26][0], fmAlpha[27][0]);
        ret +=  TLorentzVector(fmAlpha[28][0], fmAlpha[29][0], fmAlpha[30][0], fmAlpha[31][0]);
    return ret;
}

TLorentzVector  GKinFitter::GetPi0b()
{
    TLorentzVector ret(fmAlpha[32][0], fmAlpha[33][0], fmAlpha[34][0], fmAlpha[35][0]);
        ret +=  TLorentzVector(fmAlpha[36][0], fmAlpha[37][0], fmAlpha[38][0], fmAlpha[39][0]);
    return ret;
}*/


/*void    GKinFitter::Set(const GKinFitterPolarToCartesian& v)
{
    par0    = v.GetParametersH();
    V0      = v.GetCovarianceH();

    for(int i=0; i<nUnk; i++)
        unk[i]  = 0;
}

Bool_t  GKinFitter::SolveStep(const TMatrixD& x, const TMatrixD& u, TMatrixD& new_x, TMatrixD& new_u, Double_t& chiSq)
{
    TMatrixD    g_DeriPar(nCon, nPar);
    Get_g_DeriPar(x, u, g_DeriPar);
    //g_DeriPar.Print();
    TMatrixD    g_DeriPar_T(g_DeriPar);
    g_DeriPar_T.T();
    //g_Derivated_Par_T.Print();
    TMatrixD    S_Inverted(nCon, nCon);
    Get_S(g_DeriPar, S_Inverted);
    //S_Inverted.Print();
    Double_t    determinant;
    S_Inverted.Invert(&determinant);
    if(determinant==0)
    {
        std::cout << "Can not invert S in KinFitter::SolveStep" << std::endl;
        return kFALSE;
    }
    //S_Inverted.Print();

    TMatrixD    g_DeriUnk(nCon, nUnk);
    Get_g_DeriUnk(x, u, g_DeriUnk);
    TMatrixD    g_DeriUnk_T(g_DeriUnk);
    g_DeriUnk_T.T();

    TMatrixD    help(g_DeriUnk_T * S_Inverted * g_DeriUnk);
    //g_DeriUnk.Print();
    help.Invert(&determinant);
    //std::cout << determinant << std::endl;
    //g_DeriUnk.Print();
    if(determinant==0)
    {
        std::cout << "Can not invert help in KinFitter::SolveStep" << std::endl;
        return kFALSE;
    }

    TMatrixD    g(nCon, 1);
    TMatrixD    r(nCon, 1);
    Get_g(x, u, g);
    Get_r(x, u, g_DeriPar, g, r);

    new_u   = u - (help * g_DeriUnk_T * S_Inverted * r);
    lambda  = S_Inverted * (r + (g_DeriUnk * (new_u - u)));
    new_x   = par0 - (V0 * g_DeriPar_T * lambda);

    chiSq   = Calc_Chi2(new_x, g);
    std::cout << chiSq << std::endl;

    std::cout << "Beam: ";
    std::cout << GetReactionPx(new_x, new_u) << "   ";
    std::cout << GetReactionPy(new_x, new_u) << "   ";
    std::cout << GetReactionPz(new_x, new_u) << "   ";
    std::cout << GetReactionE(new_x) << "   " << std::endl;
    for(int i=0; i<6; i++)
    {
        std::cout << "Photon: ";
        std::cout << GetPhotonPx(i, new_x, new_u) << "   ";
        std::cout << GetPhotonPy(i, new_x, new_u) << "   ";
        std::cout << GetPhotonPz(i, new_x, new_u) << "   ";
        std::cout << GetPhotonE(i, new_x) << "   " << std::endl;
    }
    for(int i=0; i<nUnk; i++)
    {
        std::cout << "Vertex " << i << ": ";
        std::cout << new_u[i][0] << std::endl;
    }
    for(int i=0; i<nCon; i++)
    {
        std::cout << "Constraints " << i << ": ";
        std::cout << GetConstraint(i, new_x, new_u) << std::endl;
    }

    return kTRUE;
}

Bool_t  GKinFitter::Solve()
{
    std::cout << "Beam: ";
    std::cout << GetInitialReactionPx() << "   ";
    std::cout << GetInitialReactionPy() << "   ";
    std::cout << GetInitialReactionPz() << "   ";
    std::cout << GetInitialReactionE() << "   " << std::endl;
    for(int i=0; i<6; i++)
    {
        std::cout << "Photon: ";
        std::cout << GetInitialPhotonPx(i) << "   ";
        std::cout << GetInitialPhotonPy(i) << "   ";
        std::cout << GetInitialPhotonPz(i) << "   ";
        std::cout << GetInitialPhotonE(i) << "   " << std::endl;
    }
    for(int i=0; i<nCon; i++)
    {
        std::cout << "Constraints " << i << ": ";
        std::cout << GetConstraint(i, par0, unk0) << std::endl;
    }
    if(SolveStep(par0, unk0, par, unk, chi2)==kFALSE)
    {
        std::cout << "first SolveStep not working in KinFitter::Solve" << std::endl;
        return kFALSE;
    }
    std::cout << "ChiSq: "<< chi2 << std::endl;

    TMatrixD    new_x(nPar, 1);
    TMatrixD    new_u(nUnk, 1);
    Double_t    newChi2;

    for(int i=0; i<20; i++)
    {
        if(SolveStep(par, unk, new_x, new_u, newChi2)==kFALSE)
        {
            std::cout << i << " SolveStep not working in KinFitter::Solve" << std::endl;
            return kFALSE;
        }
        std::cout << "ChiSq: "<< newChi2 << std::endl;
        if(newChi2<chi2)
        {
            chi2    = newChi2;
            par     = new_x;
            unk     = new_u;
        }
        else
            return kTRUE;
    }

    std::cout << "No solution after 20 steps in KinFitter::Solve" << std::endl;
    return kFALSE;
}*/
