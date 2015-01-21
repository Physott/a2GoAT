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
    para(nPar,1),
    unk0(nUnk,1),
    unkn(nUnk,1),
    lambda(nCon,1),
    V0(nPar, nPar),
    V(nPar, nPar),
    chiSq(0),
    G(nCon, 1),
    GPar(nCon, nPar),
    GParT(nPar, nCon),
    GUnk(nCon, nUnk),
    GUnkT(nUnk, nCon),
    SInv(nCon, nCon),
    U(nUnk, nUnk)
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
        return ((nPart+1)*GKinFitter_ParametersPerParticle);    //+3 for the beam
    case flagRecoilAngles:
        return ((nPart+1)*GKinFitter_ParametersPerParticle) + 2;    //+1 for the beam +2 for the recoil angles
    }
}

void    GKinFitter::FillHists(GHistFit2& hists)
{
    TLorentzVector  eta(GKinFitterGamma().Particle(0, para, unkn));
                    eta     += GKinFitterGamma().Particle(1, para, unkn);
    TLorentzVector  pi0a(GKinFitterGamma().Particle(2, para, unkn));
                    pi0a    += GKinFitterGamma().Particle(3, para, unkn);
    TLorentzVector  pi0b(GKinFitterGamma().Particle(4, para, unkn));
                    pi0b    += GKinFitterGamma().Particle(5, para, unkn);
    TLorentzVector  etap(eta);
                    etap    += pi0a;
                    etap    += pi0b;

    hists.invMass.Fill(etap.M(), countIter);
    hists.invMassSub0.Fill(eta.M(), countIter);
    hists.invMassSub1.Fill(pi0a.M(), countIter);
    hists.invMassSub2.Fill(pi0b.M(), countIter);

    hists.chiSq.Fill(chiSq, countIter);
    hists.confidenceLevel.Fill(ConfidenceLevel(), countIter);
    hists.zVertex.Fill(unkn[0][0], countIter);
    for(int i=0; i<nPart; i++)
    {
        hists.photonsEnergy[i]->Fill(para[(i+1)*GKinFitter_ParametersPerParticle][0], countIter);
        hists.photonsTheta[i]->Fill(para[((i+1)*GKinFitter_ParametersPerParticle)+1][0], countIter);
        hists.photonsPhi[i]->Fill(para[((i+1)*GKinFitter_ParametersPerParticle)+2][0], countIter);
    }
    hists.beamEnergy.Fill(para[0][0], countIter);
    hists.beamTheta.Fill(para[1][0], countIter);
    hists.beamPhi.Fill(para[2][0], countIter);
    int indices[2];
    indices[0]  = 0;
    indices[1]  = 1;
    hists.con[0]->Fill(GKinFitterInvMass().Eval(para, unkn, indices, 2, MASS_ETA), countIter);
    indices[0]  = 2;
    indices[1]  = 3;
    hists.con[1]->Fill(GKinFitterInvMass().Eval(para, unkn, indices, 2, MASS_PI0), countIter);
    indices[0]  = 4;
    indices[1]  = 5;
    hists.con[2]->Fill(GKinFitterInvMass().Eval(para, unkn, indices, 2, MASS_PI0), countIter);
    for(int i=0; i<nPar; i++)
        hists.pulls[i]->Fill(Pull(i), countIter);

    switch(fitType)
    {
    case flagNoRecoil:
        hists.con[3]->Fill(MM(para, unkn, MASS_PROTON), countIter);
        break;
    case flagUnknownRecoil:
        hists.protonEnergy.Fill(unkn[1][0], countIter);
        hists.protonTheta.Fill(unkn[2][0], countIter);
        hists.protonPhi.Fill(unkn[3][0], countIter);
        {
            TLorentzVector  help(I4Vec(para, unkn, TLorentzVector(0.0, 0.0, 0.0, 0.0)));
            hists.con[3]->Fill(help.Px(), countIter);
            hists.con[4]->Fill(help.Py(), countIter);
            hists.con[5]->Fill(help.Pz(), countIter);
            hists.con[6]->Fill(help.E(), countIter);
        }
        break;
    case flagRecoilAngles:
        hists.protonEnergy.Fill(unkn[1][0], countIter);
        hists.protonTheta.Fill(para[(nPart+1)*GKinFitter_ParametersPerParticle][0], countIter);
        hists.protonPhi.Fill(para[((nPart+1)*GKinFitter_ParametersPerParticle)+1][0], countIter);
        hists.con;
        {
            TLorentzVector  help(I4Vec(para, unkn, TLorentzVector(0.0, 0.0, 0.0, 0.0)));
            hists.con[3]->Fill(help.Px(), countIter);
            hists.con[4]->Fill(help.Py(), countIter);
            hists.con[5]->Fill(help.Pz(), countIter);
            hists.con[6]->Fill(help.E(), countIter);
        }
        break;
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

Bool_t          GKinFitter::SolveStep(TMatrixD& par, TMatrixD& unk, const TMatrixD& oldPar, const TMatrixD& oldUnk, GHistFit2& hists)
{
    gDerivatePar(GPar, oldPar, oldUnk);
    GParT.Transpose(GPar);

    //GPar.Print();
    //V0.Print();
    //GParT.Print();
    SInv    = GPar*V0*GParT;
    Double_t    determinant    = SInv.Determinant();
    if(determinant==0)
        return kFALSE;
    SInv.Invert();


    gDerivateUnk(GUnk, oldPar, oldUnk);
    GUnkT.Transpose(GUnk);
    TMatrixD    mHelp(GUnkT*SInv*GUnk);
    U   = mHelp;
    determinant    = U.Determinant();
    if(determinant==0)
        return kFALSE;
    U.Invert();

    g(G, oldPar, oldUnk);
    TMatrixD    R(nCon, 1);
    r(R, oldPar, oldUnk);
    unk  = oldUnk;
    unk -= U*GUnkT*SInv*R;

    TMatrixD    diff(unk);
    diff    -= oldUnk;
    TMatrixD    help(R);
    help    += GUnk*diff;
    lambda   = SInv*help;

    par  = oldPar;
    par -= V0*GParT*lambda;

    if(CalcChiSq(par, unk)==kFALSE)
        return kFALSE;

    FillHists(hists);
    //std::cout << "SolveStep" << std::endl;
    return kTRUE;
}

Bool_t          GKinFitter::Solve(GHistFit2& hists)
{
    std::cout << "Solve" << std::endl;
    //Calc Unk0
    TLorentzVector    tot(-GKinFitterGamma().BeamAndTarget(para, unkn, targetMass));
    unk0[0][0]   = 0;
    switch(fitType)
    {
    case flagNoRecoil:
        break;
    case flagUnknownRecoil:
        for(int i=0; i<nPart; i++)
            tot += GKinFitterGamma().Particle(i, par0, unk0);
        unk0[1][0]   = tot.E();
        unk0[2][0]   = tot.Theta();
        unk0[3][0]   = tot.Phi();
        break;
    case flagRecoilAngles:
        for(int i=0; i<nPart; i++)
            tot += GKinFitterGamma().Particle(i, par0, unk0);
        unk0[1][0]   = tot.E();
        break;
    }

    if(SolveStep(para, unkn, par0, unk0, hists)==kFALSE)
        return kFALSE;
    Double_t    oldChiSq    = chiSq;
    TMatrixD    newPar(para);
    TMatrixD    newUnk(unkn);
    countIter   = 1;
    std::cout << "ChiSq: " << chiSq << std::endl;
    if(SolveStep(newPar, newUnk, para, unkn, hists)==kFALSE)
        return CalcV(para, unkn);
    countIter   = 2;
    std::cout << "ChiSq: " << chiSq << std::endl;
    while(countIter<=GKinFitter_MaxSteps)// && chiSq<oldChiSq)
    {
        std::cout << "ChiSq: " << chiSq << std::endl;
        oldChiSq    = chiSq;
        para        = newPar;
        unkn        = newUnk;
        countIter++;
        if(SolveStep(newPar, newUnk, para, unkn, hists)==kFALSE)
            return CalcV(para, unkn);
    }
    if(countIter==GKinFitter_MaxSteps)
        std::cout << "Fitter not finished. " << GKinFitter_MaxSteps << " steps are not enough." << std::endl;
    return CalcV(para, unkn);
}

Bool_t    GKinFitter::CalcChiSq(const TMatrixD& par, const TMatrixD& unk)
{
    TMatrixD    lambdaT(lambda);
    lambdaT.T();
    TMatrixD    res(lambdaT*G);
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


Bool_t    GKinFitter::CalcV(const TMatrixD& par, const TMatrixD& unk)
{
    TMatrixD    A(GParT*SInv*GPar);
    TMatrixD    B(GParT*SInv*GUnk);
    TMatrixD    BT(nUnk, nPar);
    BT.Transpose(B);

    TMatrixD    help(B*U*BT);
    A   -= help;
    A   *= V0;

    TMatrixD    I(nPar, nPar);
    I.UnitMatrix();
    I   -= A;
    V.Mult(V0, I);

    std::cout << "lambda: " << lambda[0][0] << std::endl;
    std::cout << "lambda: " << lambda[1][0] << std::endl;
    std::cout << "lambda: " << lambda[2][0] << std::endl;
    std::cout << "lambda: " << lambda[3][0] << std::endl;
    std::cout << "lambda: " << lambda[4][0] << std::endl;
    std::cout << "lambda: " << lambda[5][0] << std::endl;
    std::cout << "lambda: " << lambda[6][0] << std::endl;

    return kTRUE;
}
