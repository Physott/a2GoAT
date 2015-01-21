//////////////////////////////////////////////////////////////////
//GKinFitterBase
//////////////////////////////////////////////////////////////////


#include "GKinFitterBase.h"
#include "TDecompLU.h"
#include <iostream>


#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046






GKinFitterBase::GKinFitterBase(const Int_t nParticles, const Int_t nConstraints)    :
    nPart(nParticles),
    nCon(nConstraints),
    nPar((nPart+1)*GKinFitterBase_ParametersPerParticle),
    countPart(0),
    countCon(0),
    countIter(0),
    par0(nPar,1),
    para(nPar,1),
    unk0(0),
    unkn(0),
    lambda(nCon,1),
    V0(nPar, nPar),
    V(nPar, nPar),
    chiSq(0),
    nImConstraint(0),
    isMmConstraint(kFALSE),
    G(nCon, 1),
    GPar(nCon, nPar),
    GUnk(nCon, 1),
    SInv(nCon, nCon),
    U(0)
{
}

GKinFitterBase::~GKinFitterBase()
{

}

void    GKinFitterBase::FillHists(GHistFit2& hists)
{
    /*TLorentzVector  eta(GKinFitterBaseGamma().Particle(0, para, unkn));
                    eta     += GKinFitterBaseGamma().Particle(1, para, unkn);
    TLorentzVector  pi0a(GKinFitterBaseGamma().Particle(2, para, unkn));
                    pi0a    += GKinFitterBaseGamma().Particle(3, para, unkn);
    TLorentzVector  pi0b(GKinFitterBaseGamma().Particle(4, para, unkn));
                    pi0b    += GKinFitterBaseGamma().Particle(5, para, unkn);
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
        hists.photonsEnergy[i]->Fill(para[(i+1)*GKinFitterBase_ParametersPerParticle][0], countIter);
        hists.photonsTheta[i]->Fill(para[((i+1)*GKinFitterBase_ParametersPerParticle)+1][0] * TMath::RadToDeg(), countIter);
        hists.photonsPhi[i]->Fill(para[((i+1)*GKinFitterBase_ParametersPerParticle)+2][0] * TMath::RadToDeg(), countIter);
    }
    hists.beamEnergy.Fill(para[0][0], countIter);
    hists.beamTheta.Fill(para[1][0] * TMath::RadToDeg(), countIter);
    hists.beamPhi.Fill(para[2][0] * TMath::RadToDeg(), countIter);
    int indices[2];
    indices[0]  = 0;
    indices[1]  = 1;
    hists.con[0]->Fill(GKinFitterBaseInvMass().Eval(para, unkn, indices, 2, MASS_ETA), countIter);
    indices[0]  = 2;
    indices[1]  = 3;
    hists.con[1]->Fill(GKinFitterBaseInvMass().Eval(para, unkn, indices, 2, MASS_PI0), countIter);
    indices[0]  = 4;
    indices[1]  = 5;
    hists.con[2]->Fill(GKinFitterBaseInvMass().Eval(para, unkn, indices, 2, MASS_PI0), countIter);
    for(int i=0; i<nPar; i++)
        hists.pulls[i]->Fill(Pull(i), countIter);

   hists.protonEnergy.Fill(unkn[1][0], countIter);
   hists.protonTheta.Fill(para[(nPart+1)*GKinFitterBase_ParametersPerParticle][0] * TMath::RadToDeg(), countIter);
   hists.protonPhi.Fill(para[((nPart+1)*GKinFitterBase_ParametersPerParticle)+1][0] * TMath::RadToDeg(), countIter);

   TLorentzVector  help(I4Vec(para, unkn, TLorentzVector(0.0, 0.0, 0.0, 0.0)));
   hists.con[3]->Fill(help.Px(), countIter);
   hists.con[4]->Fill(help.Py(), countIter);
   hists.con[5]->Fill(help.Pz(), countIter);
   hists.con[6]->Fill(help.E(), countIter);*/
}



//public members

void            GKinFitterBase::AddInvMassConstraint(const int* partList, const int nPartList, const Double_t mass)
{
    imConstraint[nImConstraint].nPartList   = nPartList;
    for(int i=0; i<nPartList; i++)
        imConstraint[nImConstraint].PartList[i] = partList[i];
    imConstraint[nImConstraint].mass            = mass;
    nImConstraint++;
    countCon++;
}

void            GKinFitterBase::AddMisMassConstraint(const Double_t mass)
{
    mmConstraint    = mass;
    isMmConstraint  = kTRUE;
    countCon++;
}

void            GKinFitterBase::AddBeam(const Double_t beamEnergy, const Double_t _targetMass, const Double_t beamEnergyError, const Double_t beamSpotRadius)
{
    targetMass  = _targetMass;
    par0[0][0]  = beamEnergy;
    par0[1][0]  = 0;
    par0[2][0]  = 0;
    V0[0][0]    = beamEnergyError*beamEnergyError;
    V0[1][1]    = beamSpotRadius*beamSpotRadius/100;
    V0[2][2]    = 2*TMath::Pi();
}
void            GKinFitterBase::AddGamma(const Double_t energy, const Double_t theta, const Double_t phi, const Double_t energyError, const Double_t thetaError, const Double_t phiError)
{
    if(countPart==nPart)
    {
        std::cout << "Can not add Particle. Already " << countPart << " there." << std::endl;
        return;
    }
    par0[ (countPart+1)*GKinFitterBase_ParametersPerParticle   ][0]  = energy;
    par0[((countPart+1)*GKinFitterBase_ParametersPerParticle)+1][0]  = theta;
    par0[((countPart+1)*GKinFitterBase_ParametersPerParticle)+2][0]  = phi;
    V0[ (countPart+1)*GKinFitterBase_ParametersPerParticle   ][ (countPart+1)*GKinFitterBase_ParametersPerParticle   ]    = energyError*energyError;
    V0[((countPart+1)*GKinFitterBase_ParametersPerParticle)+1][((countPart+1)*GKinFitterBase_ParametersPerParticle)+1]    = thetaError*thetaError;
    V0[((countPart+1)*GKinFitterBase_ParametersPerParticle)+2][((countPart+1)*GKinFitterBase_ParametersPerParticle)+2]    = phiError*phiError;
    countPart++;
}
void            GKinFitterBase::AddRecoilAngles(const Double_t theta, const Double_t phi, const Double_t thetaError, const Double_t phiError)
{
    par0[ (nPart+1)*GKinFitterBase_ParametersPerParticle   ][0]  = theta;
    par0[((nPart+1)*GKinFitterBase_ParametersPerParticle)+1][0]  = phi;
    V0[ (nPart+1)*GKinFitterBase_ParametersPerParticle   ][ (nPart+1)*GKinFitterBase_ParametersPerParticle   ]    = thetaError*thetaError;
    V0[((nPart+1)*GKinFitterBase_ParametersPerParticle)+1][((nPart+1)*GKinFitterBase_ParametersPerParticle)+1]    = phiError*phiError;
}

Double_t        GKinFitterBase::GetInitialIMConstraint(const int i)    const
{
    if(i>=nImConstraint || i<0)
        return 0.0;

    TLorentzVector  tot(GetInitialPhoton(imConstraint[i].PartList[0]));
    for(int c=1; c<imConstraint[i].nPartList; c++)
        tot +=  GetInitialPhoton(imConstraint[i].PartList[c]);

    //std::cout << tot.M2() - (imConstraint[i].mass*imConstraint[i].mass) << std::endl;
    //std::cout << tot.M() - imConstraint[i].mass << std::endl;
    //std::cout << tot.M() << std::endl;
    //std::cout << imConstraint[i].mass << std::endl;
    return tot.M2() - (imConstraint[i].mass*imConstraint[i].mass);
}
void       GKinFitterBase::GetInitialIMConstraintDerPar(TMatrixD& ret, const int i)
{
    for(int c=0; c<nPar; c++)
        ret[c][0]   = 0;

    if(i>=nImConstraint || i<0)
        return;

    TLorentzVector  tot(GetInitialPhoton(imConstraint[i].PartList[0]));
    for(int c=1; c<imConstraint[i].nPartList; c++)
        tot +=  GetInitialPhoton(imConstraint[i].PartList[c]);

    for(int c=0; c<imConstraint[i].nPartList; c++)
    {
        ret[ (imConstraint[i].PartList[c]+1)*GKinFitterBase_ParametersPerParticle   ][0]   = 2 * tot * GetInitialPhotonDerivateEnergy(imConstraint[i].PartList[c]);
        ret[((imConstraint[i].PartList[c]+1)*GKinFitterBase_ParametersPerParticle)+1][0]   = 2 * tot * GetInitialPhotonDerivateTheta(imConstraint[i].PartList[c]);
        ret[((imConstraint[i].PartList[c]+1)*GKinFitterBase_ParametersPerParticle)+2][0]   = 2 * tot * GetInitialPhotonDerivatePhi(imConstraint[i].PartList[c]);
    }
}
Double_t   GKinFitterBase::GetInitialIMConstraintDerUnk(const int i)   const
{
    if(i>=nImConstraint || i<0)
        return 0.0;

    TLorentzVector  tot(GetInitialPhoton(imConstraint[i].PartList[0]));
    for(int c=1; c<imConstraint[i].nPartList; c++)
        tot +=  GetInitialPhoton(imConstraint[i].PartList[c]);

    Double_t    ret = 0;
    for(int c=0; c<imConstraint[i].nPartList; c++)
        ret += 2* tot * GetInitialPhotonDerivateZVertex(imConstraint[i].PartList[c]);

    return ret;
}
Double_t        GKinFitterBase::GetIMConstraint(const int i)    const
{
    if(i>=nImConstraint || i<0)
        return 0.0;

    TLorentzVector  tot(GetPhoton(imConstraint[i].PartList[0]));
    for(int c=1; c<imConstraint[i].nPartList; c++)
        tot +=  GetPhoton(imConstraint[i].PartList[c]);

    return tot.M2() - (imConstraint[i].mass*imConstraint[i].mass);
}
void       GKinFitterBase::GetIMConstraintDerPar(TMatrixD& ret, const int i)
{
    for(int c=0; c<nPar; c++)
        ret[c][0]   = 0;

    if(i>=nImConstraint || i<0)
        return;

    TLorentzVector  tot(GetPhoton(imConstraint[i].PartList[0]));
    for(int c=1; c<imConstraint[i].nPartList; c++)
        tot +=  GetPhoton(imConstraint[i].PartList[c]);

    for(int c=0; c<imConstraint[i].nPartList; c++)
    {
        ret[ (imConstraint[i].PartList[c]+1)*GKinFitterBase_ParametersPerParticle   ][0]   = 2 * tot * GetPhotonDerivateEnergy(imConstraint[i].PartList[c]);
        ret[((imConstraint[i].PartList[c]+1)*GKinFitterBase_ParametersPerParticle)+1][0]   = 2 * tot * GetPhotonDerivateTheta(imConstraint[i].PartList[c]);
        ret[((imConstraint[i].PartList[c]+1)*GKinFitterBase_ParametersPerParticle)+2][0]   = 2 * tot * GetPhotonDerivatePhi(imConstraint[i].PartList[c]);
    }
}
Double_t   GKinFitterBase::GetIMConstraintDerUnk(const int i)          const
{
    if(i>=nImConstraint || i<0)
        return 0.0;

    TLorentzVector  tot(GetPhoton(imConstraint[i].PartList[0]));
    for(int c=1; c<imConstraint[i].nPartList; c++)
        tot +=  GetPhoton(imConstraint[i].PartList[c]);

    Double_t    ret = 0;
    for(int c=0; c<imConstraint[i].nPartList; c++)
        ret += 2* tot * GetPhotonDerivateZVertex(imConstraint[i].PartList[c]);

    return ret;
}

Double_t        GKinFitterBase::GetInitialMMConstraint()    const
{
    TLorentzVector  tot(0.0, 0.0, 0.0, targetMass);
    tot += GetInitialBeam();

    for(int i=0; i<nPart; i++)
        tot -=  GetInitialPhoton(i);

    //std::cout << tot.M2() - (mmConstraint*mmConstraint) << std::endl;
    //std::cout << tot.M() - mmConstraint << std::endl;
    //std::cout << tot.M() << std::endl;
    //std::cout << mmConstraint << std::endl;
    return tot.M2() - (mmConstraint*mmConstraint);
}
void            GKinFitterBase::GetInitialMMConstraintDerPar(TMatrixD& ret)
{
    TLorentzVector  tot(0.0, 0.0, 0.0, targetMass);
    tot += GetInitialBeam();

    for(int i=0; i<nPart; i++)
        tot -=  GetInitialPhoton(i);

    ret[0][0]   = 2 * tot * GetInitialBeamDerivateEnergy();
    ret[1][0]   = 2 * tot * GetInitialBeamDerivateTheta();
    ret[2][0]   = 2 * tot * GetInitialBeamDerivatePhi();
    for(int c=0; c<nPart; c++)
    {
        ret[ (c+1)*GKinFitterBase_ParametersPerParticle   ][0]   = 2 * tot * GetInitialPhotonDerivateEnergy(c);
        ret[((c+1)*GKinFitterBase_ParametersPerParticle)+1][0]   = 2 * tot * GetInitialPhotonDerivateTheta(c);
        ret[((c+1)*GKinFitterBase_ParametersPerParticle)+2][0]   = 2 * tot * GetInitialPhotonDerivatePhi(c);
    }
}
Double_t        GKinFitterBase::GetInitialMMConstraintDerUnk()   const
{
    TLorentzVector  tot(0.0, 0.0, 0.0, targetMass);
    tot += GetInitialBeam();

    for(int i=0; i<nPart; i++)
        tot -=  GetInitialPhoton(i);

    Double_t    ret = 0;
    ret += 2 * tot * GetInitialBeamDerivateZVertex();
    for(int c=0; c<nPart; c++)
        ret += 2 * tot * GetInitialPhotonDerivateZVertex(c);

    return ret;
}
Double_t        GKinFitterBase::GetMMConstraint()    const
{
    TLorentzVector  tot(0.0, 0.0, 0.0, targetMass);
    tot += GetBeam();

    for(int i=0; i<nPart; i++)
        tot -=  GetPhoton(i);

    return tot.M2() - (mmConstraint*mmConstraint);
}
void            GKinFitterBase::GetMMConstraintDerPar(TMatrixD& ret)
{
    TLorentzVector  tot(0.0, 0.0, 0.0, targetMass);
    tot += GetBeam();

    for(int i=0; i<nPart; i++)
        tot -=  GetPhoton(i);

    ret[0][0]   = 2 * tot * GetBeamDerivateEnergy();
    ret[1][0]   = 2 * tot * GetBeamDerivateTheta();
    ret[2][0]   = 2 * tot * GetBeamDerivatePhi();
    for(int c=0; c<nPart; c++)
    {
        ret[ (c+1)*GKinFitterBase_ParametersPerParticle   ][0]   = 2 * tot * GetPhotonDerivateEnergy(c);
        ret[((c+1)*GKinFitterBase_ParametersPerParticle)+1][0]   = 2 * tot * GetPhotonDerivateTheta(c);
        ret[((c+1)*GKinFitterBase_ParametersPerParticle)+2][0]   = 2 * tot * GetPhotonDerivatePhi(c);
    }
}
Double_t        GKinFitterBase::GetMMConstraintDerUnk()   const
{
    TLorentzVector  tot(0.0, 0.0, 0.0, targetMass);
    tot += GetBeam();

    for(int i=0; i<nPart; i++)
        tot -=  GetPhoton(i);

    Double_t    ret = 0;
    ret += 2 * tot * GetBeamDerivateZVertex();
    for(int c=0; c<nPart; c++)
        ret += 2 * tot * GetPhotonDerivateZVertex(c);

    return ret;
}

void    GKinFitterBase::GetInitialG(TMatrixD& ret)
{
    int count   = 0;
    for(int i=0; i<nImConstraint; i++)
    {
        ret[count][0]   = GetInitialIMConstraint(i);
        count++;
    }
    if(isMmConstraint==kTRUE)
        ret[count][0]   = GetInitialMMConstraint();
}
void    GKinFitterBase::GetInitialGDerPar(TMatrixD& ret)
{
    TMatrixD    help(nPar, 1);
    int count   = 0;
    for(int i=0; i<nImConstraint; i++)
    {
        GetInitialIMConstraintDerPar(help, i);
        for(int m=0; m<nPar; m++)
            ret[count][m]   = help[m][0];
        count++;
    }
    if(isMmConstraint==kTRUE)
    {
        GetInitialMMConstraintDerPar(help);
        for(int m=0; m<nPar; m++)
            ret[count][m]   = help[m][0];
    }
}
void    GKinFitterBase::GetInitialGDerUnk(TMatrixD& ret)
{
    int count   = 0;
    for(int i=0; i<nImConstraint; i++)
    {
        ret[count][0]   = GetInitialIMConstraintDerUnk(i);
        count++;
    }
    if(isMmConstraint==kTRUE)
        ret[count][0]   = GetInitialMMConstraintDerUnk();
}
void    GKinFitterBase::GetG(TMatrixD& ret)
{
    int count   = 0;
    for(int i=0; i<nImConstraint; i++)
    {
        ret[count][0]   = GetIMConstraint(i);
        count++;
    }
    if(isMmConstraint==kTRUE)
        ret[count][0]   = GetMMConstraint();
}
void    GKinFitterBase::GetGDerPar(TMatrixD& ret)
{
    TMatrixD    help(nPar, 1);
    int count   = 0;
    for(int i=0; i<nImConstraint; i++)
    {
        GetIMConstraintDerPar(help, i);
        for(int m=0; m<nPar; m++)
            ret[count][m]   = help[m][0];
        count++;
    }
    if(isMmConstraint==kTRUE)
    {
        GetMMConstraintDerPar(help);
        for(int m=0; m<nPar; m++)
            ret[count][m]   = help[m][0];
    }
}
void    GKinFitterBase::GetGDerUnk(TMatrixD& ret)
{
    int count   = 0;
    for(int i=0; i<nImConstraint; i++)
    {
        ret[count][0]   = GetIMConstraintDerUnk(i);
        count++;
    }
    if(isMmConstraint==kTRUE)
        ret[count][0]   = GetMMConstraintDerUnk();
}

void    GKinFitterBase::GetInitialR(TMatrixD& ret)
{
    //TMatrixD    diff(par0)
}
void    GKinFitterBase::GetR(TMatrixD& ret)
{

}










void            GKinFitterBase::Print(const char* option)
{
    std::cout << "GKinFitterBase: nPart=" << nPart << ", nCon=" << nCon;

    TString str(option);
    str.ToLower();
    if(strcmp(str.Data(), "input")==0)
    {
        if(countPart<nPart)
            std::cout << "Input not completed. Only " << countPart << " Particles added yet." << std::endl;
        std::cout << "Input Parameters GKinFitterBase::par0     beam, particles, recoil angles(optional)" << std::endl;
        par0.Print();
        std::cout << "Input Covariance Matrix GKinFitterBase::V0     beam, particles, recoil angles(optional)" << std::endl;
        V0.Print();
    }
}

Bool_t          GKinFitterBase::SolveStep(TMatrixD& newPar, Double_t& newUnk,GHistFit2& hists)
{
    GetGDerPar(GPar);
    TMatrixD    help1(nCon, nPar);
    help1.MultT(V0, GPar);
    SInv.Mult(GPar, help1);
    Double_t    det = SInv.Determinant();
    if(det==0)
    {
        std::cout << "Can not invert S." << std::endl;
        return kFALSE;
    }
    SInv.InvertFast();

    GetGDerUnk(GUnk);
    TMatrixD    help2(nCon, 1);
    help2.Mult(SInv, GUnk);
    TMatrixD    help3(1, 1);
    help3.TMult(GUnk, help2);
    U   = 1/help3[0][0];



    return kTRUE;
}

Bool_t          GKinFitterBase::Solve(GHistFit2& hists)
{
    //Check Parameters
    if(countPart!=nPart)
    {
        std::cout << "Only " << countPart << " Particles filled. " << nPart << " needed." << std::endl;
        return kFALSE;
    }
    if(countCon!=nCon)
    {
        std::cout << "Only " << countCon << " Constraints filled. " << nCon << " needed." << std::endl;
        return kFALSE;
    }

    //Calc Unk0
    unk0    = 0;
    unkn    = 0;

    TMatrixD    newPar(nPar, 1);
    Double_t    newUnk;

    if(SolveStep(newPar, newUnk, hists)==kFALSE)
        return kFALSE;
    /*Double_t    oldChiSq    = chiSq;
    TMatrixD    newPar(para);
    TMatrixD    newUnk(unkn);
    countIter   = 1;
    //std::cout << "ChiSq: " << chiSq << std::endl;
    if(SolveStep(newPar, newUnk, para, unkn, hists)==kFALSE)
        return CalcV(para, unkn);
    //std::cout << "ChiSq: " << chiSq << std::endl;
    while(countIter<=5)// && chiSq<oldChiSq)
    {
        //std::cout << "ChiSq: " << chiSq << std::endl;
        oldChiSq    = chiSq;
        para        = newPar;
        unkn        = newUnk;
        countIter++;
        if(SolveStep(newPar, newUnk, para, unkn, hists)==kFALSE)
            return CalcV(para, unkn);
    }
    if(countIter==GKinFitterBase_MaxSteps)
        std::cout << "Fitter not finished. " << GKinFitterBase_MaxSteps << " steps are not enough." << std::endl;
    return CalcV(para, unkn);*/
}

Bool_t    GKinFitterBase::CalcChiSq(const TMatrixD& par, const TMatrixD& unk)
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


Bool_t    GKinFitterBase::CalcV(const TMatrixD& par, const TMatrixD& unk)
{
    /*TMatrixD    A(GParT*SInv*GPar);
    TMatrixD    B(GParT*SInv*GUnk);
    TMatrixD    BT(1, nPar);
    BT.Transpose(B);

    TMatrixD    help(B*U*BT);
    A   -= help;
    A   *= V0;

    TMatrixD    I(nPar, nPar);
    I.UnitMatrix();
    I   -= A;
    V.Mult(V0, I);

    /*std::cout << "lambda: " << lambda[0][0] << std::endl;
    std::cout << "lambda: " << lambda[1][0] << std::endl;
    std::cout << "lambda: " << lambda[2][0] << std::endl;
    std::cout << "lambda: " << lambda[3][0] << std::endl;
    std::cout << "lambda: " << lambda[4][0] << std::endl;
    std::cout << "lambda: " << lambda[5][0] << std::endl;
    std::cout << "lambda: " << lambda[6][0] << std::endl;*/

    return kTRUE;
}
