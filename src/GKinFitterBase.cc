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
    par1(nPar,1),
    par2(nPar,1),
    unk1(0),
    unk2(0),
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
    TLorentzVector  eta(GetPhoton(par2, unk2, 0));
                    eta     += GetPhoton(par2, unk2, 1);
    TLorentzVector  pi0a(GetPhoton(par2, unk2, 2));
                    pi0a    += GetPhoton(par2, unk2, 3);
    TLorentzVector  pi0b(GetPhoton(par2, unk2, 4));
                    pi0b    += GetPhoton(par2, unk2, 5);
    TLorentzVector  etap(eta);
                    etap    += pi0a;
                    etap    += pi0b;

    hists.invMass.Fill(etap.M(), countIter);
    hists.invMassSub0.Fill(eta.M(), countIter);
    hists.invMassSub1.Fill(pi0a.M(), countIter);
    hists.invMassSub2.Fill(pi0b.M(), countIter);

    hists.chiSq.Fill(chiSq, countIter);
    hists.confidenceLevel.Fill(ConfidenceLevel(), countIter);
    hists.zVertex.Fill(unk2, countIter);
    for(int i=0; i<nPart; i++)
    {
        hists.photonsEnergy[i]->Fill(par2[(i+1)*GKinFitterBase_ParametersPerParticle][0], countIter);
        hists.photonsTheta[i]->Fill(par2[((i+1)*GKinFitterBase_ParametersPerParticle)+1][0] * TMath::RadToDeg(), countIter);
        hists.photonsPhi[i]->Fill(par2[((i+1)*GKinFitterBase_ParametersPerParticle)+2][0] * TMath::RadToDeg(), countIter);
    }
    hists.beamEnergy.Fill(par2[0][0], countIter);
    hists.beamTheta.Fill(par2[1][0] * TMath::RadToDeg(), countIter);
    hists.beamPhi.Fill(par2[2][0] * TMath::RadToDeg(), countIter);
    int indices[2];
    indices[0]  = 0;
    indices[1]  = 1;
    hists.con[0]->Fill(GetIMConstraint(par2, unk2, 0), countIter);
    indices[0]  = 2;
    indices[1]  = 3;
    hists.con[1]->Fill(GetIMConstraint(par2, unk2, 1), countIter);
    indices[0]  = 4;
    indices[1]  = 5;
    hists.con[2]->Fill(GetIMConstraint(par2, unk2, 2), countIter);
    for(int i=0; i<nPar; i++)
        hists.pulls[i]->Fill(Pull(i), countIter);

   TLorentzVector  help(GetMMConstraint(par2, unk2));
   hists.con[3]->Fill(help.Px(), countIter);
   hists.con[4]->Fill(help.Py(), countIter);
   hists.con[5]->Fill(help.Pz(), countIter);
   hists.con[6]->Fill(help.E(), countIter);
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


Double_t        GKinFitterBase::GetIMConstraint(const TMatrixD& par, const Double_t unk, const int i)    const
{
    if(i>=nImConstraint || i<0)
        return 0.0;

    TLorentzVector  tot(GetPhoton(par, unk, imConstraint[i].PartList[0]));
    for(int c=1; c<imConstraint[i].nPartList; c++)
        tot +=  GetPhoton(par, unk, imConstraint[i].PartList[c]);

    return tot.M2() - (imConstraint[i].mass*imConstraint[i].mass);
}
void       GKinFitterBase::GetIMConstraintDerPar(TMatrixD& ret, const TMatrixD& par, const Double_t unk, const int i)
{
    for(int c=0; c<nPar; c++)
        ret[c][0]   = 0;

    if(i>=nImConstraint || i<0)
        return;

    TLorentzVector  tot(GetPhoton(par, unk, imConstraint[i].PartList[0]));
    for(int c=1; c<imConstraint[i].nPartList; c++)
        tot +=  GetPhoton(par, unk, imConstraint[i].PartList[c]);

    for(int c=0; c<imConstraint[i].nPartList; c++)
    {
        ret[ (imConstraint[i].PartList[c]+1)*GKinFitterBase_ParametersPerParticle   ][0]   = 2 * tot * GetPhotonDerivateEnergy(par, unk, imConstraint[i].PartList[c]);
        ret[((imConstraint[i].PartList[c]+1)*GKinFitterBase_ParametersPerParticle)+1][0]   = 2 * tot * GetPhotonDerivateTheta(par, unk, imConstraint[i].PartList[c]);
        ret[((imConstraint[i].PartList[c]+1)*GKinFitterBase_ParametersPerParticle)+2][0]   = 2 * tot * GetPhotonDerivatePhi(par, unk, imConstraint[i].PartList[c]);
    }
}
Double_t   GKinFitterBase::GetIMConstraintDerUnk(const TMatrixD& par, const Double_t unk, const int i)          const
{
    if(i>=nImConstraint || i<0)
        return 0.0;

    TLorentzVector  tot(GetPhoton(par, unk, imConstraint[i].PartList[0]));
    for(int c=1; c<imConstraint[i].nPartList; c++)
        tot +=  GetPhoton(par, unk, imConstraint[i].PartList[c]);

    Double_t    ret = 0;
    for(int c=0; c<imConstraint[i].nPartList; c++)
        ret += 2* tot * GetPhotonDerivateZVertex(par, unk, imConstraint[i].PartList[c]);

    return ret;
}

Double_t        GKinFitterBase::GetMMConstraint(const TMatrixD& par, const Double_t unk)    const
{
    TLorentzVector  tot(0.0, 0.0, 0.0, targetMass);
    tot += GetBeam(par, unk);

    for(int i=0; i<nPart; i++)
        tot -=  GetPhoton(par, unk, i);

    return tot.M2() - (mmConstraint*mmConstraint);
}
void            GKinFitterBase::GetMMConstraintDerPar(TMatrixD& ret, const TMatrixD& par, const Double_t unk)
{
    TLorentzVector  tot(0.0, 0.0, 0.0, targetMass);
    tot += GetBeam(par, unk);

    for(int i=0; i<nPart; i++)
        tot -=  GetPhoton(par, unk, i);

    ret[0][0]   = 2 * tot * GetBeamDerivateEnergy(par, unk);
    ret[1][0]   = 2 * tot * GetBeamDerivateTheta(par, unk);
    ret[2][0]   = 2 * tot * GetBeamDerivatePhi(par, unk);
    for(int c=0; c<nPart; c++)
    {
        ret[ (c+1)*GKinFitterBase_ParametersPerParticle   ][0]   = 2 * tot * GetPhotonDerivateEnergy(par, unk, c);
        ret[((c+1)*GKinFitterBase_ParametersPerParticle)+1][0]   = 2 * tot * GetPhotonDerivateTheta(par, unk, c);
        ret[((c+1)*GKinFitterBase_ParametersPerParticle)+2][0]   = 2 * tot * GetPhotonDerivatePhi(par, unk, c);
    }
}
Double_t        GKinFitterBase::GetMMConstraintDerUnk(const TMatrixD& par, const Double_t unk)   const
{
    TLorentzVector  tot(0.0, 0.0, 0.0, targetMass);
    tot += GetBeam(par, unk);

    for(int i=0; i<nPart; i++)
        tot -=  GetPhoton(par, unk, i);

    Double_t    ret = 0;
    ret += 2 * tot * GetBeamDerivateZVertex(par, unk);
    for(int c=0; c<nPart; c++)
        ret += 2 * tot * GetPhotonDerivateZVertex(par, unk, c);

    return ret;
}

void    GKinFitterBase::GetG(TMatrixD& ret, const TMatrixD& par, const Double_t unk)
{
    int count   = 0;
    for(int i=0; i<nImConstraint; i++)
    {
        ret[count][0]   = GetIMConstraint(par, unk, i);
        count++;
    }
    if(isMmConstraint==kTRUE)
        ret[count][0]   = GetMMConstraint(par, unk);
}
void    GKinFitterBase::GetGDerPar(TMatrixD& ret, const TMatrixD& par, const Double_t unk)
{
    TMatrixD    help(nPar, 1);
    int count   = 0;
    for(int i=0; i<nImConstraint; i++)
    {
        GetIMConstraintDerPar(help, par, unk, i);
        for(int m=0; m<nPar; m++)
            ret[count][m]   = help[m][0];
        count++;
    }
    if(isMmConstraint==kTRUE)
    {
        GetMMConstraintDerPar(help, par, unk);
        for(int m=0; m<nPar; m++)
            ret[count][m]   = help[m][0];
    }
}
void    GKinFitterBase::GetGDerUnk(TMatrixD& ret, const TMatrixD &par, const Double_t unk)
{
    int count   = 0;
    for(int i=0; i<nImConstraint; i++)
    {
        ret[count][0]   = GetIMConstraintDerUnk(par, unk, i);
        count++;
    }
    if(isMmConstraint==kTRUE)
        ret[count][0]   = GetMMConstraintDerUnk(par, unk);
}

void    GKinFitterBase::GetR(TMatrixD& ret, const TMatrixD &par, const Double_t unk)
{
    TMatrixD    diff(par0);
    diff    -=  par;

    ret.Mult(GPar, diff);
    GetG(G, par, unk);
    ret += G;
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

Bool_t          GKinFitterBase::SolveStep(GHistFit2& hists)
{
    par1    = par2;
    unk1    = unk2;

    GetGDerPar(GPar, par1, unk1);
    {
        TMatrixD    help(nPar, nCon);
        help.MultT(V0, GPar);
        SInv.Mult(GPar, help);
    }
    Double_t    det = SInv.Determinant();
    if(det==0)
    {
        //std::cout << "Can not invert S." << std::endl;
        return kFALSE;
    }
    SInv.InvertFast();

    GetGDerUnk(GUnk, par1, unk1);
    {
        TMatrixD    help1(nCon, 1);
        help1.Mult(SInv, GUnk);
        TMatrixD    help2(1, 1);
        help2.TMult(GUnk, help1);
        U   = 1/help2[0][0];
    }

    TMatrixD    R(nCon, 1);
    GetR(R, par1, unk1);
    //std::cout << G[0][0] << std::endl;
    //std::cout << G[1][0] << std::endl;

    {
        TMatrixD    help1(nCon, 1);
        help1.Mult(SInv, R);
        TMatrixD    help2(1, 1);
        help2.TMult(GUnk, help1);
        unk2    = unk1 -(U * help2[0][0]);
    }

    {
        TMatrixD    help(GUnk);
        help   *= unk2 - unk1;
        help   += R;
        lambda.Mult(SInv, help);
    }

    {
        TMatrixD    help1(nPar, 1);
        help1.TMult(GPar, lambda);
        TMatrixD    help2(nPar, 1);
        help2.Mult(V0, help1);
        par2     = par1;
        par2    -= help2;
    }
    //std::cout << unk2 << std::endl;

    //chiSq
    {
        TMatrixD    help1(1, 1);
        help1.TMult(lambda, G);
        chiSq    = help1[0][0];
        chiSq   *= 2;
        TMatrixD    diff(par0);
        diff    -= par2;
        TMatrixD    help2(nPar, 1);
        help2.Mult(V0, diff);
        help1.TMult(diff, help2);
        chiSq   += help1[0][0];
    }

    //V
    TMatrixD    A(nPar, nPar);
    {
        TMatrixD    help(nCon, nPar);
        help.Mult(SInv, GPar);
        A.TMult(GPar, help);
    }

    TMatrixD    B(nPar, 1);
    {
        TMatrixD    help(nCon, 1);
        help.Mult(SInv, GUnk);
        B.TMult(GPar, help);
    }
    {
        TMatrixD    help1(nPar, nPar);
        help1.MultT(B, B);
        help1   *= U;
        TMatrixD    help2(A);
        help2   -= help1;
        A.Mult(help2, V0);
        TMatrixD    help3(nPar, nPar);
        help3.UnitMatrix();
        help3   -= A;
        V.Mult(V0, help3);
    }


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

    //Calc Par2 and Unk2
    unk2    = 0;
    Double_t    fValOld;
    Double_t    fVal;
    if(nImConstraint>0)
    {
        fValOld = GetIMConstraint(par0, unk2, 0);
        unk2    = unk2 - fValOld/GetIMConstraintDerUnk(par0, unk2, 0);
        fVal    = GetIMConstraint(par0, unk2, 0);
        //std::cout << "Only " << unk2 << " Constraints filled. " << fVal << " needed." << fValOld << " needed." << std::endl;
        while(fVal<fValOld)
        {
            fValOld = fVal;
            unk2    = unk2 - fValOld/GetIMConstraintDerUnk(par0, unk2, 0);
            fVal    = GetIMConstraint(par0, unk2, 0);
            //std::cout << "Only " << unk2 << " Constraints filled. " << fVal << " needed." << fValOld << " needed." << std::endl;
        }
    }
    par2    = par0;

    //Start solving
    if(SolveStep(hists)==kFALSE)
        return kFALSE;
    FillHists(hists);
    Double_t    oldChiSq    = chiSq;
    countIter   = 1;
    //std::cout << "ChiSq: " << chiSq << std::endl;
    if(SolveStep(hists)==kFALSE)
    {
        std::cout << "countIter: " << countIter << std::endl;
        par2    = par1;
        unk2    = unk1;
        return kTRUE;
    }
    //std::cout << "ChiSq: " << chiSq << std::endl;
    while(countIter<=10 && chiSq<oldChiSq)
    {
        //std::cout << "ChiSq: " << chiSq << std::endl;
        FillHists(hists);
        oldChiSq    = chiSq;
        countIter++;
        if(SolveStep(hists)==kFALSE)
        {
            //std::cout << "countIter: " << countIter << std::endl;
            par2    = par1;
            unk2    = unk1;
            return kTRUE;
        }
    }
    if(countIter==GKinFitterBase_MaxSteps)
        std::cout << "Fitter not finished. " << GKinFitterBase_MaxSteps << " steps are not enough." << std::endl;
    //std::cout << "countIter: " << countIter << std::endl;
    par2    = par1;
    unk2    = unk1;
    return kTRUE;
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
