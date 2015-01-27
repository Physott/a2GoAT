//////////////////////////////////////////////////////////////////
//GKinFitter
//////////////////////////////////////////////////////////////////

#include "GKinFitter.h"
#include "TDecompLU.h"
#include <iostream>

//ClassImp(GKinFitter);

//-----------------------------------------------------------------------------
GKinFitter::GKinFitter(const Int_t npart, const Int_t ncon){

  fNvar=4;
  fNpart=npart;
  fNpar = fNpart*fNvar;
  fNcon=ncon;
  fNparti=0;
  fNpari=0;
  fNconi=0;
  fmUnk0    = 0;
  fmUnk     = 0;

  fmAlpha0.ResizeTo(fNpar,1);
  fmAlpha1.ResizeTo(fNpar,1);
  fmAlpha2.ResizeTo(fNpar,1);
  fmV_Alpha0.ResizeTo(fNpar,fNpar);
  fmV_Alpha.ResizeTo(fNpar,fNpar);
  fmD.ResizeTo(fNcon,fNpar);
  fmDUnk.ResizeTo(fNcon,1);
  fmd.ResizeTo(fNcon,1);
  fmlamda.ResizeTo(fNcon,1);
  fmV_D.ResizeTo(fNcon,fNcon);
}

//-----------------------------------------------------------------------------
void GKinFitter::ResetMatrices(){

  fmAlpha0.Zero();
  fmAlpha1.Zero();
  fmAlpha2.Zero();
  fmV_Alpha0.Zero();
  fmV_Alpha.Zero();
  fmD.Zero();
  fmDUnk.Zero();
  fmd.Zero();
  fmlamda.Zero();
  fmV_D.Zero();
}

//-----------------------------------------------------------------------------
Int_t GKinFitter::Solve(){


    static int    iii=0;
  //Solve according to algorithm of Paul Avery:
  //Applied Fitting Theory VI, Formulas for Kinematic Fitting
  //see www.phys.ufl.edu/~avery/fitting.html

  if(fNpart!=fNparti){
    std::cout<<"GKinFitter::Solve() Added wrong number of particles. KinFit not completed"<<std::endl;
    return -1;
  }

  //fmD.Print();
  TMatrixD mDT(fmD);
  mDT.T();
  TMatrixD mV_Dinv(fmD*fmV_Alpha0*mDT);
  fmV_D=mV_Dinv;
  TDecompLU lu(fmV_D);
  if(!lu.Decompose()){
    std::cout<<"GKinFitter::Solve() Cannot invert. KinFit not completed"<<std::endl;
    return -1;
  }
  fmV_D.Invert();


  //Derive unknowns
  TMatrixD  mDUnkT(fmDUnk);
  mDUnkT.T();
  TMatrixD  mV_DUnkinv(mDUnkT*fmV_D*fmDUnk);
  Double_t  U   = 1/mV_DUnkinv[0][0];
  TMatrixD  diff1(fmAlpha0);
            diff1   -= fmAlpha1;
  TMatrixD  R(fmD*diff1);
            R   += fmd;
  TMatrixD  help1(mDT*fmV_D*R);
            help1   *= U;
  fmUnk = fmUnk0 - help1[0][0];
  std::cout << fmUnk << "   " << U << "   " << iii << std::endl;
  if(iii>=148)
  {
      std::cout << fmUnk << "   " << U << "   " << iii << std::endl;
      mDUnkT.Print();
      fmV_D.Print();
      fmDUnk.Print();
  }
  iii++;

  //Derive langrian multipliers
  fmlamda=fmV_D*(R+(fmDUnk*(fmUnk - fmUnk0)));

  //New parameters
  fmAlpha2=fmAlpha1-fmV_Alpha0*mDT*fmlamda;

  //New Covariant matrix
  TMatrixD  A(mDT*fmV_D*fmD);
  TMatrixD  B(mDT*fmV_D*fmDUnk);
  TMatrixD  BT(B);
  BT.T();
  fmV_Alpha=fmV_Alpha0-fmV_Alpha0*(A - (U*B*BT))*fmV_Alpha0;

  //chi2
  TMatrixD mlamdaT(fmlamda);
  mlamdaT.T();
  TMatrixD mchi2(mlamdaT*fmd);
  Cchi2=mchi2[0][0];

  TMatrixD diff2(fmAlpha0);
  diff2  -= fmAlpha2;
  TMatrixD diff2T(diff2);
  diff2T.T();
  TMatrixD mvchi2(diff2T*fmV_Alpha0*diff2);
  Vchi2=mvchi2[0][0];

  return 1;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddInvMassConstraint(const Double_t Minv)
{
    TLorentzVector ptot(0.0,0.0,0.0,0.0);
    for(Int_t i=0; i<fNpart; i++)
        ptot+=GetInitialParticle(i);


  // d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=ptot.M2()-Minv*Minv;

  // D matrix (derivitives of constraint eqn)
  for(Int_t i=0; i<fNpart; i++)
  {
    Double_t    factor          = GetInitialParticleFactor(i);
    Double_t    factorDerEnergy = GetInitialParticleDerivateEnergyFactor(i);
    Double_t    factorDerPz     = GetInitialParticleDerivatePzFactor(i);
    Double_t    factorDerUnk    = GetInitialParticleDerivateUnkFactor(i);

    //[Cons Number][Var Number]
    fmD[fNconi][0+i*fNvar]=-2*ptot.X()*factor;
    fmD[fNconi][1+i*fNvar]=-2*ptot.Y()*factor;
    fmD[fNconi][2+i*fNvar] =-2*ptot.X()* factorDerPz* fmAlpha1[ i*fNvar   ][0];
    fmD[fNconi][2+i*fNvar]+=-2*ptot.Y()* factorDerPz* fmAlpha1[(i*fNvar)+1][0];
    fmD[fNconi][2+i*fNvar]+=-2*ptot.Z()* factorDerPz*(fmAlpha1[(i*fNvar)+2][0]-fmUnk0);
    fmD[fNconi][2+i*fNvar]+=-2*ptot.Z()* factor;

    fmD[fNconi][3+i*fNvar] =-2*ptot.X()*factorDerEnergy* fmAlpha1[ i*fNvar   ][0];
    fmD[fNconi][3+i*fNvar]+=-2*ptot.Y()*factorDerEnergy* fmAlpha1[(i*fNvar)+1][0];
    fmD[fNconi][3+i*fNvar]+=-2*ptot.Z()*factorDerEnergy*(fmAlpha1[(i*fNvar)+2][0]-fmUnk0);
    fmD[fNconi][3+i*fNvar]+= 2*ptot.T();

    fmDUnk[fNconi][0]      =-2*ptot.X()*factorDerUnk* fmAlpha1[ i*fNvar   ][0];
    fmDUnk[fNconi][0]     +=-2*ptot.Y()*factorDerUnk* fmAlpha1[(i*fNvar)+1][0];
    fmDUnk[fNconi][0]     +=-2*ptot.Z()*factorDerUnk*(fmAlpha1[(i*fNvar)+2][0]-fmUnk0);
    fmDUnk[fNconi][0]     += 2*ptot.Z()*factor;
  }

  //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddSubInvMassConstraint(const Int_t Np, const Int_t pid[], const Double_t Minv){

  //Add invariant mass constraint to subset of particles
  //Np is number of subs particles
  //pid[] contains the particle number i.e the order when they were added

  if(Np>fNpart){
    std::cout<<"GKinFitter::AddSubInvMassConstraint too many particles!"<<std::endl;
    return;
  }

  //Add up the particle 4 vectors
  TLorentzVector ptot(0.0,0.0,0.0,0.0);
  for(Int_t i=0; i<Np; i++){
    ptot+=GetInitialParticle(pid[i]);
  }
  //d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=ptot.M2()-Minv*Minv;

  //D matrix (derivitives of constraint eqn)
  for(Int_t i=0;i<Np;i++)
  {
    Double_t    factor          = GetInitialParticleFactor(i);
    Double_t    factorDerEnergy = GetInitialParticleDerivateEnergyFactor(i);
    Double_t    factorDerPz     = GetInitialParticleDerivatePzFactor(i);
    Double_t    factorDerUnk    = GetInitialParticleDerivateUnkFactor(i);

    //[Cons Number][Var Number]
    fmD[fNconi][0+pid[i]*fNvar] =-2*ptot.X()*factor;
    fmD[fNconi][1+pid[i]*fNvar] =-2*ptot.Y()*factor;
    fmD[fNconi][2+pid[i]*fNvar] =-2*ptot.X()* factorDerPz* fmAlpha1[ pid[i]*fNvar   ][0];
    fmD[fNconi][2+pid[i]*fNvar]+=-2*ptot.Y()* factorDerPz* fmAlpha1[(pid[i]*fNvar)+1][0];
    fmD[fNconi][2+pid[i]*fNvar]+=-2*ptot.Z()* factorDerPz*(fmAlpha1[(pid[i]*fNvar)+2][0]-fmUnk0);
    fmD[fNconi][2+pid[i]*fNvar]+= 2*ptot.Z()* factor;

    fmD[fNconi][3+pid[i]*fNvar] =-2*ptot.X()*factorDerEnergy* fmAlpha1[ pid[i]*fNvar   ][0];
    fmD[fNconi][3+pid[i]*fNvar]+=-2*ptot.Y()*factorDerEnergy* fmAlpha1[(pid[i]*fNvar)+1][0];
    fmD[fNconi][3+pid[i]*fNvar]+=-2*ptot.Z()*factorDerEnergy*(fmAlpha1[(pid[i]*fNvar)+2][0]-fmUnk0);
    fmD[fNconi][3+pid[i]*fNvar]+= 2*ptot.T();

    fmDUnk[fNconi][0]           =-2*ptot.X()*factorDerUnk* fmAlpha1[ pid[i]*fNvar   ][0];
    fmDUnk[fNconi][0]          +=-2*ptot.Y()*factorDerUnk* fmAlpha1[(pid[i]*fNvar)+1][0];
    fmDUnk[fNconi][0]          +=-2*ptot.Z()*factorDerUnk*(fmAlpha1[(pid[i]*fNvar)+2][0]-fmUnk0);
    fmDUnk[fNconi][0]          += 2*ptot.Z()*factor;
  }

  //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddTotEnergyConstraint(const Double_t Etot)
{
    TLorentzVector ptot(0.0,0.0,0.0,0.0);
    for(Int_t i=0; i<fNpart; i++)
        ptot+=GetInitialParticle(i);

  //d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=ptot.E()-Etot;

  //D matrix (derivitives of constraint eqn)
  for(Int_t i=0; i<fNpart; i++)
  {
    //[Cons Number][Var Number]
    fmD[fNconi][0+i*fNvar]=0;
    fmD[fNconi][1+i*fNvar]=0;
    fmD[fNconi][2+i*fNvar]=0;
    fmD[fNconi][3+i*fNvar]=1;

    fmDUnk[fNconi][0]      = 0;
  }

  //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddTotMomentumConstraint(TVector3 mom)
{
    TLorentzVector ptot(0.0,0.0,0.0,0.0);
    for(Int_t i=0; i<fNpart; i++)
        ptot+=GetInitialParticle(i);

  //d matrix (evaluate constraint eqn.)
  fmd[fNconi+0][0]=ptot.X()-mom.X();
  fmd[fNconi+1][0]=ptot.Y()-mom.Y();
  fmd[fNconi+2][0]=ptot.Z()-mom.Z();

  //D matrix (derivitives of constraint eqn)
  //Double_t D[3][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0}};
  for(Int_t i=0; i<fNpart; i++)
  {
    Double_t    factor          = GetInitialParticleFactor(i);
    Double_t    factorDerEnergy = GetInitialParticleDerivateEnergyFactor(i);
    Double_t    factorDerPz     = GetInitialParticleDerivatePzFactor(i);
    Double_t    factorDerUnk    = GetInitialParticleDerivateUnkFactor(i);

    //[Cons Number][Var Number]
    fmD[fNconi  ][0+i*fNvar]=factor;
    fmD[fNconi  ][1+i*fNvar]=0;
    fmD[fNconi  ][2+i*fNvar]=factorDerPz*    fmAlpha1[0+i*fNvar][0];
    fmD[fNconi  ][3+i*fNvar]=factorDerEnergy*fmAlpha1[0+i*fNvar][0];

    fmD[fNconi+1][0+i*fNvar]=0;
    fmD[fNconi+1][1+i*fNvar]=factor;
    fmD[fNconi+1][2+i*fNvar]=factorDerPz*    fmAlpha1[1+i*fNvar][0];
    fmD[fNconi+1][3+i*fNvar]=factorDerEnergy*fmAlpha1[1+i*fNvar][0];

    fmD[fNconi+2][0+i*fNvar] =0;
    fmD[fNconi+2][1+i*fNvar] =0;
    fmD[fNconi+2][2+i*fNvar] =factorDerPz*    (fmAlpha1[2+i*fNvar][0]-fmUnk0);
    fmD[fNconi+2][2+i*fNvar]-=factor;
    fmD[fNconi+2][3+i*fNvar] =factorDerEnergy*(fmAlpha1[2+i*fNvar][0]-fmUnk0);

    fmDUnk[fNconi  ][0]      = factorDerUnk* fmAlpha1[ i*fNvar   ][0];
    fmDUnk[fNconi+1][0]      = factorDerUnk* fmAlpha1[(i*fNvar)+1][0];
    fmDUnk[fNconi+2][0]      = factorDerUnk*(fmAlpha1[(i*fNvar)+2][0]-fmUnk0);
    fmDUnk[fNconi+2][0]     += -factor;
  }
  fNconi+=3;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddSubMissMassConstraint(const TLorentzVector Mom, const Int_t Np, const Int_t pid[], const Double_t MissMass){

  // Add missing mass constraint to subset of particles
  // Np is number of subset particles
  // pid[] contains the particle number i.e the order when they were added

  if(Np>fNpart){
    std::cout<<"GKinFitter::AddSubMissMassConstraint too many particles!"<<std::endl;
    return;
  }

  //Add up the particle 4 vectors
  TLorentzVector Ptot(0.0,0.0,0.0,0.0);
  for(Int_t i=0; i<Np; i++){
    Ptot += GetInitialParticle(pid[i]);
  }

  //d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=(Mom-Ptot).M2()-MissMass*MissMass;

  //D matrix (derivitives of constraint eqn)
  for(Int_t i=0 ;i<Np ;i++)
  {
    Double_t    factor          = GetInitialParticleFactor(i);
    Double_t    factorDerEnergy = GetInitialParticleDerivateEnergyFactor(i);
    Double_t    factorDerPz     = GetInitialParticleDerivatePzFactor(i);
    Double_t    factorDerUnk    = GetInitialParticleDerivateUnkFactor(i);

    //[Cons Number][Var Number]
    fmD[fNconi][0+pid[i]*fNvar] =-2*(Ptot-Mom).X()* factor;
    fmD[fNconi][1+pid[i]*fNvar] =-2*(Ptot-Mom).Y()* factor;
    fmD[fNconi][2+pid[i]*fNvar] =-2*(Ptot-Mom).X()* factorDerPz * fmAlpha1[0+pid[i]*fNvar][0];
    fmD[fNconi][2+pid[i]*fNvar]+=-2*(Ptot-Mom).Y()* factorDerPz * fmAlpha1[1+pid[i]*fNvar][0];
    fmD[fNconi][2+pid[i]*fNvar]+=-2*(Ptot-Mom).Z()* factorDerPz *(fmAlpha1[2+pid[i]*fNvar][0]-fmUnk0);
    fmD[fNconi][2+pid[i]*fNvar]+=-2*(Ptot-Mom).Z()* factor;
    fmD[fNconi][3+pid[i]*fNvar] =-2*(Ptot-Mom).X()* factorDerEnergy * fmAlpha1[0+pid[i]*fNvar][0];
    fmD[fNconi][3+pid[i]*fNvar]+=-2*(Ptot-Mom).Y()* factorDerEnergy * fmAlpha1[1+pid[i]*fNvar][0];
    fmD[fNconi][3+pid[i]*fNvar]+=-2*(Ptot-Mom).Z()* factorDerEnergy *(fmAlpha1[2+pid[i]*fNvar][0]-fmUnk0);
    fmD[fNconi][3+pid[i]*fNvar]+= 2*(Ptot-Mom).T();

    fmDUnk[fNconi][0]           =-2*(Ptot-Mom).X()*factorDerUnk* fmAlpha1[ pid[i]*fNvar   ][0];
    fmDUnk[fNconi][0]          +=-2*(Ptot-Mom).Y()*factorDerUnk* fmAlpha1[(pid[i]*fNvar)+1][0];
    fmDUnk[fNconi][0]          +=-2*(Ptot-Mom).Z()*factorDerUnk*(fmAlpha1[(pid[i]*fNvar)+2][0]-fmUnk0);
    fmDUnk[fNconi][0]          += 2*(Ptot-Mom).Z()*factor;
  }
    //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddPosKFParticle(GKinFitterParticle kfp)
{
  if(fNparti>fNpart){
    std::cout<<"GKinFitter::AddPosKFParticle already at max particles"<<std::endl;
    return;
  }

  //Add parameters to Alpha0
  fmAlpha0.SetSub(fNpari,0,kfp.GetAlpha());
  fmAlpha1.SetSub(fNpari,0,kfp.GetAlpha());
  //Add error matrix to V_Alpha0
  fmV_Alpha0.SetSub(fNpari,fNpari,kfp.GetVAlpha());

  //increment counters
  fNpari+=fNvar;
  fNparti++;
}


//-----------------------------------------------------------------------------
void GKinFitter::AddNegKFParticle(GKinFitterParticle kfp)
{
  AddPosKFParticle(kfp);
  for(int i=fNpari-fNvar; i<fNpari; i++)
  {
      fmAlpha0[i][0]    = -fmAlpha0[i][0];
      fmAlpha1[i][0]    = -fmAlpha1[i][0];
  }
}

//-----------------------------------------------------------------------------
TLorentzVector GKinFitter::GetTotalFitParticle()
{
  TLorentzVector mtot(0.0,0.0,0.0,0.0);

  //loop over the sub matrices in alpha and add to total
  for(Int_t i=0; i<fNpart; i++)
    mtot+=GetParticle(i);

  return mtot;
}

//-----------------------------------------------------------------------------
Double_t GKinFitter::GetParticleFactor(const Int_t ip)
{
    Double_t          help    = -(2*fmAlpha2[(ip*fNvar)+2][0]*fmUnk)/fmAlpha2[(ip*fNvar)+3][0];
                      help   += fmUnk*fmUnk;
                      help   += 1;
                      help    = sqrt(help);
    return  1/help;
}
//-----------------------------------------------------------------------------
TLorentzVector GKinFitter::GetParticle(Int_t ip){

  //Return the fitted particle that was added ith
  if(ip>fNpari){
    std::cout<<"GKinFitter::GetParticle particle not in fit"<<std::endl;
    return TLorentzVector(0.0,0.0,0.0,0.0);
  }

  TLorentzVector    ret(fmAlpha2[ip*fNvar][0], fmAlpha2[(ip*fNvar)+1][0], fmAlpha2[(ip*fNvar)+2][0] - fmUnk, 0);
                    ret    *= GetParticleFactor(ip);
                    ret.SetE(fmAlpha2[(ip*fNvar)+3][0]);
  return ret;
}

//-----------------------------------------------------------------------------
Double_t GKinFitter::GetInitialParticleFactor(const Int_t ip)
{
    Double_t          help    = -(2*fmAlpha1[(ip*fNvar)+2][0]*fmUnk0)/fmAlpha1[(ip*fNvar)+3][0];
                      help   += fmUnk0*fmUnk0;
                      help   += 1;
                      help    = sqrt(help);
    return  1/help;
}

//-----------------------------------------------------------------------------
TLorentzVector GKinFitter::GetInitialParticle(Int_t ip){

  //Return the unfitted particle that was added ith
  if(ip>fNpari){
    std::cout<<"GKinFitter::GetInitialParticle particle not in fit"<<std::endl;
    return TLorentzVector(0.0,0.0,0.0,0.0);
  }

  TLorentzVector    ret(fmAlpha1[ip*fNvar][0], fmAlpha1[(ip*fNvar)+1][0], fmAlpha1[(ip*fNvar)+2][0] - fmUnk0, 0);
                    ret    *= GetInitialParticleFactor(ip);
                    ret.SetE(fmAlpha1[(ip*fNvar)+3][0]);
  return ret;
}

//-----------------------------------------------------------------------------
Double_t       GKinFitter::GetInitialParticleDerivateEnergyFactor(const Int_t ip)
{
    Double_t    ret  = GetInitialParticleFactor(ip);
                ret *= ret*ret;
                ret /= fmAlpha1[(ip*fNvar)+3][0]*fmAlpha1[(ip*fNvar)+3][0];
                ret *= -fmAlpha1[(ip*fNvar)+2][0]*fmUnk0;
    return  ret;
}

//-----------------------------------------------------------------------------
Double_t       GKinFitter::GetInitialParticleDerivatePzFactor(const Int_t ip)
{
    Double_t    ret  = GetInitialParticleFactor(ip);
                ret *= ret*ret;
                ret /= fmAlpha1[(ip*fNvar)+3][0];
                ret *= fmUnk0;
    return  ret;
}

//-----------------------------------------------------------------------------
Double_t       GKinFitter::GetInitialParticleDerivateUnkFactor(const Int_t ip)
{
    Double_t    ret  = GetInitialParticleFactor(ip);
                ret *= ret*ret;
                ret *= (fmAlpha1[(ip*fNvar)+2][0]/fmAlpha1[(ip*fNvar)+3][0])-fmUnk0;
    return  ret;
}

//-----------------------------------------------------------------------------
TLorentzVector GKinFitter::GetOriginalParticle(Int_t ip){

  //Return the unfitted particle that was added ith
  if(ip>fNpari){
    std::cout<<"GKinFitter::GetInitialParticle particle not in fit"<<std::endl;
    return TLorentzVector(0.0,0.0,0.0,0.0);
  }

  TMatrixD mi(fNvar,1);
  mi=fmAlpha0.GetSub(ip*fNvar,(ip+1)*fNvar-1,0,0);
  return TLorentzVector(mi[0][0],mi[1][0],mi[2][0],mi[3][0]);
}

//-----------------------------------------------------------------------------
void GKinFitter::Debug(){

  std::cout<<"Alpha0 "<<std::endl;
  fmAlpha0.Print();
  std::cout<<"Alpha1 "<<std::endl;
  fmAlpha1.Print();
  std::cout<<"Alpha2 "<<std::endl;
  fmAlpha2.Print();
  std::cout<<"V_Alpha0 "<<std::endl;
  fmV_Alpha0.Print();
  std::cout<<"V_Alpha "<<std::endl;
  fmV_Alpha.Print();
  std::cout<<"d "<<std::endl;
  fmd.Print();
  std::cout<<"D "<<std::endl;
  fmD.Print();
  std::cout<<"V_D "<<std::endl;
  fmV_D.Print();
  std::cout<<"lamda "<<std::endl;
  fmlamda.Print();

  //Check D*deltaAlpha+d=0
  TMatrixD mdelAlpha=fmAlpha2-fmAlpha0;
  TMatrix mCheck1=fmD*mdelAlpha + fmd;
  std::cout<<"delAlpha"<<std::endl;
  mdelAlpha.Print();
  std::cout<<"Check1 "<<std::endl;
  mCheck1.Print();

}










GIterativeKinFitter::GIterativeKinFitter(const Int_t npart, const Int_t ncon) :
    GKinFitter(npart, ncon)
{
}

GIterativeKinFitter::~GIterativeKinFitter()
{

}

void GIterativeKinFitter::AddInvMassConstraint(const Double_t Minv)
{
    conType[fNconi] = ConstraintType_InvMass;
    var[fNconi]    = Minv;

    GKinFitter::AddInvMassConstraint(Minv);
}

void GIterativeKinFitter::AddSubInvMassConstraint(const Int_t Np, const Int_t pid[], const Double_t Minv)
{
    conType[fNconi] = ConstraintType_SubInvMass;
    nIndices[fNconi]      = Np;
    for(int i=0; i<Np; i++)
        indices[fNconi][i]  = pid[i];
    var[fNconi]    = Minv;

    GKinFitter::AddSubInvMassConstraint(Np, pid, Minv);
}

void GIterativeKinFitter::AddTotEnergyConstraint(const Double_t Etot)
{
    conType[fNconi] = ConstraintType_InvEnergy;
    var[fNconi]  = Etot;

    GKinFitter::AddTotEnergyConstraint(Etot);
}

void GIterativeKinFitter::AddTotMomentumConstraint(const TVector3 mom)
{
    conType[fNconi] = ConstraintType_InvMomentum;
    momentum[fNconi]     = mom;

    GKinFitter::AddTotMomentumConstraint(mom);
}

void GIterativeKinFitter::AddSubMissMassConstraint(const TLorentzVector Mom, const Int_t Np, const Int_t pid[], const Double_t MissMass)
{
    conType[fNconi] = ConstraintType_MisMass;
    beam[fNconi]    = Mom;
    nIndices[fNconi]      = Np;
    for(int i=0; i<Np; i++)
        indices[fNconi][i]  = pid[i];
    var[fNconi]    = MissMass;

    GKinFitter::AddSubMissMassConstraint(Mom, Np, pid, MissMass);
}

Int_t GIterativeKinFitter::Solve()
{
    if(nIter==0)
    {
        std::cout << "Solve Start" << std::endl;
        fmUnk0 = 0;
        if(GKinFitter::Solve()<0)
            return -1;
        nIter++;
        return 1;
    }
    if(nIter==5)
        return -1;

    fmAlpha1    = fmAlpha2;
    fmUnk0      = fmUnk;
    TMatrixD    V(fmV_Alpha);
    Double_t    CChiSq  = GetConstraintsChi2();
    Double_t    VChiSq  = GetVariablesChi2();

    ResetConstraints();
    fmD.Zero();
    fmDUnk.Zero();
    fmd.Zero();
    fmlamda.Zero();
    fmV_D.Zero();
    for(int i=0; i<fNcon; i++)
    {
        switch(conType[i])
        {
        case ConstraintType_InvMass:
            GKinFitter::AddInvMassConstraint(var[i]);
            break;
        case ConstraintType_SubInvMass:
            GKinFitter::AddSubInvMassConstraint(nIndices[i], indices[i], var[i]);
            break;
        case ConstraintType_InvEnergy:
            GKinFitter::AddTotEnergyConstraint(var[i]);
            break;
        case ConstraintType_InvMomentum:
            GKinFitter::AddTotMomentumConstraint(momentum[i]);
            i+=2;
            break;
        case ConstraintType_MisMass:
            GKinFitter::AddSubMissMassConstraint(beam[i], nIndices[i], indices[i], var[i]);
            break;
        }
    }

    if(GKinFitter::Solve()<0)
    {
        fmAlpha2    = fmAlpha1;
        fmUnk       = fmUnk0;
        fmV_Alpha   = V;
        Cchi2       = CChiSq;
        Vchi2       = VChiSq;
        return -1;
    }
    else
        nIter++;

    return 1;
}
