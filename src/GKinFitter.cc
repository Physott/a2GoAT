//////////////////////////////////////////////////////////////////
//GKinFitter
//////////////////////////////////////////////////////////////////

#include "GKinFitter.h"
#include "TDecompLU.h"
#include <iostream>

//ClassImp(GKinFitter);

//-----------------------------------------------------------------------------
GKinFitter::GKinFitter(const Int_t npart, const Int_t ncon, const Int_t unk){

  fNvar=4;
  fNpart=npart;
  fNpar = fNpart*fNvar;
  fNcon=ncon;
  fNparti=0;
  fNpari=0;
  fNconi=0;
  fNiter=0;
  fNunKnown=unk;

  fmAlpha0.ResizeTo(fNpar,1);
  fmAlpha.ResizeTo(fNpar,1);
  fmV_Alpha0.ResizeTo(fNpar,fNpar);
  fmV_Alpha.ResizeTo(fNpar,fNpar);
  fmD.ResizeTo(fNcon,fNpar);
  fmd.ResizeTo(fNcon,1);
  fmlamda.ResizeTo(fNcon,1);
  fmV_D.ResizeTo(fNcon,fNcon);
  fT.ResizeTo(fNpar,fNpar);

  fPtot.SetXYZT(0,0,0,0);

}

//-----------------------------------------------------------------------------
void GKinFitter::ResetMatrices(){

  fmAlpha0.Zero();
  fmAlpha.Zero();
  fmV_Alpha0.Zero();
  fmV_Alpha.Zero();
  fmD.Zero();
  fmd.Zero();
  fmlamda.Zero();
  fmV_D.Zero();

}

//-----------------------------------------------------------------------------
Int_t GKinFitter::Solve(){

  //Solve according to algorithm of Paul Avery:
  //Applied Fitting Theory VI, Formulas for Kinematic Fitting
  //see www.phys.ufl.edu/~avery/fitting.html

  if(fNpart!=fNparti){
    std::cout<<"GKinFitter::Solve() Added wrong number of particles. KinFit not completed"<<std::endl;
    return -1;
  }

  TMatrixD mDT=fmD;
  mDT.T();
  TMatrixD mV_Dinv=fmD*fmV_Alpha0*mDT;
  fmV_D=mV_Dinv;
  TDecompLU lu(fmV_D);
  if(!lu.Decompose()){
    std::cout<<"GKinFitter::Solve() Cannot invert. KinFit not completed"<<std::endl;
    return -1;
  }
  fmV_D.Invert();

  //Double_t det;
  //TMatrixD Einheit = fmV_D*mV_Dinv;
  //if(fNpart == 7) Einheit.Print();

  //Derive langrian multipliers
  fmlamda=fmV_D*fmd;
  //New parameters
  fmAlpha=fmAlpha0-fmV_Alpha0*mDT*fmlamda;
  //New Covariant matrix
  fmV_Alpha=fmV_Alpha0-fmV_Alpha0*mDT*fmV_D*fmD*fmV_Alpha0;
  //chi2
  TMatrixD mlamdaT=fmlamda;
  mlamdaT.T();
  TMatrixD mchi2=mlamdaT*fmd;
  fchi2=mchi2[0][0];
  fNiter++;

  return 1;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddInvMassConstraint(const Double_t Minv){

  // d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=fPtot.M2()-Minv*Minv;

  // D matrix (derivitives of constraint eqn)
  for(Int_t i=0; i<fNpart; i++){
    //[Cons Number][Var Number]
    fmD[fNconi][0+i*fNvar]=-2*fPtot.X();
    fmD[fNconi][1+i*fNvar]=-2*fPtot.Y();
    fmD[fNconi][2+i*fNvar]=-2*fPtot.Z();
    fmD[fNconi][3+i*fNvar]= 2*fPtot.T();
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
    ptot+=GetInitialParticle(pid[i]).Get4Vector();
  }

  //d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=ptot.M2()-Minv*Minv;

  //D matrix (derivitives of constraint eqn)
  for(Int_t i=0;i<Np;i++){
    //[Cons Number][Var Number]
    fmD[fNconi][0+pid[i]*fNvar]=-2*ptot.X();
    fmD[fNconi][1+pid[i]*fNvar]=-2*ptot.Y();
    fmD[fNconi][2+pid[i]*fNvar]=-2*ptot.Z();
    fmD[fNconi][3+pid[i]*fNvar]= 2*ptot.T();
  }

  //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddTotEnergyConstraint(const Double_t Etot){

  //d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=fPtot.E()-Etot;

  //D matrix (derivitives of constraint eqn)
  for(Int_t i=0; i<fNpart; i++){
    //[Cons Number][Var Number]
    fmD[fNconi][0+i*fNvar]=0;
    fmD[fNconi][1+i*fNvar]=0;
    fmD[fNconi][2+i*fNvar]=0;
    fmD[fNconi][3+i*fNvar]=1;
  }

  //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddTotMomentumConstraint(TVector3 mom){

  //d matrix (evaluate constraint eqn.)
  fmd[fNconi+0][0]=fPtot.X()-mom.X();
  fmd[fNconi+1][0]=fPtot.Y()-mom.Y();
  fmd[fNconi+2][0]=fPtot.Z()-mom.Z();

  //D matrix (derivitives of constraint eqn)
  Double_t D[3][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0}};
  for(Int_t i=0; i<fNpart; i++){
    for(Int_t j=0;j<3;j++){
      //[Cons Number][Var Number]
      fmD[fNconi+j][0+i*fNvar]=D[j][0];
      fmD[fNconi+j][1+i*fNvar]=D[j][1];
      fmD[fNconi+j][2+i*fNvar]=D[j][2];
      fmD[fNconi+j][3+i*fNvar]=D[j][3];
    }
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
    Ptot += GetInitialParticle(pid[i]).Get4Vector();
  }

  //d matrix (evaluate constraint eqn.)
  fmd[fNconi][0]=(Mom-Ptot).M2()-MissMass*MissMass;

  //D matrix (derivitives of constraint eqn)
  for(Int_t i=0 ;i<Np ;i++)
    {//[Cons Number][Var Number]
    fmD[fNconi][0+pid[i]*fNvar]=-2*(Ptot-Mom).X();
    fmD[fNconi][1+pid[i]*fNvar]=-2*(Ptot-Mom).Y();
    fmD[fNconi][2+pid[i]*fNvar]=-2*(Ptot-Mom).Z();
    fmD[fNconi][3+pid[i]*fNvar]= 2*(Ptot-Mom).T();
  }
    //increment constraint counter
  fNconi++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddPosKFParticle(GKinFitterParticle kfp){

  if(fNparti>fNpart){
    std::cout<<"GKinFitter::AddPosKFParticle already at max particles"<<std::endl;
    return;
  }

  //Add parameters to Alpha0
  fmAlpha0.SetSub(fNpari,0,kfp.GetAlpha());
  //Add error matrix to V_Alpha0
  fmV_Alpha0.SetSub(fNpari,fNpari,kfp.GetVAlpha());
  //Add transformation matrice to fT
  fT.SetSub(fNpari,fNpari,kfp.GetT());
  //ADD Lorentz Vector !!
  if(fNparti) fPtot=fPtot+kfp.Get4Vector();
  else fPtot=kfp.Get4Vector();

  //increment counters
  fNpari+=fNvar;
  fNparti++;

}

//-----------------------------------------------------------------------------
void GKinFitter::AddNegKFParticle(GKinFitterParticle kfp){

  if(fNparti>fNpart){
    std::cout<<"GKinFitter::AddPosKFParticle already at max particles"<<std::endl;
    return;
  }

  //Add parameters to Alpha0
  fmAlpha0.SetSub(fNpari,0,kfp.GetAlpha());
  //Add error matrix to V_Alpha0
  fmV_Alpha0.SetSub(fNpari,fNpari,kfp.GetVAlpha());
  //Add transformation matrice to fT
  fT.SetSub(fNpari,fNpari,kfp.GetT());
  //SUBTRACT Lorentz Vector !!
  if(fNparti) fPtot=fPtot-kfp.Get4Vector();
  else fPtot=kfp.Get4Vector();

  //increment counters
  fNpari+=fNvar;
  fNparti++;

}

//-----------------------------------------------------------------------------
GKinFitterParticle GKinFitter::GetTotalFitParticle(){

  GKinFitterParticle kfp;
  TMatrixD mtot(fNvar,1);

  //loop over the sub matrices in alpha and add to total
  for(Int_t i=0; i<fNpart; i++){
    mtot+=fmAlpha.GetSub(i*fNvar,(i+1)*fNvar-1,0,0);
  }

  //Set 4 vector, automatically sets alpha.
  kfp.Set4Vector(TLorentzVector(mtot[0][0],mtot[1][0],mtot[2][0],mtot[3][0]));
  // Add the error matrixes
  TMatrixD mV_tot(fNvar,fNvar);

  for(Int_t i=0; i<fNpart; i++){
    mV_tot+=fmV_Alpha.GetSub(i*fNvar,(i+1)*fNvar-1,i*fNvar,(i+1)*fNvar-1);
  }

  kfp.SetVAlpha(mV_tot);

  return kfp;

}

//-----------------------------------------------------------------------------
GKinFitterParticle GKinFitter::GetParticle(Int_t ip){

  //Return the fitted particle that was added ith
  if(ip>fNpari){
    std::cout<<"GKinFitter::GetParticle particle not in fit"<<std::endl;
    return GKinFitterParticle();
  }

  GKinFitterParticle kfp;
  TMatrixD mi(fNvar,1);

  mi=fmAlpha.GetSub(ip*fNvar,(ip+1)*fNvar-1,0,0);
  kfp.Set4Vector(TLorentzVector(mi[0][0],mi[1][0],mi[2][0],mi[3][0]));
  TMatrixD mVi(fNvar,fNvar);
  mVi=fmV_Alpha.GetSub(ip*fNvar,(ip+1)*fNvar-1,ip*fNvar,(ip+1)*fNvar-1);
  kfp.SetVAlpha(mVi);
  return kfp;

}

//-----------------------------------------------------------------------------
GKinFitterParticle GKinFitter::GetInitialParticle(Int_t ip){

  //Return the unfitted particle that was added ith
  if(ip>fNpari){
    std::cout<<"GKinFitter::GetInitialParticle particle not in fit"<<std::endl;
    return GKinFitterParticle();
  }

  GKinFitterParticle kfp;
  TMatrixD mi(fNvar,1);

  mi=fmAlpha0.GetSub(ip*fNvar,(ip+1)*fNvar-1,0,0);
  kfp.Set4Vector(TLorentzVector(mi[0][0],mi[1][0],mi[2][0],mi[3][0]));
  TMatrixD mVi(fNvar,fNvar);
  mVi=fmV_Alpha0.GetSub(ip*fNvar,(ip+1)*fNvar-1,ip*fNvar,(ip+1)*fNvar-1);
  kfp.SetVAlpha(mVi);
  return kfp;

}

//-----------------------------------------------------------------------------
void GKinFitter::Debug(){

  std::cout<<"Alpha0 "<<std::endl;
  fmAlpha0.Print();
  std::cout<<"Alpha "<<std::endl;
  fmAlpha.Print();
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
  std::cout<<"T "<<std::endl;
  fT.Print();

  //Check D*deltaAlpha+d=0
  TMatrixD mdelAlpha=fmAlpha-fmAlpha0;
  TMatrix mCheck1=fmD*mdelAlpha + fmd;
  std::cout<<"delAlpha"<<std::endl;
  mdelAlpha.Print();
  std::cout<<"Check1 "<<std::endl;
  mCheck1.Print();

}







GMyTrackM::GMyTrackM(const Double_t E, const Double_t Theta, const Double_t Phi, const Double_t Mass, const TVector3 &Vertex, const Double_t DE, const Double_t DTheta, const Double_t DPhi, const TVector3 &DVertex):
    e(E),
    theta(Theta),
    phi(Phi),
    mass(Mass),
    vertex(Vertex),
    de(DE),
    dtheta(DTheta),
    dphi(DPhi),
    dvertex(DVertex)
{

}

GMyTrackM::~GMyTrackM()
{
}

TMatrixD    GMyTrackM::GetParametersH()    const
{
    TMatrixD    ret(8,1);
    ret[0][0]   = GMyTrackH_nParameters * TMath::Sin(theta) * TMath::Cos(phi) - vertex.x();
    ret[1][0]   = GMyTrackH_nParameters * TMath::Sin(theta) * TMath::Sin(phi) - vertex.y();
    ret[2][0]   = GMyTrackH_nParameters * TMath::Cos(theta)                   - vertex.z();
    ret[3][0]   = e;
    ret[4][0]   = vertex.x();
    ret[5][0]   = vertex.y();
    ret[6][0]   = vertex.z();
    ret[7][0]   = mass;
}

TMatrixD    GMyTrackM::GetCovarianceH()    const
{
    TMatrixD    ret(GMyTrackH_nParameters, GMyTrackH_nParameters);
    Double_t    help    = GMyTrackM_CBRadius*TMath::Cos(theta)*TMath::Cos(phi)*dtheta;
    ret[0][0]   = help*help;
    help    = GMyTrackM_CBRadius*TMath::Sin(theta)*TMath::Sin(phi)*dphi;
    ret[0][0]   += help*help;
    ret[0][0]   += dvertex.x()*dvertex.x();

    help    = GMyTrackM_CBRadius*TMath::Cos(theta)*dtheta;
    ret[1][0]   = help*help*TMath::Sin(phi)*TMath::Cos(phi);
    help    = GMyTrackM_CBRadius*TMath::Sin(theta)*dphi;
    ret[1][0]   += help*help*TMath::Sin(phi)*TMath::Cos(phi);
    ret[0][1]   = ret[1][0];

    help    = GMyTrackM_CBRadius*TMath::Cos(theta)*TMath::Sin(phi)*dtheta;
    ret[1][1]   = help*help;
    help    = GMyTrackM_CBRadius*TMath::Sin(theta)*TMath::Cos(phi)*dphi;
    ret[1][1]   += help*help;
    ret[1][1]   += dvertex.y()*dvertex.y();

    help    = GMyTrackM_CBRadius*dtheta;
    ret[2][0]   = -help*help*TMath::Sin(theta)*TMath::Cos(theta);
    ret[0][2]   = ret[2][0];

    help    = GMyTrackM_CBRadius*dtheta;
    ret[2][1]   = help*help*TMath::Sin(theta)*TMath::Cos(theta);
    ret[1][2]   = ret[2][1];

    help    = GMyTrackM_CBRadius*TMath::Sin(theta)*dtheta;
    ret[2][2]   = help*help;
    ret[2][2]   += dvertex.z()*dvertex.z();

    ret[3][0]   = 0;    ret[0][3]   = 0;
    ret[3][1]   = 0;    ret[1][3]   = 0;
    ret[3][2]   = 0;    ret[2][3]   = 0;
    ret[3][3]   = de*de;

    ret[4][0]   = -dvertex.x()*dvertex.x(); ret[0][4]   = ret[4][0];
    ret[4][1]   = 0;                ret[1][4]   = 0;
    ret[4][2]   = 0;                ret[2][4]   = 0;
    ret[4][3]   = 0;                ret[3][4]   = 0;
    ret[4][4]   = dvertex.x()*dvertex.x();

    ret[5][0]   = 0;                ret[0][5]   = 0;
    ret[5][1]   = -dvertex.y()*dvertex.y(); ret[1][5]   = ret[5][1];
    ret[5][2]   = 0;                ret[2][5]   = 0;
    ret[5][3]   = 0;                ret[3][5]   = 0;
    ret[5][4]   = 0;                ret[4][5]   = 0;
    ret[5][5]   = dvertex.y()*dvertex.y();

    ret[6][0]   = 0;                ret[0][6]   = 0;
    ret[6][1]   = 0;                ret[1][6]   = 0;
    ret[6][2]   = -dvertex.z()*dvertex.z(); ret[2][6]   = ret[6][2];
    ret[6][3]   = 0;                ret[3][6]   = 0;
    ret[6][4]   = 0;                ret[4][6]   = 0;
    ret[6][5]   = 0;                ret[5][6]   = 0;
    ret[6][6]   = dvertex.z()*dvertex.z();

    return ret;
}




GMyTrackH::GMyTrackH(const GMyTrackM& m):
    helpVertex(),
    e(0),
    mass(0),
    vertex(),
    cm(m.GetCovarianceH())
{
    TMatrixD    p(m.GetParametersH());

    helpVertex.SetXYZ(p[0][0], p[1][0], p[2][0]);
    e   = p[3][0];
    mass= p[7][0];
    vertex.SetXYZ(p[4][0], p[5][0], p[6][0]);
}

GMyTrackH::~GMyTrackH()
{

}

TMatrixD    GMyTrackH::GetParametersW()    const
{
    TMatrixD    ret(7,1);
    Double_t    help    = TMath::Sqrt((e*e)-(mass*mass));
    help    /= helpVertex.Mag();
    ret[0][0]   = help * helpVertex.x();
    ret[1][0]   = help * helpVertex.y();
    ret[2][0]   = help * helpVertex.z();
    ret[3][0]   = e;
    ret[4][0]   = vertex.x();
    ret[5][0]   = vertex.y();
    ret[6][0]   = vertex.z();
}

TMatrixD    GMyTrackH::GetCovarianceW()    const
{
    Double_t    help    = TMath::Sqrt((e*e)-(mass*mass));
    Double_t    r       = helpVertex.Mag();

    TMatrixD    d(7,7);
    d[0][0] = e * helpVertex.X()/(help*r);
    d[0][1] = ((helpVertex.Y()*helpVertex.Y())+(helpVertex.z()*helpVertex.z()))*help/(r*r*r);
    d[0][2] = -helpVertex.x()*helpVertex.y()*help/(r*r*r);
    d[0][3] = -helpVertex.x()*helpVertex.z()*help/(r*r*r);
    d[0][4] = 0;
    d[0][5] = 0;
    d[0][6] = 0;

    d[1][0] = e * helpVertex.y()/(help*r);
    d[1][1] = -helpVertex.y()*helpVertex.x()*help/(r*r*r);
    d[1][2] = ((helpVertex.x()*helpVertex.x())+(helpVertex.z()*helpVertex.z()))*help/(r*r*r);
    d[1][3] = -helpVertex.y()*helpVertex.z()*help/(r*r*r);
    d[1][4] = 0;
    d[1][5] = 0;
    d[1][6] = 0;

    d[2][0] = e * helpVertex.z()/(help*r);
    d[2][1] = -helpVertex.z()*helpVertex.x()*help/(r*r*r);
    d[2][2] = -helpVertex.z()*helpVertex.y()*help/(r*r*r);
    d[2][3] = ((helpVertex.x()*helpVertex.x())+(helpVertex.y()*helpVertex.y()))*help/(r*r*r);
    d[2][4] = 0;
    d[2][5] = 0;
    d[2][6] = 0;

    d[3][0] = 0;
    d[3][1] = 0;
    d[3][2] = 0;
    d[3][3] = 1;
    d[3][4] = 0;
    d[3][5] = 0;
    d[3][6] = 0;

    d[4][0] = 0;
    d[4][1] = 0;
    d[4][2] = 0;
    d[4][3] = 0;
    d[4][4] = 1;
    d[4][5] = 0;
    d[4][6] = 0;

    d[5][0] = 0;
    d[5][1] = 0;
    d[5][2] = 0;
    d[5][3] = 0;
    d[5][4] = 0;
    d[5][5] = 1;
    d[5][6] = 0;

    d[6][0] = 0;
    d[6][1] = 0;
    d[6][2] = 0;
    d[6][3] = 0;
    d[6][4] = 0;
    d[6][5] = 0;
    d[6][6] = 1;

    TMatrixD    dt(d.Transpose(d));

    return d*cm*dt;
}

GMyKinFitter::GMyKinFitter(const Int_t npart, const Int_t ncon)    :
    nPart(npart),
    nPar(nPart*GMyKinFitter_fNvar),
    nCon(ncon),
    fNparti(0),
    fNpari(0),
    fNconi(0),
    fNiter(0),
    fmAlpha0(nPar,1),
    fmAlpha(nPar,1),
    fmV_Alpha0(nPar,nPar),
    fmV_Alpha(nPar,nPar),
    fmD(nCon,nPar),
    fmd(nCon,1),
    fmlamda(nCon,1),
    fmV_D(nCon,nCon),
    fT(nPar,nPar),
    fPtot(0, 0, 0, 0)
{
}

GMyKinFitter::~GMyKinFitter()
{

}

TMatrixD    GMyKinFitter::ConvertMtoH(const Double_t* pIn, const Double_t* dpIn, Double_t* pOut)
{
    TMatrixD    ret(GMyKinFitter_fNvar, GMyKinFitter_fNvar);
    Double_t    help    = GMyKinFitter_l*TMath::Cos(pIn[1])*TMath::Cos(pIn[2])*dpIn[1];
    ret[0][0]   = help*help;
    help    = GMyKinFitter_l*TMath::Sin(pIn[1])*TMath::Sin(pIn[2])*dpIn[2];
    ret[0][0]   += help*help;
    ret[0][0]   += dpIn[3]*dpIn[3];

    help    = GMyKinFitter_l*TMath::Cos(pIn[1])*dpIn[1];
    ret[1][0]   = help*help*TMath::Sin(pIn[2])*TMath::Cos(pIn[2]);
    help    = GMyKinFitter_l*TMath::Sin(pIn[1])*dpIn[2];
    ret[1][0]   += help*help*TMath::Sin(pIn[2])*TMath::Cos(pIn[2]);
    ret[0][1]   = ret[1][0];

    help    = GMyKinFitter_l*TMath::Cos(pIn[1])*TMath::Sin(pIn[2])*dpIn[1];
    ret[1][1]   = help*help;
    help    = GMyKinFitter_l*TMath::Sin(pIn[1])*TMath::Cos(pIn[2])*dpIn[2];
    ret[1][1]   += help*help;
    ret[1][1]   += dpIn[4]*dpIn[4];

    help    = GMyKinFitter_l*dpIn[1];
    ret[2][0]   = -help*help*TMath::Sin(pIn[1])*TMath::Cos(pIn[1]);
    ret[0][2]   = ret[2][0];

    help    = GMyKinFitter_l*dpIn[1];
    ret[2][1]   = help*help*TMath::Sin(pIn[1])*TMath::Cos(pIn[1]);
    ret[1][2]   = ret[2][1];

    help    = GMyKinFitter_l*TMath::Sin(pIn[1])*dpIn[1];
    ret[2][2]   = help*help;
    ret[2][2]   += dpIn[5]*dpIn[5];

    ret[3][0]   = 0;    ret[0][3]   = 0;
    ret[3][1]   = 0;    ret[1][3]   = 0;
    ret[3][2]   = 0;    ret[2][3]   = 0;
    ret[3][3]   = dpIn[0]*dpIn[0];

    ret[4][0]   = -dpIn[3]*dpIn[3]; ret[0][4]   = ret[4][0];
    ret[4][1]   = 0;                ret[1][4]   = 0;
    ret[4][2]   = 0;                ret[2][4]   = 0;
    ret[4][3]   = 0;                ret[3][4]   = 0;
    ret[4][4]   = 0;

    ret[5][0]   = 0;                ret[0][5]   = 0;
    ret[5][1]   = -dpIn[4]*dpIn[4]; ret[1][5]   = ret[5][1];
    ret[5][2]   = 0;                ret[2][5]   = 0;
    ret[5][3]   = 0;                ret[3][5]   = 0;
    ret[5][4]   = 0;                ret[4][5]   = 0;
    ret[5][5]   = 0;

    ret[6][0]   = 0;                ret[0][6]   = 0;
    ret[6][1]   = 0;                ret[1][6]   = 0;
    ret[6][2]   = -dpIn[5]*dpIn[5]; ret[2][6]   = ret[6][2];
    ret[6][3]   = 0;                ret[3][6]   = 0;
    ret[6][4]   = 0;                ret[4][6]   = 0;
    ret[6][5]   = 0;                ret[5][6]   = 0;
    ret[6][6]   = 0;

    return ret;
}

TMatrixD    GMyKinFitter::ConvertHtoW(const Double_t* pIn, const TMatrixD& dpIn, const Double_t mass, Double_t* pOut)
{
    TMatrixD    ret(GMyKinFitter_fNvar, GMyKinFitter_fNvar);
    Double_t    sqrtTerm    = TMath::Sqrt((pIn[3]*pIn[3])-(mass*mass));
    Double_t    r           = TMath::Sqrt((pIn[4]*pIn[4])+(pIn[5]*pIn[5])+(pIn[6]*pIn[6]));

    Double_t    help    = r*r*pIn[3]*pIn[0];
    ret[0][0]   = help*help*dpIn[3][3];
    help    = sqrtTerm*sqrtTerm*((pIn[1]*pIn[1])+(pIn[2]*pIn[2]));
    ret[0][0]   += help*help*dpIn[0][0];
    help    = sqrtTerm*sqrtTerm*pIn[1];
    ret[0][0]   += help*help*dpIn[1][1];
    help    = sqrtTerm*sqrtTerm*pIn[2];
    ret[0][0]   += help*help*dpIn[2][2];
    help    = sqrtTerm*r*r*r;
    ret[0][0]   /= help*help;

    help    = r*r*pIn[3];
    ret[1][0]   = help*help*pIn[0]*pIn[1]*dpIn[3][3];
    help    = sqrtTerm*sqrtTerm;
    ret[1][0]   += help*help*((pIn[1]*pIn[1])+(pIn[2]*pIn[2]))*pIn[0]*dpIn[0][0];
    help    = sqrtTerm*sqrtTerm;
    ret[1][0]   += help*help*((pIn[0]*pIn[0])+(pIn[2]*pIn[2]))*pIn[1]*dpIn[1][1];
    help    = sqrtTerm*sqrtTerm*pIn[2];
    ret[1][0]   += help*help*dpIn[2][2];
    help    = sqrtTerm*r*r*r;
    ret[1][0]   /= help*help;
    ret[0][1]   = ret[1][0];

    help    = r*r*pIn[3]*pIn[1];
    ret[1][1]   = help*help*dpIn[3][3];
    help    = sqrtTerm*sqrtTerm*pIn[0];
    ret[1][1]   += help*help*dpIn[0][0];
    help    = sqrtTerm*sqrtTerm*((pIn[0]*pIn[0])+(pIn[2]*pIn[2]));
    ret[1][1]   += help*help*dpIn[1][1];
    help    = sqrtTerm*sqrtTerm*pIn[2];
    ret[1][1]   += help*help*dpIn[2][2];
    help    = sqrtTerm*r*r*r;
    ret[1][1]   /= help*help;

    help    = r*r*pIn[3];
    ret[2][0]   = help*help*pIn[0]*pIn[2]*dpIn[3][3];
    help    = sqrtTerm*sqrtTerm;
    ret[2][0]   += help*help*((pIn[1]*pIn[1])+(pIn[2]*pIn[2]))*pIn[0]*dpIn[0][0];
    help    = sqrtTerm*sqrtTerm*pIn[1];
    ret[2][0]   += help*help*dpIn[1][1];
    help    = sqrtTerm*sqrtTerm;
    ret[2][0]   += help*help*((pIn[0]*pIn[0])+(pIn[1]*pIn[1]))*pIn[2]*dpIn[2][2];
    help    = sqrtTerm*r*r*r;
    ret[2][0]   /= help*help;
    ret[0][2]   = ret[2][0];

    help    = r*r*pIn[3];
    ret[2][1]   = help*help*pIn[1]*pIn[2]*dpIn[3][3];
    help    = sqrtTerm*sqrtTerm*pIn[0];
    ret[2][1]   += help*help*dpIn[0][0];
    help    = sqrtTerm*sqrtTerm;
    ret[2][1]   += help*help*((pIn[0]*pIn[0])+(pIn[2]*pIn[2]))*pIn[1]*dpIn[1][1];
    help    = sqrtTerm*sqrtTerm;
    ret[2][1]   += help*help*((pIn[0]*pIn[0])+(pIn[1]*pIn[1]))*pIn[2]*dpIn[2][2];
    help    = sqrtTerm*r*r*r;
    ret[2][1]   /= help*help;
    ret[1][2]   = ret[2][1];

    help    = r*r*pIn[3]*pIn[2];
    ret[2][2]   = help*help*dpIn[3][3];
    help    = sqrtTerm*sqrtTerm*pIn[0];
    ret[2][2]   += help*help*dpIn[0][0];
    help    = sqrtTerm*sqrtTerm*pIn[1];
    ret[2][2]   += help*help*dpIn[1][1];
    help    = sqrtTerm*sqrtTerm*((pIn[0]*pIn[0])+(pIn[1]*pIn[1]));
    ret[2][2]   += help*help*dpIn[2][2];
    help    = sqrtTerm*r*r*r;
    ret[2][2]   /= help*help;


    ret[3][0]   = r*r*pIn[3]*pIn[0]*dpIn[3][3];
    ret[3][0]   /= sqrtTerm*r*r*r;
    ret[0][3]   = ret[3][0];

    ret[3][1]   = r*r*pIn[3]*pIn[1]*dpIn[3][3];
    ret[3][1]   /= sqrtTerm*r*r*r;
    ret[1][3]   = ret[3][1];

    ret[3][2]   = r*r*pIn[3]*pIn[2]*dpIn[3][3];
    ret[3][2]   /= sqrtTerm*r*r*r;
    ret[2][3]   = ret[3][2];

    ret[3][3]   = dpIn[3][3];

/*
    ret[4][0]   = -dpIn[3]*dpIn[3]; ret[0][4]   = ret[4][0];
    ret[4][1]   = 0;                ret[1][4]   = 0;
    ret[4][2]   = 0;                ret[2][4]   = 0;
    ret[4][3]   = 0;                ret[3][4]   = 0;
    ret[4][4]   = 0;

    ret[5][0]   = 0;                ret[0][5]   = 0;
    ret[5][1]   = -dpIn[4]*dpIn[4]; ret[1][5]   = ret[5][1];
    ret[5][2]   = 0;                ret[2][5]   = 0;
    ret[5][3]   = 0;                ret[3][5]   = 0;
    ret[5][4]   = 0;                ret[4][5]   = 0;
    ret[5][5]   = 0;

    ret[6][0]   = 0;                ret[0][6]   = 0;
    ret[6][1]   = 0;                ret[1][6]   = 0;
    ret[6][2]   = -dpIn[5]*dpIn[5]; ret[2][6]   = ret[6][2];
    ret[6][3]   = 0;                ret[3][6]   = 0;
    ret[6][4]   = 0;                ret[4][6]   = 0;
    ret[6][5]   = 0;             ret[5][6]   = 0;
    ret[6][6]   = 0;
*/
    return ret;
}

void GMyKinFitter::AddParticle(const Double_t* p, const Double_t* dp)
{
    if(fNparti>nPart){
      std::cout<<"GKinFitter::AddPosKFParticle already at max particles"<<std::endl;
      return;
    }

    //Add parameters to Alpha0
    for(int i=0; i<GMyKinFitter_fNvar; i++)
        fmAlpha0[fNpari+i][0] = p[i];
    //Add error matrix to V_Alpha0
    for(int i=0; i<GMyKinFitter_fNvar; i++)
        fmAlpha0[fNpari+i][0] = p[i];
    //fmV_Alpha0.SetSub(fNpari,fNpari,p.GetVAlpha());
    //Add transformation matrice to fT
    //fT.SetSub(fNpari,fNpari,p.GetT());
    //ADD Lorentz Vector !!
    if(fNparti) fPtot=fPtot+p;
    else fPtot=p;

    //increment counters
    fNpari+=GMyKinFitter_fNvar;
    fNparti++;
}

void    GMyKinFitter::Reset()
{
    fNparti = 0;
    fNpari  = 0;
    fNconi  = 0;
    fNiter  = 0;
    fPtot.SetPxPyPzE(0, 0, 0, 0);
}


