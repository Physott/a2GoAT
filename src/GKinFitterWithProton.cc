//////////////////////////////////////////////////////////////////
//GKinFitter
//////////////////////////////////////////////////////////////////


#include "GKinFitterWithProton.h"
#include "TDecompLU.h"
#include <iostream>

/*

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

*/





GMyTrackWithProtonM::GMyTrackWithProtonM(const Double_t beamEnergy, const TLorentzVector& p0, const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& p3, const TLorentzVector& p4, const TLorentzVector& p5, const TLorentzVector& proton)
{
    p[0]    = 0;            p[1]    = 0;            p[2]    = 0;                //initialVertex
    p[3]    = 0;            p[4]    = 0;            p[5]    = 0;                //Distance initialVertex to eta Vertex
    p[6]    = 0;            p[7]    = 0;            p[8]    = 0;                //Distance initialVertex to pi0a Vertex
    p[9]    = 0;            p[10]   = 0;            p[11]   = 0;                //Distance initialVertex to pi0b Vertex

    p[12]   = beamEnergy;   p[13]   = 0;            p[14]   = 0;                //beam  [E, Theta, Phi]
    p[15]   = p0.E();       p[16]   = p0.Theta();   p[17]   = p0.Phi();         //gamma 0  [E, Theta, Phi]
    p[18]   = p1.E();       p[19]   = p1.Theta();   p[20]   = p1.Phi();         //gamma 0  [E, Theta, Phi]
    p[21]   = p2.E();       p[22]   = p2.Theta();   p[23]   = p2.Phi();         //gamma 0  [E, Theta, Phi]
    p[24]   = p3.E();       p[25]   = p3.Theta();   p[26]   = p3.Phi();         //gamma 0  [E, Theta, Phi]
    p[27]   = p4.E();       p[28]   = p4.Theta();   p[29]   = p4.Phi();         //gamma 0  [E, Theta, Phi]
    p[30]   = p5.E();       p[31]   = p5.Theta();   p[32]   = p5.Phi();         //gamma 0  [E, Theta, Phi]
    p[33]   = beamEnergy-p0.E()-p1.E()-p2.E()-p3.E()-p4.E()-p5.E();
    p[34]   = proton.Theta();p[35]   = proton.Phi();    //proton  [E, Theta, Phi]

    dp[0]    = 0.04;         dp[1]    = 0.04;         dp[2]    = 0.05;          //Target Abmessungen
    dp[3]    = 0.1;         dp[4]    = 0.1;         dp[5]    = 0.1;          //eta travel length
    dp[6]    = 0.1;         dp[7]    = 0.1;         dp[8]    = 0.1;          //pi0a travel length
    dp[9]    = 0.1;         dp[10]   = 0.1;         dp[11]   = 0.1;          //pi0b travel length

    dp[12]   = 1;       dp[13]   = 0.01;   dp[14]   = 0.01;
    dp[15]   = 5;      dp[16]   = 0.05;   dp[17]   = 0.05;
    dp[18]   = 5;      dp[19]   = 0.05;   dp[20]   = 0.05;
    dp[21]   = 5;      dp[22]   = 0.05;   dp[23]   = 0.05;
    dp[24]   = 5;      dp[25]   = 0.05;   dp[26]   = 0.05;
    dp[27]   = 5;      dp[28]   = 0.05;   dp[29]   = 0.05;
    dp[30]   = 5;      dp[31]   = 0.05;   dp[32]   = 0.05;
    dp[33]   = 5;      dp[34]   = 0.05;   dp[35]   = 0.05;
}

GMyTrackWithProtonM::~GMyTrackWithProtonM()
{
}

TMatrixD    GMyTrackWithProtonM::GetParametersH()    const
{
    TMatrixD    ret(44,1);
    ret[0][0]   = p[0];
    ret[1][0]   = p[1];
    ret[2][0]   = p[2];

    ret[3][0]   = p[0] + p[3];
    ret[4][0]   = p[1] + p[4];
    ret[5][0]   = p[2] + p[5];
    ret[6][0]   = p[0] + p[6];
    ret[7][0]   = p[1] + p[7];
    ret[8][0]   = p[2] + p[8];
    ret[9][0]   = p[0] + p[9];
    ret[10][0]  = p[1] + p[10];
    ret[11][0]  = p[2] + p[11];

    ret[12][0]  = GMyTrackM_CBRadius * TMath::Sin(p[13]) * TMath::Cos(p[14]) - p[0];
    ret[13][0]  = GMyTrackM_CBRadius * TMath::Sin(p[13]) * TMath::Sin(p[14]) - p[1];
    ret[14][0]  = GMyTrackM_CBRadius * TMath::Cos(p[13])                     - p[2];
    ret[15][0]  = p[12];

    ret[16][0]  = GMyTrackM_CBRadius * TMath::Sin(p[16]) * TMath::Cos(p[17]) - p[0] - p[3];
    ret[17][0]  = GMyTrackM_CBRadius * TMath::Sin(p[16]) * TMath::Sin(p[17]) - p[1] - p[4];
    ret[18][0]  = GMyTrackM_CBRadius * TMath::Cos(p[16])                     - p[2] - p[5];
    ret[19][0]  = p[15];

    ret[20][0]  = GMyTrackM_CBRadius * TMath::Sin(p[19]) * TMath::Cos(p[20]) - p[0] - p[3];
    ret[21][0]  = GMyTrackM_CBRadius * TMath::Sin(p[19]) * TMath::Sin(p[20]) - p[1] - p[4];
    ret[22][0]  = GMyTrackM_CBRadius * TMath::Cos(p[19])                     - p[2] - p[5];
    ret[23][0]  = p[18];

    ret[24][0]  = GMyTrackM_CBRadius * TMath::Sin(p[22]) * TMath::Cos(p[23]) - p[0] - p[6];
    ret[25][0]  = GMyTrackM_CBRadius * TMath::Sin(p[22]) * TMath::Sin(p[23]) - p[1] - p[7];
    ret[26][0]  = GMyTrackM_CBRadius * TMath::Cos(p[22])                     - p[2] - p[8];
    ret[27][0]  = p[21];

    ret[28][0]  = GMyTrackM_CBRadius * TMath::Sin(p[25]) * TMath::Cos(p[26]) - p[0] - p[6];
    ret[29][0]  = GMyTrackM_CBRadius * TMath::Sin(p[25]) * TMath::Sin(p[26]) - p[1] - p[7];
    ret[30][0]  = GMyTrackM_CBRadius * TMath::Cos(p[25])                     - p[2] - p[8];
    ret[31][0]  = p[24];

    ret[32][0]  = GMyTrackM_CBRadius * TMath::Sin(p[28]) * TMath::Cos(p[29]) - p[0] - p[9];
    ret[33][0]  = GMyTrackM_CBRadius * TMath::Sin(p[28]) * TMath::Sin(p[29]) - p[1] - p[10];
    ret[34][0]  = GMyTrackM_CBRadius * TMath::Cos(p[28])                     - p[2] - p[11];
    ret[35][0]  = p[27];

    ret[36][0]  = GMyTrackM_CBRadius * TMath::Sin(p[31]) * TMath::Cos(p[32]) - p[0] - p[9];
    ret[37][0]  = GMyTrackM_CBRadius * TMath::Sin(p[31]) * TMath::Sin(p[32]) - p[1] - p[10];
    ret[38][0]  = GMyTrackM_CBRadius * TMath::Cos(p[31])                     - p[2] - p[11];
    ret[39][0]  = p[30];

    ret[40][0]  = GMyTrackM_CBRadius * TMath::Sin(p[34]) * TMath::Cos(p[35]) - p[0];
    ret[41][0]  = GMyTrackM_CBRadius * TMath::Sin(p[34]) * TMath::Sin(p[35]) - p[1];
    ret[42][0]  = GMyTrackM_CBRadius * TMath::Cos(p[34])                     - p[2];
    ret[43][0]  = p[33];

    return ret;
}

TMatrixD    GMyTrackWithProtonM::GetDerivatedParametersH()    const
{
    TMatrixD    ret(44,36);
    for(int i=0; i<44; i++)
    {
        for(int j=0; j<36; j++)
            ret[i][j]   = 0;
    }

    ret[0][0]   = 1;
    ret[1][1]   = 1;
    ret[2][2]   = 1;

    ret[3][0]   = 1;
    ret[3][3]   = 1;
    ret[4][1]   = 1;
    ret[4][4]   = 1;
    ret[5][2]   = 1;
    ret[5][5]   = 1;
    ret[6][0]   = 1;
    ret[6][6]   = 1;
    ret[7][1]   = 1;
    ret[7][7]   = 1;
    ret[8][2]   = 1;
    ret[8][8]   = 1;
    ret[9][0]   = 1;
    ret[9][9]   = 1;
    ret[10][1]  = 1;
    ret[10][10] = 1;
    ret[11][2]  = 1;
    ret[11][11] = 1;

    ret[12][0]  = -1;
    ret[12][13] = GMyTrackM_CBRadius * TMath::Cos(p[13]) * TMath::Cos(p[14]);
    ret[12][14] = -GMyTrackM_CBRadius * TMath::Sin(p[13]) * TMath::Sin(p[14]);
    ret[13][1]  = -1;
    ret[13][13] = GMyTrackM_CBRadius * TMath::Cos(p[13]) * TMath::Sin(p[14]);
    ret[13][14] = GMyTrackM_CBRadius * TMath::Sin(p[13]) * TMath::Cos(p[14]);
    ret[14][2]  = -1;
    ret[14][13] = -GMyTrackM_CBRadius * TMath::Sin(p[13]);
    ret[15][12] = 1;

    ret[16][0]  = -1;
    ret[16][3]  = -1;
    ret[16][16] = GMyTrackM_CBRadius * TMath::Cos(p[16]) * TMath::Cos(p[17]);
    ret[16][17] = -GMyTrackM_CBRadius * TMath::Sin(p[16]) * TMath::Sin(p[17]);
    ret[17][1]  = -1;
    ret[17][4]  = -1;
    ret[17][16] = GMyTrackM_CBRadius * TMath::Cos(p[16]) * TMath::Sin(p[17]);
    ret[17][17] = GMyTrackM_CBRadius * TMath::Sin(p[16]) * TMath::Cos(p[17]);
    ret[18][2]  = -1;
    ret[18][5]  = -1;
    ret[18][16] = -GMyTrackM_CBRadius * TMath::Sin(p[16]);
    ret[19][15] = 1;

    ret[20][0]  = -1;
    ret[20][3]  = -1;
    ret[20][19] = GMyTrackM_CBRadius * TMath::Cos(p[19]) * TMath::Cos(p[20]);
    ret[20][20] = -GMyTrackM_CBRadius * TMath::Sin(p[19]) * TMath::Sin(p[20]);
    ret[21][1]  = -1;
    ret[21][4]  = -1;
    ret[21][19] = GMyTrackM_CBRadius * TMath::Cos(p[19]) * TMath::Sin(p[20]);
    ret[21][20] = GMyTrackM_CBRadius * TMath::Sin(p[19]) * TMath::Cos(p[20]);
    ret[22][2]  = -1;
    ret[22][5]  = -1;
    ret[22][19] = -GMyTrackM_CBRadius * TMath::Sin(p[19]);
    ret[23][18] = 1;

    ret[24][0]  = -1;
    ret[24][6]  = -1;
    ret[24][22] = GMyTrackM_CBRadius * TMath::Cos(p[22]) * TMath::Cos(p[23]);
    ret[24][23] = -GMyTrackM_CBRadius * TMath::Sin(p[22]) * TMath::Sin(p[23]);
    ret[25][1]  = -1;
    ret[25][7]  = -1;
    ret[25][22] = GMyTrackM_CBRadius * TMath::Cos(p[22]) * TMath::Sin(p[23]);
    ret[25][23] = GMyTrackM_CBRadius * TMath::Sin(p[22]) * TMath::Cos(p[23]);
    ret[26][2]  = -1;
    ret[26][8]  = -1;
    ret[26][22] = -GMyTrackM_CBRadius * TMath::Sin(p[22]);
    ret[27][21] = 1;

    ret[28][0]  = -1;
    ret[28][6]  = -1;
    ret[28][25] = GMyTrackM_CBRadius * TMath::Cos(p[25]) * TMath::Cos(p[26]);
    ret[28][26] = -GMyTrackM_CBRadius * TMath::Sin(p[25]) * TMath::Sin(p[26]);
    ret[29][1]  = -1;
    ret[29][7]  = -1;
    ret[29][25] = GMyTrackM_CBRadius * TMath::Cos(p[25]) * TMath::Sin(p[26]);
    ret[29][26] = GMyTrackM_CBRadius * TMath::Sin(p[25]) * TMath::Cos(p[26]);
    ret[30][2]  = -1;
    ret[30][8]  = -1;
    ret[30][25] = -GMyTrackM_CBRadius * TMath::Sin(p[25]);
    ret[31][24] = 1;

    ret[32][0]  = -1;
    ret[32][9]  = -1;
    ret[32][28] = GMyTrackM_CBRadius * TMath::Cos(p[28]) * TMath::Cos(p[29]);
    ret[32][29] = -GMyTrackM_CBRadius * TMath::Sin(p[28]) * TMath::Sin(p[29]);
    ret[33][1]  = -1;
    ret[33][10] = -1;
    ret[33][28] = GMyTrackM_CBRadius * TMath::Cos(p[28]) * TMath::Sin(p[29]);
    ret[33][29] = GMyTrackM_CBRadius * TMath::Sin(p[28]) * TMath::Cos(p[29]);
    ret[34][2]  = -1;
    ret[34][11] = -1;
    ret[34][28] = -GMyTrackM_CBRadius * TMath::Sin(p[28]);
    ret[35][27] = 1;

    ret[36][0]  = -1;
    ret[36][9]  = -1;
    ret[36][31] = GMyTrackM_CBRadius * TMath::Cos(p[31]) * TMath::Cos(p[32]);
    ret[36][32] = -GMyTrackM_CBRadius * TMath::Sin(p[31]) * TMath::Sin(p[32]);
    ret[37][1]  = -1;
    ret[37][10] = -1;
    ret[37][31] = GMyTrackM_CBRadius * TMath::Cos(p[31]) * TMath::Sin(p[32]);
    ret[37][32] = GMyTrackM_CBRadius * TMath::Sin(p[31]) * TMath::Cos(p[32]);
    ret[38][2]  = -1;
    ret[38][11] = -1;
    ret[38][31] = -GMyTrackM_CBRadius * TMath::Sin(p[31]);
    ret[39][30] = 1;

    ret[40][0]  = -1;
    ret[40][34] = GMyTrackM_CBRadius * TMath::Cos(p[34]) * TMath::Cos(p[35]);
    ret[40][35] = -GMyTrackM_CBRadius * TMath::Sin(p[34]) * TMath::Sin(p[35]);
    ret[41][1]  = -1;
    ret[41][34] = GMyTrackM_CBRadius * TMath::Cos(p[34]) * TMath::Sin(p[35]);
    ret[41][35] = GMyTrackM_CBRadius * TMath::Sin(p[34]) * TMath::Cos(p[35]);
    ret[42][2]  = -1;
    ret[42][34] = -GMyTrackM_CBRadius * TMath::Sin(p[34]);
    ret[43][33] = 1;

    return  ret;
}

TMatrixD    GMyTrackWithProtonM::GetCovarianceH()    const
{
    TMatrixD    cm(36, 36);
    for(int i=0; i<36; i++)
    {
        for(int j=0; j<36; j++)
        {
            if(i==j)
                cm[i][j]   = dp[i] * dp[i];
            else
                cm[i][j]   = 0;
        }
    }

    TMatrixD    dp(GetDerivatedParametersH());
    TMatrixD    dpT(dp);
    dpT.T();

    return dp * cm * dpT;
}




GMyTrackWithProtonH::GMyTrackWithProtonH(const GMyTrackWithProtonM &m):
    cm(44,44)
{
    TMatrixD    help(m.GetParametersH());

    for(int i=0; i<44; i++)
        p[i]    = help[i][0];

    cm  = m.GetCovarianceH();
}

GMyTrackWithProtonH::~GMyTrackWithProtonH()
{

}

#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046

TMatrixD    GMyTrackWithProtonH::GetParametersW()    const
{
    TMatrixD    ret(44,1);

    for(int i=0; i<12; i++)
        ret[i][0]  = p[i];

    Double_t    help;
    help    = p[15]-MASS_PROTON;
    help    /= TMath::Sqrt((p[12]*p[12])+(p[13]*p[13])+(p[14]*p[14]));
    ret[12][0]   = help * p[12];
    ret[13][0]   = help * p[13];
    ret[14][0]   = help * p[14];
    ret[15][0]   = p[15];

    for(int i=0; i<6; i++)
    {
        help    = p[4*i+19];
        help    /= TMath::Sqrt((p[4*i+16]*p[4*i+16])+(p[4*i+17]*p[4*i+17])+(p[4*i+18]*p[4*i+18]));
        ret[4*i+16][0]   = help * p[4*i+16];
        ret[4*i+17][0]   = help * p[4*i+17];
        ret[4*i+18][0]   = help * p[4*i+18];
        ret[4*i+19][0]   = p[4*i+19];
    }

    help    = TMath::Sqrt((p[43]*p[43])-(MASS_PROTON*MASS_PROTON));
    help    /= TMath::Sqrt((p[40]*p[40])+(p[41]*p[41])+(p[42]*p[42]));
    ret[40][0]   = help * p[40];
    ret[41][0]   = help * p[13];
    ret[42][0]   = help * p[42];
    ret[43][0]   = p[43];

    return ret;
}

TMatrixD    GMyTrackWithProtonH::GetDerivatedParametersW()    const
{
    TMatrixD    ret(44,44);
    for(int i=0; i<44; i++)
    {
        for(int j=0; j<44; j++)
            ret[i][j]   = 0;
    }

    for(int i=0; i<12; i++)
        ret[i][i]  = 1;

    Double_t    help;
    Double_t    r;
    help    = TMath::Sqrt((p[15]*p[15])-(MASS_PROTON*MASS_PROTON));
    r       = TMath::Sqrt((p[12]*p[12])+(p[13]*p[13])+(p[14]*p[14]));
    ret[12][12]   = help * ((p[13]*p[13])+(p[14]*p[14])) /(r*r*r);
    ret[12][13]   = -help * p[12] * p[13] /(r*r*r);
    ret[12][14]   = -help * p[12] * p[14] /(r*r*r);
    ret[12][15]   = p[15] * p[12] /(help*r);
    ret[13][12]   = -help * p[13] * p[12] /(r*r*r);
    ret[13][13]   = help * ((p[12]*p[12])+(p[14]*p[14])) /(r*r*r);
    ret[13][14]   = -help * p[13] * p[14] /(r*r*r);
    ret[13][15]   = p[15] * p[13] /(help*r);
    ret[14][12]   = -help * p[14] * p[12] /(r*r*r);
    ret[14][13]   = -help * p[14] * p[13] /(r*r*r);
    ret[14][14]   = help * ((p[12]*p[12])+(p[13]*p[13])) /(r*r*r);
    ret[14][15]   = p[15] * p[14] /(help*r);
    ret[15][15]   = 1;

    for(int i=0; i<6; i++)
    {
        r       = TMath::Sqrt((p[4*i+16]*p[4*i+16])+(p[4*i+17]*p[4*i+17])+(p[4*i+18]*p[4*i+18]));
        ret[4*i+16][4*i+16]   = p[4*i+19] * ((p[4*i+17]*p[4*i+17])+(p[4*i+18]*p[4*i+18])) /(r*r*r);
        ret[4*i+16][4*i+17]   = -p[4*i+19] * p[4*i+16] * p[4*i+17] /(r*r*r);
        ret[4*i+16][4*i+18]   = -p[4*i+19] * p[4*i+16] * p[4*i+18] /(r*r*r);
        ret[4*i+16][4*i+19]   = p[4*i+16] /r;
        ret[4*i+17][4*i+16]   = -p[4*i+19] * p[4*i+17] * p[4*i+16] /(r*r*r);
        ret[4*i+17][4*i+17]   = p[4*i+19] * ((p[4*i+16]*p[4*i+16])+(p[4*i+18]*p[4*i+18])) /(r*r*r);
        ret[4*i+17][4*i+18]   = -p[4*i+19] * p[4*i+17] * p[4*i+18] /(r*r*r);
        ret[4*i+17][4*i+19]   = p[4*i+17] /r;
        ret[4*i+18][4*i+16]   = -p[4*i+19] * p[4*i+18] * p[4*i+16] /(r*r*r);
        ret[4*i+18][4*i+17]   = -p[4*i+19] * p[4*i+18] * p[4*i+17] /(r*r*r);
        ret[4*i+18][4*i+18]   = p[4*i+19] * ((p[4*i+16]*p[4*i+16])+(p[4*i+17]*p[4*i+17])) /(r*r*r);
        ret[4*i+18][4*i+19]   = p[4*i+18] /r;
        ret[4*i+19][4*i+19]   = 1;
    }

    help    = TMath::Sqrt((p[43]*p[43])-(MASS_PROTON*MASS_PROTON));
    r       = TMath::Sqrt((p[40]*p[40])+(p[41]*p[41])+(p[42]*p[42]));
    ret[40][40]   = help * ((p[41]*p[41])+(p[42]*p[42])) /(r*r*r);
    ret[40][41]   = -help * p[40] * p[41] /(r*r*r);
    ret[40][42]   = -help * p[40] * p[42] /(r*r*r);
    ret[40][43]   = p[43] * p[40] /(help*r);
    ret[41][40]   = -help * p[41] * p[40] /(r*r*r);
    ret[41][41]   = help * ((p[40]*p[40])+(p[42]*p[42])) /(r*r*r);
    ret[41][42]   = -help * p[41] * p[42] /(r*r*r);
    ret[41][43]   = p[43] * p[41] /(help*r);
    ret[42][40]   = -help * p[42] * p[40] /(r*r*r);
    ret[42][41]   = -help * p[42] * p[41] /(r*r*r);
    ret[42][42]   = help * ((p[40]*p[40])+(p[41]*p[41])) /(r*r*r);
    ret[42][43]   = p[43] * p[42] /(help*r);
    ret[43][43]   = 1;

    return ret;
}

TMatrixD    GMyTrackWithProtonH::GetCovarianceW()    const
{
    TMatrixD    dp(GetDerivatedParametersW());
    TMatrixD    dpT(dp);
    dpT.T();

    return dp * cm * dpT;
}










GKinFitterWithProton::GKinFitterWithProton()    :
    nPar(44),
    nCon(7),
    fNiter(0),
    fmAlpha0(nPar,1),
    fmAlpha(nPar,1),
    fmV_Alpha0(nPar,nPar),
    fmV_Alpha(nPar,nPar),
    fmD(nCon,nPar),
    fmd(nCon,1),
    fmlamda(nCon,1),
    fmV_D(nCon,nCon),
    fchi2(0),
    fPtot(0, 0, 0, 0),
    solved(kFALSE)
{
}

GKinFitterWithProton::~GKinFitterWithProton()
{

}

void GKinFitterWithProton::Constraints()
{
  //Add up the particle 4 vectors
    TLorentzVector eta(fmAlpha0[16][0], fmAlpha0[17][0], fmAlpha0[18][0], fmAlpha0[19][0]);
    eta     += TLorentzVector(fmAlpha0[20][0], fmAlpha0[21][0], fmAlpha0[22][0], fmAlpha0[23][0]);
    TLorentzVector pi0a(fmAlpha0[24][0], fmAlpha0[25][0], fmAlpha0[26][0], fmAlpha0[27][0]);
    pi0a    += TLorentzVector(fmAlpha0[28][0], fmAlpha0[29][0], fmAlpha0[30][0], fmAlpha0[31][0]);
    TLorentzVector pi0b(fmAlpha0[32][0], fmAlpha0[33][0], fmAlpha0[34][0], fmAlpha0[35][0]);
    pi0b    += TLorentzVector(fmAlpha0[36][0], fmAlpha0[37][0], fmAlpha0[38][0], fmAlpha0[39][0]);
    TLorentzVector all(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]);
    all     -= TLorentzVector(fmAlpha0[40][0], fmAlpha0[41][0], fmAlpha0[42][0], fmAlpha0[43][0]);
    all     -= eta;
    all     -= pi0a;
    all     -= pi0b;
    //std::cout << all.Px() << "   " << all.Py() << "   " << all.Pz() << "   " << all.E() << "   " << std::endl;
    //TLorentzVector(fmAlpha0[40][0], fmAlpha0[41][0], fmAlpha0[42][0], fmAlpha0[43][0]).Print();
    //TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]).Print();

  //d matrix (evaluate constraint eqn.)
    fmd[0][0]   =   eta.M2() - (MASS_ETA*MASS_ETA);
    fmd[1][0]   =   pi0a.M2() - (MASS_PI0*MASS_PI0);
    fmd[2][0]   =   pi0b.M2() - (MASS_PI0*MASS_PI0);
    fmd[3][0]   =   all.Px();
    fmd[4][0]   =   all.Py();
    fmd[5][0]   =   all.Pz();
    fmd[6][0]   =   all.E();

  //D matrix (derivitives of constraint eqn)
    for(int i=0; i<7; i++)
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

    fmD[3][12] = 1;
    fmD[4][13] = 1;
    fmD[5][14] = 1;
    fmD[6][15] = 1;
    for(int i=0; i<7; i++)
    {
        fmD[3][4*i+16] = -1;
        fmD[4][4*i+17] = -1;
        fmD[5][4*i+18] = -1;
        fmD[6][4*i+19] = -1;
    }
}

TLorentzVector  GKinFitterWithProton::GetEta()
{
    TLorentzVector ret(fmAlpha[16][0], fmAlpha[17][0], fmAlpha[18][0], fmAlpha[19][0]);
        ret +=  TLorentzVector(fmAlpha[20][0], fmAlpha[21][0], fmAlpha[22][0], fmAlpha[23][0]);
    return ret;
}

TLorentzVector  GKinFitterWithProton::GetEtap()
{
    TLorentzVector ret(0.0, 0.0, 0.0, 0.0);
    for(int i=0; i<6; i++)
        ret +=  TLorentzVector(fmAlpha[4*i+16][0], fmAlpha[4*i+17][0], fmAlpha[4*i+18][0], fmAlpha[4*i+19][0]);
    return ret;
}

TLorentzVector  GKinFitterWithProton::GetPi0a()
{
    TLorentzVector ret(fmAlpha[24][0], fmAlpha[25][0], fmAlpha[26][0], fmAlpha[27][0]);
        ret +=  TLorentzVector(fmAlpha[28][0], fmAlpha[29][0], fmAlpha[30][0], fmAlpha[31][0]);
    return ret;
}

TLorentzVector  GKinFitterWithProton::GetPi0b()
{
    TLorentzVector ret(fmAlpha[32][0], fmAlpha[33][0], fmAlpha[34][0], fmAlpha[35][0]);
        ret +=  TLorentzVector(fmAlpha[36][0], fmAlpha[37][0], fmAlpha[38][0], fmAlpha[39][0]);
    return ret;
}

void    GKinFitterWithProton::Set(const Double_t beamEnergy, const TLorentzVector& p0, const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& p3, const TLorentzVector& p4, const TLorentzVector& p5, const TLorentzVector& proton)
{
    GMyTrackWithProtonH   parameters(GMyTrackWithProtonM(beamEnergy, p0, p1, p2, p3, p4, p5, proton));

    fmAlpha0    = parameters.GetParametersW();
    fmV_Alpha0  = parameters.GetCovarianceW();

    Constraints();

    solved  = kFALSE;
}

Bool_t  GKinFitterWithProton::Solve()
{

  //Solve according to algorithm of Paul Avery:
  //Applied Fitting Theory VI, Formulas for Kinematic Fitting
  //see www.phys.ufl.edu/~avery/fitting.html


  //std::cout << fmD.GetNrows() << "   " << fmD.GetNcols() << std::endl;
  //fmD.Print();
  TMatrixD mDT=fmD;
  mDT.T();
  //std::cout << mDT.GetNrows() << "   " << mDT.GetNcols() << "   " << mDT.Determinant() << std::endl;
  TMatrixD mV_Dinv=fmD*fmV_Alpha0*mDT;
  //std::cout << mV_Dinv.GetNrows() << "   " << mV_Dinv.GetNcols() << "   " << mV_Dinv.Determinant() << std::endl;
  fmV_D=mV_Dinv;
  /*TDecompLU lu(fmV_D);
  if(!lu.Decompose()){
    std::cout<<"GKinFitter::Solve() Cannot invert. KinFit not completed"<<std::endl;
    return kFALSE;
  }*/
  //std::cout << fmV_D.GetNrows() << "   " << fmV_D.GetNcols() << std::endl;
  //fmV_D.Print();
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

  solved  = kTRUE;
  return kTRUE;

}
