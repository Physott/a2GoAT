#ifndef _GKinFitter_h
#define _GKinFitter_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "GKinFitterParticle.h"

class GKinIterativeFitter;

class GKinFitter 
{
	private:
	Int_t fNvar; //Number of variables (=4 for 4 vector)
	Int_t fNpar; //Number of parameters=Npart*fNvar
	Int_t fNpart; //Number of particles
	Int_t fNcon; //Number of constraints
	Int_t fNparti; //count Number of particles added
	Int_t fNpari; //count Number of particles added
	Int_t fNconi; // countNumber of constraints added
	Int_t fNiter; // Number of times Solve has been called
	Int_t fNunKnown; // Number of unknowns
	TMatrixD fmAlpha0;	//original parameters
	TMatrixD fmAlpha;  	//fitted parameters
	TMatrixD fmV_Alpha0;//Covariance matrix for original parameters
	TMatrixD fmV_Alpha; //Covariance matrix for fitted parameters
	TMatrixD fmD;      	//Matrix of constraint derivitives
	TMatrixD fmd;      	//Vector of evaluated constraints
	TMatrixD fmlamda;  	//Vector of lagrangian multipliers
	TMatrixD fmV_D;    	//Covariance matrix of constraints (TO BE INVERTED)
	Double_t fchi2;
	TMatrixD fT;     	//Overall transforamtin matrix from Spherical->Cart
	TLorentzVector fPtot;

    Int_t SolveStep(const TMatrixD mDT); //do the least squares fit

public:
    GKinFitter(const Int_t npart, const Int_t ncon, const Int_t unk);
    virtual ~GKinFitter()							{}

    Int_t Solve(); //do the least squares fit

	//Form the D and d matrixes for the fit
    void AddInvMassConstraint(const Double_t Minv);   //based on Invariant mass of added particles
    void AddSubInvMassConstraint(const Int_t Np, const Int_t pid[], const Double_t Minv);//Add invariant mass constraint to subset of particles
    void AddTotEnergyConstraint(const Double_t Etot);  //based on total energy of added particles
    void AddTotMomentumConstraint(const TVector3 mom); //based on total 3 momentum of added particles
    void AddSubMissMassConstraint(const TLorentzVector Mom, const Int_t Np, const Int_t pid[], const Double_t MissMass); // Based on missing mass of particles in subset

    void AddPosKFParticle(const GKinFitterParticle kfp);//Add GKinFitterParticles to be fitted contribute + to fPtot
    void AddNegKFParticle(const GKinFitterParticle kfp);//Add GKinFitterParticles to be fitted contribute - to fPtot

	GKinFitterParticle GetTotalFitParticle();//returns alpha=fPtot and sum of error matrices  from each particle
    TLorentzVector Get4Vector()                             {return fPtot;}
    GKinFitterParticle GetParticle(const Int_t ip);
    GKinFitterParticle GetInitialParticle(const Int_t ip);
    Double_t GetChi2()                                      {return fchi2;}

    Double_t ConfidenceLevel()      {return TMath::Prob(fchi2,fNcon-fNunKnown);}//Note should be Ncon-Nunknowns
    Double_t Pull(const Int_t i)    {return (fmAlpha0[i][0]-fmAlpha[i][0])/sqrt(fmV_Alpha0[i][i]-fmV_Alpha[i][i]);}

    void ResetConstraints() {fNconi=0;}
    void ResetParticles()   {fNpari=0;fNparti=0;fPtot.SetXYZT(0,0,0,0);}
	void ResetMatrices();
    void Reset()            {ResetConstraints();ResetParticles();ResetMatrices();}
	void Debug();

    friend class GKinIterativeFitter;
};



class GKinIterativeFitter
{
private:
    GKinFitter          fitter;
    Int_t               nIter; // Number of times Solve has been called
    Int_t               nPar;
    GKinFitterParticle  par[20];
    Bool_t              parSign[20];
    TMatrixD            fmAlpha0;	//original parameters
    TMatrixD            fmV_Alpha0;//Covariance matrix for original parameters
    TMatrixD            result;	//fitted parameters
    TMatrixD            dresult;//Covariance matrix for fitted parameters
    Int_t               nIM;
    Double_t            IM[10];
    Int_t               IM_Np[10];
    Int_t*              IM_pid[10];
    Int_t               nMM;
    Double_t            MM[10];
    TLorentzVector      MM_mom[10];
    Int_t               nTE;
    Double_t            TE[10];
    Int_t               nTM;
    TVector3            TM[10];

    Int_t Init();
    Int_t Init(const GKinFitterParticle *p);
    Int_t SolveStep();

public:
    GKinIterativeFitter(const Int_t npart, const Int_t ncon, const Int_t unk);
    ~GKinIterativeFitter();

    Int_t Solve(); //do the least squares fit

    //Form the D and d matrixes for the fit
    void AddSubInvMassConstraint(const Int_t Np, const Int_t pid[], const Double_t Minv);
    void AddTotEnergyConstraint(const Double_t Etot)                                        {TE[nTE] = Etot; nTE++;}
    void AddTotMomentumConstraint(const TVector3 mom)                                       {TM[nTM] = mom; nTM++;}
    void AddSubMissMassConstraint(const TLorentzVector Mom, const Double_t MissMass)        {MM[nMM] = MissMass; MM_mom[nMM] = Mom; nMM++;}

    void AddPosKFParticle(const GKinFitterParticle kfp)                                     {fmAlpha0.SetSub(nPar,0,kfp.GetAlpha()); fmV_Alpha0.SetSub(nPar,nPar,kfp.GetVAlpha()); par[nPar].Set4Vector(kfp.Get4Vector()); par[nPar].SetVAlpha(kfp.GetVAlpha()); parSign[nPar]=kTRUE; nPar++;}
    void AddNegKFParticle(const GKinFitterParticle kfp)                                     {fmAlpha0.SetSub(nPar,0,kfp.GetAlpha()); fmV_Alpha0.SetSub(nPar,nPar,kfp.GetVAlpha());par[nPar].Set4Vector(kfp.Get4Vector()); par[nPar].SetVAlpha(kfp.GetVAlpha()); parSign[nPar]=kFALSE; nPar++;}

    Int_t               GetNIter()                          {return nIter;}
    GKinFitterParticle GetTotalFitParticle()                {return fitter.GetTotalFitParticle();}
    TLorentzVector Get4Vector()                             {return fitter.Get4Vector();}
    GKinFitterParticle GetParticle(const Int_t ip)          {return fitter.GetParticle(ip);}
    GKinFitterParticle GetInitialParticle(const Int_t ip)   {return fitter.GetInitialParticle(ip);}
    Double_t GetChi2()                                      {return fitter.GetChi2();}

    Double_t ConfidenceLevel()      {return fitter.ConfidenceLevel();}//Note should be Ncon-Nunknowns
    Double_t Pull(const Int_t i)    {return (fmAlpha0[i][0]-fitter.fmAlpha[i][0])/sqrt(fmV_Alpha0[i][i]-fitter.fmV_Alpha[i][i]);}

    void Reset()            {nPar=0; nIM=0; nMM=0; nTE=0; nTM=0; fitter.Reset();}
};


#endif

