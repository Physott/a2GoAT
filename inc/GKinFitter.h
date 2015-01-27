#ifndef _GKinFitter_h
#define _GKinFitter_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "iostream"

#include "GKinFitterParticle.h"

class GKinFitter 
{
protected:
    Int_t fNvar; //Number of variables (=4 for 4 vector)
    Int_t fNpar; //Number of parameters=Npart*fNvar
    Int_t fNpart; //Number of particles
    Int_t fNcon; //Number of constraints
	Int_t fNparti; //count Number of particles added
	Int_t fNpari; //count Number of particles added
    Int_t fNconi; // countNumber of constraints added
	TMatrixD fmAlpha0;	//original parameters
    TMatrixD fmAlpha1; 	//initial parameters
    TMatrixD fmAlpha2;	//fitted parameters
	TMatrixD fmV_Alpha0;//Covariance matrix for original parameters
	TMatrixD fmV_Alpha; //Covariance matrix for fitted parameters
    Double_t fmUnk0; 	//initial parameters
    Double_t fmUnk; 	//fitted parameters
	TMatrixD fmD;      	//Matrix of constraint derivitives
    TMatrixD fmDUnk;  	//Matrix of constraint derivitives
	TMatrixD fmd;      	//Vector of evaluated constraints
	TMatrixD fmlamda;  	//Vector of lagrangian multipliers
	TMatrixD fmV_D;    	//Covariance matrix of constraints (TO BE INVERTED)
    Double_t Cchi2;
    Double_t Vchi2;

    Double_t       GetParticleFactor(const Int_t ip);
    Double_t       GetInitialParticleFactor(const Int_t ip);
    Double_t       GetInitialParticleDerivateEnergyFactor(const Int_t ip);
    Double_t       GetInitialParticleDerivatePzFactor(const Int_t ip);
    Double_t       GetInitialParticleDerivateUnkFactor(const Int_t ip);

public:
    GKinFitter(const Int_t npart, const Int_t ncon);
    virtual ~GKinFitter()							{}

    virtual Int_t Solve(); //do the least squares fit

	//Form the D and d matrixes for the fit
    virtual void AddInvMassConstraint(const Double_t Minv);   //based on Invariant mass of added particles
    virtual void AddSubInvMassConstraint(const Int_t Np, const Int_t pid[], const Double_t Minv);//Add invariant mass constraint to subset of particles
    virtual void AddTotEnergyConstraint(const Double_t Etot);  //based on total energy of added particles
    virtual void AddTotMomentumConstraint(const TVector3 mom); //based on total 3 momentum of added particles
    virtual void AddSubMissMassConstraint(const TLorentzVector Mom, const Int_t Np, const Int_t pid[], const Double_t MissMass); // Based on missing mass of particles in subset

    void AddPosKFParticle(const GKinFitterParticle kfp);//Add GKinFitterParticles to be fitted
    void AddNegKFParticle(const GKinFitterParticle kfp);//Add GKinFitterParticles to be fitted

    TLorentzVector GetTotalFitParticle();//returns alpha=fPtot and sum of error matrices  from each particle
    TLorentzVector GetParticle(const Int_t ip);
    Double_t       GetUnknown()                             {return fmUnk;}
    TLorentzVector GetInitialParticle(const Int_t ip);
    Double_t       GetInitialUnknown()                      {return fmUnk0;}
    TLorentzVector GetOriginalParticle(const Int_t ip);
    Double_t GetConstraintsChi2()                           {return Cchi2;}
    Double_t GetVariablesChi2()                             {return Vchi2;}
    Double_t GetChi2()                                      {return Cchi2+Vchi2;}

    Double_t ConfidenceLevel()              {return TMath::Prob(GetChi2(),fNpar+fNcon);}//Note should be Ncon-Nunknowns
    Double_t ConstraintsConfidenceLevel()   {return TMath::Prob(GetConstraintsChi2(),fNcon);}//Note should be Ncon-Nunknowns
    Double_t VariablesConfidenceLevel()     {return TMath::Prob(GetVariablesChi2(),fNpar);}//Note should be Ncon-Nunknowns
    Double_t Pull(const Int_t i)    {return (fmAlpha0[i][0]-fmAlpha2[i][0])/sqrt(fmV_Alpha0[i][i]-fmV_Alpha[i][i]);}

    void ResetConstraints() {fNconi=0;}
    void ResetParticles()   {fNpari=0;fNparti=0;}
	void ResetMatrices();
    void Reset()            {ResetConstraints();ResetParticles();ResetMatrices();}
	void Debug();
};






class GIterativeKinFitter    : public GKinFitter
{
public:
    enum    ConstraintType
    {
        ConstraintType_InvMass,
        ConstraintType_SubInvMass,
        ConstraintType_MisMass,
        ConstraintType_InvEnergy,
        ConstraintType_InvMomentum
    };

private:
    Int_t           nIter;
    ConstraintType  conType[20];
    Double_t        var[20];
    Int_t           nIndices[20];
    Int_t           indices[20][20];
    TVector3        momentum[20];
    TLorentzVector  beam[20];

public:
    GIterativeKinFitter(const Int_t npart, const Int_t ncon);
    ~GIterativeKinFitter();

    virtual void AddInvMassConstraint(const Double_t Minv);   //based on Invariant mass of added particles
    virtual void AddSubInvMassConstraint(const Int_t Np, const Int_t pid[], const Double_t Minv);//Add invariant mass constraint to subset of particles
    virtual void AddTotEnergyConstraint(const Double_t Etot);  //based on total energy of added particles
    virtual void AddTotMomentumConstraint(const TVector3 mom); //based on total 3 momentum of added particles
    virtual void AddSubMissMassConstraint(const TLorentzVector Mom, const Int_t Np, const Int_t pid[], const Double_t MissMass); // Based on missing mass of particles in subset

    virtual Int_t   GetIterations()    {return nIter;}

    virtual void    Reset()            {nIter=0; GKinFitter::Reset();}
    virtual Int_t   Solve();
};




#endif

