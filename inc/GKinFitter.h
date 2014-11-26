#ifndef _GKinFitter_h
#define _GKinFitter_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"
/*
#include "GKinFitterParticle.h"

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
};
*/

#define GMyTrackM_CBRadius  1
#define GMyTrackM_nParameters  33
#define GMyTrackH_nParameters  40
#define GMyTrackW_nParameters  40

class GMyTrackM
{
private:
    Double_t    p[33];
    Double_t    dp[33];

public:
    GMyTrackM(const Double_t beamEnergy, const TLorentzVector& p0, const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& p3, const TLorentzVector& p4, const TLorentzVector& p5);
    ~GMyTrackM();

    TMatrixD    GetParametersH()            const;
    TMatrixD    GetDerivatedParametersH()   const;
    TMatrixD    GetCovarianceH()            const;
};


class GMyTrackH
{
private:
    Double_t    p[40];
    TMatrixD    cm;

public:
    GMyTrackH(const GMyTrackM& m);
    ~GMyTrackH();

    TMatrixD    GetParametersW()            const;
    TMatrixD    GetDerivatedParametersW()   const;
    TMatrixD    GetCovarianceW()            const;
};


class   GKinFitter
{
private:
    Int_t nPar; //Number of parameters=Npart*fNvar
    Int_t nCon; //Number of constraints
    Int_t fNiter; // Number of times Solve has been called
    TMatrixD fmAlpha0;	//original parameters
    TMatrixD fmAlpha;  	//fitted parameters
    TMatrixD fmV_Alpha0;//Covariance matrix for original parameters
    TMatrixD fmV_Alpha; //Covariance matrix for fitted parameters
    TMatrixD fmD;      	//Matrix of constraint derivitives
    TMatrixD fmd;      	//Vector of evaluated constraints
    TMatrixD fmlamda;  	//Vector of lagrangian multipliers
    TMatrixD fmV_D;    	//Covariance matrix of constraints (TO BE INVERTED)
    Double_t fchi2;
    TLorentzVector fPtot;
    Bool_t          solved;

protected:
    static  TMatrixD    ConvertMtoH(const Double_t* pIn, const Double_t* dpIn, Double_t* pOut);
    static  TMatrixD    ConvertHtoW(const Double_t* pIn, const TMatrixD& dp, const Double_t mass, Double_t *pOut);

    void    Constraints();

public:
    GKinFitter();
    ~GKinFitter();

    Double_t        ConfidenceLevel()       {return TMath::Prob(fchi2,4);}//Note should be Ncon-Nunknowns
    Double_t        GetChi2()               {return fchi2;}
    TLorentzVector  GetEta();
    TLorentzVector  GetEtap();
    TLorentzVector  GetPi0a();
    TLorentzVector  GetPi0b();
    Double_t        GetPull(const Int_t i)  {return (fmAlpha0[i][0]-fmAlpha[i][0])/sqrt(fmV_Alpha0[i][i]-fmV_Alpha[i][i]);}
    TVector3        GetVertex()             {return TVector3(fmAlpha[0][0], fmAlpha[1][0], fmAlpha[2][0]);}
    TVector3        GetVertex_sub0()        {return TVector3(fmAlpha[0][0] + fmAlpha[3][0], fmAlpha[1][0] + fmAlpha[4][0], fmAlpha[2][0] + fmAlpha[5][0]);}
    TVector3        GetVertex_sub1()        {return TVector3(fmAlpha[0][0] + fmAlpha[6][0], fmAlpha[1][0] + fmAlpha[5][0], fmAlpha[2][0] + fmAlpha[6][0]);}
    TVector3        GetVertex_sub2()        {return TVector3(fmAlpha[0][0] + fmAlpha[9][0], fmAlpha[1][0] + fmAlpha[10][0], fmAlpha[2][0] + fmAlpha[11][0]);}
    Bool_t          IsSolved()              {return solved;}
    void            Set(const Double_t beamEnergy, const TLorentzVector& p0, const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& p3, const TLorentzVector& p4, const TLorentzVector& p5);
    void            Reset();
    Bool_t          Solve();
};

#endif

