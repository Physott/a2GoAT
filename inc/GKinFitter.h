#ifndef _GKinFitter_h
#define _GKinFitter_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"


#define GKinFitter_CBRadius     0.258
#define GKinFitter_RadiatorDist 10.45

class GKinFitterPolarToCartesian
{
private:
    Double_t    p[19];
    Double_t    dp[19];

public:
    GKinFitterPolarToCartesian();
    ~GKinFitterPolarToCartesian();

    TMatrixD    GetParametersH()            const;
    TMatrixD    GetDerivatedParametersH()   const;
    TMatrixD    GetCovarianceH()            const;

    void        Set(const Double_t beamEnergy, const Double_t delta_beamEnergy);
    void        Set(const Int_t particleIndex, const Double_t E, const Double_t theta, const Double_t phi, const Double_t delta_E, const Double_t delta_theta, const Double_t delta_phi);
};

/*
class GMyTrackH
{
private:
    Double_t    p[40];
    TMatrixD    cm;

public:
    GMyTrackH(const GKinFitterPolarToCartesian& m);
    ~GMyTrackH();

    TMatrixD    GetParametersW()            const;
    TMatrixD    GetDerivatedParametersW()   const;
    TMatrixD    GetCovarianceW()            const;
};*/


class   GKinFitter
{
private:
    Int_t nPar; //Number of parameters=Npart*fNvar
    Int_t nUnk; //Number of unknown parameters
    Int_t nCon; //Number of constraints
    Int_t fNiter; // Number of times Solve has been called
    TMatrixD par0;      //vector of measured parameters
    TMatrixD par;       //vector of fitted parameters
    TMatrixD unk0;      //vector of unknowns calculated from constrain and par0
    TMatrixD unk;       //vector of fitted unknowns
    TMatrixD lambda;  	//Vector of lagrangian multipliers
    TMatrixD V0;        //Covariance matrix of measured parameters
    TMatrixD V;         //Covariance matrix of fitted parameters
    //TMatrixD fmD;      	//Matrix of constraint derivitives
    //TMatrixD fmd;      	//Vector of evaluated constraints
    //TMatrixD fmlamda;  	//Vector of lagrangian multipliers
    //TMatrixD fmV_D;    	//Covariance matrix of constraints (TO BE INVERTED)
    Double_t chi2;
    //TLorentzVector fPtot;
    //Bool_t          solved;

    /*TVector3        GetMainVertex(const TMatrixD &u)                                    {return TVector3(u[0][0], u[1][0], u[2][0]);}
    TVector3        GetMainVertex_Derivated_Unk(const Int_t index_m, const TMatrixD &u);
    TVector3        GetVertex(const Int_t index, const TMatrixD &u)                     {return TVector3(u[(index*3)+3][0], u[(index*3)+4][0], u[(index*3)+5][0]);}
    TVector3        GetVertex_Derivated_Unk(const Int_t index, const Int_t index_m, const TMatrixD &u);
    TLorentzVector  GetBeam(const TMatrixD x, const TMatrixD u);
    TLorentzVector  GetBeam_Derivated_Par(const Int_t index_n, const TMatrixD x, const TMatrixD u);
    TLorentzVector  GetBeam_Derivated_Unk(const Int_t index_m, const TMatrixD x, const TMatrixD u);
    TLorentzVector  GetParticle(const Int_t index, const TMatrixD x, const TMatrixD u);
    Double_t        Get_g(const Int_t index_k, const TMatrixD x, const TMatrixD u);
    TMatrixD        Get_g_Derivated_Par(const Int_t index_k, const TMatrixD x, const TMatrixD u);
    TMatrixD        Get_g_Derivated_Unk(const Int_t index_k, const TMatrixD x, const TMatrixD u);*/

    void        Get_lv(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret);
    void        Get_lv_Derivated_Par(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret);
    void        Get_lv_Derivated_Unk(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret);

    void        Get_g(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret);
    void        Get_g_Derivated_Par(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret);
    void        Get_g_Derivated_Unk(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret);
    void        Get_r(const TMatrixD& x, const TMatrixD& u, const TMatrixD& g_Derivated_Par, const TMatrixD& g, TMatrixD& ret);
    void        Get_S(const TMatrixD& g_Derivated_Par, TMatrixD& ret);
    Double_t    Calc_Chi2(const TMatrixD& x, const TMatrixD& g);
    Bool_t      SolveStep(const TMatrixD& x, const TMatrixD& u, TMatrixD &new_x, TMatrixD &new_u, Double_t& chiSq);

protected:
    static  TMatrixD    ConvertMtoH(const Double_t* pIn, const Double_t* dpIn, Double_t* pOut);
    static  TMatrixD    ConvertHtoW(const Double_t* pIn, const TMatrixD& dp, const Double_t mass, Double_t *pOut);

    void    Constraints();

public:
    GKinFitter();
    ~GKinFitter();

    Double_t        ConfidenceLevel()       {return TMath::Prob(chi2,4);}//Note should be Ncon-Nunknowns
    Double_t        GetChi2()               {return chi2;}
    TLorentzVector  GetEta();
    TLorentzVector  GetEtap();
    TLorentzVector  GetPi0a();
    TLorentzVector  GetPi0b();
    Double_t        GetPull(const Int_t i)  {return (par0[i][0]-par0[i][0])/sqrt(V0[i][i]-V[i][i]);}
    //Bool_t          IsSolved()              {return solved;}
    void            Set(const GKinFitterPolarToCartesian& v);
    void            Reset();
    Bool_t          Solve();
    Bool_t          ReSolve();
};

#endif

