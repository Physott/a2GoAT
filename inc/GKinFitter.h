#ifndef _GKinFitter_h
#define _GKinFitter_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "GKinFitterConverter.h"


class GKinFitter4VectorX: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

        return x[3]*helpX/r;
    }
    unsigned int NDim() const                               {return 7;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitter4VectorX();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));
        switch(ipar)
        {
        case 0:
            return x[3]*(helpY + helpZ)/(r*r*r);
        case 1:
            return -x[3]*helpX*helpY/(r*r*r);
        case 2:
            return -x[3]*helpX*helpZ/(r*r*r);
        case 3:
            return helpX/r;
        case 4:
            return x[3]*(helpY + helpZ)/(r*r*r);
        case 5:
            return x[3]*helpX*helpY/(r*r*r);
        case 6:
            return x[3]*helpX*helpZ/(r*r*r);
        default:
            return DoEval(x);
        }
    }
};
class GKinFitter4VectorY: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

        return x[3]*helpY/r;
    }
    unsigned int NDim() const                               {return 7;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitter4VectorY();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));
        switch(ipar)
        {
        case 0:
            return -x[3]*helpY*helpX/(r*r*r);
        case 1:
            return x[3]*(helpX + helpZ)/(r*r*r);
        case 2:
            return -x[3]*helpY*helpZ/(r*r*r);
        case 3:
            return helpY/r;
        case 4:
            return x[3]*helpY*helpX/(r*r*r);
        case 5:
            return x[3]*(helpX + helpZ)/(r*r*r);
        case 6:
            return x[3]*helpY*helpZ/(r*r*r);
        default:
            return DoEval(x);
        }
    }
};
class GKinFitter4VectorZ: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

        return x[3]*helpZ/r;
    }
    unsigned int NDim() const                               {return 7;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitter4VectorZ();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        Double_t    helpX = x[0] - x[4];
        Double_t    helpY = x[1] - x[5];
        Double_t    helpZ = x[2] - x[6];
        Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));
        switch(ipar)
        {
        case 0:
            return -x[3]*helpZ*helpX/(r*r*r);
        case 1:
            return -x[3]*helpZ*helpY/(r*r*r);
        case 2:
            return x[3]*(helpX + helpY)/(r*r*r);
        case 3:
            return helpZ/r;
        case 4:
            return x[3]*helpZ*helpX/(r*r*r);
        case 5:
            return x[3]*helpZ*helpY/(r*r*r);
        case 6:
            return x[3]*(helpX*helpY)/(r*r*r);
        default:
            return DoEval(x);
        }
    }
};

class GKinFitterConstrainInvMass: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        Double_t    helpX = x[0] + x[4];
        Double_t    helpY = x[1] + x[5];
        Double_t    helpZ = x[2] + x[6];
        Double_t    helpE = x[3] + x[7];

        return (helpE*helpE)-(helpX*helpX)-(helpY*helpY)-(helpZ*helpZ)-(x[8]*x[8]);
    }
    unsigned int NDim() const                               {return 9;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitterConstrainInvMass();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        if(ipar==8)
            return -2*x[8];
        switch(ipar)
        {
        case 0:
        case 4:
            return -2*(x[0] + x[4]);
        case 1:
        case 5:
            return -2*(x[1] + x[5]);
        case 2:
        case 6:
            return -2*(x[2] + x[6]);
        case 3:
        case 7:
            return 2*(x[3] + x[7]);
        }
    }
};
class GKinFitterConstrainMisMass: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        Double_t    helpX = x[0] - x[4] - x[8]  - x[12] - x[16] - x[20] - x[24];
        Double_t    helpY = x[1] - x[5] - x[9]  - x[13] - x[17] - x[21] - x[25];
        Double_t    helpZ = x[2] - x[6] - x[10] - x[14] - x[18] - x[22] - x[26];
        Double_t    helpE = x[3] - x[7] - x[11] - x[15] - x[19] - x[23] - x[27];

        return (helpE*helpE)-(helpX*helpX)-(helpY*helpY)-(helpZ*helpZ)-(x[28]*x[28]);
    }
    unsigned int NDim() const                               {return 29;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitterConstrainMisMass();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        if(ipar==28)
            return 2*x[28];
        switch(ipar)
        {
        case 0:
        case 4:
        case 8:
        case 12:
        case 16:
        case 20:
        case 24:
            return -2*(x[0] - x[4] - x[8]  - x[12] - x[16] - x[20] - x[24]);
        case 1:
        case 5:
        case 9:
        case 13:
        case 17:
        case 21:
        case 25:
            return -2*(x[1] - x[5] - x[9]  - x[13] - x[17] - x[21] - x[25]);
        case 2:
        case 6:
        case 10:
        case 14:
        case 18:
        case 22:
        case 26:
            return -2*(x[2] - x[6] - x[10] - x[14] - x[18] - x[22] - x[26]);
        case 3:
        case 7:
        case 11:
        case 15:
        case 19:
        case 23:
        case 27:
            return 2*(x[3] - x[7] - x[11] - x[15] - x[19] - x[23] - x[27]);
        }
    }
};
class GKinFitterConstrainBeam: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const
    {
        return ((x[3]-938.272046)*(x[3]-938.272046))-(x[0]*x[0])-(x[1]*x[1])-(x[2]*x[2]);
    }
    unsigned int NDim() const                               {return 4;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitterConstrainBeam();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        switch(ipar)
        {
        case 0:
            return -2*x[0];
        case 1:
            return -2*x[1];
        case 2:
            return -2*x[2];
        case 3:
            return 2*(x[3]-938.272046);
        default:
            DoEval(x);
        }
    }
};



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

    Double_t    GetReactionPx(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i=-1);
    Double_t    GetReactionPy(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i=-1);
    Double_t    GetReactionPz(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i=-1);
    Double_t    GetReactionE(const TMatrixD& x, const Int_t derivate_i=-1)                                  {if(derivate_i<0) return x[0][0]+938.272046; if(derivate_i==0) return 1; return 0;}
    Double_t    GetPhotonPx(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i=-1);
    Double_t    GetPhotonPy(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i=-1);
    Double_t    GetPhotonPz(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i=-1);
    Double_t    GetPhotonE(const Int_t i, const TMatrixD& x, const Int_t derivate_i=-1)                     {if(derivate_i<0) return x[(4*i)+4][0]; if(derivate_i==(4*i)+4) return 1; return 0;}
    Double_t    GetConstraint(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i=-1);

    Double_t    GetReactionPx_DeriPar(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);
    Double_t    GetReactionPy_DeriPar(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);
    Double_t    GetReactionPz_DeriPar(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);
    Double_t    GetReactionE_DeriPar(const TMatrixD& x, const Int_t derivate_i)                                  {if(derivate_i==0) return 1; return 0;}
    Double_t    GetPhotonPx_DeriPar(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);
    Double_t    GetPhotonPy_DeriPar(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);
    Double_t    GetPhotonPz_DeriPar(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);
    Double_t    GetPhotonE_DeriPar(const Int_t i, const TMatrixD& x, const Int_t derivate_i)                     {if(derivate_i==(4*i)+4) return 1; return 0;}
    Double_t    GetConstraint_DeriPar(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);

    Double_t    GetReactionPx_DeriUnk(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)             {return GetReactionPx(x, u, derivate_i);}
    Double_t    GetReactionPy_DeriUnk(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)             {return GetReactionPy(x, u, derivate_i);}
    Double_t    GetReactionPz_DeriUnk(const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i)             {return GetReactionPz(x, u, derivate_i);}
    Double_t    GetReactionE_DeriUnk()                                                                          {return 0;}
    Double_t    GetPhotonPx_DeriUnk(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);
    Double_t    GetPhotonPy_DeriUnk(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);
    Double_t    GetPhotonPz_DeriUnk(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);
    Double_t    GetPhotonE_DeriUnk()                                                                            {return 0;}
    Double_t    GetConstraint_DeriUnk(const Int_t i, const TMatrixD& x, const TMatrixD& u, const Int_t derivate_i);

    void        Get_g(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret);
    void        Get_g_DeriPar(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret);
    void        Get_g_DeriUnk(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret);
    void        Get_r(const TMatrixD& x, const TMatrixD& u, const TMatrixD& g_DeriPar, const TMatrixD& g, TMatrixD& ret);
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

    Double_t    GetReactionPx(const Int_t derivate_i=-1)   {return GetReactionPx(par, unk, derivate_i);}
    Double_t    GetReactionPy(const Int_t derivate_i=-1)   {return GetReactionPy(par, unk, derivate_i);}
    Double_t    GetReactionPz(const Int_t derivate_i=-1)   {return GetReactionPz(par, unk, derivate_i);}
    Double_t    GetReactionE(const Int_t derivate_i=-1)    {return GetReactionE(par, derivate_i);}
    Double_t    GetPhotonPx(const Int_t i, const Int_t derivate_i=-1)   {return GetPhotonPx(i, par, unk, derivate_i);}
    Double_t    GetPhotonPy(const Int_t i, const Int_t derivate_i=-1)   {return GetPhotonPy(i, par, unk, derivate_i);}
    Double_t    GetPhotonPz(const Int_t i, const Int_t derivate_i=-1)   {return GetPhotonPz(i, par, unk, derivate_i);}
    Double_t    GetPhotonE(const Int_t i, const Int_t derivate_i=-1)    {return GetPhotonE(i, par, derivate_i);}
    Double_t    GetInitialReactionPx(const Int_t derivate_i=-1)   {return GetReactionPx(par0, unk0, derivate_i);}
    Double_t    GetInitialReactionPy(const Int_t derivate_i=-1)   {return GetReactionPy(par0, unk0, derivate_i);}
    Double_t    GetInitialReactionPz(const Int_t derivate_i=-1)   {return GetReactionPz(par0, unk0, derivate_i);}
    Double_t    GetInitialReactionE(const Int_t derivate_i=-1)    {return GetReactionE(par0, derivate_i);}
    Double_t    GetInitialPhotonPx(const Int_t i, const Int_t derivate_i=-1)   {return GetPhotonPx(i, par0, unk0, derivate_i);}
    Double_t    GetInitialPhotonPy(const Int_t i, const Int_t derivate_i=-1)   {return GetPhotonPy(i, par0, unk0, derivate_i);}
    Double_t    GetInitialPhotonPz(const Int_t i, const Int_t derivate_i=-1)   {return GetPhotonPz(i, par0, unk0, derivate_i);}
    Double_t    GetInitialPhotonE(const Int_t i, const Int_t derivate_i=-1)    {return GetPhotonE(i, par0, derivate_i);}

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

