#ifndef _GKinFitterConverter_h
#define _GKinFitterConverter_h

#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/IFunction.h"


#define GKinFitter_CBRadius     0.258
#define GKinFitter_RadiatorDist 10.45



class GKinFitterPolarToCartesianX: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const                    {return x[0] * sin(x[1]) * cos(x[2]);}
    unsigned int NDim() const                               {return 3;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitterPolarToCartesianX();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        switch(ipar)
        {
        case 0:
            return sin(x[1]) * cos(x[2]);
        case 1:
            return x[0] * cos(x[1]) * cos(x[2]);
        case 2:
            return - x[0] * sin(x[1]) * sin(x[2]);
        }
    }
};
class GKinFitterPolarToCartesianY: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const                    {return x[0] * sin(x[1]) * sin(x[2]);}
    unsigned int NDim() const                               {return 3;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitterPolarToCartesianY();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        switch(ipar)
        {
        case 0:
            return sin(x[1]) * sin(x[2]);
        case 1:
            return x[0] * cos(x[1]) * sin(x[2]);
        case 2:
            return x[0] * sin(x[1]) * cos(x[2]);
        }
    }
};
class GKinFitterPolarToCartesianZ: public ROOT::Math::IGradientFunctionMultiDim
{
public:
    double DoEval(const double* x) const                    {return x[0] * cos(x[1]);}
    unsigned int NDim() const                               {return 2;}
    ROOT::Math::IGradientFunctionMultiDim* Clone() const    {return new GKinFitterPolarToCartesianZ();}
    double DoDerivative(const double* x, unsigned int ipar) const
    {
        if(ipar == 0)
            return cos(x[1]);
        else
            return - x[0] * sin(x[1]);
    }
};


class   GKinFitterPolarToCartesian
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


#endif

