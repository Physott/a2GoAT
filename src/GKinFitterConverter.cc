//////////////////////////////////////////////////////////////////
//GKinFitter
//////////////////////////////////////////////////////////////////


#include "GKinFitterConverter.h"
#include <iostream>


GKinFitterPolarToCartesian::GKinFitterPolarToCartesian()
{

}

GKinFitterPolarToCartesian::~GKinFitterPolarToCartesian()
{
}

void    GKinFitterPolarToCartesian::Set(const Double_t beamEnergy, const Double_t delta_beamEnergy)
{
    p[0]   = beamEnergy;

    dp[0]   = delta_beamEnergy;
    std::cout << p[0] << "          " << dp[0] << std::endl;
}

void    GKinFitterPolarToCartesian::Set(const Int_t particleIndex, const Double_t E, const Double_t theta, const Double_t phi, const Double_t delta_E, const Double_t delta_theta, const Double_t delta_phi)
{
    p[(particleIndex*3) +1]   = E;
    p[(particleIndex*3) +2]   = theta;
    p[(particleIndex*3) +3]   = phi;

    dp[(particleIndex*3) +1]   = delta_E;
    dp[(particleIndex*3) +2]   = delta_theta;
    dp[(particleIndex*3) +3]   = delta_phi;

    std::cout << E << "          " << delta_E << std::endl;
    std::cout << theta << "          " << delta_theta << std::endl;
    std::cout << phi << "          " << delta_phi << std::endl;
}

TMatrixD    GKinFitterPolarToCartesian::GetParametersH()    const
{
    TMatrixD    ret(25,1);

    ret[0][0]   = p[0];

    for(int i=0; i<6; i++)
    {
        Double_t hhh[3] = {GKinFitter_CBRadius, p[(i*3)+2], p[(i*3)+3]};
        ret[(i*4)+1][0]  = GKinFitterPolarToCartesianX().DoEval(hhh);
        ret[(i*4)+2][0]  = GKinFitterPolarToCartesianY().DoEval(hhh);
        ret[(i*4)+3][0]  = GKinFitterPolarToCartesianZ().DoEval(hhh);
        ret[(i*4)+4][0]  = p[(i*3)+1];
    }

    return ret;
}

TMatrixD    GKinFitterPolarToCartesian::GetDerivatedParametersH()    const
{
    TMatrixD    ret(25,19);
    for(int i=0; i<25; i++)
    {
        for(int j=0; j<19; j++)
            ret[i][j]   = 0;
    }

    ret[0][0]  = 1;

    for(int i=0; i<6; i++)
    {
        Double_t hhh[3] = {GKinFitter_CBRadius, p[(i*3)+2], p[(i*3)+3]};
        ret[(i*4)+1][(i*3)+2] = GKinFitterPolarToCartesianX().Derivative(hhh, 1);
        ret[(i*4)+1][(i*3)+3] = GKinFitterPolarToCartesianX().Derivative(hhh, 2);
        ret[(i*4)+2][(i*3)+2] = GKinFitterPolarToCartesianY().Derivative(hhh, 1);
        ret[(i*4)+2][(i*3)+3] = GKinFitterPolarToCartesianY().Derivative(hhh, 2);
        ret[(i*4)+3][(i*3)+2] = GKinFitterPolarToCartesianZ().Derivative(hhh, 1);
        ret[(i*4)+4][(i*3)+1] = 1;
    }

    return ret;
}

TMatrixD    GKinFitterPolarToCartesian::GetCovarianceH()    const
{
    TMatrixD    cm(19, 19);
    for(int i=0; i<19; i++)
    {
        for(int j=0; j<19; j++)
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


