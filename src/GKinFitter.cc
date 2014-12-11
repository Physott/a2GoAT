//////////////////////////////////////////////////////////////////
//GKinFitter
//////////////////////////////////////////////////////////////////


#include "GKinFitter.h"
#include "TDecompLU.h"
#include <iostream>


#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046


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
}

void    GKinFitterPolarToCartesian::Set(const Int_t particleIndex, const Double_t E, const Double_t theta, const Double_t phi, const Double_t delta_E, const Double_t delta_theta, const Double_t delta_phi)
{
    p[(particleIndex*3) +1]   = E;
    p[(particleIndex*3) +2]   = theta;
    p[(particleIndex*3) +3]   = phi;

    dp[(particleIndex*3) +1]   = delta_E;
    dp[(particleIndex*3) +2]   = delta_theta;
    dp[(particleIndex*3) +3]   = delta_phi;
}

TMatrixD    GKinFitterPolarToCartesian::GetParametersH()    const
{
    TMatrixD    ret(25,1);

    ret[0][0]   = p[0];

    for(int i=0; i<6; i++)
    {
        ret[(i*4)+1][0]  = GKinFitter_CBRadius * TMath::Sin(p[(i*3)+2]) * TMath::Cos(p[(i*3)+3]);
        ret[(i*4)+2][0]  = GKinFitter_CBRadius * TMath::Sin(p[(i*3)+2]) * TMath::Sin(p[(i*3)+3]);
        ret[(i*4)+3][0]  = GKinFitter_CBRadius * TMath::Cos(p[(i*3)+2]);
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
        ret[(i*4)+1][(i*3)+2] =   GKinFitter_CBRadius * TMath::Cos(p[(i*3)+2]) * TMath::Cos(p[(i*3)+3]);
        ret[(i*4)+1][(i*3)+3] = - GKinFitter_CBRadius * TMath::Sin(p[(i*3)+2]) * TMath::Sin(p[(i*3)+3]);
        ret[(i*4)+2][(i*3)+2] =   GKinFitter_CBRadius * TMath::Cos(p[(i*3)+2]) * TMath::Sin(p[(i*3)+3]);
        ret[(i*4)+2][(i*3)+3] =   GKinFitter_CBRadius * TMath::Sin(p[(i*3)+2]) * TMath::Cos(p[(i*3)+3]);
        ret[(i*4)+3][(i*3)+2] = - GKinFitter_CBRadius * TMath::Sin(p[(i*3)+2]);
        ret[(i*4)+4][(i*3)+2] =   p[(i*3)+1];
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



/*
GMyTrackH::GMyTrackH(const GKinFitterPolarToCartesian& m):
    cm(40,40)
{
    TMatrixD    help(m.GetParametersH());
    for(int i=0; i<40; i++)
        p[i]    = help[i][0];

    cm  = m.GetCovarianceH();
}

GMyTrackH::~GMyTrackH()
{

}

TMatrixD    GMyTrackH::GetParametersW()    const
{

    TMatrixD    ret(40,1);

    for(int i=0; i<12; i++)
        ret[i][0]  = p[i];

    Double_t    help;
    help    =  p[15]-MASS_PROTON;
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

    //std::cout << "first W  ";
    //TLorentzVector(ret[12][0], ret[13][0], ret[14][0], ret[15][0]).Print();
    return ret;
}

TMatrixD    GMyTrackH::GetDerivatedParametersW()    const
{
    TMatrixD    ret(40,40);
    for(int i=0; i<40; i++)
    {
        for(int j=0; j<40; j++)
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

    return ret;
}

TMatrixD    GMyTrackH::GetCovarianceW()    const
{
    TMatrixD    dp(GetDerivatedParametersW());
    TMatrixD    dpT(dp);
    dpT.T();

    return dp * cm * dpT;
}*/









GKinFitter::GKinFitter()    :
    nPar((6*4)+1),
    nUnk(5),
    nCon(1+3),
    fNiter(0),
    par0(nPar,1),
    par(nPar,1),
    unk0(nUnk,1),
    unk(nUnk,1),
    lambda(nCon,1),
    V0(nPar,nPar),
    V(nPar,nPar),
    //fmD(nCon,nPar),
    //fmd(nCon,1),
    //fmlamda(nCon,1),
    //fmV_D(nCon,nCon),
    chi2(0)//,
    //fPtot(0, 0, 0, 0),
    //solved(kFALSE)
{
}

GKinFitter::~GKinFitter()
{

}

void    GKinFitter::Get_lv(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret)
{
    Double_t    helpX = u[0][0]-u[3][0];
    Double_t    helpY = u[1][0]-u[4][0];
    Double_t    helpZ = u[2][0]+GKinFitter_RadiatorDist;
    Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

    ret[0][0]   = helpX * x[0][0] / r;
    ret[1][0]   = helpY * x[0][0] / r;
    ret[2][0]   = helpZ * x[0][0] / r;
    ret[3][0]   = x[0][0];

    for(int i=0; i<6; i++)
    {
        helpX = x[(4*i)+1][0]-u[0][0];
        helpY = x[(4*i)+2][0]-u[1][0];
        helpZ = x[(4*i)+3][0]-u[2][0];
        r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

        ret[(4*i)+4][0]   = helpX * x[(4*i)+4][0] / r;
        ret[(4*i)+5][0]   = helpY * x[(4*i)+4][0] / r;
        ret[(4*i)+6][0]   = helpZ * x[(4*i)+4][0] / r;
        ret[(4*i)+7][0]   = x[(4*i)+4][0];
    }
}

void        GKinFitter::Get_lv_Derivated_Par(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret)
{
    Double_t    helpX = u[0][0]-u[3][0];
    Double_t    helpY = u[1][0]-u[4][0];
    Double_t    helpZ = u[2][0]+GKinFitter_RadiatorDist;
    Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

    for(int i=0; i<4; i++)
    {
        for(int j=1; j<(6*4)+1; j++)
            ret[i][j]   = 0;
    }
    ret[0][0]   = helpX / r;
    ret[1][0]   = helpY / r;
    ret[2][0]   = helpZ / r;
    ret[3][0]   = 1;

    for(int i=0; i<6; i++)
    {
        helpX = x[(4*i)+1][0]-u[0][0];
        helpY = x[(4*i)+2][0]-u[1][0];
        helpZ = x[(4*i)+3][0]-u[2][0];
        r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

        for(int j=0; j<(4*i)+1; j++)
            ret[(4*i)+4][j]     = 0;
        ret[(4*i)+4][(4*i)+1]   = (helpY + helpZ) * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+4][(4*i)+2]   = helpX * helpY * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+4][(4*i)+3]   = helpX * helpZ * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+4][(4*i)+4]   = helpX / r;
        for(int j=(4*i)+5; j<(6*4)+1; j++)
            ret[(4*i)+4][j]     = 0;

        for(int j=0; j<(4*i)+1; j++)
            ret[(4*i)+5][j]     = 0;
        ret[(4*i)+5][(4*i)+1]   = helpY * helpX * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+5][(4*i)+2]   = (helpX + helpZ) * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+5][(4*i)+3]   = helpY * helpZ * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+5][(4*i)+4]   = helpY / r;
        for(int j=(4*i)+5; j<(6*4)+1; j++)
            ret[(4*i)+5][j]     = 0;

        for(int j=0; j<(4*i)+1; j++)
            ret[(4*i)+6][j]     = 0;
        ret[(4*i)+6][(4*i)+1]   = helpZ * helpX * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+6][(4*i)+2]   = helpZ * helpY * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+6][(4*i)+3]   = (helpX + helpY) * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+6][(4*i)+4]   = helpZ / r;
        for(int j=(4*i)+5; j<(6*4)+1; j++)
            ret[(4*i)+6][j]     = 0;

        for(int j=0; j<(4*i)+4; j++)
            ret[(4*i)+6][j]     = 0;
        ret[(4*i)+7][(4*i)+4]   = 1;
        for(int j=(4*i)+5; j<(6*4)+1; j++)
            ret[(4*i)+6][j]     = 0;
    }
}

void        GKinFitter::Get_lv_Derivated_Unk(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret)
{
    Double_t    helpX = u[0][0]-u[3][0];
    Double_t    helpY = u[1][0]-u[4][0];
    Double_t    helpZ = u[2][0]+GKinFitter_RadiatorDist;
    Double_t    r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

    ret[0][0]   = (helpY + helpZ) * x[0][0] / (r*r*r);
    ret[0][1]   = helpX * helpY * x[0][0] / (r*r*r);
    ret[0][2]   = helpX * helpZ * x[0][0] / (r*r*r);
    ret[0][3]   = -(helpY + helpZ) * x[0][0] / (r*r*r);
    ret[0][4]   = -helpX * helpY * x[0][0] / (r*r*r);

    ret[1][0]   = helpY * helpX * x[0][0] / (r*r*r);
    ret[1][1]   = (helpX + helpZ) * x[0][0] / (r*r*r);
    ret[1][2]   = helpY * helpZ * x[0][0] / (r*r*r);
    ret[1][3]   = -helpY * helpX * x[0][0] / (r*r*r);
    ret[1][4]   = -(helpX + helpZ) * x[0][0] / (r*r*r);

    ret[2][0]   = helpZ * helpX * x[0][0] / (r*r*r);
    ret[2][1]   = helpZ * helpY * x[0][0] / (r*r*r);
    ret[2][2]   = (helpX + helpY) * x[0][0] / (r*r*r);
    ret[2][3]   = -helpZ * helpX * x[0][0] / (r*r*r);
    ret[2][4]   = -helpZ * helpY * x[0][0] / (r*r*r);

    ret[3][0]   = 0;
    ret[3][1]   = 0;
    ret[3][2]   = 0;
    ret[3][3]   = 0;
    ret[3][4]   = 0;

    for(int i=0; i<6; i++)
    {
        helpX = x[(4*i)+1][0]-u[0][0];
        helpY = x[(4*i)+2][0]-u[1][0];
        helpZ = x[(4*i)+3][0]-u[2][0];
        r     = TMath::Sqrt((helpX*helpX)+(helpY*helpY)+(helpZ*helpZ));

        ret[(4*i)+4][0]   = -(helpY + helpZ) * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+4][1]   = -helpX * helpY * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+4][2]   = -helpX * helpZ * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+4][3]   = 0;
        ret[(4*i)+4][4]   = 0;

        ret[(4*i)+5][0]   = -helpY * helpX * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+5][1]   = -(helpX + helpZ) * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+5][2]   = -helpY * helpZ * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+5][3]   = 0;
        ret[(4*i)+5][4]   = 0;

        ret[(4*i)+6][0]   = -helpZ * helpX * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+6][1]   = -helpZ * helpY * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+6][2]   = -(helpX + helpY) * x[(4*i)+4][0] / (r*r*r);
        ret[(4*i)+6][3]   = 0;
        ret[(4*i)+6][4]   = 0;

        ret[(4*i)+7][0]   = 0;
        ret[(4*i)+7][1]   = 0;
        ret[(4*i)+7][2]   = 0;
        ret[(4*i)+7][3]   = 0;
        ret[(4*i)+7][4]   = 0;
    }

}

void        GKinFitter::Get_g(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret)
{
    TMatrixD    lv(7*4, 1);
    Get_lv(x, u, lv);

    Double_t    helpX = lv[4][0]+lv[8][0];
    Double_t    helpY = lv[5][0]+lv[9][0];
    Double_t    helpZ = lv[6][0]+lv[10][0];
    Double_t    helpE = lv[7][0]+lv[11][0];
    ret[0][0]   = (helpE * helpE) - (helpX * helpX) - (helpY * helpY) - (helpZ * helpZ) - (MASS_ETA * MASS_ETA);

    helpX = lv[12][0]+lv[16][0];
    helpY = lv[13][0]+lv[17][0];
    helpZ = lv[14][0]+lv[18][0];
    helpE = lv[15][0]+lv[19][0];
    ret[1][0]   = (helpE * helpE) - (helpX * helpX) - (helpY * helpY) - (helpZ * helpZ) - (MASS_PI0 * MASS_PI0);

    helpX = lv[20][0]+lv[24][0];
    helpY = lv[21][0]+lv[25][0];
    helpZ = lv[22][0]+lv[26][0];
    helpE = lv[23][0]+lv[27][0];
    ret[2][0]   = (helpE * helpE) - (helpX * helpX) - (helpY * helpY) - (helpZ * helpZ) - (MASS_PI0 * MASS_PI0);

    helpX = lv[0][0] - lv[4][0] - lv[8][0]  - lv[12][0] - lv[16][0] - lv[20][0] - lv[24][0];
    helpY = lv[1][0] - lv[5][0] - lv[9][0]  - lv[13][0] - lv[17][0] - lv[21][0] - lv[25][0];
    helpZ = lv[2][0] - lv[6][0] - lv[10][0] - lv[14][0] - lv[18][0] - lv[22][0] - lv[26][0];
    helpE = lv[3][0] - lv[7][0] - lv[11][0] - lv[15][0] - lv[19][0] - lv[23][0] - lv[27][0];
    ret[3][0]   = (helpE * helpE) - (helpX * helpX) - (helpY * helpY) - (helpZ * helpZ) - (MASS_PROTON * MASS_PROTON);
}

void        GKinFitter::Get_g_Derivated_Par(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret)
{
    TMatrixD    lv(7*4, 1);
    TMatrixD    dlv_par(7*4, 25);
    Get_lv(x, u, lv);
    Get_lv_Derivated_Par(x, u, dlv_par);

    Double_t    helpX = lv[4][0]+lv[8][0];
    Double_t    helpY = lv[5][0]+lv[9][0];
    Double_t    helpZ = lv[6][0]+lv[10][0];
    Double_t    helpE = lv[7][0]+lv[11][0];
    ret[0][0]   = 0;
    for(int i=1; i<9; i++)
    {
        Double_t    dhelpX = dlv_par[4][i]+dlv_par[8][i];
        Double_t    dhelpY = dlv_par[5][i]+dlv_par[9][i];
        Double_t    dhelpZ = dlv_par[6][i]+dlv_par[10][i];
        Double_t    dhelpE = dlv_par[7][i]+dlv_par[11][i];
        ret[0][i]   = 2 * helpE * dhelpE - 2 * helpX * dhelpX - 2 * helpY * dhelpY - 2 * helpZ * dhelpZ;
    }
    for(int i=9; i<25; i++)
        ret[0][i]   = 0;


    helpX = lv[12][0]+lv[16][0];
    helpY = lv[13][0]+lv[17][0];
    helpZ = lv[14][0]+lv[18][0];
    helpE = lv[15][0]+lv[19][0];
    for(int i=0; i<9; i++)
        ret[1][i]   = 0;
    for(int i=9; i<17; i++)
    {
        Double_t    dhelpX = dlv_par[12][i]+dlv_par[16][i];
        Double_t    dhelpY = dlv_par[13][i]+dlv_par[17][i];
        Double_t    dhelpZ = dlv_par[14][i]+dlv_par[18][i];
        Double_t    dhelpE = dlv_par[15][i]+dlv_par[19][i];
        ret[1][i]   = 2 * helpE * dhelpE - 2 * helpX * dhelpX - 2 * helpY * dhelpY - 2 * helpZ * dhelpZ;
    }
    for(int i=17; i<25; i++)
        ret[1][i]   = 0;


    helpX = lv[20][0]+lv[24][0];
    helpY = lv[21][0]+lv[25][0];
    helpZ = lv[22][0]+lv[26][0];
    helpE = lv[23][0]+lv[27][0];
    for(int i=0; i<17; i++)
        ret[2][i]   = 0;
    for(int i=17; i<25; i++)
    {
        Double_t    dhelpX = dlv_par[20][i]+dlv_par[24][i];
        Double_t    dhelpY = dlv_par[21][i]+dlv_par[25][i];
        Double_t    dhelpZ = dlv_par[22][i]+dlv_par[26][i];
        Double_t    dhelpE = dlv_par[23][i]+dlv_par[27][i];
        ret[2][i]   = 2 * helpE * dhelpE - 2 * helpX * dhelpX - 2 * helpY * dhelpY - 2 * helpZ * dhelpZ;
    }

    helpX = lv[0][0] - lv[4][0] - lv[8][0]  - lv[12][0] - lv[16][0] - lv[20][0] - lv[24][0];
    helpY = lv[1][0] - lv[5][0] - lv[9][0]  - lv[13][0] - lv[17][0] - lv[21][0] - lv[25][0];
    helpZ = lv[2][0] - lv[6][0] - lv[10][0] - lv[14][0] - lv[18][0] - lv[22][0] - lv[26][0];
    helpE = lv[3][0] - lv[7][0] - lv[11][0] - lv[15][0] - lv[19][0] - lv[23][0] - lv[27][0];
    for(int i=0; i<25; i++)
    {
        Double_t    dhelpX = dlv_par[0][i] - dlv_par[4][i] - dlv_par[8][i]  - dlv_par[12][i] - dlv_par[16][i] - dlv_par[20][i] - dlv_par[24][i];
        Double_t    dhelpY = dlv_par[1][i] - dlv_par[5][i] - dlv_par[9][i]  - dlv_par[13][i] - dlv_par[17][i] - dlv_par[21][i] - dlv_par[25][i];
        Double_t    dhelpZ = dlv_par[2][i] - dlv_par[6][i] - dlv_par[10][i] - dlv_par[14][i] - dlv_par[18][i] - dlv_par[22][i] - dlv_par[26][i];
        Double_t    dhelpE = dlv_par[3][i] - dlv_par[7][i] - dlv_par[11][i] - dlv_par[15][i] - dlv_par[19][i] - dlv_par[23][i] - dlv_par[27][i];
        ret[3][i]   = 2 * helpE * dhelpE - 2 * helpX * dhelpX - 2 * helpY * dhelpY - 2 * helpZ * dhelpZ;
    }
}

void        GKinFitter::Get_g_Derivated_Unk(const TMatrixD& x, const TMatrixD& u, TMatrixD& ret)
{
    TMatrixD    lv(7*4, 1);
    TMatrixD    dlv_unk(7*4, 5);
    Get_lv(x, u, lv);
    Get_lv_Derivated_Unk(x, u, dlv_unk);

    Double_t    helpX = lv[4][0]+lv[8][0];
    Double_t    helpY = lv[5][0]+lv[9][0];
    Double_t    helpZ = lv[6][0]+lv[10][0];
    Double_t    helpE = lv[7][0]+lv[11][0];
    for(int i=0; i<3; i++)
    {
        Double_t    dhelpX = dlv_unk[4][i]+dlv_unk[8][i];
        Double_t    dhelpY = dlv_unk[5][i]+dlv_unk[9][i];
        Double_t    dhelpZ = dlv_unk[6][i]+dlv_unk[10][i];
        Double_t    dhelpE = dlv_unk[7][i]+dlv_unk[11][i];
        ret[0][i]   = 2 * helpE * dhelpE - 2 * helpX * dhelpX - 2 * helpY * dhelpY - 2 * helpZ * dhelpZ;
    }
    for(int i=3; i<5; i++)
        ret[0][i]   = 0;


    helpX = lv[12][0]+lv[16][0];
    helpY = lv[13][0]+lv[17][0];
    helpZ = lv[14][0]+lv[18][0];
    helpE = lv[15][0]+lv[19][0];
    for(int i=0; i<3; i++)
    {
        Double_t    dhelpX = dlv_unk[12][i]+dlv_unk[16][i];
        Double_t    dhelpY = dlv_unk[13][i]+dlv_unk[17][i];
        Double_t    dhelpZ = dlv_unk[14][i]+dlv_unk[18][i];
        Double_t    dhelpE = dlv_unk[15][i]+dlv_unk[19][i];
        ret[1][i]   = 2 * helpE * dhelpE - 2 * helpX * dhelpX - 2 * helpY * dhelpY - 2 * helpZ * dhelpZ;
    }
    for(int i=3; i<5; i++)
        ret[1][i]   = 0;


    helpX = lv[20][0]+lv[24][0];
    helpY = lv[21][0]+lv[25][0];
    helpZ = lv[22][0]+lv[26][0];
    helpE = lv[23][0]+lv[27][0];
    for(int i=0; i<3; i++)
    {
        Double_t    dhelpX = dlv_unk[20][i]+dlv_unk[24][i];
        Double_t    dhelpY = dlv_unk[21][i]+dlv_unk[25][i];
        Double_t    dhelpZ = dlv_unk[22][i]+dlv_unk[26][i];
        Double_t    dhelpE = dlv_unk[23][i]+dlv_unk[27][i];
        ret[2][i]   = 2 * helpE * dhelpE - 2 * helpX * dhelpX - 2 * helpY * dhelpY - 2 * helpZ * dhelpZ;
    }
    for(int i=3; i<5; i++)
        ret[2][i]   = 0;

    helpX = lv[0][0] - lv[4][0] - lv[8][0]  - lv[12][0] - lv[16][0] - lv[20][0] - lv[24][0];
    helpY = lv[1][0] - lv[5][0] - lv[9][0]  - lv[13][0] - lv[17][0] - lv[21][0] - lv[25][0];
    helpZ = lv[2][0] - lv[6][0] - lv[10][0] - lv[14][0] - lv[18][0] - lv[22][0] - lv[26][0];
    helpE = lv[3][0] - lv[7][0] - lv[11][0] - lv[15][0] - lv[19][0] - lv[23][0] - lv[27][0];
    for(int i=0; i<5; i++)
    {
        Double_t    dhelpX = dlv_unk[0][i] - dlv_unk[4][i] - dlv_unk[8][i]  - dlv_unk[12][i] - dlv_unk[16][i] - dlv_unk[20][i] - dlv_unk[24][i];
        Double_t    dhelpY = dlv_unk[1][i] - dlv_unk[5][i] - dlv_unk[9][i]  - dlv_unk[13][i] - dlv_unk[17][i] - dlv_unk[21][i] - dlv_unk[25][i];
        Double_t    dhelpZ = dlv_unk[2][i] - dlv_unk[6][i] - dlv_unk[10][i] - dlv_unk[14][i] - dlv_unk[18][i] - dlv_unk[22][i] - dlv_unk[26][i];
        Double_t    dhelpE = dlv_unk[3][i] - dlv_unk[7][i] - dlv_unk[11][i] - dlv_unk[15][i] - dlv_unk[19][i] - dlv_unk[23][i] - dlv_unk[27][i];
        ret[3][i]   = 2 * helpE * dhelpE - 2 * helpX * dhelpX - 2 * helpY * dhelpY - 2 * helpZ * dhelpZ;
    }
}

void        GKinFitter::Get_r(const TMatrixD& x, const TMatrixD& u, const TMatrixD& g_Derivated_Par, const TMatrixD& g, TMatrixD& ret)
{
    TMatrixD    dpar(par0);
                dpar -= x;
    ret = g + (g_Derivated_Par * dpar);
}

void        GKinFitter::Get_S(const TMatrixD& g_Derivated_Par, TMatrixD& ret)
{
    TMatrixD    g_Derivated_Par_T(g_Derivated_Par);
    g_Derivated_Par_T.T();

    ret = g_Derivated_Par * V0 * g_Derivated_Par_T;

    g_Derivated_Par.Print();
    V0.Print();
    g_Derivated_Par_T.Print();
    ret.Print();
}

Double_t    GKinFitter::Calc_Chi2(const TMatrixD& x, const TMatrixD& g)
{
    TMatrixD    dpar(par0-x);
    TMatrixD    dpar_T(dpar);
    dpar_T.T();
    TMatrixD    lambda_T(lambda);
    lambda_T.T();
    TMatrixD    help0(dpar_T * V0 * dpar);
    TMatrixD    help1(lambda_T * g);
    return  help0[0][0] + (2 * help1[0][0]);
}

/*
TVector3        GKinFitter::GetMainVertex_Derivated_Unk(const Int_t index_m, const TMatrixD &u)
{
    switch(index_m)
    {
    case 0:
        return TVector3(1, 0, 0);
    case 1:
        return TVector3(0, 1, 0);
    case 2:
        return TVector3(0, 0, 1);
    default:
        return TVector3(0, 0, 0);
    }
}

TVector3        GKinFitter::GetVertex_Derivated_Unk(const Int_t index, const Int_t index_m, const TMatrixD &u)
{
    switch(index)
    {
    case 0:
        switch(index_m)
        {
        case 3:
            return TVector3(1, 0, 0);
        case 4:
            return TVector3(0, 1, 0);
        case 5:
            return TVector3(0, 0, 1);
        default:
            return TVector3(0, 0, 0);
        }
    case 1:
        switch(index_m)
        {
        case 6:
            return TVector3(1, 0, 0);
        case 7:
            return TVector3(0, 1, 0);
        case 8:
            return TVector3(0, 0, 1);
        default:
            return TVector3(0, 0, 0);
        }
    case 2:
        switch(index_m)
        {
        case 9:
            return TVector3(1, 0, 0);
        case 10:
            return TVector3(0, 1, 0);
        case 11:
            return TVector3(0, 0, 1);
        default:
            return TVector3(0, 0, 0);
        }
    default:
        return TVector3(0, 0, 0);
    }
}

TLorentzVector  GKinFitter::GetBeam(const TMatrixD x, const TMatrixD u)
{
    TVector3    beamPos(u[12][0], u[13][0], -GKinFitter_CBRadius);
    beamPos     = GetMainVertex(u)-beamPos;
    beamPos     = x[0][0] * beamPos * (1/beamPos.Mag());
    return TLorentzVector(beamPos, x[0][0]);
}

TLorentzVector  GKinFitter::GetBeam_Derivated_Par(const Int_t index_n, const TMatrixD x, const TMatrixD u)
{
    if(index_n==0 || index_n==1 || index_n==2)
    {
        TVector3    beamPos(u[12][0], u[13][0], -GKinFitter_CBRadius);
        beamPos     = GetMainVertex(u)-beamPos;
        beamPos     = beamPos * (1/beamPos.Mag());
        return TLorentzVector(beamPos, 1);
    }
    else
        return TLorentzVector(0, 0, 0, 0);
}

TLorentzVector  GKinFitter::GetBeam_Derivated_Unk(const Int_t index_m, const TMatrixD x, const TMatrixD u)
{
    if(index_n==12)
    {
        TVector3    beamPos(u[12][0], u[13][0], -GKinFitter_CBRadius);
        beamPos     = GetMainVertex(u)-beamPos;
        beamPos     = beamPos * (1/beamPos.Mag());
        return TLorentzVector(beamPos, 1);
    }
    else
        return TLorentzVector(0, 0, 0, 0);
}

TLorentzVector  GKinFitter::GetParticle(const Int_t index, const TMatrixD x, const TMatrixD u)
{
    TVector3    beamPos(x[(index*4)+1][0], x[(index*4)+2][0], x[(index*4)+3][0]);
    switch(index)
    {
    case 0:
    case 1:
        beamPos     = GetVertex(0, u)-beamPos;
        break;
    case 2:
    case 3:
        beamPos     = GetVertex(1, u)-beamPos;
        break;
    case 4:
    case 5:
        beamPos     = GetVertex(2, u)-beamPos;
        break;
    }
    beamPos     = x[(index*4)+4][0] * beamPos * (1/beamPos.Mag());
    return TLorentzVector(beamPos, x[(index*4)+4][0]);
}

Double_t    GKinFitter::Get_g(const Int_t index_k, const TMatrixD x, const TMatrixD u)
{
    switch(index_k)
    {
    case 0:                     //inv. Mass eta
        return (GetParticle(0, x, u)+GetParticle(1, x, u)).M2()-(MASS_ETA*MASS_ETA);
    case 1:                     //inv. Mass pi0 a
        return (GetParticle(2, x, u)+GetParticle(3, x, u)).M2()-(MASS_PI0*MASS_PI0);
    case 2:                     //inv. Mass pi0 b
        return (GetParticle(4, x, u)+GetParticle(5, x, u)).M2()-(MASS_PI0*MASS_PI0);
    case 3:                     //mis. Mass proton
        return (GetBeam(x, u)-GetParticle(0, x, u)-GetParticle(1, x, u)-GetParticle(2, x, u)-GetParticle(3, x, u)-GetParticle(4, x, u)-GetParticle(5, x, u)).M2()-(MASS_PROTON*MASS_PROTON);
    }
}

TMatrixD        GKinFitter::Get_g_Derivated_Par(const Int_t index_k, const TMatrixD x, const TMatrixD u)
{
    TMatrixD    ret(nPar, 1);

    switch(index_k)
    {
    case 0:                     //inv. Mass eta
        {
            TLorentzVector vec0(GetParticle(0, x, u));
            TLorentzVector vec1(GetParticle(1, x, u));
            ret[0][0]   = 0;
            ret[1][0]   = -2*(vec0.Px()+vec1.Px());
            ret[2][0]   = -2*(vec0.Py()+vec1.Py());
            ret[3][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[4][0]   = 2*(vec0.E()+vec1.E());
            ret[5][0]   = -2*(vec0.Px()+vec1.Px());
            ret[6][0]   = -2*(vec0.Py()+vec1.Py());
            ret[7][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[8][0]   = 2*(vec0.E()+vec1.E());
            for(int i=9; i<nPar; i++)
                ret[i][0]   = 0;
        }
        return ret;
    case 1:                     //inv. Mass pi0 a
        {
            TLorentzVector vec0(GetParticle(2, x, u));
            TLorentzVector vec1(GetParticle(3, x, u));
            for(int i=0; i<9; i++)
                ret[i][0]   = 0;
            ret[9][0]   = -2*(vec0.Px()+vec1.Px());
            ret[10][0]   = -2*(vec0.Py()+vec1.Py());
            ret[11][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[12][0]   = 2*(vec0.E()+vec1.E());
            ret[13][0]   = -2*(vec0.Px()+vec1.Px());
            ret[14][0]   = -2*(vec0.Py()+vec1.Py());
            ret[15][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[16][0]   = 2*(vec0.E()+vec1.E());
            for(int i=17; i<nPar; i++)
                ret[i][0]   = 0;
        }
        return ret;
    case 2:                     //inv. Mass pi0 b
        {
            TLorentzVector vec0(GetParticle(4, x, u));
            TLorentzVector vec1(GetParticle(5, x, u));
            for(int i=0; i<17; i++)
                ret[i][0]   = 0;
            ret[17][0]   = -2*(vec0.Px()+vec1.Px());
            ret[18][0]   = -2*(vec0.Py()+vec1.Py());
            ret[19][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[20][0]   = 2*(vec0.E()+vec1.E());
            ret[21][0]   = -2*(vec0.Px()+vec1.Px());
            ret[22][0]   = -2*(vec0.Py()+vec1.Py());
            ret[23][0]   = -2*(vec0.Pz()+vec1.Pz());
            ret[24][0]   = 2*(vec0.E()+vec1.E());
        }
        return ret;
    case 3:                     //mis. Mass proton
    {
        TLorentzVector beam(GetBeam(x, u));
        TLorentzVector vec0(GetParticle(0, x, u));
        TLorentzVector vec1(GetParticle(1, x, u));
        TLorentzVector vec2(GetParticle(2, x, u));
        TLorentzVector vec3(GetParticle(3, x, u));
        TLorentzVector vec4(GetParticle(4, x, u));
        TLorentzVector vec5(GetParticle(5, x, u));
        ret[0][0]   = 2*(beam.E()-vec0.E()-vec1.E()-vec2.E()-vec3.E()-vec4.E()-vec5.E())-2*(beam.Pz()-vec0.Pz()-vec1.Pz()-vec2.Pz()-vec3.Pz()-vec4.Pz()-vec5.Pz());
        ret[1][0]   = -2*(beam.Px()-vec0.Px()-vec1.Px()-vec2.Px()-vec3.Px()-vec4.Px()-vec5.Px());
        ret[2][0]   = -2*(beam.Py()-vec0.Py()-vec1.Py()-vec2.Py()-vec3.Py()-vec4.Py()-vec5.Py());
        ret[3][0]   = -2*(beam.Pz()-vec0.Pz()-vec1.Pz()-vec2.Pz()-vec3.Pz()-vec4.Pz()-vec5.Pz());
        ret[4][0]   = 2*(vec0.E()+vec1.E());
        ret[5][0]   = -2*(vec0.Px()+vec1.Px());
        ret[6][0]   = -2*(vec0.Py()+vec1.Py());
        ret[7][0]   = -2*(vec0.Pz()+vec1.Pz());
        ret[8][0]   = 2*(vec0.E()+vec1.E());
    }
    return ret;
    }
}

TMatrixD        GKinFitter::Get_g_Derivated_Unk(const Int_t index_k, const TMatrixD x, const TMatrixD u)
{

}*/




void GKinFitter::Constraints()
{
    /*
  //Add up the particle 4 vectors
    TLorentzVector eta(fmAlpha0[16][0], fmAlpha0[17][0], fmAlpha0[18][0], fmAlpha0[19][0]);
    eta     += TLorentzVector(fmAlpha0[20][0], fmAlpha0[21][0], fmAlpha0[22][0], fmAlpha0[23][0]);
    TLorentzVector pi0a(fmAlpha0[24][0], fmAlpha0[25][0], fmAlpha0[26][0], fmAlpha0[27][0]);
    pi0a    += TLorentzVector(fmAlpha0[28][0], fmAlpha0[29][0], fmAlpha0[30][0], fmAlpha0[31][0]);
    TLorentzVector pi0b(fmAlpha0[32][0], fmAlpha0[33][0], fmAlpha0[34][0], fmAlpha0[35][0]);
    pi0b    += TLorentzVector(fmAlpha0[36][0], fmAlpha0[37][0], fmAlpha0[38][0], fmAlpha0[39][0]);
    TLorentzVector all(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]);
    all     -= eta;
    all     -= pi0a;
    all     -= pi0b;*/
    /*all.Print();
    (all-TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0])).Print();
    TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]).Print();
    TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]).Print();
    std::cout << all.M() << "   " << (all-TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0])).M() << "   " << TLorentzVector(fmAlpha0[12][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]).M() << "   " << std::endl;
    std::cout << eta.M() << "   " << pi0a.M() << "   " << pi0b.M() << "   " << std::endl;
    std::cout << TLorentzVector(fmAlpha0[16][0], fmAlpha0[17][0], fmAlpha0[18][0], fmAlpha0[19][0]).M() << "   " << TLorentzVector(fmAlpha0[20][0], fmAlpha0[21][0], fmAlpha0[22][0], fmAlpha0[23][0]).M() << "   " << TLorentzVector(fmAlpha0[24][0], fmAlpha0[25][0], fmAlpha0[26][0], fmAlpha0[27][0]).M() << "   " << std::endl;
*//*
  //d matrix (evaluate constraint eqn.)
    fmd[0][0]   =   eta.M2() - (MASS_ETA*MASS_ETA);
    fmd[1][0]   =   pi0a.M2() - (MASS_PI0*MASS_PI0);
    fmd[2][0]   =   pi0b.M2() - (MASS_PI0*MASS_PI0);
    fmd[3][0]   =   all.M2() - (MASS_PROTON*MASS_PROTON);

  //D matrix (derivitives of constraint eqn)
    for(int i=0; i<4; i++)
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

    fmD[3][12] =-2*all.X();
    fmD[3][13] =-2*all.Y();
    fmD[3][14] =-2*all.Z();
    fmD[3][15] = 2*all.T();
    for(int i=0; i<6; i++)
    {
        fmD[3][4*i+16] =-2*all.X();
        fmD[3][4*i+17] =-2*all.Y();
        fmD[3][4*i+18] =-2*all.Z();
        fmD[3][4*i+19] = 2*all.T();
    }*/
}

/*TLorentzVector  GKinFitter::GetEta()
{
    TLorentzVector ret(fmAlpha[16][0], fmAlpha[17][0], fmAlpha[18][0], fmAlpha[19][0]);
        ret +=  TLorentzVector(fmAlpha[20][0], fmAlpha[21][0], fmAlpha[22][0], fmAlpha[23][0]);
    return ret;
}*/

/*TLorentzVector  GKinFitter::GetEtap()
{
    //TLorentzVector(fmAlpha[12][0], fmAlpha[13][0], fmAlpha[14][0], fmAlpha[15][0]).Print();
    //TLorentzVector(fmAlpha0[21][0], fmAlpha0[13][0], fmAlpha0[14][0], fmAlpha0[15][0]).Print();

    TLorentzVector ret(0.0, 0.0, 0.0, 0.0);
    for(int i=0; i<6; i++)
    {
        ret +=  TLorentzVector(fmAlpha[4*i+16][0], fmAlpha[4*i+17][0], fmAlpha[4*i+18][0], fmAlpha[4*i+19][0]);
        //TLorentzVector(fmAlpha[4*i+16][0], fmAlpha[4*i+17][0], fmAlpha[4*i+18][0], fmAlpha[4*i+19][0]).Print();
        //TLorentzVector(fmAlpha0[4*i+16][0], fmAlpha0[4*i+17][0], fmAlpha0[4*i+18][0], fmAlpha0[4*i+19][0]).Print();
    }
    //ret.Print();
    //std::cout << std::endl;
    return ret;
}

TLorentzVector  GKinFitter::GetPi0a()
{
    TLorentzVector ret(fmAlpha[24][0], fmAlpha[25][0], fmAlpha[26][0], fmAlpha[27][0]);
        ret +=  TLorentzVector(fmAlpha[28][0], fmAlpha[29][0], fmAlpha[30][0], fmAlpha[31][0]);
    return ret;
}

TLorentzVector  GKinFitter::GetPi0b()
{
    TLorentzVector ret(fmAlpha[32][0], fmAlpha[33][0], fmAlpha[34][0], fmAlpha[35][0]);
        ret +=  TLorentzVector(fmAlpha[36][0], fmAlpha[37][0], fmAlpha[38][0], fmAlpha[39][0]);
    return ret;
}*/

void    GKinFitter::Set(const GKinFitterPolarToCartesian& v)
{
    par0    = v.GetParametersH();
    V0      = v.GetCovarianceH();

    for(int i=0; i<nUnk; i++)
        unk[i]  = 0;
}

Bool_t  GKinFitter::SolveStep(const TMatrixD& x, const TMatrixD& u, TMatrixD& new_x, TMatrixD& new_u, Double_t& chiSq)
{
    TMatrixD    g_Derivated_Par(4, 25);
    Get_g_Derivated_Par(x, u, g_Derivated_Par);
    //g_Derivated_Par.Print();
    TMatrixD    g_Derivated_Par_T(g_Derivated_Par);
    g_Derivated_Par_T.T();
    //g_Derivated_Par_T.Print();
    TMatrixD    S_Inverted(4, 4);
    Get_S(g_Derivated_Par, S_Inverted);
    //S_Inverted.Print();
    Double_t    determinant;
    S_Inverted.Invert(&determinant);
    if(determinant==0)
    {
        std::cout << "Can not invert S in KinFitter::SolveStep" << std::endl;
        return kFALSE;
    }
    //S_Inverted.Print();

    TMatrixD    g_Derivated_Unk(4, 5);
    Get_g_Derivated_Unk(x, u, g_Derivated_Unk);
    TMatrixD    g_Derivated_Unk_T(g_Derivated_Unk);
    g_Derivated_Unk_T.T();

    TMatrixD    help(g_Derivated_Unk_T * S_Inverted * g_Derivated_Unk);
    //help.Print();
    help.Invert(&determinant);
    if(determinant==0)
    {
        std::cout << "Can not invert help in KinFitter::SolveStep" << std::endl;
        return kFALSE;
    }

    TMatrixD    g(4, 1);
    TMatrixD    r(4, 1);
    Get_g(x, u, g);
    Get_r(x, u, g_Derivated_Par, g, r);

    new_u   = u - (help * g_Derivated_Unk_T * S_Inverted * r);
    lambda  = S_Inverted * (r + (g_Derivated_Unk * (new_u - u)));
    new_x   = par0 - (V0 * g_Derivated_Par_T * lambda);

    chiSq   = Calc_Chi2(new_x, g);

    return kTRUE;
}

Bool_t  GKinFitter::Solve()
{
    if(SolveStep(par0, unk0, par, unk, chi2)==kFALSE)
    {
        std::cout << "first SolveStep not working in KinFitter::Solve" << std::endl;
        return kFALSE;
    }

    TMatrixD    new_x((6*4)+1, 1);
    TMatrixD    new_u(5, 1);
    Double_t    newChi2;

    for(int i=0; i<20; i++)
    {
        if(SolveStep(par, unk, new_x, new_u, newChi2)==kFALSE)
        {
            std::cout << i << " SolveStep not working in KinFitter::Solve" << std::endl;
            return kFALSE;
        }
        if(newChi2<chi2)
        {
            chi2    = newChi2;
            par     = new_x;
            unk     = new_u;
        }
        else
            return kTRUE;
    }

    std::cout << "No solution after 20 steps in KinFitter::Solve" << std::endl;
    return kFALSE;
}
