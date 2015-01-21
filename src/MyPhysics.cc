#include "MyPhysics.h"



MyPhysics::MyPhysics()    :
    //hist_eta("eta", "eta", kTRUE),
    //hist_eta_proton("eta_proton", "eta_proton", kTRUE),
    //hist_etap("etap", "etap", kTRUE),
    hist_etap_proton("etap_proton", "etap_proton", kTRUE)
{ 
        GHistBGSub::InitCuts(-20, 20, -55, -35);
        GHistBGSub::AddRandCut(35, 55);
}

MyPhysics::~MyPhysics()
{

}

Bool_t	MyPhysics::Start()
{
    TLorentzVector  g0(0.0, 134.98/2, 0.0, 134.98/2);
    TLorentzVector  g1(0.0, -134.98/2, 0.0, 134.98/2);
    TLorentzVector  pi(g0+g1);
    //pi.Print();

    g0.Boost(0.0, 0.3, 0.9);
    g1.Boost(0.0, 0.3, 0.9);
    pi.Boost(0.0, 0.3, 0.9);
    //g0.Print();
    //g1.Print();
    //pi.Print();

    TLorentzVector  p(-pi.Px(), -pi.Py(), -pi.Pz(), sqrt((pi.P()*pi.P())+(938.27*938.27)));
    //p.Print();

    TLorentzVector  b(pi+p);
    //b.Print();
    g0.Boost(0.40902, 0.0, 0.0);
    g1.Boost(0.40902, 0.0, 0.0);
    pi.Boost(0.40902, 0.0, 0.0);
    b.Boost(0.40902, 0.0, 0.0);
    p.Boost(0.40902, 0.0, 0.0);
    /*b.Print();
    g0.Print();
    g1.Print();
    pi.Print();
    p.Print();*/
    std::cout << b.E()-b.P() << std::endl;

    GKinFitterBase  fitter(2,2);
    GHistFit2       hfit("test", "tset", kTRUE);

    Double_t        bsmear      = b.Px()+(0.05*b.Px());
    TLorentzVector  g0smear(g0.E()-(0.05*g0.E()), 0.0, 0.0, g0.E()-(0.05*g0.E()));
                    g0smear.SetTheta(g0.Theta()+(0.05*g0.Theta()));
                    g0smear.SetPhi(g0.Phi()-(0.05*g0.Phi()));
    TLorentzVector  g1smear(g1.E()+(0.05*g1.E()), 0.0, 0.0, g1.E()+(0.05*g1.E()));
                    g1smear.SetTheta(g1.Theta()-(0.05*g1.Theta()));
                    g1smear.SetPhi(g1.Phi()+(0.05*g1.Phi()));
    std::cout << bsmear << std::endl;
    g0smear.Print();
    g1smear.Print();

    fitter.AddBeam(b.Px()+(0.05*b.Px()), 938.27, 0.05*b.Px(), 0.005);
    fitter.AddGamma(g0.E()-(0.05*g0.E()), g0.Theta()+(0.05*g0.Theta()), g0.Phi()-(0.05*g0.Phi()), 0.05*g0.E(), 0.05*g0.Theta(), 0.05*g0.Phi());
    fitter.AddGamma(g1.E()+(0.05*g1.E()), g1.Theta()-(0.05*g1.Theta()), g1.Phi()+(0.05*g1.Phi()), 0.05*g1.E(), 0.05*g1.Theta(), 0.05*g1.Phi());
    int indices[2];
    indices[0]  = 0;
    indices[1]  = 1;
    fitter.AddInvMassConstraint(indices, 2, 134.98);
    fitter.AddMisMassConstraint(938.27);

    fitter.GetInitialBeam().Print();
    fitter.GetBeam().Print();
    fitter.GetInitialPhoton(0).Print();
    fitter.GetPhoton(0).Print();
    fitter.GetInitialPhoton(1).Print();
    fitter.GetPhoton(1).Print();
    std::cout << fitter.GetInitialIMConstraint(0) << std::endl;
    std::cout << fitter.GetInitialMMConstraint() << std::endl;

    TMatrixD    gDerPar(2, 9);
    fitter.GetInitialGDerPar(gDerPar);
    gDerPar.Print();

    fitter.Solve(hfit);

    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();


    TraverseValidEvents();


	return kTRUE;
}

void	MyPhysics::ProcessEvent()
{
    /*if(eta->GetNParticles()>0)
    {
        hist_eta.Fill(*eta, *tagger, kTRUE);
        if(protons->GetNParticles()>0)
            hist_eta_proton.Fill(*eta, *protons, *tagger, kTRUE);
    }*/
    if(etap->GetNParticles()>0)
    {
        if(protons->GetNParticles()>0)
            hist_etap_proton.Fill(*etap, *protons, *tagger, kTRUE);
        //else
            //hist_etap.Fill(*etap, *tagger, kTRUE);
    }
}

void	MyPhysics::ProcessScalerRead()
{
    //hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
    //hist_eta_proton.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
    //hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
    //hist_etap_proton.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}


Bool_t	MyPhysics::Init(const char* configfile)
{
    /*SetConfigFile(configfile);
    Double_t    buf[8];
    std::string config = ReadConfig("Cut-Eta-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf %lf %lf %lf %lf\n", &buf[0], &buf[1], &buf[2], &buf[3], &buf[4], &buf[5]) == 6)
        {
            config = ReadConfig("Cut-Eta-MM");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[6], &buf[7]) == 2)
                {
                    hist_eta.SetHistMeson(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
                    cout << "Set Cuts for eta physics: ";
                    for(int i=0; i<8; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }
    config = ReadConfig("Cut-Eta-ConfidenceLevel");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            hist_eta.SetFitMeson(buf[0], buf[1]);
            cout << "Set Cuts for eta fit: " << buf[0] << "   " << buf[1] << endl;
        }
    }


    config = ReadConfig("Cut-Eta-Proton-MinAngle");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf \n", &buf[0]) == 1)
        {
            config = ReadConfig("Cut-Eta-Proton-Coplanarity");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[1], &buf[2]) == 2)
                {
                    hist_eta.SetCheckProton(buf[0], buf[1], buf[2]);
                    cout << "Set Cuts for proton checking in eta data: ";
                    for(int i=0; i<3; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }

    config = ReadConfig("Cut-Eta-Proton-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf %lf %lf %lf %lf\n", &buf[0], &buf[1], &buf[2], &buf[3], &buf[4], &buf[5]) == 6)
        {
            config = ReadConfig("Cut-Eta-Proton-MM");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[6], &buf[7]) == 2)
                {
                    hist_eta.SetHistMesonProton(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
                    cout << "Set Cuts for eta with proton physics: ";
                    for(int i=0; i<8; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }
    config = ReadConfig("Cut-Eta-Proton-ConfidenceLevel");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            hist_eta.SetFitMesonProton(buf[0], buf[1]);
            cout << "Set Cuts for eta proton fit: " << buf[0] << "   " << buf[1] << endl;
        }
    }






    config = ReadConfig("Cut-Etap-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf %lf %lf %lf %lf\n", &buf[0], &buf[1], &buf[2], &buf[3], &buf[4], &buf[5]) == 6)
        {
            config = ReadConfig("Cut-Etap-MM");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[6], &buf[7]) == 2)
                {
                    hist_etap.SetHistMeson(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
                    cout << "Set Cuts for etap physics: ";
                    for(int i=0; i<8; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }
    config = ReadConfig("Cut-Etap-ConfidenceLevel");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            hist_etap.SetFitMeson(buf[0], buf[1]);
            cout << "Set Cuts for etap fit: " << buf[0] << "   " << buf[1] << endl;
        }
    }

    config = ReadConfig("Cut-Etap-Proton-MinAngle");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf \n", &buf[0]) == 1)
        {
            config = ReadConfig("Cut-Etap-Proton-Coplanarity");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[1], &buf[2]) == 2)
                {
                    hist_etap.SetCheckProton(buf[0], buf[1], buf[2]);
                    cout << "Set Cuts for proton checking in etap data: ";
                    for(int i=0; i<3; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }

    config = ReadConfig("Cut-Etap-Proton-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf %lf %lf %lf %lf\n", &buf[0], &buf[1], &buf[2], &buf[3], &buf[4], &buf[5]) == 6)
        {
            config = ReadConfig("Cut-Etap-Proton-MM");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[6], &buf[7]) == 2)
                {
                    hist_etap.SetHistMesonProton(buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
                    cout << "Set Cuts for etap with proton physics: ";
                    for(int i=0; i<8; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }
    config = ReadConfig("Cut-Etap-Proton-ConfidenceLevel");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            hist_etap.SetFitMesonProton(buf[0], buf[1]);
            cout << "Set Cuts for etap proton fit: " << buf[0] << "   " << buf[1] << endl;
        }
    }*/

    return kTRUE;
}
