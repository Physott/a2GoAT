#include "MyPhysics.h"



MyPhysics::MyPhysics()    :
    hist_eta("eta", "eta", kFALSE),
    hist_etap("etap", "etap", kTRUE)
{ 
        GHistBGSub::InitCuts(-20, 20, -55, -35);
        GHistBGSub::AddRandCut(35, 55);
}

MyPhysics::~MyPhysics()
{

}

Bool_t	MyPhysics::Start()
{
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
    if(eta->GetNParticles()>0)
        hist_eta.Fill(*eta, *protons, *tagger, kTRUE);
    //if(etap->GetNParticles()>0)
        //hist_etap.Fill(*etap, *protons, *tagger, kTRUE);
}

void	MyPhysics::ProcessScalerRead()
{
    hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
    hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}


Bool_t	MyPhysics::Init(const char* configfile)
{
    SetConfigFile(configfile);
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
    }

    return kTRUE;
}
