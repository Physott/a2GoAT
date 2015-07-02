#include "MyPhysics.h"



MyPhysics::MyPhysics()  :
    proton("proton", "proton"),
    etap("etap", "etap"),
    etaPhotons("etaPhotons", "etaPhotons"),
    pi0Photons("pi0Photons", "pi0Photons"),
    allPhotons("allPhotons", "allPhotons")
{ 
        GHistBGSub::InitCuts(-20, 20, -535, -35);
        GHistBGSub::AddRandCut(35, 535);
}

MyPhysics::~MyPhysics()
{

}

Bool_t	MyPhysics::Start()
{
    if(!IsPhysicsFile())
    {
        cout << "ERROR: Input File is not a Physics file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TObject*    obj;
    obj =    GetInputFileRef().Get("EPT_Scaler");  if(obj){outputFile->cd();   obj->Write();} else std::cout << "Can not find object named 'EPT_Scaler'" << std::endl;
    obj =    GetInputFileRef().Get("EPT_ScalerCor");  if(obj){outputFile->cd();   obj->Write();} else std::cout << "Can not find object named 'EPT_ScalerCor'" << std::endl;
    obj =    GetInputFileRef().Get("EPT_ScalerT");  if(obj){outputFile->cd();   obj->Write();} else std::cout << "Can not find object named 'EPT_ScalerT'" << std::endl;
    obj =    GetInputFileRef().Get("EPT_ScalerCorT");  if(obj){outputFile->cd();   obj->Write();} else std::cout << "Can not find object named 'EPT_ScalerCorT'" << std::endl;
    obj =    GetInputFileRef().Get("AcceptanceTrue");  if(obj){outputFile->cd();   obj->Write();} else std::cout << "Can not find object named 'AcceptanceTrue'" << std::endl;
    obj =    GetInputFileRef().Get("AcceptanceProtonTrue");  if(obj){outputFile->cd();   obj->Write();} else std::cout << "Can not find object named 'AcceptanceProtonTrue'" << std::endl;

    TraverseValidEvents();

	return kTRUE;
}

void	MyPhysics::ProcessEvent()
{
    for(int i=0; i<GetTagger()->GetNTagged(); i++)
    {

        etap.Fill(*GetEtaPrimes(), 0, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        if(GetProtons()->GetNParticles()>0)
            proton.Fill(*GetProtons(), 0, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        etaPhotons.Fill(*GetPhotons(), 0, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        etaPhotons.Fill(*GetPhotons(), 1, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        pi0Photons.Fill(*GetPhotons(), 2, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        pi0Photons.Fill(*GetPhotons(), 3, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        pi0Photons.Fill(*GetPhotons(), 4, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        pi0Photons.Fill(*GetPhotons(), 5, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        allPhotons.Fill(*GetPhotons(), 0, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        allPhotons.Fill(*GetPhotons(), 1, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        allPhotons.Fill(*GetPhotons(), 2, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        allPhotons.Fill(*GetPhotons(), 3, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        allPhotons.Fill(*GetPhotons(), 4, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        allPhotons.Fill(*GetPhotons(), 5, GetTagger()->GetTaggedEnergy(i), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
    }
}

void	MyPhysics::ProcessScalerRead()
{

}


Bool_t	MyPhysics::Init(const char* configfile)
{
    SetConfigFile(configfile);
    //std::string config;
    //Double_t    buf[8];

    /*config = ReadConfig("Cut-Etap-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf %lf %lf %lf %lf\n", &buf[0], &buf[1], &buf[2], &buf[3], &buf[4], &buf[5]) == 6)
        {
            config = ReadConfig("Cut-Etap-MM");
            if (strcmp(config.c_str(), "nokey") != 0)
            {
                if(sscanf( config.c_str(), "%lf %lf\n", &buf[6], &buf[7]) == 2)
                {
                    hist_etap.SetCutSubIM(0, buf[0], buf[1]);
                    hist_etap.SetCutSubIM(1, buf[2], buf[3]);
                    hist_etap.SetCutSubIM(2, buf[4], buf[5]);
                    hist_etap.SetCutMM(buf[6], buf[7]);
                    cout << "Set Cuts for etap physics: ";
                    for(int i=0; i<8; i++)
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
                    hist_etap_proton.SetCutSubIM(0, buf[0], buf[1]);
                    hist_etap_proton.SetCutSubIM(1, buf[2], buf[3]);
                    hist_etap_proton.SetCutSubIM(2, buf[4], buf[5]);
                    hist_etap_proton.SetCutMM(buf[6], buf[7]);
                    cout << "Set Cuts for etap with proton physics: ";
                    for(int i=0; i<8; i++)
                        cout << buf[i] << "   ";
                    cout << endl;
                }
            }
        }
    }*/
    return kTRUE;
}
