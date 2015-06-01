#include "GoAT_Detailed.h"


GoAT_Detailed::GoAT_Detailed() :
    UseParticleReconstruction(0),
    nEvents_written(0),
    AcceptanceTrue("AcceptanceTrue", "AcceptanceTrue", 180, 0, 180, 48)
{ 
}

GoAT_Detailed::~GoAT_Detailed()
{
}

Bool_t	GoAT_Detailed::Init(const char* configfile)
{
    cout << endl << "Initialising GoAT_Detailed analysis..." << endl << endl;

    if(configfile)
        SetConfigFile(configfile);
    std::string config = ReadConfig("Period-Macro");
	if( sscanf(config.c_str(),"%d\n", &period) == 1 ) UsePeriodMacro = 1;

	cout << "==========================================================" << endl;	
	cout << "Setting up Data Checks:" << endl;	
	cout << "==========================================================" << endl;	
    if(!GDataChecks::Init())
	{
		cout << "GDataChecks Init failed!" << endl; 
		return kFALSE;
	}

	cout << "==========================================================" << endl;	
	cout << "Setting up sorting criteria:" << endl;	
	cout << "==========================================================" << endl;	
    if(!GSort::Init())
	{
		cout << "GSort Init failed!" << endl; 
		return kFALSE;
	}
	
	cout << "==========================================================" << endl;		
	cout << "Setting up analysis classes:" << endl;	
	cout << "==========================================================" << endl;	
    config = ReadConfig("DO-PARTICLE-RECONSTRUCTION");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        int buffer=0;
        sscanf( config.c_str(), "%d\n", &buffer);
        UseParticleReconstruction = (buffer==1);
    }

	if(UseParticleReconstruction) 
	{
        if(!GParticleReconstruction::Init())
		{
			cout << "GParticleReconstruction Init failed!" << endl; 
			return kFALSE;
		}
    }

    config = ReadConfig("DO-MESON-RECONSTRUCTION");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        int buffer=0;
        sscanf( config.c_str(), "%d\n", &buffer);
        UseMesonReconstruction = (buffer==1);
    }

    if(UseMesonReconstruction)
    {
        if(!GMesonReconstruction_6and7gamma::Init())
        {
            cout << "GParticleReconstruction Init failed!" << endl;
            return kFALSE;
        }
    }

	cout << endl;	

	cout << "Initialisation complete." << endl;
	cout << "==========================================================" << endl << endl;
  
	return kTRUE;
}

void	GoAT_Detailed::ProcessEvent()
{
    if(UsePeriodMacro == 1)
    {
        if(GetEventNumber() % period == 0)
            cout << "Event: " << GetEventNumber() << "  Events Accepted: " << nEvents_written << endl;
    }

    if(SortAnalyseEvent())
    {
        if(UseParticleReconstruction)
        {
            if(UseMesonReconstruction)
            {
                if(!GParticleReconstruction::ProcessEventWithoutFilling())  return;
                if(!GMesonReconstruction_6and7gamma::ProcessEventWithoutFilling())  return;
                if(!SortFillEvent())    return;
                GetElectrons()->Fill();
                GetProtons()->Fill();
                GetNeutrons()->Fill();
                GetNeutralPions()->Fill();
                GetEtas()->Fill();
                GetEtaPrimes()->Fill();
            }
            else
            {
                if(!GParticleReconstruction::ProcessEventWithoutFilling())  return;
                if(!SortFillEvent())    return;
                GetElectrons()->Fill();
                GetProtons()->Fill();
                GetNeutrons()->Fill();
            }
        }
        else if(UseMesonReconstruction)
        {
            GMesonReconstruction_6and7gamma::ProcessEventWithoutFilling();
            if(!SortFillEvent())    return;
            GetNeutralPions()->Fill();
            GetEtas()->Fill();
            GetEtaPrimes()->Fill();
        }
        GetEventParameters()->SetNReconstructed(GetNReconstructed());
        GetEventParameters()->Fill();
        GetRootinos()->Fill();
        GetPhotons()->Fill();
        GetChargedPions()->Fill();
        FillReadList();
        nEvents_written++;

        for(int i=0; i<GetTagger()->GetNTagged(); i++)
        {
            TLorentzVector  helpCM(GetGeant()->GetTrueVector(2));
            helpCM.Boost(-GetTagger()->GetVectorProtonTarget(i).BoostVector());
            AcceptanceTrue.Fill(helpCM.Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
        }
    }
}

Bool_t	GoAT_Detailed::Start()
{
    if(!IsAcquFile())
    {
        cout << "ERROR: Input File is not a Acqu file." << endl;
        return kFALSE;
    }
    SetAsGoATFile();

    if(UseParticleReconstruction)
    {
        if(UseMesonReconstruction)
        {
            GetRootinos()->CloseForInput();
            GetPhotons()->CloseForInput();
            GetElectrons()->CloseForInput();
            GetChargedPions()->CloseForInput();
            GetProtons()->CloseForInput();
            GetNeutrons()->CloseForInput();
            GetNeutralPions()->CloseForInput();
            GetEtas()->CloseForInput();
            GetEtaPrimes()->CloseForInput();
        }
        else
        {
            GetRootinos()->CloseForInput();
            GetPhotons()->CloseForInput();
            GetElectrons()->CloseForInput();
            GetChargedPions()->CloseForInput();
            GetProtons()->CloseForInput();
            GetNeutrons()->CloseForInput();
        }
    }
    else if(UseMesonReconstruction)
    {
        GetNeutralPions()->CloseForInput();
        GetEtas()->CloseForInput();
        GetEtaPrimes()->CloseForInput();
    }

    if(!TraverseValidEvents())		return kFALSE;

    return kTRUE;
}
