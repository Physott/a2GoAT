#include "ScaleMC.h"



ScaleMC::ScaleMC()  :
    CalibCB("CalibCB", "CalibCB", 800, 0, 800, 720, 0, 720),
    CalibTAPS("CalibTAPS", "CalibTAPS", 800, 0, 800, 428, 0, 428),
    CalibCBCorr("CalibCBCorr", "CalibCBCorr", 800, 0, 800, 720, 0, 720),
    CalibTAPSCorr("CalibTAPSCorr", "CalibTAPSCorr", 800, 0, 800, 428, 0, 428)
{
        GHistBGSub::InitCuts(-20, 20, -535, -35);
        GHistBGSub::AddRandCut(35, 535);
}

ScaleMC::~ScaleMC()
{

}

Bool_t	ScaleMC::Start()
{
    if(!IsAcquFile())
    {
        cout << "ERROR: Input File is not a Acqu file." << endl;
        return kFALSE;
    }

    TraverseValidEvents();

    return kTRUE;
}

void	ScaleMC::ProcessEvent()
{
    for(int i=0; i<GetTracks()->GetNTracks(); i++)
    {
        for(int j=i+1; j<GetTracks()->GetNTracks(); j++)
        {
            if(GetTracks()->HasCB(i) == kTRUE && GetTracks()->HasCB(j) == kTRUE)
            {
                CalibCB.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTracks()->GetCentralCrystal(i));
                CalibCB.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTracks()->GetCentralCrystal(j));
            }
            else
            {
                if(GetTracks()->HasCB(i) == kFALSE)
                    CalibCB.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTracks()->GetCentralCrystal(i));
                if(GetTracks()->HasCB(j) == kFALSE)
                    CalibCB.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTracks()->GetCentralCrystal(j));
            }
        }
    }


    for(int i=0; i<GetTracks()->GetNTracks(); i++)
        GetTracks()->SetClusterEnergy(i, GetTracks()->GetClusterEnergy(i)*(1.11+(-0.0001*GetTracks()->GetClusterEnergy(i))));


    for(int i=0; i<GetTracks()->GetNTracks(); i++)
    {
        for(int j=i+1; j<GetTracks()->GetNTracks(); j++)
        {
            if(GetTracks()->HasCB(i) == kTRUE && GetTracks()->HasCB(j) == kTRUE)
            {
                CalibCBCorr.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTracks()->GetCentralCrystal(i));
                CalibCBCorr.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTracks()->GetCentralCrystal(j));
            }
            else
            {
                if(GetTracks()->HasCB(i) == kFALSE)
                    CalibCBCorr.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTracks()->GetCentralCrystal(i));
                if(GetTracks()->HasCB(j) == kFALSE)
                    CalibCBCorr.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTracks()->GetCentralCrystal(j));
            }
        }
    }

    GetTracks()->Fill();
    GetTagger()->Fill();
    GetDetectorHits()->Fill();
    GetTrigger()->Fill();
}

void	ScaleMC::ProcessScalerRead()
{
}


Bool_t	ScaleMC::Init(const char* configfile)
{
    return kTRUE;
}
