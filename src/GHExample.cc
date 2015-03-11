#include "GHExample.h"

GHExample::GHExample() :
    taggerTime("taggerTime", "taggerTime", 1600, -800, 800, 48),
    IM("IM", "IM", 300, 0, 300, 48),
    TOF("TOF", "TOF", 1000, -50, 50, 700, 0, 700, 48),
    IMvsTheta("IMvsTheta", "IMvsTheta", 300, 0, 300, 180, 0, 180, 48)
{
    GHistBGSub::InitCuts(-4, 4, -530, -30);
    GHistBGSub::AddRandCut(30, 530);
}

GHExample::~GHExample()
{
}
Bool_t	GHExample::Init()
{
	return kTRUE;
}
Bool_t	GHExample::Start()
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

void	GHExample::ProcessEvent()
{
    for(int t=0; t<GetTagger()->GetNTagged(); t++)
    {
        taggerTime.Fill( GetTagger()->GetTaggedTime(t) );
        IM.Fill( GetNeutralPions()->Particle(0).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));

        for(int p=0; p<GetTracks()->GetNTracks(); p++)
        {
            if(GetTracks()->GetTheta(p)<21)
            {
                TOF.Fill( GetTagger()->GetTaggedTime(t) - GetTracks()->GetTime(p),
                          GetTracks()->GetClusterEnergy(p), GetTagger()->GetTaggedTime(t) );
            }
        }

        IMvsTheta.Fill( GetNeutralPions()->Particle(0).M(), GetNeutralPions()->Particle(0).Theta()*TMath::RadToDeg(),
                        GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t) );
    }
}

void	GHExample::ProcessScalerRead()
{
    IM.ScalerReadCorrection(GetScalers()->GetScaler(1)/GetScalers()->GetScaler(0));
}
