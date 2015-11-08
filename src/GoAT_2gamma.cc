#include "GoAT_2gamma.h"


GoAT_2gamma::GoAT_2gamma() :
    nEvents_written(0),
    im("im", "im", 1650, 0, 1650, 48),
    imCB("imCB", "imCB", 1650, 0, 1650, 48),
    imTAPS("imTAPS", "imTAPS", 1650, 0, 1650, 48),
    im2gamma("im2gamma", "im2gamma", 1650, 0, 1650, 48),
    im2gammaCB("im2gammaCB", "im2gammaCB", 1650, 0, 1650, 48),
    im2gammaTAPS("im2gammaTAPS", "im2gammaTAPS", 1650, 0, 1650, 48),
    im2gammaProton("im2gammaProton", "im2gammaProton", 1650, 0, 1650, 48),
    im2gammaProtonCB("im2gammaProtonCB", "im2gammaProtonCB", 1650, 0, 1650, 48),
    im2gammaProtonTAPS("im2gammaProtonTAPS", "im2gammaProtonTAPS", 1650, 0, 1650, 48)
{
    GHistBGSub::InitCuts(-20, 20, -535, -35);
    GHistBGSub::AddRandCut(35, 535);
}

GoAT_2gamma::~GoAT_2gamma()
{
}

Bool_t	GoAT_2gamma::Init(const char* configfile)
{
    return kTRUE;
}

void	GoAT_2gamma::ProcessEvent()
{
    if(GetTracks()->GetNTracks()==2)
    {
        for(int t=0; t<GetTagger()->GetNTagged(); t++)
        {
            im.Fill((GetTracks()->GetVector(0)+GetTracks()->GetVector(1)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            im2gamma.Fill((GetTracks()->GetVector(0)+GetTracks()->GetVector(1)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            if(GetTracks()->GetDetectors(0) & GetTracks()->DETECTOR_NaI && GetTracks()->GetDetectors(1) & GetTracks()->DETECTOR_NaI)
            {
                im.Fill((GetTracks()->GetVector(0)+GetTracks()->GetVector(1)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
                im2gammaCB.Fill((GetTracks()->GetVector(0)+GetTracks()->GetVector(1)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            }
            if(GetTracks()->GetDetectors(0) & GetTracks()->DETECTOR_BaF2 && GetTracks()->GetDetectors(1) & GetTracks()->DETECTOR_BaF2)
            {
                imTAPS.Fill((GetTracks()->GetVector(0)+GetTracks()->GetVector(1)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
                im2gammaTAPS.Fill((GetTracks()->GetVector(0)+GetTracks()->GetVector(1)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            }
        }
    }
    else if(GetTracks()->GetNTracks()==3)
    {
        int i   = 0;
        int j   = 1;
        for(int t=0; t<GetTagger()->GetNTagged(); t++)
        {
            im.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            im2gammaProton.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            if(GetTracks()->GetDetectors(i) & GetTracks()->DETECTOR_NaI && GetTracks()->GetDetectors(j) & GetTracks()->DETECTOR_NaI)
            {
                imCB.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
                im2gammaProtonCB.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            }
            if(GetTracks()->GetDetectors(i) & GetTracks()->DETECTOR_BaF2 && GetTracks()->GetDetectors(j) & GetTracks()->DETECTOR_BaF2)
            {
                imTAPS.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
                im2gammaProtonTAPS.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            }
        }
        j   = 2;
        for(int t=0; t<GetTagger()->GetNTagged(); t++)
        {
            im.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            im2gammaProton.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            if(GetTracks()->GetDetectors(i) & GetTracks()->DETECTOR_NaI && GetTracks()->GetDetectors(j) & GetTracks()->DETECTOR_NaI)
            {
                imCB.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
                im2gammaProtonCB.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            }
            if(GetTracks()->GetDetectors(i) & GetTracks()->DETECTOR_BaF2 && GetTracks()->GetDetectors(j) & GetTracks()->DETECTOR_BaF2)
            {
                imTAPS.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
                im2gammaProtonTAPS.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            }
        }
        i   = 1;
        for(int t=0; t<GetTagger()->GetNTagged(); t++)
        {
            im.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            im2gammaProton.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            if(GetTracks()->GetDetectors(i) & GetTracks()->DETECTOR_NaI && GetTracks()->GetDetectors(j) & GetTracks()->DETECTOR_NaI)
            {
                imCB.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
                im2gammaProtonCB.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            }
            if(GetTracks()->GetDetectors(i) & GetTracks()->DETECTOR_BaF2 && GetTracks()->GetDetectors(j) & GetTracks()->DETECTOR_BaF2)
            {
                imTAPS.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
                im2gammaProtonTAPS.Fill((GetTracks()->GetVector(i)+GetTracks()->GetVector(j)).M(), GetTagger()->GetTaggedTime(t), GetTagger()->GetTaggedChannel(t));
            }
        }
    }
    else
        return;

    FillReadList();
}

Bool_t	GoAT_2gamma::Start()
{
    if(!IsAcquFile())
    {
        cout << "ERROR: Input File is not a Acqu file." << endl;
        return kFALSE;
    }
    SetAsGoATFile();



    if(!TraverseValidEvents())		return kFALSE;

    return kTRUE;
}
