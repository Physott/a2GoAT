#include "PPi0Example.h"

PPi0Example::PPi0Example()
{
    GHistBGSub::InitCuts(-20, 15, -100, -40);
    GHistBGSub::AddRandCut(35, 95);

    SetTarget(938);

    time 	= new GH1("time", 	"time", 	1400, -700,	700);
    time_2g = new GH1("time_2g","time_2g", 	1400, -700, 700);

    IM 		= new GH1("IM", 	"IM", 		400,   0,	400);
    IM_2g 	= new GH1("IM_2g", 	"IM_2g", 	400,   0, 	400);

    MM		= new GH1("MM", 	"MM", 	 	400,   800, 1200);
    MM_2g	= new GH1("MM_2g", 	"MM_2g", 	400,   800, 1200);
}

PPi0Example::~PPi0Example()
{
}

Bool_t	PPi0Example::Start()
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

void	PPi0Example::ProcessEvent()
{
    // Some neutral decays
    for (Int_t i = 0; i < eta->GetNParticles(); i++)
    {
        // Fill MM for 2 photon decay
        for (Int_t t = 0; t < tagger->GetNTagged(); t++)
        {
            time->Fill(tagger->GetTagged_t(t));
            MM->Fill((tagger->GetVectorProtonTarget(t)-eta->Particle(i)).M(), tagger->GetTagged_t(t), tagger->GetTagged_ch(t), eta->Particle(i).Theta()*TMath::RadToDeg());
        }
        IM->Fill(eta->Particle(i).M(), *tagger, eta->Particle(i).Theta(), kTRUE, kTRUE);
    }

}

void	PPi0Example::ProcessScalerRead()
{

    //time.ScalerReadCorrection(5);
}
