#include "PPi0Example.h"

PPi0Example::PPi0Example()
{
    test1 	= new GH1("test1", 	"test1", 	1400, -700, 700);
    test2 	= new GH3("test2", 	"test2", 	1400, -700, 700, 48, 0, 48, 180, 0, 180);
    test3 	= new GH3("test3", 	"test3", 	1400, -700, 700, 48, 0, 48, 180, 0, 180);

    GHistBGSub::InitCuts(-20, 20, -500, -30);
    GHistBGSub::AddRandCut(30, 500);
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
    if(eta->GetNParticles()>0)
    {
        test1->Fill(eta->Particle(0).M());
        for(int t=0; t<tagger->GetNTagged(); t++)
        {
            test2->Fill(eta->Particle(0).M(), tagger->GetTagged_ch(t), eta->Particle(0).Theta(), tagger->GetTagged_t(t));
            test3->Fill(eta->Particle(0).M(), tagger->GetTagged_ch(t), eta->Particle(0).Theta(), *tagger);
        }
    }
}

void	PPi0Example::ProcessScalerRead()
{
    test1->ScalerReadCorrection(5);
    test2->ScalerReadCorrection(5);
    test3->ScalerReadCorrection(5);
}

