#include "MyPhysics.h"

MyPhysics::MyPhysics()    :
    proton_eta(TString("eta")),
    hist_eta(TString("eta")),
    hist_eta_proton(TString("eta_proton")),
    proton_etap(TString("etap")),
    hist_etap(TString("etap"), kTRUE),
    hist_etap_proton(TString("etap_proton"), kTRUE)
{
    PHist::SetCuts(-10, 5, -515, -15, 15, 510);
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

    proton_eta.Clear();
    hist_eta.Clear();
    hist_eta_proton.Clear();
    proton_etap.Clear();
    hist_etap.Clear();
    hist_etap_proton.Clear();

    TraverseEntries(0, eta->GetNEntries());

    Write();
	return kTRUE;
}

void	MyPhysics::ProcessEvent()
{
    if((GetEventNumber() % 100000 == 0) && GetEventNumber()!=0) cout << "Event: "<< GetEventNumber() << " Total Etas found: " << hist_eta.GetNFound() << endl;

    if(eta->GetNParticles()>0)
    {
        if(protons->GetNParticles()>0)
        {
            if(proton_eta.ProcessEvent(eta->Particle(0), *protons, *tagger))
            {
                hist_eta_proton.ProcessEvent(*eta, *tagger);
                return;
            }
        }
        hist_eta.ProcessEvent(*eta, *tagger);
        return;
    }

    if(etap->GetNParticles()>0)
    {
        if(protons->GetNParticles()>0)
        {
            if(proton_etap.ProcessEvent(etap->Particle(0), *protons, *tagger))
            {
                hist_etap_proton.ProcessEvent(*etap, *tagger);
                return;
            }
        }
        hist_etap.ProcessEvent(*etap, *tagger);
        return;
    }
}


Bool_t 	MyPhysics::Write()
{
    file_out->cd();
    TDirectory* curDir  = gDirectory->GetDirectory("eta");
    if(!curDir)
    {
        file_out->cd();
        gDirectory->mkdir("eta");
        curDir  = file_out->GetDirectory("eta");
    }

    curDir->cd();
    curDir  = gDirectory->GetDirectory("NoProton");
    if(!curDir)
    {
        curDir  = file_out->GetDirectory("eta");
        curDir->cd();
        gDirectory->mkdir("NoProton");
        curDir  = curDir->GetDirectory("NoProton");
    }
    hist_eta.Write(*curDir);

    curDir->cd();
    curDir  = gDirectory->GetDirectory("WithProton");
    if(!curDir)
    {
        curDir  = file_out->GetDirectory("eta");
        curDir->cd();
        gDirectory->mkdir("WithProton");
        curDir  = curDir->GetDirectory("WithProton");
    }
    proton_eta.Write(*curDir);
    hist_eta_proton.Write(*curDir);



    file_out->cd();
    curDir  = gDirectory->GetDirectory("etap");
    if(!curDir)
    {
        file_out->cd();
        gDirectory->mkdir("etap");
        curDir  = file_out->GetDirectory("etap");
    }

    curDir->cd();
    curDir  = gDirectory->GetDirectory("NoProton");
    if(!curDir)
    {
        curDir  = file_out->GetDirectory("etap");
        curDir->cd();
        gDirectory->mkdir("NoProton");
        curDir  = curDir->GetDirectory("NoProton");
    }
    hist_etap.Write(*curDir);

    curDir->cd();
    curDir  = gDirectory->GetDirectory("WithProton");
    if(!curDir)
    {
        curDir  = file_out->GetDirectory("etap");
        curDir->cd();
        gDirectory->mkdir("WithProton");
        curDir  = curDir->GetDirectory("WithProton");
    }
    proton_etap.Write(*curDir);
    hist_etap_proton.Write(*curDir);



	return kTRUE;
}
