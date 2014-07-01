#include "P3Meson.h"


P3Meson::P3Meson(const TString &_Name, const Bool_t _IsEtap)    :
    name(_Name),
    isEtap(_IsEtap),
    time_raw(TString(name).Append("_time_raw").Data(), TString(name).Append("_time_raw").Data(),		1000,-500,500),
    raw(TString(name).Append("_raw")),
    time_cutIM(TString(name).Append("_time_cutIM").Data(), TString(name).Append("_time_cutIM").Data(),		1000,-500,500),
    cutIMevent(TString(name).Append("_cutIM")),
    time_cutMM(TString(name).Append("_time_cutMM").Data(), TString(name).Append("_time_cutMM").Data(),		1000,-500,500),
    cutMMevent(TString(name).Append("_cutMM")),
    nFound(0)
{
    cutIM[0][0] = 110;
    cutIM[0][1] = 155;
    cutIM[1][0] = 110;
    cutIM[1][1] = 155;
    cutIM[2][0] = 110;
    cutIM[2][1] = 155;

    cutMM[0] = 850;
    cutMM[1] = 1050;

    if(isEtap)
    {
        cutIM[0][0] = 500;
        cutIM[0][1] = 590;
    }
}

P3Meson::~P3Meson()
{

}

/*
Bool_t	P3Meson::Start()
{
    PProtonCheck::Clear();

    TraverseEntries(0, eta->GetNEntries());

    raw.RandomSubtraction();
    cutIMevent.RandomSubtraction();
    cutMMevent.RandomSubtraction();

    PProtonCheck::RandomSubtraction();


    WriteHistograms();
	return kTRUE;
}*/

Bool_t	P3Meson::ProcessEvent(const GTreeMeson& meson, const GTreeTagger& tagger)
{
    //if(GetEventNumber() == 0) nFound = 0;
    //else if(GetEventNumber() % 100000 == 0) cout << "Event: "<< GetEventNumber() << " Total Etas found: " << nFound << endl;

    Double_t    imSub[3];
    Double_t    misMass;
    bool        passIM;

    if(meson.GetNParticles()>0)
    {

        imSub[0]    = (meson.SubParticles(0, 0)+meson.SubParticles(0, 1)).M();
        imSub[1]    = (meson.SubParticles(0, 2)+meson.SubParticles(0, 3)).M();
        imSub[2]    = (meson.SubParticles(0, 4)+meson.SubParticles(0, 5)).M();
        if((imSub[0]>cutIM[0][0] && imSub[0]<cutIM[0][1]) && (imSub[1]>cutIM[1][0] && imSub[1]<cutIM[1][1]) && (imSub[2]>cutIM[2][0] && imSub[2]<cutIM[2][1]))
            passIM  = true;
        else
            passIM  = false;

        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            time_raw.Fill(tagger.GetTagged_t(i));
            misMass = (tagger.GetVector(i)+TLorentzVector(0,0,0,MASS_PROTON) - meson.Meson(0)).M();
            raw.Fill(tagger.GetTagged_t(i), meson.Meson(0).M(), misMass);
            raw.FillSubMesons(tagger.GetTagged_t(i), imSub[0], imSub[1], imSub[2]);

            if(passIM)
            {
                time_cutIM.Fill(tagger.GetTagged_t(i));
                cutIMevent.Fill(tagger.GetTagged_t(i), meson.Meson(0).M(), (tagger.GetVector(i)+TLorentzVector(0,0,0,MASS_PROTON) - meson.Meson(0)).M());
                cutIMevent.FillSubMesons(tagger.GetTagged_t(i), imSub[0], imSub[1], imSub[2]);

                if(misMass>cutMM[0] && misMass<cutMM[1])
                {
                    time_cutMM.Fill(tagger.GetTagged_t(i));
                    cutMMevent.Fill(tagger.GetTagged_t(i), meson.Meson(0).M(), (tagger.GetVector(i)+TLorentzVector(0,0,0,MASS_PROTON) - meson.Meson(0)).M());
                    cutMMevent.FillSubMesons(tagger.GetTagged_t(i), imSub[0], imSub[1], imSub[2]);

                    nFound++;
                }
            }
        }
    }
}


Bool_t 	P3Meson::Write(TDirectory& curDir)
{
    curDir.cd();
    time_raw.Write();
    time_cutIM.Write();
    time_cutMM.Write();

    raw.Write(curDir);
    cutIMevent.Write(curDir);
    cutMMevent.Write(curDir);

	return kTRUE;
}


void    P3Meson::Clear()
{
    nFound  =0;
    time_raw.Reset();
    time_cutIM.Reset();
    time_cutMM.Reset();
    raw.Clear();
    cutIMevent.Clear();
    cutMMevent.Clear();
}

void    P3Meson::RandomSubtraction()
{
    raw.RandomSubtraction();
    cutIMevent.RandomSubtraction();
    cutMMevent.RandomSubtraction();
}

