#include "PHistEvent.h"
#include "GTreeTagger.h"


PHistEvent::PHistEvent(const TString& _Name)    :
    name(_Name),
    time(TString(name).Append("_time").Data(), TString(name).Append(" tagger time").Data(),		1000,-500,500),
    IM(TString("_IM").Prepend(name), TString("IM ").Append(name), 1500, 0, 1500),
    MM(TString("_MM").Prepend(name), TString("MM ").Append(name), 1500, 0, 1500)
{

}

PHistEvent::~PHistEvent()
{

}

void    PHistEvent::Fill(const Double_t invMass, const TLorentzVector &particle, const GTreeTagger &tagger)
{
    Double_t    misMass;
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        misMass = (tagger.GetVectorProtonTarget(i) - particle).M();
        time.Fill(tagger.GetTagged_t(i));
        IM.Fill(invMass, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        MM.Fill(misMass, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
    }
}

void    PHistEvent::Write(TDirectory& dir)
{
    TDirectory* curDir  = dir.GetDirectory(name.Data());
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir(name.Data());
        curDir  = dir.GetDirectory(name.Data());
    }
    curDir->cd();
    time.Write();
    IM.Write(*curDir);
    MM.Write(*curDir);
}






PHistEvent3Meson::PHistEvent3Meson(const TString& _Name)    :
    PHistEvent(_Name),
    sub0(TString("_sub0").Prepend(name), TString("IM sub0 ").Append(name), 1500, 0, 1500),
    sub1(TString("_sub1").Prepend(name), TString("IM sub1 ").Append(name), 1500, 0, 1500),
    sub2(TString("_sub2").Prepend(name), TString("IM sub2 ").Append(name), 1500, 0, 1500)
{

}

PHistEvent3Meson::~PHistEvent3Meson()
{

}

void    PHistEvent3Meson::FillSubMesons(const Double_t invMassSub0, const Double_t invMassSub1, const Double_t invMassSub2, const GTreeTagger& tagger)
{
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        sub0.Fill(invMassSub0, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        sub1.Fill(invMassSub1, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        sub2.Fill(invMassSub2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
    }
}

void    PHistEvent3Meson::Write(TDirectory& dir)
{
    TDirectory* curDir  = dir.GetDirectory(name.Data());
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir(name.Data());
        curDir  = dir.GetDirectory(name.Data());
    }
    curDir->cd();
    time.Write();
    IM.Write(*curDir);
    MM.Write(*curDir);
    sub0.Write(*curDir);
    sub1.Write(*curDir);
    sub2.Write(*curDir);
}






PHistEvent3MesonFit::PHistEvent3MesonFit(const TString& _Name)    :
    PHistEvent3Meson(_Name),
    ChiSq(TString("_ChiSq").Prepend(name), TString("Chi Sq ").Append(name), 100, 0, 10),
    ConfidenceLevel(TString("_ConfidenceLevel").Prepend(name), TString("Confidence Level ").Append(name), 1000, 0, 1),
    Pull({{PHistD(TString("_Pull_0_X").Prepend(name), TString(" Pull Photon 0 Px ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_0_Y").Prepend(name), TString(" Pull Photon 0 Py ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_0_Z").Prepend(name), TString(" Pull Photon 0 Pz ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_0_E").Prepend(name), TString(" Pull Photon 0 E ").Append(name), 200, -10, 10)},
        {PHistD(TString("_Pull_1_X").Prepend(name), TString(" Pull Photon 0 Px ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_1_Y").Prepend(name), TString(" Pull Photon 0 Py ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_1_Z").Prepend(name), TString(" Pull Photon 0 Pz ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_1_E").Prepend(name), TString(" Pull Photon 0 E ").Append(name), 200, -10, 10)},
        {PHistD(TString("_Pull_2_X").Prepend(name), TString(" Pull Photon 0 Px ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_2_Y").Prepend(name), TString(" Pull Photon 0 Py ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_2_Z").Prepend(name), TString(" Pull Photon 0 Pz ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_2_E").Prepend(name), TString(" Pull Photon 0 E ").Append(name), 200, -10, 10)},
        {PHistD(TString("_Pull_3_X").Prepend(name), TString(" Pull Photon 0 Px ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_3_Y").Prepend(name), TString(" Pull Photon 0 Py ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_3_Z").Prepend(name), TString(" Pull Photon 0 Pz ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_3_E").Prepend(name), TString(" Pull Photon 0 E ").Append(name), 200, -10, 10)},
        {PHistD(TString("_Pull_4_X").Prepend(name), TString(" Pull Photon 0 Px ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_4_Y").Prepend(name), TString(" Pull Photon 0 Py ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_4_Z").Prepend(name), TString(" Pull Photon 0 Pz ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_4_E").Prepend(name), TString(" Pull Photon 0 E ").Append(name), 200, -10, 10)},
        {PHistD(TString("_Pull_5_X").Prepend(name), TString(" Pull Photon 0 Px ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_5_Y").Prepend(name), TString(" Pull Photon 0 Py ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_5_Z").Prepend(name), TString(" Pull Photon 0 Pz ").Append(name), 200, -10, 10),
         PHistD(TString("_Pull_5_E").Prepend(name), TString(" Pull Photon 0 E ").Append(name), 200, -10, 10)}})
{
}

PHistEvent3MesonFit::~PHistEvent3MesonFit()
{

}

void    PHistEvent3MesonFit::FillFit(const Double_t _ChiSq, const Double_t _ConfidenceLevel, const Double_t* _Pull, const GTreeTagger& tagger)
{
    Double_t    misMass;
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        ChiSq.Fill(_ChiSq, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        ConfidenceLevel.Fill(_ConfidenceLevel, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        for(int i=0; i<6; i++)
        {
            for(int k=0; k<4; k++)
                Pull[i][k].Fill(_Pull[(4*i)+k], tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        }
    }
}

void    PHistEvent3MesonFit::Write(TDirectory& dir)
{
    TDirectory* curDir  = dir.GetDirectory(name.Data());
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir(name.Data());
        curDir  = dir.GetDirectory(name.Data());
    }
    curDir->cd();
    time.Write();
    IM.Write(*curDir);
    MM.Write(*curDir);
    sub0.Write(*curDir);
    sub1.Write(*curDir);
    sub2.Write(*curDir);
    ChiSq.Write(*curDir);
    ConfidenceLevel.Write(*curDir);
    for(int i=0; i<6; i++)
    {
        for(int k=0; k<4; k++)
            Pull[i][k].Write(*curDir);
    }
}
