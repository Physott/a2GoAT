#include "GAnalysis3Mesons.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"


GAnalysis3Mesons::GAnalysis3Mesons(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(_IsEtap),
    hist_raw(TString(name).Append("raw"), TString(title).Append("raw"), kFALSE),
    hist_SubImCut(TString(name).Append("SubImCut"), TString(title).Append("SubImCut"), kFALSE),
    hist_MMCut(TString(name).Append("MMCut"), TString(title).Append("MMCut"), kFALSE),
    fit3(kTRUE),
    fit4(kTRUE),
    fit4Beam(kTRUE),
    hist_fit3(TString(name).Append("_SubImCut_fit3"), TString(title).Append(" SubImCut fit3"), 24, 10, kFALSE),
    hist_fit4(TString(name).Append("_SubImCut_fit4"), TString(title).Append(" SubImCut fit4"), 24, 10, kFALSE),
    hist_fit4Beam(TString(name).Append("_SubImCut_fit4Beam"), TString(title).Append(" SubImCut fit4Beam"), 28, 10, kFALSE),
    hist_SubImCut_fit3(TString(name).Append("fit3"), TString(title).Append("fit3"), 24, 10, kFALSE),
    hist_SubImCut_fit4(TString(name).Append("fit4"), TString(title).Append("fit4"), 24, 10, kFALSE),
    hist_SubImCut_fit4Beam(TString(name).Append("fit4Beam"), TString(title).Append("fit4Beam"), 28, 10, kFALSE)
{
    if(_IsEtap==kTRUE)
        SetCutSubIM(0, 497, 697);
    else
        SetCutSubIM(0, 110, 160);
    SetCutSubIM(1, 110, 160);
    SetCutSubIM(2, 110, 160);

    SetCutMM(838, 1038);
}

GAnalysis3Mesons::~GAnalysis3Mesons()
{

}

void   GAnalysis3Mesons::CalcResult()
{
    hist_raw.CalcResult();
    hist_SubImCut.CalcResult();
    hist_MMCut.CalcResult();
    hist_SubImCut_fit3.CalcResult();
    hist_SubImCut_fit4.CalcResult();
    hist_SubImCut_fit4Beam.CalcResult();
    hist_fit3.CalcResult();
    hist_fit4.CalcResult();
    hist_fit4Beam.CalcResult();
}

void    GAnalysis3Mesons::Fill(const GTreeMeson& meson, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    mm;
    Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
    Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
    Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
        if(CreateHistogramsForTaggerBinning==kTRUE)
            hist_raw.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        else
            hist_raw.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));
    }

    if((sub_im_0>500 && sub_im_0<580) &&
        (sub_im_1>100 && sub_im_1<170) &&
        (sub_im_2>100 && sub_im_2<170))
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();if(CreateHistogramsForTaggerBinning==kTRUE)
            if(CreateHistogramsForTaggerBinning==kTRUE)
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            else
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

            if(mm>850 && mm<1025)
            {
                fit3.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5));
                fit4.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i));
                fit4Beam.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i));

                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

                while(fit3.Solve()>0)
                    hist_fit3.Fill(fit3);
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_fit3.FillFinal(fit3, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_fit3.FillFinal(fit3, tagger.GetTagged_t(i));

                /*while(fit4.Solve()>0)
                    hist_fit4.Fill(fit4);
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_fit4.FillFinal(fit4, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_fit4.FillFinal(fit4, tagger.GetTagged_t(i));

                while(fit4Beam.Solve()>0)
                    hist_fit4Beam.Fill(fit4Beam);
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_fit4Beam.FillFinal(fit4Beam, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_fit4Beam.FillFinal(fit4Beam, tagger.GetTagged_t(i));*/
            }
        }
    }
}

void    GAnalysis3Mesons::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* h  = arr->GetDirectory("WithoutProton");

    GHistWriteList* folder  = h->GetDirectory("Raw");
    hist_raw.PrepareWriteList(folder, TString(name).Append("_Raw").Data());

    folder  = h->GetDirectory("SubIM_Cut");
    hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
        GHistWriteList* subfolder  = folder->GetDirectory("fit3");
        hist_SubImCut_fit3.PrepareWriteList(subfolder, TString(name).Append("_fit3").Data());
        subfolder  = folder->GetDirectory("fit4");
        hist_SubImCut_fit4.PrepareWriteList(subfolder, TString(name).Append("_fit4").Data());
        subfolder  = folder->GetDirectory("fit4Beam");
        hist_SubImCut_fit4Beam.PrepareWriteList(subfolder, TString(name).Append("_fit4Beam").Data());

    folder  = h->GetDirectory("MM_Cut");
    hist_MMCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
        subfolder  = folder->GetDirectory("fit3");
        hist_fit3.PrepareWriteList(subfolder, TString(name).Append("_fit3").Data());
        subfolder  = folder->GetDirectory("fit4");
        hist_fit4.PrepareWriteList(subfolder, TString(name).Append("_fit4").Data());
        subfolder  = folder->GetDirectory("fit4Beam");
        hist_fit4Beam.PrepareWriteList(subfolder, TString(name).Append("_fit4Beam").Data());
}

void    GAnalysis3Mesons::Reset(Option_t* option)
{
    hist_raw.Reset(option);
    hist_SubImCut.Reset(option);
    hist_MMCut.Reset(option);
    hist_SubImCut_fit3.Reset(option);
    hist_SubImCut_fit4.Reset(option);
    hist_SubImCut_fit4Beam.Reset(option);
    hist_fit3.Reset(option);
    hist_fit4.Reset(option);
    hist_fit4Beam.Reset(option);
}

void    GAnalysis3Mesons::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    hist_raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MMCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit4Beam.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit4Beam.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

void    GAnalysis3Mesons::SetCutSubIM(const Int_t subNumber, const Double_t min, const Double_t max)
{
    cutSubIM[2*subNumber] = min;
    cutSubIM[(2*subNumber)+1] = max;
}

void    GAnalysis3Mesons::SetCutMM(const Double_t min, const Double_t max)
{
    cutMM[0] = min;
    cutMM[1] = max;
}











GAnalysis3MesonsProton::GAnalysis3MesonsProton(const char* name, const char* title, const Bool_t IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(IsEtap),
    checkProton(TString(name).Append("checkProton"), TString(title).Append("checkProton"), kFALSE),
    hist_raw(TString(name).Append("raw"), TString(title).Append("raw"), kFALSE),
    hist_SubImCut(TString(name).Append("SubImCut"), TString(title).Append("SubImCut"), kFALSE),
    hist_MMCut(TString(name).Append("MMCut"), TString(title).Append("MMCut"), kFALSE),
    fit3(kTRUE),
    fit4(kTRUE),
    fit4Beam(kTRUE),
    fit4Proton(kTRUE),
    fit4BeamProton(kTRUE),
    hist_SubImCut_fit3(TString(name).Append("_SubImCut_fit3"), TString(title).Append(" SubImCut fit3"), 24, 10, kFALSE),
    hist_SubImCut_fit4(TString(name).Append("_SubImCut_fit4"), TString(title).Append(" SubImCut fit4"), 24, 10, kFALSE),
    hist_SubImCut_fit4Beam(TString(name).Append("_SubImCut_fit4Beam"), TString(title).Append(" SubImCut fit4Beam"), 28, 10, kFALSE),
    hist_SubImCut_fit4Proton(TString(name).Append("_SubImCut_fit4Proton"), TString(title).Append(" SubImCut fit4Proton"), 28, 10, kFALSE),
    hist_SubImCut_fit4BeamProton(TString(name).Append("_SubImCut_fit4BeamProton"), TString(title).Append(" SubImCut fit4BeamProton"), 32, 10, kFALSE),
    hist_fit3(TString(name).Append("fit3"), TString(title).Append("fit3"), 24, 10, kFALSE),
    hist_fit4(TString(name).Append("fit4"), TString(title).Append("fit4"), 24, 10, kFALSE),
    hist_fit4Beam(TString(name).Append("fit4Beam"), TString(title).Append("fit4Beam"), 28, 10, kFALSE),
    hist_fit4Proton(TString(name).Append("fit4Proton"), TString(title).Append("fit4Proton"), 28, 10, kFALSE),
    hist_fit4BeamProton(TString(name).Append("fit4BeamProton"), TString(title).Append("fit4BeamProton"), 32, 10, kFALSE)
{

}

GAnalysis3MesonsProton::~GAnalysis3MesonsProton()
{

}

TLorentzVector  GAnalysis3MesonsProton::GetCorrectedProton(const GTreeParticle& proton)
{
    TLorentzVector  ret(proton.Particle(0));

    Double_t    E   = ret.E()-0.93827+938.27;
    Double_t    P   = sqrt((E*E)-(938.27*938.27));

    ret.SetPtEtaPhiE(P, proton.Particle(0).Theta(), proton.Particle(0).Phi(), E);
    return ret;
}

void   GAnalysis3MesonsProton::CalcResult()
{
    checkProton.CalcResult();
    hist_raw.CalcResult();
    hist_SubImCut.CalcResult();
    hist_MMCut.CalcResult();
    hist_SubImCut_fit3.CalcResult();
    hist_SubImCut_fit4.CalcResult();
    hist_SubImCut_fit4Beam.CalcResult();
    hist_SubImCut_fit4Proton.CalcResult();
    hist_SubImCut_fit4BeamProton.CalcResult();
    hist_fit3.CalcResult();
    hist_fit4.CalcResult();
    hist_fit4Beam.CalcResult();
    hist_fit4Proton.CalcResult();
    hist_fit4BeamProton.CalcResult();
}

void    GAnalysis3MesonsProton::Fill(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    theta  = meson.Particle(0).Theta();
    Double_t    phi  = meson.Particle(0).Phi();
    Double_t    mm;
    Double_t    sub_im_0    = (meson.SubPhotons(0, 0) + meson.SubPhotons(0, 1)).M();
    Double_t    sub_im_1    = (meson.SubPhotons(0, 2) + meson.SubPhotons(0, 3)).M();
    Double_t    sub_im_2    = (meson.SubPhotons(0, 4) + meson.SubPhotons(0, 5)).M();

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
        if(CreateHistogramsForTaggerBinning==kTRUE)
            hist_raw.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
        else
            hist_raw.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));
    }

    if((sub_im_0>500 && sub_im_0<580) &&
        (sub_im_1>100 && sub_im_1<170) &&
        (sub_im_2>100 && sub_im_2<170))
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            if(checkProton.Check(meson, proton, tagger.GetVectorProtonTarget(i), tagger.GetTagged_t(i))==kFALSE)
                continue;

            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();if(CreateHistogramsForTaggerBinning==kTRUE)
            if(CreateHistogramsForTaggerBinning==kTRUE)
                hist_SubImCut.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            else
                hist_SubImCut.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

            if(mm>850 && mm<1025)
            {
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MMCut.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_MMCut.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

                fit3.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5));
                fit4.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i));
                fit4Beam.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i));
                fit4Proton.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i), GetCorrectedProton(proton));
                fit4BeamProton.Set(meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), tagger.GetVectorProtonTarget(i), GetCorrectedProton(proton));

                while(fit3.Solve()>0)
                    hist_fit3.Fill(fit3);
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_fit3.FillFinal(fit3, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_fit3.FillFinal(fit3, tagger.GetTagged_t(i));

                /*while(fit4.Solve()>0)
                    hist_fit4.Fill(fit4);
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_fit4.FillFinal(fit4, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_fit4.FillFinal(fit4, tagger.GetTagged_t(i));

                while(fit4Beam.Solve()>0)
                    hist_fit4Beam.Fill(fit4Beam);
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_fit4Beam.FillFinal(fit4Beam, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_fit4Beam.FillFinal(fit4Beam, tagger.GetTagged_t(i));

                while(fit4Proton.Solve()>0)
                    hist_fit4Proton.Fill(fit4Proton);
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_fit4Proton.FillFinal(fit4Proton, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_fit4Proton.FillFinal(fit4Proton, tagger.GetTagged_t(i));

                while(fit4BeamProton.Solve()>0)
                    hist_fit4BeamProton.Fill(fit4BeamProton);
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_fit4BeamProton.FillFinal(fit4BeamProton, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_fit4BeamProton.FillFinal(fit4BeamProton, tagger.GetTagged_t(i));*/
            }
        }
    }
}

void    GAnalysis3MesonsProton::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* h  = arr->GetDirectory("WithProton");

    GHistWriteList* folder  = h->GetDirectory("CheckProton");
    checkProton.PrepareWriteList(folder, TString(name).Append("_CheckProton").Data());

    folder  = h->GetDirectory("Raw");
    hist_raw.PrepareWriteList(folder, TString(name).Append("_Raw").Data());

    folder  = h->GetDirectory("SubIM_Cut");
    hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
        GHistWriteList* subfolder  = folder->GetDirectory("fit3");
        hist_SubImCut_fit3.PrepareWriteList(subfolder, TString(name).Append("_fit3").Data());
        subfolder  = folder->GetDirectory("fit4");
        hist_SubImCut_fit4.PrepareWriteList(subfolder, TString(name).Append("_fit4").Data());
        subfolder  = folder->GetDirectory("fit4Beam");
        hist_SubImCut_fit4Beam.PrepareWriteList(subfolder, TString(name).Append("_fit4Beam").Data());
        subfolder  = folder->GetDirectory("fit4Proton");
        hist_SubImCut_fit4Proton.PrepareWriteList(subfolder, TString(name).Append("_fit4Proton").Data());
        subfolder  = folder->GetDirectory("fit4BeamProton");
        hist_SubImCut_fit4BeamProton.PrepareWriteList(subfolder, TString(name).Append("_fit4BeamProton").Data());

    folder  = h->GetDirectory("MM_Cut");
    hist_MMCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
        subfolder  = folder->GetDirectory("fit3");
        hist_fit3.PrepareWriteList(subfolder, TString(name).Append("_fit3").Data());
        subfolder  = folder->GetDirectory("fit4");
        hist_fit4.PrepareWriteList(subfolder, TString(name).Append("_fit4").Data());
        subfolder  = folder->GetDirectory("fit4Beam");
        hist_fit4Beam.PrepareWriteList(subfolder, TString(name).Append("_fit4Beam").Data());
        subfolder  = folder->GetDirectory("fit4Proton");
        hist_fit4Proton.PrepareWriteList(subfolder, TString(name).Append("_fit4Proton").Data());
        subfolder  = folder->GetDirectory("fit4BeamProton");
        hist_fit4BeamProton.PrepareWriteList(subfolder, TString(name).Append("_fit4BeamProton").Data());
}

void    GAnalysis3MesonsProton::Reset(Option_t* option)
{
    checkProton.Reset(option);
    hist_raw.Reset(option);
    hist_SubImCut.Reset(option);
    hist_MMCut.Reset(option);
    hist_SubImCut_fit3.Reset(option);
    hist_SubImCut_fit4.Reset(option);
    hist_SubImCut_fit4Beam.Reset(option);
    hist_SubImCut_fit4Proton.Reset(option);
    hist_SubImCut_fit4BeamProton.Reset(option);
    hist_fit3.Reset(option);
    hist_fit4.Reset(option);
    hist_fit4Beam.Reset(option);
    hist_fit4Proton.Reset(option);
    hist_fit4BeamProton.Reset(option);
}

void    GAnalysis3MesonsProton::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    //checkProton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MMCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit4Beam.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit4Proton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit4BeamProton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit4Beam.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit4Proton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit4BeamProton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

