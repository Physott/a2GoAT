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
    hist_fit1(TString(name).Append("_SubImCut_fit1"), TString(title).Append(" SubImCut fit1"), 24, 10, kFALSE),
    hist_fit3(TString(name).Append("_SubImCut_fit3"), TString(title).Append(" SubImCut fit3"), 24, 10, kFALSE),
    hist_fit4(TString(name).Append("_SubImCut_fit4"), TString(title).Append(" SubImCut fit4"), 24, 10, kFALSE),
    fit1(kTRUE),
    fit3(kTRUE),
    fit4(kTRUE)
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
    hist_fit1.CalcResult();
    hist_fit3.CalcResult();
    hist_fit4.CalcResult();
}

void    GAnalysis3Mesons::Fill(const GTreeMeson& meson, const GTreeParticle& photons, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    mm;
    Double_t    sub_im_0    = (photons.Particle(0) + photons.Particle(1)).M();
    Double_t    sub_im_1    = (photons.Particle(2) + photons.Particle(3)).M();
    Double_t    sub_im_2    = (photons.Particle(4) + photons.Particle(5)).M();

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
        if(CreateHistogramsForTaggerBinning==kTRUE)
            hist_raw.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
        else
            hist_raw.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i));
    }

    if((sub_im_0>500 && sub_im_0<580) &&
        (sub_im_1>100 && sub_im_1<170) &&
        (sub_im_2>100 && sub_im_2<170))
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
            if(CreateHistogramsForTaggerBinning==kTRUE)
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            else
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i));

            if(mm>850 && mm<1025)
            {
                fit1.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5), tagger.GetVectorProtonTarget(i));
                fit3.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
                fit4.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5), tagger.GetVectorProtonTarget(i));

                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                else
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i));

                while(fit1.Solve()>0)
                    hist_fit1.Fill(fit1);
                if(fit1.ConfidenceLevel()>0.05)
                {
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        hist_fit1.FillFinal(fit1, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                    else
                        hist_fit1.FillFinal(fit1, tagger.GetTaggedTime(i));
                }

                while(fit3.Solve()>0)
                    hist_fit3.Fill(fit3);
                if(fit3.ConfidenceLevel()>0.05)
                {
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        hist_fit3.FillFinal(fit3, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                    else
                        hist_fit3.FillFinal(fit3, tagger.GetTaggedTime(i));
                }

                while(fit4.Solve()>0)
                    hist_fit4.Fill(fit4);
                if(fit4.ConfidenceLevel()>0.05)
                {
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        hist_fit4.FillFinal(fit4, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                    else
                        hist_fit4.FillFinal(fit4, tagger.GetTaggedTime(i));
                }
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

    folder  = h->GetDirectory("MM_Cut");
    hist_MMCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
        GHistWriteList* subfolder  = folder->GetDirectory("fit1");
        hist_fit1.PrepareWriteList(subfolder, TString(name).Append("_fit1").Data());
        subfolder  = folder->GetDirectory("fit3");
        hist_fit3.PrepareWriteList(subfolder, TString(name).Append("_fit3").Data());
        subfolder  = folder->GetDirectory("fit4");
        hist_fit4.PrepareWriteList(subfolder, TString(name).Append("_fit4").Data());
}

void    GAnalysis3Mesons::Reset(Option_t* option)
{
    hist_raw.Reset(option);
    hist_SubImCut.Reset(option);
    hist_MMCut.Reset(option);
    hist_fit1.Reset(option);
    hist_fit3.Reset(option);
    hist_fit4.Reset(option);
}

void    GAnalysis3Mesons::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    hist_raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MMCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit1.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
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
    hist_raw_TOF(TString(name).Append("raw_TOF"), TString(title).Append("raw_TOF"), 300, -15, 15, 800, 0, 800, kFALSE),
    hist_SubImCut(TString(name).Append("SubImCut"), TString(title).Append("SubImCut"), kFALSE),
    hist_SubImCut_TOF(TString(name).Append("SubImCut_TOF"), TString(title).Append("SubImCut_TOF"), 300, -15, 15, 800, 0, 800, kFALSE),
    hist_MMCut(TString(name).Append("MMCut"), TString(title).Append("MMCut"), kFALSE),
    hist_TOF(TString(name).Append("TOF"), TString(title).Append("TOF"), 300, -15, 15, 800, 0, 800, kFALSE),
    hist_fit1(TString(name).Append("fit1"), TString(title).Append("fit1"), 24, 10, kFALSE),
    hist_fit3(TString(name).Append("fit3"), TString(title).Append("fit3"), 24, 10, kFALSE),
    hist_fit4(TString(name).Append("fit4"), TString(title).Append("fit4"), 24, 10, kFALSE),
    hist_fit7Proton(TString(name).Append("fit7Proton"), TString(title).Append("fit7Proton"), 28, 10, kFALSE),
    fit1(kTRUE),
    fit3(kTRUE),
    fit4(kTRUE),
    fit7Proton(kTRUE)
{

}

GAnalysis3MesonsProton::~GAnalysis3MesonsProton()
{

}

void   GAnalysis3MesonsProton::CalcResult()
{
    checkProton.CalcResult();
    hist_raw.CalcResult();
    hist_raw_TOF.CalcResult();
    hist_SubImCut.CalcResult();
    hist_SubImCut_TOF.CalcResult();
    hist_MMCut.CalcResult();
    hist_TOF.CalcResult();
    hist_fit1.CalcResult();
    hist_fit3.CalcResult();
    hist_fit4.CalcResult();
    hist_fit7Proton.CalcResult();
}

void    GAnalysis3MesonsProton::Fill(const GTreeMeson& meson, const GTreeParticle& photons, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    theta  = meson.Particle(0).Theta();
    Double_t    phi  = meson.Particle(0).Phi();
    Double_t    mm;
    Double_t    sub_im_0    = (photons.Particle(0) + photons.Particle(1)).M();
    Double_t    sub_im_1    = (photons.Particle(2) + photons.Particle(3)).M();
    Double_t    sub_im_2    = (photons.Particle(4) + photons.Particle(5)).M();

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
        if(CreateHistogramsForTaggerBinning==kTRUE)
            hist_raw.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
        else
            hist_raw.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i));
        hist_raw_TOF.Fill(tagger.GetTaggedTime(i)-proton.GetTime(0), proton.GetClusterEnergy(0), tagger.GetTaggedTime(i));
    }

    if((sub_im_0>500 && sub_im_0<580) &&
        (sub_im_1>100 && sub_im_1<170) &&
        (sub_im_2>100 && sub_im_2<170))
    {
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            if(checkProton.Check(meson, proton, tagger.GetVectorProtonTarget(i), tagger.GetTaggedTime(i))==kFALSE)
                continue;

            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
            if(CreateHistogramsForTaggerBinning==kTRUE)
                hist_SubImCut.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            else
                hist_SubImCut.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i));
            hist_SubImCut_TOF.Fill(tagger.GetTaggedTime(i)-proton.GetTime(0), proton.GetClusterEnergy(0), tagger.GetTaggedTime(i));

            if(mm>850 && mm<1025)
            {
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MMCut.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                else
                    hist_MMCut.Fill(im, mm, theta, phi, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i));
                hist_TOF.Fill(tagger.GetTaggedTime(i)-proton.GetTime(0), proton.GetClusterEnergy(0), tagger.GetTaggedTime(i));

                fit1.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5), tagger.GetVectorProtonTarget(i));
                fit3.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
                fit4.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5), tagger.GetVectorProtonTarget(i));
                fit7Proton.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5), tagger.GetVectorProtonTarget(i), proton.Particle(0));

                while(fit1.Solve()>0)
                    hist_fit1.Fill(fit1);
                if(fit1.ConfidenceLevel()>0.05)
                {
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        hist_fit1.FillFinal(fit1, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                    else
                        hist_fit1.FillFinal(fit1, tagger.GetTaggedTime(i));
                }

                while(fit3.Solve()>0)
                    hist_fit3.Fill(fit3);
                if(fit3.ConfidenceLevel()>0.05)
                {
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        hist_fit3.FillFinal(fit3, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                    else
                        hist_fit3.FillFinal(fit3, tagger.GetTaggedTime(i));
                }

                while(fit4.Solve()>0)
                    hist_fit4.Fill(fit4);
                if(fit4.ConfidenceLevel()>0.05)
                {
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        hist_fit4.FillFinal(fit4, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                    else
                        hist_fit4.FillFinal(fit4, tagger.GetTaggedTime(i));
                }


                while(fit7Proton.Solve()>0)
                    hist_fit7Proton.Fill(fit7Proton);
                if(fit7Proton.ConfidenceLevel()>0.05)
                {
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        hist_fit7Proton.FillFinal(fit7Proton, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                    else
                        hist_fit7Proton.FillFinal(fit7Proton, tagger.GetTaggedTime(i));
                }
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
    hist_raw_TOF.PrepareWriteList(folder, TString(name).Append("_Raw_TOF").Data());

    folder  = h->GetDirectory("SubIM_Cut");
    hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
    hist_SubImCut_TOF.PrepareWriteList(folder, TString(name).Append("_subIMCut_TOF").Data());

    folder  = h->GetDirectory("MM_Cut");
    hist_TOF.PrepareWriteList(folder, TString(name).Append("_TOF").Data());
    hist_MMCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
        GHistWriteList* subfolder  = folder->GetDirectory("fit1");
        hist_fit1.PrepareWriteList(subfolder, TString(name).Append("_fit1").Data());
        subfolder  = folder->GetDirectory("fit3");
        hist_fit3.PrepareWriteList(subfolder, TString(name).Append("_fit3").Data());
        subfolder  = folder->GetDirectory("fit4");
        hist_fit4.PrepareWriteList(subfolder, TString(name).Append("_fit4").Data());
        subfolder  = folder->GetDirectory("fit7Proton");
        hist_fit7Proton.PrepareWriteList(subfolder, TString(name).Append("_fit7Proton").Data());
}

void    GAnalysis3MesonsProton::Reset(Option_t* option)
{
    checkProton.Reset(option);
    hist_raw.Reset(option);
    hist_raw_TOF.Reset(option);
    hist_SubImCut.Reset(option);
    hist_SubImCut_TOF.Reset(option);
    hist_MMCut.Reset(option);
    hist_TOF.Reset(option);
    hist_fit1.Reset(option);
    hist_fit3.Reset(option);
    hist_fit4.Reset(option);
    hist_fit7Proton.Reset(option);
}

void    GAnalysis3MesonsProton::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    //checkProton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MMCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit1.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit7Proton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

