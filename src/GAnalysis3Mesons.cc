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
    fit(),
    hist_fit(TString(name).Append("_SubImCut_fit"), TString(title).Append(" SubImCut fit"), 40, kFALSE),
    hist_SubImCut_fit(TString(name).Append("fit"), TString(title).Append("fit"), 40, kFALSE)
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
    hist_SubImCut_fit.CalcResult();
    hist_fit.CalcResult();
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
        GKinFitterPolarToCartesian  converter;
        for(int i=0; i<6; i++)
            converter.Set(i, meson.SubPhotons(0, i).E(), meson.SubPhotons(0, i).Theta(), meson.SubPhotons(0, i).Phi(), 10, TMath::DegToRad()*3, TMath::DegToRad()*3);

        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();if(CreateHistogramsForTaggerBinning==kTRUE)
            if(CreateHistogramsForTaggerBinning==kTRUE)
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            else
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

            //std::cout << "here  "; tagger.GetVectorProtonTarget(i).Print();

            converter.Set(tagger.GetPhotonBeam_E(i), 2);
            fit.Set(converter);

            if(fit.Solve()>0)
            {
                if(fit.ConfidenceLevel()>0.1)
                {
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        hist_SubImCut_fit.Fill(fit, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        hist_SubImCut_fit.Fill(fit, tagger.GetTagged_t(i));
                }
            }

            if(mm>850 && mm<1025)
            {
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

                /*if(fit.IsSolved()==kTRUE)
                {
                    if(fit.ConfidenceLevel()>0.1)
                    {
                        if(CreateHistogramsForTaggerBinning==kTRUE)
                            hist_fit.Fill(fit, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                        else
                            hist_fit.Fill(fit, tagger.GetTagged_t(i));
                    }
                }*/
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
        GHistWriteList* subfolder  = folder->GetDirectory("fit");
        hist_SubImCut_fit.PrepareWriteList(subfolder, TString(name).Append("_fit").Data());

    folder  = h->GetDirectory("MM_Cut");
    hist_MMCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
        subfolder  = folder->GetDirectory("fit");
        hist_fit.PrepareWriteList(subfolder, TString(name).Append("_fit").Data());
}

void    GAnalysis3Mesons::Reset(Option_t* option)
{
    hist_raw.Reset(option);
    hist_SubImCut.Reset(option);
    hist_MMCut.Reset(option);
    hist_SubImCut_fit.Reset(option);
    hist_fit.Reset(option);
}

void    GAnalysis3Mesons::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    hist_raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MMCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
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
    fit(),
    hist_SubImCut_fit(TString(name).Append("_SubImCut_fit"), TString(title).Append(" SubImCut fit"), 24, kFALSE),
    hist_fit(TString(name).Append("fit"), TString(title).Append("fit"), 24, kFALSE)
{

}

GAnalysis3MesonsProton::~GAnalysis3MesonsProton()
{

}

void   GAnalysis3MesonsProton::CalcResult()
{
    checkProton.CalcResult();
    hist_raw.CalcResult();
    hist_SubImCut.CalcResult();
    hist_MMCut.CalcResult();
    hist_SubImCut_fit.CalcResult();
    hist_fit.CalcResult();
}

void    GAnalysis3MesonsProton::Fill(const GTreeMeson& meson, const GTreeParticle& proton, const GTreeTagger& tagger, const Bool_t CreateHistogramsForTaggerBinning)
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
            if(checkProton.Check(meson, proton, tagger.GetVectorProtonTarget(i), tagger.GetTagged_t(i))==kFALSE)
                continue;

            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();if(CreateHistogramsForTaggerBinning==kTRUE)
            if(CreateHistogramsForTaggerBinning==kTRUE)
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            else
                hist_SubImCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

            fit.Set(tagger.GetVectorProtonTarget(i).E(), meson.SubPhotons(0, 0), meson.SubPhotons(0, 1), meson.SubPhotons(0, 2), meson.SubPhotons(0, 3), meson.SubPhotons(0, 4), meson.SubPhotons(0, 5), proton.Particle(0));

            if(fit.Solve()>0)
            {
                if(fit.ConfidenceLevel()>0.1)
                {
                    if(CreateHistogramsForTaggerBinning==kTRUE)
                        hist_SubImCut_fit.Fill(fit, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    else
                        hist_SubImCut_fit.Fill(fit, tagger.GetTagged_t(i));
                }
            }

            if(mm>850 && mm<1025)
            {
                if(CreateHistogramsForTaggerBinning==kTRUE)
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                else
                    hist_MMCut.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTagged_t(i));

                if(fit.IsSolved()==kTRUE)
                {
                    if(fit.ConfidenceLevel()>0.1)
                    {
                        if(CreateHistogramsForTaggerBinning==kTRUE)
                            hist_fit.Fill(fit, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                        else
                            hist_fit.Fill(fit, tagger.GetTagged_t(i));
                    }
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

    folder  = h->GetDirectory("SubIM_Cut");
    hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
        GHistWriteList* subfolder  = folder->GetDirectory("fit");
        hist_SubImCut_fit.PrepareWriteList(subfolder, TString(name).Append("_fit").Data());

    folder  = h->GetDirectory("MM_Cut");
    hist_MMCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
        subfolder  = folder->GetDirectory("fit");
        hist_fit.PrepareWriteList(subfolder, TString(name).Append("_fit").Data());
}

void    GAnalysis3MesonsProton::Reset(Option_t* option)
{
    checkProton.Reset(option);
    hist_raw.Reset(option);
    hist_SubImCut.Reset(option);
    hist_MMCut.Reset(option);
    hist_SubImCut_fit.Reset(option);
    hist_fit.Reset(option);
}

void    GAnalysis3MesonsProton::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    //checkProton.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_raw.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_MMCut.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_SubImCut_fit.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    hist_fit.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}

