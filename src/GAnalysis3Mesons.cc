#include "GAnalysis3Mesons.h"
#include "GTreeTagger.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"


GAnalysis3Mesons::GAnalysis3Mesons(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(_IsEtap),
    hist_count(TString(name).Append("count"), TString(title).Append("count"), 5, 0, 5, kFALSE),
    hist_raw(TString(name).Append("raw"), TString(title).Append("raw"), kFALSE),
    hist_raw_SubAll(TString(name).Append("raw_SubAll"), TString(title).Append("raw_SubAll"), 1000, 0, 1000, kFALSE),
    hist_SubImCut(TString(name).Append("SubImCut"), TString(title).Append("SubImCut"), kFALSE),
    hist_SubImCut_SubAll(TString(name).Append("SubImCut_SubAll"), TString(title).Append("SubImCut_SubAll"), 1000, 0, 1000, kFALSE),
    hist_MMCut(TString(name).Append("MMCut"), TString(title).Append("MMCut"), kFALSE),
    hist_SubAll(TString(name).Append("SubAll"), TString(title).Append("SubAll"), 1000, 0, 1000, kFALSE),
//    fit1("fit1", kFALSE),
//    hist_fit1(TString(name).Append("hist_fit1"), TString(title).Append("hist_fit1"), kFALSE),
//    fit1Vertex("fit1Vertex", kFALSE),
//    hist_fit1Vertex(TString(name).Append("hist_fit1Vertex"), TString(title).Append("hist_fit1Vertex"), kFALSE),
//    fit3("fit3", kFALSE),
//    hist_fit3(TString(name).Append("hist_fit3"), TString(title).Append("hist_fit3"), kFALSE),
//    fit3Vertex("fit3Vertex", kFALSE),
//    hist_fit3Vertex(TString(name).Append("hist_fit3Vertex"), TString(title).Append("hist_fit3Vertex"), kFALSE),
    fit4("fit4", kFALSE),
    hist_fit4(TString(name).Append("hist_fit4"), TString(title).Append("hist_fit4"), kFALSE),
    hist_fit4_SubAll(TString(name).Append("fit4_SubAll"), TString(title).Append("fit4_SubAll"), 1000, 0, 1000, kFALSE)
//    fit4Vertex("fit4Vertex", kFALSE),
//    hist_fit4Vertex(TString(name).Append("hist_fit4Vertex"), TString(title).Append("hist_fit4Vertex"), kFALSE),
//    fitBeam1("fitBeam1", kFALSE),
//    hist_fitBeam1(TString(name).Append("hist_fitBeam1"), TString(title).Append("hist_fitBeam1"), kFALSE),
//    fitBeam1Vertex("fitBeam1Vertex", kFALSE),
//    hist_fitBeam1Vertex(TString(name).Append("hist_fitBeam1Vertex"), TString(title).Append("hist_fitBeam1Vertex"), kFALSE),
//    fitBeam3("fitBeam3", kFALSE),
//    hist_fitBeam3(TString(name).Append("hist_fitBeam3"), TString(title).Append("hist_fitBeam3"), kFALSE),
//    fitBeam3Vertex("fitBeam3Vertex", kFALSE),
//    hist_fitBeam3Vertex(TString(name).Append("hist_fitBeam3Vertex"), TString(title).Append("hist_fitBeam3Vertex"), kFALSE),
//    fitBeam4("fitBeam4", kFALSE),
//    hist_fitBeam4(TString(name).Append("hist_fitBeam4"), TString(title).Append("hist_fitBeam4"), kFALSE),
//    fitBeam4Vertex("fitBeam4Vertex", kFALSE),
//    hist_fitBeam4Vertex(TString(name).Append("hist_fitBeam4Vertex"), TString(title).Append("hist_fitBeam4Vertex"), kFALSE)
{
    if(_IsEtap==kTRUE)
        SetCutSubIM(0, 510, 600);
    else
        SetCutSubIM(0, 110, 150);
    SetCutSubIM(1, 110, 150);
    SetCutSubIM(2, 110, 150);

    SetCutMM(870, 1000);

//    fit1.AddConstraintMM();
//    fit1Vertex.AddConstraintMM();
//    fit3.AddConstraintsIM();
//    fit3Vertex.AddConstraintsIM();
    fit4.AddConstraintMM();
    fit4.AddConstraintsIM();
//    fit4Vertex.AddConstraintMM();
//    fit4Vertex.AddConstraintsIM();

//    fitBeam1.AddConstraintMM();
//    fitBeam1Vertex.AddConstraintMM();
//    fitBeam3.AddConstraintsIM();
//    fitBeam3Vertex.AddConstraintsIM();
//    fitBeam4.AddConstraintMM();
//    fitBeam4.AddConstraintsIM();
//    fitBeam4Vertex.AddConstraintMM();
//    fitBeam4Vertex.AddConstraintsIM();
}

GAnalysis3Mesons::~GAnalysis3Mesons()
{

}

void   GAnalysis3Mesons::CalcResult()
{
    hist_count.CalcResult();
    hist_raw.CalcResult();
    hist_raw_SubAll.CalcResult();
    hist_SubImCut.CalcResult();
    hist_SubImCut_SubAll.CalcResult();
    hist_MMCut.CalcResult();
    hist_SubAll.CalcResult();

//    fit1.CalcResult();
//    hist_fit1.CalcResult();
//    fit1Vertex.CalcResult();
//    hist_fit1Vertex.CalcResult();
//    fit3.CalcResult();
//    hist_fit3.CalcResult();
//    fit3Vertex.CalcResult();
//    hist_fit3Vertex.CalcResult();
    fit4.CalcResult();
    hist_fit4.CalcResult();
    hist_fit4_SubAll.CalcResult();
//    fit4Vertex.CalcResult();
//    hist_fit4Vertex.CalcResult();

//    fitBeam1.CalcResult();
//    hist_fitBeam1.CalcResult();
//    fitBeam1Vertex.CalcResult();
//    hist_fitBeam1Vertex.CalcResult();
//    fitBeam3.CalcResult();
//    hist_fitBeam3.CalcResult();
//    fitBeam3Vertex.CalcResult();
//    hist_fitBeam3Vertex.CalcResult();
//    fitBeam4.CalcResult();
//    hist_fitBeam4.CalcResult();
//    fitBeam4Vertex.CalcResult();
//    hist_fitBeam4Vertex.CalcResult();
}

void    GAnalysis3Mesons::Fill(const GTreeMeson& meson, GTreeParticle& photons, const GTreeTagger& tagger, const GTreeA2Geant& geantTree)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    theta  = meson.Particle(0).Theta()*TMath::RadToDeg();
    Double_t    phi  = meson.Particle(0).Phi()*TMath::RadToDeg();
    Double_t    mm;
    Double_t    sub_im_0    = (photons.Particle(0) + photons.Particle(1)).M();
    Double_t    sub_im_1    = (photons.Particle(2) + photons.Particle(3)).M();
    Double_t    sub_im_2    = (photons.Particle(4) + photons.Particle(5)).M();

    hist_count.Fill(0);

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
        TLorentzVector  helpCM(meson.Particle(0));
        helpCM.Boost(-tagger.GetVectorProtonTarget(i).BoostVector());

        hist_raw.Fill(im, mm, theta, phi, helpCM.Theta()*TMath::RadToDeg(), sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
        for(int k=0; k<6; k++)
        {
            for(int l=k+1; l<6; l++)
                hist_raw_SubAll.Fill((photons.Particle(k)+photons.Particle(l)).M(), tagger.GetTaggedTime(i));
        }
    }

    //cout << cutSubIM[0] << "   " << cutSubIM[1] << "   " << cutSubIM[2] << "   " << cutSubIM[3] << "   " << cutSubIM[4] << "   " << cutSubIM[5] << "   " << cutMM[6] << "   " << cutMM[7] << endl;
    if((sub_im_0>cutSubIM[0] && sub_im_0<cutSubIM[1]) &&
        (sub_im_1>cutSubIM[2] && sub_im_1<cutSubIM[3]) &&
        (sub_im_2>cutSubIM[4] && sub_im_2<cutSubIM[5]))
    {
        hist_count.Fill(1);
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
            TLorentzVector  helpCM(meson.Particle(0));
            helpCM.Boost(-tagger.GetVectorProtonTarget(i).BoostVector());

            hist_SubImCut.Fill(im, mm, theta, phi, helpCM.Theta()*TMath::RadToDeg(), sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            for(int k=0; k<6; k++)
            {
                for(int l=k+1; l<6; l++)
                    hist_SubImCut_SubAll.Fill((photons.Particle(k)+photons.Particle(l)).M(), tagger.GetTaggedTime(i));
            }

            if(mm>cutMM[0] && mm<cutMM[1])
            {
                hist_count.Fill(2);
                hist_MMCut.Fill(im, mm, theta, phi, helpCM.Theta()*TMath::RadToDeg(), sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                for(int k=0; k<6; k++)
                {
                    for(int l=k+1; l<6; l++)
                        hist_SubAll.Fill((photons.Particle(k)+photons.Particle(l)).M(), tagger.GetTaggedTime(i));
                }

//                fit1.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fit1.SetBeam(tagger.GetTaggedEnergy(i));
//                fit1Vertex.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fit1Vertex.SetBeam(tagger.GetTaggedEnergy(i));
//                fit3.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fit3Vertex.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
                fit4.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
                fit4.SetBeam(tagger.GetTaggedEnergy(i));
//                fit4Vertex.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fit4Vertex.SetBeam(tagger.GetTaggedEnergy(i));

//                fitBeam1.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitBeam1.SetBeam(tagger.GetTaggedEnergy(i));
//                fitBeam1Vertex.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitBeam1Vertex.SetBeam(tagger.GetTaggedEnergy(i));
//                fitBeam3.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitBeam3Vertex.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitBeam4.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitBeam4.SetBeam(tagger.GetTaggedEnergy(i));
//                fitBeam4Vertex.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitBeam4Vertex.SetBeam(tagger.GetTaggedEnergy(i));

//                if(fit1.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fit1.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fit1Vertex.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fit1Vertex.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fit3.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fit3.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fit3Vertex.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fit3Vertex.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                if(fit4.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i), geantTree))
                {
                    hist_count.Fill(3);
                    hist_fit4.Fill(im, mm, theta, phi, helpCM.Theta()*TMath::RadToDeg(), sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                    for(int k=0; k<6; k++)
                    {
                        for(int l=k+1; l<6; l++)
                            hist_fit4_SubAll.Fill((photons.Particle(k)+photons.Particle(l)).M(), tagger.GetTaggedTime(i));

                        photons.AddParticle(fit4.GetFittedPhoton(k).Ek, fit4.GetFittedPhoton(k).Theta, fit4.GetFittedPhoton(k).Phi,
                                            0.0, photons.GetTime(k), photons.GetClusterSize(k), photons.GetCentralCrystal(k), photons.GetCentralVeto(k), photons.GetDetectors(k),
                                            photons.GetVetoEnergy(k), photons.GetMWPC0Energy(k), photons.GetMWPC1Energy(k), photons.GetTrackIndex(k));
                    }
                }
//                if(fit4Vertex.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fit4Vertex.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));

//                if(fitBeam1.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitBeam1.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fitBeam1Vertex.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitBeam1Vertex.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fitBeam3.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitBeam3.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fitBeam3Vertex.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitBeam3Vertex.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fitBeam4.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitBeam4.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fitBeam4Vertex.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitBeam4Vertex.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            }
        }
    }
}

void    GAnalysis3Mesons::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* h  = arr->GetDirectory("WithoutProton");
    hist_count.PrepareWriteList(h, "count");

    GHistWriteList* folder  = h->GetDirectory("Raw");
    hist_raw.PrepareWriteList(folder, TString(name).Append("_Raw").Data());
    hist_raw_SubAll.PrepareWriteList(folder, "SubAll");

    folder  = h->GetDirectory("SubIM_Cut");
    hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
    hist_SubImCut_SubAll.PrepareWriteList(folder, "SubAll");

    folder  = h->GetDirectory("MM_Cut");
    hist_MMCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
    hist_SubAll.PrepareWriteList(folder, "SubAll");

//    folder  = h->GetDirectory("fit1");
//    fit1.PrepareWriteList(folder);
//    hist_fit1.PrepareWriteList(folder, "fit1");
//    folder  = folder->GetDirectory("Vertex");
//    fit1Vertex.PrepareWriteList(folder);
//    hist_fit1Vertex.PrepareWriteList(folder, "fit1Vertex");

//    folder  = h->GetDirectory("fit3");
//    fit3.PrepareWriteList(folder);
//    hist_fit3.PrepareWriteList(folder, "fit3");
//    folder  = folder->GetDirectory("Vertex");
//    fit3Vertex.PrepareWriteList(folder);
//    hist_fit3Vertex.PrepareWriteList(folder, "fit3Vertex");

    folder  = h->GetDirectory("fit4");
    fit4.PrepareWriteList(folder);
    hist_fit4.PrepareWriteList(folder, "fit4");
    hist_fit4_SubAll.PrepareWriteList(folder, "SubAll");
//    folder  = folder->GetDirectory("Vertex");
//    fit4Vertex.PrepareWriteList(folder);
//    hist_fit4Vertex.PrepareWriteList(folder, "fit4Vertex");

//    folder  = h->GetDirectory("fitBeam1");
//    fitBeam1.PrepareWriteList(folder);
//    hist_fitBeam1.PrepareWriteList(folder, "fitBeam1");
//    folder  = folder->GetDirectory("Vertex");
//    fitBeam1Vertex.PrepareWriteList(folder);
//    hist_fitBeam1Vertex.PrepareWriteList(folder, "fitBeam1Vertex");

//    folder  = h->GetDirectory("fitBeam3");
//    fitBeam3.PrepareWriteList(folder);
//    hist_fitBeam3.PrepareWriteList(folder, "fitBeam3");
//    folder  = folder->GetDirectory("Vertex");
//    fitBeam3Vertex.PrepareWriteList(folder);
//    hist_fitBeam3Vertex.PrepareWriteList(folder, "fitBeam3Vertex");

//    folder  = h->GetDirectory("fitBeam4");
//    fitBeam4.PrepareWriteList(folder);
//    hist_fitBeam4.PrepareWriteList(folder, "fitBeam4");
//    folder  = folder->GetDirectory("Vertex");
//    fitBeam4Vertex.PrepareWriteList(folder);
//    hist_fitBeam4Vertex.PrepareWriteList(folder, "fitBeam4Vertex");
}

void    GAnalysis3Mesons::Reset(Option_t* option)
{
    hist_count.Reset(option);
    hist_raw.Reset(option);
    hist_raw_SubAll.Reset(option);
    hist_SubImCut.Reset(option);
    hist_SubImCut_SubAll.Reset(option);
    hist_MMCut.Reset(option);
    hist_SubAll.Reset(option);

//    fit1.Reset(option);
//    hist_fit1.Reset(option);
//    fit1Vertex.Reset(option);
//    hist_fit1Vertex.Reset(option);
//    fit3.Reset(option);
//    hist_fit3.Reset(option);
//    fit3Vertex.Reset(option);
//    hist_fit3Vertex.Reset(option);
    fit4.Reset(option);
    hist_fit4.Reset(option);
    hist_fit4_SubAll.Reset(option);
//    fit4Vertex.Reset(option);
//    hist_fit4Vertex.Reset(option);

//    fitBeam1.Reset(option);
//    hist_fitBeam1.Reset(option);
//    fitBeam1Vertex.Reset(option);
//    hist_fitBeam1Vertex.Reset(option);
//    fitBeam3.Reset(option);
//    hist_fitBeam3.Reset(option);
//    fitBeam3Vertex.Reset(option);
//    hist_fitBeam3Vertex.Reset(option);
//    fitBeam4.Reset(option);
//    hist_fitBeam4.Reset(option);
//    fitBeam4Vertex.Reset(option);
//    hist_fitBeam4Vertex.Reset(option);
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











GAnalysis3MesonsProton::GAnalysis3MesonsProton(const char* name, const char* title, const Bool_t _IsEtap, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    isEtap(_IsEtap),
    hist_count(TString(name).Append("count"), TString(title).Append("count"), 5, 0, 5, kFALSE),
    checkProton(TString(name).Append("checkProton"), TString(title).Append("checkProton"), kFALSE),
    hist_raw(TString(name).Append("raw"), TString(title).Append("raw"), kFALSE),
    hist_raw_TOF(TString(name).Append("raw_TOF"), TString(title).Append("raw_TOF"), 300, -15, 15, 800, 0, 800, kFALSE),
    hist_raw_SubAll(TString(name).Append("raw_SubAll"), TString(title).Append("raw_SubAll"), 1000, 0, 1000, kFALSE),
    hist_SubImCut(TString(name).Append("SubImCut"), TString(title).Append("SubImCut"), kFALSE),
    hist_SubImCut_TOF(TString(name).Append("SubImCut_TOF"), TString(title).Append("SubImCut_TOF"), 300, -15, 15, 800, 0, 800, kFALSE),
    hist_SubImCut_SubAll(TString(name).Append("SubImCut_SubAll"), TString(title).Append("SubImCut_SubAll"), 1000, 0, 1000, kFALSE),
    hist_MMCut(TString(name).Append("MMCut"), TString(title).Append("MMCut"), kFALSE),
    hist_TOF(TString(name).Append("TOF"), TString(title).Append("TOF"), 300, -15, 15, 800, 0, 800, kFALSE),
    hist_SubAll(TString(name).Append("SubAll"), TString(title).Append("SubAll"), 1000, 0, 1000, kFALSE),
//    fit4("fit4", kFALSE),
//    hist_fit4(TString(name).Append("hist_fit4"), TString(title).Append("hist_fit4"), kFALSE),
//    fit4Vertex("fit4Vertex", kFALSE),
//    hist_fit4Vertex(TString(name).Append("hist_fit4Vertex"), TString(title).Append("hist_fit4Vertex"), kFALSE),
//    fitBeam4("fitBeam4", kFALSE),
//    hist_fitBeam4(TString(name).Append("hist_fitBeam4"), TString(title).Append("hist_fitBeam4"), kFALSE),
//    fitBeam4Vertex("fitBeam4Vertex", kFALSE),
//    hist_fitBeam4Vertex(TString(name).Append("hist_fitBeam4Vertex"), TString(title).Append("hist_fitBeam4Vertex"), kFALSE),
    fitProton6("fitProton6", kFALSE),
    hist_fitProton6(TString(name).Append("hist_fitProton6"), TString(title).Append("hist_fitProton6"), kFALSE),
    hist_fitProton6_SubAll(TString(name).Append("fitProton6_SubAll"), TString(title).Append("fitProton6_SubAll"), 1000, 0, 1000, kFALSE)
//    fitProton6Vertex("fitProton6Vertex", kFALSE),
//    hist_fitProton6Vertex(TString(name).Append("hist_fitProton6Vertex"), TString(title).Append("hist_fitProton6Vertex"), kFALSE),
//    fitBeamProton6("fitBeamProton6", kFALSE),
//    hist_fitBeamProton6(TString(name).Append("hist_fitBeamProton6"), TString(title).Append("hist_fitBeamProton6"), kFALSE),
//    fitBeamProton6Vertex("fitBeamProton6Vertex", kFALSE),
//    hist_fitBeamProton6Vertex(TString(name).Append("hist_fitBeamProton6Vertex"), TString(title).Append("hist_fitBeamProton6Vertex"), kFALSE)
{
    if(_IsEtap==kTRUE)
        SetCutSubIM(0, 510, 600);
    else
        SetCutSubIM(0, 110, 150);
    SetCutSubIM(1, 110, 150);
    SetCutSubIM(2, 110, 150);

    SetCutMM(870, 1000);


//    fit4.AddConstraintMM();
//    fit4.AddConstraintsIM();
//    fit4Vertex.AddConstraintMM();
//    fit4Vertex.AddConstraintsIM();

//    fitBeam4.AddConstraintMM();
//    fitBeam4.AddConstraintsIM();
//    fitBeam4Vertex.AddConstraintMM();
//    fitBeam4Vertex.AddConstraintsIM();

    fitProton6.AddConstraintsTotMomentum();
    fitProton6.AddConstraintsTotEnergy();
    fitProton6.AddConstraintsIM();
//    fitProton6Vertex.AddConstraintsTotMomentum();
//    fitProton6Vertex.AddConstraintsTotEnergy();
//    fitProton6Vertex.AddConstraintsIM();

//    fitBeamProton6.AddConstraintsTotMomentum();
//    fitBeamProton6.AddConstraintsTotEnergy();
//    fitBeamProton6.AddConstraintsIM();
//    fitBeamProton6Vertex.AddConstraintsTotMomentum();
//    fitBeamProton6Vertex.AddConstraintsTotEnergy();
//    fitBeamProton6Vertex.AddConstraintsIM();
}

GAnalysis3MesonsProton::~GAnalysis3MesonsProton()
{

}

void   GAnalysis3MesonsProton::CalcResult()
{
    hist_count.CalcResult();
    checkProton.CalcResult();
    hist_raw.CalcResult();
    hist_raw_TOF.CalcResult();
    hist_raw_SubAll.CalcResult();
    hist_SubImCut.CalcResult();
    hist_SubImCut_TOF.CalcResult();
    hist_SubImCut_SubAll.CalcResult();
    hist_MMCut.CalcResult();
    hist_TOF.CalcResult();
    hist_SubAll.CalcResult();

//    fit4.CalcResult();
//    hist_fit4.CalcResult();
//    fit4Vertex.CalcResult();
//    hist_fit4Vertex.CalcResult();

//    fitBeam4.CalcResult();
//    hist_fitBeam4.CalcResult();
//    fitBeam4Vertex.CalcResult();
//    hist_fitBeam4Vertex.CalcResult();

    fitProton6.CalcResult();
    hist_fitProton6.CalcResult();
    hist_fitProton6_SubAll.CalcResult();
//    fitProton6Vertex.CalcResult();
//    hist_fitProton6Vertex.CalcResult();

//    fitBeamProton6.CalcResult();
//    hist_fitBeamProton6.CalcResult();
//    fitBeamProton6Vertex.CalcResult();
//    hist_fitBeamProton6Vertex.CalcResult();
}

void    GAnalysis3MesonsProton::Fill(const GTreeMeson& meson, GTreeParticle& photons, GTreeParticle& proton, const GTreeTagger& tagger, const GTreeA2Geant& geantTree)
{
    Double_t    im  = meson.Particle(0).M();
    Double_t    theta  = meson.Particle(0).Theta()*TMath::RadToDeg();
    Double_t    phi  = meson.Particle(0).Phi()*TMath::RadToDeg();
    Double_t    protonTheta  = proton.Particle(0).Theta()*TMath::RadToDeg();
    Double_t    protonPhi  = proton.Particle(0).Phi()*TMath::RadToDeg();
    Double_t    mm;
    Double_t    sub_im_0    = (photons.Particle(0) + photons.Particle(1)).M();
    Double_t    sub_im_1    = (photons.Particle(2) + photons.Particle(3)).M();
    Double_t    sub_im_2    = (photons.Particle(4) + photons.Particle(5)).M();

    hist_count.Fill(0);

    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
        TLorentzVector  helpCM(meson.Particle(0));
        helpCM.Boost(-tagger.GetVectorProtonTarget(i).BoostVector());
        TLorentzVector  helpProtonCM(meson.Particle(0));
        helpProtonCM.Boost(-tagger.GetVectorProtonTarget(i).BoostVector());

        hist_raw.Fill(im, mm, theta, phi, helpCM.Theta()*TMath::RadToDeg(), proton.Particle(0).E(), protonTheta, protonPhi, helpProtonCM.Theta()*TMath::RadToDeg(), sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
        hist_raw_TOF.Fill(tagger.GetTaggedTime(i)-proton.GetTime(0), proton.GetClusterEnergy(0), tagger.GetTaggedTime(i));
        for(int k=0; k<6; k++)
        {
            for(int l=k+1; l<6; l++)
                hist_raw_SubAll.Fill((photons.Particle(k)+photons.Particle(l)).M(), tagger.GetTaggedTime(i));
        }
    }

    if((sub_im_0>cutSubIM[0] && sub_im_0<cutSubIM[1]) &&
        (sub_im_1>cutSubIM[2] && sub_im_1<cutSubIM[3]) &&
        (sub_im_2>cutSubIM[4] && sub_im_2<cutSubIM[5]))
    {
        hist_count.Fill(1);
        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            if(checkProton.Check(meson, proton, tagger.GetVectorProtonTarget(i), tagger.GetTaggedTime(i))==kFALSE)
                continue;

            mm  = (tagger.GetVectorProtonTarget(i)-meson.Particle(0)).M();
            TLorentzVector  helpCM(meson.Particle(0));
            helpCM.Boost(-tagger.GetVectorProtonTarget(i).BoostVector());
            TLorentzVector  helpProtonCM(meson.Particle(0));
            helpProtonCM.Boost(-tagger.GetVectorProtonTarget(i).BoostVector());

            hist_SubImCut.Fill(im, mm, theta, phi, helpCM.Theta()*TMath::RadToDeg(), proton.Particle(0).E(), protonTheta, protonPhi, helpProtonCM.Theta()*TMath::RadToDeg(), sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            hist_SubImCut_TOF.Fill(tagger.GetTaggedTime(i)-proton.GetTime(0), proton.GetClusterEnergy(0), tagger.GetTaggedTime(i));
            for(int k=0; k<6; k++)
            {
                for(int l=k+1; l<6; l++)
                    hist_SubImCut_SubAll.Fill((photons.Particle(k)+photons.Particle(l)).M(), tagger.GetTaggedTime(i));
            }

            if(mm>cutMM[0] && mm<cutMM[1])
            {
                hist_count.Fill(2);
                hist_MMCut.Fill(im, mm, theta, phi, helpCM.Theta()*TMath::RadToDeg(), proton.Particle(0).E(), protonTheta, protonPhi, helpProtonCM.Theta()*TMath::RadToDeg(), sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                hist_TOF.Fill(tagger.GetTaggedTime(i)-proton.GetTime(0), proton.GetClusterEnergy(0), tagger.GetTaggedTime(i));
                for(int k=0; k<6; k++)
                {
                    for(int l=k+1; l<6; l++)
                        hist_SubAll.Fill((photons.Particle(k)+photons.Particle(l)).M(), tagger.GetTaggedTime(i));
                }

//                fit4.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fit4.SetBeam(tagger.GetTaggedEnergy(i));
//                fit4Vertex.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fit4Vertex.SetBeam(tagger.GetTaggedEnergy(i));

//                fitBeam4.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitBeam4.SetBeam(tagger.GetTaggedEnergy(i));
//                fitBeam4Vertex.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitBeam4Vertex.SetBeam(tagger.GetTaggedEnergy(i));

                fitProton6.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
                fitProton6.SetBeam(tagger.GetTaggedEnergy(i));
                fitProton6.SetProton(proton.Particle(0));
//                fitProton6Vertex.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitProton6Vertex.SetBeam(tagger.GetTaggedEnergy(i));
//                fitProton6Vertex.SetProton(proton.Particle(0));

//                fitBeamProton6.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitBeamProton6.SetBeam(tagger.GetTaggedEnergy(i));
//                fitBeamProton6.SetProton(proton.Particle(0));
//                fitBeamProton6Vertex.Set(photons.Particle(0), photons.Particle(1), photons.Particle(2), photons.Particle(3), photons.Particle(4), photons.Particle(5));
//                fitBeamProton6Vertex.SetBeam(tagger.GetTaggedEnergy(i));
//                fitBeamProton6Vertex.SetProton(proton.Particle(0));

//                if(fit4.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fit4.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fit4Vertex.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fit4Vertex.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));

//                if(fitBeam4.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitBeam4.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fitBeam4Vertex.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitBeam4Vertex.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));

                if(fitProton6.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i), geantTree))
                {
                    hist_count.Fill(3);
                    hist_fitProton6.Fill(im, mm, theta, phi, helpCM.Theta()*TMath::RadToDeg(), proton.Particle(0).E(), protonTheta, protonPhi, helpProtonCM.Theta()*TMath::RadToDeg(), sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
                    for(int k=0; k<6; k++)
                    {
                        for(int l=k+1; l<6; l++)
                            hist_fitProton6_SubAll.Fill((photons.Particle(k)+photons.Particle(l)).M(), tagger.GetTaggedTime(i));

                        photons.AddParticle(fitProton6.GetFittedPhoton(k).Ek, fitProton6.GetFittedPhoton(k).Theta, fitProton6.GetFittedPhoton(k).Phi,
                                            0, photons.GetTime(k), photons.GetClusterSize(k), photons.GetCentralCrystal(k), photons.GetCentralVeto(k), photons.GetDetectors(k),
                                            photons.GetVetoEnergy(k), photons.GetMWPC0Energy(k), photons.GetMWPC1Energy(k), photons.GetTrackIndex(k));
                    }
                    proton.AddParticle(fitProton6.GetFittedProton().Ek, fitProton6.GetFittedProton().Theta, fitProton6.GetFittedProton().Phi,
                                        MASS_PROTON, proton.GetTime(0), proton.GetClusterSize(0), proton.GetCentralCrystal(0), proton.GetCentralVeto(0), proton.GetDetectors(0),
                                        proton.GetVetoEnergy(0), proton.GetMWPC0Energy(0), proton.GetMWPC1Energy(0), proton.GetTrackIndex(0));
                }
//                if(fitProton6Vertex.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitProton6Vertex.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));

//                if(fitBeamProton6.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitBeamProton6.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
//                if(fitBeamProton6Vertex.Solve(tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i)))
//                    hist_fitBeamProton6Vertex.Fill(im, mm, sub_im_0, sub_im_1, sub_im_2, tagger.GetTaggedTime(i), tagger.GetTaggedChannel(i));
            }
        }
    }
}

void    GAnalysis3MesonsProton::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistWriteList* h  = arr->GetDirectory("WithProton");
    hist_count.PrepareWriteList(h, "count");

    GHistWriteList* folder  = h->GetDirectory("CheckProton");
    checkProton.PrepareWriteList(folder, TString(name).Append("_CheckProton").Data());

    folder  = h->GetDirectory("Raw");
    hist_raw.PrepareWriteList(folder, TString(name).Append("_Raw").Data());
    hist_raw_TOF.PrepareWriteList(folder, TString(name).Append("_Raw_TOF").Data());
    hist_raw_SubAll.PrepareWriteList(folder, "SubAll");

    folder  = h->GetDirectory("SubIM_Cut");
    hist_SubImCut.PrepareWriteList(folder, TString(name).Append("_subIMCut").Data());
    hist_SubImCut_TOF.PrepareWriteList(folder, TString(name).Append("_subIMCut_TOF").Data());
    hist_SubImCut_SubAll.PrepareWriteList(folder, "SubAll");

    folder  = h->GetDirectory("MM_Cut");
    hist_TOF.PrepareWriteList(folder, TString(name).Append("_TOF").Data());
    hist_MMCut.PrepareWriteList(folder, TString(name).Append("_MMCut").Data());
    hist_SubAll.PrepareWriteList(folder, "SubAll");

//    folder  = h->GetDirectory("fit4");
//    fit4.PrepareWriteList(folder);
//    hist_fit4.PrepareWriteList(folder, "fit4");
//    folder  = folder->GetDirectory("Vertex");
//    fit4Vertex.PrepareWriteList(folder);
//    hist_fit4Vertex.PrepareWriteList(folder, "fit4Vertex");

//    folder  = h->GetDirectory("fitBeam4");
//    fitBeam4.PrepareWriteList(folder);
//    hist_fitBeam4.PrepareWriteList(folder, "fitBeam4");
//    folder  = folder->GetDirectory("Vertex");
//    fitBeam4Vertex.PrepareWriteList(folder);
//    hist_fitBeam4Vertex.PrepareWriteList(folder, "fitBeam4Vertex");

    folder  = h->GetDirectory("fitProton6");
    fitProton6.PrepareWriteList(folder);
    hist_fitProton6.PrepareWriteList(folder, "fitProton6");
    hist_fitProton6_SubAll.PrepareWriteList(folder, "SubAll");
//    folder  = folder->GetDirectory("Vertex");
//    fitProton6Vertex.PrepareWriteList(folder);
//    hist_fitProton6Vertex.PrepareWriteList(folder, "fitProton6Vertex");

//    folder  = h->GetDirectory("fitBeamProton6");
//    fitBeamProton6.PrepareWriteList(folder);
//    hist_fitBeamProton6.PrepareWriteList(folder, "fitBeamProton6");
//    folder  = folder->GetDirectory("Vertex");
//    fitBeamProton6Vertex.PrepareWriteList(folder);
//    hist_fitBeamProton6Vertex.PrepareWriteList(folder, "fitBeamProton6Vertex");
}

void    GAnalysis3MesonsProton::Reset(Option_t* option)
{
    hist_count.Reset(option);
    checkProton.Reset(option);
    hist_raw.Reset(option);
    hist_raw_TOF.Reset(option);
    hist_raw_SubAll.Reset(option);
    hist_SubImCut.Reset(option);
    hist_SubImCut_TOF.Reset(option);
    hist_SubImCut_SubAll.Reset(option);
    hist_MMCut.Reset(option);
    hist_TOF.Reset(option);
    hist_SubAll.Reset(option);

//    fit4.Reset(option);
//    hist_fit4.Reset(option);
//    fit4Vertex.Reset(option);
//    hist_fit4Vertex.Reset(option);

//    fitBeam4.Reset(option);
//    hist_fitBeam4.Reset(option);
//    fitBeam4Vertex.Reset(option);
//    hist_fitBeam4Vertex.Reset(option);

    fitProton6.Reset(option);
    hist_fitProton6.Reset(option);
    hist_fitProton6_SubAll.Reset(option);
//    fitProton6Vertex.Reset(option);
//    hist_fitProton6Vertex.Reset(option);

//    fitBeamProton6.Reset(option);
//    hist_fitBeamProton6.Reset(option);
//    fitBeamProton6Vertex.Reset(option);
//    hist_fitBeamProton6Vertex.Reset(option);
}

void    GAnalysis3MesonsProton::SetCutSubIM(const Int_t subNumber, const Double_t min, const Double_t max)
{
    cutSubIM[2*subNumber] = min;
    cutSubIM[(2*subNumber)+1] = max;
}

void    GAnalysis3MesonsProton::SetCutMM(const Double_t min, const Double_t max)
{
    cutMM[0] = min;
    cutMM[1] = max;
}
