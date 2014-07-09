#include "P3Meson.h"


P3Meson::P3Meson(const TString &_Name, const Bool_t _IsEtap)    :
    name(_Name),
    isEtap(_IsEtap),
    raw(TString(name).Append("_raw")),
    cutIM({{110,155},{110,155},{110,155}}),
    cutIMevent(TString(name).Append("_cutIM")),
    cutMM({850,1050}),
    cutMMevent(TString(name).Append("_cutMM")),
    cutFit3ConConfidenceLevel(0.1),
    fit3Con(6,3,0),
    hist_fit3Con(TString(name).Append("_fit3Con")),
    cutFit4ConConfidenceLevel(0.1),
    fit4Con(6,4,0),
    hist_fit4Con(TString(name).Append("_fit4Con")),
    nFound(0)
{
    if(isEtap)
    {
        cutIM[0][0] = 500;
        cutIM[0][1] = 590;
    }

    GammaResFile   = new TFile("configfiles/FitGammaResolution.root");
    GammaEloss     = (TH2F*)GammaResFile->Get("Eloss");
    GammaERes      = (TH2F*)GammaResFile->Get("EResIter");
    GammaThetaRes  = (TH2F*)GammaResFile->Get("ThetaRes;1");
    GammaPhiRes    = (TH2F*)GammaResFile->Get("PhiRes;1");
}

P3Meson::~P3Meson()
{

}


Bool_t	P3Meson::ProcessEvent(const GTreeMeson& meson, const GTreeTagger& tagger)
{
    //if(GetEventNumber() == 0) nFound = 0;
    //else if(GetEventNumber() % 100000 == 0) cout << "Event: "<< GetEventNumber() << " Total Etas found: " << nFound << endl;

    bool        passIM;
    bool        passFit3Con;

    if(meson.GetNParticles()>0)
    {
        Double_t    im  = meson.Particle(0).M();
        imSub[0]    = (meson.SubParticles(0, 0)+meson.SubParticles(0, 1)).M();
        imSub[1]    = (meson.SubParticles(0, 2)+meson.SubParticles(0, 3)).M();
        imSub[2]    = (meson.SubParticles(0, 4)+meson.SubParticles(0, 5)).M();
        Double_t Pull[24];

        if((imSub[0]>cutIM[0][0] && imSub[0]<cutIM[0][1]) && (imSub[1]>cutIM[1][0] && imSub[1]<cutIM[1][1]) && (imSub[2]>cutIM[2][0] && imSub[2]<cutIM[2][1]))
            passIM  = true;
        else
            passIM  = false;

        passFit3Con = DoFit3Con(meson);

        for(int i=0; i<tagger.GetNTagged(); i++)
        {
            misMass = (tagger.GetVector(i)+TLorentzVector(0,0,0,MASS_PROTON) - meson.Meson(0)).M();
            raw.Fill(im, misMass, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            raw.FillSubMesons(imSub[0], imSub[1], imSub[2], tagger.GetTagged_t(i), tagger.GetTagged_ch(i));

            if(passIM)
            {
                cutIMevent.Fill(im, misMass, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                cutIMevent.FillSubMesons(imSub[0], imSub[1], imSub[2], tagger.GetTagged_t(i), tagger.GetTagged_ch(i));

                if(passFit3Con)
                {
                    hist_fit3Con.Fill(im, misMass, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    hist_fit3Con.FillSubMesons(imSub[0], imSub[1], imSub[2], tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                    for(int i=0; i<24; i++)
                        Pull[i]   = fit3Con.Pull(i);
                    hist_fit3Con.FillFit(fit3Con.GetChi2(), conLevel, Pull, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));

                    if(misMass>cutMM[0] && misMass<cutMM[1])
                    {
                        cutMMevent.Fill(im, misMass, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                        cutMMevent.FillSubMesons(imSub[0], imSub[1], imSub[2], tagger.GetTagged_t(i), tagger.GetTagged_ch(i));

                        nFound++;
                    }
                }
            }
        }
    }
}

Bool_t    P3Meson::fitInit(const GTreeMeson& meson, GKinFitter &fitter)
{
    fitter.Reset();

    Int_t Ebin  = 0;
    Int_t Thbin = 0;
    Float_t resth = 0;
    Float_t resph = 0;
    Float_t resE  = 0;

    fitter.Reset();
    for(int i=0; i<6; i++)
    {
        Ebin  = GammaEloss->GetXaxis()->FindFixBin(meson.SubPhotons(0,i).E());
        Thbin = GammaEloss->GetYaxis()->FindFixBin(meson.SubPhotons(0,i).Theta()*TMath::RadToDeg());
        // Get resolutions
        resth = GammaThetaRes->GetBinContent(Ebin, Thbin);
        resph = GammaPhiRes->GetBinContent(Ebin, Thbin);
        resE  = GammaERes->GetBinContent(Ebin, Thbin);
        if(resth==0 || resph==0 || resE==0 ) return kFALSE; // If energy or angle is out of calibrated  range!
        // Now set particle parameters
        //                     LorentzVector
        pho[i].Set4Vector(meson.SubPhotons(0,i));
        //std::cout << "Res: " << resth << ", " << resph << ", " << photons->Particle(i).E()*resE << std::endl;
        pho[i].SetResolutions(1.5* resth, 1.5* resph, 2.2 *meson.SubPhotons(0,i).E()*resE);
        fitter.AddPosKFParticle(pho[i]);
    }

    Int_t	sub[2];
    sub[0]	= 0;
    sub[1]	= 1;
    if(isEtap)
        fitter.AddSubInvMassConstraint(2, sub, MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, sub, MASS_PI0);
    sub[0]	= 2;
    sub[1]	= 3;
    fitter.AddSubInvMassConstraint(2, sub, MASS_PI0);
    sub[0]	= 4;
    sub[1]	= 5;
    fitter.AddSubInvMassConstraint(2, sub, MASS_PI0);

    return kTRUE;
}

Bool_t    P3Meson::DoFit3Con(const GTreeMeson& meson)
{
    if(!fitInit(meson, fit3Con))
        return kFALSE;

    if(fit3Con.Solve()<0)
        return kFALSE;

    //hist_fit3Con.Fill(meson.Meson(0).M(), misMass, tagger_time, tagger_channel);
    //hist_fit3Con.FillSubMesons(imSub[0], imSub[1], imSub[2], tagger_time, tagger_channel);

    conLevel   = fit3Con.ConfidenceLevel();
    //Double_t Pull[24];

    /*for(int i=0; i<24; i++)
        Pull[i]   = fit3Con.Pull(i);

    hist_fit3Con.FillFit(chiSq, conLevel, Pull, tagger_time, tagger_channel);*/

    return kTRUE;
}


Bool_t 	P3Meson::Write(TDirectory& curDir)
{
    curDir.cd();

    raw.Write(curDir);
    cutIMevent.Write(curDir);
    hist_fit3Con.Write(curDir);
    cutMMevent.Write(curDir);
    hist_fit4Con.Write(curDir);

	return kTRUE;
}


void    P3Meson::Clear()
{
    nFound  =0;

    raw.Clear();
    cutIMevent.Clear();
    hist_fit3Con.Clear();
    cutMMevent.Clear();
    hist_fit4Con.Clear();
}

