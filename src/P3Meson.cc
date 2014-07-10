#include "P3Meson.h"


P3Meson::P3Meson(const TString &_Name, const Bool_t _IsEtap)    :
    name(_Name),
    isEtap(_IsEtap),
    raw(TString(name).Append("_raw")),
    cutIM({{110,155},{110,155},{110,155}}),
    cutIMevent(TString(name).Append("_cutIM")),
    cutMM({850,1050}),
    cutMMevent(TString(name).Append("_cutMM")),
    fit3Con(6,3,0),
    hist_fit3Con(TString(name).Append("_fit3Con")),
    cutFit3ConConfidenceLevel(0.1),
    hist_fit3Con_cutCL(TString(name).Append("_fit3ConCutCL")),
    fit4Con(6,4,0),
    hist_fit4Con(TString(name).Append("_fit4Con")),
    cutFit4ConConfidenceLevel(0.1),
    hist_fit4Con_cutCL(TString(name).Append("_fit4ConCutCL")),
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

Bool_t  P3Meson::ReconstructEvent(const GTreeMeson& meson, const GTreeTagger& tagger)
{
    if(meson.GetNParticles()>0)
    {
        bool        passIM;
        bool        passFit3Con;
        bool        found   = false;

        im          = meson.Particle(0).M();
        imSub[0]    = (meson.SubParticles(0, 0)+meson.SubParticles(0, 1)).M();
        imSub[1]    = (meson.SubParticles(0, 2)+meson.SubParticles(0, 3)).M();
        imSub[2]    = (meson.SubParticles(0, 4)+meson.SubParticles(0, 5)).M();

        if((imSub[0]>cutIM[0][0] && imSub[0]<cutIM[0][1]) && (imSub[1]>cutIM[1][0] && imSub[1]<cutIM[1][1]) && (imSub[2]>cutIM[2][0] && imSub[2]<cutIM[2][1]))
            passIM  = true;
        else
            passIM  = false;

        if(passIM)
        {
            passFit3Con = DoFit3Con(meson);
            im_fit  = fittedMeson.M();
            imSub_fit[0]    = (fittedSubParticles[0]+fittedSubParticles[1]).M();
            imSub_fit[1]    = (fittedSubParticles[2]+fittedSubParticles[3]).M();
            imSub_fit[2]    = (fittedSubParticles[4]+fittedSubParticles[5]).M();
        }


        raw.Fill(im, meson.Particle(0), tagger);
        raw.FillSubMesons(imSub[0], imSub[1], imSub[2], tagger);

        if(passIM)
        {
            cutIMevent.Fill(im, meson.Particle(0), tagger);
            cutIMevent.FillSubMesons(imSub[0], imSub[1], imSub[2], tagger);

            if(passFit3Con)
            {
                hist_fit3Con.Fill(im_fit, fittedMeson, tagger);
                hist_fit3Con.FillSubMesons(imSub_fit[0], imSub_fit[1], imSub_fit[2], tagger);
                for(int i=0; i<24; i++)
                    Pull[i]   = fit3Con.Pull(i);
                hist_fit3Con.FillFit(fit3Con.GetChi2(), conLevel, Pull, tagger);

                if(conLevel>=cutFit3ConConfidenceLevel)
                {
                    hist_fit3Con_cutCL.Fill(im_fit, fittedMeson, tagger);
                    hist_fit3Con_cutCL.FillSubMesons(imSub_fit[0], imSub_fit[1], imSub_fit[2], tagger);
                    for(int i=0; i<24; i++)
                        Pull[i]   = fit3Con.Pull(i);
                    hist_fit3Con_cutCL.FillFit(fit3Con.GetChi2(), conLevel, Pull, tagger);

                    found   = true;
                }

            }
        }

        if(found)
            return kTRUE;
    }

    return kFALSE;
}

Bool_t  P3Meson::ReconstructTagger(const GTreeMeson& meson, const GTreeTagger& tagger)
{
    Double_t    misMass;
    for(int i=0; i<tagger.GetNTagged(); i++)
    {
        misMass = (tagger.GetVectorProtonTarget(i) - fittedMeson).M();
        if(misMass>cutMM[0] && misMass<cutMM[1])
        {
            cutMMevent.Fill(im_fit, misMass, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            cutMMevent.FillSubMesons(imSub_fit[0], imSub_fit[1], imSub_fit[2], tagger.GetTagged_t(i), tagger.GetTagged_ch(i));

            DoFit4Con(meson, tagger.GetVectorProtonTarget(i));
            im_fit  = fittedMeson.M();
            imSub_fit[0]    = (fittedSubParticles[0]+fittedSubParticles[1]).M();
            imSub_fit[1]    = (fittedSubParticles[2]+fittedSubParticles[3]).M();
            imSub_fit[2]    = (fittedSubParticles[4]+fittedSubParticles[5]).M();

            hist_fit4Con.Fill(im_fit, misMass, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            hist_fit4Con.FillSubMesons(imSub_fit[0], imSub_fit[1], imSub_fit[2], tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            for(int i=0; i<24; i++)
                Pull[i]   = fit4Con.Pull(i);
            hist_fit4Con.FillFit(fit4Con.GetChi2(), conLevel, Pull, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
            if(conLevel>=cutFit4ConConfidenceLevel)
            {
                hist_fit4Con_cutCL.Fill(im_fit, misMass, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                hist_fit4Con_cutCL.FillSubMesons(imSub_fit[0], imSub_fit[1], imSub_fit[2], tagger.GetTagged_t(i), tagger.GetTagged_ch(i));
                for(int i=0; i<24; i++)
                    Pull[i]   = fit4Con.Pull(i);
                hist_fit4Con_cutCL.FillFit(fit4Con.GetChi2(), conLevel, Pull, tagger.GetTagged_t(i), tagger.GetTagged_ch(i));

                nFound++;
            }
        }
    }
}

Bool_t	P3Meson::ProcessEvent(const GTreeMeson& meson, const GTreeTagger& tagger)
{
    if(!ReconstructEvent(meson, tagger))
        return kFALSE;

    return ReconstructTagger(meson, tagger);
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

    conLevel    = fit3Con.ConfidenceLevel();
    fittedMeson = fit3Con.GetTotalFitParticle().Get4Vector();
    fittedSubParticles[0] = fit3Con.GetParticle(0).Get4Vector();
    fittedSubParticles[1] = fit3Con.GetParticle(1).Get4Vector();
    fittedSubParticles[2] = fit3Con.GetParticle(2).Get4Vector();
    fittedSubParticles[3] = fit3Con.GetParticle(3).Get4Vector();
    fittedSubParticles[4] = fit3Con.GetParticle(4).Get4Vector();
    fittedSubParticles[5] = fit3Con.GetParticle(5).Get4Vector();

    return kTRUE;
}

Bool_t  P3Meson::DoFit4Con(const GTreeMeson& meson, const TLorentzVector& beamPlusTarget)
{
    if(!fitInit(meson, fit4Con))
        return kFALSE;

    Int_t   help[6] = {0, 1, 2, 3, 4, 5};
    fit4Con.AddSubMissMassConstraint(beamPlusTarget, 6, help, 938.272046);

    if(fit4Con.Solve()<0)
        return kFALSE;

    conLevel    = fit4Con.ConfidenceLevel();
    fittedMeson = fit4Con.GetTotalFitParticle().Get4Vector();
    fittedSubParticles[0] = fit4Con.GetParticle(0).Get4Vector();
    fittedSubParticles[1] = fit4Con.GetParticle(1).Get4Vector();
    fittedSubParticles[2] = fit4Con.GetParticle(2).Get4Vector();
    fittedSubParticles[3] = fit4Con.GetParticle(3).Get4Vector();
    fittedSubParticles[4] = fit4Con.GetParticle(4).Get4Vector();
    fittedSubParticles[5] = fit4Con.GetParticle(5).Get4Vector();

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

