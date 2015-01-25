#include "GFit.h"
#include "GTreeMeson.h"


GFit3Constraints::GFit3Constraints(const Bool_t _IsEtap)    :
    solved(kFALSE),
    isEtap(_IsEtap),
    GFit(6, 3)
{

}

GFit3Constraints::~GFit3Constraints()
{

}

void    GFit3Constraints::Set(const TLorentzVector& p0,
                              const TLorentzVector& p1,
                              const TLorentzVector& p2,
                              const TLorentzVector& p3,
                              const TLorentzVector& p4,
                              const TLorentzVector& p5)
{
    solved = kFALSE;
    fitter.Reset();

    GKinFitterParticle  photons[6];

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    for(int i=0; i<6; i++)
    {
        photons[i].SetResolutions(3, 3, 10);
        fitter.AddPosKFParticle(photons[i]);
    }

    Int_t   index[6]    = {0, 1, 2, 3, 4, 5};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
}














GFit4Constraints::GFit4Constraints(const Bool_t _IsEtap)    :
    solved(kFALSE),
    isEtap(_IsEtap),
    GFit(6, 4)
{

}

GFit4Constraints::~GFit4Constraints()
{

}

void    GFit4Constraints::Set(const TLorentzVector& p0,
                              const TLorentzVector& p1,
                              const TLorentzVector& p2,
                              const TLorentzVector& p3,
                              const TLorentzVector& p4,
                              const TLorentzVector& p5,
                              const TLorentzVector& beamAndTarget)
{
    solved = kFALSE;
    fitter.Reset();

    GKinFitterParticle  photons[6];

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    for(int i=0; i<6; i++)
    {
        photons[i].SetResolutions(3, 3, 10);
        fitter.AddPosKFParticle(photons[i]);
    }

    Int_t   index[6]    = {0, 1, 2, 3, 4, 5};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    //fitter.AddSubMissMassConstraint(beamAndTarget, MASS_PROTON);
    fitter.AddSubMissMassConstraint(beamAndTarget, 6, &index[0], MASS_PROTON);
}



















GFit4ConstraintsBeam::GFit4ConstraintsBeam(const Bool_t _IsEtap)    :
    solved(kFALSE),
    isEtap(_IsEtap),
    GFit(7, 4)
{

}

GFit4ConstraintsBeam::~GFit4ConstraintsBeam()
{

}

TLorentzVector  GFit4ConstraintsBeam::GetTotalFitParticle()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i).Get4Vector();
    return ret;
}

void    GFit4ConstraintsBeam::Set(const TLorentzVector& p0,
                              const TLorentzVector& p1,
                              const TLorentzVector& p2,
                              const TLorentzVector& p3,
                              const TLorentzVector& p4,
                              const TLorentzVector& p5,
                              const TLorentzVector& beamAndTarget)
{
    solved = kFALSE;
    fitter.Reset();

    GKinFitterParticle  photons[6];
    GKinFitterParticle  beam;

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    for(int i=0; i<6; i++)
    {
        photons[i].SetResolutions(3, 3, 10);
        fitter.AddPosKFParticle(photons[i]);
    }
    beam.Set4Vector(beamAndTarget);
    beam.SetResolutions(1, 1, 2);
    fitter.AddNegKFParticle(beam);

    Int_t   index[7]    = {0, 1, 2, 3, 4, 5, 6};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    fitter.AddInvMassConstraint(MASS_PROTON);
}














GFit7ConstraintsProton::GFit7ConstraintsProton(const Bool_t _IsEtap)    :
    solved(kFALSE),
    isEtap(_IsEtap),
    GFit(7, 7)
{

}

GFit7ConstraintsProton::~GFit7ConstraintsProton()
{

}

TLorentzVector  GFit7ConstraintsProton::GetTotalFitParticle()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i).Get4Vector();
    return ret;
}

void    GFit7ConstraintsProton::Set(const TLorentzVector& p0,
                                    const TLorentzVector& p1,
                                    const TLorentzVector& p2,
                                    const TLorentzVector& p3,
                                    const TLorentzVector& p4,
                                    const TLorentzVector& p5,
                                    const TLorentzVector& beamAndTarget,
                                    const TLorentzVector& proton)
{
    solved = kFALSE;
    fitter.Reset();

    GKinFitterParticle  photons[6];
    GKinFitterParticle  pro;

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    for(int i=0; i<6; i++)
    {
        photons[i].SetResolutions(3, 3, 10);
        fitter.AddPosKFParticle(photons[i]);
    }
    TLorentzVector  help(proton);
    help.SetE(beamAndTarget.E()-p0.E()-p1.E()-p2.E()-p3.E()-p4.E()-p5.E());
    help.SetVect(proton.Vect().Unit()*(TMath::Sqrt(help.E()*help.E() - MASS_PROTON*MASS_PROTON)));
    pro.Set4Vector(help);
    pro.SetResolutions(3, 3, 25);
    fitter.AddPosKFParticle(pro);

    Int_t   index[7]    = {0, 1, 2, 3, 4, 5, 6};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    fitter.AddTotEnergyConstraint(beamAndTarget.E());
    fitter.AddTotMomentumConstraint(beamAndTarget.Vect());
}
























GFit7ConstraintsBeamProton::GFit7ConstraintsBeamProton(const Bool_t _IsEtap)    :
    solved(kFALSE),
    isEtap(_IsEtap),
    GFit(8, 7)
{

}

GFit7ConstraintsBeamProton::~GFit7ConstraintsBeamProton()
{

}

TLorentzVector  GFit7ConstraintsBeamProton::GetTotalFitParticle()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i).Get4Vector();
    return ret;
}

void    GFit7ConstraintsBeamProton::Set(const TLorentzVector& p0,
                                        const TLorentzVector& p1,
                                        const TLorentzVector& p2,
                                        const TLorentzVector& p3,
                                        const TLorentzVector& p4,
                                        const TLorentzVector& p5,
                                        const TLorentzVector& beamAndTarget,
                                        const TLorentzVector& proton)
{
    solved = kFALSE;
    fitter.Reset();

    GKinFitterParticle  photons[6];
    GKinFitterParticle  pro;
    GKinFitterParticle  beam;

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    for(int i=0; i<6; i++)
    {
        photons[i].SetResolutions(3, 3, 10);
        fitter.AddPosKFParticle(photons[i]);
    };
    beam.Set4Vector(beamAndTarget);
    beam.SetResolutions(1, 1, 2);
    fitter.AddNegKFParticle(beam);
    TLorentzVector  help(proton);
    help.SetE(beamAndTarget.E()-p0.E()-p1.E()-p2.E()-p3.E()-p4.E()-p5.E());
    help.SetVect(proton.Vect().Unit()*(TMath::Sqrt(help.E()*help.E() - MASS_PROTON*MASS_PROTON)));
    //help.Print();
    //std::cout << help.M() << std::endl;
    pro.Set4Vector(help);
    pro.SetResolutions(3, 3, 25);
    fitter.AddPosKFParticle(pro);

    Int_t   index[8]    = {0, 1, 2, 3, 4, 5, 6, 7};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    fitter.AddTotEnergyConstraint(0);
    fitter.AddTotMomentumConstraint(TVector3(0.0, 0.0, 0.0));
}












GHistFit::GHistFit(const char* name, const char* title, const Int_t _NPulls, Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    nPulls(_NPulls),
    nIter(TString(name).Append("_nIter"), TString(title).Append(" nIter"), 10, 0, 10, 48, kFALSE),
    im(TString(name).Append("_im"), TString(title).Append(" inv. Mass"), 2000, 0, 2000, 48, kFALSE),
    sub0im(TString(name).Append("_sub0im"), TString(title).Append(" sub0 inv. Mass"), 800, 0, 800, 48, kFALSE),
    sub1im(TString(name).Append("_sub1im"), TString(title).Append(" sub1 inv. Mass"), 400, 0, 400, 48, kFALSE),
    sub2im(TString(name).Append("_sub2im"), TString(title).Append(" sub2 inv. Mass"), 400, 0, 400, 48, kFALSE),
    theta(TString(name).Append("_theta"), TString(title).Append(" theta"), 180, 0, 180, 48, kFALSE),
    phi(TString(name).Append("_phi"), TString(title).Append(" phi"), 360, -180, 180, 48, kFALSE),
    Pim(TString(name).Append("_Pim"), TString(title).Append(" Proton Mass"), 2000, 0, 2000, 48, kFALSE),
    Ptheta(TString(name).Append("_Ptheta"), TString(title).Append(" Proton theta"), 180, 0, 180, 48, kFALSE),
    Pphi(TString(name).Append("_Pphi"), TString(title).Append(" Proton phi"), 360, -180, 180, 48, kFALSE),
    chiSq(TString(name).Append("_ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 100, 48, kFALSE),
    confidenceLevel(TString(name).Append("_ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, 48, kFALSE),
    pulls(TString(name).Append("_Pulls"), TString(title).Append(" Pulls"), 100, -5, 5, nPulls, 0, nPulls, kFALSE)
{

}

GHistFit::~GHistFit()
{

}

void        GHistFit::CalcResult()
{
    nIter.CalcResult();
    im.CalcResult();
    sub0im.CalcResult();
    sub1im.CalcResult();
    sub2im.CalcResult();
    theta.CalcResult();
    phi.CalcResult();
    Pim.CalcResult();
    Ptheta.CalcResult();
    Pphi.CalcResult();
    chiSq.CalcResult();
    confidenceLevel.CalcResult();
    pulls.CalcResult();
}

Int_t       GHistFit::Fill(GFit& fitter, const Double_t taggerTime)
{
    //nIter.Fill(fitter.GetNIter(), taggerTime);
    nIter.Fill(0, taggerTime);
    TLorentzVector  etap(fitter.GetTotalFitParticle());
    im.Fill(etap.M(), taggerTime);
    sub0im.Fill(fitter.GetSub(0).M(), taggerTime);
    sub1im.Fill(fitter.GetSub(1).M(), taggerTime);
    sub2im.Fill(fitter.GetSub(2).M(), taggerTime);
    theta.Fill(etap.Theta()*TMath::RadToDeg(), taggerTime);
    phi.Fill(etap.Phi()*TMath::RadToDeg(), taggerTime);
    Pim.Fill(fitter.GetRecoil().M(), taggerTime);
    Ptheta.Fill(fitter.GetRecoil().Theta()*TMath::RadToDeg(), taggerTime);
    Pphi.Fill(fitter.GetRecoil().Phi()*TMath::RadToDeg(), taggerTime);
    chiSq.Fill(fitter.GetChi2(), taggerTime);
    confidenceLevel.Fill(fitter.ConfidenceLevel(), taggerTime);
    for(int i=0; i<nPulls; i++)
        pulls.Fill(fitter.GetPull(i), i);
}

Int_t       GHistFit::Fill(GFit& fitter, const Double_t taggerTime, const Int_t taggerChannel)
{
    //nIter.Fill(fitter.GetNIter(), taggerTime, taggerChannel);
    nIter.Fill(0, taggerTime);
    TLorentzVector  etap(fitter.GetTotalFitParticle());
    im.Fill(etap.M(), taggerTime, taggerChannel);
    sub0im.Fill(fitter.GetSub(0).M(), taggerTime, taggerChannel);
    sub1im.Fill(fitter.GetSub(1).M(), taggerTime, taggerChannel);
    sub2im.Fill(fitter.GetSub(2).M(), taggerTime, taggerChannel);
    theta.Fill(etap.Theta()*TMath::RadToDeg(), taggerTime, taggerChannel);
    phi.Fill(etap.Phi()*TMath::RadToDeg(), taggerTime, taggerChannel);
    Pim.Fill(fitter.GetRecoil().M(), taggerTime, taggerChannel);
    Ptheta.Fill(fitter.GetRecoil().Theta()*TMath::RadToDeg(), taggerTime, taggerChannel);
    Pphi.Fill(fitter.GetRecoil().Phi()*TMath::RadToDeg(), taggerTime, taggerChannel);
    chiSq.Fill(fitter.GetChi2(), taggerTime, taggerChannel);
    confidenceLevel.Fill(fitter.ConfidenceLevel(), taggerTime, taggerChannel);
    for(int i=0; i<nPulls; i++)
        pulls.Fill(fitter.GetPull(i), i);
}

void    GHistFit::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        nIter.PrepareWriteList(arr, TString(name).Append("_nIter").Data());
        im.PrepareWriteList(arr, TString(name).Append("_IM").Data());
        sub0im.PrepareWriteList(arr, TString(name).Append("_Sub0IM").Data());
        sub1im.PrepareWriteList(arr, TString(name).Append("_Sub1IM").Data());
        sub2im.PrepareWriteList(arr, TString(name).Append("_Sub2IM").Data());
        theta.PrepareWriteList(arr, TString(name).Append("_Theta").Data());
        phi.PrepareWriteList(arr, TString(name).Append("_Phi").Data());
        Pim.PrepareWriteList(arr, TString(name).Append("_PrIM").Data());
        Ptheta.PrepareWriteList(arr, TString(name).Append("_PrTheta").Data());
        Pphi.PrepareWriteList(arr, TString(name).Append("_PrPhi").Data());
        chiSq.PrepareWriteList(arr, TString(name).Append("_ChiSq").Data());
        confidenceLevel.PrepareWriteList(arr, TString(name).Append("_ConfLev").Data());
        pulls.PrepareWriteList(arr, TString(name).Append("_Pulls").Data());
    }
    else
    {
        nIter.PrepareWriteList(arr);
        im.PrepareWriteList(arr);
        sub0im.PrepareWriteList(arr);
        sub1im.PrepareWriteList(arr);
        sub2im.PrepareWriteList(arr);
        theta.PrepareWriteList(arr);
        phi.PrepareWriteList(arr);
        Pim.PrepareWriteList(arr);
        Ptheta.PrepareWriteList(arr);
        Pphi.PrepareWriteList(arr);
        chiSq.PrepareWriteList(arr);
        confidenceLevel.PrepareWriteList(arr);
        pulls.PrepareWriteList(arr);
    }
}

void        GHistFit::Reset(Option_t* option)
{
    nIter.Reset(option);
    im.Reset(option);
    sub0im.Reset(option);
    sub1im.Reset(option);
    sub2im.Reset(option);
    theta.Reset(option);
    phi.Reset(option);
    Pim.Reset(option);
    Ptheta.Reset(option);
    Pphi.Reset(option);
    chiSq.Reset(option);
    confidenceLevel.Reset(option);
    pulls.Reset(option);
}

void        GHistFit::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    nIter.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub0im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub1im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub2im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    theta.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    phi.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    Pim.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    Ptheta.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    Pphi.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    chiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    confidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pulls.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}
