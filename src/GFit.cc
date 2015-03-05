#include "GFit.h"
#include "GTreeMeson.h"


GFit3Constraints::GFit3Constraints(const Bool_t _IsEtap)    :
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
    photons[0].SetResolutions(3, 3/sin(p0.Theta()), 0.02*TMath::Power(p0.E(), 0.36));
    photons[1].SetResolutions(3, 3/sin(p1.Theta()), 0.02*TMath::Power(p1.E(), 0.36));
    photons[2].SetResolutions(3, 3/sin(p2.Theta()), 0.02*TMath::Power(p2.E(), 0.36));
    photons[3].SetResolutions(3, 3/sin(p3.Theta()), 0.02*TMath::Power(p3.E(), 0.36));
    photons[4].SetResolutions(3, 3/sin(p4.Theta()), 0.02*TMath::Power(p4.E(), 0.36));
    photons[5].SetResolutions(3, 3/sin(p5.Theta()), 0.02*TMath::Power(p5.E(), 0.36));
    for(int i=0; i<6; i++)
        fitter.AddPosKFParticle(photons[i]);

    Int_t   index[6]    = {0, 1, 2, 3, 4, 5};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
}














GFit4Constraints::GFit4Constraints(const Bool_t _IsEtap)    :
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
                              const TLorentzVector& _BeamAndTarget)
{
    beamAndTarget  = _BeamAndTarget;

    solved = kFALSE;
    fitter.Reset();

    GKinFitterParticle  photons[6];

    photons[0].Set4Vector(p0);
    photons[1].Set4Vector(p1);
    photons[2].Set4Vector(p2);
    photons[3].Set4Vector(p3);
    photons[4].Set4Vector(p4);
    photons[5].Set4Vector(p5);
    photons[0].SetResolutions(3, 3/sin(p0.Theta()), 0.02*TMath::Power(p0.E(), 0.36));
    photons[1].SetResolutions(3, 3/sin(p1.Theta()), 0.02*TMath::Power(p1.E(), 0.36));
    photons[2].SetResolutions(3, 3/sin(p2.Theta()), 0.02*TMath::Power(p2.E(), 0.36));
    photons[3].SetResolutions(3, 3/sin(p3.Theta()), 0.02*TMath::Power(p3.E(), 0.36));
    photons[4].SetResolutions(3, 3/sin(p4.Theta()), 0.02*TMath::Power(p4.E(), 0.36));
    photons[5].SetResolutions(3, 3/sin(p5.Theta()), 0.02*TMath::Power(p5.E(), 0.36));
    for(int i=0; i<6; i++)
        fitter.AddPosKFParticle(photons[i]);

    Int_t   index[6]    = {0, 1, 2, 3, 4, 5};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    fitter.AddSubMissMassConstraint(beamAndTarget, 6, &index[0], MASS_PROTON);
}



















GFit4ConstraintsBeam::GFit4ConstraintsBeam(const Bool_t _IsEtap)    :
    isEtap(_IsEtap),
    GFit(7, 4)
{

}

GFit4ConstraintsBeam::~GFit4ConstraintsBeam()
{

}

TLorentzVector  GFit4ConstraintsBeam::GetMeson()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i);
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
    photons[0].SetResolutions(3, 3/sin(p0.Theta()), 0.02*TMath::Power(p0.E(), 0.36));
    photons[1].SetResolutions(3, 3/sin(p1.Theta()), 0.02*TMath::Power(p1.E(), 0.36));
    photons[2].SetResolutions(3, 3/sin(p2.Theta()), 0.02*TMath::Power(p2.E(), 0.36));
    photons[3].SetResolutions(3, 3/sin(p3.Theta()), 0.02*TMath::Power(p3.E(), 0.36));
    photons[4].SetResolutions(3, 3/sin(p4.Theta()), 0.02*TMath::Power(p4.E(), 0.36));
    photons[5].SetResolutions(3, 3/sin(p5.Theta()), 0.02*TMath::Power(p5.E(), 0.36));
    for(int i=0; i<6; i++)
        fitter.AddPosKFParticle(photons[i]);
    beam.Set4Vector(beamAndTarget);
    beam.SetResolutions(1, 1, 0.1);
    fitter.AddNegKFParticle(beam);

    Int_t   index[6]    = {0, 1, 2, 3, 4, 5};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    fitter.AddInvMassConstraint(MASS_PROTON);
}














GFit7ConstraintsProton::GFit7ConstraintsProton(const Bool_t _IsEtap)    :
    isEtap(_IsEtap),
    GFit(7, 7)
{

}

GFit7ConstraintsProton::~GFit7ConstraintsProton()
{

}

TLorentzVector  GFit7ConstraintsProton::GetMeson()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i);
    return ret;
}

void    GFit7ConstraintsProton::Set(const TLorentzVector& p0,
                                    const TLorentzVector& p1,
                                    const TLorentzVector& p2,
                                    const TLorentzVector& p3,
                                    const TLorentzVector& p4,
                                    const TLorentzVector& p5,
                                    const TLorentzVector& _BeamAndTarget,
                                    const TLorentzVector& proton)
{
    beamAndTarget = _BeamAndTarget;

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
    photons[0].SetResolutions(3, 3/sin(p0.Theta()), 0.02*TMath::Power(p0.E(), 0.36));
    photons[1].SetResolutions(3, 3/sin(p1.Theta()), 0.02*TMath::Power(p1.E(), 0.36));
    photons[2].SetResolutions(3, 3/sin(p2.Theta()), 0.02*TMath::Power(p2.E(), 0.36));
    photons[3].SetResolutions(3, 3/sin(p3.Theta()), 0.02*TMath::Power(p3.E(), 0.36));
    photons[4].SetResolutions(3, 3/sin(p4.Theta()), 0.02*TMath::Power(p4.E(), 0.36));
    photons[5].SetResolutions(3, 3/sin(p5.Theta()), 0.02*TMath::Power(p5.E(), 0.36));
    for(int i=0; i<6; i++)
        fitter.AddPosKFParticle(photons[i]);
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
    isEtap(_IsEtap),
    GFit(8, 7)
{

}

GFit7ConstraintsBeamProton::~GFit7ConstraintsBeamProton()
{

}

TLorentzVector  GFit7ConstraintsBeamProton::GetMeson()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i);
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
    photons[0].SetResolutions(3, 3/sin(p0.Theta()), 0.02*TMath::Power(p0.E(), 0.36));
    photons[1].SetResolutions(3, 3/sin(p1.Theta()), 0.02*TMath::Power(p1.E(), 0.36));
    photons[2].SetResolutions(3, 3/sin(p2.Theta()), 0.02*TMath::Power(p2.E(), 0.36));
    photons[3].SetResolutions(3, 3/sin(p3.Theta()), 0.02*TMath::Power(p3.E(), 0.36));
    photons[4].SetResolutions(3, 3/sin(p4.Theta()), 0.02*TMath::Power(p4.E(), 0.36));
    photons[5].SetResolutions(3, 3/sin(p5.Theta()), 0.02*TMath::Power(p5.E(), 0.36));
    for(int i=0; i<6; i++)
        fitter.AddPosKFParticle(photons[i]);
    beam.Set4Vector(beamAndTarget);
    beam.SetResolutions(1, 1, 0.1);
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
    im(TString(name).Append("im"), TString(title).Append(" inv. Mass"), 2000, 0, 2000, 48, kFALSE),
    sub0im(TString(name).Append("sub0im"), TString(title).Append(" sub0 inv. Mass"), 800, 0, 800, 48, kFALSE),
    sub1im(TString(name).Append("sub1im"), TString(title).Append(" sub1 inv. Mass"), 400, 0, 400, 48, kFALSE),
    sub2im(TString(name).Append("sub2im"), TString(title).Append(" sub2 inv. Mass"), 400, 0, 400, 48, kFALSE),
    theta(TString(name).Append("theta"), TString(title).Append(" theta"), 180, 0, 180, 48, kFALSE),
    phi(TString(name).Append("phi"), TString(title).Append(" phi"), 360, -180, 180, 48, kFALSE),
    Pim(TString(name).Append("Pim"), TString(title).Append(" Proton Mass"), 2000, 0, 2000, 48, kFALSE),
    Ptheta(TString(name).Append("Ptheta"), TString(title).Append(" Proton theta"), 180, 0, 180, 48, kFALSE),
    Pphi(TString(name).Append("Pphi"), TString(title).Append(" Proton phi"), 360, -180, 180, 48, kFALSE),
    chiSq(TString(name).Append("ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 100000, 48, kFALSE),
    VchiSq(TString(name).Append("VChiSq"), TString(title).Append(" VChiSq"), 1000, 0, 100000, 48, kFALSE),
    CchiSq(TString(name).Append("CChiSq"), TString(title).Append(" CChiSq"), 1000, 0, 100, 48, kFALSE),
    confidenceLevel(TString(name).Append("ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, 48, kFALSE),
    VconfidenceLevel(TString(name).Append("VConfLev"), TString(title).Append(" VConfLev"), 1000, 0, 1, 48, kFALSE),
    CconfidenceLevel(TString(name).Append("CConfLev"), TString(title).Append(" CConfLev"), 1000, 0, 1, 48, kFALSE),
    pulls(TString(name).Append("Pulls"), TString(title).Append(" Pulls"), 100, -5, 5, nPulls, 0, nPulls, kFALSE)
{

}

GHistFit::~GHistFit()
{

}

void        GHistFit::CalcResult()
{
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
    VchiSq.CalcResult();
    CchiSq.CalcResult();
    confidenceLevel.CalcResult();
    VconfidenceLevel.CalcResult();
    CconfidenceLevel.CalcResult();
    pulls.CalcResult();
}

Int_t       GHistFit::Fill(GFit& fitter, const Double_t taggerTime)
{
    TLorentzVector  etap(fitter.GetMeson());
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
    VchiSq.Fill(fitter.GetVariablesChi2(), taggerTime);
    CchiSq.Fill(fitter.GetConstraintsChi2(), taggerTime);
    confidenceLevel.Fill(fitter.ConfidenceLevel(), taggerTime);
    VconfidenceLevel.Fill(fitter.VariablesConfidenceLevel(), taggerTime);
    CconfidenceLevel.Fill(fitter.ConstraintsConfidenceLevel(), taggerTime);
    for(int i=0; i<nPulls; i++)
        pulls.Fill(fitter.GetPull(i), i);
}

Int_t       GHistFit::Fill(GFit& fitter, const Double_t taggerTime, const Int_t taggerChannel)
{

    TLorentzVector  etap(fitter.GetMeson());
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
    VchiSq.Fill(fitter.GetVariablesChi2(), taggerTime, taggerChannel);
    CchiSq.Fill(fitter.GetConstraintsChi2(), taggerTime, taggerChannel);
    confidenceLevel.Fill(fitter.ConfidenceLevel(), taggerTime, taggerChannel);
    VconfidenceLevel.Fill(fitter.VariablesConfidenceLevel(), taggerTime, taggerChannel);
    CconfidenceLevel.Fill(fitter.ConstraintsConfidenceLevel(), taggerTime, taggerChannel);
    for(int i=0; i<nPulls; i++)
        pulls.Fill(fitter.GetPull(i), i);
}

void    GHistFit::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
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
        VchiSq.PrepareWriteList(arr, TString(name).Append("_VChiSq").Data());
        CchiSq.PrepareWriteList(arr, TString(name).Append("_CChiSq").Data());
        confidenceLevel.PrepareWriteList(arr, TString(name).Append("_ConfLev").Data());
        VconfidenceLevel.PrepareWriteList(arr, TString(name).Append("_VConfLev").Data());
        CconfidenceLevel.PrepareWriteList(arr, TString(name).Append("_CConfLev").Data());
        pulls.PrepareWriteList(arr, TString(name).Append("_Pulls").Data());
    }
    else
    {
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
        VchiSq.PrepareWriteList(arr);
        CchiSq.PrepareWriteList(arr);
        confidenceLevel.PrepareWriteList(arr);
        VconfidenceLevel.PrepareWriteList(arr);
        CconfidenceLevel.PrepareWriteList(arr);
        pulls.PrepareWriteList(arr);
    }
}

void        GHistFit::Reset(Option_t* option)
{
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
    VchiSq.Reset(option);
    CchiSq.Reset(option);
    confidenceLevel.Reset(option);
    VconfidenceLevel.Reset(option);
    CconfidenceLevel.Reset(option);
    pulls.Reset(option);
}

void        GHistFit::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
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
    VchiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CchiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    confidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    VconfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CconfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pulls.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}









GHistIterativeFit::GHistIterativeFit(const char* name, const char* title, const Int_t _NPulls, const Int_t _NSteps, Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    im(TString(name).Append("im"), TString(title).Append(" inv. Mass"), 750, 500, 1250, _NSteps, 0, _NSteps, kFALSE),
    sub0im(TString(name).Append("sub0im"), TString(title).Append(" sub0 inv. Mass"), 1000, 500, 600, _NSteps, 0, _NSteps, kFALSE),
    sub1im(TString(name).Append("sub1im"), TString(title).Append(" sub1 inv. Mass"), 1000, 85, 185, _NSteps, 0, _NSteps, kFALSE),
    sub2im(TString(name).Append("sub2im"), TString(title).Append(" sub2 inv. Mass"), 1000, 85, 185, _NSteps, 0, _NSteps, kFALSE),
    mm(TString(name).Append("mm"), TString(title).Append(" mm"), 10500, -50, 1000, _NSteps, 0, _NSteps, kFALSE),
    totE(TString(name).Append("totE"), TString(title).Append(" totE"), 1000, -100, 100, _NSteps, 0, _NSteps, kFALSE),
    totPx(TString(name).Append("totPx"), TString(title).Append(" totPx"), 1000, -100, 100, _NSteps, 0, _NSteps, kFALSE),
    totPy(TString(name).Append("totPy"), TString(title).Append(" totPy"), 1000, -100, 100, _NSteps, 0, _NSteps, kFALSE),
    totPz(TString(name).Append("totPz"), TString(title).Append(" totPz"), 1000, -100, 100, _NSteps, 0, _NSteps, kFALSE),
    VchiSq(TString(name).Append("VChiSq"), TString(title).Append(" VChiSq"), 1000, 0, 1000000, _NSteps, 0, _NSteps, kFALSE),
    VconfidenceLevel(TString(name).Append("VConfLev"), TString(title).Append(" VConfLev"), 1000, 0, 1, _NSteps, 0, _NSteps, kFALSE),
    CchiSq(TString(name).Append("CChiSq"), TString(title).Append(" CChiSq"), 1000, 0, 100, _NSteps, 0, _NSteps, kFALSE),
    CconfidenceLevel(TString(name).Append("CConfLev"), TString(title).Append(" CConfLev"), 1000, 0, 1, _NSteps, 0, _NSteps, kFALSE),
    chiSq(TString(name).Append("ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 1000000, _NSteps, 0, _NSteps, kFALSE),
    confidenceLevel(TString(name).Append("ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, _NSteps, 0, _NSteps, kFALSE),
    final(TString(name).Append("Final"), TString(title).Append(" Final"), _NPulls, kFALSE)
{

}

GHistIterativeFit::~GHistIterativeFit()
{

}


void        GHistIterativeFit::CalcResult()
{
    //std::cout << im.GetName() << std::endl;
    im.CalcResult();
    sub0im.CalcResult();
    sub1im.CalcResult();
    sub2im.CalcResult();
    mm.CalcResult();
    totE.CalcResult();
    totPx.CalcResult();
    totPy.CalcResult();
    totPz.CalcResult();
    VchiSq.CalcResult();
    VconfidenceLevel.CalcResult();
    CchiSq.CalcResult();
    CconfidenceLevel.CalcResult();
    chiSq.CalcResult();
    confidenceLevel.CalcResult();

    final.CalcResult();
}

Int_t       GHistIterativeFit::Fill(GFit& fitter)
{
    TLorentzVector  etap(fitter.GetMeson());
    TLorentzVector  tot(fitter.GetTotal());
    im.Fill(etap.M(), fitter.GetIterations());
    sub0im.Fill(fitter.GetSub(0).M(), fitter.GetIterations());
    sub1im.Fill(fitter.GetSub(1).M(), fitter.GetIterations());
    sub2im.Fill(fitter.GetSub(2).M(), fitter.GetIterations());
    mm.Fill(tot.M(), fitter.GetIterations());
    totE.Fill(tot.E(), fitter.GetIterations());
    totPx.Fill(tot.Px(), fitter.GetIterations());
    totPy.Fill(tot.Py(), fitter.GetIterations());
    totPz.Fill(tot.Pz(), fitter.GetIterations());
    VchiSq.Fill(fitter.GetVariablesChi2(), fitter.GetIterations());
    VconfidenceLevel.Fill(fitter.VariablesConfidenceLevel(), fitter.GetIterations());
    CchiSq.Fill(fitter.GetConstraintsChi2(), fitter.GetIterations());
    CconfidenceLevel.Fill(fitter.ConstraintsConfidenceLevel(), fitter.GetIterations());
    chiSq.Fill(fitter.GetChi2(), fitter.GetIterations());
    confidenceLevel.Fill(fitter.ConfidenceLevel(), fitter.GetIterations());
}

void    GHistIterativeFit::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

        im.PrepareWriteList(arr, "IM");
        sub0im.PrepareWriteList(arr, "_Sub0IM");
        sub1im.PrepareWriteList(arr, "_Sub1IM");
        sub2im.PrepareWriteList(arr, "_Sub2IM");
        mm.PrepareWriteList(arr, "_MM");
        totE.PrepareWriteList(arr, "_TotE");
        totPx.PrepareWriteList(arr, "_TotPx");
        totPy.PrepareWriteList(arr, "_TotPy");
        totPz.PrepareWriteList(arr, "_TotPz");
        VchiSq.PrepareWriteList(arr, "_VChiSq");
        VconfidenceLevel.PrepareWriteList(arr, "_VConfLev");
        CchiSq.PrepareWriteList(arr, "_CChiSq");
        CconfidenceLevel.PrepareWriteList(arr, "_CConfLev");
        chiSq.PrepareWriteList(arr, "_ChiSq");
        confidenceLevel.PrepareWriteList(arr, "_ConfLev");

    GHistWriteList* finalDir    = arr->GetDirectory("Final");
    final.PrepareWriteList(finalDir, "Final");
}

void        GHistIterativeFit::Reset(Option_t* option)
{
    im.Reset(option);
    sub0im.Reset(option);
    sub1im.Reset(option);
    sub2im.Reset(option);
    mm.Reset(option);
    totE.Reset(option);
    totPx.Reset(option);
    totPy.Reset(option);
    totPz.Reset(option);
    VchiSq.Reset(option);
    VconfidenceLevel.Reset(option);
    CchiSq.Reset(option);
    CconfidenceLevel.Reset(option);
    chiSq.Reset(option);
    confidenceLevel.Reset(option);

    final.Reset(option);
}

void        GHistIterativeFit::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub0im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub1im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    sub2im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    mm.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    totE.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    totPx.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    totPy.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    totPz.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    VchiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    VconfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CchiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    CconfidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    chiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    confidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);

    final.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}
