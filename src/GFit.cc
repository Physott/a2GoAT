#include "GFit.h"
#include "GTreeMeson.h"

/*
GFit3Constraints::GFit3Constraints(const Bool_t _IsEtap)    :
    solved(kFALSE),
    isEtap(_IsEtap),
    fitter(6, 3, 0)
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
    fitter(6, 4, 0)
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
    fitter.AddSubMissMassConstraint(beamAndTarget, 6, &index[0], MASS_PROTON);
}



















GFit4ConstraintsBeam::GFit4ConstraintsBeam(const Bool_t _IsEtap)    :
    solved(kFALSE),
    isEtap(_IsEtap),
    fitter(7, 4, 0)
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

    Int_t   index[6]    = {0, 1, 2, 3, 4, 5};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    fitter.AddInvMassConstraint(MASS_PROTON);
}














GFit4ConstraintsProton::GFit4ConstraintsProton(const Bool_t _IsEtap)    :
    solved(kFALSE),
    isEtap(_IsEtap),
    fitter(7, 4, 0)
{

}

GFit4ConstraintsProton::~GFit4ConstraintsProton()
{

}

TLorentzVector  GFit4ConstraintsProton::GetTotalFitParticle()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i).Get4Vector();
    return ret;
}

void    GFit4ConstraintsProton::Set(const TLorentzVector& p0,
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
    pro.SetResolutions(3, 3, 0.1);
    fitter.AddPosKFParticle(pro);

    Int_t   index[7]    = {0, 1, 2, 3, 4, 5, 6};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    fitter.AddSubMissMassConstraint(beamAndTarget, 7, &index[0], 0);
}
























GFit4ConstraintsBeamProton::GFit4ConstraintsBeamProton(const Bool_t _IsEtap)    :
    solved(kFALSE),
    isEtap(_IsEtap),
    fitter(8, 5, 0)
{

}

GFit4ConstraintsBeamProton::~GFit4ConstraintsBeamProton()
{

}

TLorentzVector  GFit4ConstraintsBeamProton::GetTotalFitParticle()
{
    TLorentzVector ret(0, 0, 0, 0);
    for(int i=0; i<6; i++)
        ret += fitter.GetParticle(i).Get4Vector();
    return ret;
}

void    GFit4ConstraintsBeamProton::Set(const TLorentzVector& p0,
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
    pro.SetResolutions(3, 3, 0.1);
    fitter.AddPosKFParticle(pro);

    Int_t   index[8]    = {0, 1, 2, 3, 4, 5, 6, 7};
    if(isEtap==kTRUE)
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_ETA);
    else
        fitter.AddSubInvMassConstraint(2, &index[0], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[2], MASS_PI0);
    fitter.AddSubInvMassConstraint(2, &index[4], MASS_PI0);
    //fitter.AddSubInvMassConstraint(1, &index[7], MASS_PROTON);
    fitter.AddTotEnergyConstraint(0);
    //fitter.AddTotMomentumConstraint(TVector3(0, 0, 0));
    fitter.AddInvMassConstraint(0);
}



*/








GHistFit::GHistFit(const char* name, const char* title, const Int_t _NPulls, Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    nPulls(_NPulls),
    im(TString(name).Append("_im"), TString(title).Append(" inv. Mass"), 2000, 0, 2000, 48, kFALSE),
    im_sub0(TString(name).Append("_im_sub0"), TString(title).Append(" inv. Mass SubPart0"), 2000, 0, 2000, 48, kFALSE),
    im_sub1(TString(name).Append("_im_sub1"), TString(title).Append(" inv. Mass SubPart1"), 2000, 0, 2000, 48, kFALSE),
    im_sub2(TString(name).Append("_im_sub2"), TString(title).Append(" inv. Mass SubPart2"), 2000, 0, 2000, 48, kFALSE),
    chiSq(TString(name).Append("_ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 100, 48, kFALSE),
    confidenceLevel(TString(name).Append("_ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, 48, kFALSE),
    vertex_X(TString(name).Append("_Vertex_X"), TString(title).Append(" Vertex_X"), 1000, -1, 1, 48, kFALSE),
    vertex_Y(TString(name).Append("_Vertex_Y"), TString(title).Append(" Vertex_Y"), 1000, -1, 1, 48, kFALSE),
    vertex_Z(TString(name).Append("_Vertex_Z"), TString(title).Append(" Vertex_Z"), 1000, -1, 1, 48, kFALSE),
    vertex_sub0_X(TString(name).Append("_Vertex_sub0_X"), TString(title).Append(" Vertex_sub0_X"), 1000, -1, 1, 48, kFALSE),
    vertex_sub0_Y(TString(name).Append("_Vertex_sub0_Y"), TString(title).Append(" Vertex_sub0_Y"), 1000, -1, 1, 48, kFALSE),
    vertex_sub0_Z(TString(name).Append("_Vertex_sub0_Z"), TString(title).Append(" Vertex_sub0_Z"), 1000, -1, 1, 48, kFALSE),
    vertex_sub1_X(TString(name).Append("_Vertex_sub1_X"), TString(title).Append(" Vertex_sub1_X"), 1000, -1, 1, 48, kFALSE),
    vertex_sub1_Y(TString(name).Append("_Vertex_sub1_Y"), TString(title).Append(" Vertex_sub1_Y"), 1000, -1, 1, 48, kFALSE),
    vertex_sub1_Z(TString(name).Append("_Vertex_sub1_Z"), TString(title).Append(" Vertex_sub1_Z"), 1000, -1, 1, 48, kFALSE),
    vertex_sub2_X(TString(name).Append("_Vertex_sub2_X"), TString(title).Append(" Vertex_sub2_X"), 1000, -1, 1, 48, kFALSE),
    vertex_sub2_Y(TString(name).Append("_Vertex_sub2_Y"), TString(title).Append(" Vertex_sub2_Y"), 1000, -1, 1, 48, kFALSE),
    vertex_sub2_Z(TString(name).Append("_Vertex_sub2_Z"), TString(title).Append(" Vertex_sub2_Z"), 1000, -1, 1, 48, kFALSE),
    pulls(TString(name).Append("_Pulls"), TString(title).Append(" Pulls"), 100, -10, 10, nPulls, 0, nPulls, kFALSE)
{

}

GHistFit::~GHistFit()
{

}

void        GHistFit::CalcResult()
{
    im.CalcResult();
    im_sub0.CalcResult();
    im_sub1.CalcResult();
    im_sub2.CalcResult();
    chiSq.CalcResult();
    confidenceLevel.CalcResult();
    vertex_X.CalcResult();
    vertex_Y.CalcResult();
    vertex_Z.CalcResult();
    vertex_sub0_X.CalcResult();
    vertex_sub0_Y.CalcResult();
    vertex_sub0_Z.CalcResult();
    vertex_sub1_X.CalcResult();
    vertex_sub1_Y.CalcResult();
    vertex_sub1_Z.CalcResult();
    vertex_sub2_X.CalcResult();
    vertex_sub2_Y.CalcResult();
    vertex_sub2_Z.CalcResult();
    pulls.CalcResult();
}

Int_t       GHistFit::Fill(GKinFitter& fitter, const Double_t taggerTime)
{
    /*im.Fill(fitter.GetEtap().M(), taggerTime);
    im_sub0.Fill(fitter.GetEta().M(), taggerTime);
    im_sub1.Fill(fitter.GetPi0a().M(), taggerTime);
    im_sub2.Fill(fitter.GetPi0b().M(), taggerTime);*/
    //chiSq.Fill(fitter.GetChi2(), taggerTime);
    //confidenceLevel.Fill(fitter.ConfidenceLevel(), taggerTime);
    /*vertex_X.Fill(fitter.GetVertex().X(), taggerTime);
    vertex_Y.Fill(fitter.GetVertex().Y(), taggerTime);
    vertex_Z.Fill(fitter.GetVertex().Z(), taggerTime);
    vertex_sub0_X.Fill(fitter.GetVertex_sub0().X(), taggerTime);
    vertex_sub0_Y.Fill(fitter.GetVertex_sub0().Y(), taggerTime);
    vertex_sub0_Z.Fill(fitter.GetVertex_sub0().Z(), taggerTime);
    vertex_sub1_X.Fill(fitter.GetVertex_sub1().X(), taggerTime);
    vertex_sub1_Y.Fill(fitter.GetVertex_sub1().Y(), taggerTime);
    vertex_sub1_Z.Fill(fitter.GetVertex_sub1().Z(), taggerTime);
    vertex_sub2_X.Fill(fitter.GetVertex_sub2().X(), taggerTime);
    vertex_sub2_Y.Fill(fitter.GetVertex_sub2().Y(), taggerTime);
    vertex_sub2_Z.Fill(fitter.GetVertex_sub2().Z(), taggerTime);*/
    //for(int i=0; i<nPulls; i++)
        //pulls.Fill(fitter.GetPull(i), i);
}

Int_t       GHistFit::Fill(GKinFitter& fitter, const Double_t taggerTime, const Int_t taggerChannel)
{
    /*im.Fill(fitter.GetEtap().M(), taggerTime, taggerChannel);
    im_sub0.Fill(fitter.GetEta().M(), taggerTime);
    im_sub1.Fill(fitter.GetPi0a().M(), taggerTime);
    im_sub2.Fill(fitter.GetPi0b().M(), taggerTime);*/
    //chiSq.Fill(fitter.GetChi2(), taggerTime);
    //confidenceLevel.Fill(fitter.ConfidenceLevel(), taggerTime);
    /*vertex_X.Fill(fitter.GetVertex().X(), taggerTime);
    vertex_Y.Fill(fitter.GetVertex().Y(), taggerTime);
    vertex_Z.Fill(fitter.GetVertex().Z(), taggerTime);
    vertex_sub0_X.Fill(fitter.GetVertex_sub0().X(), taggerTime);
    vertex_sub0_Y.Fill(fitter.GetVertex_sub0().Y(), taggerTime);
    vertex_sub0_Z.Fill(fitter.GetVertex_sub0().Z(), taggerTime);
    vertex_sub1_X.Fill(fitter.GetVertex_sub1().X(), taggerTime);
    vertex_sub1_Y.Fill(fitter.GetVertex_sub1().Y(), taggerTime);
    vertex_sub1_Z.Fill(fitter.GetVertex_sub1().Z(), taggerTime);
    vertex_sub2_X.Fill(fitter.GetVertex_sub2().X(), taggerTime);
    vertex_sub2_Y.Fill(fitter.GetVertex_sub2().Y(), taggerTime);
    vertex_sub2_Z.Fill(fitter.GetVertex_sub2().Z(), taggerTime);*/
    //for(int i=0; i<nPulls; i++)
        //pulls.Fill(fitter.GetPull(i), i);
}

void    GHistFit::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        im.PrepareWriteList(arr, TString(name).Append("_IM").Data());
        im_sub0.PrepareWriteList(arr, TString(name).Append("_IM_sub0").Data());
        im_sub1.PrepareWriteList(arr, TString(name).Append("_IM_sub1").Data());
        im_sub2.PrepareWriteList(arr, TString(name).Append("_IM_sub2").Data());
        chiSq.PrepareWriteList(arr, TString(name).Append("_ChiSq").Data());
        confidenceLevel.PrepareWriteList(arr, TString(name).Append("_ConfLev").Data());
        vertex_X.PrepareWriteList(arr, TString(name).Append("_Vertex_X").Data());
        vertex_Y.PrepareWriteList(arr, TString(name).Append("_Vertex_Y").Data());
        vertex_Z.PrepareWriteList(arr, TString(name).Append("_Vertex_Z").Data());
        vertex_sub0_X.PrepareWriteList(arr, TString(name).Append("_Vertex_sub0_X").Data());
        vertex_sub0_Y.PrepareWriteList(arr, TString(name).Append("_Vertex_sub0_Y").Data());
        vertex_sub0_Z.PrepareWriteList(arr, TString(name).Append("_Vertex_sub0_Z").Data());
        vertex_sub1_X.PrepareWriteList(arr, TString(name).Append("_Vertex_sub1_X").Data());
        vertex_sub1_Y.PrepareWriteList(arr, TString(name).Append("_Vertex_sub1_Y").Data());
        vertex_sub1_Z.PrepareWriteList(arr, TString(name).Append("_Vertex_sub1_Z").Data());
        vertex_sub2_X.PrepareWriteList(arr, TString(name).Append("_Vertex_sub2_X").Data());
        vertex_sub2_Y.PrepareWriteList(arr, TString(name).Append("_Vertex_sub2_Y").Data());
        vertex_sub2_Z.PrepareWriteList(arr, TString(name).Append("_Vertex_sub2_Z").Data());
        pulls.PrepareWriteList(arr, TString(name).Append("_Pulls").Data());
    }
    else
    {
        im.PrepareWriteList(arr);
        im_sub0.PrepareWriteList(arr);
        im_sub1.PrepareWriteList(arr);
        im_sub2.PrepareWriteList(arr);
        chiSq.PrepareWriteList(arr);
        confidenceLevel.PrepareWriteList(arr);
        vertex_X.PrepareWriteList(arr);
        vertex_Y.PrepareWriteList(arr);
        vertex_Z.PrepareWriteList(arr);
        vertex_sub0_X.PrepareWriteList(arr);
        vertex_sub0_Y.PrepareWriteList(arr);
        vertex_sub0_Z.PrepareWriteList(arr);
        vertex_sub1_X.PrepareWriteList(arr);
        vertex_sub1_Y.PrepareWriteList(arr);
        vertex_sub1_Z.PrepareWriteList(arr);
        vertex_sub2_X.PrepareWriteList(arr);
        vertex_sub2_Y.PrepareWriteList(arr);
        vertex_sub2_Z.PrepareWriteList(arr);
        pulls.PrepareWriteList(arr);
    }
}

void        GHistFit::Reset(Option_t* option)
{
    im.Reset(option);
    im_sub0.Reset(option);
    im_sub1.Reset(option);
    im_sub2.Reset(option);
    chiSq.Reset(option);
    confidenceLevel.Reset(option);
    vertex_X.Reset(option);
    vertex_Y.Reset(option);
    vertex_Z.Reset(option);
    vertex_sub0_X.Reset(option);
    vertex_sub0_Y.Reset(option);
    vertex_sub0_Z.Reset(option);
    vertex_sub1_X.Reset(option);
    vertex_sub1_Y.Reset(option);
    vertex_sub1_Z.Reset(option);
    vertex_sub2_X.Reset(option);
    vertex_sub2_Y.Reset(option);
    vertex_sub2_Z.Reset(option);
    pulls.Reset(option);
}

void        GHistFit::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im_sub0.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im_sub1.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im_sub2.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    chiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    confidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_X.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_Y.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_Z.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub0_X.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub0_Y.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub0_Z.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub1_X.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub1_Y.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub1_Z.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub2_X.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub2_Y.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub2_Z.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pulls.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}


















GHistFit2::GHistFit2(const char* name, const char* title, const Int_t _NPulls, Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    count(0),
    nPulls(_NPulls),
    im(TString(name).Append("_im"), TString(title).Append(" inv. Mass"), 2000, 0, 2000, 5, 0, 5, kFALSE),
    im_sub0(TString(name).Append("_im_sub0"), TString(title).Append(" inv. Mass SubPart0"), 2000, 0, 2000, 5, 0, 5, kFALSE),
    im_sub1(TString(name).Append("_im_sub1"), TString(title).Append(" inv. Mass SubPart1"), 2000, 0, 2000, 5, 0, 5, kFALSE),
    im_sub2(TString(name).Append("_im_sub2"), TString(title).Append(" inv. Mass SubPart2"), 2000, 0, 2000, 5, 0, 5, kFALSE),
    chiSq(TString(name).Append("_ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 100, 5, 0, 5, kFALSE),
    confidenceLevel(TString(name).Append("_ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, 5, 0, 5, kFALSE),
    vertex_X(TString(name).Append("_Vertex_X"), TString(title).Append(" Vertex_X"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_Y(TString(name).Append("_Vertex_Y"), TString(title).Append(" Vertex_Y"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_Z(TString(name).Append("_Vertex_Z"), TString(title).Append(" Vertex_Z"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_sub0_X(TString(name).Append("_Vertex_sub0_X"), TString(title).Append(" Vertex_sub0_X"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_sub0_Y(TString(name).Append("_Vertex_sub0_Y"), TString(title).Append(" Vertex_sub0_Y"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_sub0_Z(TString(name).Append("_Vertex_sub0_Z"), TString(title).Append(" Vertex_sub0_Z"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_sub1_X(TString(name).Append("_Vertex_sub1_X"), TString(title).Append(" Vertex_sub1_X"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_sub1_Y(TString(name).Append("_Vertex_sub1_Y"), TString(title).Append(" Vertex_sub1_Y"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_sub1_Z(TString(name).Append("_Vertex_sub1_Z"), TString(title).Append(" Vertex_sub1_Z"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_sub2_X(TString(name).Append("_Vertex_sub2_X"), TString(title).Append(" Vertex_sub2_X"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_sub2_Y(TString(name).Append("_Vertex_sub2_Y"), TString(title).Append(" Vertex_sub2_Y"), 1000, -1, 1, 5, 0, 5, kFALSE),
    vertex_sub2_Z(TString(name).Append("_Vertex_sub2_Z"), TString(title).Append(" Vertex_sub2_Z"), 1000, -1, 1, 5, 0, 5, kFALSE),
    pulls0(TString(name).Append("_Pulls0"), TString(title).Append(" Pulls0"), 100, -10, 10, nPulls, 0, nPulls, kFALSE),
    pulls1(TString(name).Append("_Pulls1"), TString(title).Append(" Pulls1"), 100, -10, 10, nPulls, 0, nPulls, kFALSE),
    pulls2(TString(name).Append("_Pulls2"), TString(title).Append(" Pulls2"), 100, -10, 10, nPulls, 0, nPulls, kFALSE),
    pulls3(TString(name).Append("_Pulls3"), TString(title).Append(" Pulls3"), 100, -10, 10, nPulls, 0, nPulls, kFALSE),
    pulls4(TString(name).Append("_Pulls4"), TString(title).Append(" Pulls4"), 100, -10, 10, nPulls, 0, nPulls, kFALSE)
{

}

GHistFit2::~GHistFit2()
{

}

void        GHistFit2::CalcResult()
{
    im.CalcResult();
    im_sub0.CalcResult();
    im_sub1.CalcResult();
    im_sub2.CalcResult();
    chiSq.CalcResult();
    confidenceLevel.CalcResult();
    vertex_X.CalcResult();
    vertex_Y.CalcResult();
    vertex_Z.CalcResult();
    vertex_sub0_X.CalcResult();
    vertex_sub0_Y.CalcResult();
    vertex_sub0_Z.CalcResult();
    vertex_sub1_X.CalcResult();
    vertex_sub1_Y.CalcResult();
    vertex_sub1_Z.CalcResult();
    vertex_sub2_X.CalcResult();
    vertex_sub2_Y.CalcResult();
    vertex_sub2_Z.CalcResult();
    pulls0.CalcResult();
    pulls1.CalcResult();
    pulls2.CalcResult();
    pulls3.CalcResult();
    pulls4.CalcResult();
}

Int_t       GHistFit2::Fill(GKinFitter& fitter, const Double_t taggerTime)
{
    /*im.Fill(fitter.GetEtap().M(), count, taggerTime);
    im_sub0.Fill(fitter.GetEta().M(), count, taggerTime);
    im_sub1.Fill(fitter.GetPi0a().M(), count, taggerTime);
    im_sub2.Fill(fitter.GetPi0b().M(), count, taggerTime);*/
    //chiSq.Fill(fitter.GetChi2(), count, taggerTime);
    //confidenceLevel.Fill(fitter.ConfidenceLevel(), count, taggerTime);
    /*vertex_X.Fill(fitter.GetVertex().X(), count, taggerTime);
    vertex_Y.Fill(fitter.GetVertex().Y(), count, taggerTime);
    vertex_Z.Fill(fitter.GetVertex().Z(), count, taggerTime);
    vertex_sub0_X.Fill(fitter.GetVertex_sub0().X(), count, taggerTime);
    vertex_sub0_Y.Fill(fitter.GetVertex_sub0().Y(), count, taggerTime);
    vertex_sub0_Z.Fill(fitter.GetVertex_sub0().Z(), count, taggerTime);
    vertex_sub1_X.Fill(fitter.GetVertex_sub1().X(), count, taggerTime);
    vertex_sub1_Y.Fill(fitter.GetVertex_sub1().Y(), count, taggerTime);
    vertex_sub1_Z.Fill(fitter.GetVertex_sub1().Z(), count, taggerTime);
    vertex_sub2_X.Fill(fitter.GetVertex_sub2().X(), count, taggerTime);
    vertex_sub2_Y.Fill(fitter.GetVertex_sub2().Y(), count, taggerTime);
    vertex_sub2_Z.Fill(fitter.GetVertex_sub2().Z(), count, taggerTime);*/
    /*switch(count)
    {
    case 0:
        for(int i=0; i<nPulls; i++)
            pulls0.Fill(fitter.GetPull(i), i);
        break;
    case 1:
        for(int i=0; i<nPulls; i++)
            pulls1.Fill(fitter.GetPull(i), i);
        break;
    case 2:
        for(int i=0; i<nPulls; i++)
            pulls2.Fill(fitter.GetPull(i), i);
        break;
    case 3:
        for(int i=0; i<nPulls; i++)
            pulls3.Fill(fitter.GetPull(i), i);
        break;
    case 4:
        for(int i=0; i<nPulls; i++)
            pulls4.Fill(fitter.GetPull(i), i);
        break;
    }*/
    count++;
}

Int_t       GHistFit2::Fill(GKinFitter& fitter, const Double_t taggerTime, const Int_t taggerChannel)
{
    /*im.Fill(fitter.GetEtap().M(), count, taggerTime);
    im_sub0.Fill(fitter.GetEta().M(), count, taggerTime);
    im_sub1.Fill(fitter.GetPi0a().M(), count, taggerTime);
    im_sub2.Fill(fitter.GetPi0b().M(), count, taggerTime);*/
    //chiSq.Fill(fitter.GetChi2(), count, taggerTime);
    //confidenceLevel.Fill(fitter.ConfidenceLevel(), count, taggerTime);
    /*vertex_X.Fill(fitter.GetVertex().X(), count, taggerTime);
    vertex_Y.Fill(fitter.GetVertex().Y(), count, taggerTime);
    vertex_Z.Fill(fitter.GetVertex().Z(), count, taggerTime);
    vertex_sub0_X.Fill(fitter.GetVertex_sub0().X(), count, taggerTime);
    vertex_sub0_Y.Fill(fitter.GetVertex_sub0().Y(), count, taggerTime);
    vertex_sub0_Z.Fill(fitter.GetVertex_sub0().Z(), count, taggerTime);
    vertex_sub1_X.Fill(fitter.GetVertex_sub1().X(), count, taggerTime);
    vertex_sub1_Y.Fill(fitter.GetVertex_sub1().Y(), count, taggerTime);
    vertex_sub1_Z.Fill(fitter.GetVertex_sub1().Z(), count, taggerTime);
    vertex_sub2_X.Fill(fitter.GetVertex_sub2().X(), count, taggerTime);
    vertex_sub2_Y.Fill(fitter.GetVertex_sub2().Y(), count, taggerTime);
    vertex_sub2_Z.Fill(fitter.GetVertex_sub2().Z(), count, taggerTime);*/
    /*switch(count)
    {
    case 0:
        for(int i=0; i<nPulls; i++)
            pulls0.Fill(fitter.GetPull(i), i);
        break;
    case 1:
        for(int i=0; i<nPulls; i++)
            pulls1.Fill(fitter.GetPull(i), i);
        break;
    case 2:
        for(int i=0; i<nPulls; i++)
            pulls2.Fill(fitter.GetPull(i), i);
        break;
    case 3:
        for(int i=0; i<nPulls; i++)
            pulls3.Fill(fitter.GetPull(i), i);
        break;
    case 4:
        for(int i=0; i<nPulls; i++)
            pulls4.Fill(fitter.GetPull(i), i);
        break;
    }*/
    count++;
}

void    GHistFit2::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        im.PrepareWriteList(arr, TString(name).Append("_IM").Data());
        im_sub0.PrepareWriteList(arr, TString(name).Append("_IM_sub0").Data());
        im_sub1.PrepareWriteList(arr, TString(name).Append("_IM_sub1").Data());
        im_sub2.PrepareWriteList(arr, TString(name).Append("_IM_sub2").Data());
        chiSq.PrepareWriteList(arr, TString(name).Append("_ChiSq").Data());
        confidenceLevel.PrepareWriteList(arr, TString(name).Append("_ConfLev").Data());
        vertex_X.PrepareWriteList(arr, TString(name).Append("_Vertex_X").Data());
        vertex_Y.PrepareWriteList(arr, TString(name).Append("_Vertex_Y").Data());
        vertex_Z.PrepareWriteList(arr, TString(name).Append("_Vertex_Z").Data());
        vertex_sub0_X.PrepareWriteList(arr, TString(name).Append("_Vertex_sub0_X").Data());
        vertex_sub0_Y.PrepareWriteList(arr, TString(name).Append("_Vertex_sub0_Y").Data());
        vertex_sub0_Z.PrepareWriteList(arr, TString(name).Append("_Vertex_sub0_Z").Data());
        vertex_sub1_X.PrepareWriteList(arr, TString(name).Append("_Vertex_sub1_X").Data());
        vertex_sub1_Y.PrepareWriteList(arr, TString(name).Append("_Vertex_sub1_Y").Data());
        vertex_sub1_Z.PrepareWriteList(arr, TString(name).Append("_Vertex_sub1_Z").Data());
        vertex_sub2_X.PrepareWriteList(arr, TString(name).Append("_Vertex_sub2_X").Data());
        vertex_sub2_Y.PrepareWriteList(arr, TString(name).Append("_Vertex_sub2_Y").Data());
        vertex_sub2_Z.PrepareWriteList(arr, TString(name).Append("_Vertex_sub2_Z").Data());
        pulls0.PrepareWriteList(arr, TString(name).Append("_Pulls0").Data());
        pulls1.PrepareWriteList(arr, TString(name).Append("_Pulls1").Data());
        pulls2.PrepareWriteList(arr, TString(name).Append("_Pulls2").Data());
        pulls3.PrepareWriteList(arr, TString(name).Append("_Pulls3").Data());
        pulls4.PrepareWriteList(arr, TString(name).Append("_Pulls4").Data());
    }
    else
    {
        im.PrepareWriteList(arr);
        im_sub0.PrepareWriteList(arr);
        im_sub1.PrepareWriteList(arr);
        im_sub2.PrepareWriteList(arr);
        chiSq.PrepareWriteList(arr);
        confidenceLevel.PrepareWriteList(arr);
        vertex_X.PrepareWriteList(arr);
        vertex_Y.PrepareWriteList(arr);
        vertex_Z.PrepareWriteList(arr);
        vertex_sub0_X.PrepareWriteList(arr);
        vertex_sub0_Y.PrepareWriteList(arr);
        vertex_sub0_Z.PrepareWriteList(arr);
        vertex_sub1_X.PrepareWriteList(arr);
        vertex_sub1_Y.PrepareWriteList(arr);
        vertex_sub1_Z.PrepareWriteList(arr);
        vertex_sub2_X.PrepareWriteList(arr);
        vertex_sub2_Y.PrepareWriteList(arr);
        vertex_sub2_Z.PrepareWriteList(arr);
        pulls0.PrepareWriteList(arr);
        pulls1.PrepareWriteList(arr);
        pulls2.PrepareWriteList(arr);
        pulls3.PrepareWriteList(arr);
        pulls4.PrepareWriteList(arr);
    }
}

void        GHistFit2::Reset(Option_t* option)
{
    count   = 0;
    im.Reset(option);
    im_sub0.Reset(option);
    im_sub1.Reset(option);
    im_sub2.Reset(option);
    chiSq.Reset(option);
    confidenceLevel.Reset(option);
    vertex_X.Reset(option);
    vertex_Y.Reset(option);
    vertex_Z.Reset(option);
    vertex_sub0_X.Reset(option);
    vertex_sub0_Y.Reset(option);
    vertex_sub0_Z.Reset(option);
    vertex_sub1_X.Reset(option);
    vertex_sub1_Y.Reset(option);
    vertex_sub1_Z.Reset(option);
    vertex_sub2_X.Reset(option);
    vertex_sub2_Y.Reset(option);
    vertex_sub2_Z.Reset(option);
    pulls0.Reset(option);
    pulls1.Reset(option);
    pulls2.Reset(option);
    pulls3.Reset(option);
    pulls4.Reset(option);
}

void        GHistFit2::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    im.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im_sub0.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im_sub1.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    im_sub2.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    chiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    confidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_X.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_Y.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_Z.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub0_X.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub0_Y.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub0_Z.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub1_X.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub1_Y.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub1_Z.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub2_X.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub2_Y.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    vertex_sub2_Z.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pulls0.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pulls1.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pulls2.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pulls3.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pulls4.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}
