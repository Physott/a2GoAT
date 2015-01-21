#include "GHistFit.h"
#include "GTreeMeson.h"




GHistFit::GHistFit(const char* name, const char* title, const Int_t nPulls, Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    invMass(TString(name).Append("_im"), TString(title).Append(" inv. Mass"), 2000, 0, 2000, 48, kFALSE),
    invMassSub0(TString(name).Append("_invMassSub0"), TString(title).Append(" inv. Mass SubPart0"), 2000, 0, 2000, 48, kFALSE),
    invMassSub1(TString(name).Append("_invMassSub1"), TString(title).Append(" inv. Mass SubPart1"), 2000, 0, 2000, 48, kFALSE),
    invMassSub2(TString(name).Append("_invMassSub2"), TString(title).Append(" inv. Mass SubPart2"), 2000, 0, 2000, 48, kFALSE),
    chiSq(TString(name).Append("_ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 100, 48, kFALSE),
    confidenceLevel(TString(name).Append("_ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, 48, kFALSE),
    zVertex(TString(name).Append("_zVertex"), TString(title).Append(" zVertex"), 1000, -1, 1, 48, kFALSE),
    pulls(TString(name).Append("_Pulls"), TString(title).Append(" Pulls"), 100, -10, 10, nPulls, 0, nPulls, kFALSE)
{

}

GHistFit::~GHistFit()
{

}

void        GHistFit::CalcResult()
{
    invMass.CalcResult();
    invMassSub0.CalcResult();
    invMassSub1.CalcResult();
    invMassSub2.CalcResult();
    chiSq.CalcResult();
    confidenceLevel.CalcResult();
    zVertex.CalcResult();
    pulls.CalcResult();
}

void    GHistFit::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        invMass.PrepareWriteList(arr, TString(name).Append("_IM").Data());
        invMassSub0.PrepareWriteList(arr, TString(name).Append("_IM_sub0").Data());
        invMassSub1.PrepareWriteList(arr, TString(name).Append("_IM_sub1").Data());
        invMassSub2.PrepareWriteList(arr, TString(name).Append("_IM_sub2").Data());
        chiSq.PrepareWriteList(arr, TString(name).Append("_ChiSq").Data());
        confidenceLevel.PrepareWriteList(arr, TString(name).Append("_ConfLev").Data());
        zVertex.PrepareWriteList(arr, TString(name).Append("_Vertex_Z").Data());
    }
    else
    {
        invMass.PrepareWriteList(arr);
        invMassSub0.PrepareWriteList(arr);
        invMassSub1.PrepareWriteList(arr);
        invMassSub2.PrepareWriteList(arr);
        chiSq.PrepareWriteList(arr);
        confidenceLevel.PrepareWriteList(arr);
        zVertex.PrepareWriteList(arr);
        pulls.PrepareWriteList(arr);
    }
}

void        GHistFit::Reset(Option_t* option)
{
    invMass.Reset(option);
    invMassSub0.Reset(option);
    invMassSub1.Reset(option);
    invMassSub2.Reset(option);
    chiSq.Reset(option);
    confidenceLevel.Reset(option);
    zVertex.Reset(option);
    pulls.Reset(option);
}

void        GHistFit::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    invMass.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    invMassSub0.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    invMassSub1.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    invMassSub2.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    chiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    confidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    zVertex.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    pulls.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}


















GHistFit2::GHistFit2(const char* name, const char* title, Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    invMass(TString(name).Append("_im"), TString(title).Append(" inv. Mass"), 2000, 0, 2000, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    invMassSub0(TString(name).Append("_invMassSub0"), TString(title).Append(" inv. Mass SubPart0"), 2000, 0, 2000, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    invMassSub1(TString(name).Append("_invMassSub1"), TString(title).Append(" inv. Mass SubPart1"), 2000, 0, 2000, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    invMassSub2(TString(name).Append("_invMassSub2"), TString(title).Append(" inv. Mass SubPart2"), 2000, 0, 2000, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    chiSq(TString(name).Append("_ChiSq"), TString(title).Append(" ChiSq"), 1000, 0, 100, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    confidenceLevel(TString(name).Append("_ConfLev"), TString(title).Append(" ConfLev"), 1000, 0, 1, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    zVertex(TString(name).Append("_Vertex_X"), TString(title).Append(" Vertex_X"), 1000, -1, 1, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    protonEnergy(TString(name).Append("protonEnergy"), TString(title).Append(" protonEnergy"), 1000, -1, 1, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    protonTheta(TString(name).Append("protonTheta"), TString(title).Append(" protonTheta"), 1000, -1, 1, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    protonPhi(TString(name).Append("protonPhi"), TString(title).Append(" protonPhi"), 1000, -1, 1, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    beamEnergy(TString(name).Append("protonEnergy"), TString(title).Append(" protonEnergy"), 1000, -1, 1, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    beamTheta(TString(name).Append("protonTheta"), TString(title).Append(" protonTheta"), 1000, -1, 1, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE),
    beamPhi(TString(name).Append("protonPhi"), TString(title).Append(" protonPhi"), 1000, -1, 1, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE)
{
    for(int i=0; i<6; i++)
    {
        photonsEnergy[i]    = new GHistBGSub2(TString(name).Append("_PhotonE").Append(TString().Itoa(i,10)), TString(title).Append(" PhotonE").Append(TString().Itoa(i,10)), 1000, -1, 1, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE);
        photonsTheta[i]     = new GHistBGSub2(TString(name).Append("_PhotonTheta").Append(TString().Itoa(i,10)), TString(title).Append(" PhotonTheta").Append(TString().Itoa(i,10)), 360, 0, 180, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE);
        photonsPhi[i]       = new GHistBGSub2(TString(name).Append("_PhotonPhi").Append(TString().Itoa(i,10)), TString(title).Append(" PhotonPhi").Append(TString().Itoa(i,10)), 540, -180, 360, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE);
    }
    for(int i=0; i<7; i++)
        con[i]    = new GHistBGSub2(TString(name).Append("_Con").Append(TString().Itoa(i,10)), TString(title).Append(" Constraint ").Append(TString().Itoa(i,10)), 1000, -100000, 100000, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE);
    for(int i=0; i<23; i++)
        pulls[i]    = new GHistBGSub2(TString(name).Append("_Pull").Append(TString().Itoa(i,10)), TString(title).Append(" Pull ").Append(TString().Itoa(i,10)), 1000, -10, 10, GKinFitter_MaxSteps, 0, GKinFitter_MaxSteps, kFALSE);
}

GHistFit2::~GHistFit2()
{

}

void        GHistFit2::CalcResult()
{
    invMass.CalcResult();
    invMassSub0.CalcResult();
    invMassSub1.CalcResult();
    invMassSub2.CalcResult();
    chiSq.CalcResult();
    confidenceLevel.CalcResult();
    zVertex.CalcResult();
    for(int i=0; i<6; i++)
    {
        photonsEnergy[i]->CalcResult();
        photonsTheta[i]->CalcResult();
        photonsPhi[i]->CalcResult();
    }
    for(int i=0; i<7; i++)
        con[i]->CalcResult();
    for(int i=0; i<23; i++)
        pulls[i]->CalcResult();

}

void    GHistFit2::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        invMass.PrepareWriteList(arr, TString(name).Append("_IM").Data());
        invMassSub0.PrepareWriteList(arr, TString(name).Append("_IM_sub0").Data());
        invMassSub1.PrepareWriteList(arr, TString(name).Append("_IM_sub1").Data());
        invMassSub2.PrepareWriteList(arr, TString(name).Append("_IM_sub2").Data());
        chiSq.PrepareWriteList(arr, TString(name).Append("_ChiSq").Data());
        confidenceLevel.PrepareWriteList(arr, TString(name).Append("_ConfLev").Data());
        zVertex.PrepareWriteList(arr, TString(name).Append("_Vertex_X").Data());
        for(int i=0; i<6; i++)
        {
            photonsEnergy[i]->PrepareWriteList(arr, TString(name).Append("_photonsEnergy").Append(TString().Itoa(i,10)).Data());
            photonsTheta[i]->PrepareWriteList(arr, TString(name).Append("_photonsTheta").Append(TString().Itoa(i,10)).Data());
            photonsPhi[i]->PrepareWriteList(arr, TString(name).Append("_photonsPhi").Append(TString().Itoa(i,10)).Data());
        }
        for(int i=0; i<7; i++)
            con[i]->PrepareWriteList(arr, TString(name).Append("_Con").Append(TString().Itoa(i,10)).Data());
        for(int i=0; i<23; i++)
            pulls[i]->PrepareWriteList(arr, TString(name).Append("_Pull").Append(TString().Itoa(i,10)).Data());
    }
    else
    {
        invMass.PrepareWriteList(arr);
        invMassSub0.PrepareWriteList(arr);
        invMassSub1.PrepareWriteList(arr);
        invMassSub2.PrepareWriteList(arr);
        chiSq.PrepareWriteList(arr);
        confidenceLevel.PrepareWriteList(arr);
        zVertex.PrepareWriteList(arr);
        for(int i=0; i<6; i++)
        {
            photonsEnergy[i]->PrepareWriteList(arr);
            photonsTheta[i]->PrepareWriteList(arr);
            photonsPhi[i]->PrepareWriteList(arr);
        }
        for(int i=0; i<7; i++)
            con[i]->PrepareWriteList(arr);
        for(int i=0; i<23; i++)
            pulls[i]->PrepareWriteList(arr);
    }
}

void        GHistFit2::Reset(Option_t* option)
{
    invMass.Reset(option);
    invMassSub0.Reset(option);
    invMassSub1.Reset(option);
    invMassSub2.Reset(option);
    chiSq.Reset(option);
    confidenceLevel.Reset(option);
    zVertex.Reset(option);
    for(int i=0; i<6; i++)
    {
        photonsEnergy[i]->Reset(option);
        photonsTheta[i]->Reset(option);
        photonsPhi[i]->Reset(option);
    }
    for(int i=0; i<7; i++)
        con[i]->Reset(option);
    for(int i=0; i<23; i++)
        pulls[i]->Reset(option);
}

void        GHistFit2::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    invMass.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    invMassSub0.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    invMassSub1.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    invMassSub2.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    chiSq.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    confidenceLevel.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    zVertex.ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    for(int i=0; i<6; i++)
    {
        photonsEnergy[i]->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
        photonsTheta[i]->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
        photonsPhi[i]->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    }
    for(int i=0; i<7; i++)
        con[i]->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
    for(int i=0; i<23; i++)
        pulls[i]->ScalerReadCorrection(CorrectionFactor, CreateHistogramsForSingleScalerReads);
}
