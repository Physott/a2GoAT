#include "GHistEvent.h"
#include "GTreeTagger.h"
#include "GTreeMeson.h"




GHistEvent::GHistEvent(const char* name, const char* title, Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    im(TString(name).Append("_im"), TString(title).Append(" inv. Mass"), 1500, 0, 1500, 48, kFALSE),
    mm(TString(name).Append("_mm"), TString(title).Append(" mis. Mass"), 2000, 0, 2000, 48, kFALSE),
    EnergyTheta(TString(name).Append("_EnergyTheta"), TString(title).Append(" EnergyTheta"), 800, 0, 800, 180, 0, 180, kFALSE),
    phiTheta(TString(name).Append("_phiTheta"), TString(title).Append(" PhiTheta"), 360, -180, 180, 180, 0, 180, kFALSE)
{

}

GHistEvent::~GHistEvent()
{

}


void    GHistEvent::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    if(name)
    {
        im.PrepareWriteList(arr, TString(name).Append("_IM").Data());
        mm.PrepareWriteList(arr, TString(name).Append("_MM").Data());
        EnergyTheta.PrepareWriteList(arr, TString(name).Append("_EnergyTheta").Data());
        phiTheta.PrepareWriteList(arr, TString(name).Append("_phiTheta").Data());
    }
    else
    {
        im.PrepareWriteList(arr);
        mm.PrepareWriteList(arr);
        EnergyTheta.PrepareWriteList(arr);
        phiTheta.PrepareWriteList(arr);
    }
}













GHistEvent3Mesons::GHistEvent3Mesons(const char* name, const char* title, Bool_t linkHistogram) :
    GHistEvent(name, title, linkHistogram),
    sub0_im(TString(name).Append("_sub0im"), TString(title).Append(" sub Part. 0 inv. Mass"), 800, 0, 800, 48, kFALSE),
    sub1_im(TString(name).Append("_sub1im"), TString(title).Append(" sub Part. 1 inv. Mass"), 400, 0, 400, 48, kFALSE),
    sub2_im(TString(name).Append("_sub2im"), TString(title).Append(" sub Part. 2 inv. Mass"), 400, 0, 400, 48, kFALSE),
    subEtaPi0_im(TString(name).Append("_subEtaPi0im"), TString(title).Append(" sub Part. EtaPi0 IM"), 800, 0, 800, 400, 0, 400, kFALSE),
    subPi0Pi0_im(TString(name).Append("_subPi0Pi0im"), TString(title).Append(" sub Part. Pi0Pi0 IM"), 400, 0, 400, 400, 0, 400, kFALSE)
{

}

GHistEvent3Mesons::~GHistEvent3Mesons()
{

}


void    GHistEvent3Mesons::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistEvent::PrepareWriteList(arr, name);

    if(name)
    {
        sub0_im.PrepareWriteList(arr, TString(name).Append("_sub0IM").Data());
        sub1_im.PrepareWriteList(arr, TString(name).Append("_sub1IM").Data());
        sub2_im.PrepareWriteList(arr, TString(name).Append("_sub2IM").Data());
        subEtaPi0_im.PrepareWriteList(arr, TString(name).Append("_subEtaPi0im").Data());
        subPi0Pi0_im.PrepareWriteList(arr, TString(name).Append("_subPi0Pi0im").Data());
    }
    else
    {
        sub0_im.PrepareWriteList(arr);
        sub1_im.PrepareWriteList(arr);
        sub2_im.PrepareWriteList(arr);
        subEtaPi0_im.PrepareWriteList(arr);
        subPi0Pi0_im.PrepareWriteList(arr);
    }
}
















GHistEvent3MesonsProton::GHistEvent3MesonsProton(const char* name, const char* title, Bool_t linkHistogram) :
    GHistEvent3Mesons(name, title, linkHistogram),
    protonEnergy(TString(name).Append("_protonEnergy"), TString(title).Append(" proton energy"), 800, 0, 800, 48, kFALSE),
    protontheta(TString(name).Append("_protontheta"), TString(title).Append(" proton theta"), 180, 0, 180, 48, kFALSE),
    protonEnergyTheta(TString(name).Append("_protonEnergyTheta"), TString(title).Append(" proton EnergyTheta"), 800, 0, 800, 180, 0, 180, kFALSE),
    protonPhiTheta(TString(name).Append("_protonPhiTheta"), TString(title).Append(" proton PhiTheta"), 360, -180, 180, 180, 0, 180, kFALSE)
{

}

GHistEvent3MesonsProton::~GHistEvent3MesonsProton()
{

}

void    GHistEvent3MesonsProton::PrepareWriteList(GHistWriteList* arr, const char* name)
{
    if(!arr)
        return;

    GHistEvent3Mesons::PrepareWriteList(arr, name);

    if(name)
    {
        protonEnergy.PrepareWriteList(arr, TString(name).Append("_protonEnergy").Data());
        protontheta.PrepareWriteList(arr, TString(name).Append("_protontheta").Data());
        protonEnergyTheta.PrepareWriteList(arr, TString(name).Append("_protonEnergyTheta").Data());
        protonPhiTheta.PrepareWriteList(arr, TString(name).Append("_protonPhiTheta").Data());
    }
    else
    {
        protonEnergy.PrepareWriteList(arr);
        protontheta.PrepareWriteList(arr);
        protonEnergyTheta.PrepareWriteList(arr);
        protonPhiTheta.PrepareWriteList(arr);
    }
}
