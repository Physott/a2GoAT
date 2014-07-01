#include "PHistEvent.h"


PHistEvent::PHistEvent(const TString& _Name)    :
    name(_Name),
    IM(TString("_IM").Prepend(name), TString("IM ").Append(name), 1500, 0, 1500),
    MM(TString("_MM").Prepend(name), TString("MM ").Append(name), 1500, 0, 1500)
{

}

PHistEvent::~PHistEvent()
{

}

void    PHistEvent::Write(TDirectory& dir)
{
    TDirectory* curDir  = dir.GetDirectory(name.Data());
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir(name.Data());
        curDir  = dir.GetDirectory(name.Data());
    }
    IM.Write(*curDir);
    MM.Write(*curDir);
}






PHistEvent3Meson::PHistEvent3Meson(const TString& _Name)    :
    PHistEvent(_Name),
    sub0(TString("_sub0").Prepend(name), TString("IM sub0 ").Append(name), 1500, 0, 1500),
    sub1(TString("_sub1").Prepend(name), TString("IM sub1 ").Append(name), 1500, 0, 1500),
    sub2(TString("_sub2").Prepend(name), TString("IM sub2 ").Append(name), 1500, 0, 1500)
{

}

PHistEvent3Meson::~PHistEvent3Meson()
{

}

void    PHistEvent3Meson::Write(TDirectory& dir)
{
    TDirectory* curDir  = dir.GetDirectory(name.Data());
    if(!curDir)
    {
        dir.cd();
        gDirectory->mkdir(name.Data());
        curDir  = dir.GetDirectory(name.Data());
    }
    IM.Write(*curDir);
    MM.Write(*curDir);
    sub0.Write(*curDir);
    sub1.Write(*curDir);
    sub2.Write(*curDir);
}
