#ifndef __PHistEvent_h__
#define __PHistEvent_h__

#include "PHist.h"


class   PHistEvent
{
protected:
    TString name;

    PHistD	IM;
    PHistD	MM;

public:
    PHistEvent(const TString& _Name);
    virtual ~PHistEvent();

    inline  void    Clear();
    inline  void    Fill(const Double_t invMass, const Double_t misMass, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void    Write(TDirectory& dir);
};



void    PHistEvent::Clear()
{
    IM.Clear();
    MM.Clear();
}
void    PHistEvent::Fill(const Double_t invMass, const Double_t misMass, const Double_t taggerTime, const Int_t taggerChannel)
{
    IM.Fill(invMass, taggerTime, taggerChannel);
    MM.Fill(misMass, taggerTime, taggerChannel);
}




class   PHistEvent3Meson    : public PHistEvent
{
private:
    PHistD	sub0;
    PHistD	sub1;
    PHistD	sub2;

public:
    PHistEvent3Meson(const TString& _Name);
    virtual ~PHistEvent3Meson();

    inline  void    Clear();
    inline  void    FillSubMesons(const Double_t invMassSub0, const Double_t invMassSub1, const Double_t invMassSub2, const Double_t taggerTime, const Int_t taggerChannel);
    virtual void    Write(TDirectory& dir);
};

void    PHistEvent3Meson::Clear()
{
    PHistEvent::Clear();
    sub0.Clear();
    sub1.Clear();
    sub2.Clear();
}

void    PHistEvent3Meson::FillSubMesons(const Double_t invMassSub0, const Double_t invMassSub1, const Double_t invMassSub2, const Double_t taggerTime, const Int_t taggerChannel)
{
    sub0.Fill(invMassSub0, taggerTime, taggerChannel);
    sub1.Fill(invMassSub1, taggerTime, taggerChannel);
    sub2.Fill(invMassSub2, taggerTime, taggerChannel);
}

#endif
