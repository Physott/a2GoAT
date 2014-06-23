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

    inline  void    Fill(const Double_t taggerTime, const Double_t invMass, const Double_t misMass);
    virtual void    SetCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t Rand0Min = 0, const Double_t Rand0Max = 0, const Double_t Rand1Min = 0, const Double_t Rand1Max = 0);
    virtual void	RandomSubtraction();
    virtual void    Write(TDirectory *dir);
};





void    PHistEvent::Fill(const Double_t taggerTime, const Double_t invMass, const Double_t misMass)
{
    IM.Fill(taggerTime, invMass);
    MM.Fill(taggerTime, misMass);
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

    inline  void    FillSubMesons(const Double_t taggerTime, const Double_t invMassSub0, const Double_t invMassSub1, const Double_t invMassSub2);
    virtual void    SetCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t Rand0Min = 0, const Double_t Rand0Max = 0, const Double_t Rand1Min = 0, const Double_t Rand1Max = 0);
    virtual void	RandomSubtraction();
    virtual void    Write(TDirectory *dir);
};

void    PHistEvent3Meson::FillSubMesons(const Double_t taggerTime, const Double_t invMassSub0, const Double_t invMassSub1, const Double_t invMassSub2)
{
    sub0.Fill(taggerTime, invMassSub0);
    sub1.Fill(taggerTime, invMassSub1);
    sub2.Fill(taggerTime, invMassSub2);
}

#endif
