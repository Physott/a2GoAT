#ifndef __PHistEvent_h__
#define __PHistEvent_h__

#include "PHist.h"

class  PHistEvent
{
private:
    PHistD	IM;
    PHistD	MM;

protected:
    TString name;

public:
    PHistEvent(const TString& _Name);
    virtual ~PHistEvent();

    inline  void    Fill(const Double_t taggerTime, const Double_t invMass, const Double_t misMass);
            void    SetCuts(const Double_t PromptMin, const Double_t PromptMax, const Double_t Rand0Min = 0, const Double_t Rand0Max = 0, const Double_t Rand1Min = 0, const Double_t Rand1Max = 0);
            void	RandomSubtraction();
            void    Write(TDirectory *dir);
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
};


#endif
