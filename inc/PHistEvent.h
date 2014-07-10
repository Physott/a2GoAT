#ifndef __PHistEvent_h__
#define __PHistEvent_h__

#include <TH1D.h>
#include <TLorentzVector.h>

#include "PHist.h"

class   GTreeTagger;

class   PHistEvent
{
protected:
    TString name;

    TH1D	time;
    PHistD	IM;
    PHistD	MM;

public:
    PHistEvent(const TString& _Name);
    virtual ~PHistEvent();

    inline  void    Clear();
    inline  void    Fill(const Double_t invMass, const Double_t misMass, const Double_t taggerTime, const Int_t taggerChannel);
            void    Fill(const Double_t invMass, const TLorentzVector& particle, const GTreeTagger &tagger);
    virtual void    Write(TDirectory& dir);
};



void    PHistEvent::Clear()
{
    time.Reset();
    IM.Clear();
    MM.Clear();
}
void    PHistEvent::Fill(const Double_t invMass, const Double_t misMass, const Double_t taggerTime, const Int_t taggerChannel)
{
    time.Fill(taggerTime);
    IM.Fill(invMass, taggerTime, taggerChannel);
    MM.Fill(misMass, taggerTime, taggerChannel);
}




class   PHistEvent3Meson    : public PHistEvent
{
protected:
    PHistD	sub0;
    PHistD	sub1;
    PHistD	sub2;

public:
    PHistEvent3Meson(const TString& _Name);
    virtual ~PHistEvent3Meson();

    inline  void    Clear();
    inline  void    FillSubMesons(const Double_t invMassSub0, const Double_t invMassSub1, const Double_t invMassSub2, const Double_t taggerTime, const Int_t taggerChannel);
            void    FillSubMesons(const Double_t invMassSub0, const Double_t invMassSub1, const Double_t invMassSub2, const GTreeTagger& tagger);
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



class   PHistEvent3MesonFit    : public PHistEvent3Meson
{
private:
    PHistD	ChiSq;
    PHistD	ConfidenceLevel;
    PHistD	Pull[6][4];

public:
    PHistEvent3MesonFit(const TString& _Name);
    virtual ~PHistEvent3MesonFit();

    inline  void    Clear();
    inline  void    FillFit(const Double_t _ChiSq, const Double_t _ConfidenceLevel, const Double_t* _Pull, const Double_t taggerTime, const Int_t taggerChannel);
            void    FillFit(const Double_t _ChiSq, const Double_t _ConfidenceLevel, const Double_t* _Pull, const GTreeTagger& tagger);
    virtual void    Write(TDirectory& dir);
};

void    PHistEvent3MesonFit::Clear()
{
    PHistEvent3Meson::Clear();
    ChiSq.Clear();
    ConfidenceLevel.Clear();
    for(int i=0; i<6; i++)
    {
        for(int k=0; k<4; k++)
            Pull[i][k].Clear();
    }
}

void    PHistEvent3MesonFit::FillFit(const Double_t _ChiSq, const Double_t _ConfidenceLevel, const Double_t* _Pull, const Double_t taggerTime, const Int_t taggerChannel)
{
    ChiSq.Fill(_ChiSq, taggerTime, taggerChannel);
    ConfidenceLevel.Fill(_ConfidenceLevel, taggerTime, taggerChannel);
    for(int i=0; i<6; i++)
    {
        for(int k=0; k<4; k++)
            Pull[i][k].Fill(_Pull[(4*i)+k], taggerTime, taggerChannel);
    }
}

#endif
