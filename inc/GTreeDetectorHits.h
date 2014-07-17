#ifndef __GTreeDetectorHits_h__
#define __GTreeDetectorHits_h__


#include "GTree.h"


class  GTreeDetectorHits    : public GTree
{
private:
    Int_t		nNaI_Hits;
    Int_t		NaI_Hits[720];
    Int_t		nPID_Hits;
    Int_t		PID_Hits[24];
    Int_t		nWC_Hits;
    Int_t		WC_Hits[860];
    Int_t		nWC_Tracks;
    Int_t		nBaF2_PbWO4_Hits;
    Int_t		BaF2_PbWO4_Hits[438];
    Int_t		nVeto_Hits;
    Int_t		Veto_Hits[438];

protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeDetectorHits(GTreeManager *Manager);
    virtual ~GTreeDetectorHits();

    virtual void            Clear()                         {nNaI_Hits=0; nPID_Hits=0; nWC_Hits=0; nBaF2_PbWO4_Hits=0; nVeto_Hits=0;}

            Int_t		GetNNaI_Hits()              	const	{return nNaI_Hits;}
    const	Int_t*		GetNaI_Hits()           		const	{return NaI_Hits;}
            Int_t		GetNaI_Hits(const Int_t index)	const	{return NaI_Hits[index];}

            Int_t		GetNPID_Hits()      			const	{return nPID_Hits;}
    const	Int_t*		GetPID_Hits()               	const	{return PID_Hits;}
            Int_t		GetPID_Hits(const Int_t index)	const	{return PID_Hits[index];}

            Int_t		GetNWC_Hits()       			const	{return nWC_Hits;}
    const	Int_t*		GetWC_Hits()                	const	{return WC_Hits;}
            Int_t		GetWC_Hits(const Int_t index)	const	{return WC_Hits[index];}
            Int_t		GetNWC_Tracks()       			const	{return nWC_Tracks;}

            Int_t		GetNBaF2_PbWO4_Hits()                   const	{return nBaF2_PbWO4_Hits;}
    const	Int_t*		GetBaF2_PbWO4_Hits()                    const	{return BaF2_PbWO4_Hits;}
            Int_t		GetBaF2_PbWO4_Hits(const Int_t index)	const	{return BaF2_PbWO4_Hits[index];}

            Int_t		GetNVeto_Hits()                 const	{return nVeto_Hits;}
    const	Int_t*		GetVeto_Hits()                  const	{return Veto_Hits;}
            Int_t		GetVeto_Hits(const Int_t index)	const	{return Veto_Hits[index];}
};


#endif
