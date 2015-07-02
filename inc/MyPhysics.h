#ifndef __MyPhysics_h__
#define __MyPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GHistPhysics.h"



class	MyPhysics  : public GTreeManager
{
private:
    GHistParticle    proton;
    GHistParticle    etap;
    GHistParticle    etaPhotons;
    GHistParticle    pi0Photons;
    GHistParticle    allPhotons;

protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    MyPhysics();
    virtual ~MyPhysics();

    virtual Bool_t	Init(const char* configfile);

};
#endif
