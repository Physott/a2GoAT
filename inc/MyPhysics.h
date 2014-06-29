#ifndef __MyPhysics_h__
#define __MyPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "P3Meson.h"
#include "GTreeManager.h"

class	MyPhysics   : virtual public GTreeManager
{
private:
    PProtonCheck    proton_eta;
    P3Meson         hist_eta;
    P3Meson         hist_eta_proton;

    PProtonCheck    proton_etap;
    P3Meson         hist_etap;
    P3Meson         hist_etap_proton;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
            Bool_t	Write();
			
public:
    MyPhysics();
    virtual ~MyPhysics();

};
#endif
