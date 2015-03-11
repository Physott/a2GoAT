#ifndef __GHExample_h__
#define __GHExample_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "GH1.h"

class	GHExample  : public GTreeManager
{
private:
    GH1     taggerTime;
    GH1     IM;
    GH2     TOF;
    GH2     IMvsTheta;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    GHExample();
    virtual ~GHExample();
    virtual Bool_t  Init();

};
#endif
