#ifndef __PPi0Example_h__
#define __PPi0Example_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "GH1.h"
#include "GHistTaggerBinning.h"

class	PPi0Example : public GTreeManager
{
private:
    GHistTaggerBinning*	test1;
    GHistTaggerBinning*	test2;
    GHistTaggerBinning*	test3;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    PPi0Example();
    virtual ~PPi0Example();

};
#endif
