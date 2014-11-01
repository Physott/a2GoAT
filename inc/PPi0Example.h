#ifndef __PPi0Example_h__
#define __PPi0Example_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "GH1.h"

class	PPi0Example : public GTreeManager
{
private:
    GH1*	test1;
    GH1*	test2;
    GH1*	test3;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    PPi0Example();
    virtual ~PPi0Example();

};
#endif
