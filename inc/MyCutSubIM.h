#ifndef __MyCutSubIM_h__
#define __MyCutSubIM_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GHistEvent.h"

class	MyCutSubIM  : public GTreeManager
{
private:
    enum
    {
        isEta,
        isPi0,
        isMM
    } type;

    GH1     EPTscalers;
    GH1     EPTscalersCor;
    TH1D    EPTscalersT;
    TH1D    EPTscalersCorT;

    struct
    {
        GHistEvent3Mesons*          hist6;
        GHistEvent3MesonsProton*    hist7;
    } all, pass, fail;

    Double_t            cutSubIM[2];

protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    MyCutSubIM();
    virtual ~MyCutSubIM();

    virtual Bool_t	Init(const char* configfile);

};
#endif
