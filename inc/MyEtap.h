#ifndef __MyEtap_h__
#define __MyEtap_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "PPhysics.h"

class	MyEtap : public PPhysics
{
private:

	Double_t time;
    TH1* 	time_etap;
    TH1* 	time_etap_cuts;

    TH1* 	MM_prompt_etap;
    TH1* 	MM_random_etap;
    TH1* 	MM_etap;
	
    TH1* 	MM_prompt_etap_n;
    TH1* 	MM_random_etap_n;
    TH1* 	MM_etap_n;

    TH1* 	MM_prompt_etap_n_6g;
    TH1* 	MM_random_etap_n_6g;
    TH1* 	MM_etap_n_6g;

    TH1* 	MM_prompt_etap_n_2g;
    TH1* 	MM_random_etap_n_2g;
    TH1* 	MM_etap_n_2g;

    TH1* 	MM_prompt_etap_c;
    TH1* 	MM_random_etap_c;
    TH1* 	MM_etap_c;
	
    TH1* 	MM_prompt_etap_c_4d;
    TH1* 	MM_random_etap_c_4d;
    TH1* 	MM_etap_c_4d;

    Int_t 	N_etap;

protected:
    virtual Bool_t  Start();
    
    virtual void    ProcessEvent();
			void	PostReconstruction();	
				
			void	DefineHistograms();
			Bool_t	WriteHistograms(TFile* pfile);
			Bool_t	WriteHistograms() {return WriteHistograms(HistFile);}
			
public:
    MyEtap();
    virtual ~MyEtap();

    virtual Bool_t	Init(const char* configfile);

};
#endif
