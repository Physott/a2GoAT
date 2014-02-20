#ifndef __GTreeManager_h__
#define __GTreeManager_h__

#include <vector>

#include "GTreeRawEvent.h"
#include "GTreeTagger.h"
#include "GTreeScaler.h"
#include "GTreeParticle.h"



class  GTreeManager
{
private:
    TFile*  file_in;

    //Bool_t  EntryChecking(const GTree* tree);
    Bool_t  CreateParticle(GTreeParticle*& particleTree, const TString& _Name);
    Bool_t  OpenParticle(GTreeParticle*& particleTree, const TString& _Name);

protected:
    TFile*          file_out;

    GTreeRawEvent*  rawEvent;
    GTreeTagger*    tagger;
    GTreeScaler*    scalers;
    GTreeParticle*  photons;
    GTreeParticle*  protons;

            Bool_t  Create(const char* filename);
            Bool_t  CreatePhotons() {return CreateParticle(photons, TString("Photons"));}
            Bool_t  CreateProtons() {return CreateParticle(protons, TString("Protons"));}
            Bool_t  CreateRawEvent();
            Bool_t  CreateTagger();
            Bool_t  CreateScalers();
            Bool_t  Open(const char* filename);
            Bool_t  OpenPhotons()   {return OpenParticle(photons, TString("Photons"));}
            Bool_t  OpenProtons()   {return OpenParticle(protons, TString("Protons"));}
            Bool_t  OpenRawEvent();
            Bool_t  OpenTagger();
            Bool_t  OpenScalers();
    virtual void    ProcessEvent() = 0;
            //void    SetMinEntry(const UInt_t num)   {minEntry = num;}
            //void    SetNEntries(const UInt_t num)   {nEntries = num;}
            Bool_t  TraverseEntries(const UInt_t min, const UInt_t max);
            Bool_t  Write();

public:
    GTreeManager();
    virtual ~GTreeManager();

    virtual Bool_t  Process(const char* input_filename, const char* output_filename = 0) = 0;
};

#endif
