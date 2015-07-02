#ifndef __GTreeManager_h__
#define __GTreeManager_h__


#include <TLorentzVector.h>
#include <TH1.h>
#include <TDatabasePDG.h>

#include "GConfigFile.h"
#include "GTreeTrack.h"
#include "GTreeTagger.h"
#include "GTreeLinPol.h"
#include "GTreeScaler.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"
#include "GTreeTrigger.h"
#include "GTreeDetectorHits.h"
#include "GTreeSetupParameters.h"
#include "GTreeEventParameters.h"
#include "GHistManager.h"

#ifdef hasPluto
#include "GTreePluto.h"
#endif
#include "GTreeA2Geant.h"

#include <stdio.h>
#include <TSystem.h>


class  GTreeManager : public GHistManager, public GConfigFile
{
private:
    TFile*      inputFile;
    TObjArray   treeList;
    TObjArray   treeCorreleatedToScalerReadList;
    TObjArray   treeSingleReadList;
    TObjArray   readList;
    TObjArray   readCorreleatedToScalerReadList;
    TObjArray   readSingleReadList;
    TObjArray   writeList;
    Bool_t      isWritten;

    UInt_t  currentScalerEntry;

    Int_t   countReconstructed;

            Bool_t      TraverseValidEvents_AcquTreeFile();
            Bool_t      TraverseValidEvents_GoATTreeFile();

    //private tree variables
    GTreeTrack*         tracks;
    GTreeTagger*        tagger;
    GTreeLinPol*        linpol;
    GTreeTrigger*       trigger;
    GTreeScaler*        scalers;

    GTreeSetupParameters* setupParameters;
    GTreeEventParameters* eventParameters;
    GTreeDetectorHits*  detectorHits;
    GTreeParticle*      rootinos;
    GTreeParticle*      photons;
    GTreeParticle*      electrons;
    GTreeParticle*      chargedPions;
    GTreeParticle*      protons;
    GTreeParticle*      neutrons;
    GTreeMeson*         neutralPions;
    GTreeMeson*         etas;
    GTreeMeson*         etaPrimes;

#ifdef hasPluto
    GTreePluto*         pluto;
#endif
    GTreeA2Geant*       geant;

protected:
    TFile*          outputFile;

    //protected tree variables Getters
    GTreeTrack*   GetTracks()                 {return tracks;}
    GTreeTagger*        GetTagger()                 {return tagger;}
    GTreeLinPol*        GetLinpol()                 {return linpol;}
    GTreeTrigger*       GetTrigger()                {return trigger;}
    GTreeScaler*        GetScalers()                {return scalers;}

    GTreeSetupParameters* GetSetupParameters()      {return setupParameters;}
    GTreeEventParameters* GetEventParameters()      {return eventParameters;}
    GTreeDetectorHits*  GetDetectorHits()           {return detectorHits;}
    GTreeParticle*      GetRootinos()               {return rootinos;}
    GTreeParticle*      GetPhotons()                {return photons;}
    GTreeParticle*      GetElectrons()              {return electrons;}
    GTreeParticle*      GetChargedPions()           {return chargedPions;}
    GTreeParticle*      GetProtons()                {return protons;}
    GTreeParticle*      GetNeutrons()               {return neutrons;}
    GTreeMeson*         GetNeutralPions()           {return neutralPions;}
    GTreeMeson*         GetEtas()                   {return etas;}
    GTreeMeson*         GetEtaPrimes()              {return etaPrimes;}

#ifdef hasPluto
    GTreePluto*         GetPluto()                  {return pluto;}
#endif
    GTreeA2Geant*       GetGeant()                  {return geant;}

    //protected tree variables const Getters
    const   GTreeTrack*   GetTracks()             const       {return tracks;}
    const   GTreeTagger*        GetTagger()             const       {return tagger;}
    const   GTreeLinPol*        GetLinpol()             const       {return linpol;}
    const   GTreeTrigger*       GetTrigger()            const       {return trigger;}
    const   GTreeScaler*        GetScalers()            const       {return scalers;}

    const   GTreeSetupParameters* GetSetupParameters()  const       {return setupParameters;}
    const   GTreeEventParameters* GetEventParameters()  const       {return eventParameters;}
    const   GTreeDetectorHits*  GetDetectorHits()       const       {return detectorHits;}
    const   GTreeParticle*      GetRootinos()           const       {return rootinos;}
    const   GTreeParticle*      GetPhotons()            const       {return photons;}
    const   GTreeParticle*      GetElectrons()          const       {return electrons;}
    const   GTreeParticle*      GetChargedPions()       const       {return chargedPions;}
    const   GTreeParticle*      GetProtons()            const       {return protons;}
    const   GTreeParticle*      GetNeutrons()           const       {return neutrons;}
    const   GTreeMeson*         GetNeutralPions()       const       {return neutralPions;}
    const   GTreeMeson*         GetEtas()               const       {return etas;}
    const   GTreeMeson*         GetEtaPrimes()          const       {return etaPrimes;}

#ifdef hasPluto
    const   GTreePluto*         GetPluto()              const       {return pluto;}
#endif
    const   GTreeA2Geant*       GetGeant()              const       {return geant;}


private:
    virtual TDirectory* GetOutputDirectory()    {return outputFile;}

protected:

    TDatabasePDG *pdgDB;

            void    FillReadList()      {for(Int_t l=0; l<readList.GetEntriesFast(); l++) ((GTree*)readList[l])->Fill();}
    const   TObjArray&  GetTreeList()    const   {return treeList;}
    virtual void    ProcessEvent() = 0;
    virtual void    ProcessScalerRead() {}
            void    SetAsGoATFile();
            void    SetAsPhysicsFile();
    virtual Bool_t  Start() = 0;
            Bool_t  TraverseEntries(const UInt_t min, const UInt_t max);
            Bool_t  TraverseScalerEntries(const UInt_t min, const UInt_t max);
            Bool_t  TraverseValidEvents()   {if(IsAcquFile()) return TraverseValidEvents_AcquTreeFile(); return TraverseValidEvents_GoATTreeFile();}
    virtual Bool_t  Write();
            Bool_t  Write(const TNamed* object);
            void    GetAndUnpack(const UInt_t event);

public:
    GTreeManager();
    virtual ~GTreeManager();

    static  Int_t   CheckInput(const char* input_filename);
            Int_t   GetEventNumber()    const   {return eventParameters->GetEventNumber();}
            TFile&  GetInputFileRef()  {return *inputFile;}
            UInt_t  GetNEntries()       const;
            Int_t   GetNReconstructed() const   {return countReconstructed;}
            UInt_t  GetNScalerEntries() const;
            Bool_t  IsAcquFile()    const;
            Bool_t  IsGoATFile()    const;
            Bool_t  IsPhysicsFile()    const;
            Bool_t  TraverseFiles();
            Bool_t  StartFile(const char* inputFileName, const char* outputFileName);

    friend  class GTree;
    friend  class GTreeParticle;
    friend  class GTreeMeson;
    friend  class GTreeTagger;
};

#endif
