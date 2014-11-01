#include "GTreeManager.h"

#include <TROOT.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TFileCacheWrite.h>

using namespace std;


GTreeManager::GTreeManager()    :
    GHistManager(),
    GConfigFile(),
    file_in(0),
    file_out(0),
    treeList(),
    treeCorreleatedToScalerReadList(),
    readList(),
    readCorreleatedToScalerReadList(),
    writeList(),
    countReconstructed(0),
    rawEvent(0),
    tagger(0),
    trigger(0),
    scalers(0),
    detectorHits(0),
    rootinos(0),
    photons(0),
    electrons(0),
    chargedPi(0),
    protons(0),
    neutrons(0),
    pi0(0),
    eta(0),
    etap(0),
    linpol(0)
{
    pdgDB = TDatabasePDG::Instance();

    etap = new GTreeMeson(this, TString(pdgDB->GetParticle("eta'")->GetName()));
    eta = new GTreeMeson(this, TString(pdgDB->GetParticle("eta")->GetName()));
    pi0 = new GTreeMeson(this, TString(pdgDB->GetParticle("pi0")->GetName()));
    rootinos = new GTreeParticle(this, TString("rootino"));
    photons = new GTreeParticle(this, TString(pdgDB->GetParticle("gamma")->GetName()));
    electrons = new GTreeParticle(this, TString(pdgDB->GetParticle("e-")->GetName()));
    chargedPi = new GTreeParticle(this, TString(pdgDB->GetParticle("pi+")->GetName()));
    protons = new GTreeParticle(this, TString(pdgDB->GetParticle("proton")->GetName()));
    neutrons = new GTreeParticle(this, TString(pdgDB->GetParticle("neutron")->GetName()));
    eventParameters = new GTreeEventParameters(this);
    detectorHits = new GTreeDetectorHits(this);
    rawEvent = new GTreeRawEvent(this);
    tagger = new GTreeTagger(this);
    trigger = new GTreeTrigger(this);
    scalers = new GTreeScaler(this);
    linpol = new GTreeLinPol(this);
}

GTreeManager::~GTreeManager()
{
    while(treeList.GetEntries()>0)
    {
        if((GTree*)treeList[0])
            delete (GTree*)treeList[0];
    }
}

Bool_t  GTreeManager::TraverseEntries(const UInt_t min, const UInt_t max)
{
    if(!file_in)
        return kFALSE;

    for(UInt_t i=min; i<max; i++)
    {
        for(int l=0; l<readList.GetEntriesFast(); l++)
            ((GTree*)readList[l])->GetEntryFast(i);

        eventParameters->SetEventNumber(i);
        countReconstructed = 0;
        ProcessEvent();
    }

    return kTRUE;
}

Bool_t  GTreeManager::TraverseScalerEntries(const UInt_t min, const UInt_t max)
{
    if(!file_in)
        return kFALSE;

    if(!scalers)
        return kFALSE;
    if(!scalers->IsOpenForInput())
    {
        if(!scalers->OpenForInput())
        {
            std::cout << "Can not open treeScaler in input file." << std::endl;
            return kFALSE;
        }
    }

    for(UInt_t i=min; i<=max; i++)
    {
        scalers->GetEntryFast(i);
        countReconstructed = 0;
        ProcessEvent();
    }

    return kTRUE;
}


Bool_t  GTreeManager::StartFile(const char* input_filename, const char* output_filename)
{
    if(file_in)    file_in->Close();
    file_in = TFile::Open(input_filename);
    if(!file_in)
    {
        cout << "#ERROR: Can not open input file " << input_filename << "!" << endl;
        return kFALSE;
    }
    cout << "Opened input file " << file_in->GetName() << "!" << file_in->GetTitle() << endl;

    for(int l=0; l<treeList.GetEntries(); l++)
    {
        if(file_in->Get(((GTree*)treeList[l])->GetName()))
            ((GTree*)treeList[l])->OpenForInput();
    }
    for(int l=0; l<treeCorreleatedToScalerReadList.GetEntries(); l++)
    {
        if(file_in->Get(((GTree*)treeCorreleatedToScalerReadList[l])->GetName()))
            ((GTree*)treeCorreleatedToScalerReadList[l])->OpenForInput();
    }

    if(file_out)
        delete file_out;
    file_out = TFile::Open(output_filename, "RECREATE");
    if(!file_out)
    {
        cout << "#ERROR: Can not create output file " << output_filename << "!" << endl;
        return kFALSE;
    }
    cout << "Created output file " << file_out->GetName() << "!" << file_out->GetTitle() << endl;
    TFileCacheWrite*    cache   = new TFileCacheWrite(file_out, 104857600);

    isWritten   = kFALSE;
    ClearLinkedHistograms();

    if(!Start())
        return kFALSE;

    if(!isWritten)
        Write();
    cache->Flush();

    /*for(int l=0; l<treeList.GetEntries(); l++)
        ((GTree*)treeList[l])->Close();
    for(int l=0; l<treeCorreleatedToScalerReadList.GetEntries(); l++)
        ((GTree*)treeCorreleatedToScalerReadList[l])->Close();*/


    //if(file_in)    file_in->Close();
    //if(file_out)    file_out->Close();
    //delete  cache;

    return kTRUE;
}

Bool_t  GTreeManager::Write()
{
    if(!file_out)   return kFALSE;
    file_out->cd();

    for(int l=0; l<writeList.GetEntries(); l++)
        ((GTree*)writeList[l])->Write();

    WriteLinkedHistograms(file_out);

    isWritten   = kTRUE;

    return kTRUE;
}

Bool_t  GTreeManager::Write(const TNamed* object)
{
    if(!file_out)   return kFALSE;
    file_out->cd();
    object->Write();
    std::cout << "object " << object->GetName() << " has been written to disk." << std::endl;
    return kTRUE;
}


Bool_t  GTreeManager::TraverseValidEvents_AcquTreeFile()
{
    if(!scalers->IsOpenForInput())
    {
        cout << "No treeScaler available. Expect MC data. Loop over all events" << endl;
        cout << "\tProcess events from " << 0 << " to " << GetNEntries() << "."<< endl;
        TraverseEntries(0, GetNEntries());
        return true;
    }
    if(scalers->GetNEntries()==0)
    {
        cout << "No treeScaler available. Expect MC data. Loop over all events" << endl;
        cout << "\tProcess events from " << 0 << " to " << GetNEntries() << "."<< endl;
        TraverseEntries(0, GetNEntries());
        return true;
    }

    if(scalers->GetNEntries()<2)
    {
        std::cout << "Less than 2 scaler reads. Can not find events with correct scalers" <<std::endl;
        return kFALSE;
    }

    // find correct shift
    int shift;
    {
        double shiftMean = 0;
        for(int l=1; l<scalers->GetNEntries(); l++)
        {
            scalers->GetEntryFast(l);
            shiftMean    += scalers->GetEventNumber() - scalers->GetEventID();
        }
        shiftMean   /= scalers->GetNEntries()-1;
        int bestIndex = 0;
        scalers->GetEntryFast(0);
        double smallestDifference = shiftMean - (scalers->GetEventNumber() - scalers->GetEventID());
        for(int l=1; l<scalers->GetNEntries(); l++)
        {
            scalers->GetEntryFast(l);
            if((shiftMean - (scalers->GetEventNumber() - scalers->GetEventID())) < smallestDifference)
            {
                bestIndex = l;
                smallestDifference  = shiftMean - (scalers->GetEventNumber() - scalers->GetEventID());
            }
        }
        scalers->GetEntryFast(bestIndex);
        shift   = scalers->GetEventNumber() - scalers->GetEventID();
    }

    file_out->cd();
    TH1I*   accepted    = new TH1I("CountScalerValid", "Events with correct scalers (all=0,accepted=1,rejected=2)", 3, 0, 3);
    accepted->SetBinContent(1, GetNEntries());

    scalers->GetEntry(scalers->GetNEntries()-1);
    int start = scalers->GetEventNumber();
    scalers->GetEntry(0);
    cout << "Checking scaler reads! Valid events from " << scalers->GetEventNumber() << " to " << start << endl;
    start = scalers->GetEventNumber();

    for(int i=1; i<GetNScalerEntries(); i++)
    {
        for(int l=0; l<readCorreleatedToScalerReadList.GetEntriesFast(); l++)
            ((GTree*)readCorreleatedToScalerReadList[l])->GetEntry(i);
        if(scalers->GetEventNumber() - scalers->GetEventID() == shift)
        {
            currentScalerEntry = i;
            accepted->SetBinContent(2, accepted->GetBinContent(2) + (scalers->GetEventNumber()-start));
            TraverseEntries(start, scalers->GetEventNumber());
            ProcessScalerRead();
            for(int l=0; l<readCorreleatedToScalerReadList.GetEntriesFast(); l++)
                ((GTree*)readCorreleatedToScalerReadList[l])->Fill();
            start = scalers->GetEventNumber();
        }
    }

    accepted->SetBinContent(3, accepted->GetBinContent(1) - accepted->GetBinContent(2));

    if(!Write(accepted))  return kFALSE;

    if(accepted)    delete accepted;
    return kTRUE;
}

Bool_t  GTreeManager::TraverseValidEvents_GoATTreeFile()
{
    if(!scalers->IsOpenForInput())
    {
        cout << "No treeScaler available. Expect MC data. Loop over all events" << endl;
        cout << "\tProcess events from " << 0 << " to " << GetNEntries() << "."<< endl;
        TraverseEntries(0, GetNEntries());
        return true;
    }
    if(scalers->GetNEntries()==0)
    {
        cout << "No treeScaler available. Expect MC data. Loop over all events" << endl;
        cout << "\tProcess events from " << 0 << " to " << GetNEntries() << "."<< endl;
        TraverseEntries(0, GetNEntries());
        return true;
    }

    Int_t   event       = 0;
    Int_t   start       = 0;
    Int_t   maxEvent    = GetNEntries();
    for(int l=0; l<readList.GetEntriesFast(); l++)
        ((GTree*)readList[l])->GetEntryFast(event);

    cout << GetNScalerEntries() << " scaler reads. " << maxEvent << " events." << endl;

    for(int i=0; i<GetNScalerEntries(); i++)
    {
        for(int l=0; l<readCorreleatedToScalerReadList.GetEntriesFast(); l++)
            ((GTree*)readCorreleatedToScalerReadList[l])->GetEntry(i);
        while(eventParameters->GetEventNumber()<scalers->GetEventNumber())
        {
            event++;
            if(event>=maxEvent)
                break;
            for(int l=0; l<readList.GetEntriesFast(); l++)
                ((GTree*)readList[l])->GetEntryFast(event);
            ProcessEvent();
        }
        if(i!=0)
        {
            if(i%10 == 0)
            {
                cout << "\t" << i << " Scaler reads processed. Events from " << start << " to " << event << "." << endl;
                start   = event;
            }
        }
        ProcessScalerRead();
    }
    cout << "\t" << GetNScalerEntries() << " Scaler reads processed. Events from " << start << " to " << event << "." << endl;
}

UInt_t  GTreeManager::GetNEntries()       const
{
    for(int l=1; l<readList.GetEntriesFast(); l++)
    {
        if(((GTree*)readList[l])->GetNEntries() != ((GTree*)readList[l-1])->GetNEntries())
        {
            cout << "ERROR: input trees have different number of entries!"<< endl;
            return 0;
        }
    }

    if(readList.GetEntriesFast() == 0)
        return 0;

    return  ((GTree*)readList[0])->GetNEntries();
}

UInt_t  GTreeManager::GetNScalerEntries()       const
{
    for(int l=1; l<readCorreleatedToScalerReadList.GetEntriesFast(); l++)
    {
        if(((GTree*)readCorreleatedToScalerReadList[l])->GetNEntries() != ((GTree*)readCorreleatedToScalerReadList[l-1])->GetNEntries())
        {
            cout << "ERROR: input trees have different number of entries!"<< endl;
            return 0;
        }
    }

    if(readCorreleatedToScalerReadList.GetEntriesFast() == 0)
        return 0;

    return  ((GTree*)readCorreleatedToScalerReadList[0])->GetNEntries();
}

void    GTreeManager::SetAsGoATFile()
{
    if(!file_out)
        return;
    TNamed flag("GoAT_File", "GoAT_File");
    Write(&flag);
}

void    GTreeManager::SetAsPhysicsFile()
{
    if(!file_out)
        return;
    TNamed flag("Physics_File", "Physics_File");
    Write(&flag);
}

Bool_t  GTreeManager::IsAcquFile()    const
{
    if(IsGoATFile())
        return kFALSE;
    if(IsPhysicsFile())
        return kFALSE;
    return kTRUE;
}

Bool_t  GTreeManager::IsGoATFile()    const
{
    if(!file_in)
        return kFALSE;
    TNamed* flag    = (TNamed*)file_in->Get("GoAT_File");
    if(flag)
        return kTRUE;
    return kFALSE;
}

Bool_t  GTreeManager::IsPhysicsFile()    const
{
    if(!file_in)
        return kFALSE;
    TNamed* flag    = (TNamed*)file_in->Get("Physics_File");
    if(flag)
        return kTRUE;
    return kFALSE;
}
