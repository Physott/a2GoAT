#include "GTree.h"
#include "GTreeManager.h"


using namespace std;


GTree::GTree(GTreeManager *Manager, const TString& _Name, const Bool_t CorrelatedToScalerRead, const Bool_t SingleRead)    :
    name(_Name),
    correlatedToScalerRead(CorrelatedToScalerRead),
    singleRead(SingleRead),
    status(FLAG_CLOSED),
    inputTree(0),
    outputTree(0),
    manager(Manager),
    saveToFile(true)
{
    if(correlatedToScalerRead)
    {
        if(!manager->treeCorreleatedToScalerReadList.FindObject(this))
            manager->treeCorreleatedToScalerReadList.Add(this);
    }
    else if(singleRead)
    {
        if(!manager->treeSingleReadList.FindObject(this))
            manager->treeSingleReadList.Add(this);
    }
    else
    {
        if(!manager->treeList.FindObject(this))
            manager->treeList.Add(this);
    }
}

GTree::~GTree()
{
    if(correlatedToScalerRead)
    {
        if(manager->treeCorreleatedToScalerReadList.FindObject(this))
        {
            manager->treeCorreleatedToScalerReadList.Remove(this);
            manager->treeCorreleatedToScalerReadList.Compress();
        }
    }
    else if(singleRead)
    {
        if(manager->treeSingleReadList.FindObject(this))
        {
            manager->treeSingleReadList.Remove(this);
            manager->treeSingleReadList.Compress();
        }
    }
    else
    {
        if(manager->treeList.FindObject(this))
        {
            manager->treeList.Remove(this);
            manager->treeList.Compress();
        }
    }
}

void    GTree::Fill()
{
    if(!IsOpenForOutput())
    {
        if(!OpenForOutput())
        {
            std::cout << "Can not create " << name << " in output file." << endl;
            return;
        }
    }
    outputTree->Fill();
}

Bool_t  GTree::OpenForInput()
{
    manager->inputFile->GetObject(name.Data(),inputTree);
    if(inputTree)
    {
        SetBranchAdresses();
        status  = status | FLAG_OPENFORINPUT;
        GetEntry(0);
        if(correlatedToScalerRead)
        {
            if(!manager->readCorreleatedToScalerReadList.FindObject(this))
                manager->readCorreleatedToScalerReadList.Add(this);
        }
        else if(singleRead)
        {
            if(!manager->readSingleReadList.FindObject(this))
                manager->readSingleReadList.Add(this);
        }
        else
        {
            if(!manager->readList.FindObject(this))
                manager->readList.Add(this);
        }
        return kTRUE;
    }

    cout << "#ERROR# " << __PRETTY_FUNCTION__ << ": could not find a tree called " << name.Data() << " in input file " << manager->inputFile->GetName() << "!" << endl;
    status = status & (~FLAG_OPENFORINPUT);
    return kFALSE;
}

Bool_t  GTree::OpenForOutput()
{
    manager->outputFile->cd();
    outputTree    = new TTree(name.Data(), name.Data());
    if(outputTree)
    {
        SetBranches();
        status  = status | FLAG_OPENFOROUTPUT;
        if(!manager->writeList.FindObject(this))
            manager->writeList.Add(this);
        return kTRUE;
    }

    cout << "#ERROR# " << __PRETTY_FUNCTION__ << ": can not create output tree " << name.Data() << " in output file " << manager->inputFile->GetName() << "!" << endl;
    status = status & (~FLAG_OPENFOROUTPUT);
    return kFALSE;
}

void    GTree::Close()
{
    status = FLAG_CLOSED;

    std::cout << GetName() << std::endl;

    if(manager->writeList.FindObject(this))
    {
        manager->writeList.Remove(this);
        manager->writeList.Compress();
    }
    if(manager->readList.FindObject(this))
    {
        manager->readList.Remove(this);
        manager->readList.Compress();
    }
    if(outputTree)
        delete outputTree;

    if(!saveToFile)
    {
        manager->outputFile->Delete(TString(name).Append(";*").Data());
        saveToFile  = true;
    }
}

void    GTree::CloseForInput()
{
    status = status & ~FLAG_OPENFORINPUT;
    if(manager->readList.FindObject(this))
        manager->readList.Remove(this);
    if(inputTree)
        delete inputTree;
}

void    GTree::CloseForOutput()
{
    status = status & ~FLAG_OPENFOROUTPUT;
    if(manager->writeList.FindObject(this))
        manager->writeList.Remove(this);
    if(outputTree)
        delete outputTree;
}

void    GTree::Print() const
{
    std::cout << "GTree: Name->" << name.Data() << " Status->";
    if(!status)
        std::cout << "Closed";
    if(status & FLAG_OPENFORINPUT)
        std::cout << "Input ";
    if(status & FLAG_OPENFOROUTPUT)
        std::cout << "Output";
    std::cout << std::endl;
}

Bool_t	GTree::Write()
{
    if(!manager->outputFile)          return kFALSE;
    if(!outputTree)                   return kFALSE;
    if(!IsOpenForOutput())          return kFALSE;

    manager->outputFile->cd();
    outputTree->Write();
    std::cout << "tree " << name << " has been written to disk." << std::endl;
    return kTRUE;
}

