#ifndef __GHistManager_h__
#define __GHistManager_h__



#include <TObject.h>
#include <TObjArray.h>
#include <TH1D.h>
#include <TH1I.h>



class   GLinkedHist;

class   GHistManager
{
private:
    TObjArray   histList;

    void        AddHistogramToList(GLinkedHist* hist);
    //virtual TDirectory* GetOutputDirectory() = 0;
    void        RemoveHistogramFromList(GLinkedHist* hist);

protected:

public:
    GHistManager();
    virtual ~GHistManager();

    virtual void    Clear();
            Bool_t  WriteLinkedHistograms(TDirectory* dir);

    friend  class   GLinkedHist;
};





class   GLinkedHist : public TObject
{
public:
    GLinkedHist(const Bool_t linked = true);
    virtual ~GLinkedHist();

    virtual void	Write(TDirectory *dir) = 0;
};





/*
class   GLinkedHist1D   : public TH1D, public GLinkedHist
{
public:
    GLinkedHist1D(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max);
    virtual ~GLinkedHist1D();
};

class   GLinkedHist1I   : public TH1I, public GLinkedHist
{
public:
    GLinkedHist1I(const TString& name, const TString& title, const Int_t nBins, const Double_t min, const Double_t max);
    virtual ~GLinkedHist1I();
};
*/
#endif
