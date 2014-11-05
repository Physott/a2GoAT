#ifndef __GHistScaCor_h__
#define __GHistScaCor_h__


#include <TROOT.h>
#include <TH1.h>
#include <TObjArray.h>

#include "GHistManager.h"


#define GHSC_folderName                  "ScalerCorrection"
#define GHSC_singleScalerReadFolderName  "SingleScalerRead_"
#define GHSC_bufferNameSuffix            "_Buffer"
#define GHSC_bufferTitleSuffix           ""
#define GHSC_accumulatedNameSuffix       "_NoScaCor"
#define GHSC_accumulatedTitleSuffix      " NoScaCor"
#define GHSC_singleScalerReadNameSuffix              "_NoScaCor_ScaRead"
#define GHSC_singleScalerReadTitleSuffix             " NoScaCor ScaRead "
#define GHSC_singleScalerReadCorrectedNameSuffix     "_ScaRead"
#define GHSC_singleScalerReadCorrectedTitleSuffix    " ScaRead "





class   GTreeTagger;

class   GHistScaCor : public GHistLinked
{
private:
    Bool_t      corrected;
    Bool_t      writeUncorrected;

protected:
    TH1*        buffer;
    TH1*        accumulated;
    TH1*        accumulatedCorrected;

    TObjArray   singleScalerReads;
    TObjArray   singleScalerReadsCorrected;


    GHistScaCor(const Bool_t linkHistogram);

    virtual void    CreateSingleScalerRead();

public:
    GHistScaCor();
    GHistScaCor(const char* name, const char* title, const Int_t nbinsx, const Double_t xlow, const Double_t xup, const Bool_t linkHistogram = kTRUE);
    virtual ~GHistScaCor();

    virtual Bool_t	Add(const GHistScaCor *h, Double_t c = 1);
    virtual Bool_t	Add(const TH1* _buffer, const TH1* _accumulated, const TH1* _accumulatedCorrected, const Bool_t CorrectedInput, const Double_t c = 1);
    virtual void 	CalcResult()    {}
    virtual Int_t	Fill(Double_t x)    {return buffer->Fill(x);}
    const   TH1*    GetAccumulated()            const   {return accumulated;}
    const   TH1*    GetAccumulatedCorrected()   const   {return accumulatedCorrected;}
    const   char*   GetName()                   const   {return accumulatedCorrected->GetName();}
            Int_t   GetNScalerReadCorrections() const   {return singleScalerReads.GetEntriesFast();}
            Int_t   GetNbinsX()                 const   {return accumulatedCorrected->GetNbinsX();}
            Int_t   GetXmin()                   const   {return accumulatedCorrected->GetXaxis()->GetXmin();}
            Int_t   GetXmax()                   const   {return accumulatedCorrected->GetXaxis()->GetXmax();}
            Bool_t  IsCorrected()               const   {return corrected;}
    virtual Bool_t  IsEmpty();
    virtual void    PrepareWriteList(GHistWriteList* arr, const char* name = 0);
    virtual void    Reset(Option_t* option = "");
    virtual void	Scale(Double_t c1 = 1, Option_t* option = "");
    virtual void    ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads = kFALSE);
    virtual void	SetBins(Int_t nx, Double_t xmin, Double_t xmax);
            void    SetWriteUncorrected(const Bool_t value)                                             {writeUncorrected = value;}
    virtual Int_t   WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0);
};


#endif
