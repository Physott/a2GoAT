#include "GHistScaCor.h"

#include <iostream>




GHistScaCor::GHistScaCor(const Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    buffer(0),
    accumulated(0),
    accumulatedCorrected(0),
    singleScalerReads(),
    singleScalerReadsCorrected(),
    corrected(kFALSE),
    writeUncorrected(kTRUE)
{
    singleScalerReads.SetOwner();
    singleScalerReadsCorrected.SetOwner();
}

GHistScaCor::GHistScaCor() :
    GHistLinked(),
    buffer(new TH1D()),
    accumulated(new TH1D()),
    accumulatedCorrected(new TH1D()),
    singleScalerReads(),
    singleScalerReadsCorrected(),
    corrected(kFALSE),
    writeUncorrected(kTRUE)
{
    buffer->SetDirectory(0);
    accumulated->SetDirectory(0);
    accumulatedCorrected->SetDirectory(0);

    singleScalerReads.SetOwner();
    singleScalerReadsCorrected.SetOwner();
}

GHistScaCor::GHistScaCor(const char* name, const char* title, const Int_t nbinsx, const Double_t xlow, const Double_t xup, const Bool_t linkHistogram) :
    GHistLinked(linkHistogram),
    buffer(new TH1D(TString(name).Append(GHSC_bufferNameSuffix).Data(),
                    TString(title).Append(GHSC_bufferTitleSuffix).Data(),
                    nbinsx, xlow, xup)),
    accumulated(new TH1D(TString(name).Append(GHSC_accumulatedNameSuffix).Data(),
                        TString(title).Append(GHSC_accumulatedTitleSuffix).Data(),
                        nbinsx, xlow, xup)),
    accumulatedCorrected(new TH1D(name, title, nbinsx, xlow, xup)),
    singleScalerReads(128),
    singleScalerReadsCorrected(128),
    corrected(kFALSE),
    writeUncorrected(kTRUE)
{
    buffer->SetDirectory(0);
    accumulated->SetDirectory(0);
    accumulatedCorrected->SetDirectory(0);

    singleScalerReads.SetOwner();
    singleScalerReadsCorrected.SetOwner();
}

GHistScaCor::~GHistScaCor()
{
    if(buffer)
    {
        buffer->SetDirectory(0);
        delete buffer;
    }
    if(accumulated)
    {
        accumulated->SetDirectory(0);
        delete accumulated;
    }
    if(accumulatedCorrected)
    {
        accumulatedCorrected->SetDirectory(0);
        delete accumulatedCorrected;
    }
}

Bool_t	GHistScaCor::Add(const TH1* _buffer, const TH1* _accumulated, const TH1* _accumulatedCorrected, Double_t c)
{
    buffer->Add(_buffer, c);
    accumulated->Add(_accumulated, c);
    accumulatedCorrected->Add(_accumulatedCorrected, c);
}

Bool_t	GHistScaCor::Add(const GHistScaCor *h, Double_t c)
{
    if(corrected!=h->corrected)
    {
        if(buffer->GetEntries()>0)
        {
            std::cout << "ERROR: GHistScaCor::Add. Tried to add a corrected and an uncorrected GHistScaCor" << std::endl;
            return kFALSE;
        }
    }

    if(h->corrected==kFALSE)
    {
        corrected   = kFALSE;
        buffer->Add(h->buffer, c);
        return kTRUE;
    }
    corrected   = kTRUE;

    buffer->Add(h->buffer, c);
    accumulated->Add(h->accumulated, c);
    accumulatedCorrected->Add(h->accumulatedCorrected, c);
    buffer->SetDirectory(0);
    accumulated->SetDirectory(0);
    accumulatedCorrected->SetDirectory(0);
    for(int i=0; i<h->GetNScalerReadCorrections(); i++)
    {
        if(i>=GetNScalerReadCorrections())
            CreateSingleScalerRead();
        else
        {
            ((TH1*)singleScalerReads.At(i))->Add((TH1*)h->singleScalerReads.At(i), c);
            ((TH1*)singleScalerReadsCorrected.At(i))->Add(((TH1*)h->singleScalerReadsCorrected.At(i)), c);
            ((TH1*)singleScalerReads.At(i))->SetDirectory(0);
            ((TH1*)singleScalerReadsCorrected.At(i))->SetDirectory(0);
        }
    }
    return kTRUE;
}

void    GHistScaCor::CreateSingleScalerRead()
{
    TH1D*    uncor   = new TH1D(TString(accumulatedCorrected->GetName()).Append(GHSC_singleScalerReadNameSuffix).Append(TString::Itoa(singleScalerReads.GetEntriesFast(), 10)),
                                TString(accumulatedCorrected->GetTitle()).Append(GHSC_singleScalerReadTitleSuffix).Append(TString::Itoa(singleScalerReads.GetEntriesFast(), 10)),
                                accumulatedCorrected->GetNbinsX(),
                                accumulatedCorrected->GetXaxis()->GetXmin(),
                                accumulatedCorrected->GetXaxis()->GetXmax());
    uncor->SetDirectory(0);
    singleScalerReads.Add(uncor);

    TH1D*    cor     = new TH1D(TString(accumulatedCorrected->GetName()).Append(GHSC_singleScalerReadCorrectedNameSuffix).Append(TString::Itoa(singleScalerReadsCorrected.GetEntriesFast(), 10)),
                                TString(accumulatedCorrected->GetTitle()).Append(GHSC_singleScalerReadCorrectedTitleSuffix).Append(TString::Itoa(singleScalerReadsCorrected.GetEntriesFast(), 10)),
                                accumulatedCorrected->GetNbinsX(),
                                accumulatedCorrected->GetXaxis()->GetXmin(),
                                accumulatedCorrected->GetXaxis()->GetXmax());
    cor->SetDirectory(0);
    singleScalerReadsCorrected.Add(cor);
}

Bool_t  GHistScaCor::IsEmpty()
{
    if(accumulatedCorrected->GetEntries()!=0)
        return kFALSE;
    if(buffer->GetEntries()!=0)
        return kFALSE;
    return kTRUE;
}

void    GHistScaCor::Reset(Option_t* option)
{
    buffer->SetDirectory(0);
    accumulated->SetDirectory(0);
    accumulatedCorrected->SetDirectory(0);

    buffer->Reset(option);
    accumulated->Reset(option);
    accumulatedCorrected->Reset(option);
    singleScalerReads.Clear();
    singleScalerReadsCorrected.Clear();
    corrected   = kFALSE;
}

void	GHistScaCor::Scale(Double_t c1, Option_t* option)
{
    buffer->Scale(c1, option);
    accumulated->Scale(c1, option);
    accumulatedCorrected->Scale(c1, option);
    for(int i=0; i<singleScalerReads.GetEntriesFast(); i++)
    {
        ((TH1*)singleScalerReads.At(i))->Scale(c1, option);
        ((TH1*)singleScalerReadsCorrected.At(i))->Scale(c1, option);
    }
}

void    GHistScaCor::ScalerReadCorrection(const Double_t CorrectionFactor, const Bool_t CreateHistogramsForSingleScalerReads)
{
    accumulated->Add(buffer);

    if(CreateHistogramsForSingleScalerReads)
    {
        CreateSingleScalerRead();
        ((TH1*)singleScalerReads.Last())->Add(buffer);
        buffer->Scale(CorrectionFactor);
        ((TH1*)singleScalerReadsCorrected.Last())->Add(buffer);
    }
    else
        buffer->Scale(CorrectionFactor);

    accumulatedCorrected->Add(buffer);
    buffer->Reset();
    corrected   = kTRUE;
}

void	GHistScaCor::SetBins(Int_t nx, Double_t xmin, Double_t xmax)
{
    buffer->SetBins(nx, xmin, xmax);
    accumulated->SetBins(nx, xmin, xmax);
    accumulatedCorrected->SetBins(nx, xmin, xmax);
    for(int i=0; i<singleScalerReads.GetEntriesFast(); i++)
    {
        ((TH1*)singleScalerReads.At(i))->SetBins(nx, xmin, xmax);
        ((TH1*)singleScalerReadsCorrected.At(i))->SetBins(nx, xmin, xmax);
    }
}

void    GHistScaCor::PrepareWriteList(GHistWriteList* arr, const char *name)
{
    if(!arr)
        return;

    TString     nameBuffer;
    if(name)
    {
        if(corrected==kFALSE)
            return  arr->AddHistogram(buffer, name);
        arr->AddHistogram(accumulatedCorrected, name);
    }
    else
    {
        if(corrected==kFALSE)
            return  arr->AddHistogram(buffer, accumulatedCorrected->GetName());
        arr->AddHistogram(accumulatedCorrected, accumulatedCorrected->GetName());
    }

    if(writeUncorrected==kFALSE)
        return;

    GHistWriteList* ScalerCorrection    = arr->GetDirectory(TString(GHSC_folderName));
    if(name)
    {
        nameBuffer  = name;
        nameBuffer.Append(GHSC_accumulatedNameSuffix);
        ScalerCorrection->AddHistogram(accumulated, nameBuffer);
    }
    else
        ScalerCorrection->AddHistogram(accumulated, accumulated->GetName());

    for(int i=0; i<singleScalerReads.GetEntriesFast(); i++)
    {
        GHistWriteList* SingleScalerRead    = ScalerCorrection->GetDirectory(TString(GHSC_singleScalerReadFolderName).Append(TString::Itoa(i, 10)));
        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHSC_singleScalerReadNameSuffix).Append(TString::Itoa(i, 10)),
            SingleScalerRead->AddHistogram((TH1*)singleScalerReads.At(i), nameBuffer);
        }
        else
            SingleScalerRead->AddHistogram((TH1*)singleScalerReads.At(i), ((TH1*)singleScalerReads.At(i))->GetName());
        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHSC_singleScalerReadCorrectedNameSuffix).Append(TString::Itoa(i, 10)),
            SingleScalerRead->AddHistogram((TH1*)singleScalerReadsCorrected.At(i), nameBuffer);
        }
        else
            SingleScalerRead->AddHistogram((TH1*)singleScalerReadsCorrected.At(i), ((TH1*)singleScalerReadsCorrected.At(i))->GetName());
    }
}

Int_t   GHistScaCor::WriteWithoutCalcResult(const char* name, Int_t option, Int_t bufsize)
{
    if(corrected==kFALSE)
    {
        if(name)
            return buffer->Write(name, option, bufsize);
        else
            return buffer->Write(accumulatedCorrected->GetName(), option, bufsize);
    }

    Int_t   ret = accumulatedCorrected->Write(name, option, bufsize);

    if(writeUncorrected==kFALSE)
        return ret;

    TDirectory* parentDir   = gDirectory;
    TDirectory* dir         = GetCreateDirectory(GHSC_folderName);
    dir->cd();

    TString     nameBuffer;
    if(name)
    {
        nameBuffer  = name;
        nameBuffer.Append(GHSC_accumulatedNameSuffix);
        ret += accumulated->Write(nameBuffer.Data(), option, bufsize);
    }
    else
        ret += accumulated->Write(0, option, bufsize);

    for(int i=0; i<singleScalerReads.GetEntriesFast(); i++)
    {
        dir->cd();
        GetCreateDirectory(TString(GHSC_singleScalerReadFolderName).Append(TString::Itoa(i, 10)).Data())->cd();
        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHSC_singleScalerReadNameSuffix);
            nameBuffer.Append(TString::Itoa(i, 10));
            ret += singleScalerReads.At(i)->Write(nameBuffer.Data(), option, bufsize);
        }
        else
            ret += singleScalerReads.At(i)->Write(0, option, bufsize);
        ((TH1*)singleScalerReads.At(i))->SetDirectory(0);

        if(name)
        {
            nameBuffer  = name;
            nameBuffer.Append(GHSC_singleScalerReadCorrectedNameSuffix);
            nameBuffer.Append(TString::Itoa(i, 10));
            ret += singleScalerReadsCorrected.At(i)->Write(nameBuffer.Data(), option, bufsize);
        }
        else
            ret += singleScalerReadsCorrected.At(i)->Write(0, option, bufsize);
        ((TH1*)singleScalerReadsCorrected.At(i))->SetDirectory(0);
    }
    buffer->SetDirectory(0);
    accumulated->SetDirectory(0);
    accumulatedCorrected->SetDirectory(0);
    parentDir->cd();
    return ret;
}




