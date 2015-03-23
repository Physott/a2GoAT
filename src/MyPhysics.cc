#include "MyPhysics.h"



MyPhysics::MyPhysics()    :
    EPTscalers("EPT_Scaler", "EPT_Scaler", 1000, 0, 100000000, 48),
    EPTscalersCor("EPT_ScalerCor", "EPT_ScalerCor", 1000, 0, 100000000, 48),
    EPTscalersT("EPT_ScalerT", "EPT_ScalerT", 48, 0, 48),
    EPTscalersCorT("EPT_ScalerCorT", "EPT_ScalerCorT", 48, 0, 48)
{ 
        GHistBGSub::InitCuts(-20, 20, -535, -35);
        GHistBGSub::AddRandCut(35, 535);
}

MyPhysics::~MyPhysics()
{

}

Bool_t	MyPhysics::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    EPTscalersT.Reset();
    EPTscalersCorT.Reset();

    TraverseValidEvents();

    outputFile->cd();
    EPTscalersT.Write();
    EPTscalersCorT.Write();

	return kTRUE;
}

void	MyPhysics::ProcessEvent()
{

}

void	MyPhysics::ProcessScalerRead()
{
    /*hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
    hist_eta_proton.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
    hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
    hist_etap_proton.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));*/

    //std::cout << scalers->GetScaler(140) * Double_t(scalers->GetScaler(1)) / scalers->GetScaler(0) << "   " << scalers->GetScaler(140) << "   " << scalers->GetScaler(0) << "   " << scalers->GetScaler(1) << std::endl;
    for(int i=140; i<188; i++)
    {
        EPTscalers.Fill(Double_t(GetScalers()->GetScaler(i)), 0, i-140);
        EPTscalersCor.Fill(GetScalers()->GetScaler(i) * Double_t(GetScalers()->GetScaler(1)) / GetScalers()->GetScaler(0), 0, i-140);
        EPTscalersT.SetBinContent(i-140+1, EPTscalersT.GetBinContent(i-140+1) + Double_t(GetScalers()->GetScaler(i)));
        EPTscalersCorT.SetBinContent(i-140+1, EPTscalersCorT.GetBinContent(i-140+1) + (GetScalers()->GetScaler(i) * Double_t(GetScalers()->GetScaler(1)) / GetScalers()->GetScaler(0)));
    }
}


Bool_t	MyPhysics::Init(const char* configfile)
{
    return kTRUE;
}
