#include "MyCutSubIM.h"



MyCutSubIM::MyCutSubIM()    :
    //hist_eta("eta", "eta", kTRUE),
    //hist_eta_proton("eta_proton", "eta_proton", kTRUE),
    EPTscalers("EPT_Scaler", "EPT_Scaler", 1000, 0, 100000000, 48),
    EPTscalersCor("EPT_ScalerCor", "EPT_ScalerCor", 1000, 0, 100000000, 48),
    EPTscalersT("EPT_ScalerT", "EPT_ScalerT", 48, 0, 48),
    EPTscalersCorT("EPT_ScalerCorT", "EPT_ScalerCorT", 48, 0, 48)
{ 
        GHistBGSub::InitCuts(-20, 20, -535, -35);
        GHistBGSub::AddRandCut(35, 535);

        all.hist6   = new GHistEvent3Mesons("all_6hits", "all_6hits", kTRUE);
        all.hist7   = new GHistEvent3MesonsProton("all_7hits", "all_7hits", kTRUE);
        pass.hist6   = new GHistEvent3Mesons("pass_6hits", "pass_6hits", kTRUE);
        pass.hist7   = new GHistEvent3MesonsProton("pass_7hits", "pass_7hits", kTRUE);
        fail.hist6   = new GHistEvent3Mesons("fail_6hits", "fail_6hits", kTRUE);
        fail.hist7   = new GHistEvent3MesonsProton("fail_7hits", "fail_7hits", kTRUE);
}

MyCutSubIM::~MyCutSubIM()
{
    if(all.hist6)    delete all.hist6;
    if(all.hist7)    delete all.hist7;
    if(pass.hist6)   delete pass.hist6;
    if(pass.hist7)   delete pass.hist7;
    if(fail.hist6)   delete fail.hist6;
    if(fail.hist7)   delete fail.hist7;
}

Bool_t	MyCutSubIM::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsGoATFile();

    EPTscalersT.Reset();
    EPTscalersCorT.Reset();

    TraverseValidEvents();

    outputFile->cd();
    EPTscalersT.Write();
    EPTscalersCorT.Write();

	return kTRUE;
}

void	MyCutSubIM::ProcessEvent()
{
    if(GetEtaPrimes()->GetNParticles()>0)
    {
        TLorentzVector  helpEta  = GetPhotons()->Particle(0);
                        helpEta += GetPhotons()->Particle(1);

        Double_t    im  = GetEtaPrimes()->Particle(0).M();
        Double_t    theta  = GetEtaPrimes()->Particle(0).Theta()*TMath::RadToDeg();
        Double_t    phi  = GetEtaPrimes()->Particle(0).Phi()*TMath::RadToDeg();
        Double_t    mm;
        Double_t    sub_im[3];
        sub_im[0]   = (GetPhotons()->Particle(0) + GetPhotons()->Particle(1)).M();
        sub_im[1]   = (GetPhotons()->Particle(2) + GetPhotons()->Particle(3)).M();
        sub_im[2]   = (GetPhotons()->Particle(4) + GetPhotons()->Particle(5)).M();

        if(GetProtons()->GetNParticles()>0)
        {
            for(int i=0; i<GetTagger()->GetNTagged(); i++)
            {
                mm  = (GetTagger()->GetVectorProtonTarget(i)-GetEtaPrimes()->Particle(0)).M();
                TLorentzVector  helpMeson(GetEtaPrimes()->Particle(0));
                helpMeson.Boost(-GetTagger()->GetVectorProtonTarget(i).BoostVector());
                TLorentzVector  helpProton(GetProtons()->Particle(0));
                helpProton.Boost(-GetTagger()->GetVectorProtonTarget(i).BoostVector());

                all.hist7->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), GetProtons()->GetClusterEnergy(0), GetProtons()->GetTheta(0), GetProtons()->GetPhi(0), helpMeson.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));

                switch(type)
                {
                case isEta:
                    {
                        if((sub_im[0]>cutSubIM[0] && sub_im[0]<cutSubIM[1]))
                        {
                            pass.hist7->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), GetProtons()->GetClusterEnergy(0), GetProtons()->GetTheta(0), GetProtons()->GetPhi(0), helpProton.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                            FillReadList();
                        }
                        else
                            fail.hist7->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), GetProtons()->GetClusterEnergy(0), GetProtons()->GetTheta(0), GetProtons()->GetPhi(0), helpProton.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                    }
                    break;
                case isPi0:
                    {
                        if((sub_im[1]>cutSubIM[0] && sub_im[1]<cutSubIM[1]) &&
                           (sub_im[2]>cutSubIM[0] && sub_im[2]<cutSubIM[1]))
                        {
                            pass.hist7->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), GetProtons()->GetClusterEnergy(0), GetProtons()->GetTheta(0), GetProtons()->GetPhi(0), helpProton.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                            FillReadList();
                        }
                        else
                            fail.hist7->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), GetProtons()->GetClusterEnergy(0), GetProtons()->GetTheta(0), GetProtons()->GetPhi(0), helpProton.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                    }
                    break;
                case isMM:
                    {
                        if(mm>cutSubIM[0] && mm<cutSubIM[1])
                        {
                            pass.hist7->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), GetProtons()->GetClusterEnergy(0), GetProtons()->GetTheta(0), GetProtons()->GetPhi(0), helpProton.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                            FillReadList();
                        }
                        else
                            fail.hist7->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), GetProtons()->GetClusterEnergy(0), GetProtons()->GetTheta(0), GetProtons()->GetPhi(0), helpProton.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                    }
                    break;
                }
            }
        }
        else
        {
            for(int i=0; i<GetTagger()->GetNTagged(); i++)
            {
                mm  = (GetTagger()->GetVectorProtonTarget(i)-GetEtaPrimes()->Particle(0)).M();
                TLorentzVector  helpMeson(GetEtaPrimes()->Particle(0));
                helpMeson.Boost(-GetTagger()->GetVectorProtonTarget(i).BoostVector());
                TLorentzVector  helpProton(GetProtons()->Particle(0));
                helpProton.Boost(-GetTagger()->GetVectorProtonTarget(i).BoostVector());

                all.hist6->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));

                switch(type)
                {
                case isEta:
                    {
                        if((sub_im[0]>cutSubIM[0] && sub_im[0]<cutSubIM[1]))
                        {
                            pass.hist6->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                            FillReadList();
                        }
                        else
                            fail.hist6->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                    }
                    break;
                case isPi0:
                    {
                        if((sub_im[1]>cutSubIM[0] && sub_im[1]<cutSubIM[1]) &&
                           (sub_im[2]>cutSubIM[0] && sub_im[2]<cutSubIM[1]))
                        {
                            pass.hist6->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                            FillReadList();
                        }
                        else
                            fail.hist6->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                    }
                    break;
                case isMM:
                    {
                        if(mm>cutSubIM[0] && mm<cutSubIM[1])
                        {
                            pass.hist6->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                            FillReadList();
                        }
                        else
                            fail.hist6->Fill(im, mm, theta, phi, helpMeson.Theta()*TMath::RadToDeg(), sub_im[0], sub_im[1], sub_im[2], GetTagger()->GetTaggedTime(i), GetTagger()->GetTaggedChannel(i));
                    }
                    break;
                }
            }
        }
    }
}

void	MyCutSubIM::ProcessScalerRead()
{
    for(int i=140; i<188; i++)
    {
        EPTscalers.Fill(Double_t(GetScalers()->GetScaler(i)), 0, i-140);
        EPTscalersCor.Fill(GetScalers()->GetScaler(i) * Double_t(GetScalers()->GetScaler(1)) / GetScalers()->GetScaler(0), 0, i-140);
        EPTscalersT.SetBinContent(i-140+1, EPTscalersT.GetBinContent(i-140+1) + Double_t(GetScalers()->GetScaler(i)));
        EPTscalersCorT.SetBinContent(i-140+1, EPTscalersCorT.GetBinContent(i-140+1) + (GetScalers()->GetScaler(i) * Double_t(GetScalers()->GetScaler(1)) / GetScalers()->GetScaler(0)));
    }
}


Bool_t	MyCutSubIM::Init(const char* configfile)
{
    SetConfigFile(configfile);
    std::string config;
    Double_t    buf[2];



    config = ReadConfig("Cut-Eta-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            type        = isEta;
            cutSubIM[0] = buf[0];
            cutSubIM[1] = buf[1];
            cout << "Set Cuts for eta: " << buf[0] << "   " << buf[1] << endl;
        }
    }

    config = ReadConfig("Cut-Pi0-SubIM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            type        = isPi0;
            cutSubIM[0] = buf[0];
            cutSubIM[1] = buf[1];
            cout << "Set Cuts for pi0: " << buf[0] << "   " << buf[1] << endl;
        }
    }

    config = ReadConfig("Cut-MM");
    if (strcmp(config.c_str(), "nokey") != 0)
    {
        if(sscanf( config.c_str(), "%lf %lf\n", &buf[0], &buf[1]) == 2)
        {
            type        = isMM;
            cutSubIM[0] = buf[0];
            cutSubIM[1] = buf[1];
            cout << "Set Cuts for miss mass: " << buf[0] << "   " << buf[1] << endl;
        }
    }

    return true;
}
