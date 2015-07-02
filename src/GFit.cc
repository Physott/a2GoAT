#include "GFit.h"

using namespace std;

#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046

#include "GTreeA2Geant.h"


GFit::GFit(const char* _Name, const Bool_t linkHistogram)   :
    GHistLinked(linkHistogram),
    fitter(_Name),
    pulls(TString(_Name).Append("_pulls"), TString(_Name).Append(" pulls"), 200, -10, 10, 24, 0, 24, kFALSE),
    name(_Name),
    steps(TString(_Name).Append("_steps"), TString(_Name).Append(" steps"), 100, 0, 100, kFALSE),
    im(TString(_Name).Append("_IM"), TString(_Name).Append(" inv Mass"), 2000, 0, 2000, 48, kFALSE),
    sub0Im(TString(_Name).Append("_sub0Im"), TString(_Name).Append(" sub0 inv Mass"), 500, 300, 800, kFALSE),
    sub1Im(TString(_Name).Append("_sub1Im"), TString(_Name).Append(" sub1 inv Mass"), 300, 0, 300, kFALSE),
    sub2Im(TString(_Name).Append("_sub2Im"), TString(_Name).Append(" sub2 inv Mass"), 300, 0, 300, kFALSE),
    theta(TString(_Name).Append("_theta"), TString(_Name).Append(" theta"), 180, 0, 180, kFALSE),
    thetaCM(TString(_Name).Append("_thetaCM"), TString(_Name).Append(" thetaCM"), 180, 0, 180, 48, kFALSE),
    thetaCMTrue(TString(_Name).Append("_thetaCMTrue"), TString(_Name).Append(" thetaCMTrue"), 180, 0, 180, 48, kFALSE),
    CMtest(TString(_Name).Append("_CMtest"), TString(_Name).Append(" CMtest"), 1600, 0, 1600, kFALSE),
    phi(TString(_Name).Append("_phi"), TString(_Name).Append(" phi"), 360, -180, 180, kFALSE),
    chiSq(TString(_Name).Append("_ChiSq"), TString(_Name).Append(" ChiSq"), 1000, 0, 100, 48, kFALSE),
    confidenceLevel(TString(_Name).Append("_ConfLev"), TString(_Name).Append(" ConfLev"), 1000, 0, 1, 48, kFALSE)
{
    aplconPhotons.resize(6);

    for(size_t i=0;i<aplconPhotons.size();i++) {
        std::stringstream s;
        s << "Ph" << i;
        fitter.LinkVariable(s.str(), aplconPhotons[i].Link(), aplconPhotons[i].LinkSigma());
    }

    APLCON::Fit_Settings_t settings = fitter.GetSettings();
    settings.MaxIterations = 100;
    fitter.SetSettings(settings);
}

bool GFit::Solve(const double time, const int channel, const GTreeA2Geant& geantTree)
{
    success = false;
    result = fitter.DoFit();
    if(result.Status == APLCON::Result_Status_t::Success)
    {
        steps.Fill(result.NIterations);
        Double_t    cl  = TMath::Prob(result.ChiSquare, GetNDOF()-result.NDoF);
        confidenceLevel.Fill(cl, time);
        if(cl<0.01) return false;
        success = true;
        TLorentzVector etap(0,0,0,0);
        for(size_t t=0;t<aplconPhotons.size();t++)
            etap += FitParticle::Make(aplconPhotons[t], 0);
        im.Fill(etap.M(), time, channel);
        sub0Im.Fill((FitParticle::Make(aplconPhotons[0], 0)+FitParticle::Make(aplconPhotons[1], 0)).M(), time);
        sub1Im.Fill((FitParticle::Make(aplconPhotons[2], 0)+FitParticle::Make(aplconPhotons[3], 0)).M(), time);
        sub2Im.Fill((FitParticle::Make(aplconPhotons[4], 0)+FitParticle::Make(aplconPhotons[5], 0)).M(), time);
        theta.Fill(etap.Theta()*TMath::RadToDeg(), time);

        TLorentzVector  helpCM(etap);
        //cout << TLorentzVector(0,0,beam,beam+MASS_PROTON).BoostVector().X() << "   " << TLorentzVector(0,0,beam,beam+MASS_PROTON).BoostVector().Y() << "   " << TLorentzVector(0,0,beam,beam+MASS_PROTON).BoostVector().Z() << endl;
        TLorentzVector  helpCM3(0,0,beam,beam+MASS_PROTON);
        helpCM.Boost(-helpCM3.BoostVector());
        //cout << helpCM.Theta()*TMath::RadToDeg() << endl;
        thetaCM.Fill(helpCM.Theta()*TMath::RadToDeg(), time, channel);
        try
        {
            TLorentzVector  helpCMtrue(geantTree.GetTrueVector(2));
            helpCMtrue.Boost(-helpCM3.BoostVector());
            thetaCMTrue.Fill(helpCMtrue.Theta()*TMath::RadToDeg(), time, channel);
        }
        catch (...)
        {
        }
        helpCM3.Boost(-helpCM3.BoostVector());
        CMtest.Fill(helpCM3.P(), time);
        phi.Fill(etap.Phi()*TMath::RadToDeg(), time);
        chiSq.Fill(result.ChiSquare, time);
        for(int i=0; i<6; i++)
        {
            for(int p=0; p<3; p++)
            {
                std::stringstream s;
                s << "Ph" << i;
                s << "[" << p << "]";
                //cout << result << endl;
                //cout << s.str() << endl;
                const APLCON::Result_Variable_t& var = result.Variables.at(s.str());
                pulls.Fill(var.Pull, (3*i)+p+3, time);
            }
        }
        //cout << result << endl;
        //cout << name << " working" << endl;
        return true;
    }
    //cout << name << " not working" << endl;
    return false;
}




void    GFit::AddConstraintMM()
{
    fitter.AddConstraint("mm", {"Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5"},
                        [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5) {
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);

        return (TLorentzVector(0.0, 0.0, beam, beam + MASS_PROTON)-v0-v1-v2-v3-v4-v5).M() - MASS_PROTON;
        }
    );
}

void    GFit::AddConstraintsIM()
{
    fitter.AddConstraint("im0", {"Ph0", "Ph1"},
                        [] (vector<double>& p0, vector<double>& p1) {
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_ETA;
    }
    );

    fitter.AddConstraint("im1", {"Ph2", "Ph3"},
                        [] (vector<double>& p0, vector<double>& p1) {
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_PI0;
    }
    );

    fitter.AddConstraint("im2", {"Ph4", "Ph5"},
                        [] (vector<double>& p0, vector<double>& p1) {
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);

        return (v0+v1).M() - MASS_PI0;
    }
    );
}

    void    GFit::CalcResult()
    {
        steps.CalcResult();
        im.CalcResult();
        sub0Im.CalcResult();
        sub1Im.CalcResult();
        sub2Im.CalcResult();
        theta.CalcResult();
        thetaCM.CalcResult();
        thetaCMTrue.CalcResult();
        CMtest.CalcResult();
        phi.CalcResult();
        chiSq.CalcResult();
        confidenceLevel.CalcResult();
        pulls.CalcResult();
    }

    void    GFit::PrepareWriteList(GHistWriteList* arr, const char* name)
    {
        if(!arr)
            return;

        steps.PrepareWriteList(arr, "steps");
        im.PrepareWriteList(arr, "IM");
        sub0Im.PrepareWriteList(arr, "sub0Im");
        sub1Im.PrepareWriteList(arr, "sub1Im");
        sub2Im.PrepareWriteList(arr, "sub2Im");
        theta.PrepareWriteList(arr, "theta");
        thetaCM.PrepareWriteList(arr, "thetaCM");
        thetaCMTrue.PrepareWriteList(arr, "thetaCMTrue");
        CMtest.PrepareWriteList(arr, "CMtest");
        phi.PrepareWriteList(arr, "phi");
        chiSq.PrepareWriteList(arr, "chiSq");
        confidenceLevel.PrepareWriteList(arr, "confidenceLevel");
        pulls.PrepareWriteList(arr, "pulls");
    }

    void    GFit::Reset(Option_t* option)
    {
        steps.Reset(option);
        im.Reset(option);
        sub0Im.Reset(option);
        sub1Im.Reset(option);
        sub2Im.Reset(option);
        theta.Reset(option);
        thetaCM.Reset(option);
        thetaCMTrue.Reset(option);
        CMtest.Reset(option);
        phi.Reset(option);
        chiSq.Reset(option);
        confidenceLevel.Reset(option);
        pulls.Reset(option);
    }















    GFitVertex::GFitVertex(const char* _Name, const Bool_t linkHistogram)   :
        GFit(_Name, linkHistogram),
        vertex(TString(_Name).Append("_vertex"), TString(_Name).Append(" vertex"), 200, -10, 10, 48, kFALSE)
    {
        fitter.AddUnmeasuredVariable("Ve");
    }

    void    GFitVertex::AddConstraintMM()
    {
        fitter.AddConstraint("mm", {"Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "Ve"},
                            [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& v) {
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            ConvertTheta(p2, v);
            ConvertTheta(p3, v);
            ConvertTheta(p4, v);
            ConvertTheta(p5, v);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);
            TLorentzVector v2 = FitParticle::Make(p2, 0);
            TLorentzVector v3 = FitParticle::Make(p3, 0);
            TLorentzVector v4 = FitParticle::Make(p4, 0);
            TLorentzVector v5 = FitParticle::Make(p5, 0);

            return (TLorentzVector(0.0, 0.0, beam, beam + MASS_PROTON)-v0-v1-v2-v3-v4-v5).M() - MASS_PROTON;
            }
        );
    }

    void    GFitVertex::AddConstraintsIM()
    {
        fitter.AddConstraint("im0", {"Ph0", "Ph1", "Ve"},
                            [&] (vector<double>& p0, vector<double>& p1, vector<double>& v) {
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);

            return (v0+v1).M() - MASS_ETA;
        }
        );

        fitter.AddConstraint("im1", {"Ph2", "Ph3", "Ve"},
                            [&] (vector<double>& p0, vector<double>& p1, vector<double>& v) {
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);

            return (v0+v1).M() - MASS_PI0;
        }
        );

        fitter.AddConstraint("im2", {"Ph4", "Ph5", "Ve"},
                            [&] (vector<double>& p0, vector<double>& p1, vector<double>& v) {
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);

            return (v0+v1).M() - MASS_PI0;
        }
        );
    }

        void    GFitVertex::CalcResult()
        {
            GFit::CalcResult();
            vertex.CalcResult();
        }

        void    GFitVertex::PrepareWriteList(GHistWriteList* arr, const char* name)
        {
            if(!arr)
                return;

            GFit::PrepareWriteList(arr, "1");
            vertex.PrepareWriteList(arr, "vertex");
        }

        void    GFitVertex::Reset(Option_t* option)
        {
            GFit::Reset(option);
            vertex.Reset(option);
        }


        bool GFitVertex::Solve(const double time, const int channel, const GTreeA2Geant& geantTree)
        {
            if(GFit::Solve(time, channel, geantTree))
            {
                const APLCON::Result_Variable_t& var = result.Variables.at("Ve");
                vertex.Fill(var.Value.After, time, channel);
                return true;
            }
            return false;
        }
