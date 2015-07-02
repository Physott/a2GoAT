#include "GFitBeam.h"

using namespace std;

#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046


GFitBeam::GFitBeam(const char* _Name, const Bool_t linkHistogram)   :
    GFit(_Name, linkHistogram),
    beamTheta(TString(_Name).Append("_theta"), TString(_Name).Append(" theta"), 180, 0, 180, 48, kFALSE),
    beamPhi(TString(_Name).Append("_phi"), TString(_Name).Append(" phi"), 360, -180, 180, kFALSE)
{
    fitter.LinkVariable("Be", aplconBeam.Link(), aplconBeam.LinkSigma());
}

bool GFitBeam::Solve(const double time, const int channel, const GTreeA2Geant& geantTree)
{
    if(GFit::Solve(time, channel, geantTree))
    {
        beamTheta.Fill(aplconBeam.Theta*TMath::RadToDeg(), time, channel);
        beamPhi.Fill(aplconBeam.Phi*TMath::RadToDeg(), time);
        for(int p=0; p<3; p++)
        {
            std::stringstream s;
            s << "Be[" << p << "]";
            const APLCON::Result_Variable_t& var = result.Variables.at(s.str());
            pulls.Fill(var.Pull, p, time);
        }
        return true;
    }
    //cout << name << " not working" << endl;
    return false;
}




void    GFitBeam::AddConstraintMM()
{
    fitter.AddConstraint("mm", {"Be", "Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5"},
                        [&] (vector<double>& be, vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5) {
        TLorentzVector vb = FitParticle::Make(be, 0);
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);

        return (vb+TLorentzVector(0.0, 0.0, 0.0, MASS_PROTON)-v0-v1-v2-v3-v4-v5).M() - MASS_PROTON;
        }
    );
}

void    GFitBeam::AddConstraintsIM()
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

    void    GFitBeam::CalcResult()
    {
        GFit::CalcResult();
        beamTheta.CalcResult();
        beamPhi.CalcResult();
    }

    void    GFitBeam::PrepareWriteList(GHistWriteList* arr, const char* name)
    {
        if(!arr)
            return;

        GFit::PrepareWriteList(arr);
        beamTheta.PrepareWriteList(arr, "beamTheta");
        beamPhi.PrepareWriteList(arr, "beamPhi");
    }

    void    GFitBeam::Reset(Option_t* option)
    {
        GFit::Reset(option);
        beamTheta.Reset(option);
        beamPhi.Reset(option);
    }






    GFitBeamVertex::GFitBeamVertex(const char* _Name, const Bool_t linkHistogram)   :
        GFitBeam(_Name, linkHistogram),
        vertex(TString(_Name).Append("_vertex"), TString(_Name).Append(" vertex"), 200, -10, 10, 48, kFALSE)
    {
        fitter.AddUnmeasuredVariable("Ve");
    }

    void    GFitBeamVertex::AddConstraintMM()
    {
        fitter.AddConstraint("mm", {"Be", "Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "Ve"},
                            [&] (vector<double>& be, vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& v) {
            ConvertTheta(be, v);
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            ConvertTheta(p2, v);
            ConvertTheta(p3, v);
            ConvertTheta(p4, v);
            ConvertTheta(p5, v);
            TLorentzVector vb = FitParticle::Make(be, 0);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);
            TLorentzVector v2 = FitParticle::Make(p2, 0);
            TLorentzVector v3 = FitParticle::Make(p3, 0);
            TLorentzVector v4 = FitParticle::Make(p4, 0);
            TLorentzVector v5 = FitParticle::Make(p5, 0);

            return (vb+TLorentzVector(0.0, 0.0, 0.0, MASS_PROTON)-v0-v1-v2-v3-v4-v5).M() - MASS_PROTON;
            }
        );
    }

    void    GFitBeamVertex::AddConstraintsIM()
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

        void    GFitBeamVertex::CalcResult()
        {
            GFitBeam::CalcResult();
            vertex.CalcResult();
        }

        void    GFitBeamVertex::PrepareWriteList(GHistWriteList* arr, const char* name)
        {
            if(!arr)
                return;

            GFitBeam::PrepareWriteList(arr);
            vertex.PrepareWriteList(arr, "vertex");
        }

        void    GFitBeamVertex::Reset(Option_t* option)
        {
            GFitBeam::Reset(option);
            vertex.Reset(option);
        }


        bool GFitBeamVertex::Solve(const double time, const int channel, const GTreeA2Geant& geantTree)
        {
            if(GFitBeam::Solve(time, channel, geantTree))
            {
                const APLCON::Result_Variable_t& var = result.Variables.at("Ve");
                vertex.Fill(var.Value.After, time, channel);
                return true;
            }
            return false;
        }
