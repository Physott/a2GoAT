#include "GFitBeamProton.h"

using namespace std;

#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046




GFitBeamProton::GFitBeamProton(const char* _Name, const Bool_t linkHistogram)   :
    GFitBeam(_Name, linkHistogram),
    protonEnergy(TString(_Name).Append("_protonEnergy"), TString(_Name).Append(" protonEnergy"), 1000, 0, 1000, 48, kFALSE),
    protonTheta(TString(_Name).Append("_protonTheta"), TString(_Name).Append(" protonTheta"), 180, 0, 180, 48, kFALSE),
    protonPhi(TString(_Name).Append("_protonPhi"), TString(_Name).Append(" protonPhi"), 360, -180, 180, kFALSE)
{
    fitter.AddUnmeasuredVariable("PE", 250);
    fitter.LinkVariable("PA", GetLink(), GetLinkSigma());
}

void    GFitBeamProton::AddConstraintsTotMomentum()
{
    fitter.AddConstraint("px", {"Be", "Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE"},
                        [&] (vector<double>& be, vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe) {
        TLorentzVector vb = FitParticle::Make(be, 0);
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

        return vb.Px()-v0.Px()-v1.Px()-v2.Px()-v3.Px()-v4.Px()-v5.Px()-vp.Px();
        }
    );
    fitter.AddConstraint("py", {"Be", "Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE"},
                        [&] (vector<double>& be, vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe) {
        TLorentzVector vb = FitParticle::Make(be, 0);
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

        return vb.Py()-v0.Py()-v1.Py()-v2.Py()-v3.Py()-v4.Py()-v5.Py()-vp.Py();
        }
    );
    fitter.AddConstraint("pz", {"Be", "Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE"},
                        [&] (vector<double>& be, vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe) {
        TLorentzVector vb = FitParticle::Make(be, 0);
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

        return vb.Pz()-v0.Pz()-v1.Pz()-v2.Pz()-v3.Pz()-v4.Pz()-v5.Pz()-vp.Pz();
        }
    );
}

void    GFitBeamProton::AddConstraintsTotEnergy()
{
    fitter.AddConstraint("e", {"Be", "Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE"},
                        [&] (vector<double>& be, vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe) {
        TLorentzVector vb = FitParticle::Make(be, 0);
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

        return vb.E()+MASS_PROTON-v0.E()-v1.E()-v2.E()-v3.E()-v4.E()-v5.E()-vp.E();
        }
    );
}

bool GFitBeamProton::Solve(const double time, const int channel, const GTreeA2Geant& geantTree)
{
    if(GFitBeam::Solve(time, channel, geantTree))
    {
        const APLCON::Result_Variable_t& pe = result.Variables.at("PE");
        protonEnergy.Fill(pe.Value.After, time, channel);
        protonTheta.Fill(aplconProtonTheta*TMath::RadToDeg(), time, channel);
        protonPhi.Fill(aplconProtonPhi*TMath::RadToDeg(), time);
        for(int p=0; p<2; p++)
        {
            std::stringstream s;
            s << "PA[" << p << "]";
            const APLCON::Result_Variable_t& var = result.Variables.at(s.str());
            pulls.Fill(var.Pull, 22+p, time);
        }
        //cout << result << endl;
        //cout << name << " working" << endl;
        return true;
    }
    //cout << name << " not working" << endl;
    return false;
}


    void    GFitBeamProton::CalcResult()
    {
        GFitBeam::CalcResult();
        protonEnergy.CalcResult();
        protonTheta.CalcResult();
        protonPhi.CalcResult();
    }

    void    GFitBeamProton::PrepareWriteList(GHistWriteList* arr, const char* name)
    {
        if(!arr)
            return;

        GFitBeam::PrepareWriteList(arr);
        protonEnergy.PrepareWriteList(arr, "protonEnergy");
        protonTheta.PrepareWriteList(arr, "protonTheta");
        protonPhi.PrepareWriteList(arr, "protonPhi");
    }

    void    GFitBeamProton::Reset(Option_t* option)
    {
        GFitBeam::Reset(option);
        protonEnergy.Reset(option);
        protonTheta.Reset(option);
        protonPhi.Reset(option);
    }






    GFitBeamProtonVertex::GFitBeamProtonVertex(const char* _Name, const Bool_t linkHistogram)   :
        GFitBeamProton(_Name, linkHistogram),
        vertex(TString(_Name).Append("_vertex"), TString(_Name).Append(" vertex"), 200, -10, 10, 48, kFALSE)
    {
        fitter.AddUnmeasuredVariable("Ve");
    }

    void    GFitBeamProtonVertex::AddConstraintMM()
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

    void    GFitBeamProtonVertex::AddConstraintsIM()
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

    void    GFitBeamProtonVertex::AddConstraintsTotMomentum()
    {
        fitter.AddConstraint("px", {"Be", "Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE", "Ve"},
                            [&] (vector<double>& be, vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe, vector<double>& v) {
            ConvertTheta(be, v);
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            ConvertTheta(p2, v);
            ConvertTheta(p3, v);
            ConvertTheta(p4, v);
            ConvertTheta(p5, v);
            pa[0] = std::atan2( 25.4*sin(pa[0]), 25.4*cos(pa[0]) - v[0]);
            TLorentzVector vb = FitParticle::Make(be, 0);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);
            TLorentzVector v2 = FitParticle::Make(p2, 0);
            TLorentzVector v3 = FitParticle::Make(p3, 0);
            TLorentzVector v4 = FitParticle::Make(p4, 0);
            TLorentzVector v5 = FitParticle::Make(p5, 0);
            TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

            return vb.Px()-v0.Px()-v1.Px()-v2.Px()-v3.Px()-v4.Px()-v5.Px()-vp.Px();
            }
        );
        fitter.AddConstraint("py", {"Be", "Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE", "Ve"},
                            [&] (vector<double>& be, vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe, vector<double>& v) {
            ConvertTheta(be, v);
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            ConvertTheta(p2, v);
            ConvertTheta(p3, v);
            ConvertTheta(p4, v);
            ConvertTheta(p5, v);
            pa[0] = std::atan2( 25.4*sin(pa[0]), 25.4*cos(pa[0]) - v[0]);
            TLorentzVector vb = FitParticle::Make(be, 0);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);
            TLorentzVector v2 = FitParticle::Make(p2, 0);
            TLorentzVector v3 = FitParticle::Make(p3, 0);
            TLorentzVector v4 = FitParticle::Make(p4, 0);
            TLorentzVector v5 = FitParticle::Make(p5, 0);
            TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

            return vb.Py()-v0.Py()-v1.Py()-v2.Py()-v3.Py()-v4.Py()-v5.Py()-vp.Py();
            }
        );
        fitter.AddConstraint("pz", {"Be", "Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE", "Ve"},
                            [&] (vector<double>& be, vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe, vector<double>& v) {
            ConvertTheta(be, v);
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            ConvertTheta(p2, v);
            ConvertTheta(p3, v);
            ConvertTheta(p4, v);
            ConvertTheta(p5, v);
            pa[0] = std::atan2( 25.4*sin(pa[0]), 25.4*cos(pa[0]) - v[0]);
            TLorentzVector vb = FitParticle::Make(be, 0);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);
            TLorentzVector v2 = FitParticle::Make(p2, 0);
            TLorentzVector v3 = FitParticle::Make(p3, 0);
            TLorentzVector v4 = FitParticle::Make(p4, 0);
            TLorentzVector v5 = FitParticle::Make(p5, 0);
            TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

            return vb.Pz()-v0.Pz()-v1.Pz()-v2.Pz()-v3.Pz()-v4.Pz()-v5.Pz()-vp.Pz();
            }
        );
    }

    void    GFitBeamProtonVertex::AddConstraintsTotEnergy()
    {
        fitter.AddConstraint("e", {"Be", "Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE", "Ve"},
                            [&] (vector<double>& be, vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe, vector<double>& v) {
            ConvertTheta(be, v);
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            ConvertTheta(p2, v);
            ConvertTheta(p3, v);
            ConvertTheta(p4, v);
            ConvertTheta(p5, v);
            pa[0] = std::atan2( 25.4*sin(pa[0]), 25.4*cos(pa[0]) - v[0]);
            TLorentzVector vb = FitParticle::Make(be, 0);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);
            TLorentzVector v2 = FitParticle::Make(p2, 0);
            TLorentzVector v3 = FitParticle::Make(p3, 0);
            TLorentzVector v4 = FitParticle::Make(p4, 0);
            TLorentzVector v5 = FitParticle::Make(p5, 0);
            TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

            return vb.E()+MASS_PROTON-v0.E()-v1.E()-v2.E()-v3.E()-v4.E()-v5.E()-vp.E();
            }
        );
    }

        void    GFitBeamProtonVertex::CalcResult()
        {
            GFitBeamProton::CalcResult();
            vertex.CalcResult();
        }

        void    GFitBeamProtonVertex::PrepareWriteList(GHistWriteList* arr, const char* name)
        {
            if(!arr)
                return;

            GFitBeamProton::PrepareWriteList(arr);
            vertex.PrepareWriteList(arr, "vertex");
        }

        void    GFitBeamProtonVertex::Reset(Option_t* option)
        {
            GFitBeamProton::Reset(option);
            vertex.Reset(option);
        }


        bool GFitBeamProtonVertex::Solve(const double time, const int channel, const GTreeA2Geant& geantTree)
        {
            if(GFitBeamProton::Solve(time, channel, geantTree))
            {
                const APLCON::Result_Variable_t& var = result.Variables.at("Ve");
                vertex.Fill(var.Value.After, time, channel);
                return true;
            }
            return false;
        }
