#include "GFitProton.h"

using namespace std;

#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046




GFitProton::GFitProton(const char* _Name, const Bool_t linkHistogram)   :
    GFit(_Name, linkHistogram),
    protonEnergy(TString(_Name).Append("_protonEnergy"), TString(_Name).Append(" protonEnergy"), 1000, 0, 1000, kFALSE),
    protonTheta(TString(_Name).Append("_protonTheta"), TString(_Name).Append(" protonTheta"), 180, 0, 180, kFALSE),
    protonPhi(TString(_Name).Append("_protonPhi"), TString(_Name).Append(" protonPhi"), 360, -180, 180, kFALSE)
{
    fitter.AddUnmeasuredVariable("PE", 250);
    fitter.LinkVariable("PA", GetLink(), GetLinkSigma());
}

void    GFitProton::AddConstraintsTotMomentum()
{
    fitter.AddConstraint("px", {"Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE"},
                        [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe) {
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

        return -v0.Px()-v1.Px()-v2.Px()-v3.Px()-v4.Px()-v5.Px()-vp.Px();
        }
    );
    fitter.AddConstraint("py", {"Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE"},
                        [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe) {
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

        return -v0.Py()-v1.Py()-v2.Py()-v3.Py()-v4.Py()-v5.Py()-vp.Py();
        }
    );
    fitter.AddConstraint("pz", {"Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE"},
                        [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe) {
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

        return beam-v0.Pz()-v1.Pz()-v2.Pz()-v3.Pz()-v4.Pz()-v5.Pz()-vp.Pz();
        }
    );
}

void    GFitProton::AddConstraintsTotEnergy()
{
    fitter.AddConstraint("e", {"Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE"},
                        [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe) {
        TLorentzVector v0 = FitParticle::Make(p0, 0);
        TLorentzVector v1 = FitParticle::Make(p1, 0);
        TLorentzVector v2 = FitParticle::Make(p2, 0);
        TLorentzVector v3 = FitParticle::Make(p3, 0);
        TLorentzVector v4 = FitParticle::Make(p4, 0);
        TLorentzVector v5 = FitParticle::Make(p5, 0);
        TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

        return beam+MASS_PROTON-v0.E()-v1.E()-v2.E()-v3.E()-v4.E()-v5.E()-vp.E();
        }
    );
}

bool GFitProton::Solve(const double time, const int channel, const GTreeA2Geant& geantTree)
{
    if(GFit::Solve(time, channel, geantTree))
    {
        const APLCON::Result_Variable_t& pe = result.Variables.at("PE");
        protonEnergy.Fill(pe.Value.After, time);
        protonTheta.Fill(aplconProtonTheta*TMath::RadToDeg(), time);
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


    void    GFitProton::CalcResult()
    {
        GFit::CalcResult();
        protonEnergy.CalcResult();
        protonTheta.CalcResult();
        protonPhi.CalcResult();
    }

    void    GFitProton::PrepareWriteList(GHistWriteList* arr, const char* name)
    {
        if(!arr)
            return;

        GFit::PrepareWriteList(arr);
        protonEnergy.PrepareWriteList(arr, "protonEnergy");
        protonTheta.PrepareWriteList(arr, "protonTheta");
        protonPhi.PrepareWriteList(arr, "protonPhi");
    }

    void    GFitProton::Reset(Option_t* option)
    {
        GFit::Reset(option);
        protonEnergy.Reset(option);
        protonTheta.Reset(option);
        protonPhi.Reset(option);
    }











    GFitProtonVertex::GFitProtonVertex(const char* _Name, const Bool_t linkHistogram)   :
        GFitProton(_Name, linkHistogram),
        vertex(TString(_Name).Append("_vertex"), TString(_Name).Append(" vertex"), 200, -10, 10, 48, kFALSE)
    {
        fitter.AddUnmeasuredVariable("Ve");
    }

    void    GFitProtonVertex::AddConstraintMM()
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

    void    GFitProtonVertex::AddConstraintsIM()
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

    void    GFitProtonVertex::AddConstraintsTotMomentum()
    {
        fitter.AddConstraint("px", {"Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE", "Ve"},
                            [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe, vector<double>& v) {
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            ConvertTheta(p2, v);
            ConvertTheta(p3, v);
            ConvertTheta(p4, v);
            ConvertTheta(p5, v);
            pa[0] = std::atan2( 25.4*sin(pa[0]), 25.4*cos(pa[0]) - v[0]);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);
            TLorentzVector v2 = FitParticle::Make(p2, 0);
            TLorentzVector v3 = FitParticle::Make(p3, 0);
            TLorentzVector v4 = FitParticle::Make(p4, 0);
            TLorentzVector v5 = FitParticle::Make(p5, 0);
            TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

            return -v0.Px()-v1.Px()-v2.Px()-v3.Px()-v4.Px()-v5.Px()-vp.Px();
            }
        );
        fitter.AddConstraint("py", {"Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE", "Ve"},
                            [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe, vector<double>& v) {
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            ConvertTheta(p2, v);
            ConvertTheta(p3, v);
            ConvertTheta(p4, v);
            ConvertTheta(p5, v);
            pa[0] = std::atan2( 25.4*sin(pa[0]), 25.4*cos(pa[0]) - v[0]);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);
            TLorentzVector v2 = FitParticle::Make(p2, 0);
            TLorentzVector v3 = FitParticle::Make(p3, 0);
            TLorentzVector v4 = FitParticle::Make(p4, 0);
            TLorentzVector v5 = FitParticle::Make(p5, 0);
            TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

            return -v0.Py()-v1.Py()-v2.Py()-v3.Py()-v4.Py()-v5.Py()-vp.Py();
            }
        );
        fitter.AddConstraint("pz", {"Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE", "Ve"},
                            [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe, vector<double>& v) {
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            ConvertTheta(p2, v);
            ConvertTheta(p3, v);
            ConvertTheta(p4, v);
            ConvertTheta(p5, v);
            pa[0] = std::atan2( 25.4*sin(pa[0]), 25.4*cos(pa[0]) - v[0]);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);
            TLorentzVector v2 = FitParticle::Make(p2, 0);
            TLorentzVector v3 = FitParticle::Make(p3, 0);
            TLorentzVector v4 = FitParticle::Make(p4, 0);
            TLorentzVector v5 = FitParticle::Make(p5, 0);
            TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

            return beam-v0.Pz()-v1.Pz()-v2.Pz()-v3.Pz()-v4.Pz()-v5.Pz()-vp.Pz();
            }
        );
    }

    void    GFitProtonVertex::AddConstraintsTotEnergy()
    {
        fitter.AddConstraint("e", {"Ph0", "Ph1", "Ph2", "Ph3", "Ph4", "Ph5", "PA", "PE", "Ve"},
                            [&] (vector<double>& p0, vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p4, vector<double>& p5, vector<double>& pa, vector<double>& pe, vector<double>& v) {
            ConvertTheta(p0, v);
            ConvertTheta(p1, v);
            ConvertTheta(p2, v);
            ConvertTheta(p3, v);
            ConvertTheta(p4, v);
            ConvertTheta(p5, v);
            pa[0] = std::atan2( 25.4*sin(pa[0]), 25.4*cos(pa[0]) - v[0]);
            TLorentzVector v0 = FitParticle::Make(p0, 0);
            TLorentzVector v1 = FitParticle::Make(p1, 0);
            TLorentzVector v2 = FitParticle::Make(p2, 0);
            TLorentzVector v3 = FitParticle::Make(p3, 0);
            TLorentzVector v4 = FitParticle::Make(p4, 0);
            TLorentzVector v5 = FitParticle::Make(p5, 0);
            TLorentzVector vp = FitParticle::Make(pe[0], pa, MASS_PROTON);

            return beam+MASS_PROTON-v0.E()-v1.E()-v2.E()-v3.E()-v4.E()-v5.E()-vp.E();
            }
        );
    }

        void    GFitProtonVertex::CalcResult()
        {
            GFitProton::CalcResult();
            vertex.CalcResult();
        }

        void    GFitProtonVertex::PrepareWriteList(GHistWriteList* arr, const char* name)
        {
            if(!arr)
                return;

            GFitProton::PrepareWriteList(arr, "1");
            vertex.PrepareWriteList(arr, "vertex");
        }

        void    GFitProtonVertex::Reset(Option_t* option)
        {
            GFitProton::Reset(option);
            vertex.Reset(option);
        }


        bool GFitProtonVertex::Solve(const double time, const int channel, const GTreeA2Geant& geantTree)
        {
            if(GFitProton::Solve(time, channel, geantTree))
            {
                const APLCON::Result_Variable_t& var = result.Variables.at("Ve");
                vertex.Fill(var.Value.After, time, channel);
                return true;
            }
            return false;
        }
