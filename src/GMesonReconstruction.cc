#include "GMesonReconstruction.h"
#include "GParticleReconstruction.h"

using namespace std;


GMesonReconstruction::GMesonReconstruction()    :
    width_pi0(22),
    width_eta(40),
    width_etap(60)
{
}

GMesonReconstruction::~GMesonReconstruction()
{

}

Bool_t GMesonReconstruction::Start()
{
    pi0->CloseForInput();
    eta->CloseForInput();
    etap->CloseForInput();

    if(!TraverseEntries(0, photons->GetNEntries()))		return kFALSE;

    return kTRUE;
}

Bool_t	GMesonReconstruction::Init()
{

    string  config = ReadConfig("Do-Meson-Reconstruction");
    if (strcmp(config.c_str(), "nokey") == 0)
    {
        meson_theta_min = 0.0;
        meson_theta_max = 180.0;
    }
    else if(sscanf( config.c_str(), "%*d %lf %lf\n", &meson_theta_min, &meson_theta_max) == 2)
    {
        cout << "meson reconstruction is active over theta range [" <<
        meson_theta_min << "," << meson_theta_max <<"]" << endl;
    }
    else
    {
        meson_theta_min = 0.0;
        meson_theta_max = 180.0;
    }

    config = ReadConfig("Cut-IM-Width-Pi0");
    sscanf( config.c_str(), "%lf\n", &width_pi0);
    if(width_pi0) cout << "Pi0 IM width cut set to " << width_pi0 << " MeV" << endl;
    else
    {
        width_pi0 = DEFAULT_PI0_IM_WIDTH;
        cout << "Pi0 IM width cut set to default (" << width_pi0 << " MeV)" << endl;
    }

    config = ReadConfig("Cut-IM-Width-Eta");
    sscanf( config.c_str(), "%lf\n", &width_eta);
    if(width_pi0) cout << "Eta IM width cut set to " << width_eta << " MeV" << endl;
    else
    {
        width_eta = DEFAULT_ETA_IM_WIDTH;
        cout << "Pi0 IM width cut set to default (" << width_eta << " MeV)" << endl;
    }

    config = ReadConfig("Cut-IM-Width-Eta-Prime");
    sscanf( config.c_str(), "%lf\n", &width_etap);
    if(width_etap) cout << "Eta-Prime IM width cut set to " << width_etap << " MeV" << endl;
    else
    {
        width_etap = DEFAULT_ETAP_IM_WIDTH;
        cout << "Eta-Prime IM width cut set to default (" << width_etap << " MeV)" << endl;
    }
    cout << endl;

    return kTRUE;
}

Bool_t  GMesonReconstruction::ProcessEventWithoutFilling()
{
    pi0->Clear();
    eta->Clear();
    etap->Clear();

    if(GetNReconstructed()==6)
    {
        Reconstruct6g();
        return kTRUE;
    }
    else if(GetNReconstructed()==7)
    {
        Reconstruct7g();
        CheckProton();
        return kTRUE;
    }

    return kFALSE;
}

Bool_t    GMesonReconstruction::CheckProton()
{
    Double_t        smallestAngle;
    Double_t        help;
    if(etap->GetNParticles()>0)
        smallestAngle   = (tagger->GetVector(0) + TLorentzVector(0, 0, 0, MASS_PROTON) - etap->Particle(0)).Angle(protons->Particle(0).Vect());
    else
        smallestAngle   = (tagger->GetVector(0) + TLorentzVector(0, 0, 0, MASS_PROTON) - eta->Particle(0)).Angle(protons->Particle(0).Vect());
    Int_t           bestIndex       = 0;
    for(int i=i; i<tagger->GetNTagged(); i++)
    {
        if(etap->GetNParticles()>0)
            help   = (tagger->GetVector(i) + TLorentzVector(0, 0, 0, MASS_PROTON) - etap->Particle(0)).Angle(protons->Particle(0).Vect());
        else
            help   = (tagger->GetVector(i) + TLorentzVector(0, 0, 0, MASS_PROTON) - eta->Particle(0)).Angle(protons->Particle(0).Vect());
        if(help < smallestAngle)
        {
            smallestAngle   = help;
            bestIndex       = 0;
        }
    }
    foundTaggerHitForProton = bestIndex;
    if(smallestAngle > 4)
        return kFALSE;
    return kTRUE;
}

void    GMesonReconstruction::Reconstruct6g()
{
    TLorentzVector  meson[15][3];
    Double_t        help[2][3];
    Double_t        ChiSq[15][4];

    for(int i=0; i<15; i++)
    {
        meson[i][0] = photons->Particle(perm6g[i][0]) + photons->Particle(perm6g[i][1]);
        meson[i][1] = photons->Particle(perm6g[i][2]) + photons->Particle(perm6g[i][3]);
        meson[i][2] = photons->Particle(perm6g[i][4]) + photons->Particle(perm6g[i][5]);
        help[0][0]     = (MASS_ETA - meson[i][0].M())/width_eta;
        help[0][1]     = (MASS_ETA - meson[i][1].M())/width_eta;
        help[0][2]     = (MASS_ETA - meson[i][2].M())/width_eta;
        help[1][0]     = (MASS_PI0 - meson[i][0].M())/width_pi0;
        help[1][1]     = (MASS_PI0 - meson[i][1].M())/width_pi0;
        help[1][2]     = (MASS_PI0 - meson[i][2].M())/width_pi0;
        ChiSq[i][0] = (help[0][0]*help[0][0]) + (help[1][1]*help[1][1]) + (help[1][2]*help[1][2]);
        ChiSq[i][1] = (help[1][0]*help[1][0]) + (help[0][1]*help[0][1]) + (help[1][2]*help[1][2]);
        ChiSq[i][2] = (help[1][0]*help[1][0]) + (help[1][1]*help[1][1]) + (help[0][2]*help[0][2]);
        ChiSq[i][3] = (help[1][0]*help[1][0]) + (help[1][1]*help[1][1]) + (help[1][2]*help[1][2]);
    }

    minChiSq        = ChiSq[0][0];
    minDecayIndex   = 0;
    minIndex        = 0;

    for(int i=0; i<15; i++)
    {
        for(int d=0; d<4; d++)
        {
            if(ChiSq[i][d]<=minChiSq)
            {
                minChiSq        = ChiSq[i][d];
                minIndex        = i;
                minDecayIndex   = d;
            }
        }
    }

    if(minDecayIndex == 3)      //found 3Pi0
    {
        /*daughter_pdg[0] = 22;
        daughter_pdg[1] = 22;
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        //pi0->AddParticle(meson[minIndex][0], 2, daughter_pdg, daughter_index);
        daughter_index[0] = perm6g[minIndex][2];
        daughter_index[1] = perm6g[minIndex][3];
        //pi0->AddParticle(meson[minIndex][1], 2, daughter_pdg, daughter_index);
        daughter_index[0] = perm6g[minIndex][4];
        daughter_index[1] = perm6g[minIndex][5];
        //pi0->AddParticle(meson[minIndex][2], 2, daughter_pdg, daughter_index);*/

        reconstructedEta    = meson[minIndex][0] + meson[minIndex][1] + meson[minIndex][2];
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        daughter_index[2] = perm6g[minIndex][2];
        daughter_index[3] = perm6g[minIndex][3];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = &photons->Particle(perm6g[minIndex][0]);
        daughter[1] = &photons->Particle(perm6g[minIndex][1]);
        daughter[2] = &photons->Particle(perm6g[minIndex][2]);
        daughter[3] = &photons->Particle(perm6g[minIndex][3]);
        daughter[4] = &photons->Particle(perm6g[minIndex][4]);
        daughter[5] = &photons->Particle(perm6g[minIndex][5]);
        eta->AddParticle(6, daughter_index, daughter, 0, 0, 0);
        photons->Clear();
        return;
    }

    //found Eta2Pi0

    if(minDecayIndex == 0)  //Eta is meson[i][0]
    {
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        daughter_index[2] = perm6g[minIndex][2];
        daughter_index[3] = perm6g[minIndex][3];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = &photons->Particle(perm6g[minIndex][0]);
        daughter[1] = &photons->Particle(perm6g[minIndex][1]);
        daughter[2] = &photons->Particle(perm6g[minIndex][2]);
        daughter[3] = &photons->Particle(perm6g[minIndex][3]);
        daughter[4] = &photons->Particle(perm6g[minIndex][4]);
        daughter[5] = &photons->Particle(perm6g[minIndex][5]);
    }
    else if(minDecayIndex == 1)  //Eta is meson[i][1]
    {
        daughter_index[0] = perm6g[minIndex][2];
        daughter_index[1] = perm6g[minIndex][3];
        daughter_index[2] = perm6g[minIndex][0];
        daughter_index[3] = perm6g[minIndex][1];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = &photons->Particle(perm6g[minIndex][2]);
        daughter[1] = &photons->Particle(perm6g[minIndex][3]);
        daughter[2] = &photons->Particle(perm6g[minIndex][0]);
        daughter[3] = &photons->Particle(perm6g[minIndex][1]);
        daughter[4] = &photons->Particle(perm6g[minIndex][4]);
        daughter[5] = &photons->Particle(perm6g[minIndex][5]);
    }
    else if(minDecayIndex == 2)  //Eta is meson[i][2]
    {
        daughter_index[0] = perm6g[minIndex][4];
        daughter_index[1] = perm6g[minIndex][5];
        daughter_index[2] = perm6g[minIndex][0];
        daughter_index[3] = perm6g[minIndex][1];
        daughter_index[4] = perm6g[minIndex][2];
        daughter_index[5] = perm6g[minIndex][3];
        daughter[0] = &photons->Particle(perm6g[minIndex][4]);
        daughter[1] = &photons->Particle(perm6g[minIndex][5]);
        daughter[2] = &photons->Particle(perm6g[minIndex][0]);
        daughter[3] = &photons->Particle(perm6g[minIndex][1]);
        daughter[4] = &photons->Particle(perm6g[minIndex][2]);
        daughter[5] = &photons->Particle(perm6g[minIndex][3]);
    }

    reconstructedEtap   = meson[minIndex][0] + meson[minIndex][1] + meson[minIndex][2];
    etap->AddParticle(6, daughter_index, daughter, 0, 0, 0);
    photons->Clear();
}

void    GMesonReconstruction::Reconstruct6g(TLorentzVector** vec)
{
    TLorentzVector  meson[15][3];
    Double_t        help[2][3];
    Double_t        ChiSq[15][4];

    for(int i=0; i<15; i++)
    {
        meson[i][0] = *vec[perm6g[i][0]];
        meson[i][0] += *vec[perm6g[i][1]];
        meson[i][1] = *vec[perm6g[i][2]] + *vec[perm6g[i][3]];
        meson[i][2] = *vec[perm6g[i][4]] + *vec[perm6g[i][5]];
        help[0][0]     = (MASS_ETA - meson[i][0].M())/width_eta;
        help[0][1]     = (MASS_ETA - meson[i][1].M())/width_eta;
        help[0][2]     = (MASS_ETA - meson[i][2].M())/width_eta;
        help[1][0]     = (MASS_PI0 - meson[i][0].M())/width_pi0;
        help[1][1]     = (MASS_PI0 - meson[i][1].M())/width_pi0;
        help[1][2]     = (MASS_PI0 - meson[i][2].M())/width_pi0;
        ChiSq[i][0] = (help[0][0]*help[0][0]) + (help[1][1]*help[1][1]) + (help[1][2]*help[1][2]);
        ChiSq[i][1] = (help[1][0]*help[1][0]) + (help[0][1]*help[0][1]) + (help[1][2]*help[1][2]);
        ChiSq[i][2] = (help[1][0]*help[1][0]) + (help[1][1]*help[1][1]) + (help[0][2]*help[0][2]);
        ChiSq[i][3] = (help[1][0]*help[1][0]) + (help[1][1]*help[1][1]) + (help[1][2]*help[1][2]);
    }

    minChiSq        = ChiSq[0][0];
    minDecayIndex   = 0;
    minIndex        = 0;

    for(int i=0; i<15; i++)
    {
        for(int d=0; d<4; d++)
        {
            if(ChiSq[i][d]<=minChiSq)
            {
                minChiSq        = ChiSq[i][d];
                minIndex        = i;
                minDecayIndex   = d;
            }
        }
    }

    if(minDecayIndex == 3)      //found 3Pi0
    {
        /*daughter_pdg[0] = 22;
        daughter_pdg[1] = 22;
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        //pi0->AddParticle(meson[minIndex][0], 2, daughter_pdg, daughter_index);
        daughter_index[0] = perm6g[minIndex][2];
        daughter_index[1] = perm6g[minIndex][3];
        //pi0->AddParticle(meson[minIndex][1], 2, daughter_pdg, daughter_index);
        daughter_index[0] = perm6g[minIndex][4];
        daughter_index[1] = perm6g[minIndex][5];
        //pi0->AddParticle(meson[minIndex][2], 2, daughter_pdg, daughter_index);*/

        reconstructedEta    = meson[minIndex][0] + meson[minIndex][1] + meson[minIndex][2];
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        daughter_index[2] = perm6g[minIndex][2];
        daughter_index[3] = perm6g[minIndex][3];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = vec[perm6g[minIndex][0]];
        daughter[1] = vec[perm6g[minIndex][1]];
        daughter[2] = vec[perm6g[minIndex][2]];
        daughter[3] = vec[perm6g[minIndex][3]];
        daughter[4] = vec[perm6g[minIndex][4]];
        daughter[5] = vec[perm6g[minIndex][5]];
        return;
    }

    //found Eta2Pi0
    if(minDecayIndex == 0)  //Eta is meson[i][0]
    {
        daughter_index[0] = perm6g[minIndex][0];
        daughter_index[1] = perm6g[minIndex][1];
        daughter_index[2] = perm6g[minIndex][2];
        daughter_index[3] = perm6g[minIndex][3];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = vec[perm6g[minIndex][0]];
        daughter[1] = vec[perm6g[minIndex][1]];
        daughter[2] = vec[perm6g[minIndex][2]];
        daughter[3] = vec[perm6g[minIndex][3]];
        daughter[4] = vec[perm6g[minIndex][4]];
        daughter[5] = vec[perm6g[minIndex][5]];
    }
    else if(minDecayIndex == 1)  //Eta is meson[i][1]
    {
        daughter_index[0] = perm6g[minIndex][2];
        daughter_index[1] = perm6g[minIndex][3];
        daughter_index[2] = perm6g[minIndex][0];
        daughter_index[3] = perm6g[minIndex][1];
        daughter_index[4] = perm6g[minIndex][4];
        daughter_index[5] = perm6g[minIndex][5];
        daughter[0] = vec[perm6g[minIndex][2]];
        daughter[1] = vec[perm6g[minIndex][3]];
        daughter[2] = vec[perm6g[minIndex][0]];
        daughter[3] = vec[perm6g[minIndex][1]];
        daughter[4] = vec[perm6g[minIndex][4]];
        daughter[5] = vec[perm6g[minIndex][5]];
    }
    else if(minDecayIndex == 2)  //Eta is meson[i][2]
    {
        daughter_index[0] = perm6g[minIndex][4];
        daughter_index[1] = perm6g[minIndex][5];
        daughter_index[2] = perm6g[minIndex][0];
        daughter_index[3] = perm6g[minIndex][1];
        daughter_index[4] = perm6g[minIndex][2];
        daughter_index[5] = perm6g[minIndex][3];
        daughter[0] = vec[perm6g[minIndex][4]];
        daughter[1] = vec[perm6g[minIndex][5]];
        daughter[2] = vec[perm6g[minIndex][0]];
        daughter[3] = vec[perm6g[minIndex][1]];
        daughter[4] = vec[perm6g[minIndex][2]];
        daughter[5] = vec[perm6g[minIndex][3]];
    }

    reconstructedEtap   = meson[minIndex][0] + meson[minIndex][1] + meson[minIndex][2];
}

void    GMesonReconstruction::Reconstruct7g()
{
    Double_t    bestChiSq;
    Double_t    bestIndex;
    TLorentzVector* vec[6];
    for(int l=1; l<7; l++)
        vec[l] = &photons->Particle(l);
    Reconstruct6g(vec);
    bestChiSq   = minChiSq;
    bestIndex   = 0;
    for(int i=1; i<7; i++)
    {
        int k   = 0;
        for(int l=0; l<7; l++)
        {
            if(l!=i)
            {
                vec[k] = &photons->Particle(l);
                k++;
            }
        }
        Reconstruct6g(vec);
        if(minChiSq < bestChiSq)
        {
            bestChiSq   = minChiSq;
            bestIndex   = i;
        }
    }

    protons->AddParticle(SetMass(photons->Particle(bestIndex), pdgDB->GetParticle("proton")->Mass()), photons->GetApparatus(bestIndex), photons->Get_dE(bestIndex), photons->GetWC0_E(bestIndex), photons->GetWC1_E(bestIndex), photons->GetTime(bestIndex), photons->GetClusterSize(bestIndex));

    if(minDecayIndex == 3)      //found 3Pi0
    {
        eta->AddParticle(6, daughter_index, daughter, 0, 0, 0);
        photons->Clear();
        return;
    }

    etap->AddParticle(6, daughter_index, daughter, 0, 0, 0);
    photons->Clear();
    return;
}

void  GMesonReconstruction::ProcessEvent()
{
    if(!ProcessEventWithoutFilling())   return;

    eventParameters->SetNReconstructed(GetNReconstructed());
    eventParameters->Fill();

    pi0->Fill();
    eta->Fill();
    etap->Fill();
    FillReadList();
}

TLorentzVector  GMesonReconstruction::SetMass(const TLorentzVector& vec, const Double_t mass)
{
    Double_t P;
    Double_t E;
    Double_t T;

    T = vec.E() - vec.M();
    E = T + mass;
    P = TMath::Sqrt(E*E - mass*mass);

    return TLorentzVector(vec.Vect().Unit().x() * P, vec.Vect().Unit().y() * P, vec.Vect().Unit().z() * P, E);
}

Int_t		GMesonReconstruction::perm6g[15][6]=
{
    {0,1,2,3,4,5},
    {0,1,2,4,3,5},
    {0,1,2,5,4,3},

    {0,2,1,3,4,5},
    {0,2,1,4,3,5},
    {0,2,1,5,4,3},

    {0,3,2,1,4,5},
    {0,3,2,4,1,5},
    {0,3,2,5,4,1},

    {0,4,2,3,1,5},
    {0,4,2,1,3,5},
    {0,4,2,5,1,3},

    {0,5,2,3,4,1},
    {0,5,2,4,3,1},
    {0,5,2,1,4,3}
};
