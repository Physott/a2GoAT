#include "GTreeDetectorHits.h"



GTreeDetectorHits::GTreeDetectorHits(GTreeManager *Manager)    :
    GTree(Manager, TString("treeDetectorHits")),
    nNaI_Hits(0),
    nPID_Hits(0),
    nWC_Hits(0),
    nWC_Tracks(0),
    nBaF2_PbWO4_Hits(0),
    nVeto_Hits(0)
{
}

GTreeDetectorHits::~GTreeDetectorHits()
{
}

void    GTreeDetectorHits::SetBranchAdresses()
{
    tree_in->SetBranchAddress("nNaI_Hits", &nNaI_Hits);
    tree_in->SetBranchAddress("NaI_Hits", NaI_Hits);
    tree_in->SetBranchAddress("nPID_Hits", &nPID_Hits);
    tree_in->SetBranchAddress("PID_Hits", PID_Hits);
    tree_in->SetBranchAddress("nWC_Hits", &nWC_Hits);
    tree_in->SetBranchAddress("WC_Hits", WC_Hits);
    tree_in->SetBranchAddress("nWC_Tracks", &nWC_Tracks);
    tree_in->SetBranchAddress("nBaF2_PbWO4_Hits", &nBaF2_PbWO4_Hits);
    tree_in->SetBranchAddress("BaF2_PbWO4_Hits", BaF2_PbWO4_Hits);
    tree_in->SetBranchAddress("nVeto_Hits", &nVeto_Hits);
    tree_in->SetBranchAddress("Veto_Hits", Veto_Hits);
}

void    GTreeDetectorHits::SetBranches()
{
    tree_out->Branch("nNaI_Hits", &nNaI_Hits, "nNaI_Hits/I");
    tree_out->Branch("NaI_Hits", NaI_Hits, "NaI_Hits[nNaI_Hits]/I");
    tree_out->Branch("nPID_Hits", &nPID_Hits, "nPID_Hits/I");
    tree_out->Branch("PID_Hits", PID_Hits, "PID_Hits[nPID_Hits]/I");
    tree_out->Branch("nWC_Hits", &nWC_Hits, "nWC_Hits/I");
    tree_out->Branch("WC_Hits", WC_Hits, "WC_Hits[nWC_Hits]/I");
    tree_out->Branch("nWC_Tracks", &nWC_Tracks, "nWC_Tracks/I");
    tree_out->Branch("nBaF2_PbWO4_Hits", &nBaF2_PbWO4_Hits, "nBaF2_PbWO4_Hits/I");
    tree_out->Branch("BaF2_PbWO4_Hits", BaF2_PbWO4_Hits, "BaF2_PbWO4_Hits[nBaF2_PbWO4_Hits]/I");
    tree_out->Branch("nVeto_Hits", &nVeto_Hits, "nVeto_Hits/I");
    tree_out->Branch("Veto_Hits", Veto_Hits, "Veto_Hits[nVeto_Hits]/I");
}


