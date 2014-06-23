#ifndef __CINT__

#include "PPhysics.h"

PPhysics::PPhysics() 
{ 
}

PPhysics::~PPhysics()
{
}

Bool_t	PPhysics::Init(Char_t* configfile)
{
	
	return kTRUE;
}

void	PPhysics::Reconstruct()
{
}

void PPhysics::MissingMassPDG(const GTreeParticle& tree, PHist& hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
        particle = tree.Particle(i);
		
        for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
            FillMissingMassPair(tree, i, j, hist);
		}
	}
}	

Bool_t PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, PHist& hist)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
        FillMissingMassPair(tree, particle_index, i, hist);
	}
	return kTRUE;	
}

Bool_t PPhysics::FillMissingMassPair(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, PHist& hist)
{
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);

	Prompt = IsPrompt(time);
	Random = IsRandom(time);
	
	if ((!Prompt) && (!Random)) return kFALSE;

    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

    if (Prompt) hist.prompt->Fill(missingp4.M());
    if (Random)
    {
        if(time<Prompt_low)
            hist.rand[0]->Fill(missingp4.M());
        else
            hist.rand[1]->Fill(missingp4.M());
    }

	return kTRUE;
}

Double_t PPhysics::CalcMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.M();
}

Double_t PPhysics::CalcMissingEnergy(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree,particle_index, tagger_index);

	return missingp4.T();
}

TLorentzVector PPhysics::CalcMissingP4(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    particle	= tree.Particle(particle_index);
    beam 		= TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(tagger_index),tagger->GetPhotonBeam_E(tagger_index));
	missingp4 	= beam + target - particle;						

	return missingp4;
}

void PPhysics::FillTimePDG(const GTreeParticle& tree, TH1* Htime)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
            time = tagger->GetTagged_t(j) - tree.GetTime(i);
			Htime->Fill(time);
		}
	}
}

Bool_t PPhysics::IsPrompt(Double_t time, Double_t t_low, Double_t t_high)
{
	if ((time >= t_low) && (time <= t_high)) return kTRUE;
	
	return kFALSE;
}
	
Bool_t PPhysics::IsRandom(Double_t time, Double_t t_low1, Double_t t_high1, 
										 Double_t t_low2, Double_t t_high2 )
{

	if ((time >= t_low1) && (time <= t_high1)) return kTRUE;
	if ((time >= t_low2) && (time <= t_high2)) return kTRUE;
	
	return kFALSE;
}

Bool_t 	PPhysics::Write()
{
	return kTRUE;
}

void	PPhysics::RandomSubtraction(TH1* prompt, TH1* random, TH1* sub, Double_t ratio)
{
	sub->Add(prompt,1);
	sub->Add(random,-ratio);

}

void	PPhysics::RandomSubtraction(TH3* prompt, TH3* random, TH3* sub, Double_t ratio)
{
	sub->Add(prompt,1);
	sub->Add(random,-ratio);

}

void	PPhysics::ShowTimeCuts(TH1* timeH, TH1* cutsH, Double_t t1, Double_t t2, 
							 Double_t t3, Double_t t4, Double_t t5, Double_t t6)
{
	Double_t t;

	for (Int_t i = timeH->GetXaxis()->FindBin(t1); i < timeH->GetXaxis()->FindBin(t2); i++)
	{	
		t = timeH->GetBinContent(i);
		cutsH->SetBinContent(i,t);
	}	

	for (Int_t i = timeH->GetXaxis()->FindBin(t3); i < timeH->GetXaxis()->FindBin(t4); i++)
	{	
		t = timeH->GetBinContent(i);
		cutsH->SetBinContent(i,t);
	}	
	
	for (Int_t i = timeH->GetXaxis()->FindBin(t5); i < timeH->GetXaxis()->FindBin(t6); i++)
	{	
		t = timeH->GetBinContent(i);
		cutsH->SetBinContent(i,t);
	}	
	
}

Bool_t	PPhysics::OpenHistFile(const char* pfile, Option_t* option)
{
	HistFile = new TFile(pfile, option);
	if(!HistFile) return kFALSE;
    if(HistFile->IsZombie()) return kFALSE;
    
	cout << "HistFile " << pfile << " opened (option = "<< option << ")." << endl;
	
	return kTRUE;
}

Bool_t 	PPhysics::CloseHistFile()
{
	if(!HistFile) return kFALSE;
	HistFile->Close();
	
	return kTRUE;
}

#endif
