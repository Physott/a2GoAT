


void	ShowHistogram(const TH1D* data, const TH1D* mcSignal, const TH1D* mcBG)
{
	Double_t	max	= data->GetMaximum();
	
	TH1D*	both	= mcSignal->Clone();
	both->Add(mcBG);
	Double_t	scale	= max/both->GetMaximum();
	both->Scale(scale);
	both->SetLineColor(kYellow);
	both->Draw();
	
	mcSignal->Scale(scale);
	mcSignal->SetLineColor(kMagenta);
	mcSignal->Draw("SAME");
	
	mcBG->Scale(scale);
	mcBG->SetLineColor(kGreen);
	mcBG->Draw("SAME");
	
	data->SetLineColor(kBlack);
	data->Draw("SAME");
}


void	OpenHistogram(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TString& str)
{
	TH1D*	data		= (TH1D*)dataFile->Get(str);
	TH1D*	mcSignal	= (TH1D*)mcSignalFile->Get(str);
	TH1D*	mcBG		= (TH1D*)mcBGFile->Get(str);
	ShowHistogram(data, mcSignal, mcBG);
}

void	Fit(const TH1D* data, Double_t& nEtapFit, Double_t& dNEtapFit, Double_t& nEtapCut, Double_t& etapWidth, Double_t& BGWidth)
{
	nEtapCut	= 0;
	for(int i=1; i<=data->GetNbinsX(); i++)
	{
		if(data->GetBinCenter(i)>(957-8.83) && data->GetBinCenter(i)<(957+8.83))
			nEtapCut	+= data->GetBinContent(i);
	}
	
	Double_t	max	= data->GetMaximum();
	
	TF1*	fit = new TF1("fit","gaus(0) + gaus(3)",850,1020);
	fit->SetParameters(max, 958, etapWidth, max/10, 930, BGWidth);
	fit->SetParLimits(0, 0, max *2);
	fit->SetParLimits(1, 950, 970);
	fit->SetParLimits(2, etapWidth, etapWidth);
	fit->SetParLimits(3, 0, max/2);
	fit->SetParLimits(4, 900, 980);
	fit->SetParLimits(5, BGWidth, BGWidth);
	
	data->Fit(fit, "R");
	
	nEtapFit	= fit->GetParameter(0) * fit->GetParameter(2) * 1.7724;
	dNEtapFit	= fit->GetParError(0) * fit->GetParameter(2) * 1.7724;
}
void	FitProfile(const TH2D* data, const int i, Double_t& nEtapFit, Double_t& dNEtapFit, Double_t& nEtapCut, Double_t& etapWidth, Double_t& BGWidth)
{
	TH1D*	slice		= (TH1D*)data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1);
	slice->SetAxisRange(850, 1020);
	slice->Draw();
	
	Fit(slice, nEtapFit, dNEtapFit, nEtapCut, etapWidth, BGWidth);
}

void	FitSim(const TH1D* data, Double_t& nEtapFit, Double_t& dNEtapFit, Double_t& nEtapCut, Double_t& width)
{
	nEtapCut	= 0;
	for(int i=1; i<=data->GetNbinsX(); i++)
	{
		if(data->GetBinCenter(i)>(957-8.83) && data->GetBinCenter(i)<(957+8.83))
			nEtapCut	+= data->GetBinContent(i);
	}
	
	Double_t	max	= data->GetMaximum();
	
	TF1*	fit = new TF1("fit","gaus(0)",850,1020);
	fit->SetParameters(max, 958, 10, max/10, 930, 35);
	fit->SetParLimits(0, max/10, max *2);
	fit->SetParLimits(1, 900, 1000);
	fit->SetParLimits(2, 5, 75);
	
	data->Fit(fit, "R");
	
	nEtapFit	= fit->GetParameter(0) * fit->GetParameter(2) * 1.7724;
	dNEtapFit	= fit->GetParError(0) * fit->GetParameter(2) * 1.7724;
	
	width	= fit->GetParameter(2);
}
void	FitProfileSim(const TH2D* data, const int i, Double_t& nEtapFit, Double_t& dNEtapFit, Double_t& nEtapCut, Double_t& width)
{
	TH1D*	slice		= (TH1D*)(data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1))->Clone();
	slice->SetAxisRange(850, 1020);
	slice->Draw();
		
	FitSim(slice, nEtapFit, dNEtapFit, nEtapCut, width);
}

void	ReconstructionEff(const TH2D* mcSignal, TH1D* RecEff)
{
	RecEff->Reset();
	TH1D*	RecEffHelp	= new TH1D("help", "help", 48, 0, 48);
	TH1D*	slice;
	for(int i=0; i<48; i++)
	{
		slice		= (TH1D*)mcSignal->ProjectionX(TString("ReconstructionEff_Bin").Append(TString().Itoa(i,10)).Data(), i+1, i+1);
		RecEff->SetBinContent(i+1, slice->GetEntries());
		if(i<40)
			RecEffHelp->SetBinContent(i+1, 84910.0/40.0);
		else
		{
			RecEff->SetBinContent(i+1, 1);
			RecEffHelp->SetBinContent(i+1, 1);
		}	
	}
	RecEff->Divide(RecEffHelp);
	RecEffHelp->Scale(40.0/84910.0);
	RecEffHelp->SetLineColor(kMagenta);
	RecEffHelp->Draw();
	RecEff->Draw("Same");
}

void TaggEff(Double_t* array, Double_t* darray)
{
	array[0]	= 48.0/100.0;
	array[1]	= 47.0/100.0;
	array[2]	= 47.0/100.0;
	array[3]	= 49.0/100.0;
	array[4]	= 51.0/100.0;
	array[5]	= 53.0/100.0;
	array[6]	= 54.0/100.0;
	array[7]	= 55.0/100.0;
	array[8]	= 56.0/100.0;
	array[9]	= 56.0/100.0;
	array[10]	= 58.0/100.0;
	array[11]	= 57.0/100.0;
	array[12]	= 59.0/100.0;
	array[13]	= 59.0/100.0;
	array[14]	= 60.0/100.0;
	array[15]	= 58.0/100.0;
	array[16]	= 61.0/100.0;
	array[17]	= 60.0/100.0;
	array[18]	= 60.0/100.0;
	array[19]	= 61.0/100.0;
	array[20]	= 63.0/100.0;
	array[21]	= 62.0/100.0;
	array[22]	= 64.0/100.0;
	array[23]	= 64.0/100.0;
	array[24]	= 66.0/100.0;
	array[25]	= 62.0/100.0;
	array[26]	= 65.0/100.0;
	array[27]	= 62.0/100.0;
	array[28]	= 66.0/100.0;
	array[29]	= 63.0/100.0;
	array[30]	= 65.0/100.0;
	array[31]	= 63.0/100.0;
	array[32]	= 65.0/100.0;
	array[33]	= 65.0/100.0;
	array[34]	= 65.0/100.0;
	array[35]	= 63.0/100.0;
	array[36]	= 66.0/100.0;
	array[37]	= 62.0/100.0;
	array[38]	= 68.0/100.0;
	array[39]	= 65.0/100.0;
	array[40]	= 66.0/100.0;
	array[41]	= 65.0/100.0;
	array[42]	= 65.0/100.0;
	array[43]	= 64.0/100.0;
	array[44]	= 66.0/100.0;
	array[45]	= 62.0/100.0;
	array[46]	= 66.0/100.0;
	array[47]	= 0;
	for(int i=0; i<48; i++)
		darray[i]	= 2.0/100.0;
}
void	Result(const char* dataFileName, const char* mcSignalFileName, const char* mcBGFileName)
{
	TFile*	dataFile		= TFile::Open(dataFileName);
	if(!dataFile)
	{
		std::cout << "Can not open dataFile " << dataFileName << std::endl;
		return;
	}
	TFile*	mcSignalFile	= TFile::Open(mcSignalFileName);
	if(!mcSignalFile)
	{
		std::cout << "Can not open mcSignalFile " << mcSignalFileName << std::endl;
		return;
	}
	TFile*	mcBGFile		= TFile::Open(mcBGFileName);
	if(!mcBGFile)
	{
		std::cout << "Can not open mcBGFile " << mcBGFileName << std::endl;
		return;
	}
	TFile*	out				= TFile::Open("result.root", "RECREATE");
	if(!out)
	{
		std::cout << "Can not open output file result.root" << std::endl;
		return;
	}	
	
	
	TCanvas*		can	= new TCanvas("CanCheckProton", "CheckProton", 1500, 800);
	can->Divide(4,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/raw/_CheckProton_raw_prAng");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/raw/_CheckProton_raw_copl");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutCoplanarity/_CheckProton_cutCopl_prAng");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutCoplanarity/_CheckProton_cutCopl_copl");
	can->cd(3);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutProtonAngle/_CheckProtoncutPrAng_prAng");
	can->cd(7);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutProtonAngle/_CheckProtoncutPrAng_copl");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutBoth/_CheckProton_cutBoth_prAng");
	can->cd(8);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile,  "WithProton/CheckProton/cutBoth/_CheckProton_cutBoth_copl");
	
	out->cd();
	can->Write();
	
	
	can	= new TCanvas("CanRaw", "Raw", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_IM");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_MM");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub0IM");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub1IM");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub2IM");
	
	out->cd();
	can->Write();
	
	
	
	can	= new TCanvas("CanSubIM", "SubIM", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_IM");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_MM");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub0IM");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub1IM");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub2IM");
	
	out->cd();
	can->Write();



	can	= new TCanvas("CanMM", "MM", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_IM");
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_MM");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub0IM");
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub1IM");
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub2IM");
	
	out->cd();
	can->Write();
	
	
	
	can	= new TCanvas("CanBestFit", "BestFit", 1500, 800);
	can->Divide(3,2);
	Double_t*	BestFitWidth[5];
	can->cd(1);
	{
		BestFitWidth[0]	= new Double_t[4];
		TH1D*	help		= (TH1D*)mcSignalFile->Get("WithProton/MM_Cut/fit3/Final/Final_IM");
		BestFitWidth[0][0]	= help->GetRMS();
		BestFitWidth[0][1]	= help->GetRMSError();
		help->SetLineColor(kMagenta);
		help->Draw();
		help		= (TH1D*)mcBGFile->Get("WithProton/MM_Cut/fit3/Final/Final_IM");
		BestFitWidth[0][2]	= help->GetRMS();
		BestFitWidth[0][3]	= help->GetRMSError();
		help->SetLineColor(kGreen);
		help->Draw("SAME");
	}
	can->cd(2);
	{
		BestFitWidth[1]	= new Double_t[4];
		TH1D*	help		= (TH1D*)mcSignalFile->Get("WithProton/MM_Cut/fit4/Final/Final_IM");
		BestFitWidth[1][0]	= help->GetRMS();
		BestFitWidth[1][1]	= help->GetRMSError();
		help->SetLineColor(kMagenta);
		help->Draw();
		help		= (TH1D*)mcBGFile->Get("WithProton/MM_Cut/fit4/Final/Final_IM");
		BestFitWidth[1][2]	= help->GetRMS();
		BestFitWidth[1][3]	= help->GetRMSError();
		help->SetLineColor(kGreen);
		help->Draw("SAME");
	}
	can->cd(3);
	{
		BestFitWidth[2]	= new Double_t[4];
		TH1D*	help		= (TH1D*)mcSignalFile->Get("WithProton/MM_Cut/fit4Beam/Final/Final_IM");
		BestFitWidth[2][0]	= help->GetRMS();
		BestFitWidth[2][1]	= help->GetRMSError();
		help->SetLineColor(kMagenta);
		help->Draw();
		help		= (TH1D*)mcBGFile->Get("WithProton/MM_Cut/fit4Beam/Final/Final_IM");
		BestFitWidth[2][2]	= help->GetRMS();
		BestFitWidth[2][3]	= help->GetRMSError();
		help->SetLineColor(kGreen);
		help->Draw("SAME");
	}
	can->cd(4);
	{
		BestFitWidth[3]	= new Double_t[4];
		TH1D*	help		= (TH1D*)mcSignalFile->Get("WithProton/MM_Cut/fit4Proton/Final/Final_IM");
		BestFitWidth[3][0]	= help->GetRMS();
		BestFitWidth[3][1]	= help->GetRMSError();
		help->SetLineColor(kMagenta);
		help->Draw();
		help		= (TH1D*)mcBGFile->Get("WithProton/MM_Cut/fit4Proton/Final/Final_IM");
		BestFitWidth[3][2]	= help->GetRMS();
		BestFitWidth[3][3]	= help->GetRMSError();
		help->SetLineColor(kGreen);
		help->Draw("SAME");
	}
	can->cd(5);
	{
		BestFitWidth[4]	= new Double_t[4];
		TH1D*	help		= (TH1D*)mcSignalFile->Get("WithProton/MM_Cut/fit4BeamProton/Final/Final_IM");
		BestFitWidth[4][0]	= help->GetRMS();
		BestFitWidth[4][1]	= help->GetRMSError();
		help->SetLineColor(kMagenta);
		help->Draw();
		help		= (TH1D*)mcBGFile->Get("WithProton/MM_Cut/fit4BeamProton/Final/Final_IM");
		BestFitWidth[4][2]	= help->GetRMS();
		BestFitWidth[4][3]	= help->GetRMSError();
		help->SetLineColor(kGreen);
		help->Draw("SAME");
	}
	Double_t	x[48];
	Double_t	dx[48];
	can->cd(6);
	{
		Double_t	y[5];
		Double_t	dy[5];
		for(int i=0; i<5; i++)
		{
			x[i]	= i+1;
			dx[i]	= 0;
			y[i]	= BestFitWidth[i][2] / BestFitWidth[i][0];
			dy[i]	= TMath::Sqrt(((BestFitWidth[i][3] / BestFitWidth[i][0])*(BestFitWidth[i][3] / BestFitWidth[i][0])) + ((BestFitWidth[i][2]*BestFitWidth[i][1] / BestFitWidth[i][0])*(BestFitWidth[i][2]*BestFitWidth[i][1] / BestFitWidth[i][0])));
		}
		TGraphErrors*	bestFit = new TGraphErrors(5, x, y, dx, dy);
		bestFit->Draw();
	}
	
	can	= new TCanvas("CanFit", "Fit", 1500, 800);
	can->Divide(3, 2);
	TH2D*	data;
	{
		can->cd(1)->SetLogz();
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/IM");
		data->Draw("COL");
		can->cd(2)->SetLogz();
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/_MM");
		data->Draw("COL");
		can->cd(4)->SetLogz();
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/_Sub0IM");
		data->Draw("COL");
		can->cd(5)->SetLogz();
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/_Sub1IM");
		data->Draw("COL");
		can->cd(6)->SetLogz();
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/_Sub2IM");
		data->Draw("COL");
		//can->cd();
		//OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/fit4/_fit4_ChiSq");
		//can->cd(4);
		//OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/fit4/_fit4_ConfLev");*/
	}
	out->cd();
	can->Write();
	
	can	= new TCanvas("CanFitBinsSim", "FitBinsSim", 1500, 800);
	can->Divide(8,6);
	Double_t	nEtapFit[48];
	Double_t	dNEtapFit[48];
	Double_t	nEtapCut[48];
	Double_t	dNEtapCut[48];
	Double_t	etapWidth[48];
	Double_t	BGWidth[48];
	for(int i=0; i<48; i++)
	{
		can->cd(i+1);
		data		= (TH2D*)mcSignalFile->Get("WithProton/MM_Cut/fit4/Final/TaggerBinning/Final_IM_Bins");
		//data->Add((TH2D*)mcBGFile->Get("WithProton/MM_Cut/fit4/Final/TaggerBinning/Final_IM_Bins"));
		FitProfileSim(data, i, nEtapFit[i], dNEtapFit[i], nEtapCut[i], etapWidth[i]);
		dNEtapCut[i]	= sqrt(nEtapCut[i]);
		x[i]			= i;
		dx[i]			= 0;
	}
	/*can->cd(47);
	TGraphErrors*	graph = new TGraphErrors(47, x, nEtapFit, dx, dNEtapFit);
	graph->Draw();*/
	
	out->cd();
	can->Write();
	
	can	= new TCanvas("CanFitBins", "FitBins", 1500, 800);
	can->Divide(8,5);
	for(int i=0; i<48; i++)
	{
		can->cd(i+1);
		data		= (TH2D*)dataFile->Get("WithProton/MM_Cut/fit4/Final/TaggerBinning/Final_IM_Bins");
		FitProfile(data, i, nEtapFit[i], dNEtapFit[i], nEtapCut[i], etapWidth[i], BGWidth[i]);
		dNEtapCut[i]	= sqrt(nEtapCut[i]);
		x[i]			= i;
		dx[i]			= 0;
		//std::cout << nEtapFit[i] << "   " << dNEtapFit[i] << "   " << nEtapCut[i] << "   " << dNEtapCut[i] << std::endl;
	}
	out->cd();
	can->Write();
	return;
	
	
	can	= new TCanvas("CanResult", "Result", 1500, 800);
	can->Divide(2,2);
	TGraphErrors*	graph = new TGraphErrors(48, x, nEtapFit, dx, dNEtapFit);
	can->cd(1);
	graph->Draw();
	graph = new TGraphErrors(48, x, nEtapCut, dx, dNEtapCut);
	can->cd(2);
	graph->Draw();
	data		= (TH2D*)mcSignalFile->Get("WithProton/MM_Cut/fit4/TaggerBinning/_fit4_IM_Bins");
	TH1D*	RecEff	= new TH1D("RecEff", "RecEff", 48, 0, 48);
	can->cd(3);
	ReconstructionEff(data, RecEff);
	TH1D*	help	= (TH1D*)dataFile->Get("EPT_ScalerCorT");
	can->cd(4);
	help->Draw();
	
	can	= new TCanvas("CanTaggEff", "TaggEff", 1500, 800);
	can->Divide(2,0);
	Double_t	taggEff[48];
	Double_t	dTaggEff[48];
	TaggEff(taggEff, dTaggEff);
	graph = new TGraphErrors(48, x, taggEff, dx, dTaggEff);
	can->cd(1);
	graph->Draw();
	Double_t	sc[48];
	Double_t	dsc[48];
	for(int i=0; i<48; i++)
	{
		sc[i]			= help->GetBinContent(i+1) * taggEff[i];
		dsc[i]			= help->GetBinContent(i+1) * dTaggEff[i];
	}
	graph = new TGraphErrors(48, x, sc, dx, dsc);
	can->cd(2);
	graph->Draw();
	
	can	= new TCanvas("CanEndResult", "EndResult", 1500, 800);
	can->Divide(2,2);
	
	Double_t	y[48];
	Double_t	dy[48];
	for(int i=0; i<48; i++)
	{
		y[i]			= nEtapFit[i] / RecEff->GetBinContent(i+1);
		dy[i]			= dNEtapFit[i] / RecEff->GetBinContent(i+1);
	}
	graph = new TGraphErrors(48, x, y, dx, dy);
	can->cd(1);
	graph->Draw();
	for(int i=0; i<48; i++)
	{
		y[i]			= nEtapCut[i] / RecEff->GetBinContent(i+1);
		dy[i]			= dNEtapCut[i] / RecEff->GetBinContent(i+1);
	}
	graph = new TGraphErrors(48, x, y, dx, dy);
	can->cd(2);
	graph->Draw();
	
	
	//5.5e-29
	Double_t	y[47];
	Double_t	dy[47];
	for(int i=0; i<47; i++)
	{
		y[i]			= nEtapFit[i] / (RecEff->GetBinContent(i+1) * sc[i]);
		dy[i]			= dNEtapFit[i] / (RecEff->GetBinContent(i+1) * sc[i]);
	}
	graph = new TGraphErrors(47, x, y, dx, dy);
	can->cd(3);
	graph->Draw();
	for(int i=0; i<47; i++)
	{
		y[i]			= nEtapCut[i] * 3.8 / (RecEff->GetBinContent(i+1) * sc[i]);
		dy[i]			= dNEtapCut[i] * 3.8 / (RecEff->GetBinContent(i+1) * sc[i]);
	}
	graph = new TGraphErrors(47, x, y, dx, dy);
	can->cd(4);
	graph->Draw();
}
