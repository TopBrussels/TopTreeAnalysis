{
	
	TFile* f = new TFile("../store/mc.root");
	
	// S OVER S+B
	
	TH1F* sob = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1SoverSB_Var0");
	
	TCanvas* c1 = new TCanvas("c1","c1",800,600);
	
	gPad->SetGrid();
	
	sob->SetTitle("");
	sob->GetYaxis()->SetRangeUser(0,0.8);
	sob->GetXaxis()->SetTitle("m_{#mu j} (GeV/c^{2})");
	sob->GetYaxis()->SetTitle("b purity");
	
	sob->Draw();
	
	c1->Update();
	
	TPaveStats* sb = (TPaveStats*)sob->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	TLatex* text = new TLatex(0.13,0.87,"CMS Simulation");
	
	text->SetTextSize(0.05);
	
	text->SetNDC();
	
	text->Draw();
	
	c1->Draw();
	
	c1->SaveAs("sovern.root");
	c1->SaveAs("sovern.pdf");
	c1->SaveAs("sovern.png");
	
	// btag eff all vs SS
	
	TH1F* all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffAll");
	TH1F* ss = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffMC");
	
	TCanvas* c2 = new TCanvas("c2","c2",800,600);
	
	c2->cd();
	
	gPad->SetGrid();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("TCHE discriminant");
	all->GetYaxis()->SetTitle("Efficiency");
	
	all->Draw("hist");
	
	c2->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c2->Draw();
	
	ss->SetLineColor(kRed);
	ss->SetMarkerColor(kRed);
	ss->Draw("same");
	
	TLegend* leg = new TLegend(0.45,0.65,0.8,0.8);
	
	leg->AddEntry(all,"b candidate sample","l");
	leg->AddEntry(ss,"b enriched subsample","p");
	
	leg->Draw();
    
    text->Draw();
	
	c2->SaveAs("effcomp.pdf");
	
	// METHOD VS MC 
	
	all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll");
	ss = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasured");
	
	TCanvas* c3 = new TCanvas("c3","c3",800,600);
	
	c3->cd();
	
	gPad->SetGrid();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("TCHE discriminant");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->Scale(1./all->Integral());
	all->Draw("hist");
	
	c3->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c3->Draw();
	
	ss->SetLineColor(kRed);
	ss->SetMarkerColor(kRed);
	ss->DrawNormalized("same");
	
	leg = new TLegend(0.45,0.65,0.8,0.8);
	
	leg->AddEntry(all,"MC Truth","l");
	leg->AddEntry(ss,"Method","lp");
	
	leg->Draw();
    
    text->Draw();
	
	c3->SaveAs("measMC-truth.pdf");
	
	// LEFT btag vs right btag
	
	all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Left1DY");
	ss = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DY");
	
	TCanvas* c0 = new TCanvas("c0","c0",800,600);
	
	c0->cd();
	
	gPad->SetGrid();
	gPad->SetLogy();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("TCHE discriminant");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->Draw("");
	
	c0->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c0->Draw();
	
	ss->SetLineColor(kRed);
	ss->SetMarkerColor(kRed);
	ss->Draw("same");
	
	leg = new TLegend(0.45,0.25,0.8,0.4);
	
	leg->AddEntry(all,"b enriched","lp");
	leg->AddEntry(ss,"b depleted","lp");
	
	leg->Draw();
    
    text->Draw();
	
	c0->SaveAs("left-right-btag.pdf");
	
	// LEFT-RIGHT
	
	TH1F* leftrightrew = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_LeftRight1DX");
	
	TCanvas* c4 = new TCanvas("c4","c4",800,600);
	
	c4->cd();
	
	gPad->SetGrid();
	
	leftrightrew->SetTitle("");
	leftrightrew->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	leftrightrew->GetYaxis()->SetTitle("weight");
	leftrightrew->Draw();
	
	c4->Update();
	
	sb = (TPaveStats*)leftrightrew->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	//sb->SetOptStat(1111111);
	
	c4->Draw();
	
    text->Draw();

	c4->SaveAs("pt.pdf");
	
	// CS VS SS no rew
	
	TH1F* all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0");
	TH1F* ss = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar");
	
	TCanvas* c5 = new TCanvas("c5","c5",800,600);
	
	c5->cd();
	
	gPad->SetGrid();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("m_{#mu j} (GeV/c^{2})");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->Scale(1./all->Integral());
	
	all->GetYaxis()->SetRangeUser(0,0.1);
	
	all->Draw("hist");
	
	c5->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c5->Draw();
	
	ss->SetLineColor(kRed);
	ss->SetMarkerColor(kRed);
	ss->DrawNormalized("same");
	
	TLegend* leg = new TLegend(0.45,0.65,0.8,0.8);
	
	leg->AddEntry(all,"signal sample non-b","l");
	leg->AddEntry(ss,"control sample","p");
	
	leg->Draw();
    
    text->Draw();
	
	c5->SaveAs("cs-before.pdf");
	
	// CS VS SS no rew
	
	TH1F* all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0");
	TH1F* ss = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVarReweigh");
	
	TCanvas* c6 = new TCanvas("c6","c6",800,600);
	
	c6->cd();
	
	gPad->SetGrid();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("m_{#mu j} (GeV/c^{2})");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->Scale(1./all->Integral());
	
	all->GetYaxis()->SetRangeUser(0,0.1);
	
	all->Draw("hist");
	
	c6->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	c6->Draw();
	
	ss->SetLineColor(kRed);
	ss->SetMarkerColor(kRed);
	ss->DrawNormalized("same");
	
	TLegend* leg = new TLegend(0.45,0.65,0.8,0.8);
	
	leg->AddEntry(all,"signal sample non-b","l");
	leg->AddEntry(ss,"control sample","p");
	
	leg->Draw();
	
    text->Draw();

	c6->SaveAs("cs-after.pdf");
	
	// LEFT btag vs right btag
	
	all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Left1DY");
	
	TCanvas* c10 = new TCanvas("c10","c10",800,600);
	
	c10->cd();
	
	gPad->SetGrid();
	gPad->SetLogy();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("b-tagging discriminant value ");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->Draw("hist");
	
	c10->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	text->Draw();

	c10->Draw();
	
	c10->SaveAs("btag-left.root");
	c10->SaveAs("btag-left.pdf");
	c10->SaveAs("btag-left.png");
	
	// LEFT btag vs right btag
	
	all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DY")->Clone();
	
	TCanvas* c11 = new TCanvas("c11","c11",800,600);
	
	c11->cd();
	
	gPad->SetGrid();
	gPad->SetLogy();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("b-tagging discriminant value ");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->SetLineColor(kBlack);
	
	all->Draw("hist");
	
	c11->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	text->Draw();

	c11->Draw();
	
	c11->SaveAs("btag-right.root");
	c11->SaveAs("btag-right.pdf");
	c11->SaveAs("btag-right.png");
	
	// LEFT btag vs right btag
	
	all = (TH1F*) f->Get("Variable_1_0/Discriminator_0/nDisc_0_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR")->Clone();
	
	TCanvas* c12 = new TCanvas("c12","c12",800,600);
	
	c12->cd();
	
	gPad->SetGrid();
	gPad->SetLogy();
	
	all->SetTitle("");
	all->GetXaxis()->SetTitle("b-tagging discriminant value ");
	all->GetYaxis()->SetTitle("a.u.");
	
	all->SetLineColor(kBlack);
	
	all->Draw("hist");
	
	c12->Update();
	
	sb = (TPaveStats*)all->GetListOfFunctions()->FindObject("stats");
	sb->SetTextColor(kWhite);
	sb->SetLineColor(kWhite);
	sb->SetX1NDC(.99);
	sb->SetX2NDC(.99);
	sb->SetY1NDC(.99);
	sb->SetY2NDC(.99);
	
	text->Draw();

	c12->Draw();
	
	c12->SaveAs("btag-true.root");
	c12->SaveAs("btag-true.pdf");
	c12->SaveAs("btag-true.png");
	
	
	
}
