{
    TFile* f = new TFile("../store/data.root");
    TFile* g = new TFile("../store/mc.root");
    
    TH1F* effMC = (TH1F*) g->Get("Variable_1_0/Discriminator_6/nDisc_6_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffAll");
    TH1F* eff = (TH1F*) f->Get("Variable_1_0/Discriminator_6/nDisc_6_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredRR");

    Int_t nBin = eff->GetNbinsX();
    Float_t low = eff->GetBinLowEdge(1); 
    Float_t high = eff->GetBinLowEdge(eff->GetNbinsX()+1);
        
    double gr_MC[10000];
    double gr_x[10000];
    double gr_y[10000];
    
    double gr_ex[10000];
    double gr_ey[10000];
    
    double gr_ratio_y[10000];
    double gr_ratio_y_err_up[10000];
    double gr_ratio_y_err_down[10000];
    double gr_ratio_x[10000];
    
    double maxup = 0;
    
    double maxdown = 100;
    
    for (unsigned int i=0; i<nBin; i++) {
        
        gr_x[i] = eff->GetBinCenter(i+1);
        gr_y[i] = eff->GetBinContent(i+1);
        
        gr_ex[i] = 0;
        gr_ey[i] = eff->GetBinError(i+1);
        
        gr_MC[i] = effMC->GetBinContent(i+1);
        
        gr_ratio_x[i] = eff->GetBinCenter(i+1);
        gr_ratio_y[i] = gr_y[i]/gr_MC[i];
        ///cout << gr_ratio_y[i] << endl;
        
        float a = eff->GetBinCenter(i+1);
        float b = effMC->GetBinCenter(i+1);
        float ua = eff->GetBinError(i+1);
        float ub = effMC->GetBinError(i+1);
        
        float uSF = sqrt((pow(ua,2)/pow(b,2))+((pow(a,2)*pow(ub,2))/pow(b,4)));
        gr_ratio_y_err_up[i] =  gr_ratio_y[i]+sqrt(pow(0.048,2)+pow(uSF,2));
        gr_ratio_y_err_down[i] =  gr_ratio_y[i]-sqrt(pow(0.048,2)+pow(uSF,2));
        
        if (gr_ratio_y_err_up[i] > maxup) maxup=gr_ratio_y_err_up[i];
        if (gr_ratio_y_err_down[i] < maxdown) maxdown=gr_ratio_y_err_down[i];
        
    }
    
    cout << maxup << " " << maxdown << endl;
    
    TCanvas* canvas =  new TCanvas( "canvas", "canvas", 1000,700);
    canvas->SetBottomMargin(0.3);
    
    TMultiGraph* multi_graph = new TMultiGraph();
    TMultiGraph* multi_graph_ratio = new TMultiGraph();
    
    TGraph* graph_MC = new TGraph( nBin, gr_x, gr_MC );
    TGraphErrors* graph_nom = new TGraphErrors( nBin, gr_x, gr_y, gr_ex, gr_ey );
    TGraph* graph_ratio = new TGraph( nBin, gr_ratio_x, gr_ratio_y );
    TGraph* graph_ratio_up = new TGraph( nBin, gr_ratio_x, gr_ratio_y_err_up );
    TGraph* graph_ratio_down = new TGraph( nBin, gr_ratio_x, gr_ratio_y_err_down );
    
    graph_nom->SetMarkerStyle(20);
    graph_nom->SetMarkerColor(2);
    graph_MC->SetLineWidth(3);
    
    multi_graph->Add(graph_MC, "C");
    multi_graph->Add(graph_nom, "p");
    
    multi_graph->SetMaximum(1.1);
    multi_graph.SetMinimum(0.05);
    
    multi_graph->Draw("A");
    
    multi_graph->GetYaxis()->SetTitle("b-tag Efficiency");
    multi_graph->GetYaxis()->SetTitleOffset(0.85);
    multi_graph->GetXaxis()->SetLabelSize(0);
    multi_graph->GetXaxis()->SetTitleSize(0);
    
    float diff = eff->GetBinCenter(2)-eff->GetBinCenter(1);
    
    //multi_graph->GetXaxis()->SetLimits(eff->GetBinCenter(2) - diff , eff->GetBinCenter(eff->GetNbinsX()) + diff);
    
    TLegend* leg= new TLegend(0.156878,0.389717,0.534672,0.578527);
    
    leg->AddEntry(graph_MC,"Monte Carlo Truth","l");
    leg->AddEntry(graph_nom,"Measured Value","lp");
    
    
    leg->SetShadowColor(0);
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->Draw();
    
    latex = TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.05);
    latex->SetTextAlign(29);// align right
    latex->DrawLatex(0.15, 0.86, "CMS Preliminary");
    
    
    latex2 = TLatex();
    latex2->SetNDC();
    latex2->SetTextSize(0.05);
    latex2->SetTextAlign(30); // align right
    latex2->DrawLatex(0.87, 0.86, "2.18 fb^{-1} at #sqrt{s} = 7 TeV");
    
    TPad* pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
    pad->SetTopMargin(0.7);
    pad->SetFillColor(0);
    pad->SetFillStyle(0);
    pad->Draw();
    pad->cd(0);
    
    graph_ratio->SetMarkerStyle(20);
    graph_ratio->SetMaximum( 1.02*maxup);
    graph_ratio->SetMinimum( 0.9*maxdown);
    
    graph_ratio->GetXaxis()->SetTitle("Discriminator Value");
    graph_ratio->GetYaxis()->SetTitle("SF_{b}");
    graph_ratio->GetYaxis()->SetTitleOffset(0.85);
    //->GetYaxis()->SetNdivisions(5);
    
    //graph_ratio->GetXaxis()->SetLimits(discrim[0] - diff , discrim[len(discrim)-1] + diff)
    
    TF1* myfit_quad = new TF1("myfit_quad","[0]*x*x*x + [1]*x*x + [2]*x + [3]", low - 0.3*diff, high + 0.3*diff);
    graph_ratio_up -> Fit ("myfit_quad");
    myfit_quad->SetLineColor(4);
    myfit_quad->SetLineStyle(7);
    
    TF1* myfit_quad2 = new TF1("myfit_quad2","[0]*x*x*x + [1]*x*x + [2]*x + [3]", low - 0.3*diff, high + 0.3*diff);
    graph_ratio_down -> Fit ("myfit_quad2");
    myfit_quad2->SetLineColor(4);
    myfit_quad2->SetLineStyle(7);
    
    
    multi_graph_ratio->Add(graph_ratio, "p");  
    multi_graph_ratio->Draw("A");
    
    multi_graph_ratio->SetMaximum( 1.02*maxup);
    multi_graph_ratio->SetMinimum( 0.9*maxdown);
    
    multi_graph_ratio.GetXaxis()->SetTitle("Discriminator Value");
    multi_graph_ratio.GetYaxis()->SetTitle("SF_{b}");
    multi_graph_ratio.GetYaxis()->SetTitleOffset(0.85);
    multi_graph_ratio.GetYaxis()->SetNdivisions(5);
    //multi_graph_ratio.GetXaxis().SetLimits(discrim[0] - diff , discrim[len(discrim)-1] + diff)
        
    myfit_quad->Draw("SAME");
    myfit_quad2->Draw("SAME");
    
    //for CSV        
    Float_t line_x1_L = 0.24;
    Float_t line_x1_M = 0.68;
    Float_t line_x1_T = 0.9;
    Float_t line_y1 = 0.9*maxdown;
    Float_t line_y2 = (0.9*maxdown)+0.1;
    
    TArrow* line_L_1 = new TArrow(line_x1_L, line_y1, line_x1_L, line_y2, 0.02, "<|");
    line_L_1.SetLineColor(2);
    TArrow* line_L_2 = new TArrow(line_x1_M, line_y1, line_x1_M, line_y2, 0.02, "<|");
    line_L_2.SetLineColor(2);
    TArrow* line_L_3 = new TArrow(line_x1_T, line_y1, line_x1_T, line_y2, 0.02, "<|");
    line_L_3.SetLineColor(2);
    
    line_L_1->Draw("SAME <|");
    line_L_2->Draw("SAME <|");
    line_L_3->Draw("SAME <|");
    
    canvas->SaveAs("Plot.pdf");    
}