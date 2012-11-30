//Rebeca Gonzalez Suarez
//rebeca@cern.ch

#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "setTDRStyle.C"
using namespace std;

void dis_plotmaker(int mode = 0){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  gROOT->SetBatch(1);   // to show histos on canvas
  
  

  
  char myRootFile[300];
  double lumi = 4399;

  
  
  sprintf(myRootFile,"results/dis_an_%dpb_%d.root", lumi, mode);
  
  TFile *_file0 = TFile::Open(myRootFile);
  cout << myRootFile << endl;
  
  const int nProcess = 2;
  const int nPlots = 12;
  TString processName[nProcess] =  { "_twdr", "_tt"};
  TString processTitle[nProcess] = { "tW","t#bar{t}"};
  Color_t color[nProcess] =        { kMagenta-10, kRed+1};
  

  TString cutLabel[nPlots] =     {  "promet","met", "mll", "njets", "ptsys", "ht", "pt_leading", "eta_leading", "nvertex_2lep", "pt_max", "pt_min", "pt_all"  };
  int rebinHisto[nPlots] =       { 4, 4, 4, 1, 4, 12, 4, 4,1, 2, 2,4};
  TString cutTitle[nPlots] =     {"promet", "E_{T}^{miss}", "Inv. Mass", "# of jets",  "P_{T} system [GeV]", "H_{T} [GeV]","P_{T} of the leading jet", "eta of the leading jet","# ofvertex", "p_T of the first lepton [GeV]", "p_T  of the second lepton [GeV]", "p_T  of all jets[GeV]"};


  TString modeString[1] = {"0"};
  
  TString plotExtension = "plot_";
  

  
   //////////////////////////////////////////////////////// pt_all ///////////////////////////////////////////////////////:
  pt_all_twdr->SetLineColor(kBlue);
  pt_all_tt->SetLineColor(kRed);
  
  pt_all_twdr->SetLineWidth(2);
  pt_all_tt->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  pt_all_twdr->Rebin(4);
  pt_all_tt->Rebin(4);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(pt_all_twdr, "tW signal", "l");  
  leg ->AddEntry(pt_all_tt, "tt bckgr.", "l");

  TCanvas *c40 = new TCanvas();


  
  pt_all_tt->Draw("h");
  pt_all_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c40->SaveAs("plots/compare_tw_tt_pt_all.png");


  TCanvas *c0 = new TCanvas();
  // This is only for visualization
  pt_all_twdr->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  pt_all_tt->SetNormFactor(1);

  
  pt_all_twdr->Draw("h");
  pt_all_tt->Draw("h, sames");
  
  leg->Draw();
  
  c0->SaveAs("plots/compare_tw_tt_pt_all_norm.png");
  
  
  
  //////////////////////////////////////////////////////// nvertex_2lep ///////////////////////////////////////////////////////:
  nvertex_2lep_tt->SetLineColor(kBlue);
  nvertex_2lep_twdr->SetLineColor(kRed);
  
  nvertex_2lep_tt->SetLineWidth(2);
  nvertex_2lep_twdr->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  nvertex_2lep_tt->Rebin(1);
  nvertex_2lep_twdr->Rebin(1);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(nvertex_2lep_tt, "tW signal", "l");  
  leg ->AddEntry(nvertex_2lep_twdr, "twdr bckgr.", "l");

  TCanvas *c49 = new TCanvas();


  
  nvertex_2lep_tt->Draw("h");
  nvertex_2lep_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c49->SaveAs("plots/compare_tw_twdr_nvertex_2lep.png");


  TCanvas *c9 = new TCanvas();
  // This is only for visualization
  nvertex_2lep_tt->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  nvertex_2lep_twdr->SetNormFactor(1);

  
  nvertex_2lep_tt->Draw("h");
  nvertex_2lep_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c9->SaveAs("plots/compare_tw_twdr_nvertex_2lep_norm.png");
  
  
  //////////////////////////////////////////////////////// ht ///////////////////////////////////////////////////////:
  ht_tt->SetLineColor(kBlue);
  ht_twdr->SetLineColor(kRed);
  
  ht_tt->SetLineWidth(2);
  ht_twdr->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  ht_tt->Rebin(12);
  ht_twdr->Rebin(12);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(ht_tt, "tW signal", "l");  
  leg ->AddEntry(ht_twdr, "twdr bckgr.", "l");

  TCanvas *c48 = new TCanvas();


  
  ht_tt->Draw("h");
  ht_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c48->SaveAs("plots/compare_tw_twdr_ht.png");


  TCanvas *c8 = new TCanvas();
  // This is only for visualization
  ht_tt->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  ht_twdr->SetNormFactor(1);

  
  ht_tt->Draw("h");
  ht_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c8->SaveAs("plots/compare_tw_twdr_ht_norm.png");
  
  
  
  //////////////////////////////////////////////////////// njets ///////////////////////////////////////////////////////:
  njets_tt->SetLineColor(kBlue);
  njets_twdr->SetLineColor(kRed);
  
  njets_tt->SetLineWidth(2);
  njets_twdr->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  njets_tt->Rebin(1);
  njets_twdr->Rebin(1);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(njets_tt, "tW signal", "l");  
  leg ->AddEntry(njets_twdr, "twdr bckgr.", "l");

  TCanvas *c47 = new TCanvas();


  
  njets_tt->Draw("h");
  njets_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c47->SaveAs("plots/compare_tw_twdr_njets.png");


  TCanvas *c7 = new TCanvas();
  // This is only for visualization
  njets_tt->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  njets_twdr->SetNormFactor(1);

  
  njets_tt->Draw("h");
  njets_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c7->SaveAs("plots/compare_tw_twdr_njets_norm.png");
  
  
  
  //////////////////////////////////////////////////////// eta_leading ///////////////////////////////////////////////////////:
  eta_leading_tt->SetLineColor(kBlue);
  eta_leading_twdr->SetLineColor(kRed);
  
  eta_leading_tt->SetLineWidth(2);
  eta_leading_twdr->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  eta_leading_tt->Rebin(4);
  eta_leading_twdr->Rebin(4);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(eta_leading_tt, "tW signal", "l");  
  leg ->AddEntry(eta_leading_twdr, "twdr bckgr.", "l");

  TCanvas *c46 = new TCanvas();


  
  eta_leading_tt->Draw("h");
  eta_leading_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c46->SaveAs("plots/compare_tw_twdr_eta_leading.png");


  TCanvas *c6 = new TCanvas();
  // This is only for visualization
  eta_leading_tt->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  eta_leading_twdr->SetNormFactor(1);

  
  eta_leading_tt->Draw("h");
  eta_leading_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c6->SaveAs("plots/compare_tw_twdr_eta_leading_norm.png");
  
  //////////////////////////////////////////////////////// pt_leading ///////////////////////////////////////////////////////:
  pt_leading_tt->SetLineColor(kBlue);
  pt_leading_twdr->SetLineColor(kRed);
  
  pt_leading_tt->SetLineWidth(2);
  pt_leading_twdr->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  pt_leading_tt->Rebin(4);
  pt_leading_twdr->Rebin(4);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(pt_leading_tt, "tW signal", "l");  
  leg ->AddEntry(pt_leading_twdr, "twdr bckgr.", "l");

  TCanvas *c45 = new TCanvas();


  
  pt_leading_tt->Draw("h");
  pt_leading_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c45->SaveAs("plots/compare_tw_twdr_pt_leading.png");


  TCanvas *c5 = new TCanvas();
  // This is only for visualization
  pt_leading_tt->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  pt_leading_twdr->SetNormFactor(1);

  
  pt_leading_tt->Draw("h");
  pt_leading_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c5->SaveAs("plots/compare_tw_twdr_pt_leading_norm.png");
  
  
  
  //////////////////////////////////////////////////////// ptsys ///////////////////////////////////////////////////////:
  ptsys_tt->SetLineColor(kBlue);
  ptsys_twdr->SetLineColor(kRed);
  
  ptsys_tt->SetLineWidth(2);
  ptsys_twdr->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  ptsys_tt->Rebin(4);
  ptsys_twdr->Rebin(4);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(ptsys_tt, "tW signal", "l");  
  leg ->AddEntry(ptsys_twdr, "twdr bckgr.", "l");

  TCanvas *c44 = new TCanvas();


  
  ptsys_tt->Draw("h");
  ptsys_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c44->SaveAs("plots/compare_tw_twdr_ptsys.png");


  TCanvas *c4 = new TCanvas();
  // This is only for visualization
  ptsys_tt->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  ptsys_twdr->SetNormFactor(1);

  
  ptsys_tt->Draw("h");
  ptsys_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c4->SaveAs("plots/compare_tw_twdr_ptsys_norm.png");
  
  //////////////////////////////////////////////////////// mll ///////////////////////////////////////////////////////:
  mll_tt->SetLineColor(kBlue);
  mll_twdr->SetLineColor(kRed);
  
  mll_tt->SetLineWidth(2);
  mll_twdr->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  mll_tt->Rebin(4);
  mll_twdr->Rebin(4);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(mll_tt, "tW signal", "l");  
  leg ->AddEntry(mll_twdr, "twdr bckgr.", "l");

  TCanvas *c43 = new TCanvas();


  
  mll_tt->Draw("h");
  mll_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c43->SaveAs("plots/compare_tw_twdr_mll.png");


  TCanvas *c3 = new TCanvas();
  // This is only for visualization
  mll_tt->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  mll_twdr->SetNormFactor(1);

  
  mll_tt->Draw("h");
  mll_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c3->SaveAs("plots/compare_tw_twdr_mll_norm.png");
  
  
  //////////////////////////////////////////////////////// promet ///////////////////////////////////////////////////////:
  promet_tt->SetLineColor(kBlue);
  promet_twdr->SetLineColor(kRed);
  
  promet_tt->SetLineWidth(2);
  promet_twdr->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  promet_tt->Rebin(4);
  promet_twdr->Rebin(4);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(promet_tt, "tW signal", "l");  
  leg ->AddEntry(promet_twdr, "twdr bckgr.", "l");

  TCanvas *c42 = new TCanvas();


  
  promet_tt->Draw("h");
  promet_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c42->SaveAs("plots/compare_tw_twdr_promet.png");


  TCanvas *c2 = new TCanvas();
  // This is only for visualization
  promet_tt->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  promet_twdr->SetNormFactor(1);

  
  promet_tt->Draw("h");
  promet_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c2->SaveAs("plots/compare_tw_twdr_promet_norm.png");
  
  
  //////////////////////////////////////////////////////// met ///////////////////////////////////////////////////////:
  met_tt->SetLineColor(kBlue);
  met_twdr->SetLineColor(kRed);
  
  met_tt->SetLineWidth(2);
  met_twdr->SetLineWidth(2);

 // Rebin  // om aantal bins te kiezen 
  met_tt->Rebin(4);
  met_twdr->Rebin(4);



 
  
  leg = new TLegend(0.7,0.7,0.94,0.94);  // maak legende
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(met_tt, "tW signal", "l");  
  leg ->AddEntry(met_twdr, "twdr bckgr.", "l");

  TCanvas *c41 = new TCanvas();


  
  met_tt->Draw("h");
  met_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c41->SaveAs("plots/compare_tw_twdr_met.png");


  TCanvas *c1 = new TCanvas();
  // This is only for visualization
  met_tt->SetNormFactor(1);   // Zet alles genormaliseerd tot 1
  met_twdr->SetNormFactor(1);

  
  met_tt->Draw("h");
  met_twdr->Draw("h, sames");
  
  leg->Draw();
  
  c1->SaveAs("plots/compare_tw_twdr_met_norm.png");
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  
  /*
  //Make plots   HOW MAKE THIS WORK
   for(int iPlots = 0; iPlots< nPlots; iPlots++)
     TString Variable = cutLabel[iPlots];
    	cout << Variable << endl; 
	
       int rebins = rebinHisto[iPlots];
	
        for(int iProcess = 0; iProcess < nProcess; iProcess++)
	{
	       TString ProcessName = processName[iProcess];
    		cout << ProcessName << endl;
	       TString Naam = Variable + ProcessName ;
	       cout << Naam << endl;
	       Color_t Kleur = color[iProcess];
	       
		
		Naam->SetLineColor(Kleur); 
		Naam->SetLineWidth(2); 
		Naam->Rebin(rebins); 
		Naam->SetNormFactor(1); 
	       
	}
	
	TCanvas *c1 = new TCanvas();
	
        for(int iProcess = 0; iProcess < nProcess; nProcess++)
	{
    		cutLabel[iPlots]_processName[iProcess]->Draw("h");  
	
	}
        leg->Draw();
        c1->SaveAs("plots/" + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + "_discrimination" +".png");
      
    
    }
    
    */
    
    
    
    
    

    
  
  
}
