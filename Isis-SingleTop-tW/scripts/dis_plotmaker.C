//Isis Van Parijs


#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "setTDRStyle.C"
//#include "../Tools/interface/MultiSamplePlot.h"
//#include "../Tools/interface/PlottingTools.h"
using namespace std;

void dis_plotmaker(int mode = 0, int region = 0){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  gROOT->SetBatch(1);   // to show histos on canvas
  
  //Making it pretty
  labelcms  = new TPaveText(0.12,0.88,0.5,0.94,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.045);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS Preliminary, #sqrt{s} = 8 TeV");
  labelcms->SetBorderSize(0);

  labelcms2  = new TPaveText(0.12,0.85,0.5,0.88,"NDCBR");
  labelcms2->SetTextAlign(12);
  labelcms2->SetTextSize(0.045);
  labelcms2->SetFillColor(kWhite);
     if      (mode == 0) labelcms2->AddText("12 fb^{-1}, e#mu channel  ");	
    else if (mode == 1) labelcms2->AddText("12 fb^{-1}, #mu#mu channel  "); 	
    else if (mode == 2) labelcms2->AddText("12 fb^{-1}, ee channel  "); 
  
  labelcms2->SetBorderSize(0);
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(600);
  gStyle->SetLabelFont(18,"");
  
  gStyle->SetTitleXOffset(1.2);//1.5
  gStyle->SetTitleYOffset(1.2);//1.7
  
  char myRootFile[300];
  double lumi = 0; 
   if      (mode == 0)	 lumi = 11966.617;  	
    else if (mode == 1) lumi = 12067.294;  	
    else if (mode == 2) lumi = 12093.792; 
  

    
  sprintf(myRootFile,"results/an_%dpb_%d.root", lumi, mode); // take output from looper
  
  TFile *_file0 = TFile::Open(myRootFile);
  cout << myRootFile << endl;
  
  const int nProcess = 2;
  const int nPlots = 10;
  TString processName[nProcess] =  { "twdr", "tt"};
  TString processTitle[nProcess] = { "tW","t#bar{t}"};
  Color_t color[nProcess] =        { kBlue, kRed};
  
  if(region == 0){
  TString cutLabel[nPlots] =     {  "met_1j1t", "mll_1j1t", "ptsys_1j1t", "ht_1j1t","pt_leading_1j1t",  "nvertex_1j1t", "pt_max_1j1t","pt_min_1j1t","ht_nomet_1j1t","eta_leading_1j1t" };
  int rebinHisto[nPlots] =       {  4, 4,  4, 12, 4, 1,4,4,12,4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} 1j1t", "Inv. Mass 1j1t", "P_{T} system [GeV] 1j1t","H_{T} [GeV] 1j1t","P_{T} of the leading jet 1j1t", "# of vertex 1j1t", "p_T of the first lepton [GeV] 1j1t", "p_T  of the second lepton [GeV] 1j1t", "H_{T} no met [GeV] 1j1t","Eta of the leading jet 1j1t"};
  } else if(region == 1){
   TString cutLabel[nPlots] =     {  "met_2j1t", "mll_2j1t", "ptsys_2j1t", "ht_2j1t","pt_leading_2j1t",  "nvertex_2j1t", "pt_max_2j1t","pt_min_2j1t","ht_nomet_2j1t","eta_leading_2j1t" };
  int rebinHisto[nPlots] =       {  4, 4,  4, 12, 4, 1,4,4,12,4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} 2j1t", "Inv. Mass 2j1t", "P_{T} system [GeV] 2j1t","H_{T} [GeV] 2j1t","P_{T} of the leading jet 2j1t", "# of vertex 2j1t", "p_T of the first lepton [GeV] 2j1t", "p_T  of the second lepton [GeV] 2j1t", "H_{T} no met [GeV] 2j1t","Eta of the leading jet 2j1t"};
 
  
  
  }else if(region == 2){
   TString cutLabel[nPlots] =     {  "met_2j2t", "mll_2j2t", "ptsys_2j2t", "ht_2j2t","pt_leading_2j2t",  "nvertex_2j2t", "pt_max_2j2t","pt_min_2j2t","ht_nomet_2j2t","eta_leading_2j2t" };
  int rebinHisto[nPlots] =       {  4, 4,  4, 12, 4, 1,4,4,12,4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} 2j2t", "Inv. Mass 2j2t", "P_{T} system [GeV] 2j2t","H_{T} [GeV] 2j2t","P_{T} of the leading jet 2j2t", "# of vertex 2j2t", "p_T of the first lepton [GeV] 2j2t", "p_T  of the second lepton [GeV] 2j2t", "H_{T} no met [GeV] 2j2t","Eta of the leading jet 2j2t"};
 
  
  }

  TString modeString[3] = {"0","1","2"};
  TString regionString[3] = {"0","1","2"};
  
  TString plotExtension = "discriminatingVar_plot_"; // name of the plots
  TString plotAnalysis = "DiscriminatingVar"; // directory in plots where the plots are saved
  

  
   
  
  
  //Make plots   
  TH1F* histo_tt;
  TH1F* histo_twdr;
   

   
  
   for(int iPlots = 0; iPlots< nPlots; iPlots++)
   {
      
     leg = new TLegend(0.7,0.7,0.94,0.94);
     leg ->SetFillStyle(1001);
     leg ->SetFillColor(kWhite);
     leg ->SetBorderSize(1);
     
     
 
	
	
     histo_tt= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[1]);
	       cout << cutLabel[iPlots]+ "_" + processName[1] << endl; 
	       
	       
       histo_tt->Rebin(rebinHisto[iPlots]);
               //histo_tt->SetFillColor(color[1]);
        histo_tt->SetLineColor(color[1]);
        histo_tt->SetLineWidth(2);


       histo_twdr= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[0]);
        cout << cutLabel[iPlots]+ "_" + processName[0] << endl; 

       histo_twdr->Rebin(rebinHisto[iPlots]);
              // histo_twdr->SetFillColor(color[0]);
       histo_twdr->SetLineColor(color[0]);
       histo_twdr->SetLineWidth(2);
	      
       leg ->AddEntry(histo_tt, "tt bckgr." , "l");  
       leg ->AddEntry(histo_twdr, "twdr signal", "l");  
	
       double max = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
    
       TCanvas *c1 = new TCanvas();
       
       histo_tt->SetMaximum(max * 1.2);
       histo_tt->SetMinimum(1);
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->Draw("h");
       
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max * 1.2);
       histo_twdr->SetMinimum(1);

       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + ".png");
     // c1->SaveAs("plots/pdf/" + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");
      
      int y = 1/max; 
      
      
      TCanvas *c2 = new TCanvas();
       
       histo_tt->SetMaximum(max * 1.7);
       //histo_tt->SetMinimum(1);
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->GetYaxis()->SetTitle("events / 12 fb^{-1}");
    histo_tt->DrawNormalized("h",1);
       
       
       histo_twdr->DrawNormalized("h,sames",1);
      histo_twdr->SetMaximum(max * 1.7);
     //  histo_twdr->SetMinimum(1);

       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c2->SaveAs("plots/"  + plotAnalysis +"/"  + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + "_normalized" +  ".png");
    //  c2->SaveAs("plots/pdf/" + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + "_normalized" + ".pdf");
      

      
    } // end plots loop 
    
 
} // end constructor loop
