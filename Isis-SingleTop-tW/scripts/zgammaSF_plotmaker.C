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

void zgammaSF_plotmaker(int mode = 0, int region = 0){
  
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
  char myRootFile2[300];
  double lumi = 0; 
   if      (mode == 0)	 lumi = 11966.617;  	
    else if (mode == 1) lumi = 12067.294;  	
    else if (mode == 2) lumi = 12093.792; 
  

    
  sprintf(myRootFile,"results/an_%dpb_%d.root", lumi, mode); // take output from looper
  sprintf(myRootFile2,"results/noZjetsSF_%dpb_%d.root", lumi, mode); // take output from looper
 
  
  TFile *_file0 = TFile::Open(myRootFile);
  cout << myRootFile << endl;
  TFile *_file2 = TFile::Open(myRootFile2);
  cout << myRootFile2 << endl;
  
  
  const int nProcess = 1;
  const int nPlots = 4;
  TString processName[nProcess] =  { "zjets"};
  TString processTitle[nProcess] = { "Z+jets"};
  Color_t color[2] =        { kBlue,kRed};
  
  
  TString cutLabel[nPlots] =     {  "met_outside", "mll_outside", "met_zgamma", "mll_zgamma"};
  int rebinHisto[nPlots] =       {  2,2,2,2};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} outside Z mass window", "Inv. Mass outside Z mass window", "E_{T}^{miss} in Z mass window", "Inv. Mass in Z mass window"};
 

  TString modeString[3] = {"0","1","2"};
 
  TString plotExtension = "zgammaSF_plot_"; // name of the plots
  
  

  
   
  
  
  //Make plots   
  TH1F* histo_zjets;
  TH1F* histo_zjets_NOSF;
   

   
  
   for(int iPlots = 0; iPlots< nPlots; iPlots++)
   {
      
     leg = new TLegend(0.7,0.7,0.94,0.94);
     leg ->SetFillStyle(1001);
     leg ->SetFillColor(kWhite);
     leg ->SetBorderSize(1);
     
     
 
	
	
     histo_zjets= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[0]);
	       cout << cutLabel[iPlots]+ "_" + processName[0] << endl; 
	       
	       
       histo_zjets->Rebin(rebinHisto[iPlots]);
               //histo_zjets->SetFillColor(color[1]);
        histo_zjets->SetLineColor(color[0]);
        histo_zjets->SetLineWidth(2);


       histo_zjetsNOSF= (TH1F*) _file2->Get(cutLabel[iPlots]+ "_" + processName[0]);
	       cout << cutLabel[iPlots]+ "_" + processName[0] << endl; 
	       
	       
       histo_zjetsNOSF->Rebin(rebinHisto[iPlots]);
               //histo_zjetsNOSF->SetFillColor(color[1]);
        histo_zjetsNOSF->SetLineColor(color[1]);
        histo_zjetsNOSF->SetLineWidth(2); 
	      
       leg ->AddEntry(histo_zjets, "Z/#gamma * events" , "l");  
       leg ->AddEntry(histo_zjetsNOSF, "Z/#gamma * events - no reweighing" , "l"); 
       
      double max = TMath::Max(histo_zjets->GetMaximum(),histo_zjetsNOSF->GetMaximum());
    
       TCanvas *c1 = new TCanvas();
       histo_zjets->Draw("h");
      histo_zjets->SetMaximum(max * 1.2);
     //  histo_zjets->SetMinimum(1);
       histo_zjets->GetXaxis()->SetTitle(cutTitle[iPlots]);
     
       
         histo_zjetsNOSF->Draw("h,sames");
      histo_zjetsNOSF->SetMaximum(max * 1.2);
     //  histo_zjetsNOSF->SetMinimum(1);
       histo_zjetsNOSF->GetXaxis()->SetTitle(cutTitle[iPlots]);

       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1->SaveAs("plots/"  + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + ".png");
     // c1->SaveAs("plots/pdf/" + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");
      
     
      
    } // end plots loop 
    
 
} // end constructor loop
