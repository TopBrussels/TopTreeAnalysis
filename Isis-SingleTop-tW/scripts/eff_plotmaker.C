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

void eff_plotmaker(int mode = 0){
  
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
  labelcms->AddText("CMS simulation, #sqrt{s} = 8 TeV");
  labelcms->SetBorderSize(0);

  labelcms2  = new TPaveText(0.12,0.85,0.5,0.88,"NDCBR");
  labelcms2->SetTextAlign(12);
  labelcms2->SetTextSize(0.045);
  labelcms2->SetFillColor(kWhite);
  labelcms2->AddText("e#mu channel  ");
  labelcms2->SetBorderSize(0);
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(600);
  gStyle->SetLabelFont(18,"");
  
  gStyle->SetTitleXOffset(1.5);//1.5
  gStyle->SetTitleYOffset(1.7);//1.7
  
  char myRootFile[300];
  char myRootFilet[300];
  char myRootFilett[300];

 
  

    
  sprintf(myRootFile,"/user/ivanpari/MasterThesis_newEnv/TopBrussels/TopTreeAnalysis/Isis-SingleTop-tW/BtagFiles/btag_eff_%d_atw_dr.root",mode); // take output from singletop
  sprintf(myRootFilet,"/user/ivanpari/MasterThesis_newEnv/TopBrussels/TopTreeAnalysis/Isis-SingleTop-tW/BtagFiles/btag_eff_%d_tw_dr.root",mode); // take output from singletop
  sprintf(myRootFilett,"/user/ivanpari/MasterThesis_newEnv/TopBrussels/TopTreeAnalysis/Isis-SingleTop-tW/BtagFiles/btag_eff_%d_tt.root",mode);

  
  
   
  TFile *_file0 = TFile::Open(myRootFile);
  TFile *_file1 = TFile::Open(myRootFilet);
  TFile *_file2 = TFile::Open(myRootFilett);

  
  cout << "-------------------------------------------------------" << endl; 
  cout << " ------------ USED FILES ------------------------------" << endl; 
  cout << "Normal tw: " << myRootFilet << endl;
  cout << "Anti tw: " << myRootFile << endl;
  cout << "ttbar: " << myRootFilett << endl;
  
  cout << "-------------------------------------------------------" << endl; 
  
  const int nProcess = 4;
  const int nPlots = 12;
  const int rebins = 12;



  TString processName[nProcess] =  { "atw_dr", "tw_dr", "tt", "twdr"};
  Color_t color[nProcess] =        { kBlue, kRed, kGreen, kMagenta};
  
  


  TString cutLabel[nPlots] =     {  "ptBjet", "ptCjet", "ptLjet", "ptBTagB", "ptBTagC", "ptBTagL", "etaBjet", "etaCjet", "etaLjet", "etaBTagB", "etaBTagC", "etaBTagL" };
  TString String_Rebin[rebins]=  {  "ptBjet","ptBTagB", "etaBjet","etaBTagB","ptCjet",  "ptBTagC", "etaCjet", "etaBTagC","ptLjet","ptBTagL","etaLjet", "etaBTagL" };
  int rebinHisto[rebins] =       {2,2,2,2,5,5,4,4,4,4,2,2}; // 
  TString cutTitle[nPlots] =     { "pt jets coming from a b quark ", "pt jets coming from a c quark ", "pt jets coming from a L quark  ", "pt btagged jets coming from a b quark ","pt btagged jets coming from a c quark ", "pt btagged jets coming from a L quark ","eta jets coming from a b quark ", "eta jets coming from a c quark ", "eta jets coming from a L quark  ", "eta btagged jets coming from a b quark  ","eta btagged jets coming from a c quark ", "eta btagged jets coming from a L quark " };


  TString modeString[1] = {"0"};
  
  TString plotExtension = "final_eff_plot_"; // name of the plots
  
  TString plotAnalysis = "/user/ivanpari/MasterThesis_newEnv/TopBrussels/TopTreeAnalysis/Isis-SingleTop-tW/BtagFiles/plots/"; // directory in plots where the plots are saved
  
  
  

  
   
  
  
  //Make plots   
  TH1F* histo_atw;
  TH1F* histo_tw;
  TH1F* histo_tt;
  TH1F* histo_twdr;
  TH1F* histo_tt_eff_pt;
  TH1F* histo_twdr_eff_pt;
  TH1F* histo_tt_eff_eta;
  TH1F* histo_twdr_eff_eta;
  
  
  
   for(int iPlots = 0; iPlots< nPlots; iPlots++)
   {
      
     leg = new TLegend(0.7,0.7,0.94,0.94);
     leg ->SetFillStyle(1001);
     leg ->SetFillColor(kWhite);
     leg ->SetBorderSize(1);
     
 
	
       // Get histos  
       histo_atw= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw: " << cutLabel[iPlots]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw= (TH1F*) _file1->Get(cutLabel[iPlots]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw: " << cutLabel[iPlots]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	       
      histo_tt= (TH1F*) _file2->Get(cutLabel[iPlots]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt: " << cutLabel[iPlots]+ "_" + modeString[mode] + "_"+processName[2] << endl;        
 	
	       	       
      
        
      
      
       for(int k = 0; k<rebins; k++){
         if(cutLabel[iPlots].CompareTo(String_Rebin[k])==0 ){
	 	histo_tt->Rebin(rebinHisto[k]);
		cout << "Histo tt: " << cutLabel[iPlots] << " rebinned with " << rebinHisto[k] << endl; 
		histo_atw->Rebin(rebinHisto[k]);
		cout << "Histo atw: " << cutLabel[iPlots] << " rebinned with " << rebinHisto[k] << endl;
		histo_tw->Rebin(rebinHisto[k]);
		cout << "Histo atw: " << cutLabel[iPlots] << " rebinned with " << rebinHisto[k] << endl;
	 
	 }

       }
       
       histo_twdr = (TH1F*) histo_tw ->Clone("histo_twdr"); 
       histo_twdr->Add(histo_atw);
      
       histo_atw->SetLineColor(color[0]);
       histo_atw->SetLineWidth(2);
       histo_tw->SetLineColor(color[1]);
       histo_tw->SetLineWidth(2);
       histo_tt->SetLineColor(color[2]);
       histo_tt->SetLineWidth(2);
       histo_twdr->SetLineColor(color[3]);
       histo_twdr->SetLineWidth(2); 
       
      
     
       
       
	// Make legends      
       leg ->AddEntry(histo_tt, "tt bckgr." , "l");  
       leg ->AddEntry(histo_atw, "anti t + W" , "l"); 
       leg ->AddEntry(histo_tw, "t + W" , "l"); 
       leg ->AddEntry(histo_twdr, "twdr signal (atW + tW)" , "l");
     
       
        // determine maxima
        double max = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
      
       
    
       
       TCanvas *c1 = new TCanvas();
       
       histo_tt->Draw("h");
       histo_tt->SetMaximum(max*1.5 );
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
    
       histo_tw->Draw("h, sames");      
       histo_twdr->Draw("h, sames");
       histo_atw->Draw("h , sames");



  
       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1->SaveAs(plotAnalysis + plotExtension +  modeString[mode] + "_" + cutLabel[iPlots] + ".png");
     
      
      

      
    } // end plots loop 
    
    
    
    ///////////////////
   // pt efficiency 
    //////////////////
    TH1F* histo_atw_ptAllB;
    TH1F* histo_tw_ptAllB; 
    TH1F* histo_tt_ptAllB;
    TH1F* histo_atw_ptAllC;
    TH1F* histo_tw_ptAllC; 
    TH1F* histo_tt_ptAllC;
    TH1F* histo_atw_ptAllL;
    TH1F* histo_tw_ptAllL; 
    TH1F* histo_tt_ptAllL;
    TH1F* histo_atw_ptBtagB;
    TH1F* histo_tw_ptBtagB; 
    TH1F* histo_tt_ptBtagB;
    TH1F* histo_atw_ptBtagC;
    TH1F* histo_tw_ptBtagC; 
    TH1F* histo_tt_ptBtagC;
    TH1F* histo_atw_ptBtagL;
    TH1F* histo_tw_ptBtagL; 
    TH1F* histo_tt_ptBtagL;
   
    TH1F* histo_twdr_ptBtagB;
    TH1F* histo_twdr_ptBtagC; 
    TH1F* histo_twdr_ptBtagL;
    TH1F* histo_twdr_ptAllB;
    TH1F* histo_twdr_ptAllC; 
    TH1F* histo_twdr_ptAllL;
   
    TH1F* histo_twdr_ptEffB; 
    TH1F* histo_tt_ptEffB;
    TH1F* histo_twdr_ptEffC; 
    TH1F* histo_tt_ptEffC;
    TH1F* histo_twdr_ptEffL; 
    TH1F* histo_tt_ptEffL;
    
    

    
    legEffBpt = new TLegend(0.7,0.7,0.94,0.94);
     legEffBpt ->SetFillStyle(1001);
     legEffBpt ->SetFillColor(kWhite);
     legEffBpt ->SetBorderSize(1);
     
    legEffCpt = new TLegend(0.7,0.7,0.94,0.94);
     legEffCpt ->SetFillStyle(1001);
     legEffCpt ->SetFillColor(kWhite);
     legEffCpt ->SetBorderSize(1);
     
     legEffLpt = new TLegend(0.7,0.7,0.94,0.94);
     legEffLpt ->SetFillStyle(1001);
     legEffLpt ->SetFillColor(kWhite);
     legEffLpt ->SetBorderSize(1);
    
   
     histo_atw_ptAllB= (TH1F*) _file0->Get(cutLabel[0]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw pt all B: " << cutLabel[0]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_ptAllB= (TH1F*) _file1->Get(cutLabel[0]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw pt all B: " << cutLabel[0]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_ptAllB= (TH1F*) _file2->Get(cutLabel[0]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt pt all B: " << cutLabel[0]+ "_" + modeString[mode] + "_"+processName[2] << endl;  
	histo_atw_ptAllC= (TH1F*) _file0->Get(cutLabel[1]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw pt all C: " << cutLabel[1]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_ptAllC= (TH1F*) _file1->Get(cutLabel[1]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw pt all C: " << cutLabel[1]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_ptAllC= (TH1F*) _file2->Get(cutLabel[1]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt pt all C: " << cutLabel[1]+ "_" + modeString[mode] + "_"+processName[2] << endl;   
	histo_atw_ptAllL= (TH1F*) _file0->Get(cutLabel[2]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw pt all L: " << cutLabel[2]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_ptAllL= (TH1F*) _file1->Get(cutLabel[2]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw pt all L: " << cutLabel[2]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_ptAllL= (TH1F*) _file2->Get(cutLabel[2]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt pt all L: " << cutLabel[2]+ "_" + modeString[mode] + "_"+processName[2] << endl; 
	histo_atw_ptBtagB= (TH1F*) _file0->Get(cutLabel[3]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw pt btag B: " << cutLabel[3]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_ptBtagB= (TH1F*) _file1->Get(cutLabel[3]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw pt btag B: " << cutLabel[3]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_ptBtagB= (TH1F*) _file2->Get(cutLabel[3]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt pt btag B: " << cutLabel[3]+ "_" + modeString[mode] + "_"+processName[2] << endl;  
	histo_atw_ptBtagC= (TH1F*) _file0->Get(cutLabel[4]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw pt btag C: " << cutLabel[4]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_ptBtagC= (TH1F*) _file1->Get(cutLabel[4]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw pt btag C: " << cutLabel[4]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_ptBtagC= (TH1F*) _file2->Get(cutLabel[4]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt pt btag C: " << cutLabel[4]+ "_" + modeString[mode] + "_"+processName[2] << endl;   
	histo_atw_ptBtagL= (TH1F*) _file0->Get(cutLabel[5]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw pt btag L: " << cutLabel[5]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_ptBtagL= (TH1F*) _file1->Get(cutLabel[5]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw pt btag L: " << cutLabel[5]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_ptBtagL= (TH1F*) _file2->Get(cutLabel[5]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt pt btag L: " << cutLabel[5]+ "_" + modeString[mode] + "_"+processName[2] << endl; 
	     
	histo_twdr_ptBtagB = (TH1F*) histo_atw_ptBtagB->Clone("histo_twdr_ptBtagB");
	histo_twdr_ptBtagC = (TH1F*) histo_atw_ptBtagC->Clone("histo_twdr_ptBtagC");
	histo_twdr_ptBtagL = (TH1F*) histo_atw_ptBtagL->Clone("histo_twdr_ptBtagL");
	histo_twdr_ptAllB = (TH1F*)histo_atw_ptAllB->Clone("histo_twdr_ptAllB");
	histo_twdr_ptAllC = (TH1F*)histo_atw_ptAllC->Clone("histo_twdr_ptAllC");
	histo_twdr_ptAllL = (TH1F*)histo_atw_ptAllL->Clone("histo_twdr_ptAllL");
	histo_twdr_ptBtagB->Add(histo_tw_ptBtagB);
	histo_twdr_ptBtagC->Add(histo_tw_ptBtagC);
	histo_twdr_ptBtagL->Add(histo_tw_ptBtagL);
	histo_twdr_ptAllB->Add(histo_tw_ptAllB);
	histo_twdr_ptAllC->Add(histo_tw_ptAllC);
	histo_twdr_ptAllL->Add(histo_tw_ptAllL);
	
 
	
	histo_twdr_ptEffB = (TH1F*) histo_twdr_ptBtagB->Clone("histo_twdr_ptEffB"); 
	histo_twdr_ptEffC = (TH1F*) histo_twdr_ptBtagC->Clone("histo_twdr_ptEffC");
	histo_twdr_ptEffL = (TH1F*) histo_twdr_ptBtagL->Clone("histo_twdr_ptEffL");
	
	histo_tt_ptEffB = (TH1F*) histo_tt_ptBtagB->Clone("histo_tt_ptEffB"); 
	histo_tt_ptEffC = (TH1F*) histo_tt_ptBtagC->Clone("histo_tt_ptEffC");
	histo_tt_ptEffL = (TH1F*) histo_tt_ptBtagL->Clone("histo_tt_ptEffL");
	
	
	histo_twdr_ptEffB->Divide(histo_twdr_ptBtagB,histo_twdr_ptAllB,1.,1., "B"); 
    	histo_tt_ptEffB->Divide(histo_tt_ptBtagB,histo_tt_ptAllB,1.,1.,"B");
	 
 	histo_twdr_ptEffC->Divide(histo_twdr_ptBtagC,histo_twdr_ptAllC,1.,1., "B"); 
    	histo_tt_ptEffC->Divide(histo_tt_ptBtagC,histo_tt_ptAllC,1.,1.,"B");
	
        histo_twdr_ptEffL->Divide(histo_twdr_ptBtagL,histo_twdr_ptAllL,1.,1., "B"); 
    	histo_tt_ptEffL->Divide(histo_tt_ptBtagL,histo_tt_ptAllL,1.,1.,"B");
	
	legEffBpt ->AddEntry(histo_tt_ptEffB, "ttbar bckgr. eff" , "l");  
        legEffBpt ->AddEntry(histo_twdr_ptEffB, "twdr signal eff" , "l");
	
	legEffCpt ->AddEntry(histo_tt_ptEffC, "ttbar bckgr. eff" , "l");  
        legEffCpt ->AddEntry(histo_twdr_ptEffC, "twdr signal eff" , "l");
	
	legEffLpt ->AddEntry(histo_tt_ptEffL, "ttbar bckgr. eff" , "l");  
        legEffLpt ->AddEntry(histo_twdr_ptEffL, "twdr signal eff" , "l");
	
	double maxEffBpt = TMath::Max(histo_tt_ptEffB->GetMaximum(), histo_twdr_ptEffB->GetMaximum());
	double maxEffCpt = TMath::Max(histo_tt_ptEffC->GetMaximum(), histo_twdr_ptEffC->GetMaximum());
	double maxEffLpt = TMath::Max(histo_tt_ptEffL->GetMaximum(), histo_twdr_ptEffL->GetMaximum());
	
	TCanvas *c1effBpt = new TCanvas();
       
       histo_tt_ptEffB->Draw("h");
       histo_tt_ptEffB->SetMaximum(maxEffBpt*1.5 );
       histo_tt_ptEffB->GetYaxis()->SetTitle("Eff");
       histo_tt_ptEffB->GetXaxis()->SetTitle("pt bjets");

       histo_twdr_ptEffB->Draw("h, sames");
       histo_twdr_ptEffB->SetMaximum(maxEffBpt*1.5);
     
       
      
      legEffBpt->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1effBpt->SaveAs(plotAnalysis + plotExtension +  modeString[mode] + "_pt_eff_Bjets" + ".png");
	
	
       TCanvas *c1effCpt = new TCanvas();
       
       histo_tt_ptEffC->Draw("h");
       histo_tt_ptEffC->SetMaximum(maxEffCpt*1.5 );
       histo_tt_ptEffC->GetYaxis()->SetTitle("Eff");
       histo_tt_ptEffC->GetXaxis()->SetTitle("pt cjets");

       histo_twdr_ptEffC->Draw("h, sames");
       histo_twdr_ptEffC->SetMaximum(maxEffCpt*1.5);
 
      
      legEffCpt->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1effCpt->SaveAs(plotAnalysis + plotExtension +  modeString[mode] + "_pt_eff_Cjets" + ".png");
      
      
      	TCanvas *c1effLpt = new TCanvas();
       
       histo_tt_ptEffL->Draw("h");
       histo_tt_ptEffL->SetMaximum(maxEffLpt*1.5 );
       histo_tt_ptEffL->GetYaxis()->SetTitle("Eff");
       histo_tt_ptEffL->GetXaxis()->SetTitle("pt Ljets");

       histo_twdr_ptEffL->Draw("h, sames");
       histo_twdr_ptEffL->SetMaximum(maxEffLpt*1.5);
     
       
      
      legEffLpt->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1effLpt->SaveAs(plotAnalysis + plotExtension +  modeString[mode] + "_pt_eff_Ljets" + ".png");
      
      
      
      
     
      
        ///////////////////
   // eta efficiency 
    //////////////////
    TH1F* histo_atw_etaAllB;
    TH1F* histo_tw_etaAllB; 
    TH1F* histo_tt_etaAllB;
    TH1F* histo_atw_etaAllC;
    TH1F* histo_tw_etaAllC; 
    TH1F* histo_tt_etaAllC;
    TH1F* histo_atw_etaAllL;
    TH1F* histo_tw_etaAllL; 
    TH1F* histo_tt_etaAllL;
    TH1F* histo_atw_etaBtagB;
    TH1F* histo_tw_etaBtagB; 
    TH1F* histo_tt_etaBtagB;
    TH1F* histo_atw_etaBtagC;
    TH1F* histo_tw_etaBtagC; 
    TH1F* histo_tt_etaBtagC;
    TH1F* histo_atw_etaBtagL;
    TH1F* histo_tw_etaBtagL; 
    TH1F* histo_tt_etaBtagL;
   
    TH1F* histo_twdr_etaBtagB;
    TH1F* histo_twdr_etaBtagC; 
    TH1F* histo_twdr_etaBtagL;
    TH1F* histo_twdr_etaAllB;
    TH1F* histo_twdr_etaAllC; 
    TH1F* histo_twdr_etaAllL;
   
    TH1F* histo_twdr_etaEffB; 
    TH1F* histo_tt_etaEffB;
    TH1F* histo_twdr_etaEffC; 
    TH1F* histo_tt_etaEffC;
    TH1F* histo_twdr_etaEffL; 
    TH1F* histo_tt_etaEffL;
    
    

    
    legEffBeta = new TLegend(0.7,0.7,0.94,0.94);
     legEffBeta ->SetFillStyle(1001);
     legEffBeta ->SetFillColor(kWhite);
     legEffBeta ->SetBorderSize(1);
     
    legEffCeta = new TLegend(0.7,0.7,0.94,0.94);
     legEffCeta ->SetFillStyle(1001);
     legEffCeta ->SetFillColor(kWhite);
     legEffCeta ->SetBorderSize(1);
     
     legEffLeta = new TLegend(0.7,0.7,0.94,0.94);
     legEffLeta ->SetFillStyle(1001);
     legEffLeta ->SetFillColor(kWhite);
     legEffLeta ->SetBorderSize(1);
    
   
     histo_atw_etaAllB= (TH1F*) _file0->Get(cutLabel[6]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw eta all B: " << cutLabel[6]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_etaAllB= (TH1F*) _file1->Get(cutLabel[6]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw eta all B: " << cutLabel[6]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_etaAllB= (TH1F*) _file2->Get(cutLabel[6]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt eta all B: " << cutLabel[6]+ "_" + modeString[mode] + "_"+processName[2] << endl;  
	histo_atw_etaAllC= (TH1F*) _file0->Get(cutLabel[7]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw eta all C: " << cutLabel[7]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_etaAllC= (TH1F*) _file1->Get(cutLabel[7]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw eta all C: " << cutLabel[7]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_etaAllC= (TH1F*) _file2->Get(cutLabel[7]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt eta all C: " << cutLabel[7]+ "_" + modeString[mode] + "_"+processName[2] << endl;   
	histo_atw_etaAllL= (TH1F*) _file0->Get(cutLabel[8]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw eta all L: " << cutLabel[8]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_etaAllL= (TH1F*) _file1->Get(cutLabel[8]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw eta all L: " << cutLabel[8]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_etaAllL= (TH1F*) _file2->Get(cutLabel[8]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt eta all L: " << cutLabel[8]+ "_" + modeString[mode] + "_"+processName[2] << endl; 
	histo_atw_etaBtagB= (TH1F*) _file0->Get(cutLabel[9]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw eta btag B: " << cutLabel[9]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_etaBtagB= (TH1F*) _file1->Get(cutLabel[9]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw eta btag B: " << cutLabel[9]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_etaBtagB= (TH1F*) _file2->Get(cutLabel[9]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt eta btag B: " << cutLabel[9]+ "_" + modeString[mode] + "_"+processName[2] << endl;  
	histo_atw_etaBtagC= (TH1F*) _file0->Get(cutLabel[10]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw eta btag C: " << cutLabel[10]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_etaBtagC= (TH1F*) _file1->Get(cutLabel[10]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw eta btag C: " << cutLabel[10]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_etaBtagC= (TH1F*) _file2->Get(cutLabel[10]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt eta btag C: " << cutLabel[10]+ "_" + modeString[mode] + "_"+processName[2] << endl;   
	histo_atw_etaBtagL= (TH1F*) _file0->Get(cutLabel[11]+ "_" + modeString[mode] + "_"+processName[0]);
	     cout << "atw eta btag L: " << cutLabel[11]+ "_"  +  modeString[mode] + "_"+ processName[0] << endl; 
      histo_tw_etaBtagL= (TH1F*) _file1->Get(cutLabel[11]+ "_" + modeString[mode] + "_"+processName[1]);
	     cout << "tw eta btag L: " << cutLabel[11]+ "_" + modeString[mode] + "_"+processName[1] << endl; 	
      histo_tt_etaBtagL= (TH1F*) _file2->Get(cutLabel[11]+ "_" + modeString[mode] + "_"+processName[2]);
	     cout << "tt eta btag L: " << cutLabel[11]+ "_" + modeString[mode] + "_"+processName[2] << endl; 
	     
	histo_twdr_etaBtagB = (TH1F*) histo_atw_etaBtagB->Clone("histo_twdr_etaBtagB");
	histo_twdr_etaBtagC = (TH1F*) histo_atw_etaBtagC->Clone("histo_twdr_etaBtagC");
	histo_twdr_etaBtagL = (TH1F*) histo_atw_etaBtagL->Clone("histo_twdr_etaBtagL");
	histo_twdr_etaAllB = (TH1F*)histo_atw_etaAllB->Clone("histo_twdr_etaAllB");
	histo_twdr_etaAllC = (TH1F*)histo_atw_etaAllC->Clone("histo_twdr_etaAllC");
	histo_twdr_etaAllL = (TH1F*)histo_atw_etaAllL->Clone("histo_twdr_etaAllL");
	histo_twdr_etaBtagB->Add(histo_tw_etaBtagB);
	histo_twdr_etaBtagC->Add(histo_tw_etaBtagC);
	histo_twdr_etaBtagL->Add(histo_tw_etaBtagL);
	histo_twdr_etaAllB->Add(histo_tw_etaAllB);
	histo_twdr_etaAllC->Add(histo_tw_etaAllC);
	histo_twdr_etaAllL->Add(histo_tw_etaAllL);
	
	
	     
	histo_twdr_etaEffB = (TH1F*) histo_twdr_etaBtagB->Clone("histo_twdr_etaEffB"); 
	histo_twdr_etaEffC = (TH1F*) histo_twdr_etaBtagC->Clone("histo_twdr_etaEffC");
	histo_twdr_etaEffL = (TH1F*) histo_twdr_etaBtagL->Clone("histo_twdr_etaEffL");
	histo_tt_etaEffB = (TH1F*) histo_tt_etaBtagB->Clone("histo_tt_etaEffB"); 
	histo_tt_etaEffC = (TH1F*) histo_tt_etaBtagC->Clone("histo_tt_etaEffC");
	histo_tt_etaEffL = (TH1F*) histo_tt_etaBtagL->Clone("histo_tt_etaEffL");
	
	
	histo_twdr_etaEffB->Divide(histo_twdr_etaBtagB,histo_twdr_etaAllB,1.,1., "B"); 
    	histo_tt_etaEffB->Divide(histo_tt_etaBtagB,histo_tt_etaAllB,1.,1.,"B"); 
 	histo_twdr_etaEffC->Divide(histo_twdr_etaBtagC,histo_twdr_etaAllC,1.,1., "B"); 
    	histo_tt_etaEffC->Divide(histo_tt_etaBtagC,histo_tt_etaAllC,1.,1.,"B");
        histo_twdr_etaEffL->Divide(histo_twdr_etaBtagL,histo_twdr_etaAllL,1.,1., "B"); 
    	histo_tt_etaEffL->Divide(histo_tt_etaBtagL,histo_tt_etaAllL,1.,1.,"B");
	
	legEffBeta ->AddEntry(histo_tt_etaEffB, "ttbar bckgr. eff" , "l");  
        legEffBeta ->AddEntry(histo_twdr_etaEffB, "twdr signal eff" , "l");
	
	legEffCeta ->AddEntry(histo_tt_etaEffC, "ttbar bckgr. eff" , "l");  
        legEffCeta ->AddEntry(histo_twdr_etaEffC, "twdr signal eff" , "l");
	
	legEffLeta ->AddEntry(histo_tt_etaEffL, "ttbar bckgr. eff" , "l");  
        legEffLeta ->AddEntry(histo_twdr_etaEffL, "twdr signal eff" , "l");
	
	double maxEffBeta = TMath::Max(histo_tt_etaEffB->GetMaximum(), histo_twdr_etaEffB->GetMaximum());
	double maxEffCeta = TMath::Max(histo_tt_etaEffC->GetMaximum(), histo_twdr_etaEffC->GetMaximum());
	double maxEffLeta = TMath::Max(histo_tt_etaEffL->GetMaximum(), histo_twdr_etaEffL->GetMaximum());
	
	TCanvas *c1effBeta = new TCanvas();
       
       histo_tt_etaEffB->Draw("h");
       histo_tt_etaEffB->SetMaximum(maxEffBeta*1.5 );
       histo_tt_etaEffB->GetYaxis()->SetTitle("Eff");
       histo_tt_etaEffB->GetXaxis()->SetTitle("eta bjets");

       histo_twdr_etaEffB->Draw("h, sames");
       histo_twdr_etaEffB->SetMaximum(maxEffBeta*1.5);
    
       
      
      legEffBeta->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1effBeta->SaveAs(plotAnalysis + plotExtension +  modeString[mode] + "_eta_eff_Bjets" + ".png");
	
	
       TCanvas *c1effCeta = new TCanvas();
       
       histo_tt_etaEffC->Draw("h");
       histo_tt_etaEffC->SetMaximum(maxEffCeta*1.5 );
       histo_tt_etaEffC->GetYaxis()->SetTitle("Eff");
       histo_tt_etaEffC->GetXaxis()->SetTitle("eta cjets");

       histo_twdr_etaEffC->Draw("h, sames");
       histo_twdr_etaEffC->SetMaximum(maxEffCeta*1.5);
   
       
      
      legEffCeta->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1effCeta->SaveAs(plotAnalysis + plotExtension +  modeString[mode] + "_eta_eff_Cjets" + ".png");
      
      
      	TCanvas *c1effLeta = new TCanvas();
       
       histo_tt_etaEffL->Draw("h");
       histo_tt_etaEffL->SetMaximum(maxEffLeta*1.5 );
       histo_tt_etaEffL->GetYaxis()->SetTitle("Eff");
       histo_tt_etaEffL->GetXaxis()->SetTitle("eta Ljets");

       histo_twdr_etaEffL->Draw("h, sames");
       histo_twdr_etaEffL->SetMaximum(maxEffLeta*1.5);
    
       
      
      legEffLeta->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1effLeta->SaveAs(plotAnalysis + plotExtension +  modeString[mode] + "_eta_eff_Ljets" + ".png");
      
      //Writing in a rootfile
      
      char eff_fileB[1000];
     sprintf(eff_fileB,"/user/ivanpari/MasterThesis_newEnv/TopBrussels/TopTreeAnalysis/Isis-SingleTop-tW/BtagFiles/rootfiles/EFFbtag_%d_Bjets.root", mode);
      TFile *fileEfB = new TFile(eff_fileB,"RECREATE"); 
     histo_tt_etaEffB->SetDirectory(fileEfB) ;
      histo_twdr_etaEffB->SetDirectory(fileEfB); 
      histo_tt_ptEffB->SetDirectory(fileEfB) ;
      histo_twdr_ptEffB->SetDirectory(fileEfB) ;          
      fileEfB->Write();
      fileEfB->Close();     

       char eff_fileC[1000];
      sprintf(eff_fileC,"/user/ivanpari/MasterThesis_newEnv/TopBrussels/TopTreeAnalysis/Isis-SingleTop-tW/BtagFiles/rootfiles/EFFbtag_%d_Cjets.root", mode);
      TFile *fileEfC = new TFile(eff_fileC,"RECREATE"); 
      histo_tt_etaEffC->SetDirectory(fileEfC); 
      histo_twdr_etaEffC->SetDirectory(fileEfC) ;
      histo_tt_ptEffC->SetDirectory(fileEfC); 
      histo_twdr_ptEffC->SetDirectory(fileEfC) ;          
      fileEfC->Write();
      fileEfC->Close();
      
       char eff_fileL[1000];
      sprintf(eff_fileL,"/user/ivanpari/MasterThesis_newEnv/TopBrussels/TopTreeAnalysis/Isis-SingleTop-tW/BtagFiles/rootfiles/EFFbtag_%d_Ljets.root", mode);
      TFile *fileEfL = new TFile(eff_fileL,"RECREATE"); 
      histo_tt_etaEffL->SetDirectory(fileEfL) ;
      histo_twdr_etaEffL->SetDirectory(fileEfL) ;
      histo_tt_ptEffL->SetDirectory(fileEfL) ;
      histo_twdr_ptEffL->SetDirectory(fileEfL);           
      fileEfL->Write();
      fileEfL->Close();      
     
      
      
} // end constructor loop
