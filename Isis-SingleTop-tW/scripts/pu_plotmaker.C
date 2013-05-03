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

void pu_plotmaker(){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  gROOT->SetBatch(1);   // to show histos on canvas
  

  
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
  char myRootFile1[300];
  char myRootFile2[300];
  char myRootFile3[300];

  

    
  sprintf(myRootFile,"../pileupHistos/pileup_MC_Summer12.root"); // take output from looper
  sprintf(myRootFile1,"../pileupHistos/pileup_tW_2012Data53X_UpToRun203002/nominal.root");
  sprintf(myRootFile2,"../pileupHistos/pileup_tW_2012Data53X_UpToRun203002/sys_down.root"); // take output from looper
  sprintf(myRootFile3,"../pileupHistos/pileup_tW_2012Data53X_UpToRun203002/sys_up.root");
  
  TFile *_file0 = TFile::Open(myRootFile);
  cout << myRootFile << endl;
  TFile *_file1 = TFile::Open(myRootFile1);
  cout << myRootFile1 << endl;
  TFile *_file2 = TFile::Open(myRootFile2);
  cout << myRootFile2 << endl;
  TFile *_file3 = TFile::Open(myRootFile3);
  cout << myRootFile3<< endl;

  
  //Make plots   
  TH1F* histo_puMC;
 TH1F* histo_pu;
   TH1F* histo_puUp;
TH1F* histo_puDown;
   
     leg = new TLegend(0.7,0.7,0.94,0.94);
     leg ->SetFillStyle(1001);
     leg ->SetFillColor(kWhite);
     leg ->SetBorderSize(1);
     
     
 
	
	
     histo_puMC= (TH1F*) _file0->Get("pileup");
	 
	       
	       
       
               //histo_puMC->SetFillColor(color[1]);
        histo_puMC->SetLineColor(kRed);
        histo_puMC->SetLineWidth(2);

       histo_pu= (TH1F*) _file1->Get("pileup");
	 
	       
	       
       
               //histo_pu->SetFillColor(color[1]);
        histo_pu->SetLineColor(kBlack);
        histo_pu->SetLineWidth(2);
	
	histo_puDown= (TH1F*) _file2->Get("pileup");
	 
	       
	       
     
               //histo_puDown->SetFillColor(color[1]);
        histo_puDown->SetLineColor(kBlue);
        histo_puDown->SetLineWidth(2);
	
	
	histo_puUp= (TH1F*) _file3->Get("pileup");
	 
	       
	       
       
               //histo_puUp->SetFillColor(color[1]);
        histo_puUp->SetLineColor(kGreen);
        histo_puUp->SetLineWidth(2);
	

      
	      
       leg ->AddEntry(histo_puMC, "PU MC " , "l");  
       leg ->AddEntry(histo_pu, "PU data - normal" , "l"); 
       leg ->AddEntry(histo_puUp, "PU data - SF up" , "l"); 
       leg ->AddEntry(histo_puDown, "PU data - SF down" , "l"); 
       
       
      double max = TMath::Max(histo_puMC->GetMaximum(),histo_puUp->GetMaximum());
    
       TCanvas *c1 = new TCanvas();
       histo_pu->SetMaximum(max*1.2 );
       histo_puMC->SetMaximum(max*1.2 );
       histo_puUp->SetMaximum(max*1.2);
       histo_puDown->SetMaximum(max*1.2 );
       histo_pu->DrawNormalized("h");
      
     //  histo_puMC->SetMinimum(1);
       histo_pu->GetXaxis()->SetTitle("Pile-up");
     
       
         histo_puMC->DrawNormalized("h,sames");
      
      histo_puMC->GetXaxis()->SetTitle("Pile-up");
    histo_puUp->DrawNormalized("h,sames");
      
      histo_puUp->GetXaxis()->SetTitle("Pile-up");
      histo_puDown->DrawNormalized("h,sames");
      
histo_puDown->GetXaxis()->SetTitle("Pile-up");
       
      leg->Draw();

    
      c1->SaveAs("../pileupHistos/pileup.png");
     // c1->SaveAs("plots/pdf/" + plotExtension + modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");
      
     
      
  
    
 
} // end constructor loop
