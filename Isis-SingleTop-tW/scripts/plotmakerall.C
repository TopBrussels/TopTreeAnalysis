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

void plotmakerall(){
  
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  gROOT->SetBatch(1);
  
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
  
  labelcms2->AddText("12 fb^{-1}, ee/e#mu/#mu#mu");
    //labelcms2->AddText("2.1 fb^{-1}, ee/e#mu/#mu#mu");
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
  
  

   const int nProcess = 7;
   const int nPlots = 16;
  TString processName[nProcess] =  { "twdr", "st", "tt","di", "zjets", "wjets",  "data"};
  TString processTitle[nProcess] = { "tW", "t/s-channel", "t#bar{t}", "WW", "Z/#gamma*+jets", "W+jets",   "data"};
  Color_t color[nProcess] =        { kWhite, kMagenta-10, kRed+1, kYellow-10,  kAzure-2, kGreen-3,  kBlack};
  

  TString cutLabel[nPlots] =     { "cuts", "met", "mll", "njets", "njetsbt", "ptsys", "ht", "pt_leading", "nvertex_2lep", "nvertex", "pt_max", "pt_min","njets_final","njetsbt_final","et_jet",  "met_zgamma"};
  int rebinHisto[nPlots] =       { 1, 4, 4, 1, 1, 4, 12, 4, 1, 1,2, 2,1,1,4,5};
  TString cutTitle[nPlots] =     { "Analysis Cut", "E_{T}^{miss}", "Inv. Mass", "# of jets", "# of jets(bt)" ,  "P_{T} system [GeV]", "H_{T} [GeV]","P_{T} of the leading jet", "# of vertex after lep sel","# of  vertex ",  "p_T of the first lepton [GeV]", "p_T  of the second lepton [GeV]", "# of jets final", "# of jets(bt) final","transverse energy jet (GeV)", "E_{T}^{miss} inside Zmass window"};


 
 
 
  TH1F*  h [nPlots][nProcess];
  TH1F*   histo[3][nPlots][nProcess];
  THStack* hStack[nPlots];

  for (int iVariable = 0; iVariable < nPlots; iVariable++){
    leg = new TLegend(0.7,0.7,0.94,0.94);
    leg ->SetFillStyle(1001);
    leg ->SetFillColor(kWhite);
    leg ->SetBorderSize(1);

    hStack[iVariable] = new THStack(cutLabel[iVariable],cutLabel[iVariable]);
    for (int iProcess = 0; iProcess < nProcess; iProcess++){
      sprintf(myRootFile,"results/an_11966pb_0.root");
      
      TFile *_file0 = TFile::Open(myRootFile);
      histo[0][iVariable][iProcess] = (TH1F*) _file0->Get(cutLabel[iVariable]+ "_" + processName[iProcess]);
      histo[0][iVariable][iProcess]->Rebin(rebinHisto[iVariable]);
      histo[0][iVariable][iProcess]->SetFillColor(color[iProcess]);
      histo[0][iVariable][iProcess]->SetLineColor(kBlack);
      histo[0][iVariable][iProcess]->SetLineWidth(1);
      if (iProcess == 0) histo[0][iVariable][iProcess]->SetLineColor(kBlack);
     // if (iProcess == 4) histo[0][iVariable][iProcess]->Scale(1.77);
      
      sprintf(myRootFile,"results/an_12067pb_1.root");
      
      TFile *_file1 = TFile::Open(myRootFile);
      histo[1][iVariable][iProcess] = (TH1F*) _file1->Get(cutLabel[iVariable]+ "_" + processName[iProcess]);
      histo[1][iVariable][iProcess]->Rebin(rebinHisto[iVariable]);
      histo[1][iVariable][iProcess]->SetFillColor(color[iProcess]);
      histo[1][iVariable][iProcess]->SetLineColor(kBlack);
      histo[1][iVariable][iProcess]->SetLineWidth(1);
      if (iProcess == 0) histo[1][iVariable][iProcess]->SetLineColor(kBlack);
     // if (iProcess == 4) histo[1][iVariable][iProcess]->Scale(1.59);
      
      
      sprintf(myRootFile,"results/an_12093pb_2.root");
  
      TFile *_file2 = TFile::Open(myRootFile);
      histo[2][iVariable][iProcess] = (TH1F*) _file2->Get(cutLabel[iVariable]+ "_" + processName[iProcess]);
      histo[2][iVariable][iProcess]->Rebin(rebinHisto[iVariable]);
      histo[2][iVariable][iProcess]->SetFillColor(color[iProcess]);
      histo[2][iVariable][iProcess]->SetLineColor(kBlack);
      histo[2][iVariable][iProcess]->SetLineWidth(1);
      if (iProcess == 0) histo[2][iVariable][iProcess]->SetLineColor(kBlack);
    //  if (iProcess == 4) histo[2][iVariable][iProcess]->Scale(1.95);

      h[iVariable][iProcess] =  histo[0][iVariable][iProcess];
      h[iVariable][iProcess]->Add(histo[1][iVariable][iProcess]);
      h[iVariable][iProcess]->Add(histo[2][iVariable][iProcess]);
      
      
      
    }


    h[iVariable][5]->Add(h[iVariable][1]);
    h[iVariable][5]->Add(h[iVariable][3]);
    
    
    hStack[iVariable]->Add(h[iVariable][5]);
    hStack[iVariable]->Add(h[iVariable][4]);
    hStack[iVariable]->Add(h[iVariable][2]);
    hStack[iVariable]->Add(h[iVariable][0]);
    
    leg->AddEntry(h[iVariable][6], processTitle[6], "p");
    leg->AddEntry(h[iVariable][0], processTitle[0], "f");
    leg->AddEntry(h[iVariable][2], processTitle[2], "f");
    leg->AddEntry(h[iVariable][4], processTitle[4], "f");
    leg->AddEntry(h[iVariable][5], "Other", "f");
     
    h[iVariable][6]->SetMarkerStyle(20);
    h[iVariable][6]->SetMarkerSize(1.2);
    h[iVariable][6]->SetLineWidth(1);
    h[iVariable][6]->SetMarkerColor(kBlack);
    h[iVariable][6]->SetLineColor(kBlack);
    
    double max = TMath::Max(hStack[iVariable]->GetMaximum(), h[iVariable][6]->GetMaximum());
    TCanvas *c1 = new TCanvas();
    hStack[iVariable]->Draw("histo");
    hStack[iVariable]->SetMaximum(max * 1.2);
    hStack[iVariable]->SetMinimum(1);
    hStack[iVariable]->GetXaxis()->SetTitle(cutTitle[iVariable]);
    hStack[iVariable]->GetYaxis()->SetTitle("events / 12 fb^{-1}");
    //hStack[iVariable]->GetYaxis()->SetTitle("events / 2.1 fb^{-1}");

    hStack[iVariable]->GetYaxis()->CenterTitle(); 
    
    if (iVariable == 0){
      hStack[iVariable]->GetXaxis()->SetBinLabel(2,"Lepton Selection");
      hStack[iVariable]->GetXaxis()->SetBinLabel(3,"m_{ll}");
      hStack[iVariable]->GetXaxis()->SetBinLabel(4,"E_{T}^{miss}");
      hStack[iVariable]->GetXaxis()->SetBinLabel(5,"1 jet");
      hStack[iVariable]->GetXaxis()->SetBinLabel(6,"b-tag");
      hStack[iVariable]->GetXaxis()->SetBinLabel(7,"H_{T}");
      hStack[iVariable]->GetXaxis()->SetRangeUser(1,6);
     }
    
    if (iVariable == 5) hStack[iVariable]->GetYaxis()->SetRangeUser(1,100);
    h[iVariable][6]->Draw("e, sames");
    leg->Draw();
    labelcms->Draw();
    labelcms2->Draw();
    
    c1->SaveAs("allTogether/plots/plot_3_" + cutLabel[iVariable] + ".png");
    c1->SaveAs("allTogether/plots/pdf/plot_3_" + cutLabel[iVariable] + ".pdf");
    c1->SetLogy();
    hStack[iVariable]->SetMaximum(max * 10);
    c1->SaveAs("allTogether/plots/plot_3_" + cutLabel[iVariable] + "_log.png");
    c1->SaveAs("allTogether/plots/pdf/plot_3_" + cutLabel[iVariable] + "_log.pdf");
    
    


     delete _file0;
     delete _file1;
     delete _file2;
  }
 

}
