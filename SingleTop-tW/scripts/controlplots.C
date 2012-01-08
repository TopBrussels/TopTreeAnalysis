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

void controlplots(int mode = 0, int region = 1){
 
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
  labelcms->AddText("CMS Preliminary, #sqrt{s} = 7 TeV");
  labelcms->SetBorderSize(0);
  
    
  labelcms2  = new TPaveText(0.12,0.85,0.5,0.88,"NDCBR");
  labelcms2->SetTextAlign(12);
  labelcms2->SetTextSize(0.045);
  labelcms2->SetFillColor(kWhite);
  
  if (mode == 0) labelcms2->AddText("4.5 fb^{-1}, e#mu channel  ");
  if (mode == 1) labelcms2->AddText("4.5 fb^{-1}, #mu#mu channel  ");
  if (mode == 2) labelcms2->AddText("4.5 fb^{-1}, ee channel  ");
  
  
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
  
  
  
  if (mode == 0 && region != 1){
    cout << "WARNING: this is the emu channel, there's no DY control region ready, loading tt" << endl;
    region =1;
  } 
  char myRootFile[300];
  double lumi = 0;
  if (mode == 0 )       lumi = 4626.297;
  else if ( mode == 1)  lumi = 4534.871;
  else if ( mode == 2)  lumi = 4593.348;
  if(region == 1)  sprintf(myRootFile,"results/crtt_summer_%dpb_%d.root", lumi, mode);
  else sprintf(myRootFile,"results/crdy_%dpb_%d.root", lumi, mode);
  
  
  TFile *_file0 = TFile::Open(myRootFile);
  cout << myRootFile << endl;
  
  TString processName[8] =  { "twdr", "st", "tt","di", "zjets", "wjets",  "qcd_mu", "data"};
  TString processTitle[8] = { "tW", "t/s-channel", "t#bar{t}", "WW/WZ/ZZ", "Z/#gamma*+jets", "W+jets",  "QCD", "data"};
  Color_t color[8] =        {kWhite, kMagenta-10, kRed+1, kYellow-10,  kAzure-2, kGreen-3, 40, kBlack};
  TString modeString[3] = {"0", "1", "2"};

  TString cutLabel = "R";
  TString cutTitle = " ";
  TString nregion = "tt"; 
  if (region != 1) nregion = "dy";
  
  TH1F*  h [8];
  TH1F*  histo[8];
  TH1F*  histo2[8];
  THStack* hStack;
  THStack* hStack2;
  
  leg = new TLegend(0.7,0.7,0.94,0.94);
  leg ->SetFillStyle(1);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->SetTextSize(0.03);
  hStack = new THStack();
  hStack2 = new THStack();
  for (int iProcess = 0; iProcess < 8; iProcess++){
    h[iProcess] = (TH1F*) _file0->Get("R_" + processName[iProcess]);
    
    h[iProcess]->SetFillColor(color[iProcess]);
    h[iProcess]->SetLineColor(kBlack);
    
    histo[iProcess] = new TH1F("histo"+processName[iProcess], "", 3, 0, 3);
    histo[iProcess]->SetLineColor(kBlack);
    histo[iProcess]->SetLineWidth(1);
    histo[iProcess]->SetFillColor(color[iProcess]);
    histo[iProcess]->SetBinContent(1, h[iProcess]->GetBinContent(2));
    histo[iProcess]->SetBinContent(2, h[iProcess]->GetBinContent(7));
    histo[iProcess]->SetBinContent(3, h[iProcess]->GetBinContent(8));
    histo[iProcess]->SetBinError(1, h[iProcess]->GetBinError(2));
    histo[iProcess]->SetBinError(2, h[iProcess]->GetBinError(7));
    histo[iProcess]->SetBinError(3, h[iProcess]->GetBinError(8));
    
    histo2[iProcess] = new TH1F("histo2"+processName[iProcess], "", 5, 0, 5);
    histo2[iProcess]->SetLineColor(kBlack);
    histo2[iProcess]->SetLineWidth(1);
    histo2[iProcess]->SetFillColor(color[iProcess]);
    histo2[iProcess]->SetBinContent(1, h[iProcess]->GetBinContent(2));
    histo2[iProcess]->SetBinContent(2, h[iProcess]->GetBinContent(7));
    histo2[iProcess]->SetBinContent(3, h[iProcess]->GetBinContent(8));
    histo2[iProcess]->SetBinContent(4, h[iProcess]->GetBinContent(12));
    histo2[iProcess]->SetBinContent(5, h[iProcess]->GetBinContent(15));
    histo2[iProcess]->SetBinError(1, h[iProcess]->GetBinError(2));
    histo2[iProcess]->SetBinError(2, h[iProcess]->GetBinError(7));
    histo2[iProcess]->SetBinError(3, h[iProcess]->GetBinError(8));
    histo2[iProcess]->SetBinError(4, h[iProcess]->GetBinError(12));
    histo2[iProcess]->SetBinError(5, h[iProcess]->GetBinError(15));
    
    if (iProcess == 0){
      h[iProcess]->SetLineColor(kBlack);
      histo[iProcess]->SetLineColor(kBlack);
      histo2[iProcess]->SetLineColor(kBlack);
    }
    
    //  if (iProcess < 7) leg->AddEntry(histo[iProcess], processTitle[iProcess], "f");
    //  else leg->AddEntry(histo[iProcess], processTitle[iProcess], "p");
  }
  for (int iP = 0; iP < 7; iP ++){
    int  index = 6 - iP;
    // hStack->Add(histo[index]);
    hStack2->Add(histo2[index]);
  }
  
  
  
  histo[5]->Add(histo[1]);
  histo[5]->Add(histo[3]);
  histo[5]->Add(histo[6]);
  
  
  hStack->Add(histo[5]);
  hStack->Add(histo[4]);
  hStack->Add(histo[2]);
  hStack->Add(histo[0]);
  
  if (mode == 0) leg->AddEntry(histo[7], processTitle[7], "p");
  if (mode == 1) leg->AddEntry(histo[7], processTitle[7], "p");
  if (mode == 2) leg->AddEntry(histo[7], processTitle[7], "p");
  
  leg->AddEntry(histo[0], processTitle[0], "f");
  leg->AddEntry(histo[2], processTitle[2], "f");
  leg->AddEntry(histo[4], processTitle[4], "f");
  leg->AddEntry(histo[5], "Other", "f");
 
  
  histo[7]->SetMarkerStyle(20);
  histo[7]->SetMarkerSize(1.2);
  //histo[7]->SetLineWidth(1);
  histo[7]->SetMarkerColor(kBlack);
  histo[7]->SetLineColor(kBlack);
  
  double max = TMath::Max(hStack->GetMaximum(), histo[7]->GetMaximum());
  TCanvas *c1 = new TCanvas();
  hStack->Draw("histo");
  hStack->SetMaximum(max * 1.25);
  hStack->SetMinimum(1);
  hStack->GetXaxis()->SetTitle(cutTitle);
  //hStack->GetYaxis()->SetTitle("events / 1.1 fb^{-1}");
  //hStack->GetYaxis()->SetTitle("events / 1.4 fb^{-1}");
  hStack->GetYaxis()->SetTitle("events / 2.1 fb^{-1}");
  hStack->GetYaxis()->SetTitleOffset(1.4);
  hStack->GetYaxis()->CenterTitle(); 
  
  if (region == 1){
    hStack->GetYaxis()->SetTitleOffset(1.4);
    hStack->GetYaxis()->CenterTitle(); 
    hStack->GetXaxis()->SetBinLabel(1,"1 jet 1 tag");
    hStack->GetXaxis()->SetBinLabel(2,"2 jet 1 tag");
    hStack->GetXaxis()->SetBinLabel(3,"2 jet 2 tag");
  } else {
    hStack->GetXaxis()->SetRangeUser(1,3);
    hStack->GetXaxis()->SetBinLabel(2,"Out");
    hStack->GetXaxis()->SetBinLabel(3,"In");
    hStack->GetXaxis()->SetRangeUser(1,2);
  }
  
  for (int i = 0; i < 8; i++){
    cout << processName[i] << " ---> 	[signal]: " <<	h[i]->GetBinContent(2) << "	[2j 1t]: " << h[i]->GetBinContent(7) << " 	[2j 2t]: " << h[i]->GetBinContent(8) << " 	[>1j >1t]: " <<  h[i]->GetBinContent(12) << " 	[>1j >1t]:" <<  h[i]->GetBinContent(15) << endl;
  }
  
  
  histo[7]->Draw("e2, sames");
  leg->Draw();
  labelcms->Draw();
  labelcms2->Draw();
  
  c1->SaveAs("plots/control_summer_" + nregion + "_" + modeString[mode] + "_" + cutLabel + ".pdf");
  c1->SetLogy();
  hStack->SetMaximum(max * 10);
  c1->SaveAs("plots/control_summer_" + nregion + "_" + modeString[mode] + "_" + cutLabel + "_log.pdf");
    
  
  histo2[7]->SetMarkerStyle(20);
  histo2[7]->SetMarkerSize(1.2);
  //histo2[7]->SetLineWidth(1);
  histo2[7]->SetMarkerColor(kBlack);
  histo2[7]->SetFillColor(kWhite);
  histo2[7]->SetLineColor(kBlack);
  
  double max = TMath::Max(hStack2->GetMaximum(), histo2[7]->GetMaximum());
  TCanvas *c1 = new TCanvas();
  hStack2->Draw("histo");
  hStack2->SetMaximum(max * 1.25);
  hStack2->SetMinimum(1);
  hStack2->GetXaxis()->SetTitle(cutTitle);
  hStack2->GetYaxis()->SetTitle("events / 2.1 fb^{-1}");
  hStack2->GetYaxis()->SetTitleOffset(1.4);
  hStack2->GetYaxis()->CenterTitle(); 
  
  
  hStack2->GetYaxis()->SetTitleOffset(1.4);
  hStack2->GetYaxis()->CenterTitle(); 
  hStack2->GetXaxis()->SetBinLabel(1,"1 jet 1 tag");
  hStack2->GetXaxis()->SetBinLabel(2,"2 jet 1 tag");
  hStack2->GetXaxis()->SetBinLabel(3,"2 jet 2 tag");
  hStack2->GetXaxis()->SetBinLabel(4,">1 jet 1 tag");
  hStack2->GetXaxis()->SetBinLabel(5,">1 jet >1 tag");
  
  
  
  
  histo2[7]->Draw("e2, sames");
  leg->Draw();
  
  c1->SaveAs("plots/control2_summer_" + nregion + "_" + modeString[mode] + "_" + cutLabel + ".pdf");
  c1->SetLogy();
  hStack2->SetMaximum(max * 10);
  c1->SaveAs("plots/control2_summer_" + nregion + "_" + modeString[mode] + "_" + cutLabel + "_log.pdf");
   
  
  
}
