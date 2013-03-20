//Rebeca Gonzalez Suarez
//rebeca@cern.ch

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "setTDRStyle.C"
using namespace std;

void controlplots_main(int mode = 0){
 
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
  
  if (mode == 0) labelcms2->AddText("12 fb^{-1}, e#mu channel  ");
  if (mode == 1) labelcms2->AddText("12 fb^{-1}, #mu#mu channel  ");
  if (mode == 2) labelcms2->AddText("12 fb^{-1}, ee channel  ");
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

  sprintf(myRootFile,"outputs/an_%dpb_%d.root", lumi, mode);
  
  
  TFile *_file0 = TFile::Open(myRootFile);
  cout << myRootFile << endl;

  const int nProcess = 14;
  TString processName[nProcess] =  { "twdr", "atwdr", "s", "as","t", "at",  "tt","ww","wz","zz", "zjets","zjets_lowmll", "wjets",  "data"};
  TString processTitle[nProcess] = { "tW", "tbar W" "s-channel t", "s-channel at", "t-channel t", "t-channel at","t#bar{t}", "WW","WZ","ZZ", "Z/#gamma*+jets", "Z/#gamma*+jets low mll", "W+jets", "data"};
  Color_t color[nProcess] =        {kWhite, kMagenta-10, kRed+1, kYellow-10,   kAzure-2, kGreen-3, 40, kBlack, kBlue, kMagenta,kRed,kGreen,kYellow,kAzure};
  TString modeString[3] = {"0", "1", "2"};

  TString cutLabel = "R";
 
  
  TH3F*  h [nProcess];
  TH1F*  histo[nProcess];
  THStack* hStack;
  
  leg = new TLegend(0.7,0.7,0.94,0.94);
  leg ->SetFillStyle(1);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->SetTextSize(0.03);
  hStack = new THStack();
  hStack2 = new THStack();
  for (int iProcess = 0; iProcess <nProcess; iProcess++){
    // Get the histogram of for each process and store it as h[iProcess], so h[0] is twdr, h[1] is atwdr, etc 
    h[iProcess] = (TH3F*) _file0->Get("3d_btagged_tightjets_" + processName[iProcess]);
    // Set the fill color of the histogram and line 
    h[iProcess]->SetFillColor(color[iProcess]);
    h[iProcess]->SetLineColor(kBlack);
    
    //Define histo as histo+processname, eg histotwdr, histoatwdr
    histo[iProcess] = new TH1F("histo"+processName[iProcess], "", 3, 0, 3);
    // set line color and width and fill color
    histo[iProcess]->SetLineColor(kBlack);
    histo[iProcess]->SetLineWidth(1);
    histo[iProcess]->SetFillColor(color[iProcess]);
    // Fill histo with in bin 1 the 1j1t no loose, etc  for each process 
    histo[iProcess]->SetBinContent(1, h[iProcess]->GetBinContent(2,2,2));   // 1loose 1bt 1 jt
    histo[iProcess]->SetBinContent(2, h[iProcess]->GetBinContent(2,2,3));   // 1 loose 1 bt 2j
    histo[iProcess]->SetBinContent(3, h[iProcess]->GetBinContent(3,3,3));  // 2l 2 bt 2j
    histo[iProcess]->SetBinError(1, h[iProcess]->GetBinError(2,2,2));
    histo[iProcess]->SetBinError(2, h[iProcess]->GetBinError(2,2,3));
    histo[iProcess]->SetBinError(3, h[iProcess]->GetBinError(3,3,3));
    
 //   if (iProcess == 0 || iProcess == 1){
//      h[iProcess]->SetLineColor(kBlack);
//      histo[iProcess]->SetLineColor(kBlack);
//    }
   
  }  
  
  //Make the regions, twdr, ttbar, zjets , other
  // other = wjets, s, t, ww, zz, wz
  histo[12]->Add(histo[12]);
  histo[12]->Add(histo[9]);
  histo[12]->Add(histo[8]);
  histo[12]->Add(histo[7]);
  histo[12]->Add(histo[5]);
  histo[12]->Add(histo[4]);
  histo[12]->Add(histo[3]);
  histo[12]->Add(histo[2]);
  
  //twdr = twdr + atwdr
  histo[0]->Add(histo[1]);
  
  //zjets = zjets + zjets low
  histo[10]->Add(histo[11]);
 
  
 // put the inputs on top of eachother 
  hStack->Add(histo[12]); // other
  hStack->Add(histo[0]);  // twdr
  hStack->Add(histo[10]); // zjets
  hStack->Add(histo[6]);  // ttbar
  
  // put data in points to legend
  if (mode == 0) leg->AddEntry(histo[13], processTitle[13], "p");
  if (mode == 1) leg->AddEntry(histo[13], processTitle[13], "p");
  if (mode == 2) leg->AddEntry(histo[13], processTitle[13], "p");
  
  // make legends
  leg->AddEntry(histo[0], processTitle[0], "f"); // twdr
  leg->AddEntry(histo[10], processTitle[10], "f"); // zjets
  leg->AddEntry(histo[6], processTitle[6], "f"); // ttbar
  leg->AddEntry(histo[12], "Other", "f"); // other
 
  // make data pretty
  histo[7]->SetMarkerStyle(20);
  histo[7]->SetMarkerSize(1.2);
  histo[7]->SetMarkerColor(kBlack);
  histo[7]->SetLineColor(kBlack);
  histo[7]->SetLineWidth(1);
  
  // calculate max of mc and data for canvas drawing
  double max = TMath::Max(hStack->GetMaximum(), histo[7]->GetMaximum());
  
  // make canvas
  TCanvas *c1 = new TCanvas();
  //draw mc
  hStack->Draw("histo");
  hStack->SetMaximum(max * 1.5);
  hStack->SetMinimum(1);
  hStack->GetYaxis()->SetTitle("events / 12 fb^{-1}");
  hStack->GetYaxis()->SetTitleOffset(1.4);
  hStack->GetYaxis()->CenterTitle(); 
  
  hStack->GetYaxis()->SetTitleOffset(1.4);
  hStack->GetYaxis()->CenterTitle(); 
  hStack->GetXaxis()->SetBinLabel(1,"1 jet 1 tag");
  hStack->GetXaxis()->SetBinLabel(2,"2 jet 1 tag");
  hStack->GetXaxis()->SetBinLabel(3,"2 jet 2 tag");
  
  //draw data
  histo[7]->Draw("e, sames");
  
  //draw legend and labels
  leg->Draw();
  labelcms->Draw();
  labelcms2->Draw();
  
  c1->SaveAs("plots/main_control_tt_" + modeString[mode] + "_" + cutLabel + ".png");
  c1->SaveAs("plots/pdf/main_control_tt_" + modeString[mode] + "_" + cutLabel + ".pdf");
  c1->SetLogy();
  hStack->SetMaximum(max * 10);
  c1->SaveAs("plots/main_control_tt_" + modeString[mode] + "_" + cutLabel + "_log.png");
  c1->SaveAs("plots/pdf/main_control_tt_" + modeString[mode] + "_" + cutLabel + "_log.pdf");
  
}
