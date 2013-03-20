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
  

  double lumi = 0;
  
    if      (mode == 0)	 lumi = 11966.617;  	
    else if (mode == 1) lumi = 12067.294;  	
    else if (mode == 2) lumi = 12093.792;

/*
sprintf(name, "data");       
sprintf(name, "data1");      
sprintf(name, "data2");      
sprintf(name, "tt");         
sprintf(name, "tw_dr");      
sprintf(name, "atw_dr");     
sprintf(name, "t");          
sprintf(name, "at");         
sprintf(name, "s");          
sprintf(name, "as");         
sprintf(name, "ww");         
sprintf(name, "wz");         
sprintf(name, "zz");         
sprintf(name, "zjets");      
sprintf(name, "zjets_lowmll")
sprintf(name, "wjets");     
*/ 
  TString name;
  char myRootFile0[300];
  char myRootFile1[300];
  char myRootFile2[300];
  char myRootFile3[300];
  char myRootFile4[300];
  char myRootFile5[300];
  char myRootFile6[300];
  char myRootFile7[300];
  char myRootFile8[300];
  char myRootFile9[300];
  char myRootFile10[300];
  char myRootFile11[300];
  char myRootFile12[300];
  char myRootFile13[300];
  
  
  
  const int nProcess = 13;
  TString processName[nProcess] =  { "tw_dr", "atw_dr", "s", "as","t", "at",  "tt","ww","wz","zz", "zjets","zjets_lowmll", "wjets",  "data"};
  
  
  
  name = processName[0];
  sprintf(myRootFile0,"../outputs/out_%d_tw_dr.root", mode);
  name = processName[1];
  sprintf(myRootFile1,"../outputs/out_%d_atw_dr.root", mode);
  name = processName[2];
  sprintf(myRootFile2,"../outputs/out_%d_s.root", mode);
  name = processName[3];
  sprintf(myRootFile3,"../outputs/out_%d_as.root", mode);
  name = processName[3];
  sprintf(myRootFile3,"../outputs/out_%d_t.root", mode);
  name = processName[4];
  sprintf(myRootFile4,"../outputs/out_%d_at.root", mode);
  name = processName[5];
  sprintf(myRootFile5,"../outputs/out_%d_tt.root", mode);
  name = processName[6];
  sprintf(myRootFile6,"../outputs/out_%d_ww.root", mode);
  name = processName[7];
  sprintf(myRootFile7,"../outputs/out_%d_wz.root", mode);
  name = processName[8];
  sprintf(myRootFile8,"../outputs/out_%d_zz.root", mode);
  name = processName[9];
  sprintf(myRootFile9,"../outputs/out_%d_zjets.root", mode);
  name = processName[10];
  sprintf(myRootFile10,"../outputs/out_%d_zjets_lowmll.root", mode);
  name = processName[11];
  sprintf(myRootFile11,"../outputs/out_%d_wjets.root", mode);
  name = processName[12];
  sprintf(myRootFile12,"../outputs/out_%d_data.root", mode);
 
  
  
  TFile *_file0 = TFile::Open(myRootFile0);
  cout << myRootFile0 << endl;
    TFile *_file1 = TFile::Open(myRootFile1);
  cout << myRootFile1 << endl;
    TFile *_file2 = TFile::Open(myRootFile2);
  cout << myRootFile2 << endl;
    TFile *_file3 = TFile::Open(myRootFile3);
  cout << myRootFile3 << endl;
    TFile *_file4 = TFile::Open(myRootFile4);
  cout << myRootFile4 << endl;
    TFile *_file5 = TFile::Open(myRootFile5);
  cout << myRootFile5 << endl;
    TFile *_file6 = TFile::Open(myRootFile6);
  cout << myRootFile6 << endl;
    TFile *_file7 = TFile::Open(myRootFile7);
  cout << myRootFile7 << endl;
    TFile *_file8 = TFile::Open(myRootFile8);
  cout << myRootFile8 << endl;
    TFile *_file9 = TFile::Open(myRootFile9);
  cout << myRootFile9 << endl;
    TFile *_file10 = TFile::Open(myRootFile10);
  cout << myRootFile10 << endl;
    TFile *_file11 = TFile::Open(myRootFile11);
  cout << myRootFile11 << endl;
    TFile *_file12 = TFile::Open(myRootFile12);
  cout << myRootFile12 << endl;
 
  


  
  
  
  
  TString processTitle[nProcess] = { "tW", "tbar W" "s-channel t", "s-channel at", "t-channel t", "t-channel at","t#bar{t}", "WW","WZ","ZZ", "Z/#gamma*+jets", "Z/#gamma*+jets low mll", "W+jets", "data"};
  Color_t color[nProcess] =        {kWhite, kMagenta-10, kRed+1, kYellow-10,   kAzure-2, kGreen-3, 40, kBlack, kBlue, kMagenta,kRed,kGreen,kYellow};
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
    
    //get histo
    h[0] = (TH3F*) _file0->Get("3d_btagged_tightjets_" + processName[0]);
    // Set the fill color of the histogram and line 
    h[0]->SetFillColor(color[0]);
    h[0]->SetLineColor(kBlack);
        //get histo
    h[1] = (TH3F*) _file1->Get("3d_btagged_tightjets_" + processName[1]);
    // Set the fill color of the histogram and line 
    h[1]->SetFillColor(color[1]);
    h[1]->SetLineColor(kBlack);
   
    //get histo
    h[2] = (TH3F*) _file2->Get("3d_btagged_tightjets_" + processName[2]);
    // Set the fill color of the histogram and line 
    h[2]->SetFillColor(color[2]);
    h[2]->SetLineColor(kBlack);
        //get histo
    h[3] = (TH3F*) _file3->Get("3d_btagged_tightjets_" + processName[3]);
    // Set the fill color of the histogram and line 
    h[3]->SetFillColor(color[3]);
    h[3]->SetLineColor(kBlack);
        //get histo
    h[4] = (TH3F*) _file4->Get("3d_btagged_tightjets_" + processName[4]);
    // Set the fill color of the histogram and line 
    h[4]->SetFillColor(color[4]);
    h[4]->SetLineColor(kBlack);
        //get histo
    h[5] = (TH3F*) _file5->Get("3d_btagged_tightjets_" + processName[5]);
    // Set the fill color of the histogram and line 
    h[5]->SetFillColor(color[5]);
    h[5]->SetLineColor(kBlack);
    //get histo
    h[6] = (TH3F*) _file6->Get("3d_btagged_tightjets_" + processName[6]);
    // Set the fill color of the histogram and line 
    h[6]->SetFillColor(color[6]);
    h[6]->SetLineColor(kBlack);
    //get histo
    h[7] = (TH3F*) _file7->Get("3d_btagged_tightjets_" + processName[7]);
    // Set the fill color of the histogram and line 
    h[7]->SetFillColor(color[7]);
    h[7]->SetLineColor(kBlack);
    //get histo
    h[8] = (TH3F*) _file8->Get("3d_btagged_tightjets_" + processName[8]);
    // Set the fill color of the histogram and line 
    h[8]->SetFillColor(color[8]);
    h[8]->SetLineColor(kBlack);
    //get histo
    h[9] = (TH3F*) _file9->Get("3d_btagged_tightjets_" + processName[9]);
    // Set the fill color of the histogram and line 
    h[9]->SetFillColor(color[9]);
    h[9]->SetLineColor(kBlack);
    //get histo
    h[10] = (TH3F*) _file10->Get("3d_btagged_tightjets_" + processName[10]);
    // Set the fill color of the histogram and line 
    h[10]->SetFillColor(color[10]);
    h[10]->SetLineColor(kBlack);
    //get histo
    h[11] = (TH3F*) _file11->Get("3d_btagged_tightjets_" + processName[11]);
    // Set the fill color of the histogram and line 
    h[11]->SetFillColor(color[11]);
    h[11]->SetLineColor(kBlack);
    //get histo
    h[12] = (TH3F*) _file12->Get("3d_btagged_tightjets_" + processName[12]);
    // Set the fill color of the histogram and line 
    h[12]->SetFillColor(color[12]);
    h[12]->SetLineColor(kBlack);
    //get histo
    h[13] = (TH3F*) _file13->Get("3d_btagged_tightjets_" + processName[13]);
    // Set the fill color of the histogram and line 
    h[13]->SetFillColor(color[13]);
    h[13]->SetLineColor(kBlack);
    
    
  for (int iProcess = 0; iProcess <nProcess; iProcess++){

    
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
