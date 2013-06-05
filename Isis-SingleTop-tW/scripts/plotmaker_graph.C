

#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMarker.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "TStyle.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "setTDRStyle.C"
using namespace std;

void plotmaker_graph(int mode = 0){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  gROOT->SetBatch(1);
  
  TPaveText *labelcms  = new TPaveText(0.12,0.88,0.5,0.94,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.045);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS Preliminary, #sqrt{s} = 8 TeV");
  labelcms->SetBorderSize(0);
    
  TPaveText *labelcms2  = new TPaveText(0.12,0.85,0.5,0.88,"NDCBR");
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
  
  TPaveText *labelKStest = new TPaveText(0.12,0.80,0.6,0.84,"NDCBR");
  labelKStest->SetTextAlign(12);
  labelKStest->SetTextSize(0.035);
  labelKStest->SetFillColor(kWhite);
  labelKStest->SetBorderSize(0);

  
  char myRootFile[300];
  char myRootFileZSFsysDown[300];
  char myRootFileZSFsysUp[300];
  char myRootFileJERsysDown[300];
  char myRootFileJERsysUp[300];
  char myRootFileJESsysDown[300];
  char myRootFileJESsysUp[300];
  char myRootFilePUsysDown[300];
  char myRootFilePUsysUp[300];
  char myRootFileSFsysDown[300];
  char myRootFileSFsysUp[300];
  char myRootFileMETsysDown[300];
  char myRootFileMETsysUp[300];
  char myRootFileTopMassDown[300];
  char myRootFileTopMassUp[300];
  char myRootFileQ2Down[300];
  char myRootFileQ2Up[300];
  char myRootFileeleSFsysDown[300];
  char myRootFileeleSFsysUp[300];
  char myRootFilematchingDown[300];
  char myRootFilematchingUp[300];
  char myRootFilePDF_tt_1j1t[300];
  char myRootFilePDF_tt_1j1t_Up[300];
  char myRootFilePDF_other_1j1t[300];
  char myRootFilePDF_other_1j1t_Up[300];
  char myRootFilePDF_twdr_1j1t[300];
  char myRootFilePDF_twdr_1j1t_Up[300];
  char myRootFilePDF_tt_2j1t[300];
  char myRootFilePDF_tt_2j1t_Up[300];
  char myRootFilePDF_other_2j1t[300];
  char myRootFilePDF_other_2j1t_Up[300];
  char myRootFilePDF_twdr_2j1t[300];
  char myRootFilePDF_twdr_2j1t_Up[300];
  char myRootFilePDF_tt_2j2t[300];
  char myRootFilePDF_tt_2j2t_Up[300];
  char myRootFilePDF_other_2j2t[300];
  char myRootFilePDF_other_2j2t_Up[300];
  char myRootFilePDF_twdr_2j2t[300];
  char myRootFilePDF_twdr_2j2t_Up[300];
  double lumi = 1000;
  
    if      (mode == 0)	 lumi = 11966.617;  	
    else if (mode == 1) lumi = 12067.294;  	
    else if (mode == 2) lumi = 12093.792;  	
  
  
  
  sprintf(myRootFile,"results/an_%dpb_%d.root", (int)lumi, mode);
  sprintf(myRootFileZSFsysDown,"results/ZSFsysDown_an_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileZSFsysUp,"results/ZSFsysUp_an_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileJERsysDown,"results/JERsysDown_an_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileJERsysUp,"results/JERsysUp_an_%dpb_%d.root", (int)lumi, mode);
  
  sprintf(myRootFileJESsysDown,"results/JESsysDown_an_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileJESsysUp,"results/JESsysUp_an_%dpb_%d.root", (int)lumi, mode);
  
  sprintf(myRootFilePUsysDown,"results/PUsysDown_an_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFilePUsysUp,"results/PUsysUp_an_%dpb_%d.root", (int)lumi, mode);
  
  sprintf(myRootFileSFsysDown,"results/SFsysDown_an_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileSFsysUp,"results/SFsysUp_an_%dpb_%d.root", (int)lumi, mode);
  
  sprintf(myRootFileMETsysDown,"results/METsysDown_an_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileMETsysUp,"results/METsysUp_an_%dpb_%d.root", (int)lumi, mode);
  
   sprintf(myRootFileTopMassDown,"results/TopMassDown_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileTopMassUp,"results/TopMassUp_%dpb_%d.root", (int)lumi, mode);
  
   sprintf(myRootFileQ2Down,"results/Q2Down_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileQ2Up,"results/Q2Up_%dpb_%d.root", (int)lumi, mode);
  
   sprintf(myRootFileeleSFsysDown,"results/eleSFsysDown_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileeleSFsysUp,"results/eleSFsysUp_%dpb_%d.root", (int)lumi, mode);
  
  sprintf(myRootFilematchingDown,"results/matchingDown_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFilematchingUp,"results/matchingUp_%dpb_%d.root", (int)lumi, mode);
  
  sprintf(myRootFilePDF_tt_1j1t,"PDF_unc_output_new/pdf_signal_%d_tt_outputfile_rebinned.root", mode);
  
  sprintf(myRootFilePDF_other_1j1t,"PDF_unc_output_new/pdf_signal_%d_others_outputfile_rebinned.root", mode);
 
  sprintf(myRootFilePDF_twdr_1j1t,"PDF_unc_output_new/pdf_signal_%d_twdr_outputfile_rebinned.root", mode);
  
  sprintf(myRootFilePDF_tt_2j1t,"PDF_unc_output_new/pdf_2j1t_%d_tt_outputfile_rebinned.root", mode);
  
  sprintf(myRootFilePDF_other_2j1t,"PDF_unc_output_new/pdf_2j1t_%d_others_outputfile_rebinned.root", mode);
 
  sprintf(myRootFilePDF_twdr_2j1t,"PDF_unc_output_new/pdf_2j1t_%d_twdr_outputfile_rebinned.root", mode);
  
  sprintf(myRootFilePDF_tt_2j2t,"PDF_unc_output_new/pdf_2j2t_%d_tt_outputfile_rebinned.root", mode);
  
  sprintf(myRootFilePDF_other_2j2t,"PDF_unc_output_new/pdf_2j2t_%d_others_outputfile_rebinned.root", mode);
 
  sprintf(myRootFilePDF_twdr_2j2t,"PDF_unc_output_new/pdf_2j2t_%d_twdr_outputfile_rebinned.root", mode);
 
  
  

  
  TFile *_file0 = TFile::Open(myRootFile);
  TFile *_fileZSFsysUp = TFile::Open(myRootFileZSFsysUp);
  TFile *_fileZSFsysDown = TFile::Open(myRootFileZSFsysDown);
  TFile *_fileJERsysUp = TFile::Open(myRootFileJERsysUp);
  TFile *_fileJERsysDown = TFile::Open(myRootFileJERsysDown);
  TFile *_fileJESsysUp = TFile::Open(myRootFileJESsysUp);
  TFile *_fileJESsysDown = TFile::Open(myRootFileJESsysDown);
  TFile *_filePUsysUp = TFile::Open(myRootFilePUsysUp);
  TFile *_filePUsysDown = TFile::Open(myRootFilePUsysDown);
  TFile *_fileSFsysUp = TFile::Open(myRootFileSFsysUp);
  TFile *_fileSFsysDown = TFile::Open(myRootFileSFsysDown);
  TFile *_fileMETsysUp = TFile::Open(myRootFileMETsysUp);
  TFile *_fileMETsysDown = TFile::Open(myRootFileMETsysDown);
  TFile *_fileTopMassUp = TFile::Open(myRootFileTopMassUp);
  TFile *_fileTopMassDown = TFile::Open(myRootFileTopMassDown);
  TFile *_fileQ2Up = TFile::Open(myRootFileQ2Up);
  TFile *_fileQ2Down = TFile::Open(myRootFileQ2Down);
  TFile *_fileeleSFsysUp = TFile::Open(myRootFileeleSFsysUp);
  TFile *_fileeleSFsysDown = TFile::Open(myRootFileeleSFsysDown);
  TFile *_filematchingUp = TFile::Open(myRootFilematchingUp);
  TFile *_filematchingDown = TFile::Open(myRootFilematchingDown);
  TFile *_filePDF_tt_1j1t = TFile::Open(myRootFilePDF_tt_1j1t);
  TFile *_filePDF_other_1j1t = TFile::Open(myRootFilePDF_other_1j1t);
  TFile *_filePDF_twdr_1j1t = TFile::Open(myRootFilePDF_twdr_1j1t);
  TFile *_filePDF_tt_2j1t = TFile::Open(myRootFilePDF_tt_2j1t);
  TFile *_filePDF_other_2j1t = TFile::Open(myRootFilePDF_other_2j1t);
  TFile *_filePDF_twdr_2j1t = TFile::Open(myRootFilePDF_twdr_2j1t);
  TFile *_filePDF_tt_2j2t = TFile::Open(myRootFilePDF_tt_2j2t);
  TFile *_filePDF_other_2j2t = TFile::Open(myRootFilePDF_other_2j2t);
  TFile *_filePDF_twdr_2j2t = TFile::Open(myRootFilePDF_twdr_2j2t);

 /*  cout << "-------------------------------------------------------" << endl; 
  cout << " ------------ USED FILES ------------------------------" << endl; 
  cout << "Normal: " << myRootFile << endl;
  cout << "ZSFsysDown: " << myRootFileZSFsysDown << endl;
  cout << "ZSFsysUp: " << myRootFileZSFsysUp << endl;
  cout << "JERsysDown: " << myRootFileJERsysDown << endl;
  cout << "JERsysUp: " << myRootFileJERsysUp << endl;
  cout << "JESsysDown: " << myRootFileJESsysDown << endl;
  cout << "JESsysUp: " << myRootFileJESsysUp << endl;
  cout << "PUsysDown: " << myRootFilePUsysDown << endl;
  cout << "PUsysUp: " << myRootFilePUsysUp << endl;
  cout << "SFsysDown: " << myRootFileSFsysDown << endl;
  cout << "SFsysUp: " << myRootFileSFsysUp << endl;
  cout << "METsysDown: " << myRootFileMETsysDown << endl;
  cout << "METsysUp: " << myRootFileMETsysUp << endl;
  cout << "TopMassDown: " << myRootFileTopMassDown << endl;
  cout << "TopMassUp: " << myRootFileTopMassUp << endl;
  cout << "Q2Down: " << myRootFileQ2Down << endl;
  cout << "Q2Up: " << myRootFileQ2Up << endl;
  cout << "eleSFsysDown: " << myRootFileeleSFsysDown << endl;
  cout << "eleSFsysUp: " << myRootFileeleSFsysUp << endl;
  cout << "matchingDown: " << myRootFilematchingDown << endl;
  cout << "matchingUp: " << myRootFilematchingUp << endl;
  cout << "-------------------------------------------------------" << endl; 
  */
  const int nProcess = 10;
  const int nPlots = 28;
  const int nSys = 21;
  TString processName[nProcess] =  { "twdr", "st", "tt","di", "zjets", "wjets", "data","others","twds","mc"};
  TString processTitle[nProcess] = { "tW", "t/s-channel", "t#bar{t}", "WW", "Z/#gamma*+jets", "W+jets",   "data", "other","tW DS","mc"};
  Color_t color[nProcess] =        { kWhite, kMagenta-10, kRed+1, kYellow-10, kAzure-2, kGreen-3,  kBlack,40,kGreen,kYellow};
  
  TString FileNameProcess[3] =  { "twdr", "other", "tt"};
  TString FileNameSystematic[10] = {"JER", "JES", "BtagSF", "PDF", "ZSF", "Matching", "Q2", "TopMass", "UnclusteredMET", "DRDS"};
  TString FileNameRegion[3] = {"1j1t", "2j1", "2j2t"};

  TString cutLabel[nPlots] =     {"cuts","nvertex_1j1t","met_before_1j1t","met_before_2j1t","met_before_2j2t","met_after_1j1t","met_after_2j1t","met_after_2j2t","mll_before","mll_after","mll_after2","ht_before_1j1t","ht_before_2j1t","ht_before_2j2t","ptsys_before_1j1t","ptsys_before_2j1t","ptsys_before_2j2t","ptsys_1j1t","ptsys_2j1t","ptsys_2j2t","ht_1j1t","ht_2j1t","ht_2j2t","njets_final","njetsbt_final","njetsbt_final_loose", "met_zgamma", "mll_zgamma"};
  int rebinHisto[nPlots] =       { 1,1, 4, 4,4,4 ,4,4,4,4, 4,12, 12,12,4,4,4,4,4,4,12,12,12,1,1,1,2,2};
  TString cutTitle[nPlots] =     { "Analysis Cut","# vertices after final", "E_{T}^{miss} before Cut 1j1t", "E_{T}^{miss} before Cut 2j1t","E_{T}^{miss} before Cut 2j2t","E_{T}^{miss} after Cut 1j1t","E_{T}^{miss} after Cut 2j1t","E_{T}^{miss} after Cut 2j2t","Inv. Mass Lept. before Cut", "Inv. Mass Lept. after Cut", "Inv. Mass Lept. after Cut 2", "H_{T} [GeV] before Cut 1j1t","H_{T} [GeV] before Cut 2j1t","H_{T} [GeV] before Cut 2j2t","P_{T} [GeV] of the system before Cut 1j1t","P_{T} [GeV] of the system before Cut 2j1t","P_{T} [GeV] of the  system before Cut 2j2t","P_{T} [GeV] of the system 1j1t","P_{T}  [GeV] of the system 2j1t","P_{T} [GeV] of the system 2j2t","H_{T} [GeV] of the system 1j1t","H_{T}  [GeV] of the  system 2j1t","H_{T} [GeV] of the system 2j2t","Number of jets after full cutflow", "number of tight bjets after full cutflow", "number of loose bjets after full cutflow", "E_{T}^{miss} in Z mass window", "Inv. Mass in Z mass window"};

  TString modeString[3] = {"0", "1", "2"};
  
  TString plotExtension = "plot_graph_";
  
  /////////////////////////////////
  //  MAKE GRAPH initializers 
  /////////////////////////////////////
  TH1F* histo_tt;
  TH1F* histo_tt_ZSFsysDown;
  TH1F* histo_tt_ZSFsysUp;
  TH1F* histo_tt_JERsysDown;
  TH1F* histo_tt_JERsysUp;
  TH1F* histo_tt_JESsysDown;
  TH1F* histo_tt_JESsysUp;
  TH1F* histo_tt_PUsysDown;
  TH1F* histo_tt_PUsysUp;
  TH1F* histo_tt_SFsysDown;
  TH1F* histo_tt_SFsysUp;
  TH1F* histo_tt_METsysDown;
  TH1F* histo_tt_METsysUp;
  TH1F* histo_tt_TopMassDown;
  TH1F* histo_tt_TopMassUp;
  TH1F* histo_tt_Q2Down;
  TH1F* histo_tt_Q2Up;
  TH1F* histo_tt_eleSFsysDown;
  TH1F* histo_tt_eleSFsysUp;
  TH1F* histo_tt_matchingDown;
  TH1F* histo_tt_matchingUp;
  TH1F* histo_tt_PDFDown_1j1t;
  TH1F* histo_tt_PDFUp_1j1t;
  TH1F* histo_tt_PDFDown_2j1t;
  TH1F* histo_tt_PDFUp_2j1t;
  TH1F* histo_tt_PDFDown_2j2t;
  TH1F* histo_tt_PDFUp_2j2t;
  
  TH1F* histo_other;
  TH1F* histo_other_ZSFsysDown;
  TH1F* histo_other_ZSFsysUp;
  TH1F* histo_other_JERsysDown;
  TH1F* histo_other_JERsysUp;
  TH1F* histo_other_JESsysDown;
  TH1F* histo_other_JESsysUp;
  TH1F* histo_other_PUsysDown;
  TH1F* histo_other_PUsysUp;
  TH1F* histo_other_SFsysDown;
  TH1F* histo_other_SFsysUp;
  TH1F* histo_other_METsysDown;
  TH1F* histo_other_METsysUp;
  TH1F* histo_other_eleSFsysDown;
  TH1F* histo_other_eleSFsysUp;
  TH1F* histo_other_PDFDown_1j1t;
  TH1F* histo_other_PDFUp_1j1t;
  TH1F* histo_other_PDFDown_2j1t;
  TH1F* histo_other_PDFUp_2j1t;
  TH1F* histo_other_PDFDown_2j2t;
  TH1F* histo_other_PDFUp_2j2t;
  
  
  TH1F* histo_twdr;
  TH1F* histo_twdr_ZSFsysDown; 
  TH1F* histo_twdr_ZSFsysUp; 
  TH1F* histo_twdr_JERsysDown; 
  TH1F* histo_twdr_JERsysUp; 
  TH1F* histo_twdr_JESsysDown; 
  TH1F* histo_twdr_JESsysUp;
  TH1F* histo_twdr_PUsysDown; 
  TH1F* histo_twdr_PUsysUp;
  TH1F* histo_twdr_SFsysDown; 
  TH1F* histo_twdr_SFsysUp;
  TH1F* histo_twdr_METsysDown; 
  TH1F* histo_twdr_METsysUp;
  TH1F* histo_twdr_TopMassDown; 
  TH1F* histo_twdr_TopMassUp;
  TH1F* histo_twdr_Q2Down; 
  TH1F* histo_twdr_Q2Up;
  TH1F* histo_twdr_eleSFsysDown; 
  TH1F* histo_twdr_eleSFsysUp;
  TH1F* histo_twds; 
  TH1F* histo_twdr_PDFDown_1j1t;
  TH1F* histo_twdr_PDFUp_1j1t;
  TH1F* histo_twdr_PDFDown_2j1t;
  TH1F* histo_twdr_PDFUp_2j1t;
  TH1F* histo_twdr_PDFDown_2j2t;
  TH1F* histo_twdr_PDFUp_2j2t;
   
  
  
  
  
 

  
  //////////
  /// Make stack plot initializers 
  ///////////////////////////////
  TH1F*  h[nPlots][nProcess];
  TH1F*  hextra[nPlots];
  THStack* hStack[nPlots];
  TGraphAsymmErrors *GE[nPlots];
  TGraphAsymmErrors *graph[nPlots];
  
  
  ////////////////////////////
  ///// Loop 
  ////////////////////////////
  for ( int iVariable = 17; iVariable < 18; iVariable++){
  double graph_twdr_y = 0.; 
  double graph_twdr_ZSFsysDown_y= 0.0; 
  double graph_twdr_ZSFsysUp_y  = 0.0;
  double graph_twdr_JERsysDown_y = 0.0;
  double graph_twdr_JERsysUp_y  = 0.0;
  double graph_twdr_JESsysDown_y = 0.0;
  double graph_twdr_JESsysUp_y = 0.0;
  double graph_twdr_PUsysDown_y = 0.0;
  double graph_twdr_PUsysUp_y  = 0.0;
  double graph_twdr_SFsysDown_y = 0.0;
  double graph_twdr_SFsysUp_y = 0.0;
  double graph_twdr_METsysDown_y= 0.0;
  double graph_twdr_METsysUp_y  = 0.0;
  double graph_twdr_eleSFsysDown_y = 0.0;
  double graph_twdr_eleSFsysUp_y = 0.0;
  double graph_twdr_Q2Down_y = 0.0;
  double graph_twdr_Q2Up_y = 0.0;
  double graph_twdr_TopMassDown_y = 0.0;
  double graph_twdr_TopMassUp_y = 0.0;
  double graph_twdr_twdsDown_y = 0.0;
  double graph_twdr_twdsUp_y = 0.0;
  double graph_twdr_errorQ2_y_up = 0.0;
  double graph_twdr_errorQ2_y_down= 0.0;
  double graph_twdr_errortwds_y_up = 0.0;
  double graph_twdr_errortwds_y_down= 0.0;
  double graph_twdr_errorTopMass_y_up = 0.0;
  double graph_twdr_errorTopMass_y_down= 0.0;
  double graph_twdr_errorMET_y_up = 0.0;
  double graph_twdr_errorMET_y_down= 0.0; 
  double graph_twdr_errorZSF_y_up = 0.0;
  double graph_twdr_errorZSF_y_down= 0.0; 
  double graph_twdr_errorJER_y_up = 0.0;
  double graph_twdr_errorJER_y_down= 0.0; 
  double graph_twdr_errorJES_y_up = 0.0;
  double graph_twdr_errorJES_y_down= 0.0; 
  double graph_twdr_errorSF_y_up = 0.0;
  double graph_twdr_errorSF_y_down= 0.0; 
  double graph_twdr_erroreleSF_y_up = 0.0;
  double graph_twdr_erroreleSF_y_down= 0.0; 
  double graph_twdr_errorPU_y_up = 0.0;
  double graph_twdr_errorPU_y_down= 0.0; 
  double graph_twdr_error_up = 0.0;
  double graph_twdr_error_down= 0.0; 
  double graph_twdr_errorPDFDown_1j1t=0.0;
  double graph_twdr_errorPDFUp_1j1t=0.0;
  double graph_twdr_errorPDFDown_2j1t=0.0;
  double graph_twdr_errorPDFUp_2j1t=0.0;
  double graph_twdr_errorPDFDown_2j2t=0.0;
  double graph_twdr_errorPDFUp_2j2t=0.0;
  double graph_twdr_PDFDown_1j1t=0.0;
  double graph_twdr_PDFUp_1j1t=0.0;
  double graph_twdr_PDFDown_2j1t=0.0;
  double graph_twdr_PDFUp_2j1t=0.0;
  double graph_twdr_PDFDown_2j2t=0.0;
  double graph_twdr_PDFUp_2j2t=0.0;

  double graph_other_y = 0.; 
  double graph_other_ZSFsysDown_y= 0.0; 
  double graph_other_ZSFsysUp_y  = 0.0;
  double graph_other_JERsysDown_y = 0.0;
  double graph_other_JERsysUp_y  = 0.0;
  double graph_other_JESsysDown_y = 0.0;
  double graph_other_JESsysUp_y = 0.0;
  double graph_other_PUsysDown_y = 0.0;
  double graph_other_PUsysUp_y  = 0.0;
  double graph_other_SFsysDown_y = 0.0;
  double graph_other_SFsysUp_y = 0.0;
  double graph_other_METsysDown_y= 0.0;
  double graph_other_METsysUp_y  = 0.0;
  double graph_other_eleSFsysDown_y = 0.0;
  double graph_other_eleSFsysUp_y = 0.0;
  double graph_other_errorMET_y_up = 0.0;
  double graph_other_errorMET_y_down= 0.0; 
  double graph_other_errorZSF_y_up = 0.0;
  double graph_other_errorZSF_y_down= 0.0; 
  double graph_other_errorJER_y_up = 0.0;
  double graph_other_errorJER_y_down= 0.0; 
  double graph_other_errorJES_y_up = 0.0;
  double graph_other_errorJES_y_down= 0.0; 
  double graph_other_errorSF_y_up = 0.0;
  double graph_other_errorSF_y_down= 0.0; 
  double graph_other_erroreleSF_y_up = 0.0;
  double graph_other_erroreleSF_y_down= 0.0; 
  double graph_other_errorPU_y_up = 0.0;
  double graph_other_errorPU_y_down= 0.0; 
  double graph_other_error_up = 0.0;
  double graph_other_error_down= 0.0; 
  double graph_other_errorPDFDown_1j1t=0.0;
  double graph_other_errorPDFUp_1j1t=0.0;
  double graph_other_errorPDFDown_2j1t=0.0;
  double graph_other_errorPDFUp_2j1t=0.0;
  double graph_other_errorPDFDown_2j2t=0.0;
  double graph_other_errorPDFUp_2j2t=0.0;
  double graph_other_PDFDown_1j1t=0.0;
  double graph_other_PDFUp_1j1t=0.0;
  double graph_other_PDFDown_2j1t=0.0;
  double graph_other_PDFUp_2j1t=0.0;
  double graph_other_PDFDown_2j2t=0.0;
  double graph_other_PDFUp_2j2t=0.0;
  
  double graph_tt_y = 0.; 
  double graph_tt_JERsysDown_y = 0.0;
  double graph_tt_JERsysUp_y  = 0.0;
  double graph_tt_JESsysDown_y = 0.0;
  double graph_tt_JESsysUp_y = 0.0;
  double graph_tt_PUsysDown_y = 0.0;
  double graph_tt_PUsysUp_y  = 0.0;
  double graph_tt_SFsysDown_y = 0.0;
  double graph_tt_SFsysUp_y = 0.0;
  double graph_tt_METsysDown_y= 0.0;
  double graph_tt_METsysUp_y  = 0.0;
  double graph_tt_eleSFsysDown_y = 0.0;
  double graph_tt_eleSFsysUp_y = 0.0;
  double graph_tt_Q2Down_y = 0.0;
  double graph_tt_Q2Up_y = 0.0;
  double graph_tt_TopMassDown_y = 0.0;
  double graph_tt_TopMassUp_y = 0.0;
  double graph_tt_matchingDown_y = 0.0;
  double graph_tt_matchingUp_y = 0.0;
  double graph_tt_errorQ2_y_up = 0.0;
  double graph_tt_errorQ2_y_down= 0.0;
  double graph_tt_errormatching_y_up = 0.0;
  double graph_tt_errormatching_y_down= 0.0;
  double graph_tt_errorTopMass_y_up = 0.0;
  double graph_tt_errorTopMass_y_down= 0.0;
  double graph_tt_errorMET_y_up = 0.0;
  double graph_tt_errorMET_y_down= 0.0; 
  double graph_tt_errorJER_y_up = 0.0;
  double graph_tt_errorJER_y_down= 0.0; 
  double graph_tt_errorJES_y_up = 0.0;
  double graph_tt_errorJES_y_down= 0.0; 
  double graph_tt_errorSF_y_up = 0.0;
  double graph_tt_errorSF_y_down= 0.0; 
  double graph_tt_erroreleSF_y_up = 0.0;
  double graph_tt_erroreleSF_y_down= 0.0; 
  double graph_tt_errorPU_y_up = 0.0;
  double graph_tt_errorPU_y_down= 0.0; 
  double graph_tt_error_up = 0.0;
  double graph_tt_error_down= 0.0;
  double graph_tt_errorPDFDown_1j1t=0.0;
  double graph_tt_errorPDFUp_1j1t=0.0;
  double graph_tt_errorPDFDown_2j1t=0.0;
  double graph_tt_errorPDFUp_2j1t=0.0;
  double graph_tt_errorPDFDown_2j2t=0.0;
  double graph_tt_errorPDFUp_2j2t=0.0;
  double graph_tt_PDFDown_1j1t=0.0;
  double graph_tt_PDFUp_1j1t=0.0;
  double graph_tt_PDFDown_2j1t=0.0;
  double graph_tt_PDFUp_2j1t=0.0;
  double graph_tt_PDFDown_2j2t=0.0;
  double graph_tt_PDFUp_2j2t=0.0;
  

  

  
  
  
  
  
  
  
    TLegend *leg = new TLegend(0.7,0.7,0.94,0.94);
    leg ->SetFillStyle(1001);
    leg ->SetFillColor(kWhite);
    leg ->SetBorderSize(1);
    hStack[iVariable] = new THStack(cutLabel[iVariable],cutLabel[iVariable]);
    for (int iProcess = 0; iProcess < (nProcess-3); iProcess++){
      h[iVariable][iProcess] = (TH1F*) _file0->Get(cutLabel[iVariable]+ "_" + processName[iProcess]);
      h[iVariable][iProcess]->Rebin(rebinHisto[iVariable]);
      h[iVariable][iProcess]->SetFillColor(color[iProcess]);
      h[iVariable][iProcess]->SetLineColor(kBlack);
      h[iVariable][iProcess]->SetLineWidth(1);
   
      
    }
   
    
    h[iVariable][5]->Add(h[iVariable][1]);
    h[iVariable][5]->Add(h[iVariable][3]); // --> other
   
    
    
    hStack[iVariable]->Add(h[iVariable][5]);
    hStack[iVariable]->Add(h[iVariable][4]);
    hStack[iVariable]->Add(h[iVariable][2]);
    hStack[iVariable]->Add(h[iVariable][0]);
    
  //  if (mode == 0) leg->AddEntry(h[iVariable][6],  processTitle[6], "p");
  //  if (mode == 1) leg->AddEntry(h[iVariable][6],  processTitle[6], "p");
  //  if (mode == 2) leg->AddEntry(h[iVariable][6], processTitle[6], "p");
    
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
    
    //if (iVariable == 5) hStack[iVariable]->GetYaxis()->SetRangeUser(1,100);
   // h[iVariable][6]->Draw("e, sames");

    
    leg->Draw();
    labelcms->Draw();
    labelcms2->Draw();
    
  //  c1->SaveAs("plots/graph/" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + ".png");
   // c1->SaveAs("plots/pdf/" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + ".pdf");
    
    c1->SetLogy();
    hStack[iVariable]->SetMaximum(max * 10);
   // c1->SaveAs("plots/graph/" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + "_log.png");
   // c1->SaveAs("plots/pdf/" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + "_log.pdf");
    
    
    c1->SetLogy(0);
   
   
   
   
 
	//////////////////////////
       // Get histos tt bar 
       //////////////////////////////////
            
       histo_tt= (TH1F*) _file0->Get(cutLabel[iVariable]+ "_" + processName[2]);
	      // cout << "Normal: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 
       histo_tt_JERsysDown= (TH1F*) _fileJERsysDown->Get(cutLabel[iVariable]+ "_" + processName[2]);
	      // cout << "JERsysDown: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 	       
       histo_tt_JERsysUp= (TH1F*) _fileJERsysUp->Get(cutLabel[iVariable]+ "_" + processName[2]);
	       //cout << "JERsysUp: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 
       histo_tt_JESsysDown= (TH1F*) _fileJESsysDown->Get(cutLabel[iVariable]+ "_" + processName[2]);
	     //  cout << "JESsysDown: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 	       
       histo_tt_JESsysUp= (TH1F*) _fileJESsysUp->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "JESsysUp: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 
       histo_tt_PUsysDown= (TH1F*) _filePUsysDown->Get(cutLabel[iVariable]+ "_" + processName[2]);
	      // cout << "PUsysDown: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 	       
       histo_tt_PUsysUp= (TH1F*) _filePUsysUp->Get(cutLabel[iVariable]+ "_" + processName[2]);
	      // cout << "PUsysUp: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 
       histo_tt_SFsysDown= (TH1F*) _fileSFsysDown->Get(cutLabel[iVariable]+ "_" + processName[2]);
	     //  cout << "SFsysDown: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 	       
       histo_tt_SFsysUp= (TH1F*) _fileSFsysUp->Get(cutLabel[iVariable]+ "_" + processName[2]);
	      // cout << "SFsysUp: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 
       histo_tt_METsysDown= (TH1F*) _fileMETsysDown->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "METsysDown: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 	       
       histo_tt_METsysUp= (TH1F*) _fileMETsysUp->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "METsysUp: " << cutLabel[iVariable]+ "_" + processName[2] << endl;
       histo_tt_TopMassDown= (TH1F*) _fileTopMassDown->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "TopMassDown: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 	       
       histo_tt_TopMassUp= (TH1F*) _fileTopMassUp->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "TopMassUp: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 
	histo_tt_Q2Down= (TH1F*) _fileQ2Down->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "Q2Down: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 	       
       histo_tt_Q2Up= (TH1F*) _fileQ2Up->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "Q2Up: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 
	histo_tt_eleSFsysDown= (TH1F*) _fileeleSFsysDown->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "eleSFsysDown: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 	       
       histo_tt_eleSFsysUp= (TH1F*) _fileeleSFsysUp->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "eleSFsysUp: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 
	histo_tt_matchingDown= (TH1F*) _filematchingDown->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "matchingDown: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 	       
       histo_tt_matchingUp= (TH1F*) _filematchingUp->Get(cutLabel[iVariable]+ "_" + processName[2]);
	    //   cout << "matchingUp: " << cutLabel[iVariable]+ "_" + processName[2] << endl; 
	histo_tt_PDFDown_1j1t = (TH1F*) _filePDF_tt_1j1t->Get("histo1a");
	histo_tt_PDFUp_1j1t = (TH1F*) _filePDF_tt_1j1t->Get("histo1b");
	histo_tt_PDFDown_2j1t = (TH1F*) _filePDF_tt_2j1t->Get("histo1a");
	histo_tt_PDFUp_2j1t = (TH1F*) _filePDF_tt_2j1t->Get("histo1b");
	histo_tt_PDFDown_2j2t = (TH1F*) _filePDF_tt_2j2t->Get("histo1a");
	histo_tt_PDFUp_2j2t = (TH1F*) _filePDF_tt_2j2t->Get("histo1b");	       
	
      //////////////////////////
       // Get histos other
       //////////////////////////////////
            
       histo_other= (TH1F*) _file0->Get(cutLabel[iVariable]+ "_" + processName[7]);
	      // cout << "Normal: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 
	 histo_other_ZSFsysDown= (TH1F*) _fileZSFsysDown->Get(cutLabel[iVariable]+ "_" + processName[7]);
	     //  cout << "ZSFsysDown: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 	       
       histo_other_ZSFsysUp= (TH1F*) _fileZSFsysUp->Get(cutLabel[iVariable]+ "_" + processName[7]);
	    //   cout << "ZSFsysUp: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 
       histo_other_JERsysDown= (TH1F*) _fileJERsysDown->Get(cutLabel[iVariable]+ "_" + processName[7]);
	      // cout << "JERsysDown: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 	       
       histo_other_JERsysUp= (TH1F*) _fileJERsysUp->Get(cutLabel[iVariable]+ "_" + processName[7]);
	       //cout << "JERsysUp: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 
       histo_other_JESsysDown= (TH1F*) _fileJESsysDown->Get(cutLabel[iVariable]+ "_" + processName[7]);
	     //  cout << "JESsysDown: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 	       
       histo_other_JESsysUp= (TH1F*) _fileJESsysUp->Get(cutLabel[iVariable]+ "_" + processName[7]);
	    //   cout << "JESsysUp: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 
       histo_other_PUsysDown= (TH1F*) _filePUsysDown->Get(cutLabel[iVariable]+ "_" + processName[7]);
	      // cout << "PUsysDown: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 	       
       histo_other_PUsysUp= (TH1F*) _filePUsysUp->Get(cutLabel[iVariable]+ "_" + processName[7]);
	      // cout << "PUsysUp: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 
       histo_other_SFsysDown= (TH1F*) _fileSFsysDown->Get(cutLabel[iVariable]+ "_" + processName[7]);
	     //  cout << "SFsysDown: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 	       
       histo_other_SFsysUp= (TH1F*) _fileSFsysUp->Get(cutLabel[iVariable]+ "_" + processName[7]);
	      // cout << "SFsysUp: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 
       histo_other_METsysDown= (TH1F*) _fileMETsysDown->Get(cutLabel[iVariable]+ "_" + processName[7]);
	    //   cout << "METsysDown: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 	       
       histo_other_METsysUp= (TH1F*) _fileMETsysUp->Get(cutLabel[iVariable]+ "_" + processName[7]);
	    //   cout << "METsysUp: " << cutLabel[iVariable]+ "_" + processName[7] << endl;
	histo_other_eleSFsysDown= (TH1F*) _fileeleSFsysDown->Get(cutLabel[iVariable]+ "_" + processName[7]);
	    //   cout << "eleSFsysDown: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 	       
       histo_other_eleSFsysUp= (TH1F*) _fileeleSFsysUp->Get(cutLabel[iVariable]+ "_" + processName[7]);
	    //   cout << "eleSFsysUp: " << cutLabel[iVariable]+ "_" + processName[7] << endl; 
	histo_other_PDFDown_1j1t = (TH1F*) _filePDF_other_1j1t->Get("histo1a");
	histo_other_PDFUp_1j1t = (TH1F*) _filePDF_other_1j1t->Get("histo1b");
	histo_other_PDFDown_2j1t = (TH1F*) _filePDF_other_2j1t->Get("histo1a");
	histo_other_PDFUp_2j1t = (TH1F*) _filePDF_other_2j1t->Get("histo1b");
	histo_other_PDFDown_2j2t = (TH1F*) _filePDF_other_2j2t->Get("histo1a");
	histo_other_PDFUp_2j2t = (TH1F*) _filePDF_other_2j2t->Get("histo1b");	     
       
	//////////////////////////
       // Get histos twdr
       //////////////////////////////////
            
       histo_twdr= (TH1F*) _file0->Get(cutLabel[iVariable]+ "_" + processName[0]);
	      // cout << "Normal: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 
       histo_twdr_JERsysDown= (TH1F*) _fileJERsysDown->Get(cutLabel[iVariable]+ "_" + processName[0]);
	      // cout << "JERsysDown: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 	       
       histo_twdr_JERsysUp= (TH1F*) _fileJERsysUp->Get(cutLabel[iVariable]+ "_" + processName[0]);
	       //cout << "JERsysUp: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 
       histo_twdr_JESsysDown= (TH1F*) _fileJESsysDown->Get(cutLabel[iVariable]+ "_" + processName[0]);
	     //  cout << "JESsysDown: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 	       
       histo_twdr_JESsysUp= (TH1F*) _fileJESsysUp->Get(cutLabel[iVariable]+ "_" + processName[0]);
	    //   cout << "JESsysUp: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 
       histo_twdr_PUsysDown= (TH1F*) _filePUsysDown->Get(cutLabel[iVariable]+ "_" + processName[0]);
	      // cout << "PUsysDown: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 	       
       histo_twdr_PUsysUp= (TH1F*) _filePUsysUp->Get(cutLabel[iVariable]+ "_" + processName[0]);
	      // cout << "PUsysUp: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 
       histo_twdr_SFsysDown= (TH1F*) _fileSFsysDown->Get(cutLabel[iVariable]+ "_" + processName[0]);
	     //  cout << "SFsysDown: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 	       
       histo_twdr_SFsysUp= (TH1F*) _fileSFsysUp->Get(cutLabel[iVariable]+ "_" + processName[0]);
	      // cout << "SFsysUp: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 
       histo_twdr_METsysDown= (TH1F*) _fileMETsysDown->Get(cutLabel[iVariable]+ "_" + processName[0]);
	    //   cout << "METsysDown: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 	       
       histo_twdr_METsysUp= (TH1F*) _fileMETsysUp->Get(cutLabel[iVariable]+ "_" + processName[0]);
	    //   cout << "METsysUp: " << cutLabel[iVariable]+ "_" + processName[0] << endl;
       histo_twdr_TopMassDown= (TH1F*) _fileTopMassDown->Get(cutLabel[iVariable]+ "_" + processName[0]);
	    //   cout << "TopMaGEssDown: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 	       
       histo_twdr_TopMassUp= (TH1F*) _fileTopMassUp->Get(cutLabel[iVariable]+ "_" + processName[0]);
	    //   cout << "TopMassUp: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 
	histo_twdr_Q2Down= (TH1F*) _fileQ2Down->Get(cutLabel[iVariable]+ "_" + processName[0]);
	    //   cout << "Q2Down: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 	       
       histo_twdr_Q2Up= (TH1F*) _fileQ2Up->Get(cutLabel[iVariable]+ "_" + processName[0]);
	    //   cout << "Q2Up: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 
	histo_twdr_eleSFsysDown= (TH1F*) _fileeleSFsysDown->Get(cutLabel[iVariable]+ "_" + processName[0]);
	    //   cout << "eleSFsysDown: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 	       
       histo_twdr_eleSFsysUp= (TH1F*) _fileeleSFsysUp->Get(cutLabel[iVariable]+ "_" + processName[0]);
	    //   cout << "eleSFsysUp: " << cutLabel[iVariable]+ "_" + processName[0] << endl; 
	histo_twdr_PDFDown_1j1t = (TH1F*) _filePDF_twdr_1j1t->Get("histo1a");
	histo_twdr_PDFUp_1j1t = (TH1F*) _filePDF_twdr_1j1t->Get("histo1b");
	histo_twdr_PDFDown_2j1t = (TH1F*) _filePDF_twdr_2j1t->Get("histo1a");
	histo_twdr_PDFUp_2j1t = (TH1F*) _filePDF_twdr_2j1t->Get("histo1b");
	histo_twdr_PDFDown_2j2t = (TH1F*) _filePDF_twdr_2j2t->Get("histo1a");
	histo_twdr_PDFUp_2j2t = (TH1F*) _filePDF_twdr_2j2t->Get("histo1b");
       

   
   
   
   
   
   
   
   
   
  hextra[iVariable] = (TH1F*) h[iVariable][5]->Clone();
  hextra[iVariable]->Add(h[iVariable][4]);
  hextra[iVariable]->Add(h[iVariable][2]);
 // hextra[iVariable]->Add(h[iVariable][0]);
  
 //   TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
   // setex2->Draw();
 //   TGraphAsymmErrors *GE[iVariable] = (TGraphAsymmErrors*) hextra[iVariable];





      graph[iVariable] = (TGraphAsymmErrors*) hextra[iVariable]; // plot graph on place where ttbar stack is, so on top of ttbar/other/zjets
	
         	
      
       
       for(int i_bin = 0; i_bin < (histo_tt->GetNbinsX()) ; i_bin++){
              Double_t error = hextra[iVariable]->GetBinError(i_bin); // get stat error
	      graph_tt_y = histo_tt-> GetBinContent(i_bin);
	      graph_tt_JERsysDown_y = histo_tt_JERsysDown->GetBinContent(i_bin);
              graph_tt_JERsysUp_y = histo_tt_JERsysUp->GetBinContent(i_bin); 
	      graph_tt_JESsysDown_y = histo_tt_JESsysDown->GetBinContent(i_bin);
              graph_tt_JESsysUp_y = histo_tt_JESsysUp->GetBinContent(i_bin); 
	      graph_tt_PUsysDown_y = histo_tt_PUsysDown->GetBinContent(i_bin);
              graph_tt_PUsysUp_y = histo_tt_PUsysUp->GetBinContent(i_bin); 
	      graph_tt_SFsysDown_y = histo_tt_SFsysDown->GetBinContent(i_bin);
              graph_tt_SFsysUp_y = histo_tt_SFsysUp->GetBinContent(i_bin); 
	      graph_tt_METsysDown_y = histo_tt_METsysDown->GetBinContent(i_bin);
              graph_tt_METsysUp_y = histo_tt_METsysUp->GetBinContent(i_bin); 
	      graph_tt_eleSFsysDown_y = histo_tt_eleSFsysDown->GetBinContent(i_bin);
              graph_tt_eleSFsysUp_y = histo_tt_eleSFsysUp->GetBinContent(i_bin); 
	      graph_tt_Q2Down_y = histo_tt_Q2Down->GetBinContent(i_bin);
              graph_tt_Q2Up_y = histo_tt_Q2Up->GetBinContent(i_bin);
	      graph_tt_TopMassDown_y = histo_tt_TopMassDown->GetBinContent(i_bin);
              graph_tt_TopMassUp_y = histo_tt_TopMassUp->GetBinContent(i_bin);
	      graph_tt_matchingDown_y = histo_tt_matchingDown->GetBinContent(i_bin);
              graph_tt_matchingUp_y = histo_tt_matchingUp->GetBinContent(i_bin);
	      graph_tt_PDFDown_1j1t = histo_tt_PDFDown_1j1t->GetBinContent(i_bin);
              graph_tt_PDFUp_1j1t = histo_tt_PDFUp_1j1t->GetBinContent(i_bin);
	      
	      
	      
	      graph_twdr_y = histo_twdr-> GetBinContent(i_bin);
	      graph_twdr_JERsysDown_y = histo_twdr_JERsysDown->GetBinContent(i_bin);
              graph_twdr_JERsysUp_y = histo_twdr_JERsysUp->GetBinContent(i_bin); 
	      graph_twdr_JESsysDown_y = histo_twdr_JESsysDown->GetBinContent(i_bin);
              graph_twdr_JESsysUp_y = histo_twdr_JESsysUp->GetBinContent(i_bin); 
	      graph_twdr_PUsysDown_y = histo_twdr_PUsysDown->GetBinContent(i_bin);
              graph_twdr_PUsysUp_y = histo_twdr_PUsysUp->GetBinContent(i_bin); 
	      graph_twdr_SFsysDown_y = histo_twdr_SFsysDown->GetBinContent(i_bin);
              graph_twdr_SFsysUp_y = histo_twdr_SFsysUp->GetBinContent(i_bin); 
	      graph_twdr_METsysDown_y = histo_twdr_METsysDown->GetBinContent(i_bin);
              graph_twdr_METsysUp_y = histo_twdr_METsysUp->GetBinContent(i_bin); 
	      graph_twdr_eleSFsysDown_y = histo_twdr_eleSFsysDown->GetBinContent(i_bin);
              graph_twdr_eleSFsysUp_y = histo_twdr_eleSFsysUp->GetBinContent(i_bin); 
	      graph_twdr_Q2Down_y = histo_twdr_Q2Down->GetBinContent(i_bin);
              graph_twdr_Q2Up_y = histo_twdr_Q2Up->GetBinContent(i_bin);
	      graph_twdr_TopMassDown_y = histo_twdr_TopMassDown->GetBinContent(i_bin);
              graph_twdr_TopMassUp_y = histo_twdr_TopMassUp->GetBinContent(i_bin);
	     graph_twdr_PDFDown_1j1t = histo_twdr_PDFDown_1j1t->GetBinContent(i_bin);
              graph_twdr_PDFUp_1j1t = histo_twdr_PDFUp_1j1t->GetBinContent(i_bin);
	      
	      graph_other_y = histo_other-> GetBinContent(i_bin);
	      graph_other_ZSFsysDown_y = histo_other_ZSFsysDown->GetBinContent(i_bin);
              graph_other_ZSFsysUp_y = histo_other_ZSFsysUp->GetBinContent(i_bin); 
	      graph_other_JERsysDown_y = histo_other_JERsysDown->GetBinContent(i_bin);
              graph_other_JERsysUp_y = histo_other_JERsysUp->GetBinContent(i_bin); 
	      graph_other_JESsysDown_y = histo_other_JESsysDown->GetBinContent(i_bin);
              graph_other_JESsysUp_y = histo_other_JESsysUp->GetBinContent(i_bin); 
	      graph_other_PUsysDown_y = histo_other_PUsysDown->GetBinContent(i_bin);
              graph_other_PUsysUp_y = histo_other_PUsysUp->GetBinContent(i_bin); 
	      graph_other_SFsysDown_y = histo_other_SFsysDown->GetBinContent(i_bin);
              graph_other_SFsysUp_y = histo_other_SFsysUp->GetBinContent(i_bin); 
	      graph_other_METsysDown_y = histo_other_METsysDown->GetBinContent(i_bin);
              graph_other_METsysUp_y = histo_other_METsysUp->GetBinContent(i_bin); 
	      graph_other_eleSFsysDown_y = histo_other_eleSFsysDown->GetBinContent(i_bin);
              graph_other_eleSFsysUp_y = histo_other_eleSFsysUp->GetBinContent(i_bin); 
	      graph_other_PDFDown_1j1t = histo_other_PDFDown_1j1t->GetBinContent(i_bin);
              graph_other_PDFUp_1j1t = histo_other_PDFUp_1j1t->GetBinContent(i_bin);
	      
	      if(graph_tt_PDFDown_1j1t < graph_tt_PDFUp_1j1t){
	           graph_tt_errorPDFUp_1j1t = fabs(graph_tt_PDFUp_1j1t - graph_tt_y);
		   graph_tt_errorPDFDown_1j1t = fabs(graph_tt_PDFDown_1j1t - graph_tt_y);
              }else if (graph_tt_PDFDown_1j1t < graph_tt_PDFUp_1j1t){
	           graph_tt_errorPDFUp_1j1t = fabs(graph_tt_PDFDown_1j1t - graph_tt_y);
		   graph_tt_errorPDFDown_1j1t = fabs(graph_tt_PDFUp_1j1t - graph_tt_y);
	      }
	      cout << "graph_tt_errorPDFUp_1j1t: " << graph_tt_errorPDFUp_1j1t<< endl; 
	      cout << "graph_tt_errorPDFDown_1j1t: " << graph_tt_errorPDFDown_1j1t << endl;
	      if(graph_tt_METsysDown_y < graph_tt_METsysUp_y){
	           graph_tt_errorMET_y_up = fabs(graph_tt_METsysUp_y - graph_tt_y);
		   graph_tt_errorMET_y_down = fabs(graph_tt_METsysDown_y - graph_tt_y);
              }else if (graph_tt_METsysDown_y < graph_tt_METsysUp_y){
	           graph_tt_errorMET_y_up = fabs(graph_tt_METsysDown_y - graph_tt_y);
		   graph_tt_errorMET_y_down = fabs(graph_tt_METsysUp_y - graph_tt_y);
	      }
	     // cout << "graph_tt_errorMET_y_up: " << graph_tt_errorMET_y_up << endl; 
	     // cout << "graph_tt_errorMET_y_down: " << graph_tt_errorMET_y_down << endl;
	     if(graph_tt_JERsysDown_y < graph_tt_JERsysUp_y){
	           graph_tt_errorJER_y_up = fabs(graph_tt_JERsysUp_y - graph_tt_y);
		   graph_tt_errorJER_y_down = fabs(graph_tt_JERsysDown_y - graph_tt_y);
              }else if (graph_tt_JERsysDown_y < graph_tt_JERsysUp_y){
	           graph_tt_errorJER_y_up = fabs(graph_tt_JERsysDown_y - graph_tt_y);
		   graph_tt_errorJER_y_down = fabs(graph_tt_JERsysUp_y - graph_tt_y);
	      }
	    //  cout << "graph_tt_errorJER_y_up: " << graph_tt_errorJER_y_up << endl; 
	     // cout << "graph_tt_errorJER_y_down: " << graph_tt_errorJER_y_down << endl;
	      if(graph_tt_JESsysDown_y < graph_tt_JESsysUp_y){
	           graph_tt_errorJES_y_up = fabs(graph_tt_JESsysUp_y - graph_tt_y);
		   graph_tt_errorJES_y_down = fabs(graph_tt_JESsysDown_y - graph_tt_y);
              }else if (graph_tt_JESsysDown_y < graph_tt_JESsysUp_y){
	           graph_tt_errorJES_y_up = fabs(graph_tt_JESsysDown_y - graph_tt_y);
		   graph_tt_errorJES_y_down = fabs(graph_tt_JESsysUp_y - graph_tt_y);
	      }
	    //  cout << "graph_tt_errorJES_y_up: " << graph_tt_errorJES_y_up << endl; 
	     /// cout << "graph_tt_errorJES_y_down: " << graph_tt_errorJES_y_down << endl;
	      if(graph_tt_PUsysDown_y < graph_tt_PUsysUp_y){
	           graph_tt_errorPU_y_up = fabs(graph_tt_PUsysUp_y - graph_tt_y);
		   graph_tt_errorPU_y_down = fabs(graph_tt_PUsysDown_y - graph_tt_y);
              }else if (graph_tt_PUsysDown_y < graph_tt_PUsysUp_y){
	           graph_tt_errorPU_y_up = fabs(graph_tt_PUsysDown_y - graph_tt_y);
		   graph_tt_errorPU_y_down = fabs(graph_tt_PUsysUp_y - graph_tt_y);
	      }
	     // cout << "graph_tt_errorPU_y_up: " << graph_tt_errorPU_y_up << endl; 
	    //  cout << "graph_tt_errorPU_y_down: " << graph_tt_errorPU_y_down << endl;
	      if(graph_tt_SFsysDown_y < graph_tt_SFsysUp_y){
	           graph_tt_errorSF_y_up = fabs(graph_tt_SFsysUp_y - graph_tt_y);
		   graph_tt_errorSF_y_down = fabs(graph_tt_SFsysDown_y - graph_tt_y);
              }else if (graph_tt_SFsysDown_y < graph_tt_SFsysUp_y){
	           graph_tt_errorSF_y_up = fabs(graph_tt_SFsysDown_y - graph_tt_y);
		   graph_tt_errorSF_y_down = fabs(graph_tt_SFsysUp_y - graph_tt_y);
	      }
	//      cout << "graph_tt_errorSF_y_up: " << graph_tt_errorSF_y_up << endl; 
	//      cout << "graph_tt_errorSF_y_down: " << graph_tt_errorSF_y_down << endl;
	      if(graph_tt_eleSFsysDown_y < graph_tt_eleSFsysUp_y){
	           graph_tt_erroreleSF_y_up = fabs(graph_tt_eleSFsysUp_y - graph_tt_y);
		   graph_tt_erroreleSF_y_down = fabs(graph_tt_eleSFsysDown_y - graph_tt_y);
              }else if (graph_tt_eleSFsysDown_y < graph_tt_eleSFsysUp_y){
	           graph_tt_erroreleSF_y_up = fabs(graph_tt_eleSFsysDown_y - graph_tt_y);
		   graph_tt_erroreleSF_y_down = fabs(graph_tt_eleSFsysUp_y - graph_tt_y);
	      }
	      //cout << "graph_tt_erroreleSF_y_up: " << graph_tt_erroreleSF_y_up << endl; 
	      //cout << "graph_tt_erroreleSF_y_down: " << graph_tt_erroreleSF_y_down << endl;
	      if(graph_tt_Q2Down_y < graph_tt_Q2Up_y){
	           graph_tt_errorQ2_y_up = fabs(graph_tt_Q2Up_y - graph_tt_y);
		   graph_tt_errorQ2_y_down = fabs(graph_tt_Q2Down_y - graph_tt_y);
              }else if (graph_tt_Q2Down_y < graph_tt_Q2Up_y){
	           graph_tt_errorQ2_y_up = fabs(graph_tt_Q2Down_y - graph_tt_y);
		   graph_tt_errorQ2_y_down = fabs(graph_tt_Q2Up_y - graph_tt_y);
	      }
	     // cout << "graph_tt_errorQ2_y_up: " << graph_tt_errorQ2_y_up << endl; 
	     // cout << "graph_tt_errorQ2_y_down: " << graph_tt_errorQ2_y_down << endl;
	       if(graph_tt_matchingDown_y < graph_tt_matchingUp_y){
	           graph_tt_errormatching_y_up = fabs(graph_tt_matchingUp_y - graph_tt_y);
		   graph_tt_errormatching_y_down = fabs(graph_tt_matchingDown_y - graph_tt_y);
              }else if (graph_tt_matchingDown_y < graph_tt_matchingUp_y){
	           graph_tt_errormatching_y_up = fabs(graph_tt_matchingDown_y - graph_tt_y);
		   graph_tt_errormatching_y_down = fabs(graph_tt_matchingUp_y - graph_tt_y);
	      }
	     // cout << "graph_tt_errormatching_y_up: " << graph_tt_errormatching_y_up << endl; 
	     // cout << "graph_tt_errormatching_y_down: " << graph_tt_errormatching_y_down << endl;
	       if(graph_tt_TopMassDown_y < graph_tt_TopMassUp_y){
	           graph_tt_errorTopMass_y_up = fabs(graph_tt_TopMassUp_y - graph_tt_y);
		   graph_tt_errorTopMass_y_down = fabs(graph_tt_TopMassDown_y - graph_tt_y);
              }else if (graph_tt_TopMassDown_y < graph_tt_TopMassUp_y){
	           graph_tt_errorTopMass_y_up = fabs(graph_tt_TopMassDown_y - graph_tt_y);
		   graph_tt_errorTopMass_y_down = fabs(graph_tt_TopMassUp_y - graph_tt_y);
	      }
	     // cout << "graph_tt_errorTopMass_y_up: " << graph_tt_errorTopMass_y_up << endl; 
	    //  cout << "graph_tt_errorTopMass_y_down: " << graph_tt_errorTopMass_y_down << endl;
	    /////////////////////////////////////////////////
	    
	    if(graph_twdr_PDFDown_1j1t < graph_twdr_PDFUp_1j1t){
	           graph_twdr_errorPDFUp_1j1t = fabs(graph_twdr_PDFUp_1j1t - graph_twdr_y);
		   graph_twdr_errorPDFDown_1j1t = fabs(graph_twdr_PDFDown_1j1t - graph_twdr_y);
              }else if (graph_twdr_PDFDown_1j1t < graph_twdr_PDFUp_1j1t){
	           graph_twdr_errorPDFUp_1j1t = fabs(graph_twdr_PDFDown_1j1t - graph_twdr_y);
		   graph_twdr_errorPDFDown_1j1t = fabs(graph_twdr_PDFUp_1j1t - graph_twdr_y);
	      }
	      cout << "graph_twdr_errorPDFUp_1j1t: " << graph_twdr_errorPDFUp_1j1t<< endl; 
	      cout << "graph_twdr_errorPDFDown_1j1t: " << graph_twdr_errorPDFDown_1j1t << endl;
	    
	    if(graph_twdr_METsysDown_y < graph_twdr_METsysUp_y){
	           graph_twdr_errorMET_y_up = fabs(graph_twdr_METsysUp_y - graph_twdr_y);
		   graph_twdr_errorMET_y_down = fabs(graph_twdr_METsysDown_y - graph_twdr_y);
              }else if (graph_twdr_METsysDown_y < graph_twdr_METsysUp_y){
	           graph_twdr_errorMET_y_up = fabs(graph_twdr_METsysDown_y - graph_twdr_y);
		   graph_twdr_errorMET_y_down = fabs(graph_twdr_METsysUp_y - graph_twdr_y);
	      }
	     // cout << "graph_twdr_errorMET_y_up: " << graph_twdr_errorMET_y_up << endl; 
	     // cout << "graph_twdr_errorMET_y_down: " << graph_twdr_errorMET_y_down << endl;
	     if(graph_twdr_JERsysDown_y < graph_twdr_JERsysUp_y){
	           graph_twdr_errorJER_y_up = fabs(graph_twdr_JERsysUp_y - graph_twdr_y);
		   graph_twdr_errorJER_y_down = fabs(graph_twdr_JERsysDown_y - graph_twdr_y);
              }else if (graph_twdr_JERsysDown_y < graph_twdr_JERsysUp_y){
	           graph_twdr_errorJER_y_up = fabs(graph_twdr_JERsysDown_y - graph_twdr_y);
		   graph_twdr_errorJER_y_down = fabs(graph_twdr_JERsysUp_y - graph_twdr_y);
	      }
	    //  cout << "graph_twdr_errorJER_y_up: " << graph_twdr_errorJER_y_up << endl; 
	     // cout << "graph_twdr_errorJER_y_down: " << graph_twdr_errorJER_y_down << endl;
	      if(graph_twdr_JESsysDown_y < graph_twdr_JESsysUp_y){
	           graph_twdr_errorJES_y_up = fabs(graph_twdr_JESsysUp_y - graph_twdr_y);
		   graph_twdr_errorJES_y_down = fabs(graph_twdr_JESsysDown_y - graph_twdr_y);
              }else if (graph_twdr_JESsysDown_y < graph_twdr_JESsysUp_y){
	           graph_twdr_errorJES_y_up = fabs(graph_twdr_JESsysDown_y - graph_twdr_y);
		   graph_twdr_errorJES_y_down = fabs(graph_twdr_JESsysUp_y - graph_twdr_y);
	      }
	    //  cout << "graph_twdr_errorJES_y_up: " << graph_twdr_errorJES_y_up << endl; 
	     // cout << "graph_twdr_errorJES_y_down: " << graph_twdr_errorJES_y_down << endl;
	      if(graph_twdr_PUsysDown_y < graph_twdr_PUsysUp_y){
	           graph_twdr_errorPU_y_up = fabs(graph_twdr_PUsysUp_y - graph_twdr_y);
		   graph_twdr_errorPU_y_down = fabs(graph_twdr_PUsysDown_y - graph_twdr_y);
              }else if (graph_twdr_PUsysDown_y < graph_twdr_PUsysUp_y){
	           graph_twdr_errorPU_y_up = fabs(graph_twdr_PUsysDown_y - graph_twdr_y);
		   graph_twdr_errorPU_y_down = fabs(graph_twdr_PUsysUp_y - graph_twdr_y);
	      }
	     // cout << "graph_twdr_errorPU_y_up: " << graph_twdr_errorPU_y_up << endl; 
	    //  cout << "graph_twdr_errorPU_y_down: " << graph_twdr_errorPU_y_down << endl;
	      if(graph_twdr_SFsysDown_y < graph_twdr_SFsysUp_y){
	           graph_twdr_errorSF_y_up = fabs(graph_twdr_SFsysUp_y - graph_twdr_y);
		   graph_twdr_errorSF_y_down = fabs(graph_twdr_SFsysDown_y - graph_twdr_y);
              }else if (graph_twdr_SFsysDown_y < graph_twdr_SFsysUp_y){
	           graph_twdr_errorSF_y_up = fabs(graph_twdr_SFsysDown_y - graph_twdr_y);
		   graph_twdr_errorSF_y_down = fabs(graph_twdr_SFsysUp_y - graph_twdr_y);
	      }
	//      cout << "graph_twdr_errorSF_y_up: " << graph_twdr_errorSF_y_up << endl; 
	//      cout << "graph_twdr_errorSF_y_down: " << graph_twdr_errorSF_y_down << endl;
	      if(graph_twdr_eleSFsysDown_y < graph_twdr_eleSFsysUp_y){
	           graph_twdr_erroreleSF_y_up = fabs(graph_twdr_eleSFsysUp_y - graph_twdr_y);
		   graph_twdr_erroreleSF_y_down = fabs(graph_twdr_eleSFsysDown_y - graph_twdr_y);
              }else if (graph_twdr_eleSFsysDown_y < graph_twdr_eleSFsysUp_y){
	           graph_twdr_erroreleSF_y_up = fabs(graph_twdr_eleSFsysDown_y - graph_twdr_y);
		   graph_twdr_erroreleSF_y_down = fabs(graph_twdr_eleSFsysUp_y - graph_twdr_y);
	      }
	      //cout << "graph_twdr_erroreleSF_y_up: " << graph_twdr_erroreleSF_y_up << endl; 
	      //cout << "graph_twdr_erroreleSF_y_down: " << graph_twdr_erroreleSF_y_down << endl;
	      if(graph_twdr_Q2Down_y < graph_twdr_Q2Up_y){
	           graph_twdr_errorQ2_y_up = fabs(graph_twdr_Q2Up_y - graph_twdr_y);
		   graph_twdr_errorQ2_y_down = fabs(graph_twdr_Q2Down_y - graph_twdr_y);
              }else if (graph_twdr_Q2Down_y < graph_twdr_Q2Up_y){
	           graph_twdr_errorQ2_y_up = fabs(graph_twdr_Q2Down_y - graph_twdr_y);
		   graph_twdr_errorQ2_y_down = fabs(graph_twdr_Q2Up_y - graph_twdr_y);
	      }
	     // cout << "graph_twdr_errorQ2_y_up: " << graph_twdr_errorQ2_y_up << endl; 
	     // cout << "graph_twdr_errorQ2_y_down: " << graph_twdr_errorQ2_y_down << endl;
	       if(graph_twdr_TopMassDown_y < graph_twdr_TopMassUp_y){
	           graph_twdr_errorTopMass_y_up = fabs(graph_twdr_TopMassUp_y - graph_twdr_y);
		   graph_twdr_errorTopMass_y_down = fabs(graph_twdr_TopMassDown_y - graph_twdr_y);
              }else if (graph_twdr_TopMassDown_y < graph_twdr_TopMassUp_y){
	           graph_twdr_errorTopMass_y_up = fabs(graph_twdr_TopMassDown_y - graph_twdr_y);
		   graph_twdr_errorTopMass_y_down = fabs(graph_twdr_TopMassUp_y - graph_twdr_y);
	      }
	     // cout << "graph_twdr_errorTopMass_y_up: " << graph_twdr_errorTopMass_y_up << endl; 
	    //  cout << "graph_twdr_errorTopMass_y_down: " << graph_twdr_errorTopMass_y_down << endl;
	    
	    
	    
	     /////////////////////////////////////////////
	     
	     if(graph_other_PDFDown_1j1t < graph_other_PDFUp_1j1t){
	           graph_other_errorPDFUp_1j1t = fabs(graph_other_PDFUp_1j1t - graph_other_y);
		   graph_other_errorPDFDown_1j1t = fabs(graph_other_PDFDown_1j1t - graph_other_y);
              }else if (graph_other_PDFDown_1j1t < graph_other_PDFUp_1j1t){
	           graph_other_errorPDFUp_1j1t = fabs(graph_other_PDFDown_1j1t - graph_other_y);
		   graph_other_errorPDFDown_1j1t = fabs(graph_other_PDFUp_1j1t - graph_other_y);
	      }
	      cout << "graph_other_errorPDFUp_1j1t: " << graph_other_errorPDFUp_1j1t<< endl; 
	      cout << "graph_other_errorPDFDown_1j1t: " << graph_other_errorPDFDown_1j1t << endl;
	     
	     
	     if(graph_other_METsysDown_y < graph_other_METsysUp_y){
	           graph_other_errorMET_y_up = fabs(graph_other_METsysUp_y - graph_other_y);
		   graph_other_errorMET_y_down = fabs(graph_other_METsysDown_y - graph_other_y);
              }else if (graph_other_METsysDown_y < graph_other_METsysUp_y){
	           graph_other_errorMET_y_up = fabs(graph_other_METsysDown_y - graph_other_y);
		   graph_other_errorMET_y_down = fabs(graph_other_METsysUp_y - graph_other_y);
	      }
	     // cout << "graph_other_errorMET_y_up: " << graph_other_errorMET_y_up << endl; 
	     // cout << "graph_other_errorMET_y_down: " << graph_other_errorMET_y_down << endl;
	      if(graph_other_ZSFsysDown_y < graph_other_ZSFsysUp_y){
	           graph_other_errorZSF_y_up = fabs(graph_other_ZSFsysUp_y - graph_other_y);
		   graph_other_errorZSF_y_down = fabs(graph_other_ZSFsysDown_y - graph_other_y);
              }else if (graph_other_ZSFsysDown_y < graph_other_ZSFsysUp_y){
	           graph_other_errorZSF_y_up = fabs(graph_other_ZSFsysDown_y - graph_other_y);
		   graph_other_errorZSF_y_down = fabs(graph_other_ZSFsysUp_y - graph_other_y);
	      }
	     // cout << "graph_other_errorZSF_y_up: " << graph_other_errorZSF_y_up << endl; 
	      //cout << "graph_other_errorZSF_y_down: " << graph_other_errorZSF_y_down << endl;
	      if(graph_other_JERsysDown_y < graph_other_JERsysUp_y){
	           graph_other_errorJER_y_up = fabs(graph_other_JERsysUp_y - graph_other_y);
		   graph_other_errorJER_y_down = fabs(graph_other_JERsysDown_y - graph_other_y);
              }else if (graph_other_JERsysDown_y < graph_other_JERsysUp_y){
	           graph_other_errorJER_y_up = fabs(graph_other_JERsysDown_y - graph_other_y);
		   graph_other_errorJER_y_down = fabs(graph_other_JERsysUp_y - graph_other_y);
	      }
	    //  cout << "graph_other_errorJER_y_up: " << graph_other_errorJER_y_up << endl; 
	     // cout << "graph_other_errorJER_y_down: " << graph_other_errorJER_y_down << endl;
	      if(graph_other_JESsysDown_y < graph_other_JESsysUp_y){
	           graph_other_errorJES_y_up = fabs(graph_other_JESsysUp_y - graph_other_y);
		   graph_other_errorJES_y_down = fabs(graph_other_JESsysDown_y - graph_other_y);
              }else if (graph_other_JESsysDown_y < graph_other_JESsysUp_y){
	           graph_other_errorJES_y_up = fabs(graph_other_JESsysDown_y - graph_other_y);
		   graph_other_errorJES_y_down = fabs(graph_other_JESsysUp_y - graph_other_y);
	      }
	    //  cout << "graph_other_errorJES_y_up: " << graph_other_errorJES_y_up << endl; 
	     // cout << "graph_other_errorJES_y_down: " << graph_other_errorJES_y_down << endl;
	      if(graph_other_PUsysDown_y < graph_other_PUsysUp_y){
	           graph_other_errorPU_y_up = fabs(graph_other_PUsysUp_y - graph_other_y);
		   graph_other_errorPU_y_down = fabs(graph_other_PUsysDown_y - graph_other_y);
              }else if (graph_other_PUsysDown_y < graph_other_PUsysUp_y){
	           graph_other_errorPU_y_up = fabs(graph_other_PUsysDown_y - graph_other_y);
		   graph_other_errorPU_y_down = fabs(graph_other_PUsysUp_y - graph_other_y);
	      }
	     // cout << "graph_other_errorPU_y_up: " << graph_other_errorPU_y_up << endl; 
	    //  cout << "graph_other_errorPU_y_down: " << graph_other_errorPU_y_down << endl;
	      if(graph_other_SFsysDown_y < graph_other_SFsysUp_y){
	           graph_other_errorSF_y_up = fabs(graph_other_SFsysUp_y - graph_other_y);
		   graph_other_errorSF_y_down = fabs(graph_other_SFsysDown_y - graph_other_y);
              }else if (graph_other_SFsysDown_y < graph_other_SFsysUp_y){
	           graph_other_errorSF_y_up = fabs(graph_other_SFsysDown_y - graph_other_y);
		   graph_other_errorSF_y_down = fabs(graph_other_SFsysUp_y - graph_other_y);
	      }
	//      cout << "graph_other_errorSF_y_up: " << graph_other_errorSF_y_up << endl; 
	//      cout << "graph_other_errorSF_y_down: " << graph_other_errorSF_y_down << endl;
	      if(graph_other_eleSFsysDown_y < graph_other_eleSFsysUp_y){
	           graph_other_erroreleSF_y_up = fabs(graph_other_eleSFsysUp_y - graph_other_y);
		   graph_other_erroreleSF_y_down = fabs(graph_other_eleSFsysDown_y - graph_other_y);
              }else if (graph_other_eleSFsysDown_y < graph_other_eleSFsysUp_y){
	           graph_other_erroreleSF_y_up = fabs(graph_other_eleSFsysDown_y - graph_other_y);
		   graph_other_erroreleSF_y_down = fabs(graph_other_eleSFsysUp_y - graph_other_y);
	      }
	      //cout << "graph_other_erroreleSF_y_up: " << graph_other_erroreleSF_y_up << endl; 
	      //cout << "graph_other_erroreleSF_y_down: " << graph_other_erroreleSF_y_down << endl;
	     
	    
	      double tt_error_int_up = graph_tt_errorPDFUp_1j1t*graph_tt_errorPDFUp_1j1t+graph_tt_errorQ2_y_up*graph_tt_errorQ2_y_up+graph_tt_errorTopMass_y_up*graph_tt_errorTopMass_y_up+graph_tt_errormatching_y_up*graph_tt_errormatching_y_up+graph_tt_erroreleSF_y_up*graph_tt_erroreleSF_y_up+graph_tt_errorSF_y_up*graph_tt_errorSF_y_up+graph_tt_errorJER_y_up*graph_tt_errorJER_y_up+graph_tt_errorJES_y_up*graph_tt_errorJES_y_up+graph_tt_errorPU_y_up*graph_tt_errorPU_y_up+graph_tt_errorMET_y_up*graph_tt_errorMET_y_up;
	      double twdr_error_int_up = graph_twdr_errorPDFUp_1j1t*graph_twdr_errorPDFUp_1j1t+graph_twdr_errorQ2_y_up*graph_twdr_errorQ2_y_up+graph_twdr_errorTopMass_y_up*graph_twdr_errorTopMass_y_up+graph_twdr_erroreleSF_y_up*graph_twdr_erroreleSF_y_up+graph_twdr_errorSF_y_up*graph_twdr_errorSF_y_up+graph_twdr_errorJER_y_up*graph_twdr_errorJER_y_up+graph_twdr_errorJES_y_up*graph_twdr_errorJES_y_up+graph_twdr_errorPU_y_up*graph_twdr_errorPU_y_up+graph_twdr_errorMET_y_up*graph_twdr_errorMET_y_up;
	      double other_error_int_up = graph_other_errorPDFUp_1j1t*graph_other_errorPDFUp_1j1t+graph_other_erroreleSF_y_up*graph_other_erroreleSF_y_up+graph_other_errorSF_y_up*graph_other_errorSF_y_up+graph_other_errorJER_y_up*graph_other_errorJER_y_up+graph_other_errorJES_y_up*graph_other_errorJES_y_up+graph_other_errorZSF_y_up*graph_other_errorZSF_y_up+graph_other_errorPU_y_up*graph_other_errorPU_y_up+graph_other_errorMET_y_up*graph_other_errorMET_y_up;
	      double tt_error_int_down = graph_tt_errorPDFDown_1j1t*graph_tt_errorPDFDown_1j1t+graph_tt_errorQ2_y_down*graph_tt_errorQ2_y_down+graph_tt_errorTopMass_y_down*graph_tt_errorTopMass_y_down+graph_tt_errormatching_y_down*graph_tt_errormatching_y_down+graph_tt_erroreleSF_y_down*graph_tt_erroreleSF_y_down+graph_tt_errorSF_y_down*graph_tt_errorSF_y_down+graph_tt_errorJER_y_down*graph_tt_errorJER_y_down+graph_tt_errorJES_y_down*graph_tt_errorJES_y_down+graph_tt_errorPU_y_down*graph_tt_errorPU_y_down+graph_tt_errorMET_y_down*graph_tt_errorMET_y_down;
	      double twdr_error_int_down = graph_twdr_errorPDFDown_1j1t*graph_twdr_errorPDFDown_1j1t+graph_twdr_errorQ2_y_down*graph_twdr_errorQ2_y_down+graph_twdr_errorTopMass_y_down*graph_twdr_errorTopMass_y_down+graph_twdr_erroreleSF_y_down*graph_twdr_erroreleSF_y_down+graph_twdr_errorSF_y_down*graph_twdr_errorSF_y_down+graph_twdr_errorJER_y_down*graph_twdr_errorJER_y_down+graph_twdr_errorJES_y_down*graph_twdr_errorJES_y_down+graph_twdr_errorPU_y_down*graph_twdr_errorPU_y_down+graph_twdr_errorMET_y_down*graph_twdr_errorMET_y_down;
	      double other_error_int_down = graph_other_errorPDFDown_1j1t*graph_other_errorPDFDown_1j1t+graph_other_erroreleSF_y_down*graph_other_erroreleSF_y_down+graph_other_errorSF_y_down*graph_other_errorSF_y_down+graph_other_errorJER_y_down*graph_other_errorJER_y_down+graph_other_errorJES_y_down*graph_other_errorJES_y_down+graph_other_errorZSF_y_down*graph_other_errorZSF_y_down+graph_other_errorPU_y_down*graph_other_errorPU_y_down+graph_other_errorMET_y_down*graph_other_errorMET_y_down;
	     
	      graph_tt_error_up = sqrt(error*error + tt_error_int_up+ other_error_int_up );//+ twdr_error_int_up );
		// cout << "yhigh " << graph_tt_error_up << endl; 
              graph_tt_error_down = sqrt(error*error + tt_error_int_down   + other_error_int_down );//+ twdr_error_int_down );
	      //  cout << "ylow " << graph_tt_error_down << endl;
	     
	      graph[iVariable]->SetPointEYhigh(i_bin, graph_tt_error_up);
	      
	      graph[iVariable]->SetPointEYlow(i_bin, graph_tt_error_down);
	     

             
       }
       

















   
    leg->AddEntry(graph[iVariable], "unc stat+sys (no tW)", "f");
    hStack[iVariable]->SetMaximum(max*1.2);
    //GE[iVariable]->Draw("sames, e2"); */
    graph[iVariable]->SetFillColor(1);
    graph[iVariable]->SetFillStyle(3005);
    //graph_tt->SetMarkerSize(0);
    //graph_tt->SetLineWidth(0);
   // graph_tt->SetLineColor(kWhite);
    graph[iVariable]->Draw("sames,E3"); 
   
   // h[iVariable][6]->Draw("e, sames");
    
    c1->SaveAs("plots/graph/error_" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + ".png");
    c1->SaveAs("plots/pdf/error_" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + ".pdf");
     
    c1->SetLogy();
    hStack[iVariable]->SetMaximum(max * 10);
    c1->SaveAs("plots/graph/error_" + plotExtension+ modeString[mode] + "_" + cutLabel[iVariable] + "_log.png");
    c1->SaveAs("plots/pdf/error_" + plotExtension + modeString[mode] + "_" + cutLabel[iVariable] + "_log.pdf");
    
   
  }

}
