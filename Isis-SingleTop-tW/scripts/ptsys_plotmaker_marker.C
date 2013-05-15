//Isis Van Parijs


#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMarker.h"
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

void ptsys_plotmaker(int mode = 0, int region = 0){
 
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
  if (mode == 0) labelcms2->AddText("12 fb^{-1}, e#mu channel  ");
  if (mode == 1) labelcms2->AddText("12 fb^{-1}, #mu#mu channel  ");
  if (mode == 2) labelcms2->AddText("12 fb^{-1}, ee channel  ");
  labelcms2->SetBorderSize(0);
  
  
    alabelcms  = new TPaveText(0.12,0.88,0.5,0.94,"NDCBR");
  alabelcms->SetTextAlign(12);
  alabelcms->SetTextSize(0.045);
  alabelcms->SetFillColor(kWhite);
  alabelcms->AddText("CMS Preliminary, #sqrt{s} = 8 TeV");
  alabelcms->SetBorderSize(0);

  alabelcms2  = new TPaveText(0.12,0.85,0.5,0.88,"NDCBR");
  alabelcms2->SetTextAlign(12);
  alabelcms2->SetTextSize(0.045);
  alabelcms2->SetFillColor(kWhite);
  if (mode == 0) alabelcms2->AddText("12 fb^{-1}, e#mu channel  ");
  if (mode == 1) alabelcms2->AddText("12 fb^{-1}, #mu#mu channel  ");
  if (mode == 2) alabelcms2->AddText("12 fb^{-1}, ee channel  ");
  alabelcms2->SetBorderSize(0);
  
  
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
  
  
  double lumi = 1000;
   if      (mode == 0)	 lumi = 11966.617;  	
    else if (mode == 1) lumi = 12067.294;  	
    else if (mode == 2) lumi = 12093.792;  	
    
  sprintf(myRootFile,"results/an_%dpb_%d.root", lumi, mode); // take output from looper
  
  
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

  
  
  cout << "-------------------------------------------------------" << endl; 
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
  
  const int nProcess = 4;
  const int nPlots = 10; // 16
  const int nSys = 21;


  TString processName[nProcess] =  { "twdr", "tt","twds","others"};
  TString processTitle[nProcess] = { "tW DR","t#bar{t}","tW DS", "other"};
  Color_t color[nProcess] =        { kBlue, kRed, kGreen,kMagenta};
  
  TString systName[nSys] = { "Normal", "ZSFsysDown" , "ZSFsysUp","JERsysDown" , "JERsysUp","JESsysDown" , "JESsysUp","PUsysDown" , "PUsysUp","SFsysDown" , "SFsysUp","METsysDown" , "METsysUp","TopMassDown" , "TopMassUp","Q2Down" , "Q2Up","eleSFsysDown" , "eleSFsysUp","matchingDown" , "matchingUp"}; 

/*
  TString cutLabel[nPlots] =     {  "met", "mll", "ptsys", "ht", "pt_leading", "nvertex", "met_2j1t", "mll_2j1t", "ptsys_2j1t", "ht_2j1t", "pt_leading_2j1t","met_2j2t", "mll_2j2t", "ptsys_2j2t", "ht_2j2t", "pt_leading_2j2t" };
  int rebinHisto[nPlots] =       {  4, 4,  4, 12, 4, 1,4, 4,  4, 12, 4, 4, 4,  4, 12, 4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} ", "Inv. Mass ", "P_{T} system [GeV] ", "H_{T} [GeV] ","P_{T} of the leading jet ", "# of vertex ", "E_{T}^{miss} 2j1t", "Inv. Mass 2j1t", "P_{T} system [GeV] 2j1t", "H_{T} [GeV] 2j1t","P_{T} of the leading jet 2j1t","E_{T}^{miss} 2j2t", "Inv. Mass 2j2t", "P_{T} system [GeV] 2j2t", "H_{T} [GeV] 2j2t","P_{T} of the leading jet 2j2t" };
*/

  if(region == 0){
  TString cutLabel[nPlots] =     { "met_1j1t", "mll_1j1t", "ptsys_1j1t", "ht_1j1t","pt_leading_1j1t",  "nvertex_1j1t", "pt_max_1j1t","pt_min_1j1t","ht_nomet_1j1t","eta_leading_1j1t" };
  int rebinHisto[nPlots] =       { 4, 4,  4, 12, 4, 1,4,4,12,4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} 1j1t", "Inv. Mass 1j1t", "P_{T} system [GeV] 1j1t","H_{T} [GeV] 1j1t","P_{T} of the leading jet 1j1t", "# of vertex 1j1t", "p_T of the first lepton [GeV] 1j1t", "p_T  of the second lepton [GeV] 1j1t", "H_{T} no met [GeV] 1j1t","Eta of the leading jet 1j1t"};
  } else if(region == 1){
   TString cutLabel[nPlots] =     { "met_2j1t", "mll_2j1t", "ptsys_2j1t", "ht_2j1t","pt_leading_2j1t",  "nvertex_2j1t", "pt_max_2j1t","pt_min_2j1t","ht_nomet_2j1t","eta_leading_2j1t" };
  int rebinHisto[nPlots] =       {  4, 4,  4, 12, 4, 1,4,4,12,4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} 2j1t", "Inv. Mass 2j1t", "P_{T} system [GeV] 2j1t","H_{T} [GeV] 2j1t","P_{T} of the leading jet 2j1t", "# of vertex 2j1t", "p_T of the first lepton [GeV] 2j1t", "p_T  of the second lepton [GeV] 2j1t", "H_{T} no met [GeV] 2j1t","Eta of the leading jet 2j1t"};
  }else if(region == 2){
   TString cutLabel[nPlots] =     { "met_2j2t", "mll_2j2t", "ptsys_2j2t", "ht_2j2t","pt_leading_2j2t",  "nvertex_2j2t", "pt_max_2j2t","pt_min_2j2t","ht_nomet_2j2t","eta_leading_2j2t" };
  int rebinHisto[nPlots] =       { 4, 4,  4, 12, 4, 1,4,4,12,4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} 2j2t", "Inv. Mass 2j2t", "P_{T} system [GeV] 2j2t","H_{T} [GeV] 2j2t","P_{T} of the leading jet 2j2t", "# of vertex 2j2t", "p_T of the first lepton [GeV] 2j2t", "p_T  of the second lepton [GeV] 2j2t", "H_{T} no met [GeV] 2j2t","Eta of the leading jet 2j2t"};
 
  
  }

  TString modeString[3] = {"0","1","2"};
  TString regionString[3] = {"0","1","2"};
  
  TString plotExtension = "sys_plot_"; // name of the plots
  TString plotZSF = "ZSF_";
  TString plotJER = "JER_"; 
  TString plotJES = "JES_"; 
  TString plotPU = "PU_";
  TString plotSF = "SF_";
  TString plotMET = "MET_";
  TString plotTopMass = "TopMass_";
  TString plotQ2 = "Q2_";
  TString ploteleSF = "eleSF_";
  TString plotmatching = "matching_";
  TString plotDS = "DS_";
  
  TString plotAnalysis = "Systematics"; // directory in plots where the plots are saved
  

  
   
  
  
  //Make plots   
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
  
  
  
  TGraph* graph_tt;
  TGraph* graph_tt_ZSFsysDown;
  TGraph* graph_tt_ZSFsysUp;
  TGraph* graph_tt_JERsysDown;
  TGraph* graph_tt_JERsysUp;
  TGraph* graph_tt_JESsysDown;
  TGraph* graph_tt_JESsysUp;
  TGraph* graph_tt_PUsysDown;
  TGraph* graph_tt_PUsysUp;
  TGraph* graph_tt_SFsysDown;
  TGraph* graph_tt_SFsysUp;
  TGraph* graph_tt_METsysDown;
  TGraph* graph_tt_METsysUp;
  TGraph* graph_tt_TopMassDown;
  TGraph* graph_tt_TopMassUp;
  TGraph* graph_tt_Q2Down;
  TGraph* graph_tt_Q2Up;
  TGraph* graph_tt_eleSFsysDown;
  TGraph* graph_tt_eleSFsysUp;
  TGraph* graph_tt_matchingDown;
  TGraph* graph_tt_matchingUp;
  
  TGraph* graph_other;
  TGraph* graph_other_ZSFsysDown;
  TGraph* graph_other_ZSFsysUp;
  TGraph* graph_other_JERsysDown;
  TGraph* graph_other_JERsysUp;
  TGraph* graph_other_JESsysDown;
  TGraph* graph_other_JESsysUp;
  TGraph* graph_other_PUsysDown;
  TGraph* graph_other_PUsysUp;
  TGraph* graph_other_SFsysDown;
  TGraph* graph_other_SFsysUp;
  TGraph* graph_other_METsysDown;
  TGraph* graph_other_METsysUp;
  TGraph* graph_other_eleSFsysDown;
  TGraph* graph_other_eleSFsysUp;

  
  
  TGraph* graph_twdr;
  TGraph* graph_twdr_ZSFsysDown; 
  TGraph* graph_twdr_ZSFsysUp; 
  TGraph* graph_twdr_JERsysDown; 
  TGraph* graph_twdr_JERsysUp; 
  TGraph* graph_twdr_JESsysDown; 
  TGraph* graph_twdr_JESsysUp;
  TGraph* graph_twdr_PUsysDown; 
  TGraph* graph_twdr_PUsysUp;
  TGraph* graph_twdr_SFsysDown; 
  TGraph* graph_twdr_SFsysUp;
  TGraph* graph_twdr_METsysDown; 
  TGraph* graph_twdr_METsysUp;
  TGraph* graph_twdr_TopMassDown; 
  TGraph* graph_twdr_TopMassUp;
  TGraph* graph_twdr_Q2Down; 
  TGraph* graph_twdr_Q2Up;
  TGraph* graph_twdr_eleSFsysDown; 
  TGraph* graph_twdr_eleSFsysUp;
  TGraph* graph_twds; 
  
   for(int iPlots = 0; iPlots< nPlots; iPlots++) 
   {
      
     leg = new TLegend(0.7,0.7,0.94,0.94); 
     leg ->SetFillStyle(1001);
     leg ->SetFillColor(kWhite);
     leg ->SetBorderSize(1);
     
     legZSF = new TLegend(0.7,0.7,0.94,0.94);
     legZSF ->SetFillStyle(1001);
     legZSF ->SetFillColor(kWhite);
     legZSF ->SetBorderSize(1);
     
     legJES = new TLegend(0.7,0.7,0.94,0.94);
     legJES ->SetFillStyle(1001);
     legJES ->SetFillColor(kWhite);
     legJES ->SetBorderSize(1);
     
     legPU = new TLegend(0.7,0.7,0.94,0.94);
     legPU ->SetFillStyle(1001);
     legPU ->SetFillColor(kWhite);
     legPU ->SetBorderSize(1);
     
     legSF = new TLegend(0.7,0.7,0.94,0.94);
     legSF ->SetFillStyle(1001);
     legSF ->SetFillColor(kWhite);
     legSF ->SetBorderSize(1);
     
     legMET = new TLegend(0.7,0.7,0.94,0.94);
     legMET ->SetFillStyle(1001);
     legMET ->SetFillColor(kWhite);
     legMET ->SetBorderSize(1);
     
     legTopMass = new TLegend(0.7,0.7,0.94,0.94);
     legTopMass ->SetFillStyle(1001);
     legTopMass ->SetFillColor(kWhite);
     legTopMass ->SetBorderSize(1);
     
     legQ2 = new TLegend(0.7,0.7,0.94,0.94);
     legQ2 ->SetFillStyle(1001);
     legQ2 ->SetFillColor(kWhite);
     legQ2 ->SetBorderSize(1);
     
     legeleSF = new TLegend(0.7,0.7,0.94,0.94);
     legeleSF ->SetFillStyle(1001);
     legeleSF ->SetFillColor(kWhite);
     legeleSF ->SetBorderSize(1);
     
     legmatching = new TLegend(0.7,0.7,0.94,0.94);
     legmatching ->SetFillStyle(1001);
     legmatching ->SetFillColor(kWhite);
     legmatching ->SetBorderSize(1);
 
      legDS = new TLegend(0.7,0.7,0.94,0.94);
     legDS ->SetFillStyle(1001);
     legDS ->SetFillColor(kWhite);
     legDS ->SetBorderSize(1);
     
     
     
     

     /////////////////
	//other
     ///////////////////

       histo_other= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[3]);
	      // cout << "Normal: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 
	histo_other_ZSFsysDown= (TH1F*) _fileZSFsysDown->Get(cutLabel[iPlots]+ "_" + processName[3]);
	     //  cout << "ZSFsysDown: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 	       
       histo_other_ZSFsysUp= (TH1F*) _fileZSFsysUp->Get(cutLabel[iPlots]+ "_" + processName[3]);
	    //   cout << "ZSFsysUp: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 		
       histo_other_JERsysDown= (TH1F*) _fileJERsysDown->Get(cutLabel[iPlots]+ "_" + processName[3]);
	      // cout << "JERsysDown: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 	       
       histo_other_JERsysUp= (TH1F*) _fileJERsysUp->Get(cutLabel[iPlots]+ "_" + processName[3]);
	       //cout << "JERsysUp: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 
       histo_other_JESsysDown= (TH1F*) _fileJESsysDown->Get(cutLabel[iPlots]+ "_" + processName[3]);
	     //  cout << "JESsysDown: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 	       
       histo_other_JESsysUp= (TH1F*) _fileJESsysUp->Get(cutLabel[iPlots]+ "_" + processName[3]);
	    //   cout << "JESsysUp: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 
       histo_other_PUsysDown= (TH1F*) _filePUsysDown->Get(cutLabel[iPlots]+ "_" + processName[3]);
	      // cout << "PUsysDown: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 	       
       histo_other_PUsysUp= (TH1F*) _filePUsysUp->Get(cutLabel[iPlots]+ "_" + processName[3]);
	      // cout << "PUsysUp: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 
       histo_other_SFsysDown= (TH1F*) _fileSFsysDown->Get(cutLabel[iPlots]+ "_" + processName[3]);
	     //  cout << "SFsysDown: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 	       
       histo_other_SFsysUp= (TH1F*) _fileSFsysUp->Get(cutLabel[iPlots]+ "_" + processName[3]);
	      // cout << "SFsysUp: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 
       histo_other_METsysDown= (TH1F*) _fileMETsysDown->Get(cutLabel[iPlots]+ "_" + processName[3]);
	    //   cout << "METsysDown: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 	       
       histo_other_METsysUp= (TH1F*) _fileMETsysUp->Get(cutLabel[iPlots]+ "_" + processName[3]);
	    //   cout << "METsysUp: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 
	histo_other_eleSFsysDown= (TH1F*) _fileeleSFsysDown->Get(cutLabel[iPlots]+ "_" + processName[3]);
	    //   cout << "eleSFsysDown: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 	       
       histo_other_eleSFsysUp= (TH1F*) _fileeleSFsysUp->Get(cutLabel[iPlots]+ "_" + processName[3]);
	    //   cout << "eleSFsysUp: " << cutLabel[iPlots]+ "_" + processName[3] << endl; 
	       
	       
	       	       
       histo_other->Rebin(rebinHisto[iPlots]);
       histo_other->SetLineColor(color[3]);
       histo_other->SetLineWidth(2);
       
       histo_other_ZSFsysDown->Rebin(rebinHisto[iPlots]);
       histo_other_ZSFsysDown->SetLineColor(color[3]);
       histo_other_ZSFsysDown->SetLineStyle(3); 
       histo_other_ZSFsysDown->SetLineWidth(2);
       
       histo_other_ZSFsysUp->Rebin(rebinHisto[iPlots]);
       histo_other_ZSFsysUp->SetLineColor(color[3]);
       histo_other_ZSFsysUp->SetLineStyle(7); 
       histo_other_ZSFsysUp->SetLineWidth(2);
       
       histo_other_JERsysDown->Rebin(rebinHisto[iPlots]);
       histo_other_JERsysDown->SetLineColor(color[3]);
       histo_other_JERsysDown->SetLineStyle(3); 
       histo_other_JERsysDown->SetLineWidth(2);
       
       histo_other_JERsysUp->Rebin(rebinHisto[iPlots]);
       histo_other_JERsysUp->SetLineColor(color[3]);
       histo_other_JERsysUp->SetLineStyle(7); 
       histo_other_JERsysUp->SetLineWidth(2);     
       
       histo_other_JESsysDown->Rebin(rebinHisto[iPlots]);
       histo_other_JESsysDown->SetLineColor(color[3]);
       histo_other_JESsysDown->SetLineStyle(3); 
       histo_other_JESsysDown->SetLineWidth(2);
       
       histo_other_JESsysUp->Rebin(rebinHisto[iPlots]);
       histo_other_JESsysUp->SetLineColor(color[3]);
       histo_other_JESsysUp->SetLineStyle(7); 
       histo_other_JESsysUp->SetLineWidth(2);
       
       histo_other_PUsysDown->Rebin(rebinHisto[iPlots]);
       histo_other_PUsysDown->SetLineColor(color[3]);
       histo_other_PUsysDown->SetLineStyle(3); 
       histo_other_PUsysDown->SetLineWidth(2);
       
       histo_other_PUsysUp->Rebin(rebinHisto[iPlots]);
       histo_other_PUsysUp->SetLineColor(color[3]);
       histo_other_PUsysUp->SetLineStyle(7); 
       histo_other_PUsysUp->SetLineWidth(2); 
       
       histo_other_SFsysDown->Rebin(rebinHisto[iPlots]);
       histo_other_SFsysDown->SetLineColor(color[3]);
       histo_other_SFsysDown->SetLineStyle(3); 
       histo_other_SFsysDown->SetLineWidth(2);
       
       histo_other_SFsysUp->Rebin(rebinHisto[iPlots]);
       histo_other_SFsysUp->SetLineColor(color[3]);
       histo_other_SFsysUp->SetLineStyle(7); 
       histo_other_SFsysUp->SetLineWidth(2);
       
       histo_other_METsysDown->Rebin(rebinHisto[iPlots]);
       histo_other_METsysDown->SetLineColor(color[3]);
       histo_other_METsysDown->SetLineStyle(3); 
       histo_other_METsysDown->SetLineWidth(2);
       
       histo_other_METsysUp->Rebin(rebinHisto[iPlots]);
       histo_other_METsysUp->SetLineColor(color[3]);
       histo_other_METsysUp->SetLineStyle(7); 
       histo_other_METsysUp->SetLineWidth(2);   
  
       
       
       histo_other_eleSFsysDown->Rebin(rebinHisto[iPlots]);
       histo_other_eleSFsysDown->SetLineColor(color[3]);
       histo_other_eleSFsysDown->SetLineStyle(3); 
       histo_other_eleSFsysDown->SetLineWidth(2);
       
       histo_other_eleSFsysUp->Rebin(rebinHisto[iPlots]);
       histo_other_eleSFsysUp->SetLineColor(color[3]);
       histo_other_eleSFsysUp->SetLineStyle(7); 
       histo_other_eleSFsysUp->SetLineWidth(2);
    
       
       const int graph_size_other = histo_other->GetSize();
      // cout << "graph_size_other: " << graph_size_other << endl; 
      
      if(graph_size_other == 27){
       double graph_x_other[27];
       double graph_y_other[27];
        double graph_x_other_Q2Down[27];
       double graph_y_other_Q2Down[27];
        double graph_x_other_Q2Up[27];
       double graph_y_other_Q2Up[27];
       
       
       for(int i_bin = 0; i_bin < (int) graph_size_other ; i_bin++){
             graph_x_other[i_bin] = histo_other->GetBin(i_bin);
	      cout << "graph_x_other: " << i_bin<< " - "<< graph_x_other[i_bin] << endl;
	     graph_y_other[i_bin] = histo_other->GetBinContent(i_bin);
	      cout << "graph_y_other: " << i_bin << " - "<< graph_y_other[i_bin] << endl;
       }
       
       
       
       graph_other = new TGraph(graph_size_other,graph_x_other,graph_y_other); 
     /*  graph_other_ZSFsysDown = new TGraph( 
       graph_other_ZSFsysUp = new TGraph( 
       graph_other_JERsysDown = new TGraph( 
       graph_other_JERsysUp = new TGraph( 
       graph_other_JESsysDown = new TGraph( 
       graph_other_JESsysUp = new TGraph( 
       graph_other_PUsysDown = new TGraph( 
       graph_other_PUsysUp = new TGraph( 
       graph_other_SFsysDown = new TGraph( 
       graph_other_SFsysUp = new TGraph( 
       graph_other_METsysDown = new TGraph( 
       graph_other_METsysUp = new TGraph( 
	graph_other_eleSFsysDown
       graph_other_eleSFsysUp = new TGraph( 
       
       
      */ 
       
	//////////////////////////
       // Get histos tt bar 
       //////////////////////////////////
       histo_tt= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[1]);
	      // cout << "Normal: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
	 histo_tt_ZSFsysDown= (TH1F*) _fileZSFsysDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	     //  cout << "ZSFsysDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_ZSFsysUp= (TH1F*) _fileZSFsysUp->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "ZSFsysUp: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
       histo_tt_JERsysDown= (TH1F*) _fileJERsysDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	      // cout << "JERsysDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_JERsysUp= (TH1F*) _fileJERsysUp->Get(cutLabel[iPlots]+ "_" + processName[1]);
	       //cout << "JERsysUp: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
       histo_tt_JESsysDown= (TH1F*) _fileJESsysDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	     //  cout << "JESsysDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_JESsysUp= (TH1F*) _fileJESsysUp->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "JESsysUp: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
       histo_tt_PUsysDown= (TH1F*) _filePUsysDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	      // cout << "PUsysDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_PUsysUp= (TH1F*) _filePUsysUp->Get(cutLabel[iPlots]+ "_" + processName[1]);
	      // cout << "PUsysUp: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
       histo_tt_SFsysDown= (TH1F*) _fileSFsysDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	     //  cout << "SFsysDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_SFsysUp= (TH1F*) _fileSFsysUp->Get(cutLabel[iPlots]+ "_" + processName[1]);
	      // cout << "SFsysUp: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
       histo_tt_METsysDown= (TH1F*) _fileMETsysDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "METsysDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_METsysUp= (TH1F*) _fileMETsysUp->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "METsysUp: " << cutLabel[iPlots]+ "_" + processName[1] << endl;
       histo_tt_TopMassDown= (TH1F*) _fileTopMassDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "TopMassDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_TopMassUp= (TH1F*) _fileTopMassUp->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "TopMassUp: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
	histo_tt_Q2Down= (TH1F*) _fileQ2Down->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "Q2Down: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_Q2Up= (TH1F*) _fileQ2Up->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "Q2Up: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
	histo_tt_eleSFsysDown= (TH1F*) _fileeleSFsysDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "eleSFsysDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_eleSFsysUp= (TH1F*) _fileeleSFsysUp->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "eleSFsysUp: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
	histo_tt_matchingDown= (TH1F*) _filematchingDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "matchingDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_matchingUp= (TH1F*) _filematchingUp->Get(cutLabel[iPlots]+ "_" + processName[1]);
	    //   cout << "matchingUp: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
		       
	       
	       	       
       histo_tt->Rebin(rebinHisto[iPlots]);
       histo_tt->SetLineColor(color[0]);
       histo_tt->SetLineWidth(2);
       
       histo_tt_ZSFsysDown->Rebin(rebinHisto[iPlots]);
       histo_tt_ZSFsysDown->SetLineColor(color[0]);
       histo_tt_ZSFsysDown->SetLineStyle(3); 
       histo_tt_ZSFsysDown->SetLineWidth(2);
       
       histo_tt_ZSFsysUp->Rebin(rebinHisto[iPlots]);
       histo_tt_ZSFsysUp->SetLineColor(color[0]);
       histo_tt_ZSFsysUp->SetLineStyle(7); 
       histo_tt_ZSFsysUp->SetLineWidth(2);
       
       histo_tt_JERsysDown->Rebin(rebinHisto[iPlots]);
       histo_tt_JERsysDown->SetLineColor(color[0]);
       histo_tt_JERsysDown->SetLineStyle(3); 
       histo_tt_JERsysDown->SetLineWidth(2);
       
       histo_tt_JERsysUp->Rebin(rebinHisto[iPlots]);
       histo_tt_JERsysUp->SetLineColor(color[0]);
       histo_tt_JERsysUp->SetLineStyle(7); 
       histo_tt_JERsysUp->SetLineWidth(2);     
       
       histo_tt_JESsysDown->Rebin(rebinHisto[iPlots]);
       histo_tt_JESsysDown->SetLineColor(color[0]);
       histo_tt_JESsysDown->SetLineStyle(3); 
       histo_tt_JESsysDown->SetLineWidth(2);
       
       histo_tt_JESsysUp->Rebin(rebinHisto[iPlots]);
       histo_tt_JESsysUp->SetLineColor(color[0]);
       histo_tt_JESsysUp->SetLineStyle(7); 
       histo_tt_JESsysUp->SetLineWidth(2);
       
       histo_tt_PUsysDown->Rebin(rebinHisto[iPlots]);
       histo_tt_PUsysDown->SetLineColor(color[0]);
       histo_tt_PUsysDown->SetLineStyle(3); 
       histo_tt_PUsysDown->SetLineWidth(2);
       
       histo_tt_PUsysUp->Rebin(rebinHisto[iPlots]);
       histo_tt_PUsysUp->SetLineColor(color[0]);
       histo_tt_PUsysUp->SetLineStyle(7); 
       histo_tt_PUsysUp->SetLineWidth(2); 
       
       histo_tt_SFsysDown->Rebin(rebinHisto[iPlots]);
       histo_tt_SFsysDown->SetLineColor(color[0]);
       histo_tt_SFsysDown->SetLineStyle(3); 
       histo_tt_SFsysDown->SetLineWidth(2);
       
       histo_tt_SFsysUp->Rebin(rebinHisto[iPlots]);
       histo_tt_SFsysUp->SetLineColor(color[0]);
       histo_tt_SFsysUp->SetLineStyle(7); 
       histo_tt_SFsysUp->SetLineWidth(2);
       
       histo_tt_METsysDown->Rebin(rebinHisto[iPlots]);
       histo_tt_METsysDown->SetLineColor(color[0]);
       histo_tt_METsysDown->SetLineStyle(3); 
       histo_tt_METsysDown->SetLineWidth(2);
       
       histo_tt_METsysUp->Rebin(rebinHisto[iPlots]);
       histo_tt_METsysUp->SetLineColor(color[0]);
       histo_tt_METsysUp->SetLineStyle(7); 
       histo_tt_METsysUp->SetLineWidth(2);   
       
       histo_tt_TopMassDown->Rebin(rebinHisto[iPlots]);
       histo_tt_TopMassDown->SetLineColor(color[0]);
       histo_tt_TopMassDown->SetLineStyle(3); 
       histo_tt_TopMassDown->SetLineWidth(2);
       
       histo_tt_TopMassUp->Rebin(rebinHisto[iPlots]);
       histo_tt_TopMassUp->SetLineColor(color[0]);
       histo_tt_TopMassUp->SetLineStyle(7); 
       histo_tt_TopMassUp->SetLineWidth(2);
       
       histo_tt_Q2Down->Rebin(rebinHisto[iPlots]);
       histo_tt_Q2Down->SetLineColor(color[0]);
       histo_tt_Q2Down->SetLineStyle(3); 
       histo_tt_Q2Down->SetLineWidth(2);
       
       histo_tt_Q2Up->Rebin(rebinHisto[iPlots]);
       histo_tt_Q2Up->SetLineColor(color[0]);
       histo_tt_Q2Up->SetLineStyle(7); 
       histo_tt_Q2Up->SetLineWidth(2);
       
       
       histo_tt_eleSFsysDown->Rebin(rebinHisto[iPlots]);
       histo_tt_eleSFsysDown->SetLineColor(color[0]);
       histo_tt_eleSFsysDown->SetLineStyle(3); 
       histo_tt_eleSFsysDown->SetLineWidth(2);
       
       histo_tt_eleSFsysUp->Rebin(rebinHisto[iPlots]);
       histo_tt_eleSFsysUp->SetLineColor(color[0]);
       histo_tt_eleSFsysUp->SetLineStyle(7); 
       histo_tt_eleSFsysUp->SetLineWidth(2);
       
       histo_tt_matchingDown->Rebin(rebinHisto[iPlots]);
       histo_tt_matchingDown->SetLineColor(color[0]);
       histo_tt_matchingDown->SetLineStyle(3); 
       histo_tt_matchingDown->SetLineWidth(2);
       
       histo_tt_matchingUp->Rebin(rebinHisto[iPlots]);
       histo_tt_matchingUp->SetLineColor(color[0]);
       histo_tt_matchingUp->SetLineStyle(7); 
       histo_tt_matchingUp->SetLineWidth(2);     
       
       
       const int graph_size_tt = histo_tt->GetSize();
      // cout << "graph_size_tt: " << graph_size_tt << endl; 
       double graph_x_tt[27];
       double graph_y_tt[27];
       double graph_x_tt_Q2Up[27];
       double graph_y_tt_Q2Up[27];
       double graph_x_tt_Q2Down[27];
       double graph_y_tt_Q2Down[27];
       
       
       for(int i_bin = 0; i_bin < (int) graph_size_tt ; i_bin++){
             graph_x_tt[i_bin] = histo_tt->GetBin(i_bin);
	      cout << "graph_x_tt: " << i_bin<< " - "<< graph_x_tt[i_bin] << endl;
	     graph_y_tt[i_bin] = histo_tt->GetBinContent(i_bin);
	      cout << "graph_y_tt: " << i_bin << " - "<< graph_y_tt[i_bin] << endl;
       }
       
       
       
       graph_tt = new TGraph(graph_size_tt,graph_x_tt,graph_y_tt); 
       
       
       /*
       
  	//twds
	
	histo_twds= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[2]);
        //cout << "Normal: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
	
	       histo_twds->Rebin(rebinHisto[iPlots]);
              // histo_twdr->SetFillColor(color[0]);
       histo_twds->SetLineColor(color[2]);
       histo_twds->SetLineWidth(2);
	*/
       // get histos twdr
       histo_twdr= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[0]);
        //cout << "Normal: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_ZSFsysDown= (TH1F*) _fileZSFsysDown->Get(cutLabel[iPlots]+ "_" + processName[0]);
        //cout << "ZSFsysDown: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_ZSFsysUp= (TH1F*) _fileZSFsysUp->Get(cutLabel[iPlots]+ "_" + processName[0]);
        //cout << "ZSFsysUp: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_JERsysDown= (TH1F*) _fileJERsysDown->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "JERsysDown: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_JERsysUp= (TH1F*) _fileJERsysUp->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "JERsysUp: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_JESsysDown= (TH1F*) _fileJESsysDown->Get(cutLabel[iPlots]+ "_" + processName[0]);
        //cout << "JESsysDown: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_JESsysUp= (TH1F*) _fileJESsysUp->Get(cutLabel[iPlots]+ "_" + processName[0]);
        //cout << "JESsysUp: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_PUsysDown= (TH1F*) _filePUsysDown->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "PUsysDown: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_PUsysUp= (TH1F*) _filePUsysUp->Get(cutLabel[iPlots]+ "_" + processName[0]);
        //cout << "PUsysUp: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_SFsysDown= (TH1F*) _fileSFsysDown->Get(cutLabel[iPlots]+ "_" + processName[0]);
        //cout << "SFsysDown: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_SFsysUp= (TH1F*) _fileSFsysUp->Get(cutLabel[iPlots]+ "_" + processName[0]);
        //cout << "SFsysUp: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_METsysDown= (TH1F*) _fileMETsysDown->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "METsysDown: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_METsysUp= (TH1F*) _fileMETsysUp->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "METsysUp: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_TopMassDown= (TH1F*) _fileTopMassDown->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "TopMassDown: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_TopMassUp= (TH1F*) _fileTopMassUp->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "TopMassUp: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_Q2Down= (TH1F*) _fileQ2Down->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "Q2Down: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_Q2Up= (TH1F*) _fileQ2Up->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "Q2Up: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_eleSFsysDown= (TH1F*) _fileeleSFsysDown->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "eleSFsysDown: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
       histo_twdr_eleSFsysUp= (TH1F*) _fileeleSFsysUp->Get(cutLabel[iPlots]+ "_" + processName[0]);
       // cout << "eleSFsysUp: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 

	
	
       histo_twdr->Rebin(rebinHisto[iPlots]);
              // histo_twdr->SetFillColor(color[0]);
       histo_twdr->SetLineColor(color[1]);
       histo_twdr->SetLineWidth(2);
       
       histo_twdr_ZSFsysDown->Rebin(rebinHisto[iPlots]);       
       histo_twdr_ZSFsysDown->SetLineColor(color[1]);
       histo_twdr_ZSFsysDown->SetLineStyle(3); 
       histo_twdr_ZSFsysDown->SetLineWidth(2);
       
       histo_twdr_ZSFsysUp->Rebin(rebinHisto[iPlots]);
       histo_twdr_ZSFsysUp->SetLineColor(color[1]);
       histo_twdr_ZSFsysUp->SetLineStyle(7); 
       histo_twdr_ZSFsysUp->SetLineWidth(2);
       
       histo_twdr_JERsysDown->Rebin(rebinHisto[iPlots]);       
       histo_twdr_JERsysDown->SetLineColor(color[1]);
       histo_twdr_JERsysDown->SetLineStyle(3); 
       histo_twdr_JERsysDown->SetLineWidth(2);
       
       histo_twdr_JERsysUp->Rebin(rebinHisto[iPlots]);
       histo_twdr_JERsysUp->SetLineColor(color[1]);
       histo_twdr_JERsysUp->SetLineStyle(7); 
       histo_twdr_JERsysUp->SetLineWidth(2);

       histo_twdr_JESsysDown->Rebin(rebinHisto[iPlots]);       
       histo_twdr_JESsysDown->SetLineColor(color[1]);
       histo_twdr_JESsysDown->SetLineStyle(3); 
       histo_twdr_JESsysDown->SetLineWidth(2);
       
       histo_twdr_JESsysUp->Rebin(rebinHisto[iPlots]);
       histo_twdr_JESsysUp->SetLineColor(color[1]);
       histo_twdr_JESsysUp->SetLineStyle(7); 
       histo_twdr_JESsysUp->SetLineWidth(2);
       
       histo_twdr_PUsysDown->Rebin(rebinHisto[iPlots]);       
       histo_twdr_PUsysDown->SetLineColor(color[1]);
       histo_twdr_PUsysDown->SetLineStyle(3); 
       histo_twdr_PUsysDown->SetLineWidth(2);
       
       histo_twdr_PUsysUp->Rebin(rebinHisto[iPlots]);
       histo_twdr_PUsysUp->SetLineColor(color[1]);
       histo_twdr_PUsysUp->SetLineStyle(7); 
       histo_twdr_PUsysUp->SetLineWidth(2);
       
       histo_twdr_SFsysDown->Rebin(rebinHisto[iPlots]);       
       histo_twdr_SFsysDown->SetLineColor(color[1]);
       histo_twdr_SFsysDown->SetLineStyle(3); 
       histo_twdr_SFsysDown->SetLineWidth(2);
       
       histo_twdr_SFsysUp->Rebin(rebinHisto[iPlots]);
       histo_twdr_SFsysUp->SetLineColor(color[1]);
       histo_twdr_SFsysUp->SetLineStyle(7); 
       histo_twdr_SFsysUp->SetLineWidth(2);
       
       histo_twdr_METsysDown->Rebin(rebinHisto[iPlots]);       
       histo_twdr_METsysDown->SetLineColor(color[1]);
       histo_twdr_METsysDown->SetLineStyle(3); 
       histo_twdr_METsysDown->SetLineWidth(2);
       
       histo_twdr_METsysUp->Rebin(rebinHisto[iPlots]);
       histo_twdr_METsysUp->SetLineColor(color[1]);
       histo_twdr_METsysUp->SetLineStyle(7); 
       histo_twdr_METsysUp->SetLineWidth(2);   
       
       histo_twdr_TopMassDown->Rebin(rebinHisto[iPlots]);       
       histo_twdr_TopMassDown->SetLineColor(color[1]);
       histo_twdr_TopMassDown->SetLineStyle(3); 
       histo_twdr_TopMassDown->SetLineWidth(2);
       
       histo_twdr_TopMassUp->Rebin(rebinHisto[iPlots]);
       histo_twdr_TopMassUp->SetLineColor(color[1]);
       histo_twdr_TopMassUp->SetLineStyle(7); 
       histo_twdr_TopMassUp->SetLineWidth(2); 
       
       histo_twdr_Q2Down->Rebin(rebinHisto[iPlots]);       
       histo_twdr_Q2Down->SetLineColor(color[1]);
       histo_twdr_Q2Down->SetLineStyle(3); 
       histo_twdr_Q2Down->SetLineWidth(2);
       
       histo_twdr_Q2Up->Rebin(rebinHisto[iPlots]);
       histo_twdr_Q2Up->SetLineColor(color[1]);
       histo_twdr_Q2Up->SetLineStyle(7); 
       histo_twdr_Q2Up->SetLineWidth(2); 
       
       histo_twdr_eleSFsysDown->Rebin(rebinHisto[iPlots]);       
       histo_twdr_eleSFsysDown->SetLineColor(color[1]);
       histo_twdr_eleSFsysDown->SetLineStyle(3); 
       histo_twdr_eleSFsysDown->SetLineWidth(2);
       
       histo_twdr_eleSFsysUp->Rebin(rebinHisto[iPlots]);
       histo_twdr_eleSFsysUp->SetLineColor(color[1]);
       histo_twdr_eleSFsysUp->SetLineStyle(7); 
       histo_twdr_eleSFsysUp->SetLineWidth(2); 
       
       const int graph_size_twdr = histo_twdr->GetSize();
      // cout << "graph_size_twdr: " << graph_size_twdr << endl; 
       
       double graph_x_twdr[27];
       double graph_y_twdr[27];
       double graph_x_twdr_Q2Down[27];
       double graph_y_twdr_Q2Down[27];
       double graph_x_twdr_Q2Up[27];
       double graph_y_twdr_Q2Up[27];
       
       for(int i_bin = 0; i_bin < (int) graph_size_twdr ; i_bin++){
             graph_x_twdr[i_bin] = histo_twdr->GetBin(i_bin);
	      cout << "graph_x_twdr: " << i_bin<< " - "<< graph_x_twdr[i_bin] << endl;
	     graph_y_twdr[i_bin] = histo_twdr->GetBinContent(i_bin);
	      cout << "graph_y_twdr: " << i_bin << " - "<< graph_y_twdr[i_bin] << endl;
	  
       }
       
             
       graph_twdr = new TGraph(graph_size_twdr,graph_x_twdr,graph_y_twdr); 
       
        
	
	
	// Set color
	// ttbar
       graph_tt->SetFillColor(kRed);
       
       // twdr
       graph_twdr->SetFillColor(kBlue);
       
       
      // other
       graph_other->SetFillColor(kMagenta);
       
       
       //// Set max
       double max1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max11 = TMath::Max(histo_other->GetMaximum(), max1);
       
       
       
       
	// Make legends      
       leg ->AddEntry(histo_tt, "tt bckgr." , "l");  
       leg ->AddEntry(histo_twdr, "twdr signal", "l"); 
       leg ->AddEntry(histo_other, "other bckgr.", "l"); 
       
	
	
	
       
       
	
	TMultiGraph  *mg1  = new TMultiGraph();
	mg1->Add(graph_tt);
	mg1->Add(graph_twdr);
	mg1->Add(graph_other);
	mg1->SetMaximum(max11 * 1.5);
	
	
	
	
      TCanvas *c1 = new TCanvas();
       mg1->Draw("AB");
       mg1 -> GetXaxis()->SetTitle(cutTitle[iPlots]);
       mg1->GetYaxis()->SetTitle("#evts");
       
      for(int i_bin = 0; i_bin < (int) graph_size_twdr ; i_bin++){
             graph_x_twdr_Q2Down[i_bin] = histo_twdr_Q2Down->GetBin(i_bin);
	    //  cout << "graph_x_twdr_Q2Down: " << i_bin<< " - "<< graph_x_twdr_Q2Down[i_bin] << endl;
	     graph_y_twdr_Q2Down[i_bin] = histo_twdr_Q2Down->GetBinContent(i_bin);
	    //  cout << "graph_y_twdr_Q2Down: " << i_bin << " - "<< graph_y_twdr_Q2Down[i_bin] << endl;
	     TMarker *m_other_Q2Down = new TMarker(graph_x_twdr_Q2Down[i_bin], graph_y_twdr_Q2Down[i_bin],22);
	     m_other_Q2Down-> SetMarkerSize(1); 
	     m_other_Q2Down-> SetMarkerColor(38); 
	     m_other_Q2Down->Draw(); 
	     
	     graph_x_twdr_Q2Up[i_bin] = histo_twdr_Q2Up->GetBin(i_bin);
	    //  cout << "graph_x_twdr_Q2Up: " << i_bin<< " - "<< graph_x_twdr_Q2Up[i_bin] << endl;
	     graph_y_twdr_Q2Up[i_bin] = histo_twdr_Q2Up->GetBinContent(i_bin);
	    //  cout << "graph_y_twdr_Q2Up: " << i_bin << " - "<< graph_y_twdr_Q2Up[i_bin] << endl;
	     TMarker *m_other_Q2Up = new TMarker(graph_x_twdr_Q2Up[i_bin], graph_y_twdr_Q2Up[i_bin],22);
	     m_other_Q2Up-> SetMarkerSize(1); 
	     m_other_Q2Up-> SetMarkerColor(38); 
	     m_other_Q2Up->Draw(); 
	     
	     
	     graph_x_tt_Q2Down[i_bin] = histo_tt_Q2Down->GetBin(i_bin);
	     // cout << "graph_x_tt_Q2Down: " << i_bin<< " - "<< graph_x_tt_Q2Down[i_bin] << endl;
	     graph_y_tt_Q2Down[i_bin] = histo_tt_Q2Down->GetBinContent(i_bin);
	     // cout << "graph_y_tt_Q2Down: " << i_bin << " - "<< graph_y_tt_Q2Down[i_bin] << endl;
	     TMarker *m_tt_Q2Down = new TMarker(graph_x_tt_Q2Down[i_bin], graph_y_tt_Q2Down[i_bin],22);
	     m_tt_Q2Down-> SetMarkerSize(1); 
	     m_tt_Q2Down-> SetMarkerColor(46); 
	     m_tt_Q2Down->Draw();
	     
	     graph_x_tt_Q2Up[i_bin] = histo_tt_Q2Up->GetBin(i_bin);
	     // cout << "graph_x_tt_Q2Up: " << i_bin<< " - "<< graph_x_tt_Q2Up[i_bin] << endl;
	     graph_y_tt_Q2Up[i_bin] = histo_tt_Q2Up->GetBinContent(i_bin);
	     // cout << "graph_y_tt_Q2Up: " << i_bin << " - "<< graph_y_tt_Q2Up[i_bin] << endl;
	     TMarker *m_tt_Q2Up = new TMarker(graph_x_tt_Q2Up[i_bin], graph_y_tt_Q2Up[i_bin],22);
	     m_tt_Q2Up-> SetMarkerSize(1); 
	     m_tt_Q2Up-> SetMarkerColor(46); 
	     m_tt_Q2Up->Draw();
	     
	     
	      graph_x_twdr_Q2Down[i_bin] = histo_twdr_Q2Down->GetBin(i_bin);
	     // cout << "graph_x_twdr_Q2Down: " << i_bin<< " - "<< graph_x_twdr_Q2Down[i_bin] << endl;
	     graph_y_twdr_Q2Down[i_bin] = histo_twdr_Q2Down->GetBinContent(i_bin);
	     // cout << "graph_y_twdr_Q2Down: " << i_bin << " - "<< graph_y_twdr_Q2Down[i_bin] << endl;
	     TMarker *m_twdr_Q2Down = new TMarker(graph_x_twdr_Q2Down[i_bin], graph_y_twdr_Q2Down[i_bin],22);
	     m_twdr_Q2Down-> SetMarkerSize(1); 
	     m_twdr_Q2Down-> SetMarkerColor(6); 
	     m_twdr_Q2Down->Draw();
	     
	     graph_x_twdr_Q2Up[i_bin] = histo_twdr_Q2Up->GetBin(i_bin);
	     // cout << "graph_x_twdr_Q2Up: " << i_bin<< " - "<< graph_x_twdr_Q2Up[i_bin] << endl;
	     graph_y_twdr_Q2Up[i_bin] = histo_twdr_Q2Up->GetBinContent(i_bin);
	     // cout << "graph_y_twdr_Q2Up: " << i_bin << " - "<< graph_y_twdr_Q2Up[i_bin] << endl;
	     TMarker *m_twdr_Q2Up = new TMarker(graph_x_twdr_Q2Up[i_bin], graph_y_twdr_Q2Up[i_bin],22);
	     m_twdr_Q2Up-> SetMarkerSize(1); 
	     m_twdr_Q2Up-> SetMarkerColor(6); 
	     m_twdr_Q2Up->Draw();
	     
       }
       
      
       
       
       
       
      leg->Draw();
      labelcms->Draw();
     labelcms2->Draw();
    
      c1->SaveAs("graph/" + plotAnalysis +"/" + plotExtension + modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + ".png");
    
    
      
     
    }
    
    

    } // end plots loop 
    
 
} // end constructor loop
