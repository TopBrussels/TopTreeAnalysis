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

void syst_plotmaker(int mode = 0, int region = 0){
 
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
     
     
     
     
          aleg = new TLegend(0.7,0.7,0.94,0.94); 
     aleg ->SetFillStyle(1001);
     aleg ->SetFillColor(kWhite);
     aleg ->SetBorderSize(1);
     
     alegZSF = new TLegend(0.7,0.7,0.94,0.94);
     alegZSF ->SetFillStyle(1001);
     alegZSF ->SetFillColor(kWhite);
     alegZSF ->SetBorderSize(1);
     
     alegJES = new TLegend(0.7,0.7,0.94,0.94);
     alegJES ->SetFillStyle(1001);
     alegJES ->SetFillColor(kWhite);
     alegJES ->SetBorderSize(1);
     
     alegPU = new TLegend(0.7,0.7,0.94,0.94);
     alegPU ->SetFillStyle(1001);
     alegPU ->SetFillColor(kWhite);
     alegPU ->SetBorderSize(1);
     
     alegSF = new TLegend(0.7,0.7,0.94,0.94);
     alegSF ->SetFillStyle(1001);
     alegSF ->SetFillColor(kWhite);
     alegSF ->SetBorderSize(1);
     
     alegMET = new TLegend(0.7,0.7,0.94,0.94);
     alegMET ->SetFillStyle(1001);
     alegMET ->SetFillColor(kWhite);
     alegMET ->SetBorderSize(1);
     
     alegTopMass = new TLegend(0.7,0.7,0.94,0.94);
     alegTopMass ->SetFillStyle(1001);
     alegTopMass ->SetFillColor(kWhite);
     alegTopMass ->SetBorderSize(1);
     
     alegQ2 = new TLegend(0.7,0.7,0.94,0.94);
     alegQ2 ->SetFillStyle(1001);
     alegQ2 ->SetFillColor(kWhite);
     alegQ2 ->SetBorderSize(1);
     
     alegeleSF = new TLegend(0.7,0.7,0.94,0.94);
     alegeleSF ->SetFillStyle(1001);
     alegeleSF ->SetFillColor(kWhite);
     alegeleSF ->SetBorderSize(1);
     
     alegmatching = new TLegend(0.7,0.7,0.94,0.94);
     alegmatching ->SetFillStyle(1001);
     alegmatching ->SetFillColor(kWhite);
     alegmatching ->SetBorderSize(1);
 
      alegDS = new TLegend(0.7,0.7,0.94,0.94);
     alegDS ->SetFillStyle(1001);
     alegDS ->SetFillColor(kWhite);
     alegDS ->SetBorderSize(1);
     
     
	//other
	

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
    
       
	
       // Get histos tt bar 
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
       
  	//twds
	
	histo_twds= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[2]);
        //cout << "Normal: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
	
	       histo_twds->Rebin(rebinHisto[iPlots]);
              // histo_twdr->SetFillColor(color[0]);
       histo_twds->SetLineColor(color[2]);
       histo_twds->SetLineWidth(2);
	
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
       
       
   
        
       
       
	// Make legends      
       leg ->AddEntry(histo_tt, "tt bckgr." , "l");  
       leg ->AddEntry(histo_tt_JERsysDown, "tt bckgr. JERsysDown" , "l"); 
       leg ->AddEntry(histo_tt_JERsysUp, "tt bckgr. JERsysUp" , "l"); 
       
       leg ->AddEntry(histo_twdr, "twdr signal", "l"); 
       leg ->AddEntry(histo_twdr_JERsysDown, "twdr signal JERsysDown", "l"); 
       leg ->AddEntry(histo_twdr_JERsysUp, "twdr signal JERsysUp", "l");
       
       leg ->AddEntry(histo_other, "other bckgr.", "l"); 
       leg ->AddEntry(histo_other_JERsysDown, "other bckgr. JERsysDown", "l"); 
       leg ->AddEntry(histo_other_JERsysUp, "other bckgr. JERsysUp", "l");
       
       legZSF ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legZSF ->AddEntry(histo_tt_ZSFsysDown, "tt bckgr. ZSFsysDown" , "l"); 
       legZSF ->AddEntry(histo_tt_ZSFsysUp, "tt bckgr. ZSFsysUp" , "l");   

       legZSF ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legZSF ->AddEntry(histo_twdr_ZSFsysDown, "twdr signal ZSFsysDown", "l"); 
       legZSF ->AddEntry(histo_twdr_ZSFsysUp, "twdr signal ZSFsysUp", "l");
       
        legZSF ->AddEntry(histo_other, "other bckgr.", "l"); 
       legZSF ->AddEntry(histo_other_ZSFsysDown, "other bckgr. ZSFsysDown", "l"); 
       legZSF ->AddEntry(histo_other_ZSFsysUp, "other bckgr. ZSFsysUp", "l");
       
       legJES ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legJES ->AddEntry(histo_tt_JESsysDown, "tt bckgr. JESsysDown" , "l"); 
       legJES ->AddEntry(histo_tt_JESsysUp, "tt bckgr. JESsysUp" , "l");   

       legJES ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legJES ->AddEntry(histo_twdr_JESsysDown, "twdr signal JESsysDown", "l"); 
       legJES ->AddEntry(histo_twdr_JESsysUp, "twdr signal JESsysUp", "l");
       
        legJES ->AddEntry(histo_other, "other bckgr.", "l"); 
       legJES ->AddEntry(histo_other_JESsysDown, "other bckgr. JESsysDown", "l"); 
       legJES ->AddEntry(histo_other_JESsysUp, "other bckgr. JESsysUp", "l");
       
       legPU ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legPU ->AddEntry(histo_tt_PUsysDown, "tt bckgr. PUsysDown" , "l"); 
       legPU ->AddEntry(histo_tt_PUsysUp, "tt bckgr. PUsysUp" , "l");   

       legPU ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legPU ->AddEntry(histo_twdr_PUsysDown, "twdr signal PUsysDown", "l"); 
       legPU ->AddEntry(histo_twdr_PUsysUp, "twdr signal PUsysUp", "l");   
       
       legPU ->AddEntry(histo_other, "other bckgr.", "l"); 
       legPU ->AddEntry(histo_other_PUsysDown, "other bckgr. PUsysDown", "l"); 
       legPU ->AddEntry(histo_other_PUsysUp, "other bckgr. PUsysUp", "l");       

       legSF ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legSF ->AddEntry(histo_tt_SFsysDown, "tt bckgr. SFsysDown" , "l"); 
       legSF ->AddEntry(histo_tt_SFsysUp, "tt bckgr. SFsysUp" , "l");   

       legSF ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legSF ->AddEntry(histo_twdr_SFsysDown, "twdr signal SFsysDown", "l"); 
       legSF ->AddEntry(histo_twdr_SFsysUp, "twdr signal SFsysUp", "l");
       
       legSF ->AddEntry(histo_other, "other bckgr.", "l"); 
       legSF ->AddEntry(histo_other_SFsysDown, "other bckgr. SFsysDown", "l"); 
       legSF ->AddEntry(histo_other_SFsysUp, "other bckgr. SFsysUp", "l"); 
       
       legMET ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legMET ->AddEntry(histo_tt_METsysDown, "tt bckgr. METsysDown" , "l"); 
       legMET ->AddEntry(histo_tt_METsysUp, "tt bckgr. METsysUp" , "l");   

       legMET ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legMET ->AddEntry(histo_tt_METsysDown, "tt bckgr. METsysDown" , "l"); 
       legMET ->AddEntry(histo_tt_METsysUp, "tt bckgr. METsysUp" , "l");   
       
       legMET ->AddEntry(histo_other, "other bckgr.", "l"); 
       legMET ->AddEntry(histo_other_METsysDown, "other bckgr. METsysDown", "l"); 
       legMET ->AddEntry(histo_other_METsysUp, "other bckgr. METsysUp", "l"); 
       
       legTopMass ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legTopMass ->AddEntry(histo_twdr_TopMassDown, "twdr signal Top mass 166.5 GeV", "l"); 
       legTopMass ->AddEntry(histo_twdr_TopMassUp, "twdr signal Top mass 178.5 GeV", "l"); 
       
       legTopMass ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legTopMass ->AddEntry(histo_tt_TopMassDown, "tt bckgr. Top mass 166.5 GeV" , "l"); 
       legTopMass ->AddEntry(histo_tt_TopMassUp, "tt bckgr. Top mass 178.5 GeV" , "l");   
       
       legQ2 ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legQ2 ->AddEntry(histo_twdr_Q2Down, "twdr signal Q2 halved", "l"); 
       legQ2 ->AddEntry(histo_twdr_Q2Up, "twdr signal Q2 doubled", "l");

       legQ2 ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legQ2 ->AddEntry(histo_tt_Q2Down, "tt bckgr. Q2 halved" , "l"); 
       legQ2 ->AddEntry(histo_tt_Q2Up, "tt bckgr. Q2 doubled" , "l");   
              
       legeleSF ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legeleSF ->AddEntry(histo_twdr_eleSFsysDown, "twdr signal eleSF Down", "l"); 
       legeleSF ->AddEntry(histo_twdr_eleSFsysUp, "twdr signal eleSF Up", "l");

       legeleSF ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legeleSF ->AddEntry(histo_tt_eleSFsysDown, "tt bckgr. eleSF Down" , "l"); 
       legeleSF ->AddEntry(histo_tt_eleSFsysUp, "tt bckgr. eleSF Up" , "l");   
              
       legeleSF ->AddEntry(histo_other, "other bckgr.", "l"); 
       legeleSF ->AddEntry(histo_other_eleSFsysDown, "other bckgr. eleSFsysDown", "l"); 
       legeleSF ->AddEntry(histo_other_eleSFsysUp, "other bckgr. eleSFsysUp", "l"); 
       
        legmatching ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legmatching ->AddEntry(histo_tt_matchingDown, "tt bckgr. matching Down" , "l"); 
       legmatching ->AddEntry(histo_tt_matchingUp, "tt bckgr. matching Up" , "l");         
       
        legDS ->AddEntry(histo_twdr, "tw DR signal", "l"); 
	legDS ->AddEntry(histo_twds, "tw DS signal", "l");       
	
	
	
	aleg ->AddEntry(histo_tt, "tt bckgr." , "l");  
       aleg ->AddEntry(histo_tt_JERsysDown, "tt bckgr. JERsysDown" , "l"); 
       aleg ->AddEntry(histo_tt_JERsysUp, "tt bckgr. JERsysUp" , "l"); 
       
       aleg ->AddEntry(histo_twdr, "twdr signal", "l"); 
       aleg ->AddEntry(histo_twdr_JERsysDown, "twdr signal JERsysDown", "l"); 
       aleg ->AddEntry(histo_twdr_JERsysUp, "twdr signal JERsysUp", "l");
       
        alegZSF ->AddEntry(histo_other, "other bckgr." , "l");  
       alegZSF ->AddEntry(histo_other_ZSFsysDown, "other bckgr. ZSFsysDown" , "l"); 
       alegZSF ->AddEntry(histo_other_ZSFsysUp, "other bckgr. ZSFsysUp" , "l");   

      
       
       alegJES ->AddEntry(histo_tt, "tt bckgr." , "l");  
       alegJES ->AddEntry(histo_tt_JESsysDown, "tt bckgr. JESsysDown" , "l"); 
       alegJES ->AddEntry(histo_tt_JESsysUp, "tt bckgr. JESsysUp" , "l");   

       alegJES ->AddEntry(histo_twdr, "twdr signal", "l"); 
       alegJES ->AddEntry(histo_twdr_JESsysDown, "twdr signal JESsysDown", "l"); 
       alegJES ->AddEntry(histo_twdr_JESsysUp, "twdr signal JESsysUp", "l");
       
    
       
       alegPU ->AddEntry(histo_tt, "tt bckgr." , "l");  
       alegPU ->AddEntry(histo_tt_PUsysDown, "tt bckgr. PUsysDown" , "l"); 
       alegPU ->AddEntry(histo_tt_PUsysUp, "tt bckgr. PUsysUp" , "l");   

       alegPU ->AddEntry(histo_twdr, "twdr signal", "l"); 
       alegPU ->AddEntry(histo_twdr_PUsysDown, "twdr signal PUsysDown", "l"); 
       alegPU ->AddEntry(histo_twdr_PUsysUp, "twdr signal PUsysUp", "l");   
       
          

       alegSF ->AddEntry(histo_tt, "tt bckgr." , "l");  
       alegSF ->AddEntry(histo_tt_SFsysDown, "tt bckgr. SFsysDown" , "l"); 
       alegSF ->AddEntry(histo_tt_SFsysUp, "tt bckgr. SFsysUp" , "l");   

       alegSF ->AddEntry(histo_twdr, "twdr signal", "l"); 
       alegSF ->AddEntry(histo_twdr_SFsysDown, "twdr signal SFsysDown", "l"); 
       alegSF ->AddEntry(histo_twdr_SFsysUp, "twdr signal SFsysUp", "l");
       
       
       
       alegMET ->AddEntry(histo_tt, "tt bckgr." , "l");  
       alegMET ->AddEntry(histo_tt_METsysDown, "tt bckgr. METsysDown" , "l"); 
       alegMET ->AddEntry(histo_tt_METsysUp, "tt bckgr. METsysUp" , "l");   

       alegMET ->AddEntry(histo_tt, "tt bckgr." , "l");  
       alegMET ->AddEntry(histo_tt_METsysDown, "tt bckgr. METsysDown" , "l"); 
       alegMET ->AddEntry(histo_tt_METsysUp, "tt bckgr. METsysUp" , "l");   
       
       
       alegTopMass ->AddEntry(histo_twdr, "twdr signal", "l"); 
       alegTopMass ->AddEntry(histo_twdr_TopMassDown, "twdr signal Top mass 166.5 GeV", "l"); 
       alegTopMass ->AddEntry(histo_twdr_TopMassUp, "twdr signal Top mass 178.5 GeV", "l"); 
       
       alegTopMass ->AddEntry(histo_tt, "tt bckgr." , "l");  
       alegTopMass ->AddEntry(histo_tt_TopMassDown, "tt bckgr. Top mass 166.5 GeV" , "l"); 
       alegTopMass ->AddEntry(histo_tt_TopMassUp, "tt bckgr. Top mass 178.5 GeV" , "l");   
       
       alegQ2 ->AddEntry(histo_twdr, "twdr signal", "l"); 
       alegQ2 ->AddEntry(histo_twdr_Q2Down, "twdr signal Q2 halved", "l"); 
       alegQ2 ->AddEntry(histo_twdr_Q2Up, "twdr signal Q2 doubled", "l");

       alegQ2 ->AddEntry(histo_tt, "tt bckgr." , "l");  
       alegQ2 ->AddEntry(histo_tt_Q2Down, "tt bckgr. Q2 halved" , "l"); 
       alegQ2 ->AddEntry(histo_tt_Q2Up, "tt bckgr. Q2 doubled" , "l");   
              
       alegeleSF ->AddEntry(histo_twdr, "twdr signal", "l"); 
       alegeleSF ->AddEntry(histo_twdr_eleSFsysDown, "twdr signal eleSF Down", "l"); 
       alegeleSF ->AddEntry(histo_twdr_eleSFsysUp, "twdr signal eleSF Up", "l");

       alegeleSF ->AddEntry(histo_tt, "tt bckgr." , "l");  
       alegeleSF ->AddEntry(histo_tt_eleSFsysDown, "tt bckgr. eleSF Down" , "l"); 
       alegeleSF ->AddEntry(histo_tt_eleSFsysUp, "tt bckgr. eleSF Up" , "l");   
              

       
        alegmatching ->AddEntry(histo_tt, "tt bckgr." , "l");  
       alegmatching ->AddEntry(histo_tt_matchingDown, "tt bckgr. matching Down" , "l"); 
       alegmatching ->AddEntry(histo_tt_matchingUp, "tt bckgr. matching Up" , "l");         
       
        alegDS ->AddEntry(histo_twdr, "tw DR signal", "l"); 
	alegDS ->AddEntry(histo_twds, "tw DS signal", "l");       
       
       // determine maxima
       double max1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max11 = TMath::Max(histo_other->GetMaximum(), max1);
       double max2 = TMath::Max( histo_tt_JERsysUp->GetMaximum(), histo_twdr_JERsysUp->GetMaximum());
       double max22 = TMath::Max( histo_other_JERsysUp->GetMaximum(), max2);
       double max3 = TMath::Max(histo_tt_JERsysDown->GetMaximum(), histo_twdr_JERsysDown->GetMaximum());
       double max33 = TMath::Max(histo_other_JERsysDown->GetMaximum(), max3);
       double max4 = TMath::Max(max11,max22);
       double max = TMath::Max(max33,max4);
       
              double max_ZSF_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max_ZSF_11 = TMath::Max(histo_other->GetMaximum(), max_ZSF_1);
       double max_ZSF_2 = TMath::Max( histo_tt_ZSFsysUp->GetMaximum(), histo_twdr_ZSFsysUp->GetMaximum());
       double max_ZSF_22 = TMath::Max( histo_other_ZSFsysUp->GetMaximum(), max_ZSF_2);
       double max_ZSF_3 = TMath::Max(histo_tt_ZSFsysDown->GetMaximum(), histo_twdr_ZSFsysDown->GetMaximum());
       double max_ZSF_33 = TMath::Max(histo_other_ZSFsysDown->GetMaximum(), max_ZSF_3);
       double max_ZSF_4 = TMath::Max(max_ZSF_11,max_ZSF_22);
       double max_ZSF = TMath::Max(max_ZSF_33,max_ZSF_4);
       
       double max_JES_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max_JES_11 = TMath::Max(histo_other->GetMaximum(), max_JES_1);
       double max_JES_2 = TMath::Max( histo_tt_JESsysUp->GetMaximum(), histo_twdr_JESsysUp->GetMaximum());
       double max_JES_22 = TMath::Max( histo_other_JESsysUp->GetMaximum(), max_JES_2);
       double max_JES_3 = TMath::Max(histo_tt_JESsysDown->GetMaximum(), histo_twdr_JESsysDown->GetMaximum());
       double max_JES_33 = TMath::Max(histo_other_JESsysDown->GetMaximum(), max_JES_3);
       double max_JES_4 = TMath::Max(max_JES_11,max_JES_22);
       double max_JES = TMath::Max(max_JES_33,max_JES_4);
       
       double max_PU_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max_PU_11 = TMath::Max(histo_other->GetMaximum(), max_PU_1);
       double max_PU_2 = TMath::Max( histo_tt_PUsysUp->GetMaximum(), histo_twdr_PUsysUp->GetMaximum());
       double max_PU_22 = TMath::Max( histo_other_PUsysUp->GetMaximum(), max_PU_2);
       double max_PU_3 = TMath::Max(histo_tt_PUsysDown->GetMaximum(), histo_twdr_PUsysDown->GetMaximum());
       double max_PU_33 = TMath::Max(histo_other_PUsysDown->GetMaximum(), max_PU_3);
       double max_PU_4 = TMath::Max(max_PU_11,max_PU_22);
       double max_PU = TMath::Max(max_PU_33,max_PU_4);

       double max_SF_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max_SF_11 = TMath::Max(histo_other->GetMaximum(), max_SF_1);
       double max_SF_2 = TMath::Max( histo_tt_SFsysUp->GetMaximum(), histo_twdr_SFsysUp->GetMaximum());
       double max_SF_22 = TMath::Max( histo_other_SFsysUp->GetMaximum(), max_SF_2);
       double max_SF_3 = TMath::Max(histo_tt_SFsysDown->GetMaximum(), histo_twdr_SFsysDown->GetMaximum());
       double max_SF_33 = TMath::Max(histo_other_SFsysDown->GetMaximum(), max_SF_3);
       double max_SF_4 = TMath::Max(max_SF_11,max_SF_22);
       double max_SF = TMath::Max(max_SF_33,max_SF_4);
       
       double max_MET_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max_MET_11 = TMath::Max(histo_other->GetMaximum(), max_MET_1);
       double max_MET_2 = TMath::Max( histo_tt_METsysUp->GetMaximum(), histo_twdr_METsysUp->GetMaximum());
       double max_MET_22 = TMath::Max( histo_other_METsysUp->GetMaximum(), max_MET_2);
       double max_MET_3 = TMath::Max(histo_tt_METsysDown->GetMaximum(), histo_twdr_METsysDown->GetMaximum());
       double max_MET_33 = TMath::Max(histo_other_METsysDown->GetMaximum(), max_MET_3);
       double max_MET_4 = TMath::Max(max_MET_11,max_MET_22);
       double max_MET = TMath::Max(max_MET_33,max_MET_4)  ;  
    
       double max_TopMass_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       
       double max_TopMass_2 = TMath::Max( histo_tt_TopMassUp->GetMaximum(), histo_twdr_TopMassUp->GetMaximum());
       
       double max_TopMass_3 = TMath::Max(histo_tt_TopMassDown->GetMaximum(), histo_twdr_JERsysDown->GetMaximum());
       
       double max_TopMass_4 = TMath::Max(max_TopMass_1,max_TopMass_2);
       double max_TopMass = TMath::Max(max_TopMass_3,max_TopMass_4);
       
       double max_Q2_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
      
       double max_Q2_2 = TMath::Max( histo_tt_Q2Up->GetMaximum(), histo_twdr_Q2Up->GetMaximum());
       
       double max_Q2_3 = TMath::Max(histo_tt_Q2Down->GetMaximum(), histo_twdr_Q2Down->GetMaximum());
      
       double max_Q2_4 = TMath::Max(max_Q2_1,max_Q2_2);
       double max_Q2 = TMath::Max(max_Q2_3,max_Q2_4);
       
       double max_eleSF_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max_eleSF_11 = TMath::Max(histo_other->GetMaximum(), max_eleSF_1);
       double max_eleSF_2 = TMath::Max( histo_tt_eleSFsysUp->GetMaximum(), histo_twdr_eleSFsysUp->GetMaximum());
       double max_eleSF_22 = TMath::Max( histo_other_eleSFsysUp->GetMaximum(), max_eleSF_2);
       double max_eleSF_3 = TMath::Max(histo_tt_eleSFsysDown->GetMaximum(), histo_twdr_eleSFsysDown->GetMaximum());
       double max_eleSF_33 = TMath::Max(histo_other_eleSFsysDown->GetMaximum(), max_eleSF_3);
       double max_eleSF_4 = TMath::Max(max_eleSF_11,max_eleSF_22);
       double max_eleSF = TMath::Max(max_eleSF_33,max_eleSF_4);
       
       double max_matching_1 = TMath::Max(histo_tt->GetMaximum(), 0.);
       double max_matching_2 = TMath::Max( histo_tt_matchingUp->GetMaximum(), 0.);
       double max_matching_3 = TMath::Max(histo_tt_matchingDown->GetMaximum(),0.);
       double max_matching_4 = TMath::Max(max_matching_1,max_matching_2);
       double max_matching = TMath::Max(max_matching_3,max_matching_4);
       
        double max_DS = TMath::Max(histo_twds->GetMaximum(), histo_twdr->GetMaximum());
       
      
       TCanvas *c1 = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->Draw("h");
       
       histo_tt_JERsysUp->Draw("h, sames");
       histo_tt_JERsysUp->SetMaximum(max * 1.5);
       histo_tt_JERsysUp->SetMinimum(0);
       
       histo_tt_JERsysDown->Draw("h, sames");
       histo_tt_JERsysDown->SetMaximum(max * 1.5);
       histo_tt_JERsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max * 1.5);
       histo_twdr->SetMinimum(0);

       histo_twdr_JERsysUp->Draw("histo , sames");
       histo_twdr_JERsysUp->SetMaximum(max * 1.5);
       histo_twdr_JERsysUp->SetMinimum(0);

       histo_twdr_JERsysDown->Draw("histo , sames");
       histo_twdr_JERsysDown->SetMaximum(max * 1.5);
       histo_twdr_JERsysDown->SetMinimum(0);

       
       // other
       histo_other->Draw("histo , sames");
       histo_other->SetMaximum(max * 1.5);
       histo_other->SetMinimum(0);

       histo_other_JERsysUp->Draw("histo , sames");
       histo_other_JERsysUp->SetMaximum(max * 1.5);
       histo_other_JERsysUp->SetMinimum(0);

       histo_other_JERsysDown->Draw("histo , sames");
       histo_other_JERsysDown->SetMaximum(max * 1.5);
       histo_other_JERsysDown->SetMinimum(0);
       
       
       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotJER+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + ".png");
     // c1->SaveAs("plots/pdf/" + plotExtension + plotJER+ modeString[mode]+ "_region" + region + "_" + cutLabel[iPlots] + ".pdf");
    
    
        TCanvas *c1_ZSF = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max_ZSF * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
      histo_tt->Draw("h");
      
       
       histo_tt_ZSFsysUp->SetMaximum(max_ZSF * 1.5);
       histo_tt_ZSFsysUp->SetMinimum(0);
       histo_tt_ZSFsysUp->Draw("h, sames");
       
       
       histo_tt_ZSFsysDown->SetMaximum(max_ZSF * 1.5);
       histo_tt_ZSFsysDown->SetMinimum(0);     
       histo_tt_ZSFsysDown->Draw("h, sames");
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max_ZSF * 1.5);
       histo_twdr->SetMinimum(0);

       histo_twdr_ZSFsysUp->Draw("histo , sames");
       histo_twdr_ZSFsysUp->SetMaximum(max_ZSF * 1.5);
       histo_twdr_ZSFsysUp->SetMinimum(0);

       histo_twdr_ZSFsysDown->Draw("histo , sames");
       histo_twdr_ZSFsysDown->SetMaximum(max_ZSF * 1.5);
       histo_twdr_ZSFsysDown->SetMinimum(0);

              // other
       histo_other->Draw("histo , sames");
       histo_other->SetMaximum(max_ZSF * 1.5);
       histo_other->SetMinimum(0);

       histo_other_ZSFsysUp->Draw("histo , sames");
       histo_other_ZSFsysUp->SetMaximum(max_ZSF * 1.5);
       histo_other_ZSFsysUp->SetMinimum(0);

       histo_other_ZSFsysDown->Draw("histo , sames");
       histo_other_ZSFsysDown->SetMaximum(max_ZSF * 1.5);
       histo_other_ZSFsysDown->SetMinimum(0);
       
       
      legZSF->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_ZSF->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotZSF+ modeString[mode]+ "_region" + regionString[region] + "_" + cutLabel[iPlots] + ".png");
     // c1_ZSF->SaveAs("plots/pdf/" + plotExtension + plotZSF +modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + ".pdf");   
      
      
        TCanvas *c1_JES = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max_JES * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->Draw("h");
       
       histo_tt_JESsysUp->Draw("h, sames");
       histo_tt_JESsysUp->SetMaximum(max_JES * 1.5);
       histo_tt_JESsysUp->SetMinimum(0);
       
       histo_tt_JESsysDown->Draw("h, sames");
       histo_tt_JESsysDown->SetMaximum(max_JES * 1.5);
       histo_tt_JESsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max_JES * 1.5);
       histo_twdr->SetMinimum(0);

       histo_twdr_JESsysUp->Draw("histo , sames");
       histo_twdr_JESsysUp->SetMaximum(max_JES * 1.5);
       histo_twdr_JESsysUp->SetMinimum(0);

       histo_twdr_JESsysDown->Draw("histo , sames");
       histo_twdr_JESsysDown->SetMaximum(max_JES * 1.5);
       histo_twdr_JESsysDown->SetMinimum(0);

              // other
       histo_other->Draw("histo , sames");
       histo_other->SetMaximum(max_JES * 1.5);
       histo_other->SetMinimum(0);

       histo_other_JESsysUp->Draw("histo , sames");
       histo_other_JESsysUp->SetMaximum(max_JES * 1.5);
       histo_other_JESsysUp->SetMinimum(0);

       histo_other_JESsysDown->Draw("histo , sames");
       histo_other_JESsysDown->SetMaximum(max_JES * 1.5);
       histo_other_JESsysDown->SetMinimum(0);
       
       
      legJES->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_JES->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotJES+ modeString[mode]+ "_region" + regionString[region] + "_" + cutLabel[iPlots] + ".png");
     // c1_JES->SaveAs("plots/pdf/" + plotExtension + plotJES +modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + ".pdf");
     
     
     TCanvas *c1_PU = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max_PU * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->Draw("h");
       
       histo_tt_PUsysUp->Draw("h, sames");
       histo_tt_PUsysUp->SetMaximum(max_PU * 1.5);
       histo_tt_PUsysUp->SetMinimum(0);
       
       histo_tt_PUsysDown->Draw("h, sames");
       histo_tt_PUsysDown->SetMaximum(max_PU * 1.5);
       histo_tt_PUsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max_PU * 1.5);
       histo_twdr->SetMinimum(0);

       histo_twdr_PUsysUp->Draw("histo , sames");
       histo_twdr_PUsysUp->SetMaximum(max_PU * 1.5);
       histo_twdr_PUsysUp->SetMinimum(0);

       histo_twdr_PUsysDown->Draw("histo , sames");
       histo_twdr_PUsysDown->SetMaximum(max_PU * 1.5);
       histo_twdr_PUsysDown->SetMinimum(0);
       
       
       // other
       histo_other->Draw("histo , sames");
       histo_other->SetMaximum(max_PU * 1.5);
       histo_other->SetMinimum(0);

       histo_other_PUsysUp->Draw("histo , sames");
       histo_other_PUsysUp->SetMaximum(max_PU * 1.5);
       histo_other_PUsysUp->SetMinimum(0);

       histo_other_PUsysDown->Draw("histo , sames");
       histo_other_PUsysDown->SetMaximum(max_PU * 1.5);
       histo_other_PUsysDown->SetMinimum(0);

       
      legPU->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_PU->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotPU+ modeString[mode]+ "_region" + regionString[region] + "_" + cutLabel[iPlots] + ".png");
     // c1_PU->SaveAs("plots/pdf/" + plotExtension + plotPU +modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + ".pdf");
 
        TCanvas *c1_SF = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max_SF * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->Draw("h");
       
       histo_tt_SFsysUp->Draw("h, sames");
       histo_tt_SFsysUp->SetMaximum(max_SF * 1.5);
       histo_tt_SFsysUp->SetMinimum(0);
       
       histo_tt_SFsysDown->Draw("h, sames");
       histo_tt_SFsysDown->SetMaximum(max_SF * 1.5);
       histo_tt_SFsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max_SF * 1.5);
       histo_twdr->SetMinimum(0);

       histo_twdr_SFsysUp->Draw("histo , sames");
       histo_twdr_SFsysUp->SetMaximum(max_SF * 1.5);
       histo_twdr_SFsysUp->SetMinimum(0);

       histo_twdr_SFsysDown->Draw("histo , sames");
       histo_twdr_SFsysDown->SetMaximum(max_SF * 1.5);
       histo_twdr_SFsysDown->SetMinimum(0);
       
       // other
       histo_other->Draw("histo , sames");
       histo_other->SetMaximum(max_SF * 1.5);
       histo_other->SetMinimum(0);

       histo_other_SFsysUp->Draw("histo , sames");
       histo_other_SFsysUp->SetMaximum(max_SF * 1.5);
       histo_other_SFsysUp->SetMinimum(0);

       histo_other_SFsysDown->Draw("histo , sames");
       histo_other_SFsysDown->SetMaximum(max_SF * 1.5);
       histo_other_SFsysDown->SetMinimum(0);

       
      legSF->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_SF->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotSF+ modeString[mode]+ "_region" + regionString[region] + "_" + cutLabel[iPlots] + ".png");
     // c1_SF->SaveAs("plots/pdf/" + plotExtension + plotSF +modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + ".pdf");
     
     
     TCanvas *c1_MET = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max_MET * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->Draw("h");
       
       histo_tt_METsysUp->Draw("h, sames");
       histo_tt_METsysUp->SetMaximum(max_MET * 1.5);
       histo_tt_METsysUp->SetMinimum(0);
       
       histo_tt_METsysDown->Draw("h, sames");
       histo_tt_METsysDown->SetMaximum(max_MET * 1.5);
       histo_tt_METsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max_MET * 1.5);
       histo_twdr->SetMinimum(0);

       histo_twdr_METsysUp->Draw("histo , sames");
       histo_twdr_METsysUp->SetMaximum(max_MET * 1.5);
       histo_twdr_METsysUp->SetMinimum(0);

       histo_twdr_METsysDown->Draw("histo , sames");
       histo_twdr_METsysDown->SetMaximum(max_MET * 1.5);
       histo_twdr_METsysDown->SetMinimum(0);
       
       // other
       histo_other->Draw("histo , sames");
       histo_other->SetMaximum(max_MET * 1.5);
       histo_other->SetMinimum(0);

       histo_other_METsysUp->Draw("histo , sames");
       histo_other_METsysUp->SetMaximum(max_MET * 1.5);
       histo_other_METsysUp->SetMinimum(0);

       histo_other_METsysDown->Draw("histo , sames");
       histo_other_METsysDown->SetMaximum(max_MET * 1.5);
       histo_other_METsysDown->SetMinimum(0);

       
      legMET->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_MET->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotMET+ modeString[mode] + "_region" + regionString[region] + "_"+ cutLabel[iPlots] + ".png");
     // c1_MET->SaveAs("plots/pdf/" + plotExtension + plotMET +modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");      
      
      
    TCanvas *c1_TopMass = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max_TopMass * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->Draw("h");
    
       histo_tt_TopMassUp->Draw("h, sames");
       histo_tt_TopMassUp->SetMaximum(max_TopMass * 1.5);
       histo_tt_TopMassUp->SetMinimum(0);
       
       histo_tt_TopMassDown->Draw("h, sames");
       histo_tt_TopMassDown->SetMaximum(max_TopMass * 1.5);
       histo_tt_TopMassDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max_TopMass * 1.5);
       histo_twdr->SetMinimum(0);

       histo_twdr_TopMassUp->Draw("histo , sames");
       histo_twdr_TopMassUp->SetMaximum(max_TopMass * 1.5);
       histo_twdr_TopMassUp->SetMinimum(0);

       histo_twdr_TopMassDown->Draw("histo , sames");
       histo_twdr_TopMassDown->SetMaximum(max_TopMass * 1.5);
       histo_twdr_TopMassDown->SetMinimum(0);

       
      legTopMass->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_TopMass->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotTopMass+ modeString[mode] + "_region"  + regionString[region]+ "_" +cutLabel[iPlots] + ".png");
     // c1_TopMass->SaveAs("plots/pdf/" + plotExtension + plotTopMass +modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");         
      
  TCanvas *c1_Q2 = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max_Q2 * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->Draw("h");
    
       histo_tt_Q2Up->Draw("h, sames");
       histo_tt_Q2Up->SetMaximum(max_Q2 * 1.5);
       histo_tt_Q2Up->SetMinimum(0);
       
       histo_tt_Q2Down->Draw("h, sames");
       histo_tt_Q2Down->SetMaximum(max_Q2 * 1.5);
       histo_tt_Q2Down->SetMinimum(0);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max_Q2 * 1.5);
       histo_twdr->SetMinimum(0);

       histo_twdr_Q2Up->Draw("histo , sames");
       histo_twdr_Q2Up->SetMaximum(max_Q2 * 1.5);
       histo_twdr_Q2Up->SetMinimum(0);

       histo_twdr_Q2Down->Draw("histo , sames");
       histo_twdr_Q2Down->SetMaximum(max_Q2 * 1.5);
       histo_twdr_Q2Down->SetMinimum(0);

       
      legQ2->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_Q2->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotQ2+ modeString[mode] + "_region"  + regionString[region]+ "_" +cutLabel[iPlots] + ".png");
     // c1_Q2->SaveAs("plots/pdf/" + plotExtension + plotQ2 +modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");           
      
 TCanvas *c1_eleSF = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max_eleSF * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
      histo_tt->Draw("h");
      
       histo_tt_eleSFsysUp->Draw("h, sames");
       histo_tt_eleSFsysUp->SetMaximum(max_eleSF * 1.5);
       histo_tt_eleSFsysUp->SetMinimum(0);
       
       histo_tt_eleSFsysDown->Draw("h, sames");
       histo_tt_eleSFsysDown->SetMaximum(max_eleSF * 1.5);
       histo_tt_eleSFsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max_eleSF * 1.5);
       histo_twdr->SetMinimum(0);

       histo_twdr_eleSFsysUp->Draw("histo , sames");
       histo_twdr_eleSFsysUp->SetMaximum(max_eleSF * 1.5);
       histo_twdr_eleSFsysUp->SetMinimum(0);

       histo_twdr_eleSFsysDown->Draw("histo , sames");
       histo_twdr_eleSFsysDown->SetMaximum(max_eleSF * 1.5);
       histo_twdr_eleSFsysDown->SetMinimum(0);
       
       // other
       histo_other->Draw("histo , sames");
       histo_other->SetMaximum(max_eleSF * 1.5);
       histo_other->SetMinimum(0);

       histo_other_eleSFsysUp->Draw("histo , sames");
       histo_other_eleSFsysUp->SetMaximum(max_eleSF * 1.5);
       histo_other_eleSFsysUp->SetMinimum(0);

       histo_other_eleSFsysDown->Draw("histo , sames");
       histo_other_eleSFsysDown->SetMaximum(max_eleSF * 1.5);
       histo_other_eleSFsysDown->SetMinimum(0);

       
      legeleSF->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_eleSF->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + ploteleSF+ modeString[mode] + "_region"  + regionString[region]+ "_" +cutLabel[iPlots] + ".png");
     // c1_eleSF->SaveAs("plots/pdf/" + plotExtension + ploteleSF +modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");            
      
 TCanvas *c1_matching = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max_matching * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->Draw("h");
    
       histo_tt_matchingUp->Draw("h, sames");
       histo_tt_matchingUp->SetMaximum(max_matching * 1.5);
       histo_tt_matchingUp->SetMinimum(0);
       
       histo_tt_matchingDown->Draw("h, sames");
       histo_tt_matchingDown->SetMaximum(max_matching * 1.5);
       histo_tt_matchingDown->SetMinimum(0);     
       

       
      legmatching->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_matching->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotmatching+ modeString[mode] + "_region"  + regionString[region]+ "_" +cutLabel[iPlots] + ".png");
     // c1_matching->SaveAs("plots/pdf/" + plotExtension + plotmatching +modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");      
     
    
     TCanvas *c1_DS = new TCanvas();
       // ttbar
       
       histo_twds->SetMaximum(max_DS * 1.5);
       histo_twds->SetMinimum(0);
       histo_twds->GetYaxis()->SetTitle("#evts");
       histo_twds->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_twds->Draw("h");
    
   
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max_DS * 1.5);
       histo_twdr->SetMinimum(0);


       
      legDS->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_DS->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotDS+ modeString[mode] + "_region"  + regionString[region]+ "_" +cutLabel[iPlots] + ".png");
     // c1_DS->SaveAs("plots/pdf/" + plotExtension + plotDS +modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");  
    
    
    ///////////////////////////////////////
    /////NORMALIZED ///////////////////////
    //////////////////////////////////////
    
    
           TCanvas *c2 = new TCanvas();
       // ttbar
       
      histo_tt->SetMaximum(max * 1.5);
       //histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->DrawNormalized("h,1");
    
       histo_tt_JERsysUp->DrawNormalized("h, sames,1");
      histo_tt_JERsysUp->SetMaximum(max * 1.5);
     //  histo_tt_JERsysUp->SetMinimum(0);
       
       histo_tt_JERsysDown->DrawNormalized("h, sames,1");
      histo_tt_JERsysDown->SetMaximum(max * 1.5);
     //  histo_tt_JERsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->DrawNormalized("histo , sames,1");
      histo_twdr->SetMaximum(max * 1.5);
      // histo_twdr->SetMinimum(0);

       histo_twdr_JERsysUp->DrawNormalized("histo , sames,1");
      histo_twdr_JERsysUp->SetMaximum(max * 1.5);
       histo_twdr_JERsysUp->SetMinimum(0);

       histo_twdr_JERsysDown->DrawNormalized("histo , sames,1");
       histo_twdr_JERsysDown->SetMaximum(max * 1.5);
     //  histo_twdr_JERsysDown->SetMinimum(0);

      aleg->Draw();
      alabelcms->Draw();
      alabelcms2->Draw();
    
      c2->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotJER+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized.png");    
     
     
       // other
       histo_other->DrawNormalized("histo , sames,1");
       histo_other->SetMaximum(max * 1.5);
      // histo_other->SetMinimum(0);

       histo_other_JERsysUp->DrawNormalized("histo , sames,1");
       histo_other_JERsysUp->SetMaximum(max * 1.5);
    //   histo_other_JERsysUp->SetMinimum(0);

       histo_other_JERsysDown->DrawNormalized("histo , sames,1");
       histo_other_JERsysDown->SetMaximum(max * 1.5);
    //   histo_other_JERsysDown->SetMinimum(0);
       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c2->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotJER+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized_other.png");  
       

      TCanvas *c2_ZSF = new TCanvas();
	      
	             // other
       histo_other->GetYaxis()->SetTitle("#evts");
       histo_other->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_other->SetMaximum(max_ZSF*1.5);
      // histo_other->SetMinimum(0);
      histo_other->DrawNormalized("histo , 1");

       histo_other_ZSFsysUp->DrawNormalized("histo , sames,1");
       histo_other_ZSFsysUp->SetMaximum(max_ZSF * 1.5);
    //   histo_other_ZSFsysUp->SetMinimum(0);

       histo_other_ZSFsysDown->DrawNormalized("histo , sames,1");
       histo_other_ZSFsysDown->SetMaximum(max_ZSF * 1.5);
    //   histo_other_ZSFsysDown->SetMinimum(0);
       

      alegZSF->Draw();
      alabelcms->Draw();
      alabelcms2->Draw();
    
      c2_ZSF->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotZSF+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized_other.png");    
     
     
// ttbar
       histo_tt->DrawNormalized("h,sames,1");
      histo_tt->SetMaximum(max_ZSF* 1.5);
       //histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
    
       histo_tt_ZSFsysUp->DrawNormalized("h, sames,1");
      histo_tt_ZSFsysUp->SetMaximum(max_ZSF * 1.5);
     //  histo_tt_ZSFsysUp->SetMinimum(0);
       
       histo_tt_ZSFsysDown->DrawNormalized("h, sames,1");
      histo_tt_ZSFsysDown->SetMaximum(max_ZSF * 1.5);
     //  histo_tt_ZSFsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->DrawNormalized("histo , sames,1");
      histo_twdr->SetMaximum(max_ZSF * 1.5);
      // histo_twdr->SetMinimum(0);

       histo_twdr_ZSFsysUp->DrawNormalized("histo , sames,1");
      histo_twdr_ZSFsysUp->SetMaximum(max_ZSF * 1.5);
      // histo_twdr_ZSFsysUp->SetMinimum(0);

       histo_twdr_ZSFsysDown->DrawNormalized("histo , sames,1");
       histo_twdr_ZSFsysDown->SetMaximum(max_ZSF * 1.5);
     //  histo_twdr_ZSFsysDown->SetMinimum(0);
       
      legZSF->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c2_ZSF->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotZSF+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized.png");  
   
    
      
      
           TCanvas *c2_JES = new TCanvas();
       // ttbar
       
      histo_tt->SetMaximum(max * 1.5);
       //histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->DrawNormalized("h,1");
    
       histo_tt_JESsysUp->DrawNormalized("h, sames,1");
      histo_tt_JESsysUp->SetMaximum(max * 1.5);
     //  histo_tt_JESsysUp->SetMinimum(0);
       
       histo_tt_JESsysDown->DrawNormalized("h, sames,1");
      histo_tt_JESsysDown->SetMaximum(max * 1.5);
     //  histo_tt_JESsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->DrawNormalized("histo , sames,1");
      histo_twdr->SetMaximum(max * 1.5);
      // histo_twdr->SetMinimum(0);

       histo_twdr_JESsysUp->DrawNormalized("histo , sames,1");
      histo_twdr_JESsysUp->SetMaximum(max * 1.5);
      // histo_twdr_JESsysUp->SetMinimum(0);

       histo_twdr_JESsysDown->DrawNormalized("histo , sames,1");
       histo_twdr_JESsysDown->SetMaximum(max * 1.5);
     //  histo_twdr_JESsysDown->SetMinimum(0);

      alegJES->Draw();
      alabelcms->Draw();
      alabelcms2->Draw();
    
      c2_JES->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotJES+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized.png");    
     
     
       // other
       histo_other->DrawNormalized("histo , sames,1");
       histo_other->SetMaximum(max * 1.5);
      // histo_other->SetMinimum(0);

       histo_other_JESsysUp->DrawNormalized("histo , sames,1");
       histo_other_JESsysUp->SetMaximum(max * 1.5);
    //   histo_other_JESsysUp->SetMinimum(0);

       histo_other_JESsysDown->DrawNormalized("histo , sames,1");
       histo_other_JESsysDown->SetMaximum(max * 1.5);
    //   histo_other_JESsysDown->SetMinimum(0);
       
      legJES->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c2_JES->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotJES+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized_other.png");  
   
   
     
              TCanvas *c2_PU = new TCanvas();
       // ttbar
       
      histo_tt->SetMaximum(max * 1.5);
       //histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->DrawNormalized("h,1");
    
       histo_tt_PUsysUp->DrawNormalized("h, sames,1");
      histo_tt_PUsysUp->SetMaximum(max * 1.5);
     //  histo_tt_PUsysUp->SetMinimum(0);
       
       histo_tt_PUsysDown->DrawNormalized("h, sames,1");
      histo_tt_PUsysDown->SetMaximum(max * 1.5);
     //  histo_tt_PUsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->DrawNormalized("histo , sames,1");
      histo_twdr->SetMaximum(max * 1.5);
      // histo_twdr->SetMinimum(0);

       histo_twdr_PUsysUp->DrawNormalized("histo , sames,1");
      histo_twdr_PUsysUp->SetMaximum(max * 1.5);
       histo_twdr_PUsysUp->SetMinimum(0);

       histo_twdr_PUsysDown->DrawNormalized("histo , sames,1");
       histo_twdr_PUsysDown->SetMaximum(max * 1.5);
     //  histo_twdr_PUsysDown->SetMinimum(0);

      aleg->Draw();
      alabelcms->Draw();
      alabelcms2->Draw();
    
      c2_PU->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotPU+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized.png");    
     
     
       // other
       histo_other->DrawNormalized("histo , sames,1");
       histo_other->SetMaximum(max * 1.5);
      // histo_other->SetMinimum(0);

       histo_other_PUsysUp->DrawNormalized("histo , sames,1");
       histo_other_PUsysUp->SetMaximum(max * 1.5);
    //   histo_other_PUsysUp->SetMinimum(0);

       histo_other_PUsysDown->DrawNormalized("histo , sames,1");
       histo_other_PUsysDown->SetMaximum(max * 1.5);
    //   histo_other_PUsysDown->SetMinimum(0);
       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c2_PU->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotPU+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized_other.png");  
  
             TCanvas *c2_SF = new TCanvas();
       // ttbar
       
      histo_tt->SetMaximum(max * 1.5);
       //histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->DrawNormalized("h,1");
    
       histo_tt_SFsysUp->DrawNormalized("h, sames,1");
      histo_tt_SFsysUp->SetMaximum(max * 1.5);
     //  histo_tt_SFsysUp->SetMinimum(0);
       
       histo_tt_SFsysDown->DrawNormalized("h, sames,1");
      histo_tt_SFsysDown->SetMaximum(max * 1.5);
     //  histo_tt_SFsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->DrawNormalized("histo , sames,1");
      histo_twdr->SetMaximum(max * 1.5);
      // histo_twdr->SetMinimum(0);

       histo_twdr_SFsysUp->DrawNormalized("histo , sames,1");
      histo_twdr_SFsysUp->SetMaximum(max * 1.5);
       histo_twdr_SFsysUp->SetMinimum(0);

       histo_twdr_SFsysDown->DrawNormalized("histo , sames,1");
       histo_twdr_SFsysDown->SetMaximum(max * 1.5);
     //  histo_twdr_SFsysDown->SetMinimum(0);

      aleg->Draw();
      alabelcms->Draw();
      alabelcms2->Draw();
    
      c2_SF->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotSF+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized.png");    
     
     
       // other
       histo_other->DrawNormalized("histo , sames,1");
       histo_other->SetMaximum(max * 1.5);
      // histo_other->SetMinimum(0);

       histo_other_SFsysUp->DrawNormalized("histo , sames,1");
       histo_other_SFsysUp->SetMaximum(max * 1.5);
    //   histo_other_SFsysUp->SetMinimum(0);

       histo_other_SFsysDown->DrawNormalized("histo , sames,1");
       histo_other_SFsysDown->SetMaximum(max * 1.5);
    //   histo_other_SFsysDown->SetMinimum(0);
       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c2_SF->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotSF+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized_other.png");  
  
     
             TCanvas *c2_MET = new TCanvas();
       // ttbar
       
      histo_tt->SetMaximum(max * 1.5);
       //histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->DrawNormalized("h,1");
    
       histo_tt_METsysUp->DrawNormalized("h, sames,1");
      histo_tt_METsysUp->SetMaximum(max * 1.5);
     //  histo_tt_METsysUp->SetMinimum(0);
       
       histo_tt_METsysDown->DrawNormalized("h, sames,1");
      histo_tt_METsysDown->SetMaximum(max * 1.5);
     //  histo_tt_METsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->DrawNormalized("histo , sames,1");
      histo_twdr->SetMaximum(max * 1.5);
      // histo_twdr->SetMinimum(0);

       histo_twdr_METsysUp->DrawNormalized("histo , sames,1");
      histo_twdr_METsysUp->SetMaximum(max * 1.5);
       histo_twdr_METsysUp->SetMinimum(0);

       histo_twdr_METsysDown->DrawNormalized("histo , sames,1");
       histo_twdr_METsysDown->SetMaximum(max * 1.5);
     //  histo_twdr_METsysDown->SetMinimum(0);

      aleg->Draw();
      alabelcms->Draw();
      alabelcms2->Draw();
    
      c2_MET->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotMET+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized.png");    
     
     
       // other
       histo_other->DrawNormalized("histo , sames,1");
       histo_other->SetMaximum(max * 1.5);
      // histo_other->SetMinimum(0);

       histo_other_METsysUp->DrawNormalized("histo , sames,1");
       histo_other_METsysUp->SetMaximum(max * 1.5);
    //   histo_other_METsysUp->SetMinimum(0);

       histo_other_METsysDown->DrawNormalized("histo , sames,1");
       histo_other_METsysDown->SetMaximum(max * 1.5);
    //   histo_other_METsysDown->SetMinimum(0);
       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c2_MET->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotMET+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized_other.png");  
  
      
             TCanvas *c2_TopMass = new TCanvas();
       // ttbar
       
      histo_tt->SetMaximum(max * 1.5);
       //histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->DrawNormalized("h,1");
    
       histo_tt_TopMassUp->DrawNormalized("h, sames,1");
      histo_tt_TopMassUp->SetMaximum(max * 1.5);
     //  histo_tt_TopMassUp->SetMinimum(0);
       
       histo_tt_TopMassDown->DrawNormalized("h, sames,1");
      histo_tt_TopMassDown->SetMaximum(max * 1.5);
     //  histo_tt_TopMassDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->DrawNormalized("histo , sames,1");
      histo_twdr->SetMaximum(max * 1.5);
      // histo_twdr->SetMinimum(0);

       histo_twdr_TopMassUp->DrawNormalized("histo , sames,1");
      histo_twdr_TopMassUp->SetMaximum(max * 1.5);
       histo_twdr_TopMassUp->SetMinimum(0);

       histo_twdr_TopMassDown->DrawNormalized("histo , sames,1");
       histo_twdr_TopMassDown->SetMaximum(max * 1.5);
     //  histo_twdr_TopMassDown->SetMinimum(0);

      aleg->Draw();
      alabelcms->Draw();
      alabelcms2->Draw();
    
      c2_TopMass->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotTopMass+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized.png");    
     
     

      
  TCanvas *c2_Q2 = new TCanvas();
       // ttbar
       
       histo_tt->SetMaximum(max_Q2 * 1.5);
       histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->DrawNormalized("h");
    
       histo_tt_Q2Up->DrawNormalized("h, sames");
       histo_tt_Q2Up->SetMaximum(max_Q2 * 1.5);
       histo_tt_Q2Up->SetMinimum(0);
       
       histo_tt_Q2Down->DrawNormalized("h, sames");
       histo_tt_Q2Down->SetMaximum(max_Q2 * 1.5);
       histo_tt_Q2Down->SetMinimum(0);     
       
       // twdr
       histo_twdr->DrawNormalized("histo , sames");
       histo_twdr->SetMaximum(max_Q2 * 1.5);
       histo_twdr->SetMinimum(0);

       histo_twdr_Q2Up->DrawNormalized("histo , sames");
       histo_twdr_Q2Up->SetMaximum(max_Q2 * 1.5);
       histo_twdr_Q2Up->SetMinimum(0);

       histo_twdr_Q2Down->DrawNormalized("histo , sames");
       histo_twdr_Q2Down->SetMaximum(max_Q2 * 1.5);
       histo_twdr_Q2Down->SetMinimum(0);

       
      legQ2->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c2_Q2->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotQ2+ modeString[mode] + "_region"  + regionString[region]+ "_" +cutLabel[iPlots] + "normalized.png");
     // c2_Q2->SaveAs("plots/pdf/" + plotExtension + plotQ2 +modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");           
      
            TCanvas *c2_eleSF = new TCanvas();
       // ttbar
       
      histo_tt->SetMaximum(max * 1.5);
       //histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->DrawNormalized("h,1");
    
       histo_tt_eleSFsysUp->DrawNormalized("h, sames,1");
      histo_tt_eleSFsysUp->SetMaximum(max * 1.5);
     //  histo_tt_eleSFsysUp->SetMinimum(0);
       
       histo_tt_eleSFsysDown->DrawNormalized("h, sames,1");
      histo_tt_eleSFsysDown->SetMaximum(max * 1.5);
     //  histo_tt_eleSFsysDown->SetMinimum(0);     
       
       // twdr
       histo_twdr->DrawNormalized("histo , sames,1");
      histo_twdr->SetMaximum(max * 1.5);
      // histo_twdr->SetMinimum(0);

       histo_twdr_eleSFsysUp->DrawNormalized("histo , sames,1");
      histo_twdr_eleSFsysUp->SetMaximum(max * 1.5);
       histo_twdr_eleSFsysUp->SetMinimum(0);

       histo_twdr_eleSFsysDown->DrawNormalized("histo , sames,1");
       histo_twdr_eleSFsysDown->SetMaximum(max * 1.5);
     //  histo_twdr_eleSFsysDown->SetMinimum(0);

      aleg->Draw();
      alabelcms->Draw();
      alabelcms2->Draw();
    
      c2_eleSF->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + ploteleSF+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized.png");    
     
     
       // other
       histo_other->DrawNormalized("histo , sames,1");
       histo_other->SetMaximum(max * 1.5);
      // histo_other->SetMinimum(0);

       histo_other_eleSFsysUp->DrawNormalized("histo , sames,1");
       histo_other_eleSFsysUp->SetMaximum(max * 1.5);
    //   histo_other_eleSFsysUp->SetMinimum(0);

       histo_other_eleSFsysDown->DrawNormalized("histo , sames,1");
       histo_other_eleSFsysDown->SetMaximum(max * 1.5);
    //   histo_other_eleSFsysDown->SetMinimum(0);
       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c2_eleSF->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + ploteleSF+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized_other.png");  
      
      
            TCanvas *c2_matching = new TCanvas();
       // ttbar
       
      histo_tt->SetMaximum(max * 1.5);
       //histo_tt->SetMinimum(0);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       
    
       histo_tt_matchingUp->DrawNormalized("h, sames,1");
      histo_tt_matchingUp->SetMaximum(max * 1.5);
     //  histo_tt_matchingUp->SetMinimum(0);
       
       histo_tt_matchingDown->DrawNormalized("h, sames,1");
      histo_tt_matchingDown->SetMaximum(max * 1.5);
     //  histo_tt_matchingDown->SetMinimum(0);     
       

      aleg->Draw();
      alabelcms->Draw();
      alabelcms2->Draw();
    
      c2_matching->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotJES+ modeString[mode] + "_region" +  regionString[region]+ "_" + cutLabel[iPlots] + "normalized.png");    
     
     

    
     TCanvas *c2_DS = new TCanvas();
       // ttbar
       
       histo_twds->SetMaximum(max_DS * 1.5);
       histo_twds->SetMinimum(0);
       histo_twds->GetYaxis()->SetTitle("#evts");
       histo_twds->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_twds->DrawNormalized("h");
    
   
       // twdr
       histo_twdr->DrawNormalized("histo , sames");
       histo_twdr->SetMaximum(max_DS * 1.5);
       histo_twdr->SetMinimum(0);


       
      legDS->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c2_DS->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotDS+ modeString[mode] + "_region"  + regionString[region]+ "_" +cutLabel[iPlots] + "normalized.png");
     // c2_DS->SaveAs("plots/pdf/" + plotExtension + plotDS +modeString[mode] + "_" + cutLabel[iPlots] + ".pdf"); 
    
    
    
    
    
    
    

    } // end plots loop 
    
 
} // end constructor loop
