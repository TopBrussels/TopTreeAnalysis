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
  
  double lumi = 1000;
   if      (mode == 0)	 lumi = 11966.617;  	
    else if (mode == 1) lumi = 12067.294;  	
    else if (mode == 2) lumi = 12093.792;  	
    
  sprintf(myRootFile,"results/an_%dpb_%d.root", lumi, mode); // take output from looper
 
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
  
  
   
  TFile *_file0 = TFile::Open(myRootFile);
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
  
  cout << "-------------------------------------------------------" << endl; 
  cout << " ------------ USED FILES ------------------------------" << endl; 
  cout << "Normal: " << myRootFile << endl;
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
  cout << "-------------------------------------------------------" << endl; 
  
  const int nProcess = 2;
  const int nPlots = 10; // 16
  const int nSys = 11;


  TString processName[nProcess] =  { "twdr", "tt"};
  TString processTitle[nProcess] = { "tW","t#bar{t}"};
  Color_t color[nProcess] =        { kBlue, kRed};
  
  TString systName[nSys] = { "Normal", "JERsysDown" , "JERsysUp","JESsysDown" , "JESsysUp","PUsysDown" , "PUsysUp","SFsysDown" , "SFsysUp","METsysDown" , "METsysUp"}; 

/*
  TString cutLabel[nPlots] =     {  "met", "mll", "ptsys", "ht", "pt_leading", "nvertex", "met_2j1t", "mll_2j1t", "ptsys_2j1t", "ht_2j1t", "pt_leading_2j1t","met_2j2t", "mll_2j2t", "ptsys_2j2t", "ht_2j2t", "pt_leading_2j2t" };
  int rebinHisto[nPlots] =       {  4, 4,  4, 12, 4, 1,4, 4,  4, 12, 4, 4, 4,  4, 12, 4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} ", "Inv. Mass ", "P_{T} system [GeV] ", "H_{T} [GeV] ","P_{T} of the leading jet ", "# of vertex ", "E_{T}^{miss} 2j1t", "Inv. Mass 2j1t", "P_{T} system [GeV] 2j1t", "H_{T} [GeV] 2j1t","P_{T} of the leading jet 2j1t","E_{T}^{miss} 2j2t", "Inv. Mass 2j2t", "P_{T} system [GeV] 2j2t", "H_{T} [GeV] 2j2t","P_{T} of the leading jet 2j2t" };
*/
  if(region == 0){
  TString cutLabel[nPlots] =     {  "met_1j1t", "mll_1j1t", "ptsys_1j1t", "ht_1j1t","pt_leading_1j1t",  "nvertex_1j1t", "pt_max_1j1t","pt_min_1j1t","ht_nomet_1j1t","eta_leading_1j1t" };
  int rebinHisto[nPlots] =       {  4, 4,  4, 12, 4, 1,4,4,12,4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} 1j1t", "Inv. Mass 1j1t", "P_{T} system [GeV] 1j1t","H_{T} [GeV] 1j1t","P_{T} of the leading jet 1j1t", "# of vertex 1j1t", "p_T of the first lepton [GeV] 1j1t", "p_T  of the second lepton [GeV] 1j1t", "H_{T} no met [GeV] 1j1t","Eta of the leading jet 1j1t"};
  } else if(region == 1){
   TString cutLabel[nPlots] =     {  "met_2j1t", "mll_2j1t", "ptsys_2j1t", "ht_2j1t","pt_leading_2j1t",  "nvertex_2j1t", "pt_max_2j1t","pt_min_2j1t","ht_nomet_2j1t","eta_leading_2j1t" };
  int rebinHisto[nPlots] =       {  4, 4,  4, 12, 4, 1,4,4,12,4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} 2j1t", "Inv. Mass 2j1t", "P_{T} system [GeV] 2j1t","H_{T} [GeV] 2j1t","P_{T} of the leading jet 2j1t", "# of vertex 2j1t", "p_T of the first lepton [GeV] 2j1t", "p_T  of the second lepton [GeV] 2j1t", "H_{T} no met [GeV] 2j1t","Eta of the leading jet 2j1t"};
  }else if(region == 2){
   TString cutLabel[nPlots] =     {  "met_2j2t", "mll_2j2t", "ptsys_2j2t", "ht_2j2t","pt_leading_2j2t",  "nvertex_2j2t", "pt_max_2j2t","pt_min_2j2t","ht_nomet_2j2t","eta_leading_2j2t" };
  int rebinHisto[nPlots] =       {  4, 4,  4, 12, 4, 1,4,4,12,4};
  TString cutTitle[nPlots] =     { "E_{T}^{miss} 2j2t", "Inv. Mass 2j2t", "P_{T} system [GeV] 2j2t","H_{T} [GeV] 2j2t","P_{T} of the leading jet 2j2t", "# of vertex 2j2t", "p_T of the first lepton [GeV] 2j2t", "p_T  of the second lepton [GeV] 2j2t", "H_{T} no met [GeV] 2j2t","Eta of the leading jet 2j2t"};
 
  
  }

  TString modeString[3] = {"0","1","2"};
  
  TString plotExtension = "sys_plot_"; // name of the plots
  TString plotJER = "JER_"; 
  TString plotJES = "JES_"; 
  TString plotPU = "PU_";
  TString plotSF = "SF_";
  TString plotMET = "MET_";
  TString plotAnalysis = "Systematics"; // directory in plots where the plots are saved
  

  
   
  
  
  //Make plots   
  TH1F* histo_tt;
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
  
  TH1F* histo_twdr;
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
   
  
   for(int iPlots = 0; iPlots< nPlots; iPlots++)
   {
      
     leg = new TLegend(0.7,0.7,0.94,0.94);
     leg ->SetFillStyle(1001);
     leg ->SetFillColor(kWhite);
     leg ->SetBorderSize(1);
     
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
 
	
       // Get histos tt bar 
       histo_tt= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[1]);
	      // cout << "Normal: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
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
	       
	       	       
       histo_tt->Rebin(rebinHisto[iPlots]);
       histo_tt->SetLineColor(color[0]);
       histo_tt->SetLineWidth(2);
       
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
       
       // get histos twdr
       histo_twdr= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[0]);
        //cout << "Normal: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
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
	
	
       histo_twdr->Rebin(rebinHisto[iPlots]);
              // histo_twdr->SetFillColor(color[0]);
       histo_twdr->SetLineColor(color[1]);
       histo_twdr->SetLineWidth(2);
       
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
        
       
       
	// Make legends      
       leg ->AddEntry(histo_tt, "tt bckgr." , "l");  
       leg ->AddEntry(histo_tt_JERsysDown, "tt bckgr. JERsysDown" , "l"); 
       leg ->AddEntry(histo_tt_JERsysUp, "tt bckgr. JERsysUp" , "l"); 
       
       leg ->AddEntry(histo_twdr, "twdr signal", "l"); 
       leg ->AddEntry(histo_twdr_JERsysDown, "twdr signal JERsysDown", "l"); 
       leg ->AddEntry(histo_twdr_JERsysUp, "twdr signal JERsysUp", "l");
       
       legJES ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legJES ->AddEntry(histo_tt_JESsysDown, "tt bckgr. JESsysDown" , "l"); 
       legJES ->AddEntry(histo_tt_JESsysUp, "tt bckgr. JESsysUp" , "l");   

       legJES ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legJES ->AddEntry(histo_twdr_JESsysDown, "twdr signal JESsysDown", "l"); 
       legJES ->AddEntry(histo_twdr_JESsysUp, "twdr signal JESsysUp", "l");
       
       legPU ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legPU ->AddEntry(histo_tt_PUsysDown, "tt bckgr. PUsysDown" , "l"); 
       legPU ->AddEntry(histo_tt_PUsysUp, "tt bckgr. PUsysUp" , "l");   

       legPU ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legPU ->AddEntry(histo_twdr_PUsysDown, "twdr signal PUsysDown", "l"); 
       legPU ->AddEntry(histo_twdr_PUsysUp, "twdr signal PUsysUp", "l");          

       legSF ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legSF ->AddEntry(histo_tt_SFsysDown, "tt bckgr. SFsysDown" , "l"); 
       legSF ->AddEntry(histo_tt_SFsysUp, "tt bckgr. SFsysUp" , "l");   

       legSF ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legSF ->AddEntry(histo_twdr_SFsysDown, "twdr signal SFsysDown", "l"); 
       legSF ->AddEntry(histo_twdr_SFsysUp, "twdr signal SFsysUp", "l");
       
       legMET ->AddEntry(histo_tt, "tt bckgr." , "l");  
       legMET ->AddEntry(histo_tt_METsysDown, "tt bckgr. METsysDown" , "l"); 
       legMET ->AddEntry(histo_tt_METsysUp, "tt bckgr. METsysUp" , "l");   

       legMET ->AddEntry(histo_twdr, "twdr signal", "l"); 
       legMET ->AddEntry(histo_twdr_METsysDown, "twdr signal METsysDown", "l"); 
       legMET ->AddEntry(histo_twdr_METsysUp, "twdr signal METsysUp", "l"); 
              
       
       // determine maxima
       double max1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max2 = TMath::Max( histo_tt_JERsysUp->GetMaximum(), histo_twdr_JERsysUp->GetMaximum());
       double max3 = TMath::Max(histo_tt_JERsysDown->GetMaximum(), histo_twdr_JERsysDown->GetMaximum());
       double max4 = TMath::Max(max1,max2);
       double max = TMath::Max(max3,max4);
       
       double max_JES_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max_JES_2 = TMath::Max( histo_tt_JESsysUp->GetMaximum(), histo_twdr_JESsysUp->GetMaximum());
       double max_JES_3 = TMath::Max(histo_tt_JESsysDown->GetMaximum(), histo_twdr_JESsysDown->GetMaximum());
       double max_JES_4 = TMath::Max(max_JES_1,max_JES_2);
       double max_JES = TMath::Max(max_JES_3,max_JES_4);
       
       double max_PU_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max_PU_2 = TMath::Max( histo_tt_PUsysUp->GetMaximum(), histo_twdr_PUsysUp->GetMaximum());
       double max_PU_3 = TMath::Max(histo_tt_PUsysDown->GetMaximum(), histo_twdr_PUsysDown->GetMaximum());
       double max_PU_4 = TMath::Max(max_PU_1,max_PU_2);
       double max_PU = TMath::Max(max_PU_3,max_PU_4);

       double max_SF_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max_SF_2 = TMath::Max( histo_tt_SFsysUp->GetMaximum(), histo_twdr_SFsysUp->GetMaximum());
       double max_SF_3 = TMath::Max(histo_tt_SFsysDown->GetMaximum(), histo_twdr_SFsysDown->GetMaximum());
       double max_SF_4 = TMath::Max(max_SF_1,max_SF_2);
       double max_SF = TMath::Max(max_SF_3,max_SF_4);
       
       double max_MET_1 = TMath::Max(histo_tt->GetMaximum(), histo_twdr->GetMaximum());
       double max_MET_2 = TMath::Max( histo_tt_METsysUp->GetMaximum(), histo_twdr_METsysUp->GetMaximum());
       double max_MET_3 = TMath::Max(histo_tt_METsysDown->GetMaximum(), histo_twdr_METsysDown->GetMaximum());
       double max_MET_4 = TMath::Max(max_MET_1,max_MET_2);
       double max_MET = TMath::Max(max_MET_3,max_MET_4);    
    
       
       TCanvas *c1 = new TCanvas();
       // ttbar
       histo_tt->Draw("h");
       histo_tt->SetMaximum(max * 1.5);
       histo_tt->SetMinimum(1);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
    
       histo_tt_JERsysUp->Draw("h, sames");
       histo_tt_JERsysUp->SetMaximum(max * 1.5);
       histo_tt_JERsysUp->SetMinimum(1);
       
       histo_tt_JERsysDown->Draw("h, sames");
       histo_tt_JERsysDown->SetMaximum(max * 1.5);
       histo_tt_JERsysDown->SetMinimum(1);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max * 1.5);
       histo_twdr->SetMinimum(1);

       histo_twdr_JERsysUp->Draw("histo , sames");
       histo_twdr_JERsysUp->SetMaximum(max * 1.5);
       histo_twdr_JERsysUp->SetMinimum(1);

       histo_twdr_JERsysDown->Draw("histo , sames");
       histo_twdr_JERsysDown->SetMaximum(max * 1.5);
       histo_twdr_JERsysDown->SetMinimum(1);

       
      leg->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotJER+ modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + ".png");
     // c1->SaveAs("plots/pdf/" + plotExtension + plotJER+ modeString[mode]+ "_region" + region + "_" + cutLabel[iPlots] + ".pdf");
      
      
        TCanvas *c1_JES = new TCanvas();
       // ttbar
       histo_tt->Draw("h");
       histo_tt->SetMaximum(max * 1.5);
       histo_tt->SetMinimum(1);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
    
       histo_tt_JESsysUp->Draw("h, sames");
       histo_tt_JESsysUp->SetMaximum(max * 1.5);
       histo_tt_JESsysUp->SetMinimum(1);
       
       histo_tt_JESsysDown->Draw("h, sames");
       histo_tt_JESsysDown->SetMaximum(max * 1.5);
       histo_tt_JESsysDown->SetMinimum(1);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max * 1.5);
       histo_twdr->SetMinimum(1);

       histo_twdr_JESsysUp->Draw("histo , sames");
       histo_twdr_JESsysUp->SetMaximum(max * 1.5);
       histo_twdr_JESsysUp->SetMinimum(1);

       histo_twdr_JESsysDown->Draw("histo , sames");
       histo_twdr_JESsysDown->SetMaximum(max * 1.5);
       histo_twdr_JESsysDown->SetMinimum(1);

       
      legJES->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_JES->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotJES+ modeString[mode]+ "_region" + region + "_" + cutLabel[iPlots] + ".png");
     // c1_JES->SaveAs("plots/pdf/" + plotExtension + plotJES +modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + ".pdf");
     
     
     TCanvas *c1_PU = new TCanvas();
       // ttbar
       histo_tt->Draw("h");
       histo_tt->SetMaximum(max * 1.5);
       histo_tt->SetMinimum(1);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
    
       histo_tt_PUsysUp->Draw("h, sames");
       histo_tt_PUsysUp->SetMaximum(max * 1.5);
       histo_tt_PUsysUp->SetMinimum(1);
       
       histo_tt_PUsysDown->Draw("h, sames");
       histo_tt_PUsysDown->SetMaximum(max * 1.5);
       histo_tt_PUsysDown->SetMinimum(1);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max * 1.5);
       histo_twdr->SetMinimum(1);

       histo_twdr_PUsysUp->Draw("histo , sames");
       histo_twdr_PUsysUp->SetMaximum(max * 1.5);
       histo_twdr_PUsysUp->SetMinimum(1);

       histo_twdr_PUsysDown->Draw("histo , sames");
       histo_twdr_PUsysDown->SetMaximum(max * 1.5);
       histo_twdr_PUsysDown->SetMinimum(1);

       
      legPU->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_PU->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotPU+ modeString[mode]+ "_region" + region + "_" + cutLabel[iPlots] + ".png");
     // c1_PU->SaveAs("plots/pdf/" + plotExtension + plotPU +modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + ".pdf");
 
        TCanvas *c1_SF = new TCanvas();
       // ttbar
       histo_tt->Draw("h");
       histo_tt->SetMaximum(max * 1.5);
       histo_tt->SetMinimum(1);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
    
       histo_tt_SFsysUp->Draw("h, sames");
       histo_tt_SFsysUp->SetMaximum(max * 1.5);
       histo_tt_SFsysUp->SetMinimum(1);
       
       histo_tt_SFsysDown->Draw("h, sames");
       histo_tt_SFsysDown->SetMaximum(max * 1.5);
       histo_tt_SFsysDown->SetMinimum(1);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max * 1.5);
       histo_twdr->SetMinimum(1);

       histo_twdr_SFsysUp->Draw("histo , sames");
       histo_twdr_SFsysUp->SetMaximum(max * 1.5);
       histo_twdr_SFsysUp->SetMinimum(1);

       histo_twdr_SFsysDown->Draw("histo , sames");
       histo_twdr_SFsysDown->SetMaximum(max * 1.5);
       histo_twdr_SFsysDown->SetMinimum(1);

       
      legSF->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_SF->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotSF+ modeString[mode]+ "_region" + region + "_" + cutLabel[iPlots] + ".png");
     // c1_SF->SaveAs("plots/pdf/" + plotExtension + plotSF +modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + ".pdf");
     
     
     TCanvas *c1_MET = new TCanvas();
       // ttbar
       histo_tt->Draw("h");
       histo_tt->SetMaximum(max * 1.5);
       histo_tt->SetMinimum(1);
       histo_tt->GetYaxis()->SetTitle("#evts");
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
    
       histo_tt_METsysUp->Draw("h, sames");
       histo_tt_METsysUp->SetMaximum(max * 1.5);
       histo_tt_METsysUp->SetMinimum(1);
       
       histo_tt_METsysDown->Draw("h, sames");
       histo_tt_METsysDown->SetMaximum(max * 1.5);
       histo_tt_METsysDown->SetMinimum(1);     
       
       // twdr
       histo_twdr->Draw("histo , sames");
       histo_twdr->SetMaximum(max * 1.5);
       histo_twdr->SetMinimum(1);

       histo_twdr_METsysUp->Draw("histo , sames");
       histo_twdr_METsysUp->SetMaximum(max * 1.5);
       histo_twdr_METsysUp->SetMinimum(1);

       histo_twdr_METsysDown->Draw("histo , sames");
       histo_twdr_METsysDown->SetMaximum(max * 1.5);
       histo_twdr_METsysDown->SetMinimum(1);

       
      legMET->Draw();
      labelcms->Draw();
      labelcms2->Draw();
    
      c1_MET->SaveAs("plots/" + plotAnalysis +"/" + plotExtension + plotMET+ modeString[mode] + "_" + cutLabel[iPlots] + ".png");
     // c1_MET->SaveAs("plots/pdf/" + plotExtension + plotMET +modeString[mode] + "_" + cutLabel[iPlots] + ".pdf");      
      
      TCanvas *c2 = new TCanvas();
       histo_tt->DrawNormalized("h",1);
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->GetYaxis()->SetTitle("#evts");
       
       histo_tt_JERsysDown->DrawNormalized("h, sames",1); 
       histo_tt_JERsysUp->DrawNormalized("h, sames",1);
    
       
       
       histo_twdr->DrawNormalized("h,sames",1);
       histo_twdr_JERsysDown->DrawNormalized("h,sames",1);
       histo_twdr_JERsysUp->DrawNormalized("h,sames",1);
       
      leg->Draw();
     labelcms->Draw();
     labelcms2->Draw();
    
      c2->SaveAs("plots/"  + plotAnalysis +"/"  + plotExtension + plotJER+modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + "_normalized" +  ".png");
    //  c2->SaveAs("plots/pdf/" + plotExtension + plotJER + modeString[mode]+ "_region" + region + "_" + cutLabel[iPlots] + "_normalized" + ".pdf");
    
      TCanvas *c2_JES = new TCanvas();
       histo_tt->DrawNormalized("h",1);
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->GetYaxis()->SetTitle("#evts");
       
       histo_tt_JESsysDown->DrawNormalized("h, sames",1); 
       histo_tt_JESsysUp->DrawNormalized("h, sames",1);
    
       
       
       histo_twdr->DrawNormalized("h,sames",1);
       histo_twdr_JESsysDown->DrawNormalized("h,sames",1);
       histo_twdr_JESsysUp->DrawNormalized("h,sames",1);
       
      legJES->Draw();
     labelcms->Draw();
     labelcms2->Draw();
    
      c2_JES->SaveAs("plots/"  + plotAnalysis +"/"  + plotExtension +plotJES+ modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + "_normalized" +  ".png");
    //  c2_JES->SaveAs("plots/pdf/" + plotExtension + plotJES + modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + "_normalized" + ".pdf");
    
    TCanvas *c2_PU = new TCanvas();
       histo_tt->DrawNormalized("h",1);
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->GetYaxis()->SetTitle("#evts");
       
       histo_tt_PUsysDown->DrawNormalized("h, sames",1); 
       histo_tt_PUsysUp->DrawNormalized("h, sames",1);
    
       
       
       histo_twdr->DrawNormalized("h,sames",1);
       histo_twdr_PUsysDown->DrawNormalized("h,sames",1);
       histo_twdr_PUsysUp->DrawNormalized("h,sames",1);
       
      legPU->Draw();
     labelcms->Draw();
     labelcms2->Draw();
    
      c2_PU->SaveAs("plots/"  + plotAnalysis +"/"  + plotExtension +plotPU+ modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + "_normalized" +  ".png");
    //  c2_PU->SaveAs("plots/pdf/" + plotExtension + plotPU + modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + "_normalized" + ".pdf");
      
      TCanvas *c2_SF = new TCanvas();
       histo_tt->DrawNormalized("h",1);
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->GetYaxis()->SetTitle("#evts");
       
       histo_tt_SFsysDown->DrawNormalized("h, sames",1); 
       histo_tt_SFsysUp->DrawNormalized("h, sames",1);
    
       
       
       histo_twdr->DrawNormalized("h,sames",1);
       histo_twdr_SFsysDown->DrawNormalized("h,sames",1);
       histo_twdr_SFsysUp->DrawNormalized("h,sames",1);
       
      legSF->Draw();
     labelcms->Draw();
     labelcms2->Draw();
    
      c2_SF->SaveAs("plots/"  + plotAnalysis +"/"  + plotExtension +plotSF+ modeString[mode]+ "_region" + region + "_" + cutLabel[iPlots] + "_normalized" +  ".png");
    //  c2_SF->SaveAs("plots/pdf/" + plotExtension + plotSF + modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + "_normalized" + ".pdf");
    
    TCanvas *c2_MET = new TCanvas();
       histo_tt->DrawNormalized("h",1);
       histo_tt->GetXaxis()->SetTitle(cutTitle[iPlots]);
       histo_tt->GetYaxis()->SetTitle("#evts");
       
       histo_tt_METsysDown->DrawNormalized("h, sames",1); 
       histo_tt_METsysUp->DrawNormalized("h, sames",1);
    
       
       
       histo_twdr->DrawNormalized("h,sames",1);
       histo_twdr_METsysDown->DrawNormalized("h,sames",1);
       histo_twdr_METsysUp->DrawNormalized("h,sames",1);
       
      legMET->Draw();
     labelcms->Draw();
     labelcms2->Draw();
    
      c2_MET->SaveAs("plots/"  + plotAnalysis +"/"  + plotExtension +plotMET+ modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + "_normalized" +  ".png");
    //  c2_MET->SaveAs("plots/pdf/" + plotExtension + plotMET + modeString[mode] + "_region" + region+ "_" + cutLabel[iPlots] + "_normalized" + ".pdf");      

      
    } // end plots loop 
    
 
} // end constructor loop
