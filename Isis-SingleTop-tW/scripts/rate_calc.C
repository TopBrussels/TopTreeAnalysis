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

void rate_calc(int mode = 2, int region = 0){
 


  
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
  char myRootFilePDFtt[300];
  char myRootFilePDFtwdr[300];
  char myRootFilePDFother[300];
  
  
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
  
   if(region == 0){
 	  sprintf(myRootFilePDFtt,"renamingROOT/PDF_unc_results/pdf_1j1t_%d_tt_outputfile_rebinned.root", mode);
 	 sprintf(myRootFilePDFtwdr,"renamingROOT/PDF_unc_results/pdf_1j1t_%d_twdr_outputfile_rebinned.root", mode);
	  sprintf(myRootFilePDFother,"renamingROOT/PDF_unc_results/pdf_1j1t_%d_others_outputfile_rebinned.root", mode);
	  
   }else if(region == 1){
      sprintf(myRootFilePDFtt,"renamingROOT/PDF_unc_results/pdf_2j1t_%d_tt_outputfile_rebinned.root", mode);
	  sprintf(myRootFilePDFtwdr,"renamingROOT/PDF_unc_results/pdf_2j1t_%d_twdr_outputfile_rebinned.root", mode);
	  sprintf(myRootFilePDFother,"renamingROOT/PDF_unc_results/pdf_2j1t_%d_others_outputfile_rebinned.root", mode);
	  
   
   } else if (region == 2){http://mon.iihe.ac.be/~ivanpari/23may_rates/analysisJochenSyst8TeVRunTheta/
     sprintf(myRootFilePDFtt,"renamingROOT/PDF_unc_results/pdf_2j2t_%d_tt_outputfile_rebinned.root", mode);
	  sprintf(myRootFilePDFtwdr,"renamingROOT/PDF_unc_results/pdf_2j2t_%d_twdr_outputfile_rebinned.root", mode);
	  sprintf(myRootFilePDFother,"renamingROOT/PDF_unc_results/pdf_2j2t_%d_others_outputfile_rebinned.root", mode);
   
   
   }
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
  TFile *_filePDFtt = TFile::Open(myRootFilePDFtt);
  TFile *_filePDFtwdr = TFile::Open(myRootFilePDFtwdr);
  TFile *_filePDFother = TFile::Open(myRootFilePDFother);
 
  
  
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
  const int nPlots = 1; // 16
  const int nSys = 21;


  TString processName[nProcess] =  { "twdr", "tt","twds","others"};
  TString processTitle[nProcess] = { "tW DR","t#bar{t}","tW DS", "other"};
  Color_t color[nProcess] =        { kBlue, kRed, kGreen,kMagenta};
  
  TString systName[nSys] = { "Normal", "ZSFsysDown" , "ZSFsysUp","JERsysDown" , "JERsysUp","JESsysDown" , "JESsysUp","PUsysDown" , "PUsysUp","SFsysDown" , "SFsysUp","METsysDown" , "METsysUp","TopMassDown" , "TopMassUp","Q2Down" , "Q2Up","eleSFsysDown" , "eleSFsysUp","matchingDown" , "matchingUp"}; 
 

  if(region == 0){
  TString cutLabel[nPlots] =     { "ptsys_1j1t" };
  
   } else if(region == 1){
   TString cutLabel[nPlots] =     { "ptsys_2j1t"};
  
   }else if(region == 2){
   TString cutLabel[nPlots] =     { "ptsys_2j2t" };
  

  
  }

  TString modeString[3] = {"0","1","2"};
  TString regionString[3] = {"0","1","2"};
  
  TString plotExtension = "sys_plot_"; // name of the plots
  
  
  TString plotAnalysis = "Systematics"; // directory in plots where the plots are saved
  

  
   
  
  
  //Make plots   
  TH1F* histo_tt;
  
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
   TH1F* histo_tt_PDFDown;
  TH1F* histo_tt_PDFUp;
  
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
   TH1F* histo_other_PDFDown;
  TH1F* histo_other_PDFUp;
  
  
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
  TH1F* histo_twdr_TopMassDown; 
  TH1F* histo_twdr_TopMassUp;
  TH1F* histo_twdr_Q2Down; 
  TH1F* histo_twdr_Q2Up;
  TH1F* histo_twdr_eleSFsysDown; 
  TH1F* histo_twdr_eleSFsysUp;
  TH1F* histo_twdr_PDFDown; 
TH1F* histo_twdr_PDFUp;
  TH1F* histo_twds; 
  
  
  
   //Make double  
  double double_tt= 0.0 ; 
  double double_tt_ZSFsysDown= 0.0 ; 
  double double_tt_ZSFsysUp= 0.0 ; 
  double double_tt_JERsysDown= 0.0 ; 
  double double_tt_JERsysUp= 0.0 ; 
  double double_tt_JESsysDown= 0.0 ; 
  double double_tt_JESsysUp= 0.0 ; 
  double double_tt_PUsysDown= 0.0 ; 
  double double_tt_PUsysUp= 0.0 ; 
  double double_tt_SFsysDown= 0.0 ; 
  double double_tt_SFsysUp= 0.0 ; 
  double double_tt_METsysDown= 0.0 ; 
  double double_tt_METsysUp= 0.0 ; 
  double double_tt_TopMassDown= 0.0 ; 
  double double_tt_TopMassUp= 0.0 ; 
  double double_tt_Q2Down= 0.0 ; 
  double double_tt_Q2Up= 0.0 ; 
  double double_tt_eleSFsysDown= 0.0 ; 
  double double_tt_eleSFsysUp= 0.0 ; 
  double double_tt_matchingDown= 0.0 ; 
  double double_tt_matchingUp= 0.0 ; 
  double double_tt_PDFDown= 0.0 ; 
  double double_tt_PDFUp= 0.0 ; 
  
  double double_other= 0.0 ; 
  double double_other_ZSFsysDown= 0.0 ; 
  double double_other_ZSFsysUp= 0.0 ; 
  double double_other_JERsysDown= 0.0 ; 
  double double_other_JERsysUp= 0.0 ; 
  double double_other_JESsysDown= 0.0 ; 
  double double_other_JESsysUp= 0.0 ; 
  double double_other_PUsysDown= 0.0 ; 
  double double_other_PUsysUp= 0.0 ; 
  double double_other_SFsysDown= 0.0 ; 
  double double_other_SFsysUp= 0.0 ; 
  double double_other_METsysDown= 0.0 ; 
  double double_other_METsysUp= 0.0 ; 
  double double_other_eleSFsysDown= 0.0 ; 
  double double_other_eleSFsysUp= 0.0 ; 
  double double_other_PDFDown= 0.0 ; 
  double double_other_PDFUp= 0.0 ; 

  
  
  double double_twdr= 0.0 ; 
  double double_twdr_ZSFsysDown= 0.0 ;  
  double double_twdr_ZSFsysUp= 0.0 ;  
  double double_twdr_JERsysDown= 0.0 ;  
  double double_twdr_JERsysUp= 0.0 ;  
  double double_twdr_JESsysDown= 0.0 ;  
  double double_twdr_JESsysUp= 0.0 ; 
  double double_twdr_PUsysDown= 0.0 ;  
  double double_twdr_PUsysUp= 0.0 ; 
  double double_twdr_SFsysDown= 0.0 ;  
  double double_twdr_SFsysUp= 0.0 ; 
  double double_twdr_METsysDown= 0.0 ;  
  double double_twdr_METsysUp= 0.0 ; 
  double double_twdr_TopMassDown= 0.0 ;  
  double double_twdr_TopMassUp= 0.0 ; 
  double double_twdr_Q2Down= 0.0 ;  
  double double_twdr_Q2Up= 0.0 ; 
  double double_twdr_PDFDown= 0.0 ;  
  double double_twdr_PDFUp= 0.0 ; 
  double double_twdr_eleSFsysDown= 0.0 ;  
  double double_twdr_eleSFsysUp= 0.0 ; 
  double double_twds= 0.0 ;  
  
   for(int iPlots = 0; iPlots< nPlots; iPlots++)
   {
     
     
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
	histo_other_PDFDown= (TH1F*) _filePDFother->Get("histo1a");
	histo_other_PDFUp= (TH1F*) _filePDFother->Get("histo1b");       
	       
	
	
       // Get histos tt bar 
       histo_tt= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[1]);
	      // cout << "Normal: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
	histo_tt_PDFDown= (TH1F*) _filePDFtt->Get("histo1a");
	histo_tt_PDFUp= (TH1F*) _filePDFtt->Get("histo1b");
	  
       
       histo_tt_JERsysDown= (TH1F*) _fileJERsysDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	      // cout << "JERsysDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
       histo_tt_JERsysUp= (TH1F*) _fileJERsysUp->Get(cutLabel[iPlots]+ "_" + processName[1]);
	       //cout << "JERsysUp: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 
       histo_tt_JESsysDown= (TH1F*) _fileJESsysDown->Get(cutLabel[iPlots]+ "_" + processName[1]);
	     //  cout <<"JESsysDown: " << cutLabel[iPlots]+ "_" + processName[1] << endl; 	       
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
		       
	       
	       	       
      
       
  	//twds
	
	histo_twds= (TH1F*) _file0->Get(cutLabel[iPlots]+ "_" + processName[2]);
        //cout << "Normal: " << cutLabel[iPlots]+ "_" + processName[0] << endl; 
	
	
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
       histo_twdr_PDFDown= (TH1F*) _filePDFtwdr->Get("histo1a");
	histo_twdr_PDFUp= (TH1F*) _filePDFtwdr->Get("histo1b");
	




       //Get integrals 
      double_tt = histo_tt->Integral(1,histo_tt->GetNbinsX()); 
      //cout << double_tt << endl; 
      
         double_tt_JERsysDown=histo_tt_JERsysDown ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_JERsysUp=histo_tt_JERsysUp ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_JESsysDown=histo_tt_JESsysDown ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_JESsysUp=histo_tt_JESsysUp ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_PUsysDown=histo_tt_PUsysDown ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_PUsysUp=histo_tt_PUsysUp ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_SFsysDown=histo_tt_SFsysDown ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_SFsysUp=histo_tt_SFsysUp ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_METsysDown=histo_tt_METsysDown ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_METsysUp=histo_tt_METsysUp ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_TopMassDown=histo_tt_TopMassDown ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_TopMassUp=histo_tt_TopMassUp ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_Q2Down=histo_tt_Q2Down ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_Q2Up=histo_tt_Q2Up ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_eleSFsysDown=histo_tt_eleSFsysDown ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_eleSFsysUp=histo_tt_eleSFsysUp ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_matchingDown=histo_tt_matchingDown ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_matchingUp=histo_tt_matchingUp ->Integral(1,histo_tt->GetNbinsX()) ; 
	  double_tt_PDFDown=histo_tt_PDFDown ->Integral(1,histo_tt->GetNbinsX()) ; 
         double_tt_PDFUp=histo_tt_PDFUp ->Integral(1,histo_tt->GetNbinsX()) ; 

         double_other=histo_other ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_ZSFsysDown=histo_other_ZSFsysDown->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_ZSFsysUp=histo_other_ZSFsysUp ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_JERsysDown=histo_other_JERsysDown ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_JERsysUp=histo_other_JERsysUp ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_JESsysDown=histo_other_JESsysDown ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_JESsysUp=histo_other_JESsysUp ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_PUsysDown=histo_other_PUsysDown ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_PUsysUp=histo_other_PUsysUp ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_SFsysDown=histo_other_SFsysDown ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_SFsysUp=histo_other_SFsysUp ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_METsysDown=histo_other_METsysDown ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_METsysUp=histo_other_METsysUp ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_eleSFsysDown=histo_other_eleSFsysDown ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_eleSFsysUp=histo_other_eleSFsysUp ->Integral(1,histo_other->GetNbinsX()) ; 
	  double_other_PDFDown=histo_other_PDFDown ->Integral(1,histo_other->GetNbinsX()) ; 
         double_other_PDFUp=histo_other_PDFUp ->Integral(1,histo_other->GetNbinsX()) ; 



         double_twdr=histo_twdr ->Integral(1,histo_twdr->GetNbinsX()) ;   
         double_twdr_JERsysDown=histo_twdr_JERsysDown ->Integral(1,histo_twdr->GetNbinsX()) ;  
         double_twdr_JERsysUp=histo_twdr_JERsysUp ->Integral(1,histo_twdr->GetNbinsX()) ;  
         double_twdr_JESsysDown=histo_twdr_JESsysDown ->Integral(1,histo_twdr->GetNbinsX()) ;  
         double_twdr_JESsysUp=histo_twdr_JESsysUp ->Integral(1,histo_twdr->GetNbinsX()) ; 
         double_twdr_PUsysDown=histo_twdr_PUsysDown ->Integral(1,histo_twdr->GetNbinsX()) ;  
         double_twdr_PUsysUp=histo_twdr_PUsysUp ->Integral(1,histo_twdr->GetNbinsX()) ; 
         double_twdr_SFsysDown=histo_twdr_SFsysDown ->Integral(1,histo_twdr->GetNbinsX()) ;  
         double_twdr_SFsysUp=histo_twdr_SFsysUp ->Integral(1,histo_twdr->GetNbinsX()) ; 
         double_twdr_METsysDown=histo_twdr_METsysDown ->Integral(1,histo_twdr->GetNbinsX()) ;  
         double_twdr_METsysUp=histo_twdr_METsysUp ->Integral(1,histo_twdr->GetNbinsX()) ; 
         double_twdr_TopMassDown=histo_twdr_TopMassDown ->Integral(1,histo_twdr->GetNbinsX()) ;  
         double_twdr_TopMassUp=histo_twdr_TopMassUp ->Integral(1,histo_twdr->GetNbinsX()) ; 
         double_twdr_Q2Down=histo_twdr_Q2Down ->Integral(1,histo_twdr->GetNbinsX()) ;  
         double_twdr_Q2Up=histo_twdr_Q2Up ->Integral(1,histo_twdr->GetNbinsX()) ; 
         double_twdr_eleSFsysDown=histo_twdr_eleSFsysDown ->Integral(1,histo_twdr->GetNbinsX()) ;  
         double_twdr_eleSFsysUp=histo_twdr_eleSFsysUp ->Integral(1,histo_twdr->GetNbinsX()) ; 
         double_twds=histo_twds ->Integral(1,histo_twdr->GetNbinsX()) ;  
	 double_twdr_PDFDown=histo_twdr_PDFDown ->Integral(1,histo_twdr->GetNbinsX()) ;  
         double_twdr_PDFUp=histo_twdr_PDFUp ->Integral(1,histo_twdr->GetNbinsX()) ;
	 
	 cout << "--------------------------------------------------"<< endl; 
	 cout << " Rate influences for ttbar " << endl; 
	 cout << "--------------------------------------------------"<< endl; 
	 double rate_tt_JERsysUp = (double_tt_JERsysUp - double_tt)/double_tt;
	 cout << "JERsysUp: " << rate_tt_JERsysUp << endl;  
	 double rate_tt_JERsysDown = (double_tt_JERsysDown - double_tt)/double_tt; 
	 cout << "JERsysDown: " << rate_tt_JERsysDown << endl;
	 double rate_tt_JESsysUp = (double_tt_JESsysUp - double_tt)/double_tt;
	 cout << "JESsysUp: " << rate_tt_JESsysUp << endl;  
	 double rate_tt_JESsysDown = (double_tt_JESsysDown - double_tt)/double_tt; 
	 cout << "JESsysDown: " << rate_tt_JESsysDown << endl;  
	 double rate_tt_PUsysUp = (double_tt_PUsysUp - double_tt)/double_tt;
	 cout << "PUsysUp: " << rate_tt_PUsysUp << endl;  
	 double rate_tt_PUsysDown = (double_tt_PUsysDown - double_tt)/double_tt; 
	 cout << "PUsysDown: " << rate_tt_PUsysDown << endl; 
	 double rate_tt_METsysUp = (double_tt_METsysUp - double_tt)/double_tt;
	 cout << "METsysUp: " << rate_tt_METsysUp << endl;  
	 double rate_tt_METsysDown = (double_tt_METsysDown - double_tt)/double_tt; 
	 cout << "METsysDown: " << rate_tt_METsysDown << endl; 
	 double rate_tt_SFsysUp = (double_tt_SFsysUp - double_tt)/double_tt;
	 cout << "SFsysUp: " << rate_tt_SFsysUp << endl;  
	 double rate_tt_SFsysDown = (double_tt_SFsysDown - double_tt)/double_tt; 
	 cout << "SFsysDown: " << rate_tt_SFsysDown << endl; 
	 double rate_tt_eleSFsysUp = (double_tt_eleSFsysUp - double_tt)/double_tt;
	 cout << "eleSFsysUp: " << rate_tt_eleSFsysUp << endl;  
	 double rate_tt_eleSFsysDown = (double_tt_eleSFsysDown - double_tt)/double_tt; 
	 cout << "eleSFsysDown: " << rate_tt_eleSFsysDown << endl; 
	 double rate_tt_Q2Up = (double_tt_Q2Up - double_tt)/double_tt;
	 cout << "Q2Up: " << rate_tt_Q2Up << endl;  
	 double rate_tt_Q2Down = (double_tt_Q2Down - double_tt)/double_tt; 
	 cout << "Q2Down: " << rate_tt_Q2Down << endl; 
	 double rate_tt_matchingUp = (double_tt_matchingUp - double_tt)/double_tt;
	 cout << "matchingUp: " << rate_tt_matchingUp << endl;  
	 double rate_tt_matchingDown = (double_tt_matchingDown - double_tt)/double_tt; 
	 cout << "matchingDown: " << rate_tt_matchingDown << endl; 
        double rate_tt_TopMassUp = (double_tt_TopMassUp - double_tt)/double_tt;
	 cout << "TopMassUp: " << rate_tt_TopMassUp << endl;  
	 double rate_tt_TopMassDown = (double_tt_TopMassDown - double_tt)/double_tt; 
	 cout << "TopMassDown: " << rate_tt_TopMassDown << endl; 
	 double rate_tt_PDFUp = (double_tt_PDFUp - double_tt)/double_tt;
	 cout << "PDFUp: " << rate_tt_PDFUp << endl;  
	 double rate_tt_PDFDown = (double_tt_PDFDown - double_tt)/double_tt; 
	 cout << "PDFDown: " << rate_tt_PDFDown << endl; 
	 
	 cout << "--------------------------------------------------"<< endl; 
	 cout << " Rate influences for twdr " << endl; 
	 cout << "--------------------------------------------------"<< endl; 
	 double rate_twdr_JERsysUp = (double_twdr_JERsysUp - double_twdr)/double_twdr;
	 cout << "JERsysUp: " << rate_twdr_JERsysUp << endl;  
	 double rate_twdr_JERsysDown = (double_twdr_JERsysDown - double_twdr)/double_twdr; 
	 cout << "JERsysDown: " << rate_twdr_JERsysDown << endl;
	 double rate_twdr_JESsysUp = (double_twdr_JESsysUp - double_twdr)/double_twdr;
	 cout << "JESsysUp: " << rate_twdr_JESsysUp << endl;  
	 double rate_twdr_JESsysDown = (double_twdr_JESsysDown - double_twdr)/double_twdr; 
	 cout << "JESsysDown: " << rate_twdr_JESsysDown << endl;  
	 double rate_twdr_PUsysUp = (double_twdr_PUsysUp - double_twdr)/double_twdr;
	 cout << "PUsysUp: " << rate_twdr_PUsysUp << endl;  
	 double rate_twdr_PUsysDown = (double_twdr_PUsysDown - double_twdr)/double_twdr; 
	 cout << "PUsysDown: " << rate_twdr_PUsysDown << endl; 
	 double rate_twdr_METsysUp = (double_twdr_METsysUp - double_twdr)/double_twdr;
	 cout << "METsysUp: " << rate_twdr_METsysUp << endl;  
	 double rate_twdr_METsysDown = (double_twdr_METsysDown - double_twdr)/double_twdr; 
	 cout << "METsysDown: " << rate_twdr_METsysDown << endl; 
	 double rate_twdr_SFsysUp = (double_twdr_SFsysUp - double_twdr)/double_twdr;
	 cout << "SFsysUp: " << rate_twdr_SFsysUp << endl;  
	 double rate_twdr_SFsysDown = (double_twdr_SFsysDown - double_twdr)/double_twdr; 
	 cout << "SFsysDown: " << rate_twdr_SFsysDown << endl; 
	 double rate_twdr_eleSFsysUp = (double_twdr_eleSFsysUp - double_twdr)/double_twdr;
	 cout << "eleSFsysUp: " << rate_twdr_eleSFsysUp << endl;  
	 double rate_twdr_eleSFsysDown = (double_twdr_eleSFsysDown - double_twdr)/double_twdr; 
	 cout << "eleSFsysDown: " << rate_twdr_eleSFsysDown << endl; 
	 double rate_twdr_Q2Up = (double_twdr_Q2Up - double_twdr)/double_twdr;
	 cout << "Q2Up: " << rate_twdr_Q2Up << endl;  
	 double rate_twdr_Q2Down = (double_twdr_Q2Down - double_twdr)/double_twdr; 
	 cout << "Q2Down: " << rate_twdr_Q2Down << endl; 
        double rate_twdr_TopMassUp = (double_twdr_TopMassUp - double_twdr)/double_twdr;
	 cout << "TopMassUp: " << rate_twdr_TopMassUp << endl;  
	 double rate_twdr_TopMassDown = (double_twdr_TopMassDown - double_twdr)/double_twdr; 
	 cout << "TopMassDown: " << rate_twdr_TopMassDown << endl; 
	 double rate_twdr_PDFUp = (double_twdr_PDFUp - double_twdr)/double_twdr;
	 cout << "PDFUp: " << rate_twdr_PDFUp << endl;  
	 double rate_twdr_PDFDown = (double_twdr_PDFDown - double_twdr)/double_twdr; 
	 cout << "PDFDown: " << rate_twdr_PDFDown << endl; 
       
	 cout << "--------------------------------------------------"<< endl; 
	 cout << " Rate influences for other " << endl; 
	 cout << "--------------------------------------------------"<< endl; 
	 double rate_other_JERsysUp = (double_other_JERsysUp - double_other)/double_other;
	 cout << "JERsysUp: " << rate_other_JERsysUp << endl;  
	 double rate_other_JERsysDown = (double_other_JERsysDown - double_other)/double_other; 
	 cout << "JERsysDown: " << rate_other_JERsysDown << endl;
	 double rate_other_JESsysUp = (double_other_JESsysUp - double_other)/double_other;
	 cout << "JESsysUp: " << rate_other_JESsysUp << endl;  
	 double rate_other_JESsysDown = (double_other_JESsysDown - double_other)/double_other; 
	 cout << "JESsysDown: " << rate_other_JESsysDown << endl;  
	 double rate_other_PUsysUp = (double_other_PUsysUp - double_other)/double_other;
	 cout << "PUsysUp: " << rate_other_PUsysUp << endl;  
	 double rate_other_PUsysDown = (double_other_PUsysDown - double_other)/double_other; 
	 cout << "PUsysDown: " << rate_other_PUsysDown << endl; 
	 double rate_other_METsysUp = (double_other_METsysUp - double_other)/double_other;
	 cout << "METsysUp: " << rate_other_METsysUp << endl;  
	 double rate_other_METsysDown = (double_other_METsysDown - double_other)/double_other; 
	 cout << "METsysDown: " << rate_other_METsysDown << endl; 
	 double rate_other_SFsysUp = (double_other_SFsysUp - double_other)/double_other;
	 cout << "SFsysUp: " << rate_other_SFsysUp << endl;  
	 double rate_other_SFsysDown = (double_other_SFsysDown - double_other)/double_other; 
	 cout << "SFsysDown: " << rate_other_SFsysDown << endl; 
	 double rate_other_eleSFsysUp = (double_other_eleSFsysUp - double_other)/double_other;
	 cout << "eleSFsysUp: " << rate_other_eleSFsysUp << endl;  
	 double rate_other_eleSFsysDown = (double_other_eleSFsysDown - double_other)/double_other; 
	 cout << "eleSFsysDown: " << rate_other_eleSFsysDown << endl; 
	 double rate_other_ZSFsysUp = (double_other_ZSFsysUp - double_other)/double_other;
	 cout << "ZSFsysUp: " << rate_other_ZSFsysUp << endl;  
	 double rate_other_ZSFsysDown = (double_other_ZSFsysDown - double_other)/double_other; 
	 cout << "ZSFsysDown: " << rate_other_ZSFsysDown << endl; 
	 double rate_other_PDFUp = (double_other_PDFUp - double_other)/double_other;
	 cout << "PDFUp: " << rate_other_PDFUp << endl;  
	 double rate_other_PDFDown = (double_other_PDFDown - double_other)/double_other; 
	 cout << "PDFDown: " << rate_other_PDFDown << endl; 
	 

    } // end plots loop 
    
 
} // end constructor loop
