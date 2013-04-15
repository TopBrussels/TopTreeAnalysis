

#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "inputs.h"


using namespace std;
void sys_tables(int mode = 0){

  double lumi = luminosity;
  
  if (mode == 0 )        lumi =  11966.617;
  else if ( mode == 1)   lumi = 12067.294;
  else if ( mode == 2)   lumi = 12093.792;
  
  
  char myTexFile[300];
  sprintf(myTexFile,"tables/syst_table_%d_%dpb.tex", mode, lumi);
  ofstream salida(myTexFile); 
  
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
  char myRootFileTopMassDown[300];
  char myRootFileTopMassUp[300];
  char myRootFileQ2Down[300];
  char myRootFileQ2Up[300];
  char myRootFileeleSFsysDown[300];
  char myRootFileeleSFsysUp[300];
  char myRootFilematchingDown[300];
  char myRootFilematchingUp[300];
 
  
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
  
   sprintf(myRootFileTopMassDown,"results/TopMassDown_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileTopMassUp,"results/TopMassUp_%dpb_%d.root", (int)lumi, mode);
  
   sprintf(myRootFileQ2Down,"results/Q2Down_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileQ2Up,"results/Q2Up_%dpb_%d.root", (int)lumi, mode);
  
   sprintf(myRootFileeleSFsysDown,"results/eleSFsysDown_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFileeleSFsysUp,"results/eleSFsysUp_%dpb_%d.root", (int)lumi, mode);
  
   sprintf(myRootFilematchingDown,"results/matchingDown_%dpb_%d.root", (int)lumi, mode);
 
  sprintf(myRootFilematchingUp,"results/matchingUp_%dpb_%d.root", (int)lumi, mode);
  
   
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
  
  
  const int nProcess = 2;
  const int nSys = 20;
  TString processName[nProcess] = { "twdr", "tt"};
  TString processLabel[nProcess] = { "\\textbf{$tW$}", "\\textbf{$t \\bar{t}$}"};
   TString systName[nSys] = { "Normal", "JERsysDown" , "JERsysUp","JESsysDown" , "JESsysUp","PUsysDown" , "PUsysUp","SFsysDown" , "SFsysUp","METsysDown" , "METsysUp","TopMassDown" ,"TopMassUp","Q2Down" , "Q2Up","eleSFsysDown" , "eleSFsysUp","matchingDown" , "matchingUp", "tW DS"}; 

  
  
 // if (mode == 1)        processLabel[4] = "\\textbf{$DoubleMu$}";
 // else if (mode == 2)   processLabel[4] = "\\textbf{$DoubleElectron$}";
  
  TString cutLabel[8] = { "blank0", "blank1", "Lepton Sel.", "Inv. Mass", "$E_{T}^{miss}$", "1 Jet", "b-tagging", "$H_{T}$"};
  
  TH1F*  h [nProcess];
  TH1F*  h_JERsysUp[nProcess];
  TH1F*  h_JERsysDown[nProcess];
  TH1F*  h_JESsysUp[nProcess];
  TH1F*  h_JESsysDown[nProcess];
  TH1F*  h_PUsysUp[nProcess];
  TH1F*  h_PUsysDown[nProcess];
  TH1F*  h_SFsysUp[nProcess];
  TH1F*  h_SFsysDown[nProcess];
  TH1F*  h_METsysUp[nProcess];
  TH1F*  h_METsysDown[nProcess];
  TH1F*  h_eleSFsysUp[nProcess];
  TH1F*  h_eleSFsysDown[nProcess];
  TH1F*  h_TopMassUp[nProcess];
  TH1F*  h_TopMassDown[nProcess];
  TH1F*  h_Q2Up[nProcess];
  TH1F*  h_Q2Down[nProcess];
  TH1F*  h_matchingUp[nProcess];
  TH1F*  h_matchingDown[nProcess];
  TH1F*  h_twds;
  
  
  
  
  for(int i=0; i<nProcess; i++){
     h[i] = (TH1F*) _file0->Get("cuts_"+processName[i]);
     h_JERsysUp[i]=(TH1F*) _fileJERsysUp->Get("cuts_"+processName[i]);
     h_JERsysDown[i]= (TH1F*) _fileJERsysDown->Get("cuts_"+processName[i]);
     h_JESsysUp[i]=(TH1F*) _fileJESsysUp->Get("cuts_"+processName[i]);
     h_JESsysDown[i]= (TH1F*) _fileJESsysDown->Get("cuts_"+processName[i]);
     h_PUsysUp[i]=(TH1F*) _filePUsysUp->Get("cuts_"+processName[i]);
     h_PUsysDown[i]=(TH1F*) _filePUsysDown->Get("cuts_"+processName[i]);
     h_SFsysUp[i]=(TH1F*) _fileSFsysUp->Get("cuts_"+processName[i]);
     h_SFsysDown[i]=(TH1F*) _fileSFsysDown->Get("cuts_"+processName[i]);
     h_METsysUp[i]=(TH1F*) _fileMETsysUp->Get("cuts_"+processName[i]);
     h_METsysDown[i]= (TH1F*) _fileMETsysDown->Get("cuts_"+processName[i]);
     h_eleSFsysUp[i]= (TH1F*) _fileeleSFsysUp->Get("cuts_"+processName[i]);
     h_eleSFsysDown[i] = (TH1F*) _fileeleSFsysDown->Get("cuts_"+processName[i]);
     h_TopMassUp[i]=(TH1F*) _fileTopMassUp->Get("cuts_"+processName[i]);
     h_TopMassDown[i]= (TH1F*) _fileTopMassDown->Get("cuts_"+processName[i]);
     h_Q2Up[i]=(TH1F*) _fileQ2Up->Get("cuts_"+processName[i]);
     h_Q2Down[i]=(TH1F*) _fileQ2Down->Get("cuts_"+processName[i]); 
     h_matchingUp[i]= (TH1F*) _filematchingUp->Get("cuts_"+processName[i]);
     h_matchingDown[i] = (TH1F*) _filematchingDown->Get("cuts_"+processName[i]);
    
  }
  
  h_twds = (TH1F*) _file0->Get("cuts_twds");
  

  
  double vectorValue[nProcess][17][4];
  double vectorValue_JERsysUp[nProcess][17][4];
  double vectorValue_JERsysDown[nProcess][17][4];
  double vectorValue_JESsysUp[nProcess][17][4];
  double vectorValue_JESsysDown[nProcess][17][4];
  double vectorValue_PUsysUp[nProcess][17][4];
  double vectorValue_PUsysDown[nProcess][17][4];
  double vectorValue_SFsysUp[nProcess][17][4];
  double vectorValue_SFsysDown[nProcess][17][4];
  double vectorValue_METsysUp[nProcess][17][4];
  double vectorValue_METsysDown[nProcess][17][4];
  double vectorValue_eleSFsysUp[nProcess][17][4];
  double vectorValue_eleSFsysDown[nProcess][17][4];
  double vectorValue_TopMassUp[nProcess][17][4];
  double vectorValue_TopMassDown[nProcess][17][4];
  double vectorValue_Q2Up[nProcess][17][4];
  double vectorValue_Q2Down[nProcess][17][4];
  double vectorValue_matchingUp[nProcess][17][4];
  double vectorValue_matchingDown[nProcess][17][4];
  double vectorValue_twds[1][17][4];
  
  
  for (int i = 0; i < 16; i++){
    for (int j = 0; j < nProcess; j++){
      vectorValue[j][i][0] = h[j]->GetBinContent(i);
      vectorValue[j][i][1] = precision(h[j]->GetBinError(i));
      vectorValue[j][i][2] = h[j]->GetBinError(i);
    
      vectorValue_JERsysUp[j][i][0] =  h_JERsysUp[j]->GetBinContent(i);
      vectorValue_JERsysUp[j][i][1] = precision(h_JERsysUp[j]->GetBinError(i));
      vectorValue_JERsysUp[j][i][2] = h_JERsysUp[j]->GetBinError(i);

      vectorValue_JERsysDown[j][i][0] = h_JERsysDown[j]->GetBinContent(i);
      vectorValue_JERsysDown[j][i][1] = precision(h_JERsysDown[j]->GetBinError(i));
      vectorValue_JERsysDown[j][i][2] = h_JERsysDown[j]->GetBinError(i);
            
      vectorValue_JESsysUp[j][i][0] = h_JESsysUp[j]->GetBinContent(i);
      vectorValue_JESsysUp[j][i][1] = precision(h_JESsysUp[j]->GetBinError(i));
      vectorValue_JESsysUp[j][i][2] = h_JESsysUp[j]->GetBinError(i);

      vectorValue_JESsysDown[j][i][0] = h_JESsysDown[j]->GetBinContent(i);
      vectorValue_JESsysDown[j][i][1] = precision(h_JESsysDown[j]->GetBinError(i));
      vectorValue_JESsysDown[j][i][2] = h_JESsysDown[j]->GetBinError(i);	    

      
      vectorValue_PUsysUp[j][i][0] = h_PUsysUp[j]->GetBinContent(i);
      vectorValue_PUsysUp[j][i][1] = precision(h_PUsysUp[j]->GetBinError(i));
      vectorValue_PUsysUp[j][i][2] = h_PUsysUp[j]->GetBinError(i);

      vectorValue_PUsysDown[j][i][0] = h_PUsysDown[j]->GetBinContent(i);
      vectorValue_PUsysDown[j][i][1] = precision(h_PUsysDown[j]->GetBinError(i));
      vectorValue_PUsysDown[j][i][2] = h_PUsysDown[j]->GetBinError(i);
      
      vectorValue_SFsysUp[j][i][0] = h_SFsysUp[j]->GetBinContent(i);
      vectorValue_SFsysUp[j][i][1] = precision(h_SFsysUp[j]->GetBinError(i));
      vectorValue_SFsysUp[j][i][2] = h_SFsysUp[j]->GetBinError(i);

      vectorValue_SFsysDown[j][i][0] = h_SFsysDown[j]->GetBinContent(i);
      vectorValue_SFsysDown[j][i][1] = precision(h_SFsysDown[j]->GetBinError(i));
      vectorValue_SFsysDown[j][i][2] = h_SFsysDown[j]->GetBinError(i);
      
      vectorValue_METsysUp[j][i][0] = h_METsysUp[j]->GetBinContent(i);
      vectorValue_METsysUp[j][i][1] = precision(h_METsysUp[j]->GetBinError(i));
      vectorValue_METsysUp[j][i][2] = h_METsysUp[j]->GetBinError(i);

      vectorValue_METsysDown[j][i][0] = h_METsysDown[j]->GetBinContent(i);
      vectorValue_METsysDown[j][i][1] = precision(h_METsysDown[j]->GetBinError(i));
      vectorValue_METsysDown[j][i][2] = h_METsysDown[j]->GetBinError(i);
      
      vectorValue_eleSFsysUp[j][i][0] = h_eleSFsysUp[j]->GetBinContent(i);
      vectorValue_eleSFsysUp[j][i][1] = precision(h_eleSFsysUp[j]->GetBinError(i));
      vectorValue_eleSFsysUp[j][i][2] = h_eleSFsysUp[j]->GetBinError(i);

      vectorValue_eleSFsysDown[j][i][0] = h_eleSFsysDown[j]->GetBinContent(i);
      vectorValue_eleSFsysDown[j][i][1] = precision(h_eleSFsysDown[j]->GetBinError(i));
      vectorValue_eleSFsysDown[j][i][2] = h_eleSFsysDown[j]->GetBinError(i);
      
      vectorValue_TopMassUp[j][i][0] = h_TopMassUp[j]->GetBinContent(i);
      vectorValue_TopMassUp[j][i][1] = precision(h_TopMassUp[j]->GetBinError(i));
      vectorValue_TopMassUp[j][i][2] = h_TopMassUp[j]->GetBinError(i);

      vectorValue_TopMassDown[j][i][0] = h_TopMassDown[j]->GetBinContent(i);
      vectorValue_TopMassDown[j][i][1] = precision(h_TopMassDown[j]->GetBinError(i));
      vectorValue_TopMassDown[j][i][2] = h_TopMassDown[j]->GetBinError(i);
      
      vectorValue_Q2Up[j][i][0] = h_Q2Up[j]->GetBinContent(i);
      vectorValue_Q2Up[j][i][1] = precision(h_Q2Up[j]->GetBinError(i));
      vectorValue_Q2Up[j][i][2] = h_Q2Up[j]->GetBinError(i);

      vectorValue_Q2Down[j][i][0] = h_Q2Down[j]->GetBinContent(i);
      vectorValue_Q2Down[j][i][1] = precision(h_Q2Down[j]->GetBinError(i));
      vectorValue_Q2Down[j][i][2] = h_Q2Down[j]->GetBinError(i);
      
      vectorValue_matchingUp[j][i][0] = h_matchingUp[j]->GetBinContent(i);
      vectorValue_matchingUp[j][i][1] = precision(h_matchingUp[j]->GetBinError(i));
      vectorValue_matchingUp[j][i][2] = h_matchingUp[j]->GetBinError(i);

      vectorValue_matchingDown[j][i][0] = h_matchingDown[j]->GetBinContent(i);
      vectorValue_matchingDown[j][i][1] = precision(h_matchingDown[j]->GetBinError(i));
      vectorValue_matchingDown[j][i][2] = h_matchingDown[j]->GetBinError(i);
 
    }  
    
      vectorValue_twds[0][i][0] = h_twds->GetBinContent(i);
      vectorValue_twds[0][i][1] = precision(h_twds->GetBinError(i));
      vectorValue_twds[0][i][2] = h_twds->GetBinError(i);
  }
  
  salida << "\\documentclass[a4paper,12pt]{article}" << endl;
  salida << "\\begin{document}" << endl;
  salida << endl;
  salida << endl;
  
       
  /////////////////////
  /// JER      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 0; l < 1; l++){

   for (int k = 0; k < 3; k++){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 1; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_JERsysDown[j][i][1]) << vectorValue_JERsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_JERsysDown[j][i][1])<< vectorValue_JERsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_JERsysUp[j][i][1]) << vectorValue_JERsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_JERsysUp[j][i][1])<< vectorValue_JERsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;

  /////////////////////
  /// JES      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 0; l < 1; l++){

   for (int k = 0; k < 5; k++){
      if(k == 0|| k == 3 || k == 4){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 1; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_JESsysDown[j][i][1]) << vectorValue_JESsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_JESsysDown[j][i][1])<< vectorValue_JESsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_JESsysUp[j][i][1]) << vectorValue_JESsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_JESsysUp[j][i][1])<< vectorValue_JESsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
  
    /////////////////////
  /// PU      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 0; l < 1; l++){

   for (int k = 0; k < 7; k++){
       if(k == 0|| k == 5 || k == 6){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 1; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_PUsysDown[j][i][1]) << vectorValue_PUsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_PUsysDown[j][i][1])<< vectorValue_PUsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_PUsysUp[j][i][1]) << vectorValue_PUsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_PUsysUp[j][i][1])<< vectorValue_PUsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
    /////////////////////
  /// SF      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 0; l < 1; l++){

   for (int k = 0; k < 9; k++){
       if(k == 0|| k == 7 || k == 8){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 1; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_SFsysDown[j][i][1]) << vectorValue_SFsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_SFsysDown[j][i][1])<< vectorValue_SFsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_SFsysUp[j][i][1]) << vectorValue_SFsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_SFsysUp[j][i][1])<< vectorValue_SFsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
  
  
    /////////////////////
  /// MET      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 0; l < 1; l++){

   for (int k = 0; k < 11; k++){
      if(k == 0|| k == 9 || k == 10){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	
    }
  }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 1; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_METsysDown[j][i][1]) << vectorValue_METsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_METsysDown[j][i][1])<< vectorValue_METsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_METsysUp[j][i][1]) << vectorValue_METsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_METsysUp[j][i][1])<< vectorValue_METsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
  
    /////////////////////
  /// eleSF      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 0; l < 1; l++){

   for (int k = 0; k < 17; k++){
       if(k == 0|| k == 15 || k == 16){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 1; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_eleSFsysDown[j][i][1]) << vectorValue_eleSFsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_eleSFsysDown[j][i][1])<< vectorValue_eleSFsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_eleSFsysUp[j][i][1]) << vectorValue_eleSFsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_eleSFsysUp[j][i][1])<< vectorValue_eleSFsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
  
  
    /////////////////////
  /// TopMass      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 0; l < 1; l++){

   for (int k = 0; k < 13; k++){
       if(k == 0|| k == 11 || k == 12){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 1; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_TopMassDown[j][i][1]) << vectorValue_TopMassDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_TopMassDown[j][i][1])<< vectorValue_TopMassDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_TopMassUp[j][i][1]) << vectorValue_TopMassUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_TopMassUp[j][i][1])<< vectorValue_TopMassUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
    
    /////////////////////
  /// Q2      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 0; l < 1; l++){

   for (int k = 0; k < 15; k++){
      if(k == 0|| k == 13 || k == 14){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 1; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_Q2Down[j][i][1]) << vectorValue_Q2Down[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_Q2Down[j][i][1])<< vectorValue_Q2Down[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_Q2Up[j][i][1]) << vectorValue_Q2Up[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_Q2Up[j][i][1])<< vectorValue_Q2Up[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;

  
    /////////////////////
  /// matching      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 0; l < 1; l++){

   for (int k = 0; k < 19; k++){
     if(k == 0|| k ==  17|| k == 18){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 1; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_matchingDown[j][i][1]) << vectorValue_matchingDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_matchingDown[j][i][1])<< vectorValue_matchingDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_matchingUp[j][i][1]) << vectorValue_matchingUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_matchingUp[j][i][1])<< vectorValue_matchingUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
  
    salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
 
     /////////////////////
  /// twds     ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 0; l < 1; l++){

   for (int k = 0; k < 20; k++){
     if(k == 0|| k ==  19){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 0; j < 1; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_twds[0][i][1]) << vectorValue_twds[0][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_twds[0][i][1])<< vectorValue_twds[0][i][2];
	
	
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
 


///// FOR TTBAR /////////

  /////////////////////
  /// JER      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 1; l < nProcess; l++){

   for (int k = 0; k < 3; k++){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 1; j < nProcess; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_JERsysDown[j][i][1]) << vectorValue_JERsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_JERsysDown[j][i][1])<< vectorValue_JERsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_JERsysUp[j][i][1]) << vectorValue_JERsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_JERsysUp[j][i][1])<< vectorValue_JERsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;

  /////////////////////
  /// JES      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 1; l < nProcess; l++){

   for (int k = 0; k < 5; k++){
      if(k == 0|| k == 3 || k == 4){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 1; j < nProcess; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_JESsysDown[j][i][1]) << vectorValue_JESsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_JESsysDown[j][i][1])<< vectorValue_JESsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_JESsysUp[j][i][1]) << vectorValue_JESsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_JESsysUp[j][i][1])<< vectorValue_JESsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
  
    /////////////////////
  /// PU      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 1; l < nProcess; l++){

   for (int k = 0; k < 7; k++){
       if(k == 0|| k == 5 || k == 6){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 1; j < nProcess; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_PUsysDown[j][i][1]) << vectorValue_PUsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_PUsysDown[j][i][1])<< vectorValue_PUsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_PUsysUp[j][i][1]) << vectorValue_PUsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_PUsysUp[j][i][1])<< vectorValue_PUsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
    /////////////////////
  /// SF      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 1; l < nProcess; l++){

   for (int k = 0; k < 9; k++){
       if(k == 0|| k == 7 || k == 8){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 1; j < nProcess; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_SFsysDown[j][i][1]) << vectorValue_SFsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_SFsysDown[j][i][1])<< vectorValue_SFsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_SFsysUp[j][i][1]) << vectorValue_SFsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_SFsysUp[j][i][1])<< vectorValue_SFsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
  
  
    /////////////////////
  /// MET      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 1; l < nProcess; l++){

   for (int k = 0; k < 11; k++){
      if(k == 0|| k == 9 || k == 10){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	
    }
  }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 1; j < nProcess; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_METsysDown[j][i][1]) << vectorValue_METsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_METsysDown[j][i][1])<< vectorValue_METsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_METsysUp[j][i][1]) << vectorValue_METsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_METsysUp[j][i][1])<< vectorValue_METsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
  
    /////////////////////
  /// eleSF      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 1; l < nProcess; l++){

   for (int k = 0; k < 17; k++){
       if(k == 0|| k == 15 || k == 16){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 1; j < nProcess; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_eleSFsysDown[j][i][1]) << vectorValue_eleSFsysDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_eleSFsysDown[j][i][1])<< vectorValue_eleSFsysDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_eleSFsysUp[j][i][1]) << vectorValue_eleSFsysUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_eleSFsysUp[j][i][1])<< vectorValue_eleSFsysUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
  
  
    /////////////////////
  /// TopMass      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 1; l < nProcess; l++){

   for (int k = 0; k < 13; k++){
       if(k == 0|| k == 11 || k == 12){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 1; j < nProcess; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_TopMassDown[j][i][1]) << vectorValue_TopMassDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_TopMassDown[j][i][1])<< vectorValue_TopMassDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_TopMassUp[j][i][1]) << vectorValue_TopMassUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_TopMassUp[j][i][1])<< vectorValue_TopMassUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
    
    /////////////////////
  /// Q2      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 1; l < nProcess; l++){

   for (int k = 0; k < 15; k++){
      if(k == 0|| k == 13 || k == 14){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 1; j < nProcess; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_Q2Down[j][i][1]) << vectorValue_Q2Down[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_Q2Down[j][i][1])<< vectorValue_Q2Down[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_Q2Up[j][i][1]) << vectorValue_Q2Up[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_Q2Up[j][i][1])<< vectorValue_Q2Up[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;

  
    /////////////////////
  /// matching      ///////
  ////////////////////
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  
  for (int l = 1; l < nProcess; l++){

   for (int k = 0; k < 19; k++){
     if(k == 0|| k ==  17|| k == 18){
         //   cout << " & " <<  processLabel[l] << " " << systName[k]; 
            salida << " & " << processLabel[l] << " " << systName[k];
	}
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
 
  for (int i=2; i < 8; i++){
    salida << cutLabel[i];
    for (int j = 1; j < nProcess; j++){ // only for twdr
      
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_matchingDown[j][i][1]) << vectorValue_matchingDown[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_matchingDown[j][i][1])<< vectorValue_matchingDown[j][i][2];
	
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue_matchingUp[j][i][1]) << vectorValue_matchingUp[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue_matchingUp[j][i][1])<< vectorValue_matchingUp[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
 
 
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  






  
  salida << "\\end{document}" << endl;
 
}


double normalization(double nevents, double xsec, double lumi){
  
  if (nevents !=0) return lumi*xsec/nevents;
  else return 1;
  
}

double precision(double error){
  
  int precisionValue;
  double factErr = 0; 
  int iN = 0;
  if (error == 0 || error >= 1) precisionValue = 0;
  else if (error < 1) {
    iN = 0;
    factErr = 0; 
    while (factErr < 1){
      factErr = error*(10**iN);
      iN++;  
    }
    precisionValue = iN-1;
  }
  
  if (factErr > 9.5) precisionValue-=1;
  
  return precisionValue;
  
}

double errorEfi (double efi, double errTotal, double totalEvents, double error, double number){
  
  if (totalEvents !=0 && number != 0) return efi*((errTotal/totalEvents)+(error/number));
  else return 0;
  
}

double efficiency (double finalevents, double totalevents){
  
  if (totalevents !=0) return finalevents*100/totalevents;
  else return 0;
  
}
