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
#include "inputs.h"

using namespace std;
void datacardmaker_v02(){

  ofstream datacard("singletop_tW_5fb_dycr.txt"); 
 
  TString processName[3] =  { "twdr", "tt","others_2"};
  char myRootFile[300]; 
  TH1F*  hdata [3];
  TH1F*  hnominal [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    sprintf(myRootFile,"outputs/out_%d_data.root", i);
    TFile *_file0 = TFile::Open(myRootFile);
    hdata[mode] = (TH1F*) _file0->Get("R");
    for (int j = 0; j < 3; j++){
      sprintf(myRootFile,"outputs/out_%d_", i);
      TFile *_file1 = TFile::Open(myRootFile + processName[j] + ".root");
      hnominal[mode][j] = (TH1F*) _file1->Get("R");
    
    }
  } 
  
  TH1F*  hdatady [3];
  TH1F*  hnominaldy [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    sprintf(myRootFile,"outputs/out_%d_data.root", i);
    TFile *_file2 = TFile::Open(myRootFile);
    hdatady[mode] = (TH1F*) _file2->Get("R_dy");
    for (int j = 0; j < 3; j++){
      sprintf(myRootFile,"outputs/out_%d_", i);
      TFile *_file3 = TFile::Open(myRootFile + processName[j] + ".root");
      hnominaldy[mode][j] = (TH1F*) _file3->Get("R_dy");
    
    }
  } 
  
  
  datacard << setprecision(4) << "# this is the version with *exclusive* jet / tag bins " << endl;
  datacard << setprecision(4) << "# based on 4.9/fb " << endl;
  datacard << setprecision(4) << "imax 11 # number of bins " << endl;
  datacard << setprecision(4) << "jmax 2 # number of processes - 1 " << endl;
  datacard << setprecision(4) << "kmax * # number of uncertainties " << endl;
  datacard << setprecision(4) << "------------ " << endl;
  datacard << setprecision(4) << "bin \t\tee1j1t\t emu1j1t\t mumu1j1t\t ee2j1t\t emu2j1t\t mumu2j1t\t ee2j2t\t emu2j2t\t mumu2j2t\t  eein\t  mumuin " << endl;
  datacard << setprecision(4) << "observation" ;
  
  for (int i = 0; i < 3; i++){
    datacard << setprecision(4) << "\t " << hdata[i]->GetBinContent(2) ;
  }
  for (int i = 0; i < 3; i++){
    datacard << setprecision(4) << "\t " << hdata[i]->GetBinContent(7) ;
  }
  for (int i = 0; i < 3; i++){
    datacard << setprecision(4) << "\t " << hdata[i]->GetBinContent(8) ;
  }
  for (int i = 0; i < 3; i++){
    if (i!=1) datacard << setprecision(4) << "\t " << hdatady[i]->GetBinContent(3) ;
   
  }
  
  datacard << setprecision(4) << endl;
  datacard << setprecision(4) << "------------ " << endl;
  datacard << setprecision(4) << "bin        \tee1j1t\t ee1j1t\t ee1j1t\t emu1j1t\t emu1j1t\t emu1j1t\t mumu1j1t\t mumu1j1t\t mumu1j1t\t ee2j1t\t ee2j1t \t ee2j1t\t emu2j1t\t emu2j1t\t emu2j1t\t mumu2j1t\t mumu2j1t\t mumu2j1t\t ee2j2t\t ee2j2t\t ee2j2t\t emu2j2t\t emu2j2t\t emu2j2t\t mumu2j2t\t mumu2j2t\t mumu2j2t\t eein\t eein\t eein\t mumuin\t mumuin\t mumuin" << endl;
  datacard << setprecision(4) << "process    \tst\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t " << endl;
  datacard << setprecision(4) << "process    \t0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t " << endl;
  datacard << setprecision(4) << "rate       \t" ;
  
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(4) << hnominal[i][j]->GetBinContent(2) << "\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(4) << hnominal[i][j]->GetBinContent(7) << "\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(4) << hnominal[i][j]->GetBinContent(8) << "\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if (i != 1) datacard << setprecision(4) << hnominaldy[i][j]->GetBinContent(3) << "\t ";
 
    }
  }
  datacard << setprecision(4) << endl;

  datacard << setprecision(4) << "---------------------- " << endl;
  datacard << setprecision(4) <<  "lumi      lnN\t";
  for (int i = 0; i < 4; i++){
    for (int j = 0; j< 3; j++){
      for (int k = 0; k < 3; k++){
        if (i < 3 || j != 1)  datacard << setprecision(4) << "1.022\t" ;
	else if (j! = 1 && i !=3)  datacard << setprecision(4) << "-\t" ;
      }
    }
  }
  datacard << endl;
  datacard << setprecision(4) <<  "hlt       lnN\t";
  for (int i = 0; i < 4; i++){
    for (int j = 0; j< 3; j++){
      for (int k = 0; k < 3; k++){
        if (i < 3 || j != 1)  datacard << setprecision(4) << "1.015\t" ;
	else if (j! = 1 || i !=3)  datacard << setprecision(4) << "-\t" ;
      }
    }
  }
  datacard << endl;
  
  datacard << setprecision(4) <<  "ele       lnN\t";
  for (int i = 0; i < 4; i++){
    for (int j = 0; j< 3; j++){
      for (int k = 0; k < 3; k++){
        if ((i < 3 && j != 2) || (i == 3 && j == 0))  datacard << setprecision(4) << "1.02\t" ;
	else if (j! = 1 || i !=3)  datacard << setprecision(4) << "-\t" ;
      }
    }
  }
  datacard << endl;
  

  datacard << setprecision(4) <<  "mu        lnN\t";
  for (int i = 0; i < 4; i++){
    for (int j = 0; j< 3; j++){
      for (int k = 0; k < 3; k++){
        if ((i < 3 && j != 0) || (i == 3 && j == 2))  datacard << setprecision(4) << "1.01\t" ;
	else if (j != 1 || i !=3)  datacard << setprecision(4) << "-\t" ;
      }
    }
  }
  datacard << endl;
  
  

  
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 3; j++){
      if (j == 0){
        sprintf(myRootFile,"outputs/out_%d_tw_sup.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hup[mode][j] = (TH1F*) _file1->Get("R");
	sprintf(myRootFile,"outputs/out_%d_tw_sdo.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hdown[mode][j] = (TH1F*) _file1->Get("R");
      } else if (j == 1) {
	sprintf(myRootFile,"outputs/out_%d_tt_scaleup.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hup[mode][j] = (TH1F*) _file1->Get("R");
	sprintf(myRootFile,"outputs/out_%d_tt_scaledown.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hdown[mode][j] = (TH1F*) _file1->Get("R");
      }
    }
  } 
  

  
   TH1F*  hupdy [3][3];
  TH1F*  hdowndy [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
       for (int j = 0; j < 3; j++){
      if (j == 0){
        sprintf(myRootFile,"outputs/out_%d_tw_sup.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hupdy[mode][j] = (TH1F*) _file1->Get("R_dy");
	sprintf(myRootFile,"outputs/out_%d_tw_sdo.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hdowndy[mode][j] = (TH1F*) _file1->Get("R_dy");
      } else if (j == 1) {
	sprintf(myRootFile,"outputs/out_%d_tt_scaleup.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hupdy[mode][j] = (TH1F*) _file1->Get("R_dy");
	sprintf(myRootFile,"outputs/out_%d_tt_scaledown.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hdowndy[mode][j] = (TH1F*) _file1->Get("R_dy");
      }
    }
  } 
  
  datacard << setprecision(4) << "ttscale   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(2) + hdown[i][j]->GetBinContent(2))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hup[i][j]->GetBinContent(2) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      } 
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(7) + hdown[i][j]->GetBinContent(7))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hup[i][j]->GetBinContent(7) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(8) + hdown[i][j]->GetBinContent(8))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hup[i][j]->GetBinContent(8) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i!=1 && j == 1){
	double average = (hupdy[i][j]->GetBinContent(3) + hdowndy[i][j]->GetBinContent(3))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hupdy[i][j]->GetBinContent(3) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else if (i! = 1) datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(4) << endl;
  
  datacard << setprecision(4) << "twscale   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = (hup[i][j]->GetBinContent(2) + hdown[i][j]->GetBinContent(2))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hup[i][j]->GetBinContent(2) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      } 
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = (hup[i][j]->GetBinContent(7) + hdown[i][j]->GetBinContent(7))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hup[i][j]->GetBinContent(7) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = (hup[i][j]->GetBinContent(8) + hdown[i][j]->GetBinContent(8))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hup[i][j]->GetBinContent(8) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i!=1 && j == 0){
	double average = (hupdy[i][j]->GetBinContent(3) + hdowndy[i][j]->GetBinContent(3))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hupdy[i][j]->GetBinContent(3) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else if (i! = 1) datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(4) << endl;
  
  
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 3; j++){
      if (j == 1) {
	sprintf(myRootFile,"outputs/out_%d_tt_matchingup.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hup[mode][j] = (TH1F*) _file1->Get("R");
	sprintf(myRootFile,"outputs/out_%d_tt_matchingdown.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hdown[mode][j] = (TH1F*) _file1->Get("R");
      }
    }
  } 
  TH1F*  hupdy [3][3];
  TH1F*  hdowndy [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 3; j++){
      if (j == 1) {
	sprintf(myRootFile,"outputs/out_%d_tt_matchingup.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hupdy[mode][j] = (TH1F*) _file1->Get("R_dy");
	sprintf(myRootFile,"outputs/out_%d_tt_matchingdown.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hdowndy[mode][j] = (TH1F*) _file1->Get("R_dy");
      }
    }
  } 
  
  
  
  datacard << setprecision(4) << "ttmatch   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(2) + hdown[i][j]->GetBinContent(2))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hup[i][j]->GetBinContent(2) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      } 
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(7) + hdown[i][j]->GetBinContent(7))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hup[i][j]->GetBinContent(7) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(8) + hdown[i][j]->GetBinContent(8))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hup[i][j]->GetBinContent(8) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i!=1 && j == 1){
	double average = (hupdy[i][j]->GetBinContent(8) + hdowndy[i][j]->GetBinContent(3))/2;
	if (average !=0) datacard << setprecision(4) << 1 + fabs((hupdy[i][j]->GetBinContent(3) - average)/average) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else if (i! = 1) datacard << setprecision(4) << "-\t ";
    }
  }
  
  datacard << setprecision(4) << endl;
  
  TH1F*  h [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 3; j++){
      if (j == 0) {
	sprintf(myRootFile,"outputs/out_%d_twds.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	h[mode][j] = (TH1F*) _file1->Get("R");
      }
    }
  } 
  
  TH1F*  hdy [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 3; j++){
      if (j == 0) {
	sprintf(myRootFile,"outputs/out_%d_twds.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hdy[mode][j] = (TH1F*) _file1->Get("R_dy");
      }
    }
  } 
  
  datacard << setprecision(4) << "twdrds    lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(4) << 1 + fabs((h[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(4) << 1 + fabs((h[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(4) << 1 + fabs((h[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0 && i != 1){
	if (hnominaldy[i][j]->GetBinContent(3) !=0) datacard << setprecision(4) << 1 + fabs((hdy[i][j]->GetBinContent(3) - hnominaldy[i][j]->GetBinContent(3))/hnominaldy[i][j]->GetBinContent(3)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else if (i! = 1) datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(4) << endl;
  
  TString SystName = "PU";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 3; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hup[mode][j] = (TH1F*) _file1->Get("R");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdown[mode][j] = (TH1F*) _file1->Get("R"); 
    }
  } 
  
    TH1F*  hupdy [3][3];
  TH1F*  hdowndy [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 3; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hupdy[mode][j] = (TH1F*) _file1->Get("R_dy");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdowndy[mode][j] = (TH1F*) _file1->Get("R_dy"); 
    }
  } 
  
  datacard << setprecision(4) << "pu        lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
   //   if(j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(4) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)), fabs((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) ) << "\t ";
	else datacard << setprecision(4) << "-\t ";
     // }
     // else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
   //   if(j == 0){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(4) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)), fabs((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) ) << "\t ";
	else datacard << setprecision(4) << "-\t ";
   //   }
     // else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
   //   if(j == 0){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(4) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)), fabs((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) ) << "\t ";
	else datacard << setprecision(4) << "-\t ";
    // }
     // else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
     if(i != 1){
	if (hnominaldy[i][j]->GetBinContent(3) !=0) datacard << setprecision(4) << 1 + TMath::Max(fabs((hupdy[i][j]->GetBinContent(3) - hnominaldy[i][j]->GetBinContent(3))/hnominaldy[i][j]->GetBinContent(3)), fabs((hdowndy[i][j]->GetBinContent(3) - hnominaldy[i][j]->GetBinContent(3))/hnominaldy[i][j]->GetBinContent(3)) ) << "\t ";
	else datacard << setprecision(4) << "-\t ";
     }
     
    }
  }
  
  datacard << setprecision(4) << endl;
  
  datacard << setprecision(4) <<  "ttxs      lnN\t";
  for (int i = 0; i < 4; i++){
    for (int j = 0; j< 3; j++){
      for (int k = 0; k < 3; k++){
        if ((i < 3 || j != 1) && k == 1)  datacard << setprecision(4) << "1.06\t" ;
	else if (i < 3 || j != 1)  datacard << setprecision(4) << "-\t" ;
      }
    }
  }
  datacard << endl;
   datacard << setprecision(4) <<  "dyxs      lnN\t";
  for (int i = 0; i < 4; i++){
    for (int j = 0; j< 3; j++){
      for (int k = 0; k < 3; k++){
        if ((i < 3 || j != 1) && k == 3)  datacard << setprecision(4) << "1.05\t" ;
	else if (i < 3 || j != 1)  datacard << setprecision(4) << "-\t" ;
      }
    }
  }
  datacard << endl;

  
  TString SystName = "JES";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 2; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hup[mode][j] = (TH1F*) _file1->Get("R");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdown[mode][j] = (TH1F*) _file1->Get("R"); 
    }
  } 
  
  TH1F*  hupdy [3][3];
  TH1F*  hdowndy [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 2; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hupdy[mode][j] = (TH1F*) _file1->Get("R_dy");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdowndy[mode][j] = (TH1F*) _file1->Get("R_dy"); 
    }
  } 
  
  datacard << setprecision(4) << "jes       lnN\t";
 for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
      	if (hnominal[i][j]->GetBinContent(2) !=0){
	  double upv = (hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
          double downv =  (hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
          double maxv = TMath::Max(fabs(upv), fabs(downv));
          if (upv*downv >= 0)  datacard << setprecision(4) << 1 + maxv << "\t ";
	  else if (upv < downv)  datacard << setprecision(4) << 1 + downv << "/" << 1+ upv<< "\t ";
	  else  datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	} else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
      	if (hnominal[i][j]->GetBinContent(7) !=0){
	  double upv = (hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7);
          double downv =  (hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7);
          double maxv = TMath::Max(fabs(upv), fabs(downv));
          if (upv*downv >= 0)  datacard << setprecision(4) << 1 + maxv << "\t ";
	  else if (upv < downv)  datacard << setprecision(4) << 1 + downv << "/" << 1+ upv<< "\t ";
	  else  datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	} else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
      	if (hnominal[i][j]->GetBinContent(8) !=0){
	  double upv = (hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8);
          double downv =  (hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8);
          double maxv = TMath::Max(fabs(upv), fabs(downv));
          if (upv*downv >= 0)  datacard << setprecision(4) << 1 + maxv << "\t ";
	  else if (upv < downv)  datacard << setprecision(4) << 1 + downv << "/" << 1+ upv<< "\t ";
	  else  datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	} else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  
    for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3 && i != 1){
      	if (hnominaldy[i][j]->GetBinContent(8) !=0){
	  double upv = (hupdy[i][j]->GetBinContent(8) - hnominaldy[i][j]->GetBinContent(8))/hnominaldy[i][j]->GetBinContent(8);
          double downv =  (hdowndy[i][j]->GetBinContent(8) - hnominaldy[i][j]->GetBinContent(8))/hnominaldy[i][j]->GetBinContent(8);
          double maxv = TMath::Max(fabs(upv), fabs(downv));
          if (upv*downv >= 0)  datacard << setprecision(4) << 1 + maxv << "\t ";
	  else if (upv < downv)  datacard << setprecision(4) << 1 + downv << "/" << 1+ upv<< "\t ";
	  else  datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	} else datacard << setprecision(4) << "-\t ";
      }
    
    }
  }
  
  datacard << setprecision(4) << endl;
  
  TString SystName = "SF";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 2; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hup[mode][j] = (TH1F*) _file1->Get("R");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hdown[mode][j] = (TH1F*) _file1->Get("R"); 
    }
  } 
  
  TH1F*  hupdy [3][3];
  TH1F*  hdowndy [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 2; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hupdy[mode][j] = (TH1F*) _file1->Get("R_dy");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hdowndy[mode][j] = (TH1F*) _file1->Get("R_dy"); 
    }
  } 
  
  datacard << setprecision(4) << "btag      lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
    //  if(j < 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(4) << 1 + ((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) << "/" << 1+ ((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
//      }
  //    else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
     // if(j < 2){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(4) << 1 + ((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) << "/" << 1+ ((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
//      }
  //    else datacard << setprecision(4) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
   //   if(j < 2){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(4) << 1 + ((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) << "/" << 1+ ((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
    //  }
     // else datacard << setprecision(4) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
     if(i != 1){
	if (hnominaldy[i][j]->GetBinContent(3) !=0) datacard << setprecision(4) << 1 + ((hupdy[i][j]->GetBinContent(3) - hnominaldy[i][j]->GetBinContent(3))/hnominaldy[i][j]->GetBinContent(3)) << "/" << 1+ ((hdowndy[i][j]->GetBinContent(3) - hnominaldy[i][j]->GetBinContent(3))/hnominaldy[i][j]->GetBinContent(3)) << "\t ";
	elsedatacard << setprecision(4) << "-\t ";
     }
      
    }
  }
  datacard << setprecision(4) << endl;  
  
  TString SystName = "JER";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 2; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hup[mode][j] = (TH1F*) _file1->Get("R");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdown[mode][j] = (TH1F*) _file1->Get("R"); 
    }
  } 
  TH1F*  hupdy [3][3];
  TH1F*  hdowndy [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 2; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hupdy[mode][j] = (TH1F*) _file1->Get("R_dy");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdowndy[mode][j] = (TH1F*) _file1->Get("R_dy"); 
    }
  } 
  
  datacard << setprecision(4) << "jer       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
     // if(j < 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(5) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)), fabs((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) ) << "\t ";
	else datacard << setprecision(5) << "-\t ";
    //  }
    //  else datacard << setprecision(5) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
 //     if(j < 2){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(5) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)), fabs((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) ) << "\t ";
	else datacard << setprecision(5) << "-\t ";
  //    }
   //   else datacard << setprecision(5) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
  //    if(j < 2){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(5) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)), fabs((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) ) << "\t ";
	else datacard << setprecision(5) << "-\t ";
    //  }
     // else datacard << setprecision(5) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
     if(i != 1){
	if (hnominaldy[i][j]->GetBinContent(3) !=0) datacard << setprecision(5) << 1 + TMath::Max(fabs((hupdy[i][j]->GetBinContent(3) - hnominaldy[i][j]->GetBinContent(3))/hnominaldy[i][j]->GetBinContent(3)), fabs((hdowndy[i][j]->GetBinContent(3) - hnominaldy[i][j]->GetBinContent(3))/hnominaldy[i][j]->GetBinContent(3))) << "\t ";
	else datacard << setprecision(5) << "-\t ";
     }
     
    }
  }
  datacard << setprecision(5) << endl;
  
  
  TString SystName = "MET";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 2; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hup[mode][j] = (TH1F*) _file1->Get("R");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdown[mode][j] = (TH1F*) _file1->Get("R"); 
    }
  } 
  
  
  TH1F*  hupdy [3][3];
  TH1F*  hdowndy [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 2; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hupdy[mode][j] = (TH1F*) _file1->Get("R_dy");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdowndy[mode][j] = (TH1F*) _file1->Get("R_dy"); 
    }
  } 
  
  
  datacard << setprecision(4) << "met       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 

	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(4) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)), fabs((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) ) << "\t ";
	else datacard << setprecision(4) << "-\t ";


    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
   
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(4) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)), fabs((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) ) << "\t ";
	else datacard << setprecision(4) << "-\t ";
    
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
   
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(4) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)), fabs((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) ) << "\t ";
	else datacard << setprecision(4) << "-\t ";
    
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
     if(i != 1 ){
	if (hnominaldy[i][j]->GetBinContent(3) !=0) datacard << setprecision(4) << 1 + ((hupdy[i][j]->GetBinContent(3) - hnominaldy[i][j]->GetBinContent(3))/hnominaldy[i][j]->GetBinContent(3)) << "/" << 1+ ((hdowndy[i][j]->GetBinContent(3) - hnominaldy[i][j]->GetBinContent(3))/hnominaldy[i][j]->GetBinContent(3)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
     }
     
    }
  }
  datacard << setprecision(4) << endl;
  
  datacard << setprecision(4) << "pdf       lnN\t1.020/0.978     1.024/0.975     -     1.018/0.980      1.024/0.975     -     1.018/0.98     1.024/0.975     -    1.025/0.974    1.024/0.975    -     1.021/0.977     1.023/0.976    -     1.019/0.978     1.023/0.976    -    1.019/0.979     1.022/0.976     -     1.021/0.978     1.024/0.975    -     1.02/0.978      1.023/0.976     -";
  datacard << "\t1.020/0.975\t1.025/0.975\t-  \t1.020/0.975\t1.025/0.975\t-" << endl;

  
  
  datacard << setprecision(4) << "# mc statistics for signal:" << endl;
  datacard << setprecision(4) << "mcstatst1 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(4) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(4) << 1 + (hnominal[i][j]->GetBinError(7)/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(4) << 1 + (hnominal[i][j]->GetBinError(8)/hnominal[i][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0 && i != 1){
	if (hnominaldy[i][j]->GetBinContent(3) !=0) datacard << setprecision(4) << 1 + (hnominaldy[i][j]->GetBinError(3)/hnominaldy[i][j]->GetBinContent(3)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else if (i! = 1) datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(4) << endl;
  
  datacard << setprecision(4) << "# mc statitics for background" << endl;
  datacard << setprecision(4) << "mcstatot2 lnN\t";
   for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(4) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(4) << 1 + (hnominal[i][j]->GetBinError(7)/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(4) << 1 + (hnominal[i][j]->GetBinError(8)/hnominal[i][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1 && i != 1){
	if (hnominaldy[i][j]->GetBinContent(3) !=0) datacard << setprecision(4) << 1 + (hnominaldy[i][j]->GetBinError(3)/hnominaldy[i][j]->GetBinContent(3)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else if (i! = 1) datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(4) << endl;
  
  datacard << setprecision(4) << "mcstatot3 lnN\t";
   for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(4) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(4) << 1 + (hnominal[i][j]->GetBinError(7)/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(4) << 1 + (hnominal[i][j]->GetBinError(8)/hnominal[i][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2 && i != 1){
	if (hnominaldy[i][j]->GetBinContent(3) !=0) datacard << setprecision(4) << 1 + (hnominaldy[i][j]->GetBinError(3)/hnominaldy[i][j]->GetBinContent(3)) << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else if (i! = 1)  datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(4) << endl;
  
}



