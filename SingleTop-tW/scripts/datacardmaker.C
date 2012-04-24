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
void datacardmaker(){

  ofstream datacard("singletop_tW_5fb.txt"); 
 
  TString processName[3] =  { "twdr", "tt","others"};
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
      if (j == 2 && mode == 0){
       hnominal[mode][j]->SetBinContent(2,  hnominal[mode][j]->GetBinContent(2) + 21.5);
       hnominal[mode][j]->SetBinContent(7,  hnominal[mode][j]->GetBinContent(7) + 1.26);
       hnominal[mode][j]->SetBinContent(8,  hnominal[mode][j]->GetBinContent(8) + 0);
      }
      if (j == 2 && mode == 1){
       hnominal[mode][j]->SetBinContent(2,  hnominal[mode][j]->GetBinContent(2) + 62);
       hnominal[mode][j]->SetBinContent(7,  hnominal[mode][j]->GetBinContent(7) + 10.6);
       hnominal[mode][j]->SetBinContent(8,  hnominal[mode][j]->GetBinContent(8) + 0.6);
      }
      if (j == 2 && mode == 2){
       hnominal[mode][j]->SetBinContent(2,  hnominal[mode][j]->GetBinContent(2) + 26);
       hnominal[mode][j]->SetBinContent(7,  hnominal[mode][j]->GetBinContent(7) + 5.9);
       hnominal[mode][j]->SetBinContent(8,  hnominal[mode][j]->GetBinContent(8) + 1.2);
      }
    
    }
  } 
  
  
  datacard << setprecision(3) << "# this is the version with *exclusive* jet / tag bins " << endl;
  datacard << setprecision(3) << "# based on 4.9/fb " << endl;
  datacard << setprecision(3) << "imax 9 # number of bins " << endl;
  datacard << setprecision(3) << "jmax 2 # number of processes - 1 " << endl;
  datacard << setprecision(3) << "kmax * # number of uncertainties " << endl;
  datacard << setprecision(3) << "------------ " << endl;
  datacard << setprecision(3) << "bin \t\tee1j1t\t emu1j1t\t mumu1j1t\t ee2j1t\t emu2j1t\t mumu2j1t\t ee2j2t\t emu2j2t\t mumu2j2t " << endl;
  datacard << setprecision(3) << "observation" ;
  
  for (int i = 0; i < 3; i++){
    datacard << setprecision(3) << "\t " << hdata[i]->GetBinContent(2) ;
  }
  for (int i = 0; i < 3; i++){
    datacard << setprecision(3) << "\t " << hdata[i]->GetBinContent(7) ;
  }
  for (int i = 0; i < 3; i++){
    datacard << setprecision(3) << "\t " << hdata[i]->GetBinContent(8) ;
  }
  
  datacard << setprecision(3) << endl;
  datacard << setprecision(3) << "------------ " << endl;
  datacard << setprecision(3) << "bin        \tee1j1t\t ee1j1t\t ee1j1t\t emu1j1t\t emu1j1t\t emu1j1t\t mumu1j1t\t mumu1j1t\t mumu1j1t\t ee2j1t\t ee2j1t \t ee2j1t\t emu2j1t\t emu2j1t\t emu2j1t\t mumu2j1t\t mumu2j1t\t mumu2j1t\t ee2j2t\t ee2j2t\t ee2j2t\t emu2j2t\t emu2j2t\t emu2j2t\t mumu2j2t\t mumu2j2t\t mumu2j2t" << endl;
  datacard << setprecision(3) << "process    \tst\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt\t other\t st\t tt    \t other\t st\t tt\t other" << endl;
  datacard << setprecision(3) << "process    \t0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2\t 0\t 1\t 2" << endl;
  datacard << setprecision(3) << "rate       \t" ;
  
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << hnominal[i][j]->GetBinContent(2) << "\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << hnominal[i][j]->GetBinContent(7) << "\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << hnominal[i][j]->GetBinContent(8) << "\t ";
    }
  }
  datacard << setprecision(3) << endl;

  datacard << setprecision(3) << "---------------------- " << endl;
  datacard << setprecision(3) <<  "lumi      lnN\t1.022\t 1.022\t -\t1.022\t 1.022\t 1.022\t 1.022\t 1.022\t -\t 1.022\t 1.022\t 1.022\t 1.022\t 1.022\t 1.022\t1.022\t 1.022\t 1.022	1.022\t1.022\t 1.022    1.022\t 1.022\t 1.022    1.022\t 1.022\t 1.022 " << endl;
  datacard << setprecision(3) <<  "hlte      lnN\t1.015\t 1.015\t -\t1.011\t 1.011\t 1.011\t -\t  -\t  -\t 1.015\t 1.015\t 1.015\t 1.011\t 1.011\t 1.011\t-\t -\t -\t1.015\t1.015\t 1.015    1.011\t 1.011\t 1.011	-\t	-\t	- " << endl;
  datacard << setprecision(3) <<  "hltmu     lnN\t-\t  -\t  -\t1.011\t 1.011\t 1.011\t 1.015\t 1.015\t -\t -\t	-\t	-\t1.011\t 1.011\t 1.011	1.015\t 1.015\t 1.015\t -\t	-\t  -\t 1.011\t 1.011\t 1.011    1.015\t 1.015\t 1.015 " << endl;
  datacard << setprecision(3) <<  "ele       lnN\t1.02\t1.02\t-\t1.02\t1.02\t 1.02\t -\t  -\t -\t1.02\t 1.02\t1.02\t1.02\t1.02\t 1.02\t -\t -\t -\t1.02\t 1.02\t 1.02\t 1.02\t 1.02\t1.02\t -\t -\t - " << endl;
  datacard << setprecision(3) <<  "mu        lnN\t-\t  -\t  -\t1.01\t1.01\t 1.01\t 1.01\t1.01\t-\t  -\t  -\t -\t1.01\t1.01\t 1.01\t 1.01\t1.01\t1.01	-\t -\t  -\t 1.01\t 1.01\t1.01	 1.01\t1.01\t1.01 " << endl;
	

  
  
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    for (int j = 0; j < 3; j++){
      if (j == 0){
	sprintf(myRootFile,"outputs/out_%d_tbar_sup.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hup[mode][j] = (TH1F*) _file1->Get("R");
	sprintf(myRootFile,"outputs/out_%d_tbar_sdo.root", i);
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
  
  datacard << setprecision(3) << "ttscale   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(2) + hdown[i][j]->GetBinContent(2))/2;
	if (average !=0) datacard << setprecision(3) << 1 + fabs((hup[i][j]->GetBinContent(2) - average)/average) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      } 
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(7) + hdown[i][j]->GetBinContent(7))/2;
	if (average !=0) datacard << setprecision(3) << 1 + fabs((hup[i][j]->GetBinContent(7) - average)/average) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(8) + hdown[i][j]->GetBinContent(8))/2;
	if (average !=0) datacard << setprecision(3) << 1 + fabs((hup[i][j]->GetBinContent(8) - average)/average) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "twscale   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = (hup[i][j]->GetBinContent(2) + hdown[i][j]->GetBinContent(2))/2;
	if (average !=0) datacard << setprecision(3) << 1 + fabs((hup[i][j]->GetBinContent(2) - average)/average) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      } 
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = (hup[i][j]->GetBinContent(7) + hdown[i][j]->GetBinContent(7))/2;
	if (average !=0) datacard << setprecision(3) << 1 + fabs((hup[i][j]->GetBinContent(7) - average)/average) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = (hup[i][j]->GetBinContent(8) + hdown[i][j]->GetBinContent(8))/2;
	if (average !=0) datacard << setprecision(3) << 1 + fabs((hup[i][j]->GetBinContent(8) - average)/average) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  
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
  
  
  
  datacard << setprecision(3) << "ttmatch   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(2) + hdown[i][j]->GetBinContent(2))/2;
	if (average !=0) datacard << setprecision(3) << 1 + fabs((hup[i][j]->GetBinContent(2) - average)/average) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      } 
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(7) + hdown[i][j]->GetBinContent(7))/2;
	if (average !=0) datacard << setprecision(3) << 1 + fabs((hup[i][j]->GetBinContent(7) - average)/average) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = (hup[i][j]->GetBinContent(8) + hdown[i][j]->GetBinContent(8))/2;
	if (average !=0) datacard << setprecision(3) << 1 + fabs((hup[i][j]->GetBinContent(8) - average)/average) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  
  datacard << setprecision(3) << endl;
  
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
  
  datacard << setprecision(3) << "twdrds    lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + fabs((h[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + fabs((h[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(3) << 1 + fabs((h[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
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
  
  datacard << setprecision(3) << "pu        lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)), fabs((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) ) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)), fabs((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) ) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(3) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)), fabs((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) ) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << "ttxs      lnN\t -\t 1.06\t -\t -\t1.06\t -\t -\t 1.06\t -\t -\t 1.06\t -\t -\t 1.06\t -\t -\t 1.06\t - \t -\t 1.06 \t - \t -\t 1.06 \t - \t- \t1.06 \t -" << endl;
  
  
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
  
  datacard << setprecision(3) << "jes       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + ((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) << "/" << 1+ ((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 2){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + ((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) << "/" << 1+ ((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 2){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(3) << 1 + ((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) << "/" << 1+ ((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
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
  
  datacard << setprecision(3) << "btag      lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + ((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) << "/" << 1+ ((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 2){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + ((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) << "/" << 1+ ((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  } 
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 2){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(3) << 1 + ((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) << "/" << 1+ ((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;  
  
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
  
  datacard << setprecision(3) << "jer       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)), fabs((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) ) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 2){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)), fabs((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) ) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 2){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(3) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)), fabs((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) ) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  
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
  
  datacard << setprecision(3) << "met       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)), fabs((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) ) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)), fabs((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) ) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(3) << 1 + TMath::Max(fabs((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)), fabs((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) ) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "pdf       lnN\t1.020/0.978     1.024/0.975     -     1.018/0.980      1.024/0.975     -     1.018/0.98     1.024/0.975     -    1.025/0.974    1.024/0.975    -     1.021/0.977     1.023/0.976    -     1.019/0.978     1.023/0.976    -    1.019/0.979     1.022/0.976     -     1.021/0.978     1.024/0.975    -     1.02/0.978      1.023/0.976     -" << endl;
  
  datacard << setprecision(3) << "# mc statistics for signal:" << endl;
  datacard << setprecision(3) << "mcstatst1 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(7)/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	if (hnominal[i][j]->GetBinContent(8) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(8)/hnominal[i][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "# mc statitics for background" << endl;
  datacard << setprecision(3) << "mcstatot1 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <9; j++){ 
      if(j == 2){
	if (i == 0) datacard << setprecision(3) << "1.5      " ;
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot2 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <9; j++){ 
      if(j == 5 && i == 0){
	datacard << setprecision(3) << "1.2      ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot3 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <9; j++){ 
      if(j == 8 && i == 0){
	datacard << setprecision(3) << "1.5      ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  
  datacard << setprecision(3) << "mcstatot4 lnN\t";
  for (int i = 0; i < 9; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2 && i == 3){
	if (hnominal[0][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + (hnominal[0][j]->GetBinError(7)/hnominal[0][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot5 lnN\t";
  for (int i = 0; i < 9; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2 && i == 4){
	if (hnominal[1][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + (hnominal[1][j]->GetBinError(7)/hnominal[1][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
     else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  
  datacard << setprecision(3) << "mcstatot6 lnN\t";
  for (int i = 0; i < 9; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2 && i == 5){
	if (hnominal[2][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + (hnominal[2][j]->GetBinError(7)/hnominal[2][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot7 lnN\t";
  for (int i = 0; i < 9; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2 && i == 6){
	if (hnominal[0][j]->GetBinContent(8) !=0) datacard << setprecision(3) << 1 + (hnominal[0][j]->GetBinError(8)/hnominal[0][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot8 lnN\t";
  for (int i = 0; i < 9; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2 && i == 7){
	if (hnominal[1][j]->GetBinContent(8) !=0) datacard << setprecision(3) << 1 + (hnominal[1][j]->GetBinError(8)/hnominal[1][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
     }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  
  datacard << setprecision(3) << "mcstatot9 lnN\t";
  for (int i = 0; i < 9; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2 && i == 8){
	if (hnominal[2][j]->GetBinContent(8) !=0) datacard << setprecision(3) << 1 + (hnominal[2][j]->GetBinError(8)/hnominal[2][j]->GetBinContent(8)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
}


