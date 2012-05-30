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
void datacardmaker_evolution(){

  ofstream datacard("singletop_tW_5fb_ultimate.txt"); 
 
  char myRootFile[300];
  sprintf(myRootFile,"results/Histos_cutbased_CR.root");
  
  TFile *_file0 = TFile::Open(myRootFile);
  
  const int nProcess = 4;
  TString processName[nProcess+1] =  { "data", "twdr", "tt", "zjets", "others"};
  
  TH1F*  hdata [3];
  TH1F*  hnominal [3][3];
  
  for (int process = 0; process < nProcess; process++){
    for (int i = 0; i < 3; i++){
      
      int mode = 0;
      if (i < 2) mode = i+1;
      char modeName[300];
      if (i == 0) sprintf(modeName,"emu");
      else if (i == 1) sprintf(modeName,"mumu");
      else sprintf(modeName,"ee");
      char title[300];
      sprintf(title,"R_%s_",modeName);
      
      if (process == 0) hdata[mode] = (TH1F*) (_file0->Get(title + processName[process]))->Clone();
      else {
	hnominal[mode][process-1] = (TH1F*) (_file0->Get(title + processName[process]))->Clone();
	if (process == 3){
	  TH1F* h1 = (TH1F*) (_file0->Get(title + processName[process+1]))->Clone();
	  hnominal[mode][process-1]->Add(h1);
	}
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
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
    
    for (int j = 0; j < 3; j++){
      if (j == 0){ 
        char title[300];
        sprintf(title,"R_%s_tw_sup",modeName);
        hup[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
        sprintf(title,"R_%s_tw_sdo",modeName);
	hdown[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
      } else if (j == 1) {
        char title[300];
        sprintf(title,"R_%s_tt_scaleup",modeName);
        hup[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
        sprintf(title,"R_%s_tt_scaledown",modeName);
	hdown[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
      }
    }
  } 
  
  datacard << setprecision(3) << "ttscale   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = ((hup[0][j]->GetBinContent(2) + hup[1][j]->GetBinContent(2) + hup[2][j]->GetBinContent(2)) + (hdown[0][j]->GetBinContent(2) +  hdown[1][j]->GetBinContent(2) + hdown[2][j]->GetBinContent(2)))/2;
	if (average !=0) datacard << setprecision(3) << 1 + ((hup[0][j]->GetBinContent(2) + hup[1][j]->GetBinContent(2) + hup[2][j]->GetBinContent(2)) - average)/average << "\t ";
	else datacard << setprecision(3) << "-\t ";
      } 
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = ((hup[0][j]->GetBinContent(7) + hup[1][j]->GetBinContent(7) + hup[2][j]->GetBinContent(7)) + (hdown[0][j]->GetBinContent(7) +  hdown[1][j]->GetBinContent(7) + hdown[2][j]->GetBinContent(7)))/2;
	if (average !=0) datacard << setprecision(3) << 1 + ((hup[0][j]->GetBinContent(7) + hup[1][j]->GetBinContent(7) + hup[2][j]->GetBinContent(7)) - average)/average << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = ((hup[0][j]->GetBinContent(8) + hup[1][j]->GetBinContent(8) + hup[2][j]->GetBinContent(8)) + (hdown[0][j]->GetBinContent(8) +  hdown[1][j]->GetBinContent(8) + hdown[2][j]->GetBinContent(8)))/2;
	if (average !=0) datacard << setprecision(3) << 1 + ((hup[0][j]->GetBinContent(8) + hup[1][j]->GetBinContent(8) + hup[2][j]->GetBinContent(8)) - average)/average << "\t ";
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
	double average = ((hup[0][j]->GetBinContent(2) + hup[1][j]->GetBinContent(2) + hup[2][j]->GetBinContent(2)) + (hdown[0][j]->GetBinContent(2) +  hdown[1][j]->GetBinContent(2) + hdown[2][j]->GetBinContent(2)))/2;
	if (average !=0) datacard << setprecision(3) << 1 + ((hup[0][j]->GetBinContent(2) + hup[1][j]->GetBinContent(2) + hup[2][j]->GetBinContent(2)) - average)/average << "\t ";
	else datacard << setprecision(3) << "-\t ";
      } 
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = ((hup[0][j]->GetBinContent(7) + hup[1][j]->GetBinContent(7) + hup[2][j]->GetBinContent(7)) + (hdown[0][j]->GetBinContent(7) +  hdown[1][j]->GetBinContent(7) + hdown[2][j]->GetBinContent(7)))/2;
	if (average !=0) datacard << setprecision(3) << 1 + ((hup[0][j]->GetBinContent(7) + hup[1][j]->GetBinContent(7) + hup[2][j]->GetBinContent(7)) - average)/average << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = ((hup[0][j]->GetBinContent(8) + hup[1][j]->GetBinContent(8) + hup[2][j]->GetBinContent(8)) + (hdown[0][j]->GetBinContent(8) +  hdown[1][j]->GetBinContent(8) + hdown[2][j]->GetBinContent(8)))/2;
	if (average !=0) datacard << setprecision(3) << 1 + ((hup[0][j]->GetBinContent(8) + hup[1][j]->GetBinContent(8) + hup[2][j]->GetBinContent(8)) - average)/average << "\t ";
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
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else  sprintf(modeName,"ee");
    
    for (int j = 0; j < 3; j++){
      if (j == 1) {
        char title[300];
        sprintf(title,"R_%s_tt_matchingup",modeName);
        hup[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
        sprintf(title,"R_%s_tt_matchingdown",modeName);
	hdown[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
      }
    }
  } 
 
  
  datacard << setprecision(3) << "ttmatch   lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = ((hup[0][j]->GetBinContent(2) + hup[1][j]->GetBinContent(2) + hup[2][j]->GetBinContent(2)) + (hdown[0][j]->GetBinContent(2) +  hdown[1][j]->GetBinContent(2) + hdown[2][j]->GetBinContent(2)))/2;
	if (average !=0) datacard << setprecision(4) << 1 + ((hup[0][j]->GetBinContent(2) + hup[1][j]->GetBinContent(2) + hup[2][j]->GetBinContent(2)) - average)/average << "\t ";
	else datacard << setprecision(4) << "-\t ";
      } 
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = ((hup[0][j]->GetBinContent(7) + hup[1][j]->GetBinContent(7) + hup[2][j]->GetBinContent(7)) + (hdown[0][j]->GetBinContent(7) +  hdown[1][j]->GetBinContent(7) + hdown[2][j]->GetBinContent(7)))/2;
	if (average !=0) datacard << setprecision(4) << 1 + ((hup[0][j]->GetBinContent(7) + hup[1][j]->GetBinContent(7) + hup[2][j]->GetBinContent(7)) - average)/average << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 1){
	double average = ((hup[0][j]->GetBinContent(8) + hup[1][j]->GetBinContent(8) + hup[2][j]->GetBinContent(8)) + (hdown[0][j]->GetBinContent(8) +  hdown[1][j]->GetBinContent(8) + hdown[2][j]->GetBinContent(8)))/2;
	if (average !=0) datacard << setprecision(4) << 1 + ((hup[0][j]->GetBinContent(8) + hup[1][j]->GetBinContent(8) + hup[2][j]->GetBinContent(8)) - average)/average << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  
  datacard << setprecision(3) << endl;
  
  TH1F*  h [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;  
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
     
    for (int j = 0; j < 3; j++){
      if (j == 0){ 
	char title[300];
	sprintf(title,"R_%s_twds",modeName);
	h[mode][j] = (TH1F*) (_file0->Get(title))->Clone();
      }
    }
  } 
   
   
  datacard << setprecision(4) << "twdrds    lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = hnominal[0][j]->GetBinContent(2) + hnominal[1][j]->GetBinContent(2) + hnominal[2][j]->GetBinContent(2);
	if (average !=0) datacard << setprecision(4) << 1 + ((h[0][j]->GetBinContent(2) + h[1][j]->GetBinContent(2) + h[2][j]->GetBinContent(2)) - average)/average << "\t ";
	else datacard << setprecision(4) << "-\t ";
      } 
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = hnominal[0][j]->GetBinContent(7) + hnominal[1][j]->GetBinContent(7) + hnominal[2][j]->GetBinContent(7);
	if (average !=0) datacard << setprecision(4) << 1 + ((h[0][j]->GetBinContent(7) + h[1][j]->GetBinContent(7) + h[2][j]->GetBinContent(7)) - average)/average << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 0){
	double average = hnominal[0][j]->GetBinContent(8) + hnominal[1][j]->GetBinContent(8) + hnominal[2][j]->GetBinContent(8);
	if (average !=0) datacard << setprecision(4) << 1 + ((h[0][j]->GetBinContent(8) + h[1][j]->GetBinContent(8) + h[2][j]->GetBinContent(8)) - average)/average << "\t ";
	else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(4) << endl;
   
  TString SystName = "PU";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
     
    for (int j = 0; j < 3; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      if (j == 2){
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__plus"))->Clone();
	hup[mode][j]->Add(h1);
      }
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
      if (j == 2){
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__minus"))->Clone();
	hdown[mode][j]->Add(h1);
      } 
    }
  } 
   
  datacard << setprecision(3) << "pu        lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
             
	if (hnominal[i][j]->GetBinContent(2) !=0){
	  double upv = (hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
	  double downv =  (hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
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
	 
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
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
	 
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
	} else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(4) << endl;
   
   
  TString SystName = "JES";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
     
    for (int j = 0; j < 3; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      if (j == 2){
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__plus"))->Clone();
	hup[mode][j]->Add(h1);
      }
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
      if (j == 2){
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__minus"))->Clone();
	hdown[mode][j]->Add(h1);
      }  
    }
  } 
   
  datacard << setprecision(3) << "jes       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnominal[i][j]->GetBinContent(2) !=0){
	  double upv = (hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
	  double downv =  (hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
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
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
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
	  
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
	} else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(4) << endl;
   
   
  TString SystName = "btag";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
     
    for (int j = 0; j < 3; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      if (j == 2){
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__plus"))->Clone();
	hup[mode][j]->Add(h1);
      }
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
      if (j == 2){
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__minus"))->Clone();
	hdown[mode][j]->Add(h1);
      } 
    }
  } 
   
  datacard << setprecision(4) << "btag      lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnominal[i][j]->GetBinContent(2) !=0){
	  double upv = (hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
	  double downv =  (hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	  datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
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
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
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
	 
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
	} else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  } datacard << setprecision(4) << endl;  
  
  TString SystName = "JER";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
     
    for (int j = 0; j < 3; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      if (j == 2){
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__plus"))->Clone();
	hup[mode][j]->Add(h1);
      }
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
      if (j == 2){
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__minus"))->Clone();
	hdown[mode][j]->Add(h1);
      } 
    }
  } 
  datacard << setprecision(4) << "jer       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnominal[i][j]->GetBinContent(2) !=0){
	  double upv = (hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
	  double downv =  (hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	 
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
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
	  
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
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
	 
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
	} else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(4) << endl;
   
   
  TString SystName = "MET";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
    int mode = 0;
    if (i < 2) mode = i+1;
    char modeName[300];
    if (i == 0) sprintf(modeName,"emu");
    else if (i == 1) sprintf(modeName,"mumu");
    else sprintf(modeName,"ee");
    
    for (int j = 0; j < 3; j++){
      char title[300];
      sprintf(title,"R_%s_",modeName);
      hup[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__plus" ))->Clone();
      if (j == 2){
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__plus"))->Clone();
	hup[mode][j]->Add(h1);
      }
      hdown[mode][j] = (TH1F*) (_file0->Get(title + processName[j+1]+ "__" + SystName + "__minus"))->Clone();
      if (j == 2){
	TH1F* h1 = (TH1F*) (_file0->Get(title + processName[j+2]+ "__" + SystName + "__minus"))->Clone();
	hdown[mode][j]->Add(h1);
      }
    }
  } 
   
  datacard << setprecision(3) << "met       lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j < 3){
	if (hnominal[i][j]->GetBinContent(2) !=0){
	  double upv = (hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
	  double downv =  (hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2);
	  double maxv = TMath::Max(fabs(upv), fabs(downv));
	
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
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
	
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
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
	  if (fabs(upv) >= 0.001 && fabs(downv) >= 0.001)datacard << setprecision(4) << 1 + upv << "/" << 1+ downv<< "\t ";
	  else if (fabs(upv) < 0.001 && fabs(downv) >= 0.001 && downv > -1)datacard << setprecision(4) << 1+ downv<< "\t ";
	  else if (fabs(downv) < 0.001 && fabs(upv) >= 0.001 && upv > -1)datacard << setprecision(4) << 1+ upv<< "\t ";
	  else datacard << setprecision(4) << "-\t ";
	} else datacard << setprecision(4) << "-\t ";
      }
      else datacard << setprecision(4) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
  datacard << setprecision(3) << "pdf       lnN\t1.020/0.978     1.024/0.975     -     1.018/0.980      1.024/0.975     -     1.018/0.98     1.024/0.975     -    1.025/0.974    1.024/0.975    -     1.021/0.977     1.023/0.976    -     1.019/0.978     1.023/0.976    -    1.019/0.979     1.022/0.976     -     1.021/0.978     1.024/0.975    -     1.02/0.978      1.023/0.976     -" << endl;

   
  datacard << setprecision(3) << "dynorm    lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==0 && j == 2)datacard << setprecision(3) << "1.33\t ";
      if(i ==1 && j == 2)datacard << setprecision(3) << "1.30\t ";
      if(i ==2 && j == 2)datacard << setprecision(3) << "1.27\t ";
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2)datacard << setprecision(3) << "1.3\t ";
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(j == 2)datacard << setprecision(3) << "-\t ";
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
  datacard << "ttxs      lnN\t -\t 1.06\t -\t -\t1.06\t -\t -\t 1.06\t -\t -\t 1.06\t -\t -\t 1.06\t -\t -\t 1.06\t - \t -\t 1.06 \t - \t -\t 1.06 \t - \t- \t1.06 \t -" << endl;
  datacard << "dyxs      lnN\t -\t -\t 1.05\t  -\t -\t1.05\t -\t -\t 1.05\t -\t -\t 1.05\t -\t -\t 1.05\t -\t -\t 1.05\t - \t -\t 1.05 \t - \t -\t 1.05 \t - \t- \t1.05 " << endl;
   
  datacard << setprecision(3) << "# mc statistics for signal:" << endl;
  datacard << setprecision(3) << "mcstatst1 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==0 && j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatst2 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==0 && j == 0){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(7)/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatst3 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==0 && j == 0){
	if (hnominal[i][j]->GetBinContent(7) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(7)/hnominal[i][j]->GetBinContent(7)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatst4 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==1 && j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){  
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
  datacard << setprecision(3) << "mcstatst5 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==1 && j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){  
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
  datacard << setprecision(3) << "mcstatst6 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){  
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==1 && j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
   
  datacard << setprecision(3) << "mcstatst7 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==2 && j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
   
  datacard << setprecision(3) << "mcstatst8 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==2 && j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
   
  datacard << setprecision(3) << "mcstatst9 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i ==2 && j == 0){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
  datacard << setprecision(3) << "# mc statitics for background" << endl;
  datacard << setprecision(3) << "mcstatot1 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 1){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
  datacard << setprecision(3) << "mcstatot2 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 1){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot3 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 1){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot4 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 1){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot5 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 1){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot6 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 1){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot7 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 1){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot8 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 1){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot9 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 1){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot10 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot11 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot12 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 0 && j == 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot13 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot14 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  datacard << setprecision(3) << "mcstatot15 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 1 && j == 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
  datacard << setprecision(3) << "mcstatot16 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
  datacard << setprecision(3) << "mcstatot17 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
   
  datacard << setprecision(3) << "mcstatot18 lnN\t";
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      datacard << setprecision(3) << "-\t ";
    }
  }
  for (int i = 0; i < 3; i++){
    for (int j = 0; j <3; j++){ 
      if(i == 2 && j == 2){
	if (hnominal[i][j]->GetBinContent(2) !=0) datacard << setprecision(3) << 1 + (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2)) << "\t ";
	else datacard << setprecision(3) << "-\t ";
      }
      else datacard << setprecision(3) << "-\t ";
    }
  }
  datacard << setprecision(3) << endl;
  
  
  
}



