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
void syst(){

  TString processName[4] =  { "twdr", "tt", "zjets","others"};
  char myRootFile[300]; 
  TH1F*  hnominal [3][4];
  for (int mode = 0; mode < 3; mode++){
   for (int j = 0; j < 3; j++){
     sprintf(myRootFile,"outputs/out_%d_", mode);
     TFile *_file1 = TFile::Open(myRootFile + processName[j] + ".root");
     hnominal[mode][j] = (TH1F*) _file1->Get("R");
     /*
     if (j == 2 && mode == 0){
       hnominal[mode][j]->SetBinContent(2,  hnominal[mode][j]->GetBinContent(2) + 21.5);
       hnominal[mode][j]->SetBinContent(7,  hnominal[mode][j]->GetBinContent(7) + 2.5);
      }
      if (j == 2 && mode == 1){
       hnominal[mode][j]->SetBinContent(2,  hnominal[mode][j]->GetBinContent(2) + 62);
       hnominal[mode][j]->SetBinContent(7,  hnominal[mode][j]->GetBinContent(7) + 18.7);
       hnominal[mode][j]->SetBinContent(8,  hnominal[mode][j]->GetBinContent(8) + 1.2);
      }
      if (j == 2 && mode == 2){
       hnominal[mode][j]->SetBinContent(2,  hnominal[mode][j]->GetBinContent(2) + 26);
       hnominal[mode][j]->SetBinContent(7,  hnominal[mode][j]->GetBinContent(7) + 13.3);
       hnominal[mode][j]->SetBinContent(8,  hnominal[mode][j]->GetBinContent(8) + 3);
      }*/
     
   }
 } 
  cout << "Breakdown of the systematics " << endl;
  cout << "* Lumi " << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      2.2%" << endl;
  }
  cout << "* HLT (electron):" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      1.1%" << endl;
    else if (i == 1) cout << "mode:" << i << "      -" << endl;
    else if (i == 2) cout << "mode:" << i << "      1.5%" << endl;
  }
  cout << "* HLT (muon):" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      1.1%" << endl;
    else if (i == 1) cout << "mode:" << i << "      1.15%" << endl;
    else if (i == 2) cout << "mode:" << i << "      -" << endl;
  }
  cout << "* Electron ID:" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      2%" << endl;
    else if (i == 1) cout << "mode:" << i << "      -" << endl;
    else if (i == 2) cout << "mode:" << i << "      2%" << endl;
  }
  cout << "* Muon ID:" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      1%" << endl;
    else if (i == 1) cout << "mode:" << i << "      1%" << endl;
    else if (i == 2) cout << "mode:" << i << "      -" << endl;
  }
  
  
  cout.precision(2);
  
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
  
    for (int j = 0; j < 3; j++){
      if (j == 0){
	sprintf(myRootFile,"outputs/out_%d_tw_sup.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hup[i][j] = (TH1F*) _file1->Get("R");
	sprintf(myRootFile,"outputs/out_%d_tw_sdo.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hdown[i][j] = (TH1F*) _file1->Get("R");
      } else if (j == 1) {
	sprintf(myRootFile,"outputs/out_%d_tt_scaleup.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hup[i][j] = (TH1F*) _file1->Get("R");
	sprintf(myRootFile,"outputs/out_%d_tt_scaledown.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hdown[i][j] = (TH1F*) _file1->Get("R");
      }
    }
  } 
  
  
  cout << "* Factorization/Normalization Scale " << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j < 3; j++){ 
      if(j < 2){
	cout << processName[j] << ":" ;
        double average = ((hup[0][j]->GetBinContent(2) + hup[1][j]->GetBinContent(2) + hup[2][j]->GetBinContent(2)) + (hdown[0][j]->GetBinContent(2) +  hdown[1][j]->GetBinContent(2) + hdown[2][j]->GetBinContent(2)))/2;
	if (average !=0) cout << fabs(((hup[0][j]->GetBinContent(2) + hup[1][j]->GetBinContent(2) + hup[2][j]->GetBinContent(2)) - average)/average)*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j1t]" ;
     for (int j = 0; j < 3; j++){ 
      if(j < 2){
	cout << processName[j] << ":" ;
        double average = ((hup[0][j]->GetBinContent(7) + hup[1][j]->GetBinContent(7) + hup[2][j]->GetBinContent(7)) + (hdown[0][j]->GetBinContent(7) +  hdown[1][j]->GetBinContent(7) + hdown[2][j]->GetBinContent(7)))/2;
	if (average !=0) cout << fabs(((hup[0][j]->GetBinContent(7) + hup[1][j]->GetBinContent(7) + hup[2][j]->GetBinContent(7)) - average)/average)*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j2t]" ;
     for (int j = 0; j < 3; j++){ 
      if(j < 2){
	cout << processName[j] << ":" ;
        double average = ((hup[0][j]->GetBinContent(8) + hup[1][j]->GetBinContent(8) + hup[2][j]->GetBinContent(8)) + (hdown[0][j]->GetBinContent(8) +  hdown[1][j]->GetBinContent(8) + hdown[2][j]->GetBinContent(8)))/2;
	if (average !=0) cout << fabs(((hup[0][j]->GetBinContent(8) + hup[1][j]->GetBinContent(8) + hup[2][j]->GetBinContent(8)) - average)/average)*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << endl;
    
  }
  
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
  
    for (int j = 0; j < 3; j++){
      if (j == 1) {
	sprintf(myRootFile,"outputs/out_%d_tt_matchingup.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hup[i][j] = (TH1F*) _file1->Get("R");
	sprintf(myRootFile,"outputs/out_%d_tt_matchingdown.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	hdown[i][j] = (TH1F*) _file1->Get("R");
      }
    }
  } 
  
  
  cout << "* ME/PS matching thresholds" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j < 3; j++){ 
      if(j == 1){
	cout << processName[j] << ":" ;
        double average = ((hup[0][j]->GetBinContent(2) + hup[1][j]->GetBinContent(2) + hup[2][j]->GetBinContent(2)) + (hdown[0][j]->GetBinContent(2) +  hdown[1][j]->GetBinContent(2) + hdown[2][j]->GetBinContent(2)))/2;
	if (average !=0) cout << fabs(((hup[0][j]->GetBinContent(2) + hup[1][j]->GetBinContent(2) + hup[2][j]->GetBinContent(2)) - average)/average)*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j1t]" ;
     for (int j = 0; j < 3; j++){ 
      if(j == 1){
	cout << processName[j] << ":" ;
        double average = ((hup[0][j]->GetBinContent(7) + hup[1][j]->GetBinContent(7) + hup[2][j]->GetBinContent(7)) + (hdown[0][j]->GetBinContent(7) +  hdown[1][j]->GetBinContent(7) + hdown[2][j]->GetBinContent(7)))/2;
	if (average !=0) cout << fabs(((hup[0][j]->GetBinContent(7) + hup[1][j]->GetBinContent(7) + hup[2][j]->GetBinContent(7)) - average)/average)*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j2t]" ;
     for (int j = 0; j < 3; j++){ 
      if(j == 1){
	cout << processName[j] << ":" ;
        double average = ((hup[0][j]->GetBinContent(8) + hup[1][j]->GetBinContent(8) + hup[2][j]->GetBinContent(8)) + (hdown[0][j]->GetBinContent(8) +  hdown[1][j]->GetBinContent(8) + hdown[2][j]->GetBinContent(8)))/2;
	if (average !=0) cout << fabs(((hup[0][j]->GetBinContent(8) + hup[1][j]->GetBinContent(8) + hup[2][j]->GetBinContent(8)) - average)/average)*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << endl;
  }
 
  
  TH1F*  h [3][4];
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      if (j == 0) {
	sprintf(myRootFile,"outputs/out_%d_twds.root", i);
	TFile* _file1 = TFile::Open(myRootFile);
	h[i][j] = (TH1F*) _file1->Get("R");
      }
    }
  } 
  
  cout << "* DRDS scheme" << endl;
   for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j < 3; j++){ 
      if(j == 0){
	cout << processName[j] << ":" ;
        double average = hnominal[0][j]->GetBinContent(2) +  hnominal[1][j]->GetBinContent(2) + hnominal[2][j]->GetBinContent(2);
	if (average !=0) cout << fabs(((h[0][j]->GetBinContent(2) + h[1][j]->GetBinContent(2) + h[2][j]->GetBinContent(2)) - average)/average)*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j1t]" ;
     for (int j = 0; j < 3; j++){ 
      if(j == 0){
	cout << processName[j] << ":" ;
        double average = hnominal[0][j]->GetBinContent(7) +  hnominal[1][j]->GetBinContent(7) + hnominal[2][j]->GetBinContent(7);
	if (average !=0) cout << fabs(((h[0][j]->GetBinContent(7) + h[1][j]->GetBinContent(7) + h[2][j]->GetBinContent(7)) - average)/average)*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j2t]" ;
     for (int j = 0; j < 3; j++){ 
      if(j == 0){
	cout << processName[j] << ":" ;
        double average = hnominal[0][j]->GetBinContent(8) +  hnominal[1][j]->GetBinContent(8) + hnominal[2][j]->GetBinContent(8);
	if (average !=0) cout << fabs(((h[0][j]->GetBinContent(8) + h[1][j]->GetBinContent(8) + h[2][j]->GetBinContent(8)) - average)/average)*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << endl;
  }
  
  
  TString SystName = "PU";
  TH1F*  hup [3][4];
  TH1F*  hdown [3][4];
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hup[i][j] = (TH1F*) _file1->Get("R");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdown[i][j] = (TH1F*) _file1->Get("R"); 
    }
  } 
  
  cout << "* PU" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
       for (int j = 0; j <3; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(2) !=0) cout << TMath::Max(fabs((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2))*100 , fabs((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2))*100 )<< "%   " ;
	else cout << " -      ";
	
	//if (j == 2) cout << hup[i][j]->GetBinContent(2) << ", " << hnominal[i][j]->GetBinContent(2) << ", " <<hdown[i][j]->GetBinContent(2) << endl;
      }
    }cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(7) !=0) cout << TMath::Max(fabs((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7))*100 , fabs((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7))*100 )<< "%   " ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(8) !=0) cout << TMath::Max(fabs((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8))*100 , fabs((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8))*100 )<< "%   " ;
	else cout << " -\t";
      }
    }
    cout << endl;
  }
  
  cout << "* tt cross-section:" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      6%" << endl;
    else if (i == 1) cout << "mode:" << i << "      6%" << endl;
    else if (i == 2) cout << "mode:" << i << "      6" << endl;
  }
  
  TString SystName = "JES";
  TH1F*  hup [3][4];
  TH1F*  hdown [3][4];
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hup[i][j] = (TH1F*) _file1->Get("R");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdown[i][j] = (TH1F*) _file1->Get("R"); 
    }
  } 
  
  cout << "* JES" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <3; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(2) !=0) cout << ((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2))*100 << "%/" <<  ((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2))*100 << "%   " ;
	else cout << " -      ";
      }
    }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(7) !=0) cout << ((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7))*100 << "%/" << 
	((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7))*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(8) !=0) cout << ((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8))*100 << "%/" << 
	((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8))*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << endl;
  }

 TString SystName = "SF";
  TH1F*  hup [3][4];
  TH1F*  hdown [3][4];
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hup[i][j] = (TH1F*) _file1->Get("R");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdown[i][j] = (TH1F*) _file1->Get("R"); 
    }
  } 
  
  cout << "* B-tagging" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <3; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(2) !=0) cout << ((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2))*100 << "%/" <<  ((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2))*100 << "%   " ;
	else cout << " -      ";
      }
    }cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(7) !=0) cout << ((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7))*100 << "%/" << 
	((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7))*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(8) !=0) cout << ((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8))*100 << "%/" << 
	((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8))*100 << "%\t" ;
	else cout << " -\t";
      }
    }
    cout << endl;
  }

 TString SystName = "JER";
  TH1F*  hup [3][4];
  TH1F*  hdown [3][4];
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hup[i][j] = (TH1F*) _file1->Get("R");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdown[i][j] = (TH1F*) _file1->Get("R"); 
    }
  } 
  
  cout << "* JER" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
    for (int j = 0; j <3; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(2) !=0) cout << TMath::Max(fabs((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2))*100 , fabs((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2))*100 )<< "%   " ;
	else cout << " -      ";
      }
    }cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(7) !=0) cout << TMath::Max(fabs((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7))*100 , fabs((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7))*100 )<< "%   " ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(8) !=0) cout << TMath::Max(fabs((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8))*100 , fabs((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8))*100 )<< "%   " ;
	else cout << " -\t";
      }
    }
    cout << endl;
  }

 TString SystName = "MET";
  TH1F*  hup [3][4];
  TH1F*  hdown [3][4];
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      sprintf(myRootFile,"_%d_", i);
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
      hup[i][j] = (TH1F*) _file1->Get("R");
      TFile* _file1 = TFile::Open("outputs/" + SystName + "sysDown"+ myRootFile + processName[j] + ".root");
      hdown[i][j] = (TH1F*) _file1->Get("R"); 
    }
  } 
  
  cout << "* Missing ET (unclustered)" << endl;
  for (int i = 0; i < 3; i++){
    cout << "mode:" << i << "      " ;
      for (int j = 0; j <3; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(2) !=0) cout << TMath::Max(fabs((hup[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2))*100 , fabs((hdown[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2))*100 )<< "%   " ;
	else cout << " -      ";
      }
    }cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(7) !=0) cout << TMath::Max(fabs((hup[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7))*100 , fabs((hdown[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7))*100 )<< "%   " ;
	else cout << " -\t";
      }
    }
    cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if (j < 3){
	cout << processName[j] << ":" ;
	if (hnominal[i][j]->GetBinContent(8) !=0) cout << TMath::Max(fabs((hup[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8))*100 , fabs((hdown[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8))*100 )<< "%   " ;
	else cout << " -\t";
      }
    }
    cout << endl;
  }

  cout << "* pdf (signal only, for all of them check https://fblekman.web.cern.ch/fblekman/documents/pdf_systs_cut-and-count.txt):" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      -2.0/1.8%" << endl;
    else if (i == 1) cout << "mode:" << i << "      -2.0/1.8%" << endl;
    else if (i == 2) cout << "mode:" << i << "      -2.2/2.0%" << endl;
  }
  

 cout << "* MC statistics for signal:" << endl;
 for (int i = 0; i < 3; i++){
   cout << "mode:" << i << "      " ;
   for (int j = 0; j <4; j++){ 
     if(j < 3){
       cout << processName[j] << ":" ;
       if (hnominal[i][j]->GetBinContent(2) !=0) cout << (hnominal[i][j]->GetBinError(2)/hnominal[i][j]->GetBinContent(2))*100 << "%\t" ;
       else cout << "-\t";
     }
   }
    cout << "\t[2j1t]" ;
    for (int j = 0; j <4; j++){ 
      if(j == 0 || j == 2){
	cout << processName[j] << ":" ;
       if (hnominal[i][j]->GetBinContent(7) !=0) cout << (hnominal[i][j]->GetBinError(7)/hnominal[i][j]->GetBinContent(7))*100 << "%\t" ;
       else cout << "-\t";
      }
    }
    cout << "\t[2j2t]" ;
    for (int j = 0; j <4; j++){ 
      if(j == 0 || j == 2){
	cout << processName[j] << ":" ;
       if (hnominal[i][j]->GetBinContent(8) !=0) cout << (hnominal[i][j]->GetBinError(8)/hnominal[i][j]->GetBinContent(8))*100 << "%\t" ;
       else cout << "-\t";
      }
    }
   cout << endl;
 }

 cout << "* z+jet background normalization:" << endl;
  for (int i = 0; i < 3; i++){
    if (i == 0) cout << "mode:" << i << "      50%\t[2j1t]: 50%\t" << endl;
    else if (i == 1) cout << "mode:" << i << "      50%\t[2j1t]: 50%\t" << endl;
    else if (i == 2) cout << "mode:" << i << "      50%\t[2j1t]: 50%\t" << endl;
  } 
   
}


