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
#include "setTDRStyle.C"
using namespace std;

void syst(int mode = 0){
  
  char myRootFile[300];
  double lumi = 1000;
  if (mode == 0 )        lumi = 4626.297;
  else if ( mode == 1)   lumi = 4534.871;
  else if ( mode == 2)   lumi = 4593.348;
  
  TString processName[4] =  { "twdr", "tt", "zjets","others"};
  
  TH1F*  hnominal [4];
  cout << "Nominal values" << endl;
  
  for (int iProcess = 0; iProcess < 4; iProcess++){
    sprintf(myRootFile,"outputs/out_%d_", mode);
    
    TFile *_file0 = TFile::Open(myRootFile + processName[iProcess]+".root");
    
    hnominal[iProcess] = (TH1F*) _file0->Get("cutflow");
    if (iProcess<2)  cout << "[" << processName[iProcess] << "]:		" ;
    else  cout << "[" << processName[iProcess] << "]:	" ;
    cout  << hnominal[iProcess]->GetBinContent(12) << "+-" << hnominal[iProcess]->GetBinError(12) << "		stats: "<< hnominal[iProcess]->GetBinError(12)*100/hnominal[iProcess]->GetBinContent(12) << "%" << endl;
  }
  
  
  TString SystName = "JESsysUp";
  cout << endl;
  
  cout << "Systematic calculated: " << SystName << endl;
  TH1F*  h [4];
  for (int iProcess = 0; iProcess < 4; iProcess++){
    sprintf(myRootFile,"_%d_", mode);
    
    TFile *_file0 = TFile::Open("outputs/sys" + SystName + myRootFile + processName[iProcess]+".root");
    
    h[iProcess] = (TH1F*) _file0->Get("cutflow");
    if (iProcess<2)  cout << "[" << processName[iProcess] << "]:		" ;
    else  cout << "[" << processName[iProcess] << "]:	" ;
    cout  << (h[iProcess]->GetBinContent(12) -  hnominal[iProcess]->GetBinContent(12))*100/hnominal[iProcess]->GetBinContent(12) << "%		(" <<  h[iProcess]->GetBinContent(12) << "+-" << h[iProcess]->GetBinError(12) << ")" << endl;
    
  }
  
  
  
  
}
