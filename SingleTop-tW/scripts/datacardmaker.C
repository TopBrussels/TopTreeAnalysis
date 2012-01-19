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

 ofstream datacard("singletop_tW45.txt"); 

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
   }
 } 

 datacard << "# this is the version with *exclusive* jet / tag bins " << endl;
 datacard << "# based on 4.5/fb " << endl;
 datacard << "imax 9 # number of bins " << endl;
 datacard << "jmax 2 # number of processes - 1 " << endl;
 datacard << "kmax * # number of uncertainties " << endl;
 datacard << "------------ " << endl;
 datacard << "bin            ee1j1t  emu1j1t  mumu1j1t ee2j1t    emu2j1t  mumu2j1t  ee2j2t  emu2j2t  mumu2j2t " << endl;
 datacard << "observation   " ;
 
 for (int i = 0; i < 3; i++){
 datacard << "	" << hdata[i]->GetBinContent(2) ;
 }
 for (int i = 0; i < 3; i++){
 datacard << "	" << hdata[i]->GetBinContent(7) ;
 }
 for (int i = 0; i < 3; i++){
 datacard << "	" << hdata[i]->GetBinContent(8) ;
 }
 
 datacard << endl;
 datacard << "------------ " << endl;
 datacard << "bin            ee1j1t       ee1j1t       ee1j1t    emu1j1t      emu1j1t     emu1j1t    mumu1j1t     mumu1j1t     mumu1j1t    ee2j1t      ee2j1t       ee2j1t    emu2j1t      emu2j1t     emu2j1t   mumu2j1t     mumu2j1t     mumu2j1t  ee2j2t     ee2j2t      ee2j2t   emu2j2t     emu2j2t      emu2j2t  mumu2j2t     mumu2j2t     mumu2j2t" << endl;
 datacard << "process        st           tt           other     st           tt          other      st           tt           other       st          tt           other     st           tt          other      st           tt           other     st         tt          other    st          tt           other    st           tt           other" << endl;
 datacard << "process        0            1            2         0            1           2          0            1            2           0           1            2         0            1           2          0            1            2         0          1           2        0           1            2        0            1            2" << endl;
 datacard << "rate           " ;

 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     datacard << hnominal[i][j]->GetBinContent(2) << "          ";
   }
 }
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     datacard << hnominal[i][j]->GetBinContent(7) << "          ";
   }
 }
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     datacard << hnominal[i][j]->GetBinContent(8) << "          ";
   }
 }
 datacard << endl;
 datacard << "---------------------- " << endl;
 datacard << "lumi      lnN  1.045        1.045        -         1.045        1.045       1.045      1.045        1.045        -           1.045       1.045        1.045     1.045        1.045       1.045      1.045        1.045        1.045     1.045      1.045       1.045    1.045       1.045        1.045    1.045        1.045        1.045 " << endl;
 datacard << "hlte      lnN  1.015        1.015        -         1.011        1.011       1.011      -            -            -           1.015       1.015        1.015     1.011        1.011       1.011      -            -            -         1.015      1.015       1.015    1.011       1.011        1.011    -            -            - " << endl;
 datacard << "hltmu     lnN  -            -            -         1.011        1.011       1.011      1.015        1.015        -           -           -            -         1.011        1.011       1.011      1.015        1.015        1.015     -          -           -        1.011       1.011        1.011    1.015        1.015        1.015 " << endl;
 datacard << "ele       lnN  1.02         1.02         -         1.02         1.02        1.02       -            -            -           1.02        1.02         1.02      1.02         1.02        1.02       -            -            -         1.02       1.02        1.02     1.02        1.02         1.02     -            -            - " << endl;
 datacard << "mu        lnN  -            -            -         1.01         1.01        1.01       1.01         1.01         -           -           -            -         1.01         1.01        1.01       1.01         1.01         1.01      -          -           -        1.01        1.01         1.01     1.01         1.01         1.01 " << endl;
 datacard.precision(3);
 
 TH1F*  hup [3][3];
 TH1F*  hdown [3][3];
 for (int i = 0; i < 3; i++){
   int mode = 0;
   if (i < 2) mode = i+1;
   for (int j = 0; j < 3; j++){
     if (j == 0){
       sprintf(myRootFile,"outputs/noPU_%d_tbar_sup.root", i);
       TFile* _file1 = TFile::Open(myRootFile);
       hup[mode][j] = (TH1F*) _file1->Get("R");
       sprintf(myRootFile,"outputs/noPU_%d_tbar_sdo.root", i);
       TFile* _file1 = TFile::Open(myRootFile);
       hdown[mode][j] = (TH1F*) _file1->Get("R");
     } else if (j == 1) {
       sprintf(myRootFile,"outputs/noPU_%d_tt_scaleup.root", i);
       TFile* _file1 = TFile::Open(myRootFile);
       hup[mode][j] = (TH1F*) _file1->Get("R");
       sprintf(myRootFile,"outputs/noPU_%d_tt_scaledown.root", i);
       TFile* _file1 = TFile::Open(myRootFile);
       hdown[mode][j] = (TH1F*) _file1->Get("R");
     }
   }
 } 
 
 datacard << "ttscale   lnN  ";
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j < 2){
       double average = (hup[i][j]->GetBinContent(2) + hdown[i][j]->GetBinContent(2))/2;
       if (average !=0) datacard << 1 + fabs((hup[i][j]->GetBinContent(2) - average)/average) << "         ";
       else datacard << "-         ";
     }
     else datacard << "-         ";
   }
 }
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j < 2){
       double average = (hup[i][j]->GetBinContent(7) + hdown[i][j]->GetBinContent(7))/2;
       if (average !=0) datacard << 1 + fabs((hup[i][j]->GetBinContent(7) - average)/average) << "         ";
       else datacard << "-         ";
     }
     else datacard << "-         ";
   }
 }
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j < 2){
       double average = (hup[i][j]->GetBinContent(8) + hdown[i][j]->GetBinContent(8))/2;
       if (average !=0) datacard << 1 + fabs((hup[i][j]->GetBinContent(8) - average)/average) << "         ";
       else datacard << "-         ";
     }
     else datacard << "-         ";
   }
 }
 datacard << endl;
 
 TH1F*  hup [3][3];
 TH1F*  hdown [3][3];
 for (int i = 0; i < 3; i++){
   int mode = 0;
   if (i < 2) mode = i+1;
   for (int j = 0; j < 3; j++){
     if (j == 1) {
       sprintf(myRootFile,"outputs/noPU_%d_tt_matchingup.root", i);
       TFile* _file1 = TFile::Open(myRootFile);
       hup[mode][j] = (TH1F*) _file1->Get("R");
       sprintf(myRootFile,"outputs/noPU_%d_tt_matchingdown.root", i);
       TFile* _file1 = TFile::Open(myRootFile);
       hdown[mode][j] = (TH1F*) _file1->Get("R");
     }
   }
 } 

 datacard << "ttmatch   lnN  ";
  for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j == 1){
       double average = (hup[i][j]->GetBinContent(2) + hdown[i][j]->GetBinContent(2))/2;
       if (average !=0) datacard << 1 + fabs((hup[i][j]->GetBinContent(2) - average)/average) << "         ";
       else datacard << "-            ";
     }
     else datacard << "-            ";
   }
 }
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j == 1){
       double average = (hup[i][j]->GetBinContent(7) + hdown[i][j]->GetBinContent(7))/2;
       if (average !=0) datacard << 1 + fabs((hup[i][j]->GetBinContent(7) - average)/average) << "         ";
       else datacard << "-            ";
     }
     else datacard << "-            ";
   }
 }
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j == 1){
       double average = (hup[i][j]->GetBinContent(8) + hdown[i][j]->GetBinContent(8))/2;
       if (average !=0) datacard << 1 + fabs((hup[i][j]->GetBinContent(8) - average)/average) << "         ";
       else datacard << "-            ";
     }
     else datacard << "-            ";
   }
 }
 datacard << endl;

 TH1F*  hup [3][3];
 TH1F*  hdown [3][3];
 for (int i = 0; i < 3; i++){
   int mode = 0;
   if (i < 2) mode = i+1;
   for (int j = 0; j < 3; j++){
     if (j == 1) {
       sprintf(myRootFile,"outputs/noPU_%d_tt_largeISR.root", i);
       TFile* _file1 = TFile::Open(myRootFile);
       hup[mode][j] = (TH1F*) _file1->Get("R");
       sprintf(myRootFile,"outputs/noPU_%d_tt_smallISR.root", i);
       TFile* _file1 = TFile::Open(myRootFile);
       hdown[mode][j] = (TH1F*) _file1->Get("R");
     }
   }
 } 

 datacard << "ttifsr    lnN  ";
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j == 1){
       double average = (hup[i][j]->GetBinContent(2) + hdown[i][j]->GetBinContent(2))/2;
       if (average !=0) datacard << 1 + fabs((hup[i][j]->GetBinContent(2) - average)/average) << "         ";
       else datacard << "-            ";
     }
     else datacard << "-            ";
   }
 }
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j == 1){
       double average = (hup[i][j]->GetBinContent(7) + hdown[i][j]->GetBinContent(7))/2;
       if (average !=0) datacard << 1 + fabs((hup[i][j]->GetBinContent(7) - average)/average) << "         ";
       else datacard << "-            ";
     }
     else datacard << "-            ";
   }
 }
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j == 1){
       double average = (hup[i][j]->GetBinContent(8) + hdown[i][j]->GetBinContent(8))/2;
       if (average !=0) datacard << 1 + fabs((hup[i][j]->GetBinContent(8) - average)/average) << "         ";
       else datacard << "-            ";
     }
     else datacard << "-            ";
   }
 }
 datacard << endl;
 
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
 
 datacard << "twdrds    lnN  ";
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j == 0){
       if (hnominal[i][j]->GetBinContent(2) !=0) datacard << 1 + fabs((h[i][j]->GetBinContent(2) - hnominal[i][j]->GetBinContent(2))/hnominal[i][j]->GetBinContent(2)) << "         ";
       else datacard << "-            ";
     }
     else datacard << "-            ";
   }
 }
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j == 0){
       if (hnominal[i][j]->GetBinContent(7) !=0) datacard << 1 + fabs((h[i][j]->GetBinContent(7) - hnominal[i][j]->GetBinContent(7))/hnominal[i][j]->GetBinContent(7)) << "         ";
       else datacard << "-            ";
     }
     else datacard << "-            ";
   }
 }
 for (int i = 0; i < 3; i++){
   for (int j = 0; j <3; j++){ 
     if(j == 0){
       if (hnominal[i][j]->GetBinContent(8) !=0) datacard << 1 + fabs((h[i][j]->GetBinContent(8) - hnominal[i][j]->GetBinContent(8))/hnominal[i][j]->GetBinContent(8)) << "         ";
       else datacard << "-            ";
     }
     else datacard << "-            ";
   }
 }
 datacard << endl;

  TString SystName = "PU";
  TH1F*  hup [3][3];
  TH1F*  hdown [3][3];
  for (int i = 0; i < 3; i++){
   int mode = 0;
   if (i < 2) mode = i+1;
   for (int j = 0; j < 3; j++){
     
       sprintf(myRootFile,"_%d_", i);
       TFile* _file1 = TFile::Open("outputs/" + SystName + "sysUp"+ myRootFile + processName[j] + ".root");
       /*
       
       hup[mode][j] = (TH1F*) _file1->Get("R");
       sprintf(myRootFile,"outputs/noPU_%d_tt_smallISR.root", i);
       TFile* _file1 = TFile::Open(myRootFile);
       hdown[mode][j] = (TH1F*) _file1->Get("R");
     */
   }
 } 




}


