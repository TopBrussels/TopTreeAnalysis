//Rebeca Gonzalez Suarez
//rebeca@cern.ch

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "inputs.h"


using namespace std;
void tables_region(int mode = 0){

  double lumi = luminosity;
  
  if (mode == 0 )        lumi =  11966.617;
  else if ( mode == 1)   lumi = 12067.294;
  else if ( mode == 2)   lumi = 12093.792;
  
  
  char myTexFile[300];
  sprintf(myTexFile,"tables/main_table_%d_%dpb.tex", mode, lumi);
  ofstream salida(myTexFile); 
  
  char myRootFile[300];
  sprintf(myRootFile,"../outputs/out_%dpb_%d.root", lumi, mode);
  TFile *_file0 = TFile::Open(myRootFile);
  
  const int nProcess = 19;
  const int nProcess_real = 14;
  TString processName[nProcess] = {  "twdr", "atwdr", "s", "as","t", "at", "tt","ww","wz","zz", "zjets","zjets_lowmll", "wjets",  "data","di","zjets","twdr","other","mc"};
  TString processLabel[nProcess] = { "\\textbf{$tW$}","\\textbf{$\bar{t}W$}","\\textbf{s-channel t}", "\\textbf{s-channel $\bar{t}$}", "\\textbf{t-channel t}", "\\textbf{t-channel $\bar{t}$}", "\\textbf{$t \\bar{t}$}",  "\\textbf{$WW$}", "\\textbf{$WZ$}", "\\textbf{$ZZ$}", , "\\textbf{$Z+jets$}",  "\\textbf{$Z+jets low$}","\\textbf{$W+jets$}", "\\textbf{$MuEG$}","\\textbf{$WW$, $WZ$, $ZZ$}",  "\\textbf{$Z+jets$}","\\textbf{$single$ $top$}", "\\textbf{$Other$}", "\\textbf{$MC$}"};
  
  if (mode == 1)        processLabel[14] = "\\textbf{$DoubleMu$}";
  else if (mode == 2)   processLabel[14] = "\\textbf{$DoubleElectron$}";
  
  const in nCutLabels = 10;
  TString cutLabel[nCutLabels] = {"blank0", "blank1", "All","HLT","Lepton Sel.", "Inv. Mass", "$E_{T}^{miss}$", "1 Jet", "b-tagging", "$H_{T}$"};
  
  TH1F*  h [nProcess];
  for(int i=0; i<nProcess; i++){
    if(i<nProcess_real){
    	h[i] = (TH1F*) _file0->Get("cutflow_"+processName[i]);
    }
    else if (i == 14) { //di
     // clone ww
     h[i] =  (TH1F*)h[7]->Clone();
     // add wz
     h[i]->Add(h[8]);
    // add zz
     h[i]->Add(h[9]);
    }
    else if (i == 15) { //zjets
     // clone zjets
     h[i] =  (TH1F*)h[10]->Clone();
     // add zjets low
     h[i]->Add(h[11]);
    }    
    else if (i == 15) { //twdr
     // clone twdr
     h[i] =  (TH1F*)h[0]->Clone();
     // add atwdr
     h[i]->Add(h[1]);
   
    }
    else if (i == 16) { //other
     // clone s
     h[i] =  (TH1F*)h[2]->Clone();
     // add as
     h[i]->Add(h[3]);
    // add t
     h[i]->Add(h[4]);
     // add at
     h[i]->Add(h[5]);
     // add di 
     h[i]->Add(h[14]);
     // add wjets
     h[i]->Add(h[12]);
    }    
  }
  
  int nCuts = 10; 
  int nBins = nCuts +1;
  double vectorValue[nProcess][nBins][4];

  for (int i = 2; i < nBins; i++){
     // loop over processes
    for (int j = 0; j < nProcess; j++){
      vectorValue[j][i][0] = h[j]->GetBinContent(i);
      vectorValue[j][i][1] = precision(h[j]->GetBinError(i));
      vectorValue[j][i][2] = h[j]->GetBinError(i);
    }  
  }

  salida << "\\documentclass{article}" << endl;
  salida << "\\begin{document}" << endl;
  salida << endl;
  salida << endl;
  
  // table with everything in detail
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|c|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  // titles all samples
  for (int i = 0; i < nProcess_real; i++){
    salida << " & " << processLabel[i] ;
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
  
  for (int i=2; i < nBins; i++){
    salida << cutLabel[i-1];
    for (int j = 0; j < nProcess_real; j++){
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
      
    }
    salida <<  " \\\\  " << endl; 
  }
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;


   // table with tw ttbar zgamma other
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 6; i < nProcess; i++){
    if (i == 15 || i == 16 || i ==6|| i == 17) salida << " & " << processLabel[i] ;
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
  
  for (int i=2; i < nBins; i++){
    salida << cutLabel[i-1];
    for (int j = 0; j < nProcess; j++){
       if(j == 6 || j == 15 || j == 16 || j == 17){
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][2];
	}
    }
    salida <<  " \\\\  " << endl;
  }
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  salida << endl;
  salida << endl;
  
 
 
  // table for data - MC 
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 12; i < nProcess; i++){
    if( i == 13 || i == 18){
     salida << " & " << processLabel[i] ;
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
  
  for (int i=2; i < nBins; i++){
    salida << cutLabel[i-1] ;
    for (int j = 12; j < nProcess; j++){
      if (j == 13){

	  salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][0] ; 
 
      } else if (j == 18){
	  salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][0] ; 
	  salida << " $\\pm $"  << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][2];
	} 
      }
      
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
