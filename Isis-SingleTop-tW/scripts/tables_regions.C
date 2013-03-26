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
void tables_regions(int mode = 0){

  double lumi = luminosity;
  
  if (mode == 0 )        lumi =  11966.617;
  else if ( mode == 1)   lumi = 12067.294;
  else if ( mode == 2)   lumi = 12093.792;
  
  
  char myTexFile[300];
  sprintf(myTexFile,"tables/region_table_%d_%dpb.tex", mode, lumi);
  ofstream salida(myTexFile); 
 
  char myRootFile[300];
  sprintf(myRootFile,"results/an_%dpb_%d.root", lumi, mode);
  TFile *_file0 = TFile::Open(myRootFile);
  
  const int nProcess = 9;
  TString processName[nProcess] = { "twdr", "tt", "di", "wjets", "zjets", "st",  "others", "data", "mc"};
  TString processLabel[nProcess] = { "\\textbf{$tW$}", "\\textbf{$t \\bar{t}$}", "\\textbf{$WW$, $WZ$, $ZZ$}", "\\textbf{$W+jets$}", "\\textbf{$Z+jets$}",
			       "\\textbf{$single$ $top$}", "\\textbf{$Other$}", "\\textbf{$MuEG$}", "\\textbf{$MC$}"};
  
  if (mode == 1)        processLabel[7] = "\\textbf{$DoubleMu$}";
  else if (mode == 2)   processLabel[7] = "\\textbf{$DoubleElectron$}";
  
  TString cutLabel[8] = { "blank0", "blank1", "1jet 1tag 1 loose", "2jets 1 tag 1 loose", "2jets 2tags 2 loose", "1 jet 1 tag no Ht cut" , "2jets 1tag no Ht cut", "2 jets 2 tags no Ht cut"};
  
  TH1F*  h [nProcess];
  for(int i=0; i<nProcess; i++){
    h[i] = (TH1F*) _file0->Get("R_"+processName[i]);
    cout << "Get for " << i << " th loop " << "R_" << processName[i] << endl;
  }
  
  double vectorValue[nProcess][10][3];
  for (int i = 1; i < 10; i++){
    for (int j = 0; j < nProcess; j++){
      vectorValue[j][i][0] = h[j]->GetBinContent(i);
      vectorValue[j][i][1] = precision(h[j]->GetBinError(i));
      vectorValue[j][i][2] = h[j]->GetBinError(i);
    }  
  }

  salida << "\\documentclass[a4paper,12pt]{article}" << endl;
  salida << "\\begin{document}" << endl;
  salida << endl;
  salida << endl;
  
  
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 0; i < 6; i++){
    salida << " & " << processLabel[i] ;
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
  
  for (int i=1; i <4; i++){
    
    	salida << cutLabel[i+1];
    
    for (int j = 0; j < 6; j++){
       if(i < 4 ){
	salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
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
  
  
    salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 0; i < 6; i++){
    salida << " & " << processLabel[i] ;
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
  
  for (int i=7; i <10; i++){
        salida << cutLabel[i-2];
    
    for (int j = 0; j < 6; j++){
       if(i > 6 ){
         salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
      
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

// tw tt zjets other
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 0; i < 7; i++){
    if( i ==0 || i == 1 || i == 4 || i ==6 ){
       salida << " & " << processLabel[i] ;
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
    
  for (int i=1; i < 4; i++){
    
    	salida << cutLabel[i+1];
    
    for (int j = 0; j < 7; j++){
        if(j < 2 || j == 4 || j ==6){
	 salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	 salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
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
  
  
    salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|c|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 0; i < 7; i++){
    if(i < 2 || i == 4 || i ==6){
      salida << " & " << processLabel[i] ;
    }
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
  
  for (int i=7; i <10; i++){
        salida << cutLabel[i-2];
    
    for (int j = 0; j < 7; j++){
       if(j < 2 || j == 4 || j ==6){
         salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ; 
	salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
      
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
  
  
  // data MC
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 7; i < 9; i++){
    
       salida << " & " << processLabel[i] ;
    
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
    
  for (int i=1; i < 4; i++){
    
    	salida << cutLabel[i+1];
    
    for (int j = 7; j < 9; j++){
        
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
  
  
    salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|c|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 7; i < 9; i++){
    
      salida << " & " << processLabel[i] ;
    
  }
  salida << "  \\\\ " << endl; 
  salida << "  \\hline " << endl;
  
  for (int i=7; i <10; i++){
        salida << cutLabel[i-2];
    
    for (int j = 7; j <9; j++){
       
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
