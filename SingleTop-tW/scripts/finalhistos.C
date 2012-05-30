#include "shapelooper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "inputs.h"
#include "Riostream.h"
#include <vector>
#include <string>
#include "TFile.h"
#include "TChain.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

void finalhistos(){

  char myRootFile[300];
  sprintf(myRootFile,"results/Histos_cutbased.root");
  
  TFile *_file0 = TFile::Open(myRootFile);
  cout << myRootFile << " read" << endl;
  
  const int nProcess = 4;
  TString processName[nProcess+1] =  { "data", "twdr", "tt", "zjets", "others"};
  TString newProcessName[nProcess] =  { "DATA", "st", "tt", "other"};
 
  const int nSyst = 6;
  TString systName[nSyst] =  { "null", "PU", "JES", "JER", "MET", "btag"};
  TString modeSyst[2] =  { "plus", "minus"};
    
  
  char newRootFile[300];
  sprintf(newRootFile,"results/Compact_histos.root");
  
  TFile f_var(newRootFile, "UPDATE");
  cout << newRootFile << " created" << endl;
  
  for (int mode =0; mode < 3; mode ++){
    for (int syst = 0; syst < nSyst; syst++){
      for (int process = 0; process < nProcess; process++){
        for (int i =0; i < 2; i++){
          
	  char modeName[300];
          if (mode == 0) sprintf(modeName,"emu");
          else if (mode == 1) sprintf(modeName,"mumu");
          else sprintf(modeName,"ee");
	  
	  char title[300];
	  char oldtitle[300];
	  
	  if (syst == 0){
	    sprintf(oldtitle,"%s1j1t__",modeName);
	    sprintf(title,"%s__",modeName);
	    
	    // cout <<oldtitle + processName[process] ;
	    
	    TH1F* h = (TH1F*) _file0->Get(oldtitle + processName[process]);
	    if (process == nProcess -1){
	      cout <<oldtitle << "  ";
	      TH1F* h1 = (TH1F*) _file0->Get(oldtitle + processName[process+1]);
	      cout << "check: " << h->Integral() <<" + " << h1->Integral()<< h->Integral() + h1->Integral();
	      h->Add(h1);
	      cout << " --> " << h->Integral() << endl; 
	     
	    }
	    TH1F* histo_1j1t = h->Clone();
	    histo_1j1t->SetName(title + newProcessName[process]);
	    //   cout << " saved as " << title + newProcessName[process] << endl;
	    i ++;
          } else if (process != 0){
	  
	    sprintf(oldtitle,"%s1j1t__",modeName);
	    sprintf(title,"%s__",modeName);
	    
	    // cout <<oldtitle + processName[process] + "__" + systName[syst] + "__" + modeSyst[i] ;
	    
	    TH1F* h = (TH1F*) _file0->Get(oldtitle + processName[process] + "__" + systName[syst] + "__" + modeSyst[i]);
	   
	    if (process == nProcess -1){
	      cout <<oldtitle  + systName[syst] + "__" + modeSyst[i]<< "  ";
	      TH1F* h1 = (TH1F*) _file0->Get(oldtitle + processName[process+1]+ "__" + systName[syst] + "__" + modeSyst[i]);
	      cout << "check: " << h->Integral() <<" + " << h1->Integral()<< h->Integral() + h1->Integral();
	      h->Add(h1);
	      cout << " --> " << h->Integral() << endl; 
	    }
	    TH1F* histo_1j1t = h->Clone();
	    histo_1j1t->SetName(title + newProcessName[process]+ "__" + systName[syst] + "__" + modeSyst[i]);
	    
	    //  cout << " saved as " << title + newProcessName[process]+ "__" + systName[syst] + "__" + modeSyst[i] << endl;
	 
	    
	  }
	
	  
  
  
        }
      }
    }
  } 
  
  
  for (int mode =0; mode < 3; mode ++){
    for (int syst = 0; syst < nSyst; syst++){
      for (int process = 0; process < nProcess; process++){
        for (int i =0; i < 2; i++){
          
	  char modeName[300];
          if (mode == 0) sprintf(modeName,"emu");
          else if (mode == 1) sprintf(modeName,"mumu");
          else sprintf(modeName,"ee");
	  
	  char title[300];
	  char oldtitle[300];
	  
	  if (syst == 0){
	    sprintf(oldtitle,"%s2j1t__",modeName);
	    sprintf(title,"%s2j1t__",modeName);
	    
	    //  cout <<oldtitle + processName[process] ;
	    
	    TH1F* h = (TH1F*) _file0->Get(oldtitle + processName[process]);
	    if (process == nProcess -1){
	      cout <<oldtitle << "  ";
	      TH1F* h1 = (TH1F*) _file0->Get(oldtitle + processName[process+1]);
	      cout << "check: " << h->Integral() <<" + " << h1->Integral()<< h->Integral() + h1->Integral();
	      h->Add(h1);
	      cout << " --> " << h->Integral() << endl; 
	    }
	    TH1F* histo_2j1t = h->Clone();
	    histo_2j1t->SetName(title + newProcessName[process]);
	    //    cout << " saved as " << title + newProcessName[process] << endl;
	    i ++;
          } else if (process != 0){
	  
	    sprintf(oldtitle,"%s2j1t__",modeName);
	    sprintf(title,"%s2j1t__",modeName);
	    
	    //   cout <<oldtitle + processName[process] + "__" + systName[syst] + "__" + modeSyst[i] ;
	    
	    TH1F* h = (TH1F*) _file0->Get(oldtitle + processName[process] + "__" + systName[syst] + "__" + modeSyst[i]);
	   
	    if (process == nProcess -1){
	      cout <<oldtitle  + systName[syst] + "__" + modeSyst[i]<< "  ";
	      TH1F* h1 = (TH1F*) _file0->Get(oldtitle + processName[process+1]+ "__" + systName[syst] + "__" + modeSyst[i]);
	      cout << "check: " << h->Integral() <<" + " << h1->Integral()<< h->Integral() + h1->Integral();
	      h->Add(h1);
	      cout << " --> " << h->Integral() << endl; 
	    }
	    TH1F* histo_2j1t = h->Clone();
	    histo_2j1t->SetName(title + newProcessName[process]+ "__" + systName[syst] + "__" + modeSyst[i]);
	    
	    //    cout << " saved as " << title + newProcessName[process]+ "__" + systName[syst] + "__" + modeSyst[i] << endl;
	 
	    
	  }
	
	  
  
  
        }
      }
    }
  } 
  
  
  
  for (int mode =0; mode < 3; mode ++){
    for (int syst = 0; syst < nSyst; syst++){
      for (int process = 0; process < nProcess; process++){
        for (int i =0; i < 2; i++){
          
	  char modeName[300];
          if (mode == 0) sprintf(modeName,"emu");
          else if (mode == 1) sprintf(modeName,"mumu");
          else sprintf(modeName,"ee");
	  
	  char title[300];
	  char oldtitle[300];
	  
	  if (syst == 0){
	    sprintf(oldtitle,"%s2j2t__",modeName);
	    sprintf(title,"%s2j2t__",modeName);
	    
	    // cout <<oldtitle + processName[process] ;
	    
	    TH1F* h = (TH1F*) _file0->Get(oldtitle + processName[process]);
	    if (process == nProcess -1){
	      cout <<oldtitle << "  ";
	      TH1F* h1 = (TH1F*) _file0->Get(oldtitle + processName[process+1]);
	      cout << "check: " << h->Integral() <<" + " << h1->Integral()<< h->Integral() + h1->Integral();
	      h->Add(h1);
	      cout << " --> " << h->Integral() << endl; 
	    }
	    TH1F* histo_2j2t = h->Clone();
	    histo_2j2t->SetName(title + newProcessName[process]);
	    //   cout << " saved as " << title + newProcessName[process] << endl;
	    i ++;
          } else if (process != 0){
	  
	    sprintf(oldtitle,"%s2j2t__",modeName);
	    sprintf(title,"%s2j2t__",modeName);
	    
	    //   cout <<oldtitle + processName[process] + "__" + systName[syst] + "__" + modeSyst[i] ;
	    
	    TH1F* h = (TH1F*) _file0->Get(oldtitle + processName[process] + "__" + systName[syst] + "__" + modeSyst[i]);
	   
	    if (process == nProcess -1){
	      cout <<oldtitle  + systName[syst] + "__" + modeSyst[i] << "  ";
	      TH1F* h1 = (TH1F*) _file0->Get(oldtitle + processName[process+1]+ "__" + systName[syst] + "__" + modeSyst[i]);
	      cout << "---> check: " << h->Integral() <<" + " << h1->Integral()<< h->Integral() + h1->Integral();
	      h->Add(h1);
	      cout << " --> " << h->Integral() << endl; 
	    }
	    TH1F* histo_2j2t = h->Clone();
	    histo_2j2t->SetName(title + newProcessName[process]+ "__" + systName[syst] + "__" + modeSyst[i]);
	    
	    //   cout << " saved as " << title + newProcessName[process]+ "__" + systName[syst] + "__" + modeSyst[i] << endl;
	 
	    
	  }
	
	  
  
  
        }
      }
    }
  } 
  
  
  f_var.Write();
  f_var.Close();

}
