#include "TFile.h"
#include "TH1D.h"

#include <iostream>

using namespace std;

void MakePileUpRootFiles()
{
  
  const int NbOfFiles = 4;
  string comment = "";
  //string comment = "_Systematic_Up_5perc";
  //string comment = "_Systematic_Down_5perc";
  string filenames[NbOfFiles] = {
    "DoubleMu_Run2012A_22Jan2013_TopTreeID_3495_1372688083_PileupHistogram"+comment+".root",
    "DoubleMu_Run2012B_22Jan2013_TopTreeID_3494_1372688071_PileupHistogram"+comment+".root",
    "DoubleMu_Run2012C_22Jan2013_TopTreeID_3493_1372689333_PileupHistogram"+comment+".root",
    "DoubleMu_Run2012D_22Jan2013_TopTreeID_3490_1372689340_PileupHistogram"+comment+".root"
    /*
    "MuEG_Run2012A_22Jan2013_TopTreeID_3497_1372750900_PileupHistogram"+comment+".root",
    "MuEG_Run2012B_22Jan2013_TopTreeID_3496_1372687900_PileupHistogram"+comment+".root",
    "MuEG_Run2012C_22Jan2013_TopTreeID_3498_1372751025_PileupHistogram"+comment+".root",
    "MuEG_Run2012D_22Jan2013_TopTreeID_3499_1372751043_PileupHistogram"+comment+".root"
    */
  };
  const double IntLumi[NbOfFiles] = {
    866.6,//863.8,
    4380.7,//4375,
    7102,//7146,
    7277//7292
  };

  TFile *files[NbOfFiles];
  TH1D* pileup = 0;
  
  for(int i=0;i<NbOfFiles;i++){
    files[i] = TFile::Open(filenames[i].c_str());
    TH1D* tmp = (TH1D*)files[i]->Get("pileup");
    if(i==0){
      pileup = (TH1D*)tmp->Clone();
      pileup->Scale(IntLumi[i]/pileup->Integral());
      cout << "Integral = " << pileup->Integral() << endl;
    }
    else{
      pileup->Add(tmp,IntLumi[i]/tmp->Integral());
    }
  }
  cout << "Integral = " << pileup->Integral() << endl;
  pileup->Draw();
  
  TFile *fout = new TFile("fout.root","RECREATE");
  fout->cd();
  pileup->Write();
  fout->Close();
}
