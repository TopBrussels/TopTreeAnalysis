#include "TFile.h"
#include "TH1D.h"

#include <iostream>

using namespace std;

void MakePileUpRootFiles()
{
  
  const int NbOfFiles = 3;
  string filenames[NbOfFiles] = {
    "DoubleMu_Run2012A_13Jul_TopTreeID_2219_PileupHistogram.root",
    "DoubleMu_Run2012A_06Aug_TopTreeID_2204_PileupHistogram.root",
    "DoubleMu_Run2012B_13Jul_TopTreeID_2235_PileupHistogram.root"//,
    //"DoubleMu_Run2012C_24Aug_TopTreeID_2203_PileupHistogram.root",
    //"DoubleMu_Run2012C_Promptv2_TopTreeID_2226_PileupHistogram.root",
    //"DoubleMu_Run2012D_Promptv1_TopTreeID_2218_PileupHistogram.root"
  };
  const double IntLumi[NbOfFiles] = {
    800.342,
    80.074,
    4315//,
    //495.003,
    //6392,
    //6767
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
