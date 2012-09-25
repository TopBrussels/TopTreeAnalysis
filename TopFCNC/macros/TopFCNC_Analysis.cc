#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
/*
// Root headers
#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPaveStats.h"
#include "TRandom3.h"
*/
//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../../Content/interface/Dataset.h"
#include "../../Tools/interface/MultiSamplePlot.h"
#include "../interface/TopFCNC_Evt.h"
#include "../interface/TopFCNC_KinFit.h"

#include "../../macros/Style.C"

using namespace std;
using namespace TopTree;

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

int main (int argc, char *argv[])
{

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FCNC analysis ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
  //setMyStyle();


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Output ROOT file
  string postfix = "_Analysis";
  string channelpostfix = "_DiMuonTrigger";
  string comments = "_Run2011AB";
  string rootFileName ("TopFCNC"+postfix+channelpostfix+comments+".root");

  Float_t Luminosity = -1.;
  Float_t EventWeight = 1.;

  vector<Dataset*> dataSets; // needed for MSPlots
  vector<string>   inputFiles;
  
  inputFiles.push_back("../rootfiles/TopFCNC_EventSelection_DiMu_Run2011AB_TTree_Data.root");
  inputFiles.push_back("../rootfiles/TopFCNC_EventSelection_DiMu_Run2011AB_TTree_ST_tbar_tWch_DR.root");
  inputFiles.push_back("../rootfiles/TopFCNC_EventSelection_DiMu_Run2011AB_TTree_ST_t_tWch_DR.root");
  inputFiles.push_back("../rootfiles/TopFCNC_EventSelection_DiMu_Run2011AB_TTree_ww.root");
  inputFiles.push_back("../rootfiles/TopFCNC_EventSelection_DiMu_Run2011AB_TTree_wz.root");
  inputFiles.push_back("../rootfiles/TopFCNC_EventSelection_DiMu_Run2011AB_TTree_zz.root");
  inputFiles.push_back("../rootfiles/TopFCNC_EventSelection_DiMu_Run2011AB_TTree_TTjets.root");
  inputFiles.push_back("../rootfiles/TopFCNC_EventSelection_DiMu_Run2011AB_TTree_Z_Jets.root");
  inputFiles.push_back("../rootfiles/TopFCNC_EventSelection_DiMu_Run2011AB_TTree_ttbar_fcnc.root");

  for(unsigned int iDataSet=0; iDataSet<inputFiles.size(); iDataSet++)
  {
    TFile* inputFile = new TFile(inputFiles[iDataSet].c_str(),"READ");
    TTree* inConfigTree = (TTree*) inputFile->Get("configTree");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    dataSets.push_back( (Dataset*) dataSet->Clone() );
    if(dataSet->Name().find("Data") == 0 || dataSet->Name().find("data") == 0 || dataSet->Name().find("DATA") == 0 )
      Luminosity = dataSet->EquivalentLumi();
    delete tc_dataset;
    inputFile->Close();
    delete inputFile;
  }
  
  if( Luminosity<0 ) Luminosity = 5000.;
  cout<<"Executing analysis for an integrated luminosity of " << Luminosity << " pb^-1" << endl;
    

  TFile      *fout  = new TFile (rootFileName.c_str(), "RECREATE");
  TDirectory *myDir = 0;

  ResolutionFit *resFitLeptons = new ResolutionFit("Muon");
  resFitLeptons->LoadResolutions("../rootfiles/leptonReso_FromZJets.root");

  ResolutionFit *resFitLightJets = new ResolutionFit("LightJet");
  resFitLightJets->LoadResolutions("../rootfiles/lightJetReso.root");
  
  ResolutionFit *resFitBJets = new ResolutionFit("BJet");
  resFitBJets->LoadResolutions("../rootfiles/bJetReso.root");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ///////////////////// MS plots /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  MSPlot["KinFit_Prob"]        = new MultiSamplePlot(dataSets,"KinFit_Prob",100,0,1,"Prob.");
  MSPlot["KinFit_Chi2"]        = new MultiSamplePlot(dataSets,"KinFit_Chi2",500,0,100,"#chi^{2}");
  MSPlot["KinFit_ReducedChi2"] = new MultiSamplePlot(dataSets,"KinFit_ReducedChi2",100,0,1,"#chi^{2}/Ndf");

  MSPlot["KinFit_HadWMass"]    = new MultiSamplePlot(dataSets,"KinFit_HadWMass",120,50,110,"m_{W} [Gev/c^{2}]");
  MSPlot["KinFit_HadTopMass"]  = new MultiSamplePlot(dataSets,"KinFit_HadTopMass",280,100,240,"m^{SM}_{top} [Gev/c^{2}]");
  MSPlot["KinFit_FcncTopMass"] = new MultiSamplePlot(dataSets,"KinFit_FcncTopMass",280,100,240,"m^{FCNC}_{top} [Gev/c^{2}]");

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "TopFCNC"+postfix+channelpostfix+comments;
  pathPNG += "_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////// Analysis //////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(unsigned int iDataSet=0; iDataSet<inputFiles.size(); iDataSet++)
  {
    string dataSetName = dataSets[iDataSet]->Name();
    cout << "Processing DataSet: " << dataSetName << endl;

    TFile* inFile = new TFile(inputFiles[iDataSet].c_str(),"READ");
    
    TTree* inTree = (TTree*) inFile->Get("Tree");
    TBranch* m_br = (TBranch*) inTree->GetBranch("TheTopFCNC_Evt");
    
    int nEvent = 100;//inTree->GetEntries();
    
    TTree* inConfigTree = (TTree*) inFile->Get("configTree");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    
    TopFCNC_Evt* topFCNC_Evt = 0;
    m_br->SetAddress(&topFCNC_Evt);
    m_br->SetAutoDelete(kTRUE);
    
    TopFCNC_KinFit *topFCNC_KinFit = new TopFCNC_KinFit(dataSet, resFitLeptons, resFitLightJets, resFitBJets);
    topFCNC_KinFit->SetMaxNbIter(60);
    topFCNC_KinFit->SetVerbosity(false);
    topFCNC_KinFit->SetFitVerbosity(false);

    for(int iEvt=0; iEvt<nEvent; iEvt++)
    {
      inTree->GetEvent(iEvt);
//      cout << "event: " << iEvt << endl;
      if(iEvt%1000 == 0)
		    std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock()-start)/CLOCKS_PER_SEC << " ("<<100*(iEvt)/(nEvent)<<"%)"<<flush<<"\r";

      if(!topFCNC_Evt->isDiLeptonic()) continue;
      
      topFCNC_KinFit->FitEvent(topFCNC_Evt);//, float WMass = 80.4, float Zmass = 91.2, float topMass = 172.5/173.3 MC vs data);

      topFCNC_Evt->ReconstructEvt();

      EventWeight = Luminosity*topFCNC_Evt->eventWeight();
      
      MSPlot["KinFit_Prob"]->Fill(topFCNC_KinFit->GetProb(), dataSet, true, EventWeight);
      MSPlot["KinFit_Chi2"]->Fill(topFCNC_KinFit->GetChi2(), dataSet, true, EventWeight);
      if(topFCNC_KinFit->GetNdof()!=0)
        MSPlot["KinFit_ReducedChi2"]->Fill(topFCNC_KinFit->GetChi2()/topFCNC_KinFit->GetNdof(), dataSet, true, EventWeight);

      MSPlot["KinFit_HadWMass"]   ->Fill(topFCNC_Evt->W().M(), dataSet, true, EventWeight);
      MSPlot["KinFit_HadTopMass"] ->Fill(topFCNC_Evt->smDecayTop().M(), dataSet, true, EventWeight);
      MSPlot["KinFit_FcncTopMass"]->Fill(topFCNC_Evt->fcncDecayTop().M(), dataSet, true, EventWeight);

    } // loop on events
  } // loop on datasets

  fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
	  MultiSamplePlot *temp = it->second;
	  //temp->addText("CMS preliminary");
	  string name = it->first;
	  temp->Draw(false, name, true, true, true, true, true,1,true); // merge TT/QCD/W/Z/ST/
	  //Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
	  temp->Write(fout, name, true, pathPNG, "pdf");
  }
  
  //delete  
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

