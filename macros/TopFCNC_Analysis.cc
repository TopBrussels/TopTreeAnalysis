#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// Root headers
#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include  "TGraphErrors.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPaveStats.h"
#include "TRandom3.h"

// RooFit headers
#include "RooAddition.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooCatType.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormula.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"

#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"

#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"

#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooTable.h"
#include "Roo1DTable.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ModelConfig.h"

#include "RooStats/HypoTestInverterOriginal.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/HybridCalculatorOriginal.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
//#include "../Selection/interface/SelectionTable.h"
//#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
//#include "../Tools/interface/JetTools.h"
//#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
//#include "../Tools/interface/TTreeLoader.h"
//#include "../Tools/interface/AnalysisEnvironmentLoader.h"
//#include "../Reconstruction/interface/JetCorrectorParameters.h"
//#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
//#include "../Reconstruction/interface/MakeBinning.h"
//#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../TopFCNC/interface/TopFCNC_Evt.h"
#include "../TopFCNC/interface/TopFCNC_KinFit.h"

#include "Style.C"

using namespace RooFit;
using namespace RooStats ;
using namespace std;
//using namespace TopTree;

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
  string channelpostfix = "_diMuon";
  string postfix = "_Analysis";
  string rootFileName ("TopFCNC"+postfix+channelpostfix+".root");

  Float_t Luminosity = -1.;

  vector<Dataset*> dataSets; // needed for MSPlots
  vector<string>   inputFiles;
  
  //inputFiles.push_back("TopFCNC_EventSelection_DiMu_TTree_data.root");
  inputFiles.push_back("TopFCNC_EventSelection_DiMu_TTree_ttjets.root");
  
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

  ResolutionFit *resFitLightJets = new ResolutionFit("LightJet");
  resFitLightJets->LoadResolutions("lightJetReso.root");
  
  ResolutionFit *resFitBJets = new ResolutionFit("BJet");
  resFitBJets->LoadResolutions("bJetReso.root");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ///////////////////// MS plots /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  MSPlot["KinFit_Prob"] = new MultiSamplePlot(dataSets,"KinFit_Prob",100,0,1,"Prob.");
  MSPlot["KinFit_Chi2"] = new MultiSamplePlot(dataSets,"KinFit_Chi2",100,0,1,"#chi^{2}");

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "TopFCNC_Analysis"+postfix+channelpostfix;
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
    
    int nEvent = inTree->GetEntries();
    
    TTree* inConfigTree = (TTree*) inFile->Get("configTree");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    
    TopFCNC_Evt* topFCNC_Evt = 0;
    m_br->SetAddress(&topFCNC_Evt);
    
    TopFCNC_KinFit *topFCNC_KinFit = new TopFCNC_KinFit(dataSet, resFitLightJets, resFitBJets);
    topFCNC_KinFit->SetVerbosity(false);

    for(int iEvt=0; iEvt<nEvent; iEvt++)
    {
      inTree->GetEvent(iEvt);
//      cout << "event: " << iEvt << endl;
      if(iEvt%1000 == 0)
        std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC <<endl;

      if(!topFCNC_Evt->isDiLeptonic()) continue;
      
      topFCNC_KinFit->FitEvent(topFCNC_Evt);//, float WMass = 80.4, float Zmass = 91.2, float topMass = 172.5);

      MSPlot["KinFit_Prob"]->Fill(topFCNC_KinFit->GetProb(), dataSet, true, Luminosity*topFCNC_Evt->eventWeight());
      MSPlot["KinFit_Chi2"]->Fill(topFCNC_KinFit->GetChi2(), dataSet, true, Luminosity*topFCNC_Evt->eventWeight());
    }
  }

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

