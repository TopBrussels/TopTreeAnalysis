#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// Root headers
//#include "XXX.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"

#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/BtagScaleFactor.h"

#include "../interface/TopFCNC_Evt.h"
#include "../interface/TopFCNC_KinFit.h"

//#include "TopTreeAnalysisBase/macros/Style.C"
#include "Style.C"

using namespace std;
using namespace TopTree;

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

int main (int argc, char *argv[])
{

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FCNC analysis ! " << endl;
  cout << "*************************************************************" << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Choose leptonic channel
  bool UseMuChannel = true;
  TopFCNC_Evt::LeptonType leptChannel = (UseMuChannel ? TopFCNC_Evt::kMuon : TopFCNC_Evt::kElec);
  
  bool DEBUG = false;
  bool IsData = false;

  //Output ROOT file
  string postfix = "_KinFitTreeMaker";
  string channelpostfix = (UseMuChannel ? "_DiMuTrigger" : "_DiElecTrigger");
//  string comments = "_Run2012ABCD_Skim2Mu";
  string comments = "_Run2012ABCD_Skim2Mu3Jets";
  string resoprefix = "Stijns_";
  //string resoprefix = "TTbar_FCNC_";
  //string postfix = "_TopGenEvt";

  string treepath = "$HOME/AnalysisCode/CMSSW_53X/TopBrussels/TopTreeAnalysis/TopFCNC/macros/TopFCNC_EventSelection"+channelpostfix+comments+"_Trees/";
  //string treepath = "$HOME/CMSSW_5_3_3_patch3/src/TopBrussels/TopTreeAnalysis/TopFCNC/macros/TopFCNC_EventSelection_DiMuTrigger_Run2012A_Trees/";
  string resopath = "$HOME/AnalysisCode/CMSSW_53X/TopBrussels/TopTreeAnalysis/TopFCNC/macros/ResolutionFiles/";
  //string resopath = "$HOME/CMSSW_5_3_3_patch3/src/TopBrussels/TopTreeAnalysis/TopFCNC/macros/ResolutionFiles/";

  //vector<Dataset*> dataSets; // needed for MSPlots
  string bckgdName("Z_4Jets");
  string inputFile(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_"+bckgdName+".root");
  string rootFileName("TopFCNC_Analysis"+channelpostfix+comments+"_TTree_"+bckgdName+".root");

/*
  "Data"
  "ST_TBarToDilepton_tW-ch"
  "ST_TBarToThadWlep_tW-ch"
  "ST_TBarToTlepWhad_tW-ch"
  "ST_TToDilepton_tW-ch"
  "ST_TToThadWlep_tW-ch"
  "ST_TToTlepWhad_tW-ch"
  "ST_T_t-ch"
  "ST_Tbar_t-ch"
  "TTJetsTocZbW"
  "TT_FullLeptMGDecays"
  "TT_HadronicMGDecays"
  "TT_SemiLeptMGDecays"
  "WW_To2L2Nu"
  "WZ_To2L2Q"
  "WZ_To3LNu"
  "W_1Jets"
  "W_2Jets"
  "W_3Jets"
  "W_4Jets"
  "ZZ_To2L2Nu"
  "ZZ_To2L2Q"
  "ZZ_To4L"
  "Z_M-10To50"
  "Z_M-50"
  "Z_1Jets"
  "Z_2Jets"
  "Z_3Jets"
  "Z_4Jets"
*/

  ResolutionFit *resFitLeptons = 0;
  if(UseMuChannel){
    resFitLeptons = new ResolutionFit("Muon");
    resFitLeptons->LoadResolutions(resopath+resoprefix+"muonReso.root");
  }
  else{
    resFitLeptons = new ResolutionFit("Electron");
    resFitLeptons->LoadResolutions(resopath+resoprefix+"electronReso.root");
  }
  
  ResolutionFit *resFitBJets = new ResolutionFit("BJet");
  resFitBJets->LoadResolutions(resopath+resoprefix+"bJetReso.root");

  ResolutionFit *resFitQJets = new ResolutionFit("QJet");
  resFitQJets->LoadResolutions(resopath+resoprefix+"qJetReso.root");

  ResolutionFit *resFitLightJets = new ResolutionFit("LightJet");
  resFitLightJets->LoadResolutions(resopath+resoprefix+"lightJetReso.root");
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////// Histograms //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  histo2D["KinFit_HadWMass_"+bckgdName]    = new TH2F(("KinFit_HadWMass_"+bckgdName).c_str(),";m_{jj} [GeV/c^2];#chi^{2}/Ndf",120,50,110,200,0,10);
  histo2D["KinFit_LepZMass_"+bckgdName]    = new TH2F(("KinFit_LepZMass_"+bckgdName).c_str(),";m_{ll} [GeV/c^2];#chi^{2}/Ndf",160,50,130,200,0,10);

  histo2D["KinFit_HadTopMass_"+bckgdName]  = new TH2F(("KinFit_HadTopMass_"+bckgdName).c_str(),";m_{jjj} [GeV/c^2];#chi^{2}/Ndf",280,100,240,200,0,10);
  histo2D["KinFit_FcncTopMass_"+bckgdName] = new TH2F(("KinFit_FcncTopMass_"+bckgdName).c_str(),";m_{llj} [GeV/c^2];#chi^{2}/Ndf",280,100,240,200,0,10);
  
  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "TopFCNC"+postfix+channelpostfix+comments;
  pathPNG += "_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  ////////////////////////////////////////////////////////////////////
  //////////////////// Counter ///////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  double Counter = 0;

  ////////////////////////////////////////////////////////////////////
  //////////////////// TTrees ////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  TFile* inFile = TFile::Open(inputFile.c_str());
  
  // Config tree
  TTree* inConfigTree = (TTree*) inFile->Get("configTree");
  TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");

  TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
  d_br->SetAddress(&tc_dataset);
  
  inConfigTree->GetEvent(0);
  Dataset* dataSet = (Dataset*) tc_dataset->At(0);

  string dataSetName = dataSet->Name();
  
  float EventWeight = 1.;
  float NormFactor = dataSet->NormFactor();
  float Luminosity = -1.;
  
  cout << "Processing DataSet: " << dataSetName << endl;
  
  if(dataSetName.find("Data") == string::npos && dataSetName.find("data") == string::npos && dataSetName.find("DATA") == string::npos )
    IsData = false;
  else
    IsData = true;
  
  if(IsData)
    Luminosity = dataSet->EquivalentLumi();
  else
    Luminosity = 19626;

  cout<<"Executing analysis for an integrated luminosity of " << Luminosity << "/pb" << endl;
  
  // Input tree
  TTree* inTree = (TTree*)   inFile->Get("Tree");
  TBranch* m_br = (TBranch*) inTree->GetBranch("TheTopFCNC_Evt");
  
  TopFCNC_Evt* topFCNC_Evt = 0;
  m_br->SetAddress(&topFCNC_Evt);
  m_br->SetAutoDelete(kTRUE);
  
  // Output trees
  TFile *outFile  = new TFile (rootFileName.c_str(), "RECREATE");
  TTree *oTree = inTree->CloneTree(0);
  TTree *oConfigTree = inConfigTree->CloneTree();
  
  ////////////////////////////////////////////////////////////////////
  ////////////////// B-tagging SF ////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  string value = "mean"; // "mean","max","min"
  string tagger = "CSV";
  string period = "ABCD";
  vector<string> bckgdName_vec;
  bckgdName_vec.push_back(bckgdName);

  string btagfilename("TopFCNC_BtaggingEff_DiMuTrigger_Run2012ABCD_Doc/TopFCNC_BtaggingEff_DiMuTrigger"+comments+".root");
  BtagScaleFactor* myBtagScaleFactor = new BtagScaleFactor(btagfilename, bckgdName_vec, value, tagger, period);

  double SFb_CSVL = 0;
  double SFb_CSVM = 0;
  double SFb_CSVT = 0;
  
  // Add branches to the output tree
  oTree->Branch("SFb_CSVL",&SFb_CSVL,"SFb_CSVL/D");
  oTree->Branch("SFb_CSVM",&SFb_CSVM,"SFb_CSVM/D");
  oTree->Branch("SFb_CSVT",&SFb_CSVT,"SFb_CSVT/D");

  ////////////////////////////////////////////////////////////////////
  /////////////////// KinFitter //////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  double topMass = 172.5;
  if(IsData)
    topMass = 173.3;
  // Top quark mass LHC average = 172.6 GeV/c²
  
  TopFCNC_KinFit *topFCNC_KinFit = new TopFCNC_KinFit(dataSet, resFitLeptons, resFitBJets, resFitQJets, resFitLightJets,80.4,91.2,topMass);
  
  // Add constraints to the kinfitter
  vector<string> constraints;
  constraints.push_back("kHadWMass");
  constraints.push_back("kLepZMass");
  constraints.push_back("kHadTopMass");
  constraints.push_back("kFcncTopMass");
  //constraints.push_back("kEqualTopMasses");
  //constraints.push_back("kSumPx");
  //constraints.push_back("kSumPy");
  
  topFCNC_KinFit->SetConstraints(constraints);
  
  topFCNC_KinFit->SetMaxNbIter(50);
  topFCNC_KinFit->SetMaxDeltaS(0.00005);
  topFCNC_KinFit->SetMaxF(0.0001);
  
  topFCNC_KinFit->SetVerbosity(false);
  topFCNC_KinFit->SetFitVerbosity(0);
  
  double kin_prob        = -1.;
  double kin_chi2        = -1.;
  double kin_ndf         = -1.;
  double kin_chi2ByNdf   = -1.;
  double kin_hadWmass    = -1.;
  double kin_lepZmass    = -1.;
  double kin_hadtopmass  = -1.;
  double kin_fcnctopmass = -1.;

  // Add branches to the output tree
  oTree->Branch("kin_prob",&kin_prob,"kin_prob/D");
  oTree->Branch("kin_chi2",&kin_chi2,"kin_chi2/D");
  oTree->Branch("kin_ndf",&kin_ndf,"kin_ndf/D");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////// Analysis ////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //int nEvent = inTree->GetEntries();
  int nEvent = 1000;
  cout<< " - number of entries: "<<nEvent<<endl;
  cout<< " - number of events: "<<NormFactor*Luminosity*nEvent<<endl;
    
  for(int iEvt=0; iEvt<nEvent; iEvt++)
  {
    inTree->GetEvent(iEvt);
    if(iEvt%1000 == 0)
      std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock()-start)/CLOCKS_PER_SEC << " ("<<100*(iEvt)/(nEvent)<<"%)"<<flush<<"\r";
    
    if(!topFCNC_Evt->isDiLeptonic(leptChannel)) continue;

    topFCNC_KinFit->FitEvent(topFCNC_Evt);
    
    topFCNC_Evt->ReconstructEvt();
    
    EventWeight = topFCNC_Evt->eventWeight();
    
    Counter += NormFactor*Luminosity*EventWeight;
    
    kin_prob        = topFCNC_KinFit->GetProb();
    kin_chi2        = topFCNC_KinFit->GetChi2();
    kin_ndf         = topFCNC_KinFit->GetNdof();
    if(kin_ndf!=0)
      kin_chi2ByNdf = kin_chi2/topFCNC_KinFit->GetNdof();
    else
      kin_chi2ByNdf = 999;
    
    kin_hadWmass    = topFCNC_Evt->W().M();
    kin_lepZmass    = topFCNC_Evt->Z().M();
    kin_hadtopmass  = topFCNC_Evt->smDecayTop().M();
    kin_fcnctopmass = topFCNC_Evt->fcncDecayTop().M();
    
    histo2D[("KinFit_HadWMass_"+dataSet->Name()).c_str()]->Fill(kin_hadWmass,kin_chi2ByNdf);
    histo2D[("KinFit_LepZMass_"+dataSet->Name()).c_str()]->Fill(kin_lepZmass,kin_chi2ByNdf);

    histo2D[("KinFit_HadTopMass_"+dataSet->Name()).c_str()]->Fill(kin_hadtopmass,kin_chi2ByNdf);
    histo2D[("KinFit_FcncTopMass_"+dataSet->Name()).c_str()]->Fill(kin_fcnctopmass,kin_chi2ByNdf);
    
    vector<TRootJet> selectedJets = topFCNC_Evt->selectedJets();

    SFb_CSVL = (IsData ? 1 : myBtagScaleFactor->GetBtagSF(0, selectedJets, "L", 0.244));
    SFb_CSVM = (IsData ? 1 : myBtagScaleFactor->GetBtagSF(0, selectedJets, "M", 0.679));
    SFb_CSVT = (IsData ? 1 : myBtagScaleFactor->GetBtagSF(0, selectedJets, "T", 0.898));
    
    if(DEBUG){
      cout << fixed << setprecision(3);
      cout << "**********************************************" << endl;
      cout << "Btag weight (CSVL) = " << SFb_CSVL << endl;
      cout << "Btag weight (CSVM) = " << SFb_CSVM << endl;
      cout << "Btag weight (CSVT) = " << SFb_CSVT << endl;
    }
    oTree->Fill();
  } // loop on events
  //oTree->Print();
  delete topFCNC_KinFit;
  
  cout << "/************** Number of selected events ***************/" <<endl;
  cout << " - Dataset: " << dataSetName << endl;
  cout << " - inclusive: " << Counter << endl;
  cout << "/********************************************************/" <<endl;

  //MultiSample plots
  outFile->cd();
  oConfigTree->Write();
  oTree->Write();

  //TH1X histograms
  TDirectory* th1dir = outFile->mkdir("Histos1D");
  th1dir->cd();
  for(map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
    TH1F *temp = it->second;
    temp->Write();
  }
  //TH2X histograms
  TDirectory* th2dir = outFile->mkdir("Histos2D");
  th2dir->cd();
  for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {
    TH2F *temp = it->second;
    temp->Write();
  }
  
  //delete
  delete inFile;
  delete outFile;
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}
