#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
// Root headers
#include "TEfficiency.h"

/*
// Root headers
#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TEfficiency.h"
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

  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
  //setMyStyle();


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
  
  Float_t Luminosity = -1.;
  Float_t NormFactor = 1.;
  Float_t EventWeight = 1.;

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

  if( Luminosity<0 ) Luminosity = 19626;
  cout<<"Executing analysis for an integrated luminosity of " << Luminosity << " pb^-1" << endl;
  
  TFile      *fout  = new TFile (rootFileName.c_str(), "RECREATE");
  TDirectory *myDir = 0;

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

  ////////////////////////////////////////////////////////////////////
  ///////////////////// MS plots /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
/*
  MSPlot["KinFit_NbOfBtaggedJets_CSVL"] = new MultiSamplePlot(dataSets,"KinFit_NbOfBtaggedJets_CSVL",10,0,10,"Nb. of b-tagged jets (CSVL)");
  MSPlot["KinFit_NbOfBtaggedJets_CSVM"] = new MultiSamplePlot(dataSets,"KinFit_NbOfBtaggedJets_CSVM",10,0,10,"Nb. of b-tagged jets (CSVM)");
  MSPlot["KinFit_NbOfBtaggedJets_CSVT"] = new MultiSamplePlot(dataSets,"KinFit_NbOfBtaggedJets_CSVT",10,0,10,"Nb. of b-tagged jets (CSVT)");
  MSPlot["KinFit_BtaggingSF_CSVL"] = new MultiSamplePlot(dataSets,"KinFit_BtaggingSF_CSVL",80,0,4,"B-tagging SF");
  MSPlot["KinFit_BtaggingSF_CSVM"] = new MultiSamplePlot(dataSets,"KinFit_BtaggingSF_CSVM",80,0,4,"B-tagging SF");
  MSPlot["KinFit_BtaggingSF_CSVT"] = new MultiSamplePlot(dataSets,"KinFit_BtaggingSF_CSVT",80,0,4,"B-tagging SF");
  
  // NO BTAG CUT________________________________________________________________________________________
  MSPlot["KinFit_Prob"]        = new MultiSamplePlot(dataSets,"KinFit_Prob",100,0,1,"Prob.");
  MSPlot["KinFit_Chi2"]        = new MultiSamplePlot(dataSets,"KinFit_Chi2",500,0,100,"#chi^{2}");
  MSPlot["KinFit_ReducedChi2"] = new MultiSamplePlot(dataSets,"KinFit_ReducedChi2",200,0,10,"#chi^{2}/Ndf");

  MSPlot["KinFit_HadWMass"]    = new MultiSamplePlot(dataSets,"KinFit_HadWMass",120,50,110,"m_{W} [Gev/c^{2}]");
  MSPlot["KinFit_HadTopMass"]  = new MultiSamplePlot(dataSets,"KinFit_HadTopMass",280,100,240,"m^{SM}_{top} [Gev/c^{2}]");
  MSPlot["KinFit_FcncTopMass"] = new MultiSamplePlot(dataSets,"KinFit_FcncTopMass",280,100,240,"m^{FCNC}_{top} [Gev/c^{2}]");

  MSPlot["KinFit_LepZ_Pt"]     = new MultiSamplePlot(dataSets,"KinFit_LepZ_Pt",400,0,200,"p^{ll}_{T} [Gev/c]");
  MSPlot["KinFit_Lep_DR"]      = new MultiSamplePlot(dataSets,"KinFit_Lep_DR",50,0,5,"#Delta R(l^{+}l^{-})");
  MSPlot["KinFit_Met"]         = new MultiSamplePlot(dataSets,"KinFit_Met",300,0,150,"\\slashE_{T} [Gev/c]");
  MSPlot["KinFit_Njets"]       = new MultiSamplePlot(dataSets,"KinFit_Njets",8,4,12,"Nb. of jets");
  MSPlot["KinFit_LepZ_Pt_Over_JetSyst_Pt"] = new MultiSamplePlot(dataSets, "KinFit_LepZ_Pt_Over_JetSyst_Pt", 500, 0, 10, "p_{T}^{ll}/p_{T}^{#Sigma jets}");
  MSPlot["KinFit_LepZ_DPhi_JetSyst"]       = new MultiSamplePlot(dataSets, "KinFit_LepZ_DPhi_JetSyst", 70, 0, 3.5, "#Delta #Phi(l^{+}+l^{-},#Sigma jets))");
  // CSV Tight Working point_____________________________________________________________________________
  MSPlot["KinFit_Prob_AtLeast1Btag_CSVT"]        = new MultiSamplePlot(dataSets,"KinFit_Prob_AtLeast1Btag_CSVT",100,0,1,"Prob.");
  MSPlot["KinFit_Chi2_AtLeast1Btag_CSVT"]        = new MultiSamplePlot(dataSets,"KinFit_Chi2_AtLeast1Btag_CSVT",500,0,100,"#chi^{2}");
  MSPlot["KinFit_ReducedChi2_AtLeast1Btag_CSVT"] = new MultiSamplePlot(dataSets,"KinFit_ReducedChi2_AtLeast1Btag_CSVT",200,0,10,"#chi^{2}/Ndf");
  
  MSPlot["KinFit_HadWMass_AtLeast1Btag_CSVT"]    = new MultiSamplePlot(dataSets,"KinFit_HadWMass_AtLeast1Btag_CSVT",120,50,110,"m_{W} [Gev/c^{2}]");
  MSPlot["KinFit_HadTopMass_AtLeast1Btag_CSVT"]  = new MultiSamplePlot(dataSets,"KinFit_HadTopMass_AtLeast1Btag_CSVT",280,100,240,"m^{SM}_{top} [Gev/c^{2}]");
  MSPlot["KinFit_FcncTopMass_AtLeast1Btag_CSVT"] = new MultiSamplePlot(dataSets,"KinFit_FcncTopMass_AtLeast1Btag_CSVT",280,100,240,"m^{FCNC}_{top} [Gev/c^{2}]");
  
  MSPlot["KinFit_LepZ_Pt_AtLeast1Btag_CSVT"]     = new MultiSamplePlot(dataSets,"KinFit_LepZ_Pt_AtLeast1Btag_CSVT",400,0,200,"p^{ll}_{T} [Gev/c]");
  MSPlot["KinFit_Lep_DR_AtLeast1Btag_CSVT"]      = new MultiSamplePlot(dataSets,"KinFit_Lep_DR_AtLeast1Btag_CSVT",50,0,5,"#Delta R(l^{+}l^{-})");
  MSPlot["KinFit_Met_AtLeast1Btag_CSVT"]         = new MultiSamplePlot(dataSets,"KinFit_Met_AtLeast1Btag_CSVT",300,0,150,"\\slashE_{T} [Gev/c]");
  MSPlot["KinFit_Njets_AtLeast1Btag_CSVT"]       = new MultiSamplePlot(dataSets,"KinFit_Njets_AtLeast1Btag_CSVT",8,4,12,"Nb. of jets");
  MSPlot["KinFit_LepZ_Pt_Over_JetSyst_Pt_AtLeast1Btag_CSVT"] = new MultiSamplePlot(dataSets, "KinFit_LepZ_Pt_Over_JetSyst_Pt_AtLeast1Btag_CSVT", 500, 0, 10, "p_{T}^{ll}/p_{T}^{#Sigma jets}");
  MSPlot["KinFit_LepZ_DPhi_JetSyst_AtLeast1Btag_CSVT"]       = new MultiSamplePlot(dataSets, "KinFit_LepZ_DPhi_JetSyst_AtLeast1Btag_CSVT", 70, 0, 3.5, "#Delta #Phi(l^{+}+l^{-},#Sigma jets))");

  // CSV Medium Working point____________________________________________________________________________
  MSPlot["KinFit_Prob_AtLeast1Btag_CSVM"]        = new MultiSamplePlot(dataSets,"KinFit_Prob_AtLeast1Btag_CSVM",100,0,1,"Prob.");
  MSPlot["KinFit_Chi2_AtLeast1Btag_CSVM"]        = new MultiSamplePlot(dataSets,"KinFit_Chi2_AtLeast1Btag_CSVM",500,0,100,"#chi^{2}");
  MSPlot["KinFit_ReducedChi2_AtLeast1Btag_CSVM"] = new MultiSamplePlot(dataSets,"KinFit_ReducedChi2_AtLeast1Btag_CSVM",200,0,10,"#chi^{2}/Ndf");

  MSPlot["KinFit_HadWMass_AtLeast1Btag_CSVM"]    = new MultiSamplePlot(dataSets,"KinFit_HadWMass_AtLeast1Btag_CSVM",120,50,110,"m_{W} [Gev/c^{2}]");
  MSPlot["KinFit_HadTopMass_AtLeast1Btag_CSVM"]  = new MultiSamplePlot(dataSets,"KinFit_HadTopMass_AtLeast1Btag_CSVM",280,100,240,"m^{SM}_{top} [Gev/c^{2}]");
  MSPlot["KinFit_FcncTopMass_AtLeast1Btag_CSVM"] = new MultiSamplePlot(dataSets,"KinFit_FcncTopMass_AtLeast1Btag_CSVM",280,100,240,"m^{FCNC}_{top} [Gev/c^{2}]");

  MSPlot["KinFit_LepZ_Pt_AtLeast1Btag_CSVM"]     = new MultiSamplePlot(dataSets,"KinFit_LepZ_Pt_AtLeast1Btag_CSVM",400,0,200,"p^{ll}_{T} [Gev/c]");
  MSPlot["KinFit_Lep_DR_AtLeast1Btag_CSVM"]      = new MultiSamplePlot(dataSets,"KinFit_Lep_DR_AtLeast1Btag_CSVM",50,0,5,"#Delta R(l^{+}l^{-})");
  MSPlot["KinFit_Met_AtLeast1Btag_CSVM"]         = new MultiSamplePlot(dataSets,"KinFit_Met_AtLeast1Btag_CSVM",300,0,150,"\\slashE_{T} [Gev/c]");
  MSPlot["KinFit_Njets_AtLeast1Btag_CSVM"]       = new MultiSamplePlot(dataSets,"KinFit_Njets_AtLeast1Btag_CSVM",8,4,12,"Nb. of jets");
  MSPlot["KinFit_LepZ_Pt_Over_JetSyst_Pt_AtLeast1Btag_CSVM"] = new MultiSamplePlot(dataSets, "KinFit_LepZ_Pt_Over_JetSyst_Pt_AtLeast1Btag_CSVM", 500, 0, 10, "p_{T}^{ll}/p_{T}^{#Sigma jets}");
  MSPlot["KinFit_LepZ_DPhi_JetSyst_AtLeast1Btag_CSVM"]       = new MultiSamplePlot(dataSets, "KinFit_LepZ_DPhi_JetSyst_AtLeast1Btag_CSVM", 70, 0, 3.5, "#Delta #Phi(l^{+}+l^{-},#Sigma jets))");
  
  // CSV Loose Working point (Veto) _____________________________________________________________________
  MSPlot["KinFit_Prob_NoBtag_CSVL"]        = new MultiSamplePlot(dataSets,"KinFit_Prob_NoBtag_CSVL",100,0,1,"Prob.");
  MSPlot["KinFit_Chi2_NoBtag_CSVL"]        = new MultiSamplePlot(dataSets,"KinFit_Chi2_NoBtag_CSVL",500,0,100,"#chi^{2}");
  MSPlot["KinFit_ReducedChi2_NoBtag_CSVL"] = new MultiSamplePlot(dataSets,"KinFit_ReducedChi2_NoBtag_CSVL",200,0,10,"#chi^{2}/Ndf");

  MSPlot["KinFit_HadWMass_NoBtag_CSVL"]    = new MultiSamplePlot(dataSets,"KinFit_HadWMass_NoBtag_CSVL",120,50,110,"m_{W} [Gev/c^{2}]");
  MSPlot["KinFit_HadTopMass_NoBtag_CSVL"]  = new MultiSamplePlot(dataSets,"KinFit_HadTopMass_NoBtag_CSVL",240,100,240,"m^{SM}_{top} [Gev/c^{2}]");
  MSPlot["KinFit_FcncTopMass_NoBtag_CSVL"] = new MultiSamplePlot(dataSets,"KinFit_FcncTopMass_NoBtag_CSVL",240,100,240,"m^{FCNC}_{top} [Gev/c^{2}]");
  
  MSPlot["KinFit_LepZ_Pt_NoBtag_CSVL"]     = new MultiSamplePlot(dataSets,"KinFit_LepZ_Pt_NoBtag_CSVL",400,0,200,"p^{ll}_{T} [Gev/c]");
  MSPlot["KinFit_Lep_DR_NoBtag_CSVL"]      = new MultiSamplePlot(dataSets,"KinFit_Lep_DR_NoBtag_CSVL",50,0,5,"#Delta R(l^{+}l^{-})");
  MSPlot["KinFit_Met_NoBtag_CSVL"]         = new MultiSamplePlot(dataSets,"KinFit_Met_NoBtag_CSVL",300,0,150,"\\slashE_{T} [Gev/c]");
  MSPlot["KinFit_Njets_NoBtag_CSVL"]       = new MultiSamplePlot(dataSets,"KinFit_Njets_NoBtag_CSVL",8,4,12,"Nb. of jets");
  MSPlot["KinFit_LepZ_Pt_Over_JetSyst_Pt_NoBtag_CSVL"] = new MultiSamplePlot(dataSets, "KinFit_LepZ_Pt_Over_JetSyst_Pt_NoBtag_CSVL", 500, 0, 10, "p_{T}^{ll}/p_{T}^{#Sigma jets}");
  MSPlot["KinFit_LepZ_DPhi_JetSyst_NoBtag_CSVL"]       = new MultiSamplePlot(dataSets, "KinFit_LepZ_DPhi_JetSyst_NoBtag_CSVL", 70, 0, 3.5, "#Delta #Phi(l^{+}+l^{-},#Sigma jets))");
*/

  histo2D["KinFit_HadWMass_"+bckgdName]   = new TH2F(("KinFit_HadWMass_"+bckgdName).c_str(),";W mass [GeV/c^2];#chi^{2}/Ndf",120,50,110,200,0,10);
  histo2D["KinFit_HadTopMass_"+bckgdName] = new TH2F(("KinFit_HadTopMass_"+bckgdName).c_str(),";Top mass [GeV/c^2];#chi^{2}/Ndf",280,100,240,200,0,10);
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
  
  cout << "Processing DataSet: " << dataSetName << endl;
  
  if(dataSetName.find("Data") == string::npos && dataSetName.find("data") == string::npos && dataSetName.find("DATA") == string::npos )
    IsData = false;
  else
    IsData = true;
  
  NormFactor = dataSet->NormFactor();
  
  // Input tree
  TTree* inTree = (TTree*)   inFile->Get("Tree");
  TBranch* m_br = (TBranch*) inTree->GetBranch("TheTopFCNC_Evt");
  
  TopFCNC_Evt* topFCNC_Evt = 0;
  m_br->SetAddress(&topFCNC_Evt);
  m_br->SetAutoDelete(kTRUE);
  
  // Output tree
  TTree *oTree = inTree->CloneTree(0);
  
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
  int nEvent = 200;
  cout<< " - number of entries: "<<nEvent<<endl;
  cout<< " - number of events: "<<NormFactor*Luminosity*nEvent<<endl;
    
/*
  int kin_nbofbtag_csvt = 0;
  int kin_nbofbtag_csvm = 0;
  int kin_nbofbtag_csvl = 0;
*/  
  for(int iEvt=0; iEvt<nEvent; iEvt++)
  {
    inTree->GetEvent(iEvt);
    if(iEvt%1000 == 0)
      std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock()-start)/CLOCKS_PER_SEC << " ("<<100*(iEvt)/(nEvent)<<"%)"<<flush<<"\r";
    
    if(!topFCNC_Evt->isDiLeptonic(leptChannel)) continue;

    topFCNC_KinFit->FitEvent(topFCNC_Evt);
    
//    cout << "After KinFit: Light jet 1: " << topFCNC_Evt->quark1FromW() << endl;
    
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
    kin_hadtopmass  = topFCNC_Evt->smDecayTop().M();
    kin_fcnctopmass = topFCNC_Evt->fcncDecayTop().M();
    
    histo2D[("KinFit_HadWMass_"+dataSet->Name()).c_str()]->Fill(kin_hadWmass,kin_chi2ByNdf);
    histo2D[("KinFit_HadTopMass_"+dataSet->Name()).c_str()]->Fill(kin_hadtopmass,kin_chi2ByNdf);
    
    
    vector<TRootJet> selectedJets = topFCNC_Evt->selectedJets();

    SFb_CSVL = (IsData ? 1 : myBtagScaleFactor->GetBtagSF(0, selectedJets, "L", 0.244));
    SFb_CSVM = (IsData ? 1 : myBtagScaleFactor->GetBtagSF(0, selectedJets, "M", 0.679));
    SFb_CSVT = (IsData ? 1 : myBtagScaleFactor->GetBtagSF(0, selectedJets, "T", 0.898));
    
    if(DEBUG){
      cout << fixed << setprecision(3);
      cout << "**********************************************" << endl;
      cout << "Btag weight (CSVL) = " << SFb_CSVL << endl;//" / Nb of b-tagged jets = " << kin_nbofbtag_csvl << endl;
      cout << "Btag weight (CSVM) = " << SFb_CSVM << endl;//" / Nb of b-tagged jets = " << kin_nbofbtag_csvm << endl;
      cout << "Btag weight (CSVT) = " << SFb_CSVT << endl;//" / Nb of b-tagged jets = " << kin_nbofbtag_csvt << endl;
    }
    oTree->Fill();
/*
    if(kin_nbofbtag_csvl==0) MSPlot["KinFit_BtaggingSF_CSVL"]->Fill(SFb_CSVL, dataSets[iDataSet], true, Luminosity*EventWeight);
    if(kin_nbofbtag_csvm>=1) MSPlot["KinFit_BtaggingSF_CSVM"]->Fill(SFb_CSVM, dataSets[iDataSet], true, Luminosity*EventWeight);
    if(kin_nbofbtag_csvt>=1) MSPlot["KinFit_BtaggingSF_CSVT"]->Fill(SFb_CSVT, dataSets[iDataSet], true, Luminosity*EventWeight);
    
    MSPlot["KinFit_NbOfBtaggedJets_CSVL"]->Fill(kin_nbofbtag_csvl, dataSets[iDataSet], true, Luminosity*EventWeight);
    MSPlot["KinFit_NbOfBtaggedJets_CSVM"]->Fill(kin_nbofbtag_csvm, dataSets[iDataSet], true, Luminosity*EventWeight);
    MSPlot["KinFit_NbOfBtaggedJets_CSVT"]->Fill(kin_nbofbtag_csvt, dataSets[iDataSet], true, Luminosity*EventWeight);
    
    // CSV Medium Working point ___________________________________________________________________________
    if(bdisc>0.679){ // DO NOT FORGET THE B-TAGGING SF WHEN MC !!!!!!
      // CSV Tight Working point ___________________________________________________________________________
      if(bdisc>0.898){ // DO NOT FORGET THE B-TAGGING SF WHEN MC !!!!!!
        Counters[1][iDataSet] += NormFactor*Luminosity*EventWeight*SFb_CSVT;
        selecTableDiMu.Fill(iDataSet,4,EventWeight);
        
        MSPlot["KinFit_Prob_AtLeast1Btag_CSVT"]->Fill(kin_prob, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVT);
        MSPlot["KinFit_Chi2_AtLeast1Btag_CSVT"]->Fill(kin_chi2, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVT);
        
        if(topFCNC_KinFit->GetNdof()!=0)
          MSPlot["KinFit_ReducedChi2_AtLeast1Btag_CSVT"]->Fill(kin_chi2ByNdf, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVT);
        
        MSPlot["KinFit_HadWMass_AtLeast1Btag_CSVT"]   ->Fill(kin_hadWmass, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVT);
        MSPlot["KinFit_HadTopMass_AtLeast1Btag_CSVT"] ->Fill(kin_hadtopmass, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVT);
        MSPlot["KinFit_FcncTopMass_AtLeast1Btag_CSVT"]->Fill(kin_fcnctopmass, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVT);
        MSPlot["KinFit_LepZ_Pt_AtLeast1Btag_CSVT"]    ->Fill(kin_lepZpt, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVT);
        MSPlot["KinFit_Lep_DR_AtLeast1Btag_CSVT"]     ->Fill(kin_lepDR, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVT);
        MSPlot["KinFit_Met_AtLeast1Btag_CSVT"]        ->Fill(kin_met, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVT);
        MSPlot["KinFit_Njets_AtLeast1Btag_CSVT"]->Fill(selectedJets.size(), dataSets[iDataSet], true, Luminosity*EventWeight);
        MSPlot["KinFit_LepZ_Pt_Over_JetSyst_Pt_AtLeast1Btag_CSVT"]->Fill(kin_lepZpt/kin_jetSystpt, dataSets[iDataSet], true, Luminosity*EventWeight);
        MSPlot["KinFit_LepZ_DPhi_JetSyst_AtLeast1Btag_CSVT"]->Fill(kin_lepZdphi_jetSyst, dataSets[iDataSet], true, Luminosity*EventWeight);
      }
      Counters[2][iDataSet] += NormFactor*Luminosity*EventWeight*SFb_CSVM;
      selecTableDiMu.Fill(iDataSet,3,EventWeight);
      
      MSPlot["KinFit_Prob_AtLeast1Btag_CSVM"]->Fill(kin_prob, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVM);
      MSPlot["KinFit_Chi2_AtLeast1Btag_CSVM"]->Fill(kin_chi2, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVM);
      
      if(topFCNC_KinFit->GetNdof()!=0)
        MSPlot["KinFit_ReducedChi2_AtLeast1Btag_CSVM"]->Fill(kin_chi2ByNdf, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVM);
      
      MSPlot["KinFit_HadWMass_AtLeast1Btag_CSVM"]   ->Fill(kin_hadWmass, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVM);
      MSPlot["KinFit_HadTopMass_AtLeast1Btag_CSVM"] ->Fill(kin_hadtopmass, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVM);
      MSPlot["KinFit_FcncTopMass_AtLeast1Btag_CSVM"]->Fill(kin_fcnctopmass, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVM);
      MSPlot["KinFit_LepZ_Pt_AtLeast1Btag_CSVM"]    ->Fill(kin_lepZpt, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVM);
      MSPlot["KinFit_Lep_DR_AtLeast1Btag_CSVM"]     ->Fill(kin_lepDR, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVM);
      MSPlot["KinFit_Met_AtLeast1Btag_CSVM"]        ->Fill(kin_met, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVM);
      MSPlot["KinFit_Njets_AtLeast1Btag_CSVM"]->Fill(selectedJets.size(), dataSets[iDataSet], true, Luminosity*EventWeight);
      MSPlot["KinFit_LepZ_Pt_Over_JetSyst_Pt_AtLeast1Btag_CSVM"]->Fill(kin_lepZpt/kin_jetSystpt, dataSets[iDataSet], true, Luminosity*EventWeight);
      MSPlot["KinFit_LepZ_DPhi_JetSyst_AtLeast1Btag_CSVM"]->Fill(kin_lepZdphi_jetSyst, dataSets[iDataSet], true, Luminosity*EventWeight);
    }
    // CSV Loose Working point (Veto) _____________________________________________________________________
    else if(bdisc<0.244){ // DO NOT FORGET THE B-TAGGING SF WHEN MC !!!!!!!!!!
      Counters[3][iDataSet] += NormFactor*Luminosity*EventWeight*SFb_CSVL;
      selecTableDiMu.Fill(iDataSet,2,EventWeight);
      
      
      MSPlot["KinFit_Prob_NoBtag_CSVL"]->Fill(kin_prob, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVL);
      MSPlot["KinFit_Chi2_NoBtag_CSVL"]->Fill(kin_chi2, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVL);
      
      if(topFCNC_KinFit->GetNdof()!=0)
        MSPlot["KinFit_ReducedChi2_NoBtag_CSVL"]->Fill(kin_chi2ByNdf, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVL);
      
      MSPlot["KinFit_HadWMass_NoBtag_CSVL"]   ->Fill(kin_hadWmass, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVL);
      MSPlot["KinFit_HadTopMass_NoBtag_CSVL"] ->Fill(kin_hadtopmass, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVL);
      MSPlot["KinFit_FcncTopMass_NoBtag_CSVL"]->Fill(kin_fcnctopmass, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVL);
      MSPlot["KinFit_LepZ_Pt_NoBtag_CSVL"]    ->Fill(kin_lepZpt, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVL);
      MSPlot["KinFit_Lep_DR_NoBtag_CSVL"]     ->Fill(kin_lepDR, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVL);
      MSPlot["KinFit_Met_NoBtag_CSVL"]        ->Fill(kin_met, dataSets[iDataSet], true, Luminosity*EventWeight*SFb_CSVL);
      MSPlot["KinFit_Njets_NoBtag_CSVL"]->Fill(selectedJets.size(), dataSets[iDataSet], true, Luminosity*EventWeight);
      MSPlot["KinFit_LepZ_Pt_Over_JetSyst_Pt_NoBtag_CSVL"]->Fill(kin_lepZpt/kin_jetSystpt, dataSets[iDataSet], true, Luminosity*EventWeight);
      MSPlot["KinFit_LepZ_DPhi_JetSyst_NoBtag_CSVL"]->Fill(kin_lepZdphi_jetSyst, dataSets[iDataSet], true, Luminosity*EventWeight);
    }
*/
  } // loop on events
  //oTree->Print();
  delete topFCNC_KinFit;
  
  cout << "/************** Number of selected events ***************/" <<endl;
  cout << " - Dataset: " << dataSetName << endl;
  cout << " - inclusive: " << Counter << endl;//s[0][iDataSet] << endl;
  cout << "/********************************************************/" <<endl;

  /*
   //Selection tables
   if(UseMuChannel){
   //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergeTTV, bool NP_mass)
   selecTableDiMu.TableCalculator(false, true, true, true, true, true);
   //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
   selecTableDiMu.Write("TopFCNC"+postfix+channelpostfix+comments+"_SelectionTable_DiMu.tex",true,true,true,true,false,false,true);
   }
   */
  //MultiSample plots
  fout->cd();
  inConfigTree->Write();
  oTree->Write();
/*
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    //temp->addText("CMS preliminary");
    string name = it->first;
    name += comments;
    cout<<"Booking MS :"<<name<<endl;
    temp->Draw(false, name, true, true, true, true, true,1,true); // merge TT/QCD/W/Z/ST/
    //Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
    temp->Write(fout, name, false, pathPNG, "pdf");
  }
*/
  //TH1X histograms
  TDirectory* th1dir = fout->mkdir("Histos1D");
  th1dir->cd();
  for(map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
    TH1F *temp = it->second;
    temp->Write();
	  //TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	  //tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  //TH2X histograms
  TDirectory* th2dir = fout->mkdir("Histos2D");
  th2dir->cd();
  for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {
    TH2F *temp = it->second;
    temp->Write();
	  //TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	  //tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  
  //delete
  delete fout;
  delete inFile;
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}
