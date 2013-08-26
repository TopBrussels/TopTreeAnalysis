///////////////////////////
///// TODO & COMMENTS /////
/////////////////////////// 

#include "TStyle.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

#include "TopTreeAnalysis/JESMeasurement/interface/FullKinFit.h"
#include "TopTreeAnalysis/JESMeasurement/interface/JetCombiner.h"
#include "TopTreeAnalysis/JESMeasurement/interface/LightMonster.h"

#include "Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

map<string,TGraphAsymmErrors*> graphAsymmErr;
map<string,TGraphErrors*> graphErr;
map<string,TGraph*> graph;
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;
map<string,TGraph2D*> graph2D;
map<string,TProfile*> profile;

int main (int argc, char *argv[])
{
  bool doPF2PAT = false;

  clock_t start = clock();

  cout << "******************************************" << endl;
  cout << " Beginning of the program for TTbar JES ! " << endl;
  cout << "******************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

  //which systematic to run?
  string systematic = "Nominal";
  if (argc >= 2)
		systematic = string(argv[1]);
  cout << "Systematic to be used:  " << systematic << endl;
  if( ! (systematic == "Nominal" || systematic == "InvertedIso" || systematic == "JESPlus" || systematic == "JESMinus" || systematic == "JERPlus" || systematic == "JERMinus" || systematic == "AlignPlus" || systematic == "AlignMinus" || systematic == "bVSbbarJES" || systematic == "bVSbbarJESPtDep" || systematic == "bTagPlus" || systematic == "bTagMinus" || systematic == "bVSbbarTag" || systematic == "mWdiff" || systematic == "mWdiffSmall" || systematic == "mWdiffSmaller") )
  {
    cout << "Unknown systematic!!!" << endl;
    cout << "Possible options are: Nominal, InvertedIso, JESPlus, JESMinus, JERPlus, JERMinus" << endl;
    exit(-1);
  }
  
  //xml file
  string xmlFileName = "../config/myJESconfig.xml";
  if (argc >= 3)
    xmlFileName=string(argv[2]);
  const char *xmlfile = xmlFileName.c_str();
  cout << "Used config file:  " << xmlfile << endl;
  
  //Output ROOT file
  string rootFileName ("TTbarJES.root");
 	const char *rootfile = rootFileName.c_str();
  
  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  
  /////////////////////////
  // Which decay channel //
  /////////////////////////
  
  bool semiElectron = true; // use semiElectron channel?
  bool semiMuon = true; // use semiMuon channel?
  if(semiElectron && semiMuon) cout << "  --> Using semiMuon and semiElectron channel..." << endl;
  else
  {
    if(semiMuon) cout << " --> Using the semiMuon channel..." << endl;
    else cout << " --> Using the semiElectron channel..." << endl;
  }
  
  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<"Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = anaEnv.Verbose;
  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
 
 	cout << "analysis environment luminosity for rescaling "<< oldLuminosity << endl;
 	
  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) // Rename datasets, taking systematic name into account
  {
    string dataSetName = datasets[i]->Name();
    datasets[i]->SetName( dataSetName+"_"+systematic );
  }
  
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;
	for (unsigned int d = 0; d < datasets.size (); d++)
  {
		//cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
		if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
    string dataSetName = datasets[d]->Name();
		if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
		{
		  Luminosity = datasets[d]->EquivalentLumi();
		  break;
	  }
	}
	if(Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
	
  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootJet* > init_jets_corrected;
  vector < TRootMET* > mets;

  TFile *fout = new TFile (rootfile, "RECREATE");
  //Global variable
  TRootEvent* event = 0;
  
  string pathPNG = "PlotsJES/"; 	
  mkdir(pathPNG.c_str(),0777);
	
	//nof selected events
  float NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
  TH2F* KinFitPNegDeriv = 0;
  
  // Define the plots
  MSPlot["AllMuonsRelPFIso"] = new MultiSamplePlot(datasets, "AllMuonsRelPFIso", 100, 0, 5, "Muons PF Iso");
  MSPlot["AllElectronsRelPFIso"] = new MultiSamplePlot(datasets, "AllElectronsRelPFIso", 100, 0, 5, "Electrons PF Iso");
  MSPlot["SelectedEventsMuonsRelPFIso"] = new MultiSamplePlot(datasets, "SelectedEventsMuonsRelPFIso", 100, 0, 5, "Muons PF Iso");
  MSPlot["SelectedEventsElectronsRelPFIso"] = new MultiSamplePlot(datasets, "SelectedEventsElectronsRelPFIso", 100, 0, 5, "Electrons PF Iso");
  MSPlot["nEventsAfterCutsSemiMu"] = new MultiSamplePlot(datasets, "nEventsAfterCutsSemiMu", 23, -0.5, 22.5, "Nr. of events after each cut, SemiMu");
  MSPlot["nEventsAfterCutsSemiEl"] = new MultiSamplePlot(datasets, "nEventsAfterCutsSemiEl", 23, -0.5, 22.5, "Nr. of events after each cut, SemiEl");
	MSPlot["nPrimaryVertices"] = new MultiSamplePlot(datasets, "nPrimaryVertices", 12, -0.5, 11.5, "Nr. of primary vertices");
	MSPlot["nGoodPrimaryVertices"] = new MultiSamplePlot(datasets, "nGoodPrimaryVertices", 12, -0.5, 11.5, "Nr. of good primary vertices");
  MSPlot["nSelectedMuons"] = new MultiSamplePlot(datasets, "nSelectedMuons", 5, -0.5, 4.5, "Nr. of selected Muons");
  MSPlot["nSelectedElectrons"] = new MultiSamplePlot(datasets, "nSelectedElectrons", 5, -0.5, 4.5, "Nr. of selected Electrons");
  MSPlot["nLooseOtherMuons"] = new MultiSamplePlot(datasets, "nLooseOtherMuons", 5, -0.5, 4.5, "Nr. of Loose Other Muons");
  MSPlot["nLooseElectrons"] = new MultiSamplePlot(datasets, "nLooseElectrons", 5, -0.5, 4.5, "Nr. of Loose Electrons");
  MSPlot["nSelectedJets"] = new MultiSamplePlot(datasets, "nSelectedJets", 10, -0.5, 9.5, "Nr. of Selected Jets");
  
  histo1D["JetPtSelEv"] = new TH1F("JetPtSelEv","JetPtSelEv",100,30,330);
  histo1D["METSelEv"] = new TH1F("METSelEv","METSelEv",50,0,250);
  
  histo1D["FourthJetPt"] = new TH1F("FourthJetPt","FourthJetPt",100,0,100);
  histo1D["FourthJetPtTriggered"] = new TH1F("FourthJetPtTriggered","FourthJetPtTriggered",100,0,100);
  histo1D["AlignSystSF"] = new TH1F("AlignSystSF","AlignSystSF",200,-.001,.001);
  
  histo1D["mWGen"] = new TH1F("mWGen","mWGen",100,60,100);
  histo1D["mTop"] = new TH1F("mTop","mTop",100,100,250);
  histo1D["mTop_genJet"] = new TH1F("mTop_genJet","mTop_genJet",100,100,250);
  histo1D["nJets"] = new TH1F("nJets","nJets",7,3.5,10.5);
  histo2D["nJets_VS_mTop"] = new TH2F("nJets_VS_mTop","nJets_VS_mTop",7,3.5,10.5,100,100,250);
  histo2D["nJets_VS_mTopGenJet"] = new TH2F("nJets_VS_mTopGenJet","nJets_VS_mTopGenJet",7,3.5,10.5,100,100,250);
  histo1D["dRMin"] = new TH1F("dRMin","dRMin",25,0.5,3.);
  histo2D["dRMin_VS_mTop"] = new TH2F("dRMin_VS_mTop","dRMin_VS_mTop",25,0.5,3.,100,100,250);
  histo2D["dRMin_VS_mTopGenJet"] = new TH2F("dRMin_VS_mTopGenJet","dRMin_VS_mTopGenJet",25,0.5,3.,100,100,250);
  histo1D["dRLights"] = new TH1F("dRLights","dRLights",25,0.5,3.5);
  histo2D["dRLights_VS_mTop"] = new TH2F("dRLights_VS_mTop","dRLights_VS_mTop",25,0.5,3.5,100,100,250);
  histo2D["dRLights_VS_mTopGenJet"] = new TH2F("dRLights_VS_mTopGenJet","dRLights_VS_mTopGenJet",25,0.5,3.5,100,100,250);
  histo1D["MinDRLightB"] = new TH1F("MinDRLightB","MinDRLightB",25,0.5,3.5);
  histo2D["MinDRLightB_VS_mTop"] = new TH2F("MinDRLightB_VS_mTop","MinDRLightB_VS_mTop",25,0.5,3.5,100,100,250);
  histo2D["MinDRLightB_VS_mTopGenJet"] = new TH2F("MinDRLightB_VS_mTopGenJet","MinDRLightB_VS_mTopGenJet",25,0.5,3.5,100,100,250);
  histo1D["MET"] = new TH1F("MET","MET",25,0,250);
  histo2D["MET_VS_mTop"] = new TH2F("MET_VS_mTop","MET_VS_mTop",25,0,250,100,100,250);
  histo2D["MET_VS_mTopGenJet"] = new TH2F("MET_VS_mTopGenJet","MET_VS_mTopGenJet",25,0,250,100,100,250);
  histo1D["HT"] = new TH1F("HT","HT",25,100,800);
  histo2D["HT_VS_mTop"] = new TH2F("HT_VS_mTop","HT_VS_mTop",25,100,800,100,100,250);
  histo2D["HT_VS_mTopGenJet"] = new TH2F("HT_VS_mTopGenJet","HT_VS_mTopGenJet",25,100,800,100,100,250);
  histo1D["mTTbar"] = new TH1F("mTTbar","mTTbar",25,200,1000);
  histo2D["mTTbar_VS_mTop"] = new TH2F("mTTbar_VS_mTop","mTTbar_VS_mTop",25,200,1000,100,100,250);
  histo2D["mTTbar_VS_mTopGenJet"] = new TH2F("mTTbar_VS_mTopGenJet","mTTbar_VS_mTopGenJet",25,200,1000,100,100,250);
  histo1D["PtTTbar"] = new TH1F("PtTTbar","PtTTbar",25,0,200);
  histo2D["PtTTbar_VS_mTop"] = new TH2F("PtTTbar_VS_mTop","PtTTbar_VS_mTop",25,0,200,100,100,250);
  histo2D["PtTTbar_VS_mTopGenJet"] = new TH2F("PtTTbar_VS_mTopGenJet","PtTTbar_VS_mTopGenJet",25,0,200,100,100,250);
  histo1D["nBtags"] = new TH1F("nBtags","nBtags",5,-0.5,4.5);
  histo2D["nBtags_VS_mTop"] = new TH2F("nBtags_VS_mTop","nBtags_VS_mTop",5,-0.5,4.5,100,100,250);
  histo2D["nBtags_VS_mTopGenJet"] = new TH2F("nBtags_VS_mTopGenJet","nBtags_VS_mTopGenJet",5,-0.5,4.5,100,100,250);
  histo1D["PtTop"] = new TH1F("PtTop","PtTop",25,0,400);
  histo2D["PtTop_VS_mTop"] = new TH2F("PtTop_VS_mTop","PtTop_VS_mTop",25,0,400,100,100,250);
  histo2D["PtTop_VS_mTopGenJet"] = new TH2F("PtTop_VS_mTopGenJet","PtTop_VS_mTopGenJet",25,0,400,100,100,250);
  histo1D["EtaTop"] = new TH1F("EtaTop","EtaTop",25,-5,5);
  histo2D["EtaTop_VS_mTop"] = new TH2F("EtaTop_VS_mTop","EtaTop_VS_mTop",25,-5,5,100,100,250);
  histo2D["EtaTop_VS_mTopGenJet"] = new TH2F("EtaTop_VS_mTopGenJet","EtaTop_VS_mTopGenJet",25,-5,5,100,100,250);
  histo1D["PtBjet"] = new TH1F("PtBjet","PtBjet",25,30,280);
  histo2D["PtBjet_VS_mTop"] = new TH2F("PtBjet_VS_mTop","PtBjet_VS_mTop",25,30,280,100,100,250);
  histo2D["PtBjet_VS_mTopGenJet"] = new TH2F("PtBjet_VS_mTopGenJet","PtBjet_VS_mTopGenJet",25,30,280,100,100,250);
  histo1D["EtaBjet"] = new TH1F("EtaBjet","EtaBjet",25,-2.5,2.5);
  histo2D["EtaBjet_VS_mTop"] = new TH2F("EtaBjet_VS_mTop","EtaBjet_VS_mTop",25,-2.5,2.5,100,100,250);
  histo2D["EtaBjet_VS_mTopGenJet"] = new TH2F("EtaBjet_VS_mTopGenJet","EtaBjet_VS_mTopGenJet",25,-2.5,2.5,100,100,250);
  histo1D["jetAreaBjet"] = new TH1F("jetAreaBjet","jetAreaBjet",25,0.5,1.);
  histo2D["jetAreaBjet_VS_mTop"] = new TH2F("jetAreaBjet_VS_mTop","jetAreaBjet_VS_mTop",25,0.5,1.,100,100,250);
  histo2D["jetAreaBjet_VS_mTopGenJet"] = new TH2F("jetAreaBjet_VS_mTopGenJet","jetAreaBjet_VS_mTopGenJet",25,0.5,1.,100,100,250);
  histo1D["jetAreaTotal"] = new TH1F("jetAreaTotal","jetAreaTotal",25,1.5,3.);
  histo2D["jetAreaTotal_VS_mTop"] = new TH2F("jetAreaTotal_VS_mTop","jetAreaTotal_VS_mTop",25,1.5,3.,100,100,250);
  histo2D["jetAreaTotal_VS_mTopGenJet"] = new TH2F("jetAreaTotal_VS_mTopGenJet","jetAreaTotal_VS_mTopGenJet",25,1.5,3.,100,100,250);
  histo1D["nParticlesBJet"] = new TH1F("nParticlesBJet","nParticlesBJet",30,-0.5,59.5);
  histo2D["nParticlesBJet_VS_mTop"] = new TH2F("nParticlesBJet_VS_mTop","nParticlesBJet_VS_mTop",30,-0.5,59.5,100,100,250);
  histo2D["nParticlesBJet_VS_mTopGenJet"] = new TH2F("nParticlesBJet_VS_mTopGenJet","nParticlesBJet_VS_mTopGenJet",30,-0.5,59.5,100,100,250);
  histo1D["nParticlesTotal"] = new TH1F("nParticlesTotal","nParticlesTotal",25,-29.5,129.5);
  histo2D["nParticlesTotal_VS_mTop"] = new TH2F("nParticlesTotal_VS_mTop","nParticlesTotal_VS_mTop",25,-29.5,129.5,100,100,250);
  histo2D["nParticlesTotal_VS_mTopGenJet"] = new TH2F("nParticlesTotal_VS_mTopGenJet","nParticlesTotal_VS_mTopGenJet",25,-29.5,129.5,100,100,250);
  histo1D["nChargedBJet"] = new TH1F("nChargedBJet","nChargedBJet",18,-0.5,35.5);
  histo2D["nChargedBJet_VS_mTop"] = new TH2F("nChargedBJet_VS_mTop","nChargedBJet_VS_mTop",18,-0.5,35.5,100,100,250);
  histo2D["nChargedBJet_VS_mTopGenJet"] = new TH2F("nChargedBJet_VS_mTopGenJet","nChargedBJet_VS_mTopGenJet",18,-0.5,35.5,100,100,250);
  histo1D["nChargedTotal"] = new TH1F("nChargedTotal","nChargedTotal",30,5.5,65.5);
  histo2D["nChargedTotal_VS_mTop"] = new TH2F("nChargedTotal_VS_mTop","nChargedTotal_VS_mTop",30,5.5,65.5,100,100,250);
  histo2D["nChargedTotal_VS_mTopGenJet"] = new TH2F("nChargedTotal_VS_mTopGenJet","nChargedTotal_VS_mTopGenJet",30,5.5,65.5,100,100,250);
  
  histo1D["nParticlesInJet"] = new TH1F("nParticlesInJet","nParticlesInJet",60,-0.5,59.5);
  histo1D["nChargedParticles"] = new TH1F("nChargedParticles","nChargedParticles",35,-0.5,34.5);
  histo1D["jetArea"] = new TH1F("jetArea","jetArea",50,0.5,1.);
  histo1D["mW_dRsmall"] = new TH1F("mW_dRsmall","mW_dRsmall",50,0,200);
  histo1D["mW_dRlarge"] = new TH1F("mW_dRlarge","mW_dRlarge",50,0,200);
  histo1D["mTop_dRsmall"] = new TH1F("mTop_dRsmall","mTop_dRsmall",50,0,300);
  histo1D["mTop_dRlarge"] = new TH1F("mTop_dRlarge","mTop_dRlarge",50,0,300);
  
  int nBinsJet = 12;
  float binningDR[13] = {0.5, 0.65, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.5};
  float binningJetArea[13] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.9, 0.95};
  float binningNpart[13] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65};
  for(size_t i=0; i<nBinsJet; i++)
  {
    stringstream ss; ss << i;
    histo1D["lightJetResp_DR_"+ss.str()] = new TH1F(("lightJetResp_DR_"+ss.str()).c_str(),("lightJetResp_DR_"+ss.str()).c_str(),100,0.,2.);
    histo1D["bJetResp_DR_"+ss.str()] = new TH1F(("bJetResp_DR_"+ss.str()).c_str(),("bJetResp_DR_"+ss.str()).c_str(),100,0.,2.);
    histo1D["lightJetResp_Area_"+ss.str()] = new TH1F(("lightJetResp_Area_"+ss.str()).c_str(),("lightJetResp_Area_"+ss.str()).c_str(),100,0.,2.);
    histo1D["bJetResp_Area_"+ss.str()] = new TH1F(("bJetResp_Area_"+ss.str()).c_str(),("bJetResp_Area_"+ss.str()).c_str(),100,0.,2.);
    histo1D["lightJetResp_nPart_"+ss.str()] = new TH1F(("lightJetResp_nPart_"+ss.str()).c_str(),("lightJetResp_nPart_"+ss.str()).c_str(),100,0.,2.);
    histo1D["bJetResp_nPart_"+ss.str()] = new TH1F(("bJetResp_nPart_"+ss.str()).c_str(),("bJetResp_nPart_"+ss.str()).c_str(),100,0.,2.);
  }
  
  // JES plots
  profile["JetResp_Raw_0to10PV"] = new TProfile("JetResp_Raw_0to10PV","JetResp_Raw_0to10PV",200,5,505);
  profile["JetResp_Raw_10to15PV"] = new TProfile("JetResp_Raw_10to15PV","JetResp_Raw_10to15PV",200,5,505);
  profile["JetResp_Raw_15to20PV"] = new TProfile("JetResp_Raw_15to20PV","JetResp_Raw_15to20PV",200,5,505);
  profile["JetResp_Raw_20to25PV"] = new TProfile("JetResp_Raw_20to25PV","JetResp_Raw_20to25PV",200,5,505);
  profile["JetResp_Raw_25to30PV"] = new TProfile("JetResp_Raw_25to30PV","JetResp_Raw_25to30PV",200,5,505);
  profile["JetResp_Raw_30PlusPV"] = new TProfile("JetResp_Raw_30PlusPV","JetResp_Raw_30PlusPV",200,5,505);
  
  profile["JetResp_L1_0to10PV"] = new TProfile("JetResp_L1_0to10PV","JetResp_L1_0to10PV",200,5,505);
  profile["JetResp_L1_10to15PV"] = new TProfile("JetResp_L1_10to15PV","JetResp_L1_10to15PV",200,5,505);
  profile["JetResp_L1_15to20PV"] = new TProfile("JetResp_L1_15to20PV","JetResp_L1_15to20PV",200,5,505);
  profile["JetResp_L1_20to25PV"] = new TProfile("JetResp_L1_20to25PV","JetResp_L1_20to25PV",200,5,505);
  profile["JetResp_L1_25to30PV"] = new TProfile("JetResp_L1_25to30PV","JetResp_L1_25to30PV",200,5,505);
  profile["JetResp_L1_30PlusPV"] = new TProfile("JetResp_L1_30PlusPV","JetResp_L1_30PlusPV",200,5,505);
  
  profile["JetResp_L1L2L3_0to10PV"] = new TProfile("JetResp_L1L2L3_0to10PV","JetResp_L1L2L3_0to10PV",200,5,505);
  profile["JetResp_L1L2L3_10to15PV"] = new TProfile("JetResp_L1L2L3_10to15PV","JetResp_L1L2L3_10to15PV",200,5,505);
  profile["JetResp_L1L2L3_15to20PV"] = new TProfile("JetResp_L1L2L3_15to20PV","JetResp_L1L2L3_15to20PV",200,5,505);
  profile["JetResp_L1L2L3_20to25PV"] = new TProfile("JetResp_L1L2L3_20to25PV","JetResp_L1L2L3_20to25PV",200,5,505);
  profile["JetResp_L1L2L3_25to30PV"] = new TProfile("JetResp_L1L2L3_25to30PV","JetResp_L1L2L3_25to30PV",200,5,505);
  profile["JetResp_L1L2L3_30PlusPV"] = new TProfile("JetResp_L1L2L3_30PlusPV","JetResp_L1L2L3_30PlusPV",200,5,505);
  
  profile["JetResp_L1L2L3_Incl"] = new TProfile("JetResp_L1L2L3_Incl","JetResp_L1L2L3_Incl",200,5,505);
  profile["JetResp_L1L2L3_smallEta_Incl"] = new TProfile("JetResp_L1L2L3_smallEta_Incl","JetResp_smallEta_L1L2L3_Incl",200,5,505);
  profile["JetResp_L1L2L3_mediumEta_Incl"] = new TProfile("JetResp_L1L2L3_mediumEta_Incl","JetResp_L1L2L3_mediumEta_Incl",200,5,505);
  profile["JetResp_L1L2L3_largeEta_Incl"] = new TProfile("JetResp_L1L2L3_largeEta_Incl","JetResp_L1L2L3_largeEta_Incl",200,5,505);
  
  profile["JetResp_L1L2L3_gluon"] = new TProfile("JetResp_L1L2L3_gluon","JetResp_L1L2L3_gluon",200,5,505);
  profile["JetResp_L1L2L3_uds"] = new TProfile("JetResp_L1L2L3_uds","JetResp_L1L2L3_uds",200,5,505);
  profile["JetResp_L1L2L3_c"] = new TProfile("JetResp_L1L2L3_c","JetResp_L1L2L3_c",200,5,505);
  profile["JetResp_L1L2L3_b"] = new TProfile("JetResp_L1L2L3_b","JetResp_L1L2L3_b",200,5,505);
  profile["JetResp_L1L2L3_unm"] = new TProfile("JetResp_L1L2L3_unm","JetResp_L1L2L3_unm",200,5,505);
  
  profile["JetResp_L1L2L3_smallEta_gluon"] = new TProfile("JetResp_L1L2L3_smallEta_gluon","JetResp_smallEta_L1L2L3_gluon",200,5,505);
  profile["JetResp_L1L2L3_smallEta_uds"] = new TProfile("JetResp_L1L2L3_smallEta_uds","JetResp_smallEta_L1L2L3_uds",200,5,505);
  profile["JetResp_L1L2L3_smallEta_c"] = new TProfile("JetResp_L1L2L3_smallEta_c","JetResp_L1L2L3_smallEta_c",200,5,505);
  profile["JetResp_L1L2L3_smallEta_b"] = new TProfile("JetResp_L1L2L3_smallEta_b","JetResp_L1L2L3_smallEta_b",200,5,505);
  profile["JetResp_L1L2L3_smallEta_unm"] = new TProfile("JetResp_L1L2L3_smallEta_unm","JetResp_L1L2L3_smallEta_unm",200,5,505);
  
  // JES unc plots
  graph["JecUnc_Pt30_PileUp"] = new TGraph(5000);
  graph["JecUnc_Pt30_PileUp"]->SetNameTitle("JecUnc_Pt30_PileUp","JecUnc_Pt30_PileUp");
  graph["JecUnc_Pt30_Relative"] = new TGraph(5000);
  graph["JecUnc_Pt30_Relative"]->SetNameTitle("JecUnc_Pt30_Relative","JecUnc_Pt30_Relative");
  graph["JecUnc_Pt30_Time"] = new TGraph(5000);
  graph["JecUnc_Pt30_Time"]->SetNameTitle("JecUnc_Pt30_Time","JecUnc_Pt30_Time");
  graph["JecUnc_Pt30_Flavour"] = new TGraph(5000);
  graph["JecUnc_Pt30_Flavour"]->SetNameTitle("JecUnc_Pt30_Flavour","JecUnc_Pt30_Flavour");
  graph["JecUnc_Pt30_Absolute"] = new TGraph(5000);
  graph["JecUnc_Pt30_Absolute"]->SetNameTitle("JecUnc_Pt30_Absolute","JecUnc_Pt30_Absolute");
  graph["JecUnc_Pt30_Total"] = new TGraph(5000);
  graph["JecUnc_Pt30_Total"]->SetNameTitle("JecUnc_Pt30_Total","JecUnc_Pt30_Total");
  
  graph["JecUnc_Pt50_PileUp"] = new TGraph(5000);
  graph["JecUnc_Pt50_PileUp"]->SetNameTitle("JecUnc_Pt50_PileUp","JecUnc_Pt50_PileUp");
  graph["JecUnc_Pt50_Relative"] = new TGraph(5000);
  graph["JecUnc_Pt50_Relative"]->SetNameTitle("JecUnc_Pt50_Relative","JecUnc_Pt50_Relative");
  graph["JecUnc_Pt50_Time"] = new TGraph(5000);
  graph["JecUnc_Pt50_Time"]->SetNameTitle("JecUnc_Pt50_Time","JecUnc_Pt50_Time");
  graph["JecUnc_Pt50_Flavour"] = new TGraph(5000);
  graph["JecUnc_Pt50_Flavour"]->SetNameTitle("JecUnc_Pt50_Flavour","JecUnc_Pt50_Flavour");
  graph["JecUnc_Pt50_Absolute"] = new TGraph(5000);
  graph["JecUnc_Pt50_Absolute"]->SetNameTitle("JecUnc_Pt50_Absolute","JecUnc_Pt50_Absolute");
  graph["JecUnc_Pt50_Total"] = new TGraph(5000);
  graph["JecUnc_Pt50_Total"]->SetNameTitle("JecUnc_Pt50_Total","JecUnc_Pt50_Total");
  
  graph["JecUnc_Pt70_PileUp"] = new TGraph(5000);
  graph["JecUnc_Pt70_PileUp"]->SetNameTitle("JecUnc_Pt70_PileUp","JecUnc_Pt70_PileUp");
  graph["JecUnc_Pt70_Relative"] = new TGraph(5000);
  graph["JecUnc_Pt70_Relative"]->SetNameTitle("JecUnc_Pt70_Relative","JecUnc_Pt70_Relative");
  graph["JecUnc_Pt70_Time"] = new TGraph(5000);
  graph["JecUnc_Pt70_Time"]->SetNameTitle("JecUnc_Pt70_Time","JecUnc_Pt70_Time");
  graph["JecUnc_Pt70_Flavour"] = new TGraph(5000);
  graph["JecUnc_Pt70_Flavour"]->SetNameTitle("JecUnc_Pt70_Flavour","JecUnc_Pt70_Flavour");
  graph["JecUnc_Pt70_Absolute"] = new TGraph(5000);
  graph["JecUnc_Pt70_Absolute"]->SetNameTitle("JecUnc_Pt70_Absolute","JecUnc_Pt70_Absolute");
  graph["JecUnc_Pt70_Total"] = new TGraph(5000);
  graph["JecUnc_Pt70_Total"]->SetNameTitle("JecUnc_Pt70_Total","JecUnc_Pt70_Total");
  
  graph["JecUnc_Eta0_PileUp"] = new TGraph(1921);
  graph["JecUnc_Eta0_PileUp"]->SetNameTitle("JecUnc_Eta0_PileUp","JecUnc_Eta0_PileUp");
  graph["JecUnc_Eta0_Relative"] = new TGraph(1921);
  graph["JecUnc_Eta0_Relative"]->SetNameTitle("JecUnc_Eta0_Relative","JecUnc_Eta0_Relative");
  graph["JecUnc_Eta0_Time"] = new TGraph(1921);
  graph["JecUnc_Eta0_Time"]->SetNameTitle("JecUnc_Eta0_Time","JecUnc_Eta0_Time");
  graph["JecUnc_Eta0_Flavour"] = new TGraph(1921);
  graph["JecUnc_Eta0_Flavour"]->SetNameTitle("JecUnc_Eta0_Flavour","JecUnc_Eta0_Flavour");
  graph["JecUnc_Eta0_Absolute"] = new TGraph(1921);
  graph["JecUnc_Eta0_Absolute"]->SetNameTitle("JecUnc_Eta0_Absolute","JecUnc_Eta0_Absolute");
  graph["JecUnc_Eta0_Total"] = new TGraph(1921);
  graph["JecUnc_Eta0_Total"]->SetNameTitle("JecUnc_Eta0_Total","JecUnc_Eta0_Total");
  
  graph["JecUnc_Eta1p2_PileUp"] = new TGraph(1921);
  graph["JecUnc_Eta1p2_PileUp"]->SetNameTitle("JecUnc_Eta1p2_PileUp","JecUnc_Eta1p2_PileUp");
  graph["JecUnc_Eta1p2_Relative"] = new TGraph(1921);
  graph["JecUnc_Eta1p2_Relative"]->SetNameTitle("JecUnc_Eta1p2_Relative","JecUnc_Eta1p2_Relative");
  graph["JecUnc_Eta1p2_Time"] = new TGraph(1921);
  graph["JecUnc_Eta1p2_Time"]->SetNameTitle("JecUnc_Eta1p2_Time","JecUnc_Eta1p2_Time");
  graph["JecUnc_Eta1p2_Flavour"] = new TGraph(1921);
  graph["JecUnc_Eta1p2_Flavour"]->SetNameTitle("JecUnc_Eta1p2_Flavour","JecUnc_Eta1p2_Flavour");
  graph["JecUnc_Eta1p2_Absolute"] = new TGraph(1921);
  graph["JecUnc_Eta1p2_Absolute"]->SetNameTitle("JecUnc_Eta1p2_Absolute","JecUnc_Eta1p2_Absolute");
  graph["JecUnc_Eta1p2_Total"] = new TGraph(1921);
  graph["JecUnc_Eta1p2_Total"]->SetNameTitle("JecUnc_Eta1p2_Total","JecUnc_Eta1p2_Total");
  
  graph["JecUnc_Eta2p3_PileUp"] = new TGraph(1921);
  graph["JecUnc_Eta2p3_PileUp"]->SetNameTitle("JecUnc_Eta2p3_PileUp","JecUnc_Eta2p3_PileUp");
  graph["JecUnc_Eta2p3_Relative"] = new TGraph(1921);
  graph["JecUnc_Eta2p3_Relative"]->SetNameTitle("JecUnc_Eta2p3_Relative","JecUnc_Eta2p3_Relative");
  graph["JecUnc_Eta2p3_Time"] = new TGraph(1921);
  graph["JecUnc_Eta2p3_Time"]->SetNameTitle("JecUnc_Eta2p3_Time","JecUnc_Eta2p3_Time");
  graph["JecUnc_Eta2p3_Flavour"] = new TGraph(1921);
  graph["JecUnc_Eta2p3_Flavour"]->SetNameTitle("JecUnc_Eta2p3_Flavour","JecUnc_Eta2p3_Flavour");
  graph["JecUnc_Eta2p3_Absolute"] = new TGraph(1921);
  graph["JecUnc_Eta2p3_Absolute"]->SetNameTitle("JecUnc_Eta2p3_Absolute","JecUnc_Eta2p3_Absolute");
  graph["JecUnc_Eta2p3_Total"] = new TGraph(1921);
  graph["JecUnc_Eta2p3_Total"]->SetNameTitle("JecUnc_Eta2p3_Total","JecUnc_Eta2p3_Total");
  
  graph["FlavUnc_Eta0_gluon"] = new TGraph(1921);
  graph["FlavUnc_Eta0_gluon"]->SetNameTitle("FlavUnc_Eta0_gluon","FlavUnc_Eta0_gluon");
  graph["FlavUnc_Eta0_uds"] = new TGraph(1921);
  graph["FlavUnc_Eta0_uds"]->SetNameTitle("FlavUnc_Eta0_uds","FlavUnc_Eta0_uds");
  graph["FlavUnc_Eta0_c"] = new TGraph(1921);
  graph["FlavUnc_Eta0_c"]->SetNameTitle("FlavUnc_Eta0_c","FlavUnc_Eta0_c");
  graph["FlavUnc_Eta0_b"] = new TGraph(1921);
  graph["FlavUnc_Eta0_b"]->SetNameTitle("FlavUnc_Eta0_b","FlavUnc_Eta0_b");
  
  graph["FlavUnc_Eta1p9_gluon"] = new TGraph(1921);
  graph["FlavUnc_Eta1p9_gluon"]->SetNameTitle("FlavUnc_Eta1p9_gluon","FlavUnc_Eta1p9_gluon");
  graph["FlavUnc_Eta1p9_uds"] = new TGraph(1921);
  graph["FlavUnc_Eta1p9_uds"]->SetNameTitle("FlavUnc_Eta1p9_uds","FlavUnc_Eta1p9_uds");
  graph["FlavUnc_Eta1p9_c"] = new TGraph(1921);
  graph["FlavUnc_Eta1p9_c"]->SetNameTitle("FlavUnc_Eta1p9_c","FlavUnc_Eta1p9_c");
  graph["FlavUnc_Eta1p9_b"] = new TGraph(1921);
  graph["FlavUnc_Eta1p9_b"]->SetNameTitle("FlavUnc_Eta1p9_b","FlavUnc_Eta1p9_b");
  
  graph["FlavUnc_Eta2p3_gluon"] = new TGraph(1921);
  graph["FlavUnc_Eta2p3_gluon"]->SetNameTitle("FlavUnc_Eta2p3_gluon","FlavUnc_Eta2p3_gluon");
  graph["FlavUnc_Eta2p3_uds"] = new TGraph(1921);
  graph["FlavUnc_Eta2p3_uds"]->SetNameTitle("FlavUnc_Eta2p3_uds","FlavUnc_Eta2p3_uds");
  graph["FlavUnc_Eta2p3_c"] = new TGraph(1921);
  graph["FlavUnc_Eta2p3_c"]->SetNameTitle("FlavUnc_Eta2p3_c","FlavUnc_Eta2p3_c");
  graph["FlavUnc_Eta2p3_b"] = new TGraph(1921);
  graph["FlavUnc_Eta2p3_b"]->SetNameTitle("FlavUnc_Eta2p3_b","FlavUnc_Eta2p3_b");
  
  // some b-tag plots
  histo1D["CSV_b"] = new TH1F("CSV_b","CSV_b",100,0.,1.);
  histo1D["CSV_c"] = new TH1F("CSV_c","CSV_c",100,0.,1.);
  histo1D["CSV_uds"] = new TH1F("CSV_uds","CSV_uds",100,0.,1.);
  histo1D["CSV_g"] = new TH1F("CSV_g","CSV_g",100,0.,1.);
  histo1D["CSV_unm"] = new TH1F("CSV_unm","CSV_unm",100,0.,1.);
  
  // some EleId plots
  histo1D["ElIdMVA_matched"] = new TH1F("ElIdMVA_matched","ElIdMVA_matched",100,-1,1);
  histo1D["ElIdMVA_unmatched"] = new TH1F("ElIdMVA_unmatched","ElIdMVA_unmatched",100,-1,1);
  histo1D["ElIso_matched"] = new TH1F("ElIso_matched","ElIso_matched",100,0,0.4);
  histo1D["ElIso_unmatched"] = new TH1F("ElIso_unmatched","ElIso_unmatched",100,0,0.4);
  
  // some MuId plots
  histo1D["MuIso_matched"] = new TH1F("MuIso_matched","MuIso_matched",100,0,0.2);
  histo1D["MuIso_unmatched"] = new TH1F("MuIso_unmatched","MuIso_unmatched",100,0,0.2);
  histo1D["MuNPixelHits_matched"] = new TH1F("MuNPixelHits_matched","MuNPixelHits_matched",10,-0.5,9.5);
  histo1D["MuNPixelHits_unmatched"] = new TH1F("MuNPixelHits_unmatched","MuNPixelHits_unmatched",10,-0.5,9.5);
  histo1D["Mud0_matched"] = new TH1F("Mud0_matched","Mud0_matched",100,-0.05,0.05);
  histo1D["Mud0_unmatched"] = new TH1F("Mud0_unmatched","Mud0_unmatched",100,-0.05,0.05);
  histo1D["MuChi2_matched"] = new TH1F("MuChi2_matched","MuChi2_matched",100,0,50);
  histo1D["MuChi2_unmatched"] = new TH1F("MuChi2_unmatched","MuChi2_unmatched",100,0,50);
  histo1D["MuDZ_matched"] = new TH1F("MuDZ_matched","MuDZ_matched",100,-2,2);
  histo1D["MuDZ_unmatched"] = new TH1F("MuDZ_unmatched","MuDZ_unmatched",100,-2,2);
  histo1D["MuNsegm_matched"] = new TH1F("MuNsegm_matched","MuNsegm_matched",20,-0.5,19.5);
  histo1D["MuNsegm_unmatched"] = new TH1F("MuNsegm_unmatched","MuNsegm_unmatched",20,-0.5,19.5);
  histo1D["MuNTrackLayers_matched"] = new  TH1F("MuNTrackLayers_matched","MuNTrackLayers_matched",30,-0.5,29.5);
  histo1D["MuNTrackLayers_unmatched"] = new  TH1F("MuNTrackLayers_unmatched","MuNTrackLayers_unmatched",30,-0.5,29.5);
  
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
  
  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("preselected"));
  CutsSelecTableSemiMu.push_back(string("trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
  CutsSelecTableSemiMu.push_back(string("Veto electron"));
  
  vector<string> CutsSelecTableSemiEl;
  CutsSelecTableSemiEl.push_back(string("preselected"));
  CutsSelecTableSemiEl.push_back(string("trigged"));
  CutsSelecTableSemiEl.push_back(string("Good PV"));
  CutsSelecTableSemiEl.push_back(string("1 selected electron"));
  CutsSelecTableSemiEl.push_back(string("Veto muon"));
  CutsSelecTableSemiEl.push_back(string("Veto 2nd electron"));
  
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableSemiEl.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 2 jets"));
  CutsSelecTableSemiEl.push_back(string("$\\geq$ 2 jets"));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 3 jets"));
  CutsSelecTableSemiEl.push_back(string("$\\geq$ 3 jets"));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 4 jets"));
  CutsSelecTableSemiEl.push_back(string("$\\geq$ 4 jets"));
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  selecTableSemiMu.SetLuminosity(Luminosity);
  SelectionTable selecTableSemiEl(CutsSelecTableSemiEl, datasets);
  selecTableSemiEl.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  
  ////////////////////////////////////////
  /// initialize LumiReWeighting stuff ///
  ////////////////////////////////////////
  
  cout << "Initializing LumiReWeighting stuff" << endl;
  LumiReWeighting LumiWeights = LumiReWeighting("PileUpReweighting/pileup_MC_Summer12_S10.root", "PileUpReweighting/pileup_2012Data53X_UpToRun208686_Mu/nominal.root", "pileup", "pileup");
  
  /////////////////////////////
  /// ResolutionFit Stuff
  /////////////////////////////
  
  bool CalculateResolutions = false; // If false, the resolutions will be loaded from a previous calculation
  bool ResolutionsClosure = false;
  
  ResolutionFit *resFitLightJets = 0, *resFitLightJets_PtTopRew = 0, *resFitBJets = 0, *resFitBJets_PtTopRew = 0, *resFitLightJetsL7 = 0, *resFitBJetsL7 = 0, *resFitBJets_B = 0, *resFitBJets_Bbar = 0;
  if(CalculateResolutions)
  {
    resFitLightJets = new ResolutionFit("LightJet");
    resFitLightJets_PtTopRew = new ResolutionFit("LightJet_PtTopRew");
    resFitBJets = new ResolutionFit("BJet");
    resFitBJets_PtTopRew = new ResolutionFit("BJet_PtTopRew");
    resFitBJets_B = new ResolutionFit("BJet_B");
    resFitBJets_Bbar = new ResolutionFit("BJet_Bbar");
    if(ResolutionsClosure)
    {
      resFitLightJets->LoadResolutions("resolutions/lightJetReso.root");
      resFitLightJets_PtTopRew->LoadResolutions("resolutions/lightJetReso_PtTopRew.root");
      resFitBJets->LoadResolutions("resolutions/bJetReso.root");
      resFitBJets_PtTopRew->LoadResolutions("resolutions/bJetReso_PtTopRew.root");
      resFitBJets_B->LoadResolutions("resolutions/bJetReso_PtTopRew.root");
      resFitBJets_Bbar->LoadResolutions("resolutions/bJetReso_PtTopRew.root");
    }
  }
  else
  {
    resFitLightJets = new ResolutionFit("LightJet");
    resFitLightJets->LoadResolutions("resolutions/lightJetReso.root");
    resFitLightJetsL7 = new ResolutionFit("LightJetL7");
    resFitLightJetsL7->LoadResolutions("resolutions/lightJetReso_AfterL7.root");
    resFitBJets = new ResolutionFit("BJet");
    resFitBJets->LoadResolutions("resolutions/bJetReso.root");
    resFitBJetsL7 = new ResolutionFit("BJetL7");
    resFitBJetsL7->LoadResolutions("resolutions/bJetReso_AfterL7.root");
  }
  if (verbose > 0)
    cout << " - ResolutionFit instantiated ..." << endl;
  
  /////////////////////////////////////
  // Initialize JetCombination Stuff //
  /////////////////////////////////////
	bool Tprime = false; // If false, regular variables are used in MVA

  bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff
  string MVAmethod = "Likelihood"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)
    
  JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, Tprime);
  JetCombiner* jetCombinerGenPartGenJets = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, Tprime);
  JetCombiner* jetCombinerGenJetsPFJets = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, Tprime);
  if (verbose > 0)
    cout << " - JetCombiner instantiated ..." << endl;
  
  /////////////////////////////////////////////////
  // Which analysis to execute? Top Mass or JES? //
  /////////////////////////////////////////////////
  
  bool measureTopMass = false;
  bool measureTopMassDifference = true;
  
  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
		
		//open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;
    
    string dataSetName = datasets[d]->Name();
    string previousFilename = "";
    int iFile = -1;
    
    /////////////////////////////////////
   	/// Initialize JEC factors
   	/////////////////////////////////////
   	
    vector<JetCorrectorParameters> vCorrParam;
    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) // Data!
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("JECFiles/FT_53_V21_AN4_Summer13_Data_L1FastJet_AK5PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("JECFiles/FT_53_V21_AN4_Summer13_Data_L2Relative_AK5PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("JECFiles/FT_53_V21_AN4_Summer13_Data_L3Absolute_AK5PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("JECFiles/FT_53_V21_AN4_Summer13_Data_L2L3Residual_AK5PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
    }
    else
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("JECFiles/START53_V23_Summer13_L1FastJet_AK5PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("JECFiles/START53_V23_Summer13_L2Relative_AK5PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("JECFiles/START53_V23_Summer13_L3Absolute_AK5PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
    }
//    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START53_V23_Summer13_Uncertainty_AK5PFchs.txt");
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "SubTotalMC")));
    
    JetTools *jetTools = new JetTools(vCorrParam, "JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "JECFiles/flavourFractions_4OrMoreJets.root", true, true); // JetTools(vCorrParam, jecUnc, true);
    
    JetTools *jetTools_Exact4 = new JetTools(vCorrParam, "JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "JECFiles/flavourFractions_Exact4Jets.root", true, true); // JetTools(vCorrParam, jecUnc, true);
    
    histo2D["jesUnc"] = new TH2F("jesUnc", "jesUnc", 101, 29.5, 130.5, 25, -0.05, 2.45);
    graph2D["jesUnc_graph2D"] = new TGraph2D(101*25);
    histo1D["jesUncPt30"] = new TH1F("jesUncPt30", "jesUncPt30", 25, -0.05, 2.45);
    histo1D["jesUncPt40"] = new TH1F("jesUncPt40", "jesUncPt40", 25, -0.05, 2.45);
    histo1D["jesUncPt50"] = new TH1F("jesUncPt50", "jesUncPt50", 25, -0.05, 2.45);
    histo1D["jesUncPt60"] = new TH1F("jesUncPt60", "jesUncPt60", 25, -0.05, 2.45);
    histo1D["jesUncPt70"] = new TH1F("jesUncPt70", "jesUncPt70", 25, -0.05, 2.45);
    
    histo2D["jesUnc_NEW"] = new TH2F("jesUnc_NEW", "jesUnc_NEW", 101, 29.5, 130.5, 25, -0.05, 2.45);
    graph2D["jesUnc_graph2D_NEW"] = new TGraph2D(101*25);
    histo1D["jesUncPt30_NEW"] = new TH1F("jesUncPt30_NEW", "jesUncPt30_NEW", 25, -0.05, 2.45);
    histo1D["jesUncPt40_NEW"] = new TH1F("jesUncPt40_NEW", "jesUncPt40_NEW", 25, -0.05, 2.45);
    histo1D["jesUncPt50_NEW"] = new TH1F("jesUncPt50_NEW", "jesUncPt50_NEW", 25, -0.05, 2.45);
    histo1D["jesUncPt60_NEW"] = new TH1F("jesUncPt60_NEW", "jesUncPt60_NEW", 25, -0.05, 2.45);
    histo1D["jesUncPt70_NEW"] = new TH1F("jesUncPt70_NEW", "jesUncPt70_NEW", 25, -0.05, 2.45);
    
    histo2D["jesUnc_NEW_Exact4"] = new TH2F("jesUnc_NEW_Exact4", "jesUnc_NEW_Exact4", 101, 29.5, 130.5, 25, -0.05, 2.45);
    graph2D["jesUnc_graph2D_NEW_Exact4"] = new TGraph2D(101*25);
    histo1D["jesUncPt30_NEW_Exact4"] = new TH1F("jesUncPt30_NEW_Exact4", "jesUncPt30_NEW_Exact4", 25, -0.05, 2.45);
    histo1D["jesUncPt40_NEW_Exact4"] = new TH1F("jesUncPt40_NEW_Exact4", "jesUncPt40_NEW_Exact4", 25, -0.05, 2.45);
    histo1D["jesUncPt50_NEW_Exact4"] = new TH1F("jesUncPt50_NEW_Exact4", "jesUncPt50_NEW_Exact4", 25, -0.05, 2.45);
    histo1D["jesUncPt60_NEW_Exact4"] = new TH1F("jesUncPt60_NEW_Exact4", "jesUncPt60_NEW_Exact4", 25, -0.05, 2.45);
    histo1D["jesUncPt70_NEW_Exact4"] = new TH1F("jesUncPt70_NEW_Exact4", "jesUncPt70_NEW_Exact4", 25, -0.05, 2.45);
    
    int iP=0;
    for(float jetPt=30; jetPt<=130; jetPt++)
    {
      for(float jetEta=0; jetEta<2.5; jetEta+=0.1)
      {
        jecUnc->setJetEta(jetEta);
        jecUnc->setJetPt(jetPt);
        float uncPos = jecUnc->getUncertainty(true);
        jecUnc->setJetEta(jetEta);
        jecUnc->setJetPt(jetPt);
        float uncNeg = jecUnc->getUncertainty(false);
        histo2D["jesUnc"]->Fill(jetPt, jetEta, 0.5*(uncPos+uncNeg));
        graph2D["jesUnc_graph2D"]->SetPoint(iP, jetPt, jetEta, 0.5*(uncPos+uncNeg));
        if(jetPt == 30) histo1D["jesUncPt30"]->Fill(jetEta, 0.5*(uncPos+uncNeg));
        else if(jetPt == 40) histo1D["jesUncPt40"]->Fill(jetEta, 0.5*(uncPos+uncNeg));
        else if(jetPt == 50) histo1D["jesUncPt50"]->Fill(jetEta, 0.5*(uncPos+uncNeg));
        else if(jetPt == 60) histo1D["jesUncPt60"]->Fill(jetEta, 0.5*(uncPos+uncNeg));
        else if(jetPt == 70) histo1D["jesUncPt70"]->Fill(jetEta, 0.5*(uncPos+uncNeg));
        
        //new JESunc stuff
        float unc_NEW = 0.5*(jetTools->calculateJESUnc(jetEta, jetPt, "plus") + jetTools->calculateJESUnc(jetEta, jetPt, "minus"));
        histo2D["jesUnc_NEW"]->Fill(jetPt, jetEta, unc_NEW);
        graph2D["jesUnc_graph2D_NEW"]->SetPoint(iP, jetPt, jetEta, unc_NEW);
        if(jetPt == 30) histo1D["jesUncPt30_NEW"]->Fill(jetEta, unc_NEW);
        else if(jetPt == 40) histo1D["jesUncPt40_NEW"]->Fill(jetEta, unc_NEW);
        else if(jetPt == 50) histo1D["jesUncPt50_NEW"]->Fill(jetEta, unc_NEW);
        else if(jetPt == 60) histo1D["jesUncPt60_NEW"]->Fill(jetEta, unc_NEW);
        else if(jetPt == 70) histo1D["jesUncPt70_NEW"]->Fill(jetEta, unc_NEW);
        
        float unc_NEW_Exact4 = 0.5*(jetTools_Exact4->calculateJESUnc(jetEta, jetPt, "plus") + jetTools_Exact4->calculateJESUnc(jetEta, jetPt, "minus"));
        histo2D["jesUnc_NEW_Exact4"]->Fill(jetPt, jetEta, unc_NEW_Exact4);
        graph2D["jesUnc_graph2D_NEW_Exact4"]->SetPoint(iP, jetPt, jetEta, unc_NEW_Exact4);
        if(jetPt == 30) histo1D["jesUncPt30_NEW_Exact4"]->Fill(jetEta, unc_NEW_Exact4);
        else if(jetPt == 40) histo1D["jesUncPt40_NEW_Exact4"]->Fill(jetEta, unc_NEW_Exact4);
        else if(jetPt == 50) histo1D["jesUncPt50_NEW_Exact4"]->Fill(jetEta, unc_NEW_Exact4);
        else if(jetPt == 60) histo1D["jesUncPt60_NEW_Exact4"]->Fill(jetEta, unc_NEW_Exact4);
        else if(jetPt == 70) histo1D["jesUncPt70_NEW_Exact4"]->Fill(jetEta, unc_NEW_Exact4);
        iP++;
      }
    }
    
    if(true) // Make JES unc plots for 8 TeV
    {
      JetCorrectionUncertainty *jecUnc_PuDataMC = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "PileUpDataMC")));
      JetCorrectionUncertainty *jecUnc_PuBias = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "PileUpBias")));
      
      JetCorrectionUncertainty *jecUnc_Rel = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "SubTotalRelative")));
      
      JetCorrectionUncertainty *jecUnc_Time = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "Time")));
      
      JetCorrectionUncertainty *jecUnc_Abs = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "Absolute")));
      JetCorrectionUncertainty *jecUnc_HighPtExtra = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "HighPtExtra")));
      JetCorrectionUncertainty *jecUnc_SinglePiECAL = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "SinglePionECAL")));
      JetCorrectionUncertainty *jecUnc_SinglePiHCAL = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "SinglePionHCAL")));
      
      JetCorrectionUncertainty *jecUnc_gluon = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "FlavorPureGluon")));
      JetCorrectionUncertainty *jecUnc_uds = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "FlavorPureQuark")));
      JetCorrectionUncertainty *jecUnc_c = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "FlavorPureCharm")));
      JetCorrectionUncertainty *jecUnc_b = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt", "FlavorPureBottom")));
      
//      jetTools->getJECunc(eta, pt, jesUnc);
      int i;
      for(float eta = -2.5, i = 0; eta <= 2.5; eta += 0.001, i++)
      {
        float uncPU = sqrt( pow( jetTools->getJECunc(eta, 30, jecUnc_PuDataMC) ,2) + pow( jetTools->getJECunc(eta, 30, jecUnc_PuBias) ,2) );
        float uncRel = jetTools->getJECunc(eta, 30, jecUnc_Rel);
        float uncTime = jetTools->getJECunc(eta, 30, jecUnc_Time);
        float uncAbs = sqrt( pow( jetTools->getJECunc(eta, 30, jecUnc_Abs) ,2) + pow( jetTools->getJECunc(eta, 30, jecUnc_HighPtExtra) ,2) + pow( jetTools->getJECunc(eta, 30, jecUnc_SinglePiECAL) ,2) + pow( jetTools->getJECunc(eta, 30, jecUnc_SinglePiHCAL) ,2) );
        float uncFlav = jetTools->calculateJESUnc(eta, 30, "plus", "onlyFlavour");
        float uncTotal = sqrt( uncPU*uncPU + uncRel*uncRel + uncTime*uncTime + uncAbs*uncAbs + uncFlav*uncFlav );
        
        graph["JecUnc_Pt30_PileUp"]->SetPoint(i, eta, uncPU);
        graph["JecUnc_Pt30_Relative"]->SetPoint(i, eta, uncRel);
        graph["JecUnc_Pt30_Time"]->SetPoint(i, eta, uncTime);
        graph["JecUnc_Pt30_Flavour"]->SetPoint(i, eta, uncFlav);
        graph["JecUnc_Pt30_Absolute"]->SetPoint(i, eta, uncAbs);
        graph["JecUnc_Pt30_Total"]->SetPoint(i, eta, uncTotal);
        
        uncPU = sqrt( pow( jetTools->getJECunc(eta, 50, jecUnc_PuDataMC) ,2) + pow( jetTools->getJECunc(eta, 50, jecUnc_PuBias) ,2) );
        uncRel = jetTools->getJECunc(eta, 50, jecUnc_Rel);
        uncTime = jetTools->getJECunc(eta, 50, jecUnc_Time);
        uncAbs = sqrt( pow( jetTools->getJECunc(eta, 50, jecUnc_Abs) ,2) + pow( jetTools->getJECunc(eta, 50, jecUnc_HighPtExtra) ,2) + pow( jetTools->getJECunc(eta, 50, jecUnc_SinglePiECAL) ,2) + pow( jetTools->getJECunc(eta, 50, jecUnc_SinglePiHCAL) ,2) );
        uncFlav = jetTools->calculateJESUnc(eta, 50, "plus", "onlyFlavour");
        uncTotal = sqrt( uncPU*uncPU + uncRel*uncRel + uncTime*uncTime + uncAbs*uncAbs + uncFlav*uncFlav );
        
        graph["JecUnc_Pt50_PileUp"]->SetPoint(i, eta, uncPU);
        graph["JecUnc_Pt50_Relative"]->SetPoint(i, eta, uncRel);
        graph["JecUnc_Pt50_Time"]->SetPoint(i, eta, uncTime);
        graph["JecUnc_Pt50_Flavour"]->SetPoint(i, eta, uncFlav);
        graph["JecUnc_Pt50_Absolute"]->SetPoint(i, eta, uncAbs);
        graph["JecUnc_Pt50_Total"]->SetPoint(i, eta, uncTotal);
        
        uncPU = sqrt( pow( jetTools->getJECunc(eta, 70, jecUnc_PuDataMC) ,2) + pow( jetTools->getJECunc(eta, 70, jecUnc_PuBias) ,2) );
        uncRel = jetTools->getJECunc(eta, 70, jecUnc_Rel);
        uncTime = jetTools->getJECunc(eta, 70, jecUnc_Time);
        uncAbs = sqrt( pow( jetTools->getJECunc(eta, 70, jecUnc_Abs) ,2) + pow( jetTools->getJECunc(eta, 70, jecUnc_HighPtExtra) ,2) + pow( jetTools->getJECunc(eta, 70, jecUnc_SinglePiECAL) ,2) + pow( jetTools->getJECunc(eta, 70, jecUnc_SinglePiHCAL) ,2) );
        uncFlav = jetTools->calculateJESUnc(eta, 70, "plus", "onlyFlavour");
        uncTotal = sqrt( uncPU*uncPU + uncRel*uncRel + uncTime*uncTime + uncAbs*uncAbs + uncFlav*uncFlav );
        
        graph["JecUnc_Pt70_PileUp"]->SetPoint(i, eta, uncPU);
        graph["JecUnc_Pt70_Relative"]->SetPoint(i, eta, uncRel);
        graph["JecUnc_Pt70_Time"]->SetPoint(i, eta, uncTime);
        graph["JecUnc_Pt70_Flavour"]->SetPoint(i, eta, uncFlav);
        graph["JecUnc_Pt70_Absolute"]->SetPoint(i, eta, uncAbs);
        graph["JecUnc_Pt70_Total"]->SetPoint(i, eta, uncTotal);
      }
      
      for(float pt = 20, i=0; pt <= 500; pt += 0.25, i++)
      {
        float uncPU = sqrt( pow( jetTools->getJECunc(0, pt, jecUnc_PuDataMC) ,2) + pow( jetTools->getJECunc(0, pt, jecUnc_PuBias) ,2) );
        float uncRel = jetTools->getJECunc(0, pt, jecUnc_Rel);
        float uncTime = jetTools->getJECunc(0, pt, jecUnc_Time);
        float uncAbs = sqrt( pow( jetTools->getJECunc(0, pt, jecUnc_Abs) ,2) + pow( jetTools->getJECunc(0, pt, jecUnc_HighPtExtra) ,2) + pow( jetTools->getJECunc(0, pt, jecUnc_SinglePiECAL) ,2) + pow( jetTools->getJECunc(0, pt, jecUnc_SinglePiHCAL) ,2) );
        float uncFlav = jetTools->calculateJESUnc(0, pt, "plus", "onlyFlavour");
        float uncTotal = sqrt( uncPU*uncPU + uncRel*uncRel + uncTime*uncTime + uncAbs*uncAbs + uncFlav*uncFlav );
        
        graph["JecUnc_Eta0_PileUp"]->SetPoint(i, pt, uncPU);
        graph["JecUnc_Eta0_Relative"]->SetPoint(i, pt, uncRel);
        graph["JecUnc_Eta0_Time"]->SetPoint(i, pt, uncTime);
        graph["JecUnc_Eta0_Flavour"]->SetPoint(i, pt, uncFlav);
        graph["JecUnc_Eta0_Absolute"]->SetPoint(i, pt, uncAbs);
        graph["JecUnc_Eta0_Total"]->SetPoint(i, pt, uncTotal);
        
        uncPU = sqrt( pow( jetTools->getJECunc(1.2, pt, jecUnc_PuDataMC) ,2) + pow( jetTools->getJECunc(1.2, pt, jecUnc_PuBias) ,2) );
        uncRel = jetTools->getJECunc(1.2, pt, jecUnc_Rel);
        uncTime = jetTools->getJECunc(1.2, pt, jecUnc_Time);
        uncAbs = sqrt( pow( jetTools->getJECunc(1.2, pt, jecUnc_Abs) ,2) + pow( jetTools->getJECunc(1.2, pt, jecUnc_HighPtExtra) ,2) + pow( jetTools->getJECunc(1.2, pt, jecUnc_SinglePiECAL) ,2) + pow( jetTools->getJECunc(1.2, pt, jecUnc_SinglePiHCAL) ,2) );
        uncFlav = jetTools->calculateJESUnc(1.2, pt, "plus", "onlyFlavour");
        uncTotal = sqrt( uncPU*uncPU + uncRel*uncRel + uncTime*uncTime + uncAbs*uncAbs + uncFlav*uncFlav );
        
        graph["JecUnc_Eta1p2_PileUp"]->SetPoint(i, pt, uncPU);
        graph["JecUnc_Eta1p2_Relative"]->SetPoint(i, pt, uncRel);
        graph["JecUnc_Eta1p2_Time"]->SetPoint(i, pt, uncTime);
        graph["JecUnc_Eta1p2_Flavour"]->SetPoint(i, pt, uncFlav);
        graph["JecUnc_Eta1p2_Absolute"]->SetPoint(i, pt, uncAbs);
        graph["JecUnc_Eta1p2_Total"]->SetPoint(i, pt, uncTotal);
        
        uncPU = sqrt( pow( jetTools->getJECunc(2.3, pt, jecUnc_PuDataMC) ,2) + pow( jetTools->getJECunc(2.3, pt, jecUnc_PuBias) ,2) );
        uncRel = jetTools->getJECunc(2.3, pt, jecUnc_Rel);
        uncTime = jetTools->getJECunc(2.3, pt, jecUnc_Time);
        uncAbs = sqrt( pow( jetTools->getJECunc(2.3, pt, jecUnc_Abs) ,2) + pow( jetTools->getJECunc(2.3, pt, jecUnc_HighPtExtra) ,2) + pow( jetTools->getJECunc(2.3, pt, jecUnc_SinglePiECAL) ,2) + pow( jetTools->getJECunc(2.3, pt, jecUnc_SinglePiHCAL) ,2) );
        uncFlav = jetTools->calculateJESUnc(2.3, pt, "plus", "onlyFlavour");
        uncTotal = sqrt( uncPU*uncPU + uncRel*uncRel + uncTime*uncTime + uncAbs*uncAbs + uncFlav*uncFlav );
        
        graph["JecUnc_Eta2p3_PileUp"]->SetPoint(i, pt, uncPU);
        graph["JecUnc_Eta2p3_Relative"]->SetPoint(i, pt, uncRel);
        graph["JecUnc_Eta2p3_Time"]->SetPoint(i, pt, uncTime);
        graph["JecUnc_Eta2p3_Flavour"]->SetPoint(i, pt, uncFlav);
        graph["JecUnc_Eta2p3_Absolute"]->SetPoint(i, pt, uncAbs);
        graph["JecUnc_Eta2p3_Total"]->SetPoint(i, pt, uncTotal);
        
        graph["FlavUnc_Eta0_gluon"]->SetPoint(i, pt, fabs(jetTools->getJECunc(0, pt, jecUnc_gluon)));
        graph["FlavUnc_Eta0_uds"]->SetPoint(i, pt, fabs(jetTools->getJECunc(0, pt, jecUnc_uds)));
        graph["FlavUnc_Eta0_c"]->SetPoint(i, pt, fabs(jetTools->getJECunc(0, pt, jecUnc_c)));
        graph["FlavUnc_Eta0_b"]->SetPoint(i, pt, fabs(jetTools->getJECunc(0, pt, jecUnc_b)));
        
        graph["FlavUnc_Eta1p9_gluon"]->SetPoint(i, pt, fabs(jetTools->getJECunc(1.9, pt, jecUnc_gluon)));
        graph["FlavUnc_Eta1p9_uds"]->SetPoint(i, pt, fabs(jetTools->getJECunc(1.9, pt, jecUnc_uds)));
        graph["FlavUnc_Eta1p9_c"]->SetPoint(i, pt, fabs(jetTools->getJECunc(1.9, pt, jecUnc_c)));
        graph["FlavUnc_Eta1p9_b"]->SetPoint(i, pt, fabs(jetTools->getJECunc(1.9, pt, jecUnc_b)));
        
//        cout << "i:  " << i << "  pt: " << pt << "  unc: " << fabs(jetTools->getJECunc(2.3, pt, jecUnc_gluon)) << endl;
        
        graph["FlavUnc_Eta2p3_gluon"]->SetPoint(i, pt, fabs(jetTools->getJECunc(2.3, pt, jecUnc_gluon)));
        graph["FlavUnc_Eta2p3_uds"]->SetPoint(i, pt, fabs(jetTools->getJECunc(2.3, pt, jecUnc_uds)));
        graph["FlavUnc_Eta2p3_c"]->SetPoint(i, pt, fabs(jetTools->getJECunc(2.3, pt, jecUnc_c)));
        graph["FlavUnc_Eta2p3_b"]->SetPoint(i, pt, fabs(jetTools->getJECunc(2.3, pt, jecUnc_b)));
      }
    }
    
    ////////////////////////////////
    // LOAD THE FULLKINFIT OBJECT //
    ////////////////////////////////
    
    FullKinFit* kinFit = NULL;
    if( ! CalculateResolutions )
    {
      kinFit = new FullKinFit(datasets[d], resFitLightJets, resFitBJets, measureTopMass, measureTopMassDifference);
      kinFit->SetResFitL7(resFitLightJetsL7, resFitBJetsL7);
    }
    
    ////////////////////////////////////////////////////////////
    // CREATE OUTPUT FILE AND TTREE FOR STORAGE OF THE NTUPLE //
    ////////////////////////////////////////////////////////////
    
    string decayChannel;
    if(semiElectron && semiMuon) decayChannel = "SemiLep";
    else
    {
      if(semiMuon) decayChannel = "SemiMu";
      if(semiElectron) decayChannel = "SemiEl";
    }
    
    string monsterFileTitle = "Monsters/KinFit_Monsters_"+dataSetName+"_"+decayChannel+".root";
    if(measureTopMass)
      monsterFileTitle = "Monsters/KinFit_Monsters_TopMass_"+dataSetName+"_"+decayChannel+".root";
    else if(measureTopMassDifference)
      monsterFileTitle = "Monsters/KinFit_LightMonsters_TopMassDiff_"+dataSetName+"_"+decayChannel+".root";
    
    cout << "INFO: creating Monsters file "+monsterFileTitle << endl;
    
    TFile* MonsterFile = new TFile(monsterFileTitle.c_str(),"RECREATE");
    
    LightMonster* lightMonster = 0;
    
    TTree* MonsterTree = new TTree("MonsterTree","Tree containing the monsters");
    
    if(measureTopMassDifference)
      MonsterTree->Branch("TheLightMonster","LightMonster",&lightMonster,256000);

    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1, itriggerSemiEl = -1, previousRun = -1;
    float nBjets = 0, nBjetsBtag = 0, nNonBjets = 0, nNonBjetsBtag = 0;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
//    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    for (unsigned int ievt = 0; ievt < 10; ievt++)
    {
      nEvents[d]++;
      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
      
      //load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      vector<TRootGenJet*> genjets;
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
      {
        genjets = treeLoader.LoadGenJet(ievt);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
//      cout << "run: " << event->runId() << "  lumi: " << event->lumiBlockId() << "  event: " << event->eventId() << endl;
      
      // scale factors for the event
      float scaleFactor = 1.;
      float lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
        lumiWeight = 1;
      
      TRootGenEvent* genEvt = 0;
      // Load the GenEvent
      if(dataSetName.find("TT") == 0) genEvt = treeLoader.LoadGenEvent(ievt);
      
      // Clone the init_jets vector, otherwise the corrections will be removed
      for(unsigned int i=0; i<init_jets_corrected.size(); i++)
        if(init_jets_corrected[i]) delete init_jets_corrected[i];
      init_jets_corrected.clear();

      for(unsigned int i=0; i<init_jets.size(); i++)
        init_jets_corrected.push_back( (TRootJet*) init_jets[i]->Clone() );
      
      // check which file in the dataset it is to have the HLTInfo right
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename)
      {
      	previousFilename = currentFilename;
	      iFile++;
	      cout<<"File changed!!! => iFile = "<<iFile<<endl;
      }
      
      int currentRun = event->runId();
      if(previousRun != currentRun)
      {
        previousRun = currentRun;
        
        // semi-muon
        if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
        {
          if( event->runId() <= 190738 )
      	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);
      	  else if( event->runId() >= 190782 && event->runId() <= 193621 )
      	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v12"), currentRun, iFile);
      	  else if( event->runId() >= 193834 && event->runId() <= 196531 )
      	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
      	  else if( event->runId() >= 198049 && event->runId() <= 199608)
      	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v14"), currentRun, iFile);
      	  else if( event->runId() >= 199698 && event->runId() <= 208686)
      	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v15"), currentRun, iFile);
          else
            cout << "Unknown run for SemiMu HLTpath selection: " << event->runId() << endl;
          if( itriggerSemiMu == 9999 )
          {
            cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
            exit(-1);
          }
        }
        else itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun);
        
        // semi-electron
        if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
        {
          if( event->runId() <= 190738 )
      	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v8"), currentRun, iFile);
      	  else if( event->runId() >= 190782 && event->runId() <= 191411 )
	          itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v9"), currentRun, iFile);
      	  else if( event->runId() >= 191695 && event->runId() <= 196531)
      	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun, iFile);
      	  else if( event->runId() >= 198049 && event->runId() <= 208686)
	          itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v11"), currentRun, iFile);
          else
            cout << "Unknown run for SemiEl HLTpath selection: " << event->runId() << endl;
          if( itriggerSemiEl == 9999 )
          {
            cout << "itriggerSemiEl == 9999 for SemiEl HLTpath selection: " << event->runId() << endl;
            exit(-1);
          }
        }
        else itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun);
      }
      // Apply Jet Corrections on-the-fly
      if( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
      {
        jetTools->unCorrectMETTypeOne(init_jets_corrected, mets[0], true);
        jetTools->correctJets(init_jets_corrected, event->kt6PFJets_rho(), true);
        jetTools->correctMETTypeOne(init_jets_corrected, mets[0], true);
      }
      else
      {
        jetTools->unCorrectMETTypeOne(init_jets_corrected, mets[0], false);
        jetTools->correctJets(init_jets_corrected, event->kt6PFJets_rho(), false);
        jetTools->correctMETTypeOne(init_jets_corrected, mets[0], false);
      }
      
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
      {
        if(systematic == "JERPlus")
          jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus", false);
        else if(systematic == "JERMinus")
          jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus", false);
        else
          jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal", false);
        
        // Correct jets for JES uncertainy systematics
        if(systematic == "JESPlus")
          jetTools->correctJetJESUnc(init_jets_corrected, mets[0], "plus");
        else if(systematic == "JESMinus")
          jetTools->correctJetJESUnc(init_jets_corrected, mets[0], "minus");
        
        // Make JES correction plots, first match PFJets with genJets
        vector< pair<size_t, size_t> > indexVector; //first index = RecoJet, second index = GenJet
        vector<bool> mLock(genjets.size(),false);   // when locked, genJet is already matched to a recoJet
        for(size_t i=0; i<init_jets_corrected.size(); i++)
        {
          pair<size_t, size_t> tmpIndex;
          float minDR = 9999.;
          for(size_t j=0; j<genjets.size(); j++)
          {
            if( ! mLock[j] )
            {
              if( init_jets_corrected[i]->DeltaR(*genjets[j]) < 0.4 && init_jets_corrected[i]->DeltaR(*genjets[j]) < minDR )
              {
                minDR = init_jets_corrected[i]->DeltaR(*genjets[j]);
                tmpIndex = pair<size_t, size_t>(i,j);
              }
            }
          }
          if(minDR < 999.)
          {
            mLock[tmpIndex.second] = true;
            indexVector.push_back(tmpIndex);
          }
        }
        
        // Nr. of good PV
        int nPV = 0;
        for(size_t i=0; i<vertex.size(); i++)
        {
          float Rho = sqrt(vertex[i]->x()*vertex[i]->x()+vertex[i]->y()*vertex[i]->y());
          if(!vertex[i]->isFake() && vertex[i]->ndof()>4 && abs(vertex[0]->z())<24 && Rho<2)
            nPV++;
        }
        
        // Now loop over matched combis and fill the plots
        for(size_t i=0; i<indexVector.size(); i++)
        {
          // First some preselection
          if( fabs(init_jets_corrected[indexVector[i].first]->Eta()) > 2.4 ) continue;
          
          // inJets[indexVector[i].first], inGenJets[indexVector[i].second]
          float PtRaw = init_jets_corrected[indexVector[i].first]->Pt() / init_jets_corrected[indexVector[i].first]->getJetCorrFactor("L1FastJetL2L3");
          float PtL1 = PtRaw * init_jets_corrected[indexVector[i].first]->getJetCorrFactor("L1FastJet");
          float responseRaw = PtRaw / genjets[indexVector[i].second]->Pt();
          float responseL1 = PtL1 / genjets[indexVector[i].second]->Pt();
          float response = init_jets_corrected[indexVector[i].first]->Pt() / genjets[indexVector[i].second]->Pt();
          
          int partFlav = init_jets_corrected[indexVector[i].first]->partonFlavour();
          
          if(response > 0 && response < 2)  //TODO: Change to 50 bins within 0 and 2 ??
          {
            profile["JetResp_L1L2L3_Incl"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            float eta = init_jets_corrected[indexVector[i].first]->Eta();
            if( fabs(eta) < 0.8 )
              profile["JetResp_L1L2L3_smallEta_Incl"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            else if( fabs(eta) < 1.6 )
              profile["JetResp_L1L2L3_mediumEta_Incl"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            else if( fabs(eta) < 2.4 )
              profile["JetResp_L1L2L3_largeEta_Incl"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            
            if(nPV <= 10)
            {
              profile["JetResp_Raw_0to10PV"]->Fill(genjets[indexVector[i].second]->Pt(), responseRaw);
              profile["JetResp_L1_0to10PV"]->Fill(genjets[indexVector[i].second]->Pt(), responseL1);
              profile["JetResp_L1L2L3_0to10PV"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
            else if(nPV <= 15)
            {
              profile["JetResp_Raw_10to15PV"]->Fill(genjets[indexVector[i].second]->Pt(), responseRaw);
              profile["JetResp_L1_10to15PV"]->Fill(genjets[indexVector[i].second]->Pt(), responseL1);
              profile["JetResp_L1L2L3_10to15PV"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
            else if(nPV <= 20)
            {
              profile["JetResp_Raw_15to20PV"]->Fill(genjets[indexVector[i].second]->Pt(), responseRaw);
              profile["JetResp_L1_15to20PV"]->Fill(genjets[indexVector[i].second]->Pt(), responseL1);
              profile["JetResp_L1L2L3_15to20PV"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
            else if(nPV <= 25)
            {
              profile["JetResp_Raw_20to25PV"]->Fill(genjets[indexVector[i].second]->Pt(), responseRaw);
              profile["JetResp_L1_20to25PV"]->Fill(genjets[indexVector[i].second]->Pt(), responseL1);
              profile["JetResp_L1L2L3_20to25PV"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
            else if(nPV <= 30)
            {
              profile["JetResp_Raw_25to30PV"]->Fill(genjets[indexVector[i].second]->Pt(), responseRaw);
              profile["JetResp_L1_25to30PV"]->Fill(genjets[indexVector[i].second]->Pt(), responseL1);
              profile["JetResp_L1L2L3_25to30PV"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
            else
            {
              profile["JetResp_Raw_30PlusPV"]->Fill(genjets[indexVector[i].second]->Pt(), responseRaw);
              profile["JetResp_L1_30PlusPV"]->Fill(genjets[indexVector[i].second]->Pt(), responseL1);
              profile["JetResp_L1L2L3_30PlusPV"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
            
            if( partFlav == 21 )
            {
              profile["JetResp_L1L2L3_gluon"]->Fill(genjets[indexVector[i].second]->Pt(), response);
              if( fabs(init_jets_corrected[indexVector[i].first]->Eta()) < 1.3 )
                profile["JetResp_L1L2L3_smallEta_gluon"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
            else if( fabs(partFlav) != 0 && fabs(partFlav) < 4 )
            {
              profile["JetResp_L1L2L3_uds"]->Fill(genjets[indexVector[i].second]->Pt(), response);
              if( fabs(init_jets_corrected[indexVector[i].first]->Eta()) < 1.3 )
                profile["JetResp_L1L2L3_smallEta_uds"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
            else if( fabs(partFlav) == 4 )
            {
              profile["JetResp_L1L2L3_c"]->Fill(genjets[indexVector[i].second]->Pt(), response);
              if( fabs(init_jets_corrected[indexVector[i].first]->Eta()) < 1.3 )
                profile["JetResp_L1L2L3_smallEta_c"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
            else if( fabs(partFlav) == 5 )
            {
              profile["JetResp_L1L2L3_b"]->Fill(genjets[indexVector[i].second]->Pt(), response);
              if( fabs(init_jets_corrected[indexVector[i].first]->Eta()) < 1.3 )
                profile["JetResp_L1L2L3_smallEta_b"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
            else
            {
              profile["JetResp_L1L2L3_unm"]->Fill(genjets[indexVector[i].second]->Pt(), response);
              if( fabs(init_jets_corrected[indexVector[i].first]->Eta()) < 1.3 )
                profile["JetResp_L1L2L3_smallEta_unm"]->Fill(genjets[indexVector[i].second]->Pt(), response);
            }
          }
        }
        
        // Make some b-tagging plots
        for(size_t i=0; i<init_jets_corrected.size(); i++)
        {
          if( fabs(init_jets_corrected[i]->Eta()) >= 2.4 || init_jets_corrected[i]->Pt() <= 30 )
            continue;
          
          int partFlav = init_jets_corrected[i]->partonFlavour();
          if(fabs(partFlav) == 5)
            histo1D["CSV_b"]->Fill( init_jets_corrected[i]->btag_combinedSecondaryVertexBJetTags() );
          else if(fabs(partFlav) == 4)
            histo1D["CSV_c"]->Fill( init_jets_corrected[i]->btag_combinedSecondaryVertexBJetTags() );
          else if(fabs(partFlav) == 3 || fabs(partFlav) == 2 || fabs(partFlav) == 1)
            histo1D["CSV_uds"]->Fill( init_jets_corrected[i]->btag_combinedSecondaryVertexBJetTags() );
          else if(fabs(partFlav) == 21)
            histo1D["CSV_g"]->Fill( init_jets_corrected[i]->btag_combinedSecondaryVertexBJetTags() );
          else
            histo1D["CSV_unm"]->Fill( init_jets_corrected[i]->btag_combinedSecondaryVertexBJetTags() );
        }
        
        
        // Do matching with leptons
        int muonMatchIndex = -1, electronMatchIndex = -1;
        vector<TRootMCParticle*> mcParticles;
        mcParticles = treeLoader.LoadMCPart(ievt);
        for(size_t i=0; i<mcParticles.size(); i++)
        {
          if(mcParticles[i]->status() == 3 && fabs(mcParticles[i]->type()) == 11 && fabs(mcParticles[i]->motherType()) == 24)
          {
            // Look for the closest electron
            float minDr = 9999;
            for(size_t j=0; j<init_electrons.size(); j++)
            {
              float DR = init_electrons[j]->DeltaR(*mcParticles[i]);
              if(DR < 0.1 && DR < minDr)
              {
                minDr = DR;
                electronMatchIndex = j;
              }
            }
          }
          else if(mcParticles[i]->status() == 3 && fabs(mcParticles[i]->type()) == 13 && fabs(mcParticles[i]->motherType()) == 24)
          {
            // Look for the closes muon
            float minDr = 9999;
            for(size_t j=0; j<init_muons.size(); j++)
            {
              float DR = init_muons[j]->DeltaR(*mcParticles[i]);
              if(DR < 0.1 && DR < minDr)
              {
                minDr = DR;
                muonMatchIndex = j;
              }
            }
          }
        }
        
        // Make electron Id plots
        for(size_t j=0; j<init_electrons.size(); j++)
        {
          double EffectiveArea = 0.;
          double SCeta = fabs(init_electrons[j]->superClusterEta());
          // HCP 2012 updated for electron conesize = 0.3, taken from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h?revision=1.4&view=markup
          if      (SCeta >= 0.0   && SCeta < 1.0   ) EffectiveArea = 0.130;
          else if (SCeta >= 1.0   && SCeta < 1.479 ) EffectiveArea = 0.137;
          else if (SCeta >= 1.479 && SCeta < 2.0   ) EffectiveArea = 0.067;
          else if (SCeta >= 2.0   && SCeta < 2.2   ) EffectiveArea = 0.089;
          else if (SCeta >= 2.2   && SCeta < 2.3   ) EffectiveArea = 0.107;
          else if (SCeta >= 2.3   && SCeta < 2.4   ) EffectiveArea = 0.110;
          else if (SCeta >= 2.4) EffectiveArea = 0.138;
          float reliso = (init_electrons[j]->chargedHadronIso() + max( init_electrons[j]->neutralHadronIso() + init_electrons[j]->photonIso()  - event->kt6PFJets_rho()*EffectiveArea, 0.) )/ init_electrons[j]->Pt(); // EffectiveArea corrected
          
          if(j == electronMatchIndex)
          {
            histo1D["ElIdMVA_matched"]->Fill(init_electrons[j]->mvaTrigId());
            histo1D["ElIso_matched"]->Fill(reliso);
          }
          else
          {
            histo1D["ElIdMVA_unmatched"]->Fill(init_electrons[j]->mvaTrigId());
            histo1D["ElIso_unmatched"]->Fill(reliso);
          }
        }
        
        // Make mu Id plots
        for(size_t j=0; j<init_muons.size(); j++)
        {
          float reliso = (init_muons[j]->chargedHadronIso() + max( 0.0, init_muons[j]->neutralHadronIso() + init_muons[j]->photonIso() - 0.5*init_muons[j]->puChargedHadronIso() ) ) / init_muons[j]->Pt(); // dBeta corrected
          
          if(j == muonMatchIndex)
          {
            histo1D["MuIso_matched"]->Fill( reliso );
            histo1D["MuNPixelHits_matched"]->Fill( init_muons[j]->nofValidPixelHits() );
            histo1D["Mud0_matched"]->Fill( init_muons[j]->d0() );
            
            histo1D["MuChi2_matched"]->Fill( init_muons[j]->chi2() );
            histo1D["MuDZ_matched"]->Fill( init_muons[j]->dz() );
            histo1D["MuNsegm_matched"]->Fill( init_muons[j]->nofMatchedStations() );
            histo1D["MuNTrackLayers_matched"]->Fill( init_muons[j]->nofTrackerLayersWithMeasurement() );
          }
          else
          {
            histo1D["MuIso_unmatched"]->Fill( reliso );
            histo1D["MuNPixelHits_unmatched"]->Fill( init_muons[j]->nofValidPixelHits() );
            histo1D["Mud0_unmatched"]->Fill( init_muons[j]->d0() );
            
            histo1D["MuChi2_unmatched"]->Fill( init_muons[j]->chi2() );
            histo1D["MuDZ_unmatched"]->Fill( init_muons[j]->dz() );
            histo1D["MuNsegm_unmatched"]->Fill( init_muons[j]->nofMatchedStations() );
            histo1D["MuNTrackLayers_unmatched"]->Fill( init_muons[j]->nofTrackerLayersWithMeasurement() );
          }
        }
      }
      
//      continue; // For plotting
      
      for(unsigned i=0; i<init_muons.size(); i++)
        MSPlot["AllMuonsRelPFIso"]->Fill((init_muons[i]->chargedHadronIso()+init_muons[i]->neutralHadronIso()+init_muons[i]->photonIso())/init_muons[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      for(unsigned i=0; i<init_electrons.size(); i++)
        MSPlot["AllElectronsRelPFIso"]->Fill((init_electrons[i]->chargedHadronIso()+init_electrons[i]->neutralHadronIso()+init_electrons[i]->photonIso())/init_electrons[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      
      if(systematic == "AlignPlus" || systematic == "AlignMinus")
      {
        // Scale PFJets up/down according to their charge (for tracker misalignment studies)
        for(unsigned int iJet=0; iJet<init_jets_corrected.size(); iJet++)
        {
          TRootPFJet* jet = jetTools->convertToPFJets(init_jets_corrected[iJet]);
          if(jet->chargedMultiplicity() > 0)
          {
            float chargedFraction = jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction() + jet->chargedMuEnergyFraction();
            float chargedAveragePt = jet->Pt() * chargedFraction / jet->chargedMultiplicity();
            float charge = jet->charge();
            float deltaPtFraction = ( chargedAveragePt * chargedAveragePt * charge * 0.0001 ) / jet->Pt();
            if( systematic == "AlignMinus" ) deltaPtFraction *= -1.;
            init_jets_corrected[iJet]->SetPxPyPzE(jet->Px()*(1+deltaPtFraction), jet->Py()*(1+deltaPtFraction), jet->Pz()*(1+deltaPtFraction), jet->E()*(1+deltaPtFraction));
            histo1D["AlignSystSF"]->Fill(deltaPtFraction);
          }
        }
      }
      
      if(systematic == "bVSbbarJES" || systematic == "bVSbbarJESPtDep")
      {
        // Scale PFJets up/down according to their flavour (b or bbar)
        for(unsigned int iJet=0; iJet<init_jets_corrected.size(); iJet++)
        {
          float factor = 0.000779 / 2.; // divide by two, to scale b-jets up by half of the effect, and bbar ets down by half of the effect
          if(systematic == "bVSbbarJESPtDep") // Pt dependent!
          {
            float jetPt = init_jets_corrected[iJet]->Pt();
            if(jetPt > 30 && jetPt < 40) factor = 0.00139575 / 2.;
            else if(jetPt > 40 && jetPt < 50) factor = 0.00319358 / 2.;
            else if(jetPt > 50 && jetPt < 60) factor = 0.001213 / 2.;
            else if(jetPt > 60 && jetPt < 70) factor = 0.00148661 / 2.;
            else if(jetPt > 70 && jetPt < 80) factor = 0.00125843 / 2.;
            else if(jetPt > 80 && jetPt < 90) factor = 0.00124397 / 2.;
            else if(jetPt > 90 && jetPt < 100) factor = -0.00253308 / 2.;
            else if(jetPt > 100 && jetPt < 120) factor = 0.000311997 / 2.;
            else if(jetPt > 120 && jetPt < 150) factor = -0.000796185 / 2.;
            else if(jetPt > 150 && jetPt < 250) factor = -0.000396841 / 2.;
            else if(jetPt > 250) factor = -0.00129253 / 2.;
          }
          
          if(init_jets_corrected[iJet]->partonFlavour() == 5)
            init_jets_corrected[iJet]->SetPxPyPzE(init_jets_corrected[iJet]->Px()*(1+factor), init_jets_corrected[iJet]->Py()*(1+factor), init_jets_corrected[iJet]->Pz()*(1+factor), init_jets_corrected[iJet]->E()*(1+factor));
          else if(init_jets_corrected[iJet]->partonFlavour() == -5)
            init_jets_corrected[iJet]->SetPxPyPzE(init_jets_corrected[iJet]->Px()*(1-factor), init_jets_corrected[iJet]->Py()*(1-factor), init_jets_corrected[iJet]->Pz()*(1-factor), init_jets_corrected[iJet]->E()*(1-factor));
        }
      }
      
      /////////////////////////////
      //   Selection
      /////////////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho());
      selection.setJetCuts(30,2.4,0.01,1.,0.98,0.3,0.1);
      selection.setMuonCuts(25,2.1,0.12,0.2,0.3,1,0.5,5,0); // DR mu-jets cleaning still needed?
      selection.setElectronCuts(32,2.5,0.1,0.02,0.5,0.3,0); // 32 GeV too low? DR el-jets cleaning still needed?
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setLooseElectronCuts(20,2.5,0.15,0.);
      
      bool triggedSemiMu = false;
      if( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("Data_SingleEl") == 0 ) )
        triggedSemiMu = treeLoader.EventTrigged (itriggerSemiMu);
      bool triggedSemiEl = false;
      if( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("Data_SingleMu") == 0 ) )
        triggedSemiEl = treeLoader.EventTrigged (itriggerSemiEl);
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);
      
      vector<TRootJet*> selectedJets, selectedJetsNoMu, selectedJetsNoEl;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons(10,2.5,0.2);
      vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(20,2.5,0.15);
      if( systematic == "InvertedIso" )
      {
        vector<TRootMuon*> overlapMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0]);
        vector<TRootElectron*> overlapElectrons = selection.GetSelectedElectronsInvIso(0.2);
        selectedJetsNoMu = selection.GetSelectedJets(overlapMuons,true);
        selectedJetsNoEl = selection.GetSelectedJets(overlapElectrons,true);
        selectedMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0], selectedJetsNoMu);
        selectedElectrons = selection.GetSelectedElectronsInvIso(0.2, selectedJetsNoEl); // Also mvaId < 0.5
      }
      else // Normal selection
      {
        selectedJets = selection.GetSelectedJets(true);
        selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
        selectedElectrons = selection.GetSelectedElectrons(selectedJets);
      }
      
      MSPlot["nPrimaryVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
   		MSPlot["nGoodPrimaryVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nSelectedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nSelectedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLooseOtherMuons"]->Fill(vetoMuons.size()-1, datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLooseElectrons"]->Fill(vetoElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      
      bool eventSelectedSemiMu = false;
      bool eventSelectedSemiEl = false;
      
      MSPlot["nEventsAfterCutsSemiMu"]->Fill(0, datasets[d], true, Luminosity*scaleFactor);
      selecTableSemiMu.Fill(d,0,scaleFactor*lumiWeight);
      if( triggedSemiMu && semiMuon )
      {
        MSPlot["nEventsAfterCutsSemiMu"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
        selecTableSemiMu.Fill(d,1,scaleFactor*lumiWeight);
        if( isGoodPV )
        {
          MSPlot["nEventsAfterCutsSemiMu"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
          selecTableSemiMu.Fill(d,2,scaleFactor*lumiWeight);
          if( selectedMuons.size() == 1 )
          {
            MSPlot["nEventsAfterCutsSemiMu"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
  		      selecTableSemiMu.Fill(d,3,scaleFactor*lumiWeight);
        		if( ( systematic != "InvertedIso" && vetoMuons.size() == 1 ) || ( systematic == "InvertedIso" && vetoMuons.size() == 0 ) ) // if InvertedIso, selected muon not part of vetoMuons vector!
        		{
        		  MSPlot["nEventsAfterCutsSemiMu"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
              selecTableSemiMu.Fill(d,4,scaleFactor*lumiWeight);
              if( vetoElectrons.size() == 0 )
              {
                MSPlot["nEventsAfterCutsSemiMu"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
                selecTableSemiMu.Fill(d,5,scaleFactor*lumiWeight);
                if(systematic == "InvertedIso") selectedJets = selectedJetsNoMu;
                if(selectedJets.size()>=1)
                {
                  MSPlot["nEventsAfterCutsSemiMu"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
                  selecTableSemiMu.Fill(d,6,scaleFactor*lumiWeight);
                  if(selectedJets.size()>=2)
                  {
                    MSPlot["nEventsAfterCutsSemiMu"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
                    selecTableSemiMu.Fill(d,7,scaleFactor*lumiWeight);
                    if(selectedJets.size()>=3)
                    {
                      MSPlot["nEventsAfterCutsSemiMu"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
                      selecTableSemiMu.Fill(d,8,scaleFactor*lumiWeight);
                      if(selectedJets.size()>=4)
                      {
                        histo1D["FourthJetPt"]->Fill(selectedJets[3]->Pt());
                        if(triggedSemiMu) histo1D["FourthJetPtTriggered"]->Fill(selectedJets[3]->Pt());
                        MSPlot["nEventsAfterCutsSemiMu"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
                        selecTableSemiMu.Fill(d,9,scaleFactor*lumiWeight);
                        eventSelectedSemiMu = true;
                        float reliso = (selectedMuons[0]->chargedHadronIso()+selectedMuons[0]->neutralHadronIso()+selectedMuons[0]->photonIso())/selectedMuons[0]->Pt();
                        MSPlot["SelectedEventsMuonsRelPFIso"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      
      MSPlot["nEventsAfterCutsSemiEl"]->Fill(0, datasets[d], true, Luminosity*scaleFactor);
      selecTableSemiEl.Fill(d,0,scaleFactor*lumiWeight);
      if( semiElectron && triggedSemiEl )
      {
        MSPlot["nEventsAfterCutsSemiEl"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
        selecTableSemiEl.Fill(d,1,scaleFactor*lumiWeight);
        if( isGoodPV )
        {
          MSPlot["nEventsAfterCutsSemiEl"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
          selecTableSemiEl.Fill(d,2,scaleFactor*lumiWeight);
          if( selectedElectrons.size() == 1 )
          {
            MSPlot["nEventsAfterCutsSemiEl"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
            selecTableSemiEl.Fill(d,3,scaleFactor*lumiWeight);
            if( vetoMuons.size() == 0 )
            {
              MSPlot["nEventsAfterCutsSemiEl"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
              selecTableSemiEl.Fill(d,4,scaleFactor*lumiWeight);
              if( ( systematic != "InvertedIso" && vetoElectrons.size() == 1 ) || ( systematic == "InvertedIso" && vetoElectrons.size() == 0 ) ) // if InvertedIso, selected electron not part of vetoElectrons vector!
              {
                MSPlot["nEventsAfterCutsSemiEl"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
                selecTableSemiEl.Fill(d,5,scaleFactor*lumiWeight);
                if(systematic == "InvertedIso") selectedJets = selectedJetsNoEl;
                if( selectedJets.size()>=1 )
                {
                  MSPlot["nEventsAfterCutsSemiEl"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
                  selecTableSemiEl.Fill(d,6,scaleFactor*lumiWeight);
                  if( selectedJets.size()>=2 )
                  {
                    MSPlot["nEventsAfterCutsSemiEl"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
                    selecTableSemiEl.Fill(d,7,scaleFactor*lumiWeight);
                    if( selectedJets.size()>=3 )
                    {
                      MSPlot["nEventsAfterCutsSemiEl"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
                      selecTableSemiEl.Fill(d,8,scaleFactor*lumiWeight);
                      if( selectedJets.size()>=4 )
                      {
                        histo1D["FourthJetPt"]->Fill(selectedJets[3]->Pt());
                        if(triggedSemiEl) histo1D["FourthJetPtTriggered"]->Fill(selectedJets[3]->Pt());
                        MSPlot["nEventsAfterCutsSemiEl"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
                        selecTableSemiEl.Fill(d,9,scaleFactor*lumiWeight);
                        eventSelectedSemiEl = true;
                        float reliso = (selectedElectrons[0]->chargedHadronIso()+selectedElectrons[0]->neutralHadronIso()+selectedElectrons[0]->photonIso())/selectedElectrons[0]->Pt();
                        MSPlot["SelectedEventsElectronsRelPFIso"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      
      if( ! ( eventSelectedSemiEl || eventSelectedSemiMu ) ) continue;
      
      for(size_t i=0; i<selectedJets.size(); i++)
        histo1D["JetPtSelEv"]->Fill(selectedJets[i]->Pt());
      histo1D["METSelEv"]->Fill(mets[0]->Pt());
      
      if( eventSelectedSemiEl && eventSelectedSemiMu )
        cout << "Event selected in semiEl and semiMu channel???" << endl;
      
      vector<TRootMCParticle*> mcParticles;
      if( dataSetName.find("TTbarJets") == 0 || dataSetName.find("TT_") == 0 || dataSetName.find("ST_") == 0 || dataSetName.find("SingleTop") == 0 )
      {
        mcParticles = treeLoader.LoadMCPart(ievt);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      if( eventSelectedSemiMu )
        jetCombiner->ProcessEvent(datasets[d], mcParticles, selectedJets, selectedMuons[0], init_electrons, init_muons, genEvt, scaleFactor, false);
      else
        jetCombiner->ProcessEvent(datasets[d], mcParticles, selectedJets, selectedElectrons[0], init_electrons, init_muons, genEvt, scaleFactor, false);
      
      if(CalculateResolutions && dataSetName.find("TTbarJets") == 0)
      {
        jetCombiner->FillResolutions(resFitLightJets, resFitBJets, resFitBJets_B, resFitBJets_Bbar, scaleFactor*lumiWeight);
        float pTtop = 0, pTantiTop = 0;
        for(unsigned int iPart=0; iPart<mcParticles.size(); iPart++)
        {
          if( mcParticles[iPart]->status() != 3) continue;
          if( mcParticles[iPart]->type() == 6 )
          {
            pTtop = mcParticles[iPart]->Pt();
          }
          else if( mcParticles[iPart]->type() == -6 )
          {
            pTantiTop = mcParticles[iPart]->Pt();
          }
        }
        float pTtopWeight = sqrt( exp(0.156 - 0.00137*pTtop) * exp(0.156 - 0.00137*pTantiTop) );
        jetCombiner->FillResolutions(resFitLightJets_PtTopRew, resFitBJets_PtTopRew, 0, 0, scaleFactor*lumiWeight*pTtopWeight);
      }
      
      if ( !TrainMVA && !CalculateResolutions )
      {
        int nBtags = 0;
        vector<float> bTagCSV;
        vector<TLorentzVector> otherSelectedJets;
        vector<int> partonFlavour;
        double bTagCutValue = 0.679; // nominal: 0.679
        if(systematic == "bTagPlus") bTagCutValue = 0.6427; // absolute change in eff: 1.001 % 
        else if(systematic == "bTagMinus") bTagCutValue = 0.7155; // absolute change in eff: -1.003 %
        else if(systematic == "bVSbbarTag")
        {
          if(eventSelectedSemiMu)
          {
            if(selectedMuons[0]->charge() > 0) bTagCutValue = 0.6000; // abs change: 2.2727 %
            else bTagCutValue = 0.7532; // abs change:  -2.2725%
          }
          else
          {
            if(selectedElectrons[0]->charge() > 0) bTagCutValue = 0.6000;
            else bTagCutValue = 0.7532;
          }
        }
        
        for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++)
        {
          otherSelectedJets.push_back( *selectedJets[iJet] );
          bTagCSV.push_back(selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags());
          partonFlavour.push_back( selectedJets[iJet]->partonFlavour() );
          if( selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags() > bTagCutValue )
          {
            nBtags++;
            if( fabs(selectedJets[iJet]->partonFlavour()) == 5 ) nBjetsBtag++;
            else nNonBjetsBtag++;
          }
          if( fabs(selectedJets[iJet]->partonFlavour()) == 5 ) nBjets++;
          else nNonBjets++;
        }
        
        if(nBtags < 1) continue;
        
//        continue;
        
        //get the MC matched jet combination, not the MVA best matched
	      vector<unsigned int> mcJetCombi = jetCombiner->GetGoodJetCombination();
        int hadrBJetIndex = mcJetCombi[2], lightJet1Index = mcJetCombi[0], lightJet2Index = mcJetCombi[1], leptBJetIndex = mcJetCombi[3];
        
        if( mcJetCombi[0] < 9999 && mcJetCombi[1] < 9999 && mcJetCombi[2] < 9999 && mcJetCombi[3] < 9999 && (genEvt->isSemiLeptonic( TRootGenEvent::kMuon ) || genEvt->isSemiLeptonic( TRootGenEvent::kElec )) )
        {
          vector<TLorentzVector> matchedQuarks, genJetsTLV, matchedGenJets, selectedJetsTLV, matchedPFJetsGenJets;
          matchedQuarks.push_back(jetCombiner->GetLightQuark1());
          matchedQuarks.push_back(jetCombiner->GetLightQuark2());
          matchedQuarks.push_back(jetCombiner->GetHadrBQuark());
          matchedQuarks.push_back(jetCombiner->GetLeptBQuark());
          vector<float> dummyBtag (4, 0.);
          for(size_t i=0; i<genjets.size(); i++)
            genJetsTLV.push_back( *genjets[i] );
          jetCombinerGenPartGenJets->ProcessEvent(datasets[d], matchedQuarks, genJetsTLV, dummyBtag, TLorentzVector(), true, scaleFactor, false);
          vector<unsigned int> mcGenJetCombi = jetCombinerGenPartGenJets->GetGoodJetCombination(); // matching of genParticles with genJets
          
          if( mcGenJetCombi[0] < 9999 && mcGenJetCombi[1] < 9999 && mcGenJetCombi[2] < 9999 && mcGenJetCombi[3] < 9999 )
          {
            matchedGenJets.push_back(*genjets[mcGenJetCombi[0]]);
            matchedGenJets.push_back(*genjets[mcGenJetCombi[1]]);
            matchedGenJets.push_back(*genjets[mcGenJetCombi[2]]);
            matchedGenJets.push_back(*genjets[mcGenJetCombi[3]]);
            for(size_t i=0; i<selectedJets.size(); i++)
              selectedJetsTLV.push_back( *selectedJets[i] );
            jetCombinerGenJetsPFJets->ProcessEvent(datasets[d], matchedGenJets, selectedJetsTLV, dummyBtag, TLorentzVector(), true, scaleFactor, false);
            vector<unsigned int> mcGenJetPFJetCombi = jetCombinerGenJetsPFJets->GetGoodJetCombination(); // matching of PFJets with genJets
            
            if( mcGenJetPFJetCombi[0] < 9999 && mcGenJetPFJetCombi[1] < 9999 && mcGenJetPFJetCombi[2] < 9999 && mcGenJetPFJetCombi[3] < 9999 )
            {
              matchedPFJetsGenJets.push_back(selectedJetsTLV[mcGenJetPFJetCombi[0]]);
              matchedPFJetsGenJets.push_back(selectedJetsTLV[mcGenJetPFJetCombi[1]]);
              matchedPFJetsGenJets.push_back(selectedJetsTLV[mcGenJetPFJetCombi[2]]);
              matchedPFJetsGenJets.push_back(selectedJetsTLV[mcGenJetPFJetCombi[3]]);
              
              vector<TRootJet*> matchedJets;
              for(size_t i=0; i<4; i++)
                for(size_t j=0; j<selectedJets.size(); j++)
                  if( fabs(matchedPFJetsGenJets[i].Pt() - selectedJets[j]->Pt()) < 0.01 && fabs(matchedPFJetsGenJets[i].Eta() - selectedJets[j]->Eta()) < 0.01 )
                    matchedJets.push_back( selectedJets[j] );
              vector<TRootPFJet*> pfJets = jetTools->convertToPFJets(matchedJets);
              
              float lightCorr1 = 1 - resFitLightJets->EtCorrection(&matchedPFJetsGenJets[0]);
              float lightCorr2 = 1 - resFitLightJets->EtCorrection(&matchedPFJetsGenJets[1]);
              float bCorr = 1 - resFitBJets->EtCorrection(&matchedPFJetsGenJets[2]);
              
              if(matchedPFJetsGenJets[0].DeltaR(matchedQuarks[0]) < 0.3 && matchedPFJetsGenJets[1].DeltaR(matchedQuarks[1]) < 0.3 && matchedPFJetsGenJets[2].DeltaR(matchedQuarks[2]) < 0.3 && matchedPFJetsGenJets[3].DeltaR(matchedQuarks[3]) < 0.3)
              {
                // Now do something with the fully matched events.....
                float mTop = (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]).M();
                histo1D["mTop"]->Fill(mTop);
                float mTopGenJet = (matchedGenJets[0]+matchedGenJets[1]+matchedGenJets[2]).M();
                histo1D["mTop_genJet"]->Fill(mTopGenJet);
                
                float dRMin = matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[1] );
                if( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] ) < dRMin ) dRMin = matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] );
                if( matchedPFJetsGenJets[2].DeltaR( matchedPFJetsGenJets[1] ) < dRMin ) dRMin = matchedPFJetsGenJets[2].DeltaR( matchedPFJetsGenJets[1] );
                
                histo1D["dRMin"]->Fill(dRMin);
                histo2D["dRMin_VS_mTop"]->Fill(dRMin, mTop);
                histo2D["dRMin_VS_mTopGenJet"]->Fill(dRMin, mTopGenJet);
                histo1D["nJets"]->Fill(selectedJets.size());
                histo2D["nJets_VS_mTop"]->Fill(selectedJets.size(), mTop);
                histo2D["nJets_VS_mTopGenJet"]->Fill(selectedJets.size(), mTopGenJet);
                for(size_t i=0; i<matchedJets.size(); i++)
                  histo1D["jetArea"]->Fill(matchedJets[i]->jetArea());
                histo1D["jetAreaBjet"]->Fill( matchedJets[2]->jetArea() );
                histo2D["jetAreaBjet_VS_mTop"]->Fill( matchedJets[2]->jetArea(), mTop);
                histo2D["jetAreaBjet_VS_mTopGenJet"]->Fill( matchedJets[2]->jetArea(), mTopGenJet);
                histo1D["jetAreaTotal"]->Fill( matchedJets[0]->jetArea()+matchedJets[1]->jetArea()+matchedJets[2]->jetArea() );
                histo2D["jetAreaTotal_VS_mTop"]->Fill( matchedJets[0]->jetArea()+matchedJets[1]->jetArea()+matchedJets[2]->jetArea(), mTop);
                histo2D["jetAreaTotal_VS_mTopGenJet"]->Fill( matchedJets[0]->jetArea()+matchedJets[1]->jetArea()+matchedJets[2]->jetArea(), mTopGenJet);
                
                if( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[1] ) < 1. )
                  histo1D["mW_dRsmall"]->Fill( (matchedPFJetsGenJets[0]*lightCorr1+matchedPFJetsGenJets[1]*lightCorr2).M() );
                else
                  histo1D["mW_dRlarge"]->Fill( (matchedPFJetsGenJets[0]*lightCorr1+matchedPFJetsGenJets[1]*lightCorr2).M() );
                
                if(dRMin < 1.)
                  histo1D["mTop_dRsmall"]->Fill( (matchedPFJetsGenJets[0]*lightCorr1+matchedPFJetsGenJets[1]*lightCorr2+matchedPFJetsGenJets[2]*bCorr).M() );
                else
                  histo1D["mTop_dRlarge"]->Fill( (matchedPFJetsGenJets[0]*lightCorr1+matchedPFJetsGenJets[1]*lightCorr2+matchedPFJetsGenJets[2]*bCorr).M() );
                
  //              float mTopNoFitL7 = ( selectedJets[combi[0]]*lightCorr1 + selectedJets[combi[1]]*lightCorr2 + selectedJets[combi[2]]*bCorr ).M();
                
                float dRMinLight1 = matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[1] );
                if(matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] ) < dRMinLight1)
                  dRMinLight1 = matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] );
                
                float dRMinLight2 = matchedPFJetsGenJets[2].DeltaR( matchedPFJetsGenJets[1] );
                if(matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[1] ) < dRMinLight2)
                  dRMinLight2 = matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[1] );
                
                float dRMinHadrB = matchedPFJetsGenJets[2].DeltaR( matchedPFJetsGenJets[1] );
                if(matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] ) < dRMinHadrB)
                  dRMinHadrB = matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] );
                
                for(size_t i=0; i<nBinsJet; i++)
                {
                  stringstream ss; ss << i;
                  if( dRMinLight1 > binningDR[i] && dRMinLight1 <= binningDR[i+1] )
                    histo1D["lightJetResp_DR_"+ss.str()]->Fill( matchedPFJetsGenJets[0].Pt() / matchedQuarks[0].Pt() );
                  if( dRMinLight2 > binningDR[i] && dRMinLight2 <= binningDR[i+1] )
                    histo1D["lightJetResp_DR_"+ss.str()]->Fill( matchedPFJetsGenJets[1].Pt() / matchedQuarks[1].Pt() );
                  if( dRMinHadrB > binningDR[i] && dRMinHadrB <= binningDR[i+1] )
                    histo1D["bJetResp_DR_"+ss.str()]->Fill( matchedPFJetsGenJets[2].Pt() / matchedQuarks[2].Pt() );
                  
                  for(size_t j=0; j<3; j++)
                  {
                    if( matchedJets[j]->jetArea() > binningJetArea[i] && matchedJets[j]->jetArea() <= binningJetArea[i+1] )
                    {
                      if(j<2) histo1D["lightJetResp_Area_"+ss.str()]->Fill( matchedPFJetsGenJets[j].Pt() / matchedQuarks[j].Pt() );
                      else histo1D["bJetResp_Area_"+ss.str()]->Fill( matchedPFJetsGenJets[j].Pt() / matchedQuarks[j].Pt() );
                    }
                    float nPart = pfJets[j]->chargedMultiplicity()+pfJets[j]->neutralMultiplicity()+pfJets[j]->muonMultiplicity();
                    if( nPart > binningNpart[i] && nPart <= binningNpart[i+1] )
                    {
                      if(j<2) histo1D["lightJetResp_nPart_"+ss.str()]->Fill( matchedPFJetsGenJets[j].Pt() / matchedQuarks[j].Pt() );
                      else  histo1D["bJetResp_nPart_"+ss.str()]->Fill( matchedPFJetsGenJets[j].Pt() / matchedQuarks[j].Pt() );
                    }
                  }
                }
                
                for(size_t i=0; i<pfJets.size(); i++)
                {
                  histo1D["nParticlesInJet"]->Fill( pfJets[i]->nConstituents() );
                  histo1D["nChargedParticles"]->Fill(pfJets[i]->chargedMultiplicity());
                }
                histo1D["nParticlesBJet"]->Fill( pfJets[2]->nConstituents() );
                histo2D["nParticlesBJet_VS_mTop"]->Fill( pfJets[2]->nConstituents(), mTop);
                histo2D["nParticlesBJet_VS_mTopGenJet"]->Fill( pfJets[2]->nConstituents(), mTopGenJet);
                histo1D["nParticlesTotal"]->Fill( pfJets[0]->nConstituents()+pfJets[1]->nConstituents()+pfJets[2]->nConstituents() );
                histo2D["nParticlesTotal_VS_mTop"]->Fill( pfJets[0]->nConstituents()+pfJets[1]->nConstituents()+pfJets[2]->nConstituents(), mTop);
                histo2D["nParticlesTotal_VS_mTopGenJet"]->Fill( pfJets[0]->nConstituents()+pfJets[1]->nConstituents()+pfJets[2]->nConstituents(), mTopGenJet);
                histo1D["nChargedBJet"]->Fill( pfJets[2]->chargedMultiplicity() );
                histo2D["nChargedBJet_VS_mTop"]->Fill( pfJets[2]->chargedMultiplicity(), mTop);
                histo2D["nChargedBJet_VS_mTopGenJet"]->Fill( pfJets[2]->chargedMultiplicity(), mTopGenJet);
                histo1D["nChargedTotal"]->Fill( pfJets[0]->chargedMultiplicity()+pfJets[1]->chargedMultiplicity()+pfJets[2]->chargedMultiplicity() );
                histo2D["nChargedTotal_VS_mTop"]->Fill( pfJets[0]->chargedMultiplicity()+pfJets[1]->chargedMultiplicity()+pfJets[2]->chargedMultiplicity(), mTop);
                histo2D["nChargedTotal_VS_mTopGenJet"]->Fill( pfJets[0]->chargedMultiplicity()+pfJets[1]->chargedMultiplicity()+pfJets[2]->chargedMultiplicity(), mTopGenJet);
                
                histo1D["dRLights"]->Fill( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[1] ) );
                histo2D["dRLights_VS_mTop"]->Fill( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[1] ), mTop);
                histo2D["dRLights_VS_mTopGenJet"]->Fill( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[1] ), mTopGenJet);
                if( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] ) < matchedPFJetsGenJets[2].DeltaR( matchedPFJetsGenJets[1] ) )
                {
                  histo1D["MinDRLightB"]->Fill( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] ) );
                  histo2D["MinDRLightB_VS_mTop"]->Fill( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] ), mTop);
                  histo2D["MinDRLightB_VS_mTopGenJet"]->Fill( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] ), mTopGenJet);
                }
                else
                {
                  histo1D["MinDRLightB"]->Fill( matchedPFJetsGenJets[2].DeltaR( matchedPFJetsGenJets[1] ) );
                  histo2D["MinDRLightB_VS_mTop"]->Fill( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] ), mTop);
                  histo2D["MinDRLightB_VS_mTopGenJet"]->Fill( matchedPFJetsGenJets[0].DeltaR( matchedPFJetsGenJets[2] ), mTopGenJet);
                }
                histo1D["MET"]->Fill(mets[0]->Pt());
                histo2D["MET_VS_mTop"]->Fill(mets[0]->Pt(), mTop);
                histo2D["MET_VS_mTopGenJet"]->Fill(mets[0]->Pt(), mTopGenJet);
                float ht = 0;
                for(size_t i=0; i<selectedJets.size(); i++)
                  ht += selectedJets[i]->Pt();
                histo1D["HT"]->Fill(ht);
                histo2D["HT_VS_mTop"]->Fill(ht, mTop);
                histo2D["HT_VS_mTopGenJet"]->Fill(ht, mTopGenJet);
                
                //calculate pz neutrino (for mttbar)
                TLorentzVector lepton;
                if(eventSelectedSemiMu) lepton = *selectedMuons[0];
                else lepton = *selectedElectrons[0];
                
                double M_W  = 80.4;
                double M_mu =  0.10566;
                double pznu = 0.;
                
                double a = M_W*M_W - M_mu*M_mu + 2.0*lepton.Px()*mets[0]->Px() + 2.0*lepton.Py()*mets[0]->Py();
                double A = 4.0*(pow(lepton.E(),2)- pow(lepton.Pz(),2));
                double B = -4.0*a*lepton.Pz();
                double C = 4.0*pow(lepton.E(),2)*(pow(mets[0]->Px(),2) + pow(mets[0]->Py(),2)) - a*a;
                double tmproot = B*B - 4.0*A*C;
                
                if(tmproot < 0) pznu = - B/(2*A); // take real part for complex roots
                else
                {
                  double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
                  double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);
                  pznu = TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ? tmpsol1 : tmpsol2;
                }
                TLorentzVector met(mets[0]->Px(), mets[0]->Py(), pznu, sqrt(mets[0]->Px()*mets[0]->Px() + mets[0]->Py()*mets[0]->Py() + pznu*pznu ) );
                
                histo1D["mTTbar"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]+matchedPFJetsGenJets[3]+lepton+met).M() );
                histo2D["mTTbar_VS_mTop"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]+matchedPFJetsGenJets[3]+lepton+met).M(), mTop);
                histo2D["mTTbar_VS_mTopGenJet"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]+matchedPFJetsGenJets[3]+lepton+met).M(), mTopGenJet);
                histo1D["PtTTbar"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]+matchedPFJetsGenJets[3]+lepton+met).Pt() );
                histo2D["PtTTbar_VS_mTop"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]+matchedPFJetsGenJets[3]+lepton+met).Pt(), mTop);
                histo2D["PtTTbar_VS_mTopGenJet"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]+matchedPFJetsGenJets[3]+lepton+met).Pt(), mTopGenJet);
                histo1D["nBtags"]->Fill(nBtags);
                histo2D["nBtags_VS_mTop"]->Fill(nBtags, mTop);
                histo2D["nBtags_VS_mTopGenJet"]->Fill(nBtags, mTopGenJet);
                histo1D["PtTop"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]).Pt() );
                histo2D["PtTop_VS_mTop"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]).Pt(), mTop);
                histo2D["PtTop_VS_mTopGenJet"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]).Pt(), mTopGenJet);
                histo1D["EtaTop"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]).Eta() );
                histo2D["EtaTop_VS_mTop"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]).Eta(), mTop);
                histo2D["EtaTop_VS_mTopGenJet"]->Fill( (matchedPFJetsGenJets[0]+matchedPFJetsGenJets[1]+matchedPFJetsGenJets[2]).Eta(), mTopGenJet);
                histo1D["PtBjet"]->Fill( matchedPFJetsGenJets[2].Pt() );
                histo2D["PtBjet_VS_mTop"]->Fill( matchedPFJetsGenJets[2].Pt(), mTop);
                histo2D["PtBjet_VS_mTopGenJet"]->Fill( matchedPFJetsGenJets[2].Pt(), mTopGenJet);
                histo1D["EtaBjet"]->Fill( matchedPFJetsGenJets[2].Eta() );
                histo2D["EtaBjet_VS_mTop"]->Fill( matchedPFJetsGenJets[2].Eta(), mTop);
                histo2D["EtaBjet_VS_mTopGenJet"]->Fill( matchedPFJetsGenJets[2].Eta(), mTopGenJet);
              }
            }
          }
        }
        
				///////////////
	      // KINFITTER //
	      ///////////////

   	    kinFit->SetJets(selectedJets);
   	    
   	    vector<TH2F> monsterVector;
   	    vector< float > mvaValsVector;
   	    vector< vector<unsigned int> > mvaResultsVector;
   	    vector< vector< float > > topMassVector; // mTopFit, sigmaMTopFit, chi2MTopFit
   	    for(unsigned int iCombi=0; iCombi<12; iCombi++)
   	    {
   	      pair<float, vector<unsigned int> > tmpMvaVals = jetCombiner->getMVAValue(MVAmethod, iCombi+1);
   	      mvaResultsVector.push_back(tmpMvaVals.second);
   	      mvaValsVector.push_back(tmpMvaVals.first);
          kinFit->SetMVAStuff(tmpMvaVals);
          
          bool correctCombi = false;
          if( ( (tmpMvaVals.second[0] == lightJet1Index && tmpMvaVals.second[1] == lightJet2Index) || (tmpMvaVals.second[1] == lightJet1Index && tmpMvaVals.second[0] == lightJet2Index) ) && tmpMvaVals.second[2] == hadrBJetIndex ) correctCombi = true;
          
          if(measureTopMassDifference)
          {
            vector<float> tmp;
            float* res = 0;
            if(systematic == "mWdiff" || systematic == "mWdiffSmall" || systematic == "mWdiffSmaller")
            {
              double dMw = 0.3;
              if(systematic == "mWdiffSmall") dMw = 0.2;
              else if(systematic == "mWdiffSmaller") dMw = 0.1;
              if(eventSelectedSemiMu)
              {
                if(selectedMuons[0]->charge() > 0) res = kinFit->EstimateTopMass(event, 80.4-dMw, false, iCombi, correctCombi);
                else res = kinFit->EstimateTopMass(event, 80.4+dMw, false, iCombi, correctCombi);
              }
              else
              {
                if(selectedElectrons[0]->charge() > 0) res = kinFit->EstimateTopMass(event, 80.4-dMw, false, iCombi, correctCombi);
                else res = kinFit->EstimateTopMass(event, 80.4+dMw, false, iCombi, correctCombi);
              }
            }
            else // Nominal !!
            {
              if(correctCombi) res = kinFit->EstimateTopMass(event, 80.4, true, iCombi, correctCombi);
              else res = kinFit->EstimateTopMass(event, 80.4, false, iCombi, correctCombi);
            }
            
            tmp.push_back(res[0]);
            tmp.push_back(res[1]);
            tmp.push_back(res[2]);
            topMassVector.push_back(tmp);
            delete res;
          }
          else
          {
        	  TH2F* histo = 0;
        	  if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
        	    histo = kinFit->FitEvent(event, 80.4, 173.3, false, iCombi); // As measured by the Tevatron //switch first boolean to true to save all monsters
            else
              histo = kinFit->FitEvent(event, 80.4, 172.5, false, iCombi); // As used in the MC
            monsterVector.push_back(*histo);
            delete histo;
          }
        }
        
        if(measureTopMassDifference)
        {
          // check top quark masses and which top is decaying hadronically
          float topMass = -1., antiTopMass = -1.;
          bool topDecayedLept = false;
          TLorentzVector top, antiTop;
          
          if(dataSetName.find("TTbarJets") == 0 || dataSetName.find("TT_") == 0 || dataSetName.find("ST_") == 0 || dataSetName.find("SingleTop") == 0)
          {
            for(unsigned int iPart=0; iPart<mcParticles.size(); iPart++)
            {
              if( mcParticles[iPart]->status() != 3) continue;
//                cout << "type: " << mcParticles[iPart]->type() << endl;
              if( mcParticles[iPart]->type() == 6 )
              {
                topMass = mcParticles[iPart]->M();
                top = *mcParticles[iPart];
              }
              else if( mcParticles[iPart]->type() == -6 )
              {
                antiTopMass = mcParticles[iPart]->M();
                antiTop = *mcParticles[iPart];
              }
              else if( mcParticles[iPart]->type() == -13 && mcParticles[iPart]->motherType() == 24 && mcParticles[iPart]->grannyType() == 6 )
                topDecayedLept = true;
              else if( mcParticles[iPart]->type() == -11 && mcParticles[iPart]->motherType() == 24 && mcParticles[iPart]->grannyType() == 6 )
                topDecayedLept = true;
              if( fabs( mcParticles[iPart]->type() ) == 24 )
                histo1D["mWGen"]->Fill( mcParticles[iPart]->M() );
            }
          }
          
          lightMonster = new LightMonster();
          lightMonster->setEventID( event->eventId() );
          lightMonster->setRunID( event->runId() );
          lightMonster->setLumiBlockID( event->lumiBlockId() );
          lightMonster->setIdParton1( event->idParton1() );
          lightMonster->setXParton1( event->xParton1() );
          lightMonster->setIdParton2( event->idParton2() );
          lightMonster->setXParton2( event->xParton2() );
          lightMonster->setFactorizationScale( event->factorizationScale() );
          lightMonster->setFlavHistPath( event->flavorHistoryPath() );
          lightMonster->setNPV(vertex.size());
          lightMonster->setNTruePU(event->nTruePU());
          lightMonster->setTopMass(topMass);
          lightMonster->setAntiTopMass(antiTopMass);
          lightMonster->setSelectedSemiMu(eventSelectedSemiMu);
          if(dataSetName.find("TTbarJets") == 0 || dataSetName.find("TT_") == 0)
          {
            lightMonster->setSemiMuDecay(genEvt->isSemiLeptonic( TRootGenEvent::kMuon ));
            lightMonster->setSemiElDecay(genEvt->isSemiLeptonic( TRootGenEvent::kElec ));
          }
          lightMonster->setTopDecayedLept(topDecayedLept);
          lightMonster->setAll4JetsMCMatched( jetCombiner->All4JetsMatched_MCdef() );
          lightMonster->setAllHadronicJetsMCMatched( jetCombiner->HadronicTopJetsMatched_MCdef() );
          lightMonster->setMvaVals(mvaValsVector);
          lightMonster->setMvaResults(mvaResultsVector);
          lightMonster->setEventWeight(scaleFactor);
          lightMonster->setMTopFitResults(topMassVector);
          lightMonster->setHadrBJet( hadrBJetIndex );
          lightMonster->setHadrLJet1( lightJet1Index );
          lightMonster->setHadrLJet2( lightJet2Index );
          lightMonster->setLeptBJet( leptBJetIndex );
          lightMonster->setMET( *mets[0] );
          lightMonster->setSelectedJets( otherSelectedJets );
          lightMonster->setBTagCSV(bTagCSV);
          lightMonster->setPartonFlavour(partonFlavour);
          vector<TLorentzVector> triggerLeptons;
          std::map<std::string, std::vector<TopTree::triggeredObject> > trigFilters = event->getTriggerFilters();
          if( eventSelectedSemiMu )
          {
            lightMonster->setLepton( *selectedMuons[0] );
            lightMonster->setLeptonCharge( selectedMuons[0]->charge() );
            float reliso = (selectedMuons[0]->chargedHadronIso() + max( 0.0, selectedMuons[0]->neutralHadronIso() + selectedMuons[0]->photonIso() - 0.5*selectedMuons[0]->puChargedHadronIso() ) ) / selectedMuons[0]->Pt(); // dBeta corrected
            lightMonster->setLeptonPFRelIso( reliso );
            if(trigFilters.find("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15") != trigFilters.end())
            {
              std::vector<TopTree::triggeredObject> trigObj = trigFilters["hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15"];
              for(int i=0; i<trigObj.size(); i++)
              {
                TLorentzVector tmp;
                tmp.SetPtEtaPhiM(trigObj[i].pt, trigObj[i].eta, trigObj[i].phi, 0.);
                triggerLeptons.push_back( tmp );
              }
            }
            else if(trigFilters.find("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10") != trigFilters.end())
            {
              std::vector<TopTree::triggeredObject> trigObj = trigFilters["hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10"];
              for(int i=0; i<trigObj.size(); i++)
              {
                TLorentzVector tmp;
                tmp.SetPtEtaPhiM(trigObj[i].pt, trigObj[i].eta, trigObj[i].phi, 0.);
                triggerLeptons.push_back( tmp );
              }
            }
            else cout << "No known triggerFilter found for run:  " << event->runId() << endl;
          }
          else
          {
            lightMonster->setLepton( *selectedElectrons[0] );
            lightMonster->setLeptonCharge( selectedElectrons[0]->charge() );
            double EffectiveArea = 0.;
            double SCeta = fabs(selectedElectrons[0]->superClusterEta());
            // HCP 2012 updated for electron conesize = 0.3, taken from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h?revision=1.4&view=markup
            if      (SCeta >= 0.0   && SCeta < 1.0   ) EffectiveArea = 0.130;
            else if (SCeta >= 1.0   && SCeta < 1.479 ) EffectiveArea = 0.137;
            else if (SCeta >= 1.479 && SCeta < 2.0   ) EffectiveArea = 0.067;
            else if (SCeta >= 2.0   && SCeta < 2.2   ) EffectiveArea = 0.089;
            else if (SCeta >= 2.2   && SCeta < 2.3   ) EffectiveArea = 0.107;
            else if (SCeta >= 2.3   && SCeta < 2.4   ) EffectiveArea = 0.110;
            else if (SCeta >= 2.4) EffectiveArea = 0.138;
            float reliso = (selectedElectrons[0]->chargedHadronIso() + max( selectedElectrons[0]->neutralHadronIso() + selectedElectrons[0]->photonIso()  - event->kt6PFJets_rho()*EffectiveArea, 0.) )/ selectedElectrons[0]->Pt(); // EffectiveArea corrected
            lightMonster->setLeptonPFRelIso( reliso );
            if(trigFilters.find("hltEle27WP80TrackIsoFilter") != trigFilters.end())
            {
              std::vector<TopTree::triggeredObject> trigObj = trigFilters["hltEle27WP80TrackIsoFilter"];
              for(int i=0; i<trigObj.size(); i++)
              {
                TLorentzVector tmp;
                tmp.SetPtEtaPhiM(trigObj[i].pt, trigObj[i].eta, trigObj[i].phi, 0.);
                triggerLeptons.push_back( tmp );
              }
            }
            else cout << "No known triggerFilter found for run:  " << event->runId() << endl;
          }
          lightMonster->setTriggerLeptons( triggerLeptons );
          lightMonster->setHadrBQuark( jetCombiner->GetHadrBQuark() );
          lightMonster->setHadrLQuark1( jetCombiner->GetLightQuark1() );
          lightMonster->setHadrLQuark2( jetCombiner->GetLightQuark2() );
          lightMonster->setLeptBQuark( jetCombiner->GetLeptBQuark() );
          lightMonster->setGenTop( top );
          lightMonster->setGenAntiTop( antiTop );
          
          MonsterTree->Fill();
          delete lightMonster;
	      }
      }// end !TrainMVA && eventSelected
    }				//loop on events
    cout<<endl;
    
    cout << "bTagEff = " << nBjetsBtag / nBjets << "   missTagRate = " << nNonBjetsBtag / nNonBjets << endl;

    if( !CalculateResolutions ) kinFit->Write(fout, true, pathPNG+"FullKinFit/");
    delete kinFit;
    
    MonsterFile->cd();
    
    TTree *configTreeMonsterFile = new TTree("configTreeMonsterFile","configuration Tree in Monster File");
    TClonesArray* tcdatasetmonsterfile = new TClonesArray("Dataset",1);
    configTreeMonsterFile->Branch("Dataset","TClonesArray",&tcdatasetmonsterfile);
    new ((*tcdatasetmonsterfile)[0]) Dataset(*datasets[d]);
    
    configTreeMonsterFile->Fill();
    configTreeMonsterFile->Write();
    MonsterTree->Write();
    MonsterFile->Close();
    delete MonsterFile;
    
    if(jetTools) delete jetTools;
    
    treeLoader.UnLoadDataset(); //important: free memory
  }				//loop on datasets
  
  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;
  
  //Selection tables
  selecTableSemiMu.TableCalculator(false, true, true, true, true);
  selecTableSemiEl.TableCalculator(false, true, true, true, true);
  string selectiontableSemiMu = "SelectionTable_SemiMu_JES.tex";
  string selectiontableSemiEl = "SelectionTable_SemiEl_JES.tex";
	selecTableSemiMu.Write(selectiontableSemiMu.c_str());
	selecTableSemiEl.Write(selectiontableSemiEl.c_str());

  // Do some special things with certain plots (normalize, BayesDivide, ... )
  if (verbose > 0)
    cout << "Treating the special plots." << endl;
  
  double respLightDR[nBinsJet], respLightDRErr[nBinsJet], respBDR[nBinsJet], respBDRErr[nBinsJet], respLightArea[nBinsJet], respLightAreaErr[nBinsJet], respBArea[nBinsJet], respBAreaErr[nBinsJet], respLightNpart[nBinsJet], respLightNpartErr[nBinsJet], respBNpart[nBinsJet], respBNpartErr[nBinsJet];
  double dR[nBinsJet], dRErr[nBinsJet], area[nBinsJet], areaErr[nBinsJet], nPart[nBinsJet], nPartErr[nBinsJet];
  
  for(size_t i=0; i<nBinsJet; i++)
  {
    dR[i] = (binningDR[i]+binningDR[i+1])/2.;
    dRErr[i] = fabs(binningDR[i+1]-binningDR[i])/2.;
    area[i] = (binningJetArea[i]+binningJetArea[i+1])/2.;
    areaErr[i] = fabs(binningJetArea[i+1]-binningJetArea[i])/2.;
    nPart[i] = (binningNpart[i]+binningNpart[i+1])/2.;
    nPartErr[i] = fabs(binningNpart[i+1]-binningNpart[i])/2.;
    stringstream ss; ss << i;
    vector<double> lightDR = resFitLightJets->ExtractSigmaMean(histo1D["lightJetResp_DR_"+ss.str()]);
    respLightDR[i] = lightDR[2];
    respLightDRErr[i] = lightDR[3];
    vector<double> bDR = resFitLightJets->ExtractSigmaMean(histo1D["bJetResp_DR_"+ss.str()]);
    respBDR[i] = bDR[2];
    respBDRErr[i] = bDR[3];
    vector<double> lightArea = resFitLightJets->ExtractSigmaMean(histo1D["lightJetResp_Area_"+ss.str()]);
    respLightArea[i] = lightArea[2];
    respLightAreaErr[i] = lightArea[3];
    vector<double> bArea = resFitLightJets->ExtractSigmaMean(histo1D["bJetResp_Area_"+ss.str()]);
    respBArea[i] = bArea[2];
    respBAreaErr[i] = bArea[3];
    vector<double> lightNpart = resFitLightJets->ExtractSigmaMean(histo1D["lightJetResp_nPart_"+ss.str()]);
    respLightNpart[i] = lightNpart[2];
    respLightNpartErr[i] = lightNpart[3];
    vector<double> bNpart = resFitLightJets->ExtractSigmaMean(histo1D["bJetResp_nPart_"+ss.str()]);
    respBNpart[i] = bNpart[2];
    respBNpartErr[i] = bNpart[3];
  }
  
  graphErr["lightResp_DR"] = new TGraphErrors(nBinsJet, dR, respLightDR, dRErr, respLightDRErr);
  graphErr["lightResp_DR"]->SetNameTitle("lightResp_DR","lightResp_DR");
  graphErr["bResp_DR"] = new TGraphErrors(nBinsJet, dR, respBDR, dRErr, respBDRErr);
  graphErr["bResp_DR"]->SetNameTitle("bResp_DR","bResp_DR");
  graphErr["lightResp_Area"] = new TGraphErrors(nBinsJet, area, respLightArea, areaErr, respLightAreaErr);
  graphErr["lightResp_Area"]->SetNameTitle("lightResp_Area","lightResp_Area");
  graphErr["bResp_Area"] = new TGraphErrors(nBinsJet, area, respBArea, areaErr, respBAreaErr);
  graphErr["bResp_Area"]->SetNameTitle("bResp_Area","bResp_Area");
  graphErr["lightResp_nPart"] = new TGraphErrors(nBinsJet, nPart, respLightNpart, nPartErr, respLightNpartErr);
  graphErr["lightResp_nPart"]->SetNameTitle("lightResp_nPart","lightResp_nPart");
  graphErr["bResp_nPart"] = new TGraphErrors(nBinsJet, nPart, respBNpart, nPartErr, respBNpartErr);
  graphErr["bResp_nPart"]->SetNameTitle("bResp_nPart","bResp_nPart");
  
  ///////////////////
  // Writing
  //////////////////
  if (verbose > 1)
  	cout << " - Writing outputs to the files ..." << endl;

	string pathPNGJetCombi = pathPNG+"JetCombination/";

  mkdir((pathPNG+"MSPlot/").c_str(),0777);
  mkdir(pathPNGJetCombi.c_str(),0777);
  
  jetCombiner->Write(fout, true, pathPNGJetCombi, false);
  
  // Fill the resolution histograms and calculate the resolutions
  if(CalculateResolutions)
  {
    cout << "Fit the resolution stuff..." << endl;
    mkdir((pathPNG+"resFit_LightJet/").c_str(),0777);
    mkdir((pathPNG+"resFit_LightJet_PtTopRew/").c_str(),0777);
    mkdir((pathPNG+"resFit_BJet/").c_str(),0777);
    mkdir((pathPNG+"resFit_BJet_PtTopRew/").c_str(),0777);
    mkdir((pathPNG+"resFit_BJet_B/").c_str(),0777);
    mkdir((pathPNG+"resFit_BJet_Bbar/").c_str(),0777);
    
    resFitLightJets->WritePlots(fout, true, pathPNG+"resFit_LightJet/");
    resFitLightJets->WriteResolutions("lightJetReso.root");
    resFitLightJets_PtTopRew->WritePlots(fout, true, pathPNG+"resFit_LightJet_PtTopRew/");
    resFitLightJets_PtTopRew->WriteResolutions("lightJetReso_PtTopRew.root");
    resFitBJets->WritePlots(fout, true, pathPNG+"resFit_BJet/");
    resFitBJets->WriteResolutions("bJetReso.root");
    resFitBJets_PtTopRew->WritePlots(fout, true, pathPNG+"resFit_BJet_PtTopRew/");
    resFitBJets_PtTopRew->WriteResolutions("bJetReso_PtTopRew.root");
    resFitBJets_B->WritePlots(fout, true, pathPNG+"resFit_BJet_B/");
    resFitBJets_B->WriteResolutions("bJetReso_B.root");
    resFitBJets_Bbar->WritePlots(fout, true, pathPNG+"resFit_BJet_Bbar/");
    resFitBJets_Bbar->WriteResolutions("bJetReso_Bbar.root");
  }
  
  cout << "Writing the histograms..." << endl;
  // 1D 
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(false, name, true, true, true, true, true);
    temp->Write(fout, name, true, pathPNG+"MSPlot/");
  }
  
//void Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false, bool mergeVVV = false, bool mergeSameSignWW = false);
  
  //Write histograms
  fout->cd();
  th1dir->cd();
	for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
	{
		TH1F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}

  // 2D
  TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
  th2dir->cd();
	for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
	{
		TH2F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}

	//Write TGraphAsymmErrors
	for(map<string,TGraphAsymmErrors*>::const_iterator it = graphAsymmErr.begin(); it != graphAsymmErr.end(); it++)
	{
	  TGraphAsymmErrors *temp = it->second;
	  temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}
	
	fout->cd();
	//Write TGraph2D
	for(map<string,TGraph2D*>::const_iterator it = graph2D.begin(); it != graph2D.end(); it++)
	{
	  TGraph2D *temp = it->second;
	  temp->SetNameTitle( (it->first).c_str(), (it->first).c_str() );
	  temp->Write();
	}
	
	fout->cd();
	//Write TProfile
	for(map<string,TProfile*>::const_iterator it = profile.begin(); it != profile.end(); it++)
	{
	  TProfile *temp = it->second;
	  temp->Write();
	  TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}
	
  fout->cd();
  //add configuration info
  fout->cd();
  configTree->Fill();
  configTree->Write();

  //Write TGraphErrors
  fout->cd();
  for(map<string,TGraphErrors*>::const_iterator it = graphErr.begin(); it != graphErr.end(); it++)
  {
    TGraphErrors *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  
  //Write TGraph
  fout->cd();
  for(map<string,TGraph*>::const_iterator it = graph.begin(); it != graph.end(); it++)
  {
    TGraph *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }

  //
  if (verbose > 1)
    cout << " - Done with writing the module outputs in the ouput file ..." << endl;
  cout << " - Closing the output file now..." << endl;
//  fout->Write();
  fout->Close();

  //delete
//  if(resFitLightJets) delete resFitLightJets;
//  if(resFitBJets) delete resFitBJets;
//  if(resFitBJets_B) delete resFitBJets_B;
//  if(resFitBJets_Bbar) delete resFitBJets_Bbar;
  if(jetCombiner) delete jetCombiner;
  
  delete fout;
  delete tcdatasets;
  delete configTree;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;

  return 0;
}
