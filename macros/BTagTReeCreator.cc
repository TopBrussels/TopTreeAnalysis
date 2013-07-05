////////////////////////////
///// TODO & COMMENTS /////
/////////////////////////// 

#include "TStyle.h"
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
#include "TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../BtagEffAnalysis/interface/TRootNTuple.h"
#include "Style.C"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/MCInformation/interface/Lumi3DReWeighting.h"

using namespace std;
using namespace reweight;
using namespace TopTree;

int main (int argc, char *argv[])
{

  /////////////////////////////////////
  // SWITCHES FOR SYSTAMATIC SAMPLES //
  /////////////////////////////////////
  
  int systematic = 0; // 0: off 1: minus 2: plus
  
  if (argc > 1) // JES shift is command line parameter 1
    if (atoi(argv[1]) == 0 || atoi(argv[1]) == 1 || atoi(argv[1]) == 2 || atoi(argv[1]) == 3 || atoi(argv[1]) == 4  || atoi(argv[1]) == 5 || atoi(argv[1]) == 6  || atoi(argv[1]) == 7 || atoi(argv[1]) == 8 || atoi(argv[1]) == 9 || atoi(argv[1]) == 10 || atoi(argv[1]) == 11 || atoi(argv[1]) == 12 || atoi(argv[1]) == 13 || atoi(argv[1]) == 14 || atoi(argv[1]) == 15 || atoi(argv[1]) == 16 || atoi(argv[1]) == 17 || atoi(argv[1]) == 18 ||atoi(argv[1]) == -1) // -1 == do the invertediso setup
      systematic=atoi(argv[1]);
  
  bool runSpecificSample=false;
  
  int nSpecSample=-19;

  if (argc > 2) {

    runSpecificSample=true;

    nSpecSample=atoi(argv[2]);

  }

  bool correctJets=false;
  double jetScale_mu = 1;
  double jetScale_el = 1;

  std::map<std::string,double> jetScale_const;

  // MadGraph
  /*jetScale_const["mu"]=1.00858;
  jetScale_const["mu_JESMinus"]=1.02365;
  jetScale_const["mu_JESPlus"]=0.99427;
  jetScale_const["mu_JERMinus"]=1.01093;
  jetScale_const["mu_JERPlus"]=1.00418;
  jetScale_const["mu_LessPU"]=1.00986;
  jetScale_const["mu_MorePU"]=1.00759;
  jetScale_const["mu_METUnClusteredEnergyMinus"]=1.00858;
  jetScale_const["mu_METUnClusteredEnergyPlus"]=1.00858;
  jetScale_const["mu_ttjetsscaleup"]=1.00869;
  jetScale_const["mu_ttjetsscaledown"]=1.00421;
  jetScale_const["mu_ttjetsmatchingup"]=1.00958;
  jetScale_const["mu_ttjetsmatchingdown"]=1.00285;
  jetScale_const["mu_ttjetsmassup"]=1.00182;
  jetScale_const["mu_ttjetsmassdown"]=1.00358;
  jetScale_const["el"]=1.01116;
  jetScale_const["el_JESMinus"]=1.02926;
  jetScale_const["el_JESPlus"]=0.99878;
  jetScale_const["el_JERMinus"]=1.01688;
  jetScale_const["el_JERPlus"]=1.01361;
  jetScale_const["el_LessPU"]=1.01663;
  jetScale_const["el_MorePU"]=1.01523;
  jetScale_const["el_METUnClusteredEnergyMinus"]=1.01116;
  jetScale_const["el_METUnClusteredEnergyPlus"]=1.01116;
  jetScale_const["el_ttjetsscaleup"]=1.01791;
  jetScale_const["el_ttjetsscaledown"]=1.02034;
  jetScale_const["el_ttjetsmatchingup"]=1.01806;
  jetScale_const["el_ttjetsmatchingdown"]=1.00836;
  jetScale_const["el_ttjetsmassup"]=1.02497;
  jetScale_const["el_ttjetsmassdown"]=1.02085;*/

  // MC@NLO
  /* jetScale_const["mu"]=1.01102;
  jetScale_const["mu_JESMinus"]=1.02892;
  jetScale_const["mu_JESPlus"]=0.992269;
  jetScale_const["mu_JERMinus"]=1.01835;
  jetScale_const["mu_JERPlus"]=1.01137;
  jetScale_const["mu_LessPU"]=1.01587;
  jetScale_const["mu_MorePU"]=1.01395;
  jetScale_const["mu_METUnClusteredEnergyMinus"]=1.01102;
  jetScale_const["mu_METUnClusteredEnergyPlus"]=1.01102;
  jetScale_const["mu_ttjetsscaleup"]=1.00869;
  jetScale_const["mu_ttjetsscaledown"]=1.00421;
  jetScale_const["mu_ttjetsmatchingup"]=1.00958;
  jetScale_const["mu_ttjetsmatchingdown"]=1.00285;
  jetScale_const["mu_ttjetsmassup"]=1.00182;
  jetScale_const["mu_ttjetsmassdown"]=1.00358;
  jetScale_const["el"]=1.01759;
  jetScale_const["el_JESMinus"]=1.03879;
  jetScale_const["el_JESPlus"]=1.00362;
  jetScale_const["el_JERMinus"]=1.02954;
  jetScale_const["el_JERPlus"]=1.02301;
  jetScale_const["el_LessPU"]=1.02611;
  jetScale_const["el_MorePU"]=1.02517;
  jetScale_const["el_METUnClusteredEnergyMinus"]=1.01759;
  jetScale_const["el_METUnClusteredEnergyPlus"]=1.01759;
  jetScale_const["el_ttjetsscaleup"]=1.01791;
  jetScale_const["el_ttjetsscaledown"]=1.02034;
  jetScale_const["el_ttjetsmatchingup"]=1.01806;
  jetScale_const["el_ttjetsmatchingdown"]=1.00836;
  jetScale_const["el_ttjetsmassup"]=1.02508;
  jetScale_const["el_ttjetsmassdown"]=1.02085;

  jetScale_const["mu_madgraph"]=1.00857;
  jetScale_const["el_madgraph"]=1.01116; */
  
  /*  jetScale_const["mu"]=1.01102;
  jetScale_const["mu_JESMinus"]=1.02077;
  jetScale_const["mu_JESPlus"]=1.0007;
  jetScale_const["mu_JERMinus"]=1.02325;
  jetScale_const["mu_JERPlus"]=0.997945;
  jetScale_const["mu_LessPU"]=1.01587;
  jetScale_const["mu_MorePU"]=1.01395;
  jetScale_const["mu_METUnClusteredEnergyMinus"]=1.01102;
  jetScale_const["mu_METUnClusteredEnergyPlus"]=1.01102;
  jetScale_const["mu_ttjetsscaleup"]=1.00869;
  jetScale_const["mu_ttjetsscaledown"]=1.00421;
  jetScale_const["mu_ttjetsmatchingup"]=1.00958;
  jetScale_const["mu_ttjetsmatchingdown"]=1.00285;
  jetScale_const["mu_ttjetsmassup"]=1.00182;
  jetScale_const["mu_ttjetsmassdown"]=1.00358;
  jetScale_const["el"]=1.0251;
  jetScale_const["el_JESMinus"]=1.04199;
  jetScale_const["el_JESPlus"]=1.01186;
  jetScale_const["el_JERMinus"]=1.02954;
  jetScale_const["el_JERPlus"]=1.02241;
  jetScale_const["el_LessPU"]=1.02611;
  jetScale_const["el_MorePU"]=1.0232;
  jetScale_const["el_METUnClusteredEnergyMinus"]=1.0251;
  jetScale_const["el_METUnClusteredEnergyPlus"]=1.0251;
  jetScale_const["el_ttjetsscaleup"]=1.01791;
  jetScale_const["el_ttjetsscaledown"]=1.02034;
  jetScale_const["el_ttjetsmatchingup"]=1.01806;
  jetScale_const["el_ttjetsmatchingdown"]=1.00836;
  jetScale_const["el_ttjetsmassup"]=1.02497;
  jetScale_const["el_ttjetsmassdown"]=1.02085;
  */

  // reweighed madgraph

  /*jetScale_const["mu"]=1.01149;
jetScale_const["mu_JESMinus"]=1.02657;
jetScale_const["mu_JESPlus"]=0.997401;
jetScale_const["mu_JERMinus"]=1.00384;
jetScale_const["mu_JERPlus"]=1.00663;
jetScale_const["mu_LessPU"]=1.00104;
jetScale_const["mu_MorePU"]=1.01083;
jetScale_const["mu_METUnClusteredEnergyMinus"]=1.01149;
jetScale_const["mu_METUnClusteredEnergyPlus"]=1.01149;
jetScale_const["mu_madgraph"]=1.01199;
jetScale_const["mu_ttjetsscaleup"]=1.01253;
jetScale_const["mu_ttjetsscaledown"]=1.00429;
jetScale_const["mu_ttjetsmatchingup"]=1.01332;
jetScale_const["mu_ttjetsmatchingdown"]=1.01162;
jetScale_const["mu_ttjetsmassup"]=1.00516;
jetScale_const["mu_ttjetsmassdown"]=0.990899;
jetScale_const["el"]=1.01108;
jetScale_const["el_JESMinus"]=1.02488;
jetScale_const["el_JESPlus"]=0.991653;
jetScale_const["el_JERMinus"]=1.01761;
jetScale_const["el_JERPlus"]=1.00267;
jetScale_const["el_LessPU"]=1.01205;
jetScale_const["el_MorePU"]=1.01005;
jetScale_const["el_METUnClusteredEnergyMinus"]=1.01108;
jetScale_const["el_METUnClusteredEnergyPlus"]=1.01108;
jetScale_const["el_madgraph"]=1.01585;
jetScale_const["el_ttjetsscaleup"]=1.01368;
jetScale_const["el_ttjetsscaledown"]=1.00298;
jetScale_const["el_ttjetsmatchingup"]=1.00357;
jetScale_const["el_ttjetsmatchingdown"]=1.01378;
jetScale_const["el_ttjetsmassup"]=1.02164;
jetScale_const["el_ttjetsmassdown"]=1.00997;

jetScale_const["mu_mcnlo"]=1.00911;
jetScale_const["el_mcnlo"]=1.02049;*/

jetScale_const["mu"]=1.01473;
jetScale_const["mu_JESMinus"]=1.03564;
jetScale_const["mu_JESPlus"]=0.996635;
jetScale_const["mu_JERMinus"]=1.02016;
jetScale_const["mu_JERPlus"]=1.01176;
jetScale_const["mu_LessPU"]=1.0153;
jetScale_const["mu_MorePU"]=1.01368;
jetScale_const["mu_METUnClusteredEnergyMinus"]=1.01473;
jetScale_const["mu_METUnClusteredEnergyPlus"]=1.01473;
jetScale_const["mu_ttjetsscaleup"]=1.01653;
jetScale_const["mu_ttjetsscaledown"]=0.999601;
jetScale_const["mu_ttjetsmatchingup"]=1.01584;
jetScale_const["mu_ttjetsmatchingdown"]=1.00853;
jetScale_const["mu_ttjetsmassup"]=1.00443;
jetScale_const["mu_ttjetsmassdown"]=1.00431;
jetScale_const["el"]=1.01235;
jetScale_const["el_JESMinus"]=1.03126;
jetScale_const["el_JESPlus"]=0.992497;
jetScale_const["el_JERMinus"]=1.01841;
jetScale_const["el_JERPlus"]=1.00999;
jetScale_const["el_LessPU"]=1.01327;
jetScale_const["el_MorePU"]=1.0113;
jetScale_const["el_METUnClusteredEnergyMinus"]=1.01686;
jetScale_const["el_METUnClusteredEnergyPlus"]=1.01686;
jetScale_const["el_ttjetsscaleup"]=1.01214;
jetScale_const["el_ttjetsscaledown"]=1.01635;
jetScale_const["el_ttjetsmatchingup"]=1.01188;
jetScale_const["el_ttjetsmatchingdown"]=1.01203;
jetScale_const["el_ttjetsmassup"]=1.03155;
jetScale_const["el_ttjetsmassdown"]=1.00407;

//jetScale_const["mu_mcnlo"]=1.01049;
//jetScale_const["el_mcnlo"]=1.02111;

jetScale_const["mu_mcnlo"]=1.02181;
jetScale_const["el_mcnlo"]=1.01462;

  if (argc > 1) {
    if (atoi(argv[1]) >= 40) {
      systematic=atoi(argv[1])-40;
      correctJets=true;
    }
  }

  //jetScale_el=jetScale_mu; // temp test

  //cout << systematic << " " << correctJets << endl;

  //exit(1);

  string postfix = ""; // to relabel the names of systematic samples

  if (systematic == 1)
    postfix="_JESMinus";
  if (systematic == 2)
    postfix="_JESPlus";
  if (systematic == 3)
    postfix="_METUnClusteredEnergyMinus";
  if (systematic == 4)
    postfix="_METUnClusteredEnergyPlus";
  if (systematic == 5)
    postfix="_ZeroPUevents";
  if (systematic == 6)
    postfix="_JERMinus";
  if (systematic == 7)
    postfix="_JERPlus";
  if (systematic == 10)
    postfix="_MorePU";
  if (systematic == 11)
    postfix="_LessPU";
  if (systematic == 12)
    postfix="_NomPU";

  if (systematic == 13)
    postfix="_ttjetsscaleup";
  if (systematic == 14)
    postfix="_ttjetsscaledown";
  if (systematic == 15)
    postfix="_ttjetsmatchingup";
  if (systematic == 16)
    postfix="_ttjetsmatchingdown";
  if (systematic == 17)
    postfix="_ttjetsmassup";
  if (systematic == 18)
    postfix="_ttjetsmassdown";

  if (correctJets) {

    jetScale_mu = jetScale_const["mu"+postfix];
    jetScale_el = jetScale_const["el"+postfix];
    
  } else {

    jetScale_mu=1;
    jetScale_el=1;

  }
   
  if (argc > 1) {
    if (atoi(argv[1]) == 60) {
      systematic=0;
      correctJets=true;
      jetScale_mu = jetScale_const["mu_madgraph"];
      jetScale_el = jetScale_const["el_madgraph"];
    }
  }

  cout << "systematic: " << systematic << " -> Postfix = " << postfix << " -> JetScale mu/el = " << jetScale_mu << "/" << jetScale_el << endl;

  // return 0;

  //exit(1);

  bool doPF2PAT = false;

  bool useMassesAndResolutions = false;

  //useMassesAndResolutions = true;

  if (systematic != 0 || runSpecificSample)
    useMassesAndResolutions = true;

  double btagcut=3;
  //btagcut=9999.;

  clock_t start = clock();

  cout << "********************************************************" << endl;
  cout << " Beginning of the program for creating the BTag Trees ! " << endl;
  cout << "********************************************************" << endl;

  //SetsStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

  //xml file
  string xmlFileName ="../config/myBTAGconfig.xml";
  //xmlFileName ="../config/myBTAGconfig_newcalib.xml";
  //string xmlFileName ="../config/myBTAGconfig_fall11.xml";
  //xmlFileName="../config/testje.xml";

  if (systematic == 8)
    xmlFileName ="../config/myBTAGconfig_ttjetsyst.xml";
    //xmlFileName ="../config/myBTAGconfig_ttjetsyst_newcalib.xml";
  else if (systematic == 9)
    xmlFileName ="../config/myBTAGconfig_wjetsyst.xml";

  if (argc > 3)
    xmlFileName = (string)argv[3];
  
  const char *xmlfile = xmlFileName.c_str();

  cout << "used config file: " << xmlfile << endl;

  //Output ROOT file
  string rootFileName ("BtaggingOutput"+postfix+".root");

  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);

  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<"Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
  
  cout << "analysis environment luminosity for rescaling "<< oldLuminosity << endl;
  
  /////////////////////
  // Load Datasets
  /////////////////////

	double yieldQCD_MC_SemiEl = 0;

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  vector < Dataset* > datasetsMu;
  vector < Dataset* > datasetsEl;
  vector < Dataset* > datasetsPlot;

  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;

  float LuminosityMu = oldLuminosity;
  float LuminosityEl = oldLuminosity;

  bool isSemiMu = false;
  bool isSemiE = false;

  bool foundMu = false;
  bool foundEl = false;

  for (unsigned int d = 0; d < datasets.size (); d++) {
    //cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
    if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
    string dataSetName = datasets[d]->Name();
    
    if(dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
      LuminosityMu = datasets[d]->EquivalentLumi();
      foundMu=true;
    }  
    if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) {
      LuminosityEl = datasets[d]->EquivalentLumi();
      foundEl=true;  
    }  

    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // when only using one channel
      Luminosity = datasets[d]->EquivalentLumi();
    }   

    if( dataSetName.find("QCD") == 0 ) datasets[d]->SetColor(kYellow);
    if( dataSetName.find("TT") == 0 ) datasets[d]->SetColor(kRed+1);
    if( dataSetName.find("TTbarJets_Other") == 0 ) datasets[d]->SetColor(kRed-7);
    if( dataSetName.find("WJets") == 0 )
    {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      datasets[d]->SetColor(kGreen-3);
    }
    if( dataSetName.find("ZJets") == 0 )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      datasets[d]->SetColor(kAzure-2);
    }
    if( dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") ==0 )
      datasets[d]->SetColor(kMagenta);
    
  }

  if(!foundMu && !foundEl && Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  else {
    if(LuminosityMu != oldLuminosity) cout << "Muon PD: changed analysis environment luminosity to "<< LuminosityMu << endl;
    if(LuminosityEl != oldLuminosity) cout << "Electron PD: changed analysis environment luminosity to "<< LuminosityEl << endl;
  }

  // make a datasets vector only for 
  if (foundMu) {
      for (unsigned int d = 0; d < datasets.size (); d++) {
	string dataSetName = datasets[d]->Name();
	if ( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) ) 
	  datasetsMu.push_back(datasets[d]);
      }
  }
      
  if (foundEl) {
    for (unsigned int d = 0; d < datasets.size (); d++) {
      string dataSetName = datasets[d]->Name();
      if ( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) ) 
	datasetsEl.push_back(datasets[d]);
    }
  }

  //vector of objects
  //cout << " - Variable declaration ..." << endl;
  /* vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootMET* > mets;
  vector < TRootGenJet* > genjet;*/
  
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //Global variable
  //TRootEvent* event = 0;
  
  //nof selected events
  double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
  double nTTbar = 0;
  double nTTbar_w = 0;

  //for chi2 performance
  int nAll=0;
  int nGoodCombi=0;
  int nGoodChi2=0;

  ////////////////////////////////////
  /// Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////
  
  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;
  
  histo1D["hadronicPartonTopMass"] = new TH1F("hadronicPartonTopMass","Hadronic Top Mass, using the Partons", 200, 150, 200);
  histo1D["hadronicPartonWMass"] = new TH1F("hadronicPartonWMass","Hadronic W Mass, using the Partons",100,0,200);
  histo1D["hadronicRecoTopMass"] = new TH1F("hadronicRecoTopMass","Hadronic Top Mass, using the RecoJets", 100, 100, 300);
  histo1D["hadronicRecoWMass"] = new TH1F("hadronicRecoWMass","Hadronic W Mass, using the RecoJets",100,0,200);

  histo1D["hadronicRecoTopMass-REW"] = new TH1F("hadronicRecoTopMass-REW","Hadronic Top Mass, using the RecoJets", 100, 100, 300);
  histo1D["hadronicRecoWMass-REW"] = new TH1F("hadronicRecoWMass-REW","Hadronic W Mass, using the RecoJets",100,0,200);
  histo1D["hadronicRecoTopMass-REW3D"] = new TH1F("hadronicRecoTopMass-REW3D","Hadronic Top Mass, using the RecoJets", 100, 100, 300);
  histo1D["hadronicRecoWMass-REW3D"] = new TH1F("hadronicRecoWMass-REW3D","Hadronic W Mass, using the RecoJets",100,0,200);
  
  ////////////////////////////////////
  /// MultiSamplePlot
  ////////////////////////////////////

  map<string,MultiSamplePlot*> MSPlot;

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
  CutsSelecTableSemiEl.push_back(string("Veto 2nd electron from Z-decay"));
  CutsSelecTableSemiEl.push_back(string("Conversion veto"));
  
  char LabelNJets[100];
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));

  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  selecTableSemiMu.SetLuminosity(LuminosityMu);
  SelectionTable selecTableSemiEl(CutsSelecTableSemiEl, datasets);
  selecTableSemiEl.SetLuminosity(LuminosityEl);

  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////
  // PileUp Reweighting //
  ////////////////////////

  //cout << Luminosity << endl;

  LumiReWeighting LumiWeights_mu, LumiWeightsUp_mu, LumiWeightsDown_mu;
  LumiReWeighting LumiWeights_el, LumiWeightsUp_el, LumiWeightsDown_el;
  
  LumiWeights_mu = LumiReWeighting("PileUpReweighting/pileup_MC_Summer12_S10.root", "PileUpReweighting/pileup_2012Data53X_UpToRun208686_Mu/nominal.root", "pileup", "pileup");
  LumiWeightsUp_mu = LumiReWeighting("PileUpReweighting/pileup_MC_Summer12_S10.root", "PileUpReweighting/pileup_2012Data53X_UpToRun208686_Mu/sys_up.root", "pileup", "pileup");
  LumiWeightsDown_mu = LumiReWeighting("PileUpReweighting/pileup_MC_Summer12_S10.root", "PileUpReweighting/pileup_2012Data53X_UpToRun208686_Mu/sys_down.root", "pileup", "pileup");

  LumiWeights_el = LumiReWeighting("PileUpReweighting/pileup_MC_Summer12_S10.root", "PileUpReweighting/pileup_2012Data53X_UpToRun208686_El/nominal.root", "pileup", "pileup");
  LumiWeightsUp_el = LumiReWeighting("PileUpReweighting/pileup_MC_Summer12_S10.root", "PileUpReweighting/pileup_2012Data53X_UpToRun208686_El/sys_up.root", "pileup", "pileup");
  LumiWeightsDown_el = LumiReWeighting("PileUpReweighting/pileup_MC_Summer12_S10.root", "PileUpReweighting/pileup_2012Data53X_UpToRun208686_El/sys_down.root", "pileup", "pileup");
  

  cout << " Initialized LumiReWeighting stuff" << endl;
  
  //exit(1);

  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////

  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++) {
    //   for (unsigned int d = 2; d < 3; d++) {

    // for the data driven QCD estimate, we count the number of data and inverted iso data in 3 bins of the runrange to make shure all runranges are equaly represented. We rescale the bins to data and then to the yield obtained with MC

    vector<int> QCDEstcounts; QCDEstcounts.push_back(0); QCDEstcounts.push_back(0); QCDEstcounts.push_back(0);

    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();

    if (runSpecificSample && d != nSpecSample) continue;

    //if (!(dataSetName.find("QCD_") == 0)) continue;

    if (systematic == -1 && !(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) ) continue;

    // temp fix for pfmet name in stijn's ElectronHad sample

    //if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0 )
    //anaEnv.METCollection = "PFMET_patMETsPF2PAT";
    
    if(dataSetName.find("InvIso") != string::npos) {
	anaEnv.METCollection = anaEnv.METCollection+"NoLeptonCleaning";
	anaEnv.MuonCollection = anaEnv.MuonCollection+"NoLeptonCleaning";
	anaEnv.ElectronCollection = anaEnv.ElectronCollection+"NoLeptonCleaning";
	anaEnv.JetCollection = anaEnv.JetCollection+"NoLeptonCleaning";
    }

    //cout << dataSetName << endl;
    //cout << anaEnv.METCollection << endl;
    //cout << anaEnv.ElectronCollection << endl;

    //exit(1);

    int nSelectedMu=0;
    int nSelectedEl=0;
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    //open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;
    
    
    
    /////////////////////////////////////
    /// Initialize JEC factors
    /////////////////////////////////////
   	    
    /*vector<JetCorrectorParameters> vCorrParam;

    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt", "Total")));

    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);*/

    // new JEC

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
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START53_V23_Summer13_Uncertainty_AK5PFchs.txt");
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
    
    ////////////////////////////////////////////////////////////
    // CREATE OUTPUT FILE AND TTREE FOR STORAGE OF THE NTUPLE //
    ////////////////////////////////////////////////////////////

    string dir = "BtagTrees/";

    mkdir(dir.c_str(),0777);

    //if (argc > 2)
    //  dir += string(argv[2])+"/";
    
    mkdir(dir.c_str(),0777);


    string bTitle = dir+"BtagTree_"+dataSetName+postfix+".root";

    cout << "INFO: creating BtagTree "+bTitle << endl;
      
    TFile* BTreeFile = new TFile(bTitle.c_str(),"RECREATE");
      
    TTree* BtagTTree= new TTree("tree","Tree containing my NTuple");
    // the ntuple container
    TRootNTuple* NTuple = new TRootNTuple();

    BtagTTree->Branch("TheNTuple","TRootNTuple",&NTuple);

    // create a histos map to store histos in the btagfile

    map<string,TH1F*> histo1D_btag;
    map<string,TH1D*> histo1D_syst;
    map<string,TH2F*> histo2D_btag;
    map<string,TGraphErrors*> graphErr_btag;
  
    // JET CALIBRATION

    if (correctJets && systematic == 0 & dataSetName.find("TTbarJets_mcnlo") != string::npos) {
	jetScale_mu = jetScale_const["mu_mcnlo"];
	jetScale_el = jetScale_const["el_mcnlo"];
    
	cout << "CHANGING JETSCALE FOR MC@NLO ==> systematic: " << systematic << " -> JetScale mu/el = " << jetScale_mu << "/" << jetScale_el << endl;
    }

    ////////////////////////////////////////////////////////
    // GET MASSES AND RESOLUTIONS FOR CHI2 JETCOMBINATION //
    ////////////////////////////////////////////////////////

    // load the mass values and resolutions for chi2 jetcombination

    float Chi2Wmass = -9999;
    float SigmaChi2Wmass = -9999;
    float Chi2Topmass = -9999;
    float SigmaChi2Topmass = -9999;
    
    if (useMassesAndResolutions) {

      //cout << "Test" << endl;
      
      string filename = "BtagMassPlots";

      //if (argc > 2)
      //filename += "_"+string(argv[2]);
      
      filename += ".root";
      
      cout << "INFO: Using masses and widths for the Chi2 jetcombiner from: " << filename << endl;
	    
      TFile* res = new TFile(filename.c_str(),"READ");
      
      res->cd();
      
      TF1* WmassFit = (TF1*)res->Get("hadronicRecoWMass_Fitted");
      //TF1* WmassFit = (TF1*)res->Get("hadronicRecoWMass-REW_Fitted");
      
      if (WmassFit) {
	
	//cout << "Fitted Wmass: " << WmassFit->GetParameter(1) << "+-" << WmassFit->GetParameter(2) << endl;
	Chi2Wmass = WmassFit->GetParameter(1); 
	SigmaChi2Wmass = WmassFit->GetParameter(2); 	

      }

      TF1* TopmassFit = (TF1*)res->Get("hadronicRecoTopMass_Fitted");
      //TF1* TopmassFit = (TF1*)res->Get("hadronicRecoTopMass-REW_Fitted");
      
      if (TopmassFit) {
	
	//cout << "Fitted Topmass: " << TopmassFit->GetParameter(1) << "+-" << TopmassFit->GetParameter(2) << endl;
	Chi2Topmass = TopmassFit->GetParameter(1); 
	SigmaChi2Topmass = TopmassFit->GetParameter(2);
      
      }
      res->Close();
      
      delete res;

    }

    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    // btag eff check

    int nBjets=0;
    int nJets=0;
    int nTagged=0;
    int nTagged_trueb=0;
    int nTagged_nonb=0;

    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
      //      for (unsigned int ievt = 0; ievt < 10000; ievt++)
    {
      
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;

      //cout << mets[400]->Pt();
      
      nEvents[d]++;
            
      if(ievt%1000 == 0)
	std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected: " << nSelectedMu << " (mu+jets) " << nSelectedEl << " (e+jets)" << flush<<"\r";

      ////////////////
      // LOAD EVENT //
      ////////////////

      //TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);  

      // do not use events with true pileup < 6 (not enough MC there)

      /*if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
	if ((int)event->nTruePU() < 6) {

	  //cout << "Skipping event " << ievt << " with # true PU = " << (int)event->nTruePU()  << endl;

	  continue;

	}
	}*/

      if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
        genjets = treeLoader.LoadGenJet(ievt,false);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }

      // check with genEvent which ttbar channel it is
      if(dataSetName.find("TTbarJets") == 0)  {
	//cout << "LOADING GenEvent" << endl;
	TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
	if( genEvt->isSemiLeptonic(TRootGenEvent::kMuon) ) {
	  isSemiMu=true;
	  isSemiE=false;
	}
	else if( genEvt->isSemiLeptonic(TRootGenEvent::kElec) ) {
	  isSemiMu=false;
	  isSemiE=true;
	}
	else {
	  isSemiMu=false;
	  isSemiE=false;
	}
      }


      /////////////////////////////////
      // DETERMINE EVENT SCALEFACTOR //
      /////////////////////////////////
      
      // scale factor for the event
      float scaleFactor = 1.;
      float scaleFactor_mu = 1;
      float scaleFactor_el = 1;
      
      // Load the GenEvent and calculate the branching ratio correction
      /*if(dataSetName.find("TTbarJets") == 0)
	{
	  //cout << "LOADING GenEvent" << endl;
	  TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
	  if( genEvt->isSemiLeptonic() )
	    scaleFactor *= (0.108*9.)*(0.676*1.5);
	  else if( genEvt->isFullHadronic() )
	    scaleFactor *= (0.676*1.5)*(0.676*1.5);
	  else if( genEvt->isFullLeptonic() )
	    scaleFactor *= (0.108*9.)*(0.108*9.);

      }*/

      // JES CORRECTION

      // before correction
      //for (unsigned int j=0;j<init_jets.size();j++)
      //	cout << "jet " << j << " pt " << init_jets[j]->Pt() << " eta " << init_jets[j]->Eta() << endl;

      //for (unsigned int j=0;j<init_jets_corrected.size();j++)
      //cout << "jet " << j << " pt " << init_jets_corrected[j]->Pt() << " eta " << init_jets_corrected[j]->Eta() << endl;
      
      // Clone the init_jets vector, otherwise the corrections will be removed
      //for(unsigned int i=0; i<init_jets_corrected.size(); i++)
      //if(init_jets_corrected[i]) delete init_jets_corrected[i];
      //init_jets_corrected.clear();
      
      //for(unsigned int i=0; i<init_jets.size(); i++)
      //  init_jets_corrected.push_back( (TRootJet*) init_jets[i]->Clone() );
      
      //////////////////////////////////////
      // Apply Jet Corrections on-the-fly //   
      //////////////////////////////////////

      // not needed for now, GT contains good stuff
      /*if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) {
      	//jetTools->correctJets(init_jets_corrected,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
      } else {
	jetTools->correctJets(init_jets_corrected,event->kt6PFJets_rho(),false); //last boolean: isData (needed for L2L3Residual...)
	}*/

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
      
      
      // after correction
      //for (unsigned int j=0;j<init_jets_corrected.size();j++)
      //cout << "after new JES jet " << j << " pt " << init_jets_corrected[j]->Pt() << " eta " << init_jets_corrected[j]->Eta() << endl;
      //cout << "****** end event *******"<< endl;
      //cout << endl;
      //exit(1);

      // PU reweighting

      // old method

      double lumiWeight_mu = LumiWeights_mu.ITweight( (int)event->nTruePU() );
      double lumiWeight_el = LumiWeights_el.ITweight( (int)event->nTruePU() );

      //cout << lumiWeight_mu << " " << lumiWeight_el << endl;

      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) {
	lumiWeight_mu=1;
	lumiWeight_el=1;
      }
      if (systematic == 10) {
	lumiWeight_mu = LumiWeightsUp_mu.ITweight( (int)event->nTruePU() );
	lumiWeight_el = LumiWeightsUp_el.ITweight( (int)event->nTruePU() );
      }
      else if (systematic == 11) {
	lumiWeight_mu = LumiWeightsDown_mu.ITweight( (int)event->nTruePU() );
	lumiWeight_el = LumiWeightsDown_el.ITweight( (int)event->nTruePU() );
      }

      scaleFactor_mu = scaleFactor_mu*lumiWeight_mu;
      scaleFactor_el = scaleFactor_el*lumiWeight_el;

      ///////////////////
      // TRIGGER SETUP //
      ///////////////////

      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename){
	previousFilename = currentFilename;
	iFile++;
	cout<<"File changed!!! => iFile = "<<iFile<<endl;
      }

      int currentRun = event->runId();

      if(previousRun != currentRun) {
        previousRun = currentRun;
	
	//semi-mu
	if(dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
	  
	  // 2.7/fb recalib 
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

	  else {
	    cout << "Unknown run for HLTpath selection: " << event->runId() << endl;
	    exit(1);
	  }

	  if( itriggerSemiMu == 9999 )
          {
            cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
            exit(-1);
          }

	} else {
	  itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun); // Summer12 DR53X
	  if( itriggerSemiMu == 9999 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun); // Summer12 DR53X
	}

	// semi-el
	// semi-electron
        if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0 ) {
         
	  
	  // 2.7/fb recalib 
	  if( event->runId() <= 190738 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v8"), currentRun, iFile);
	  
	  else if( event->runId() >= 190782 && event->runId() <= 191411 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v9"), currentRun, iFile);

	  else if( event->runId() >= 191695 && event->runId() <= 196531)
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun, iFile);
	  
	  else if( event->runId() >= 198049 && event->runId() <= 208686)
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v11"), currentRun, iFile);
  
	  else { 
            cout << "Unknown run for SemiEl HLTpath selection: " << event->runId() << endl;
	    exit(1);
	  }
	  if( itriggerSemiEl == 9999 )
	    {
	      cout << "itriggerSemiEl == 9999 for SemiEl HLTpath selection: " << event->runId() << endl;
	      exit(-1);
	    }
        }
        else
        {
	  itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun); // Summer12 DR53X
	  if( itriggerSemiEl == 9999 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v8"), currentRun); // Summer12 DR53X  
	}
	
      }


      if (itriggerSemiMu == 9999 && itriggerSemiEl == 9999) {
	  
	cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT IN RUN " << event->runId() << endl;

	exit(1);
	  
      }

      /////////////////////////////////////////////////////////////////////////////
      // JES SYSTEMATICS && SMEAR JET RESOLUTION TO MIMIC THE RESOLUTION IN DATA //
      /////////////////////////////////////////////////////////////////////////////

      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {

	/*if (event->eventId() == 5942766) {

	  cout << "bla" << endl;
	  
	  for (int j = 0; j < init_jets_corrected.size(); j++) {
	    
	    cout << "init jet " << j << " after JER -> Pt: " << init_jets_corrected[j]->Pt() << endl;
	    
	  } cout << endl;
	 
	  }*/

	if(systematic == 6)
	  jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus",false); //false means don't use old numbers but newer ones...
	else if(systematic == 7)
	  jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus",false);
	else
	  jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal",false);

	// Correct jets for JES uncertainy systematics
	
	/*for (unsigned int j=0;j<init_jets_corrected.size();j++)
	  cout << "jet " << j << " pt " << init_jets_corrected[j]->Pt() << endl;
	
	  cout << "MET " << mets[0]->Pt() << endl;*/
	
	if (systematic == 1)
	  jetTools->correctJetJESUnc(init_jets_corrected, mets[0], "minus",1);
	else if (systematic == 2)
	  jetTools->correctJetJESUnc(init_jets_corrected, mets[0], "plus",1);
	
      }

      /////////////////////
      // EVENT SELECTION //
      /////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho());
      selection.setJetCuts(40,2.5,0.01,1.,0.98,0.3,0.1);
      selection.setMuonCuts(25,2.1,0.12,0.2,0.3,1,0.5,5,0); // DR mu-jets cleaning still needed?
      
      //selection.setElectronCuts(32,2.5,0.1,0.02,0.5,0.3,0); // Ele25 uses 30 so Ele27 32? DR el-jets cleaning still needed?
      selection.setElectronCuts(32,2.5,0.1,0.02,0.9,0.3,0); // mva > 0.9

      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setLooseElectronCuts(20,2.5,0.15,0.);

      bool triggedSemiMu = false;
      bool triggedSemiEl = false;

      if( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) )
        triggedSemiMu = treeLoader.EventTrigged (itriggerSemiMu);
      if( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) )
        triggedSemiEl = treeLoader.EventTrigged (itriggerSemiEl);

      bool isGoodPV = false;
        
      //isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);
      isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);

      vector<TRootJet*> selectedJets, selectedJetsNoMu, selectedJetsNoEl;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons(10,2.5,0.2);
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(20,2.5,0.15);
      vector<TRootElectron*> vetoElectronsSemiEl = selection.GetSelectedLooseElectrons(20,2.5,0.15);

      vector<TRootMCParticle*> mcParticles;
      
      if(dataSetName.find("TTbarJets") == 0)
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }

      // reweigh pttop for madgraph

      double w_pttop = 1; 
      double w_pttop_top=1;
      double w_pttop_atop=1;
	
      //if (dataSetName.find("TTbarJets_madgraph_rew") != string::npos) {
      //if (dataSetName.find("TTbarJets") != string::npos && dataSetName.find("_norew") == string::npos && dataSetName.find("_mcnlo") == string::npos && dataSetName.find("_powheg") == string::npos) { // && dataSetName.find("_ttjets") == string::npos) {

	if (dataSetName.find("TTbarJets") != string::npos && dataSetName.find("_norew") == string::npos && dataSetName.find("_mcnlo") == string::npos){ // && dataSetName.find("_ttjets") == string::npos) {
	
	for (int p=0; p<mcParticles.size(); p++) {
	  
	  if (mcParticles[p]->status() != 3 ) continue;
	  if (mcParticles[p]->Pt() == 0 ) continue;
	  	  
	  double pt = mcParticles[p]->Pt();

	  //if (mcParticles[p]->type() == 6) w_pttop = (1.4*pow(10.0,-6.0)*pt*pt)-(2.0*pow(10.0,-3.0)*pt)+1.2;
	  
	  if (mcParticles[p]->type() == 6) w_pttop_top = exp(0.156-(0.00137*pt));
	  if (mcParticles[p]->type() == -6) w_pttop_atop = exp(0.156-(0.00137*pt));

	}

	w_pttop = sqrt(w_pttop_top*w_pttop_atop);

	//cout << w_pttop_top << " " << w_pttop_atop << " " << w_pttop << endl;

	scaleFactor_mu = scaleFactor_mu*w_pttop;
	scaleFactor_el = scaleFactor_el*w_pttop;
      }

      bool eventselectedSemiMu = false;
      bool eventselectedSemiEl = false;

      // semi-mu selection

      /*for (unsigned int j=0;j<selectedJets.size();j++)
	cout << "1) jet " << j << " pt " << selectedJets[j]->Pt() << endl;
	cout << "# jets selected " << selectedJets.size() << endl << endl;*/
      
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) && correctJets) {
	
	jetTools->scaleJets(init_jets_corrected,mets[0],jetScale_mu);
	
      }
	
      selectedJets = selection.GetSelectedJets(true);

      /*if (selectedJets.size() >= 4) {
	  if (selectedJets[0]->Pt() < 40) selectedJets.clear();
	  if (selectedJets[1]->Pt() < 40) selectedJets.clear();
	  if (selectedJets[2]->Pt() < 30) selectedJets.clear();
	  if (selectedJets[3]->Pt() < 30) selectedJets.clear();
	  }*/
      
      selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);

      if (triggedSemiMu) {
	if (isGoodPV) {
	  if(dataSetName.find("InvIso") != string::npos) selectedJets = selectedJetsNoMu;
	  if (selectedMuons.size() == 1) {
	    //cout << "selectedLooseMuons: " << selectedLooseMuons.size() << endl;
	    if( vetoMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && vetoMuons.size() == 0 ) ) { // if InvertedIso, selected muon not part of vetoMuons vector!
	      if (vetoElectronsSemiMu.size() == 0) {
		if (selectedJets.size() >= 4) {
		  //if (selectedJets.size() == 4) {
		  //if (selectedJets[0]->Pt() > 45 && selectedJets[1]->Pt() > 45 && selectedJets[2]->Pt() > 45 && selectedJets[3]->Pt() > 20) {
		  eventselectedSemiMu = true;
		  // }
		}
	      }
	    }
	  }
	}
      }


      // semi-el selection

      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) && correctJets && !eventselectedSemiMu) {
	
	jetTools->scaleJets(init_jets_corrected,mets[0],jetScale_el/jetScale_mu);
             
	selectedJets = selection.GetSelectedJets(true);

	/*if (selectedJets.size() >= 4) {
	  if (selectedJets[0]->Pt() < 40) selectedJets.clear();
	  if (selectedJets[1]->Pt() < 40) selectedJets.clear();
	  if (selectedJets[2]->Pt() < 30) selectedJets.clear();
	  if (selectedJets[3]->Pt() < 30) selectedJets.clear();
	  }*/

      }

      selectedElectrons = selection.GetSelectedElectrons(selectedJets);

      double lumiWeight = lumiWeight_el;
      scaleFactor = scaleFactor_el;

      selecTableSemiEl.Fill(d,0,scaleFactor*lumiWeight);

      if( triggedSemiEl) {
	if(dataSetName.find("InvIso") != string::npos) selectedJets = selectedJetsNoEl;
       selecTableSemiEl.Fill(d,1,scaleFactor*lumiWeight);
       if (isGoodPV ) {
	 selecTableSemiEl.Fill(d,2,scaleFactor*lumiWeight);
	 if( selectedElectrons.size() == 1 ) {
	   selecTableSemiEl.Fill(d,3,scaleFactor*lumiWeight);
	   if( vetoMuons.size() == 0 ) {
	     selecTableSemiEl.Fill(d,4,scaleFactor*lumiWeight);
	     if (vetoElectronsSemiEl.size() == 1) {
	       //if( !selection.foundZCandidate(selectedElectrons[0], vetoElectronsSemiEl) ) {
	       selecTableSemiEl.Fill(d,5,scaleFactor*lumiWeight);
	       //if( selection.passConversionRejection(selectedElectrons[0]) ) {
		 selecTableSemiEl.Fill(d,6,scaleFactor*lumiWeight);
		 if( selectedJets.size()>=1 ) {
		   selecTableSemiEl.Fill(d,7,scaleFactor*lumiWeight);
		   if( selectedJets.size()>=2 ) {
		     selecTableSemiEl.Fill(d,8,scaleFactor*lumiWeight);
		     if( selectedJets.size()>=3 ) {
		       selecTableSemiEl.Fill(d,9,scaleFactor*lumiWeight);
		       if( selectedJets.size()>=4 ) {
			 selecTableSemiEl.Fill(d,10,scaleFactor*lumiWeight);

			 //if (selectedJets[0]->Pt() > 45 && selectedJets[1]->Pt() > 45 && selectedJets[2]->Pt() > 45 && selectedJets[3]->Pt() > 45) {
			   eventselectedSemiEl=true;
			 //}
		       }
		     }
		   }
		 //}
	       }
	     }
	   }
	 }
       }
     }
     //eventselectedSemiEl=false;
     /*if (eventselected) {
	if (selectedJets[0]->Pt() < 70 || selectedJets[1]->Pt() < 50) 
	eventselected=false;
	}*/
     
     /*vector<TLorentzVector> triggerLeptons;
     std::map<std::string, std::vector<TopTree::triggeredObject> > trigFilters = event->getTriggerFilters();
     if( eventselectedSemiMu ) {
       if(trigFilters.find("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15") != trigFilters.end()) {
	 std::vector<TopTree::triggeredObject> trigObj = trigFilters["hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15"];
	 for(int i=0; i<trigObj.size(); i++) {
	   TLorentzVector tmp;
	   tmp.SetPtEtaPhiM(trigObj[i].pt, trigObj[i].eta, trigObj[i].phi, 0.);
	   triggerLeptons.push_back( tmp );
	 }
       }
       else if(trigFilters.find("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10") != trigFilters.end()) {
	 std::vector<TopTree::triggeredObject> trigObj = trigFilters["hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10"];
	 for(int i=0; i<trigObj.size(); i++) {
	   TLorentzVector tmp;
	   tmp.SetPtEtaPhiM(trigObj[i].pt, trigObj[i].eta, trigObj[i].phi, 0.);
	   triggerLeptons.push_back( tmp );
	 }
       }
       else cout << "No known triggerFilter found for run:  " << event->runId() << endl;

       bool matched=false;
       for(int i=0; i<triggerLeptons.size(); i++) {

	 //cout << i << " dR(lepton,triglepton) = " << selectedMuons[0]->DeltaR(triggerLeptons[i]) << endl << endl << endl;

	 if (selectedMuons[0]->DeltaR(triggerLeptons[i]) < 0.1) matched=true;

	 }

	 if (!matched) eventselectedSemiMu = false;

     }
     else if (eventselectedSemiEl) {
       if(trigFilters.find("hltEle27WP80TrackIsoFilter") != trigFilters.end()) {
	 std::vector<TopTree::triggeredObject> trigObj = trigFilters["hltEle27WP80TrackIsoFilter"];
	 for(int i=0; i<trigObj.size(); i++) {
	   TLorentzVector tmp;
	   tmp.SetPtEtaPhiM(trigObj[i].pt, trigObj[i].eta, trigObj[i].phi, 0.);
	   triggerLeptons.push_back( tmp );
	 }
       }
       else cout << "No known triggerFilter found for run:  " << event->runId() << endl;

       bool matched=false;
       for(int i=0; i<triggerLeptons.size(); i++) {
	 
	 //cout << i << " dR(lepton,triglepton) = " << selectedMuons[0]->DeltaR(triggerLeptons[i]) << endl << endl << endl;
	 
	 if (selectedElectrons[0]->DeltaR(triggerLeptons[i]) < 0.1) matched=true;
	 
       }
       
       if (!matched) eventselectedSemiEl = false;
     }
     */

     if (!eventselectedSemiMu && !eventselectedSemiEl) continue;

     /*
     /*int nBtags = 0;
     vector<float> bTagCSV;
     vector<TLorentzVector> otherSelectedJets;
     double bTagCutValue = 0.679; 
        
     for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++)
       {
	 otherSelectedJets.push_back( *selectedJets[iJet] );
	 bTagCSV.push_back(selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags());
	 if( selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags() > bTagCutValue )
	   {
	     nBtags++;
	   }
       }
     
       if(nBtags < 1) continue;*/

     if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
       
       if (systematic == 3)
	 jetTools->correctMETUnclusteredEnergy(mets[0],init_jets_corrected, init_muons, init_electrons, "minus");
       else if (systematic == 4)
	 jetTools->correctMETUnclusteredEnergy(mets[0],init_jets_corrected, init_muons, init_electrons, "plus"); 

     }

     if (eventselectedSemiMu && eventselectedSemiEl) cout << "selected in both channels are we?????" << endl;
     
     if (eventselectedSemiMu) scaleFactor = scaleFactor_mu;
     if (eventselectedSemiEl) scaleFactor = scaleFactor_el;

     /*if (scaleFactor_mu != scaleFactor_el) {
       cout << endl << "Mu: " << eventselectedSemiMu << " El: " << eventselectedSemiEl << endl;
       cout << "Scalefactor mu " << scaleFactor_mu << " el " << scaleFactor_el << endl;
       cout << "Scalefactor sel " << scaleFactor<< endl;
       }*/

     /*if (event->eventId() == 5942766) {

       cout << "bla" << endl;

       for (int j = 0; j < init_jets_corrected.size(); j++) {
	 
	 cout << "init jet " << j << " after JER -> Pt: " << init_jets_corrected[j]->Pt() << endl;
	 
       } cout << endl;

       for (int j = 0; j < selectedJets.size(); j++) {
	 
	 cout << "selected jet " << j << " after JER -> Pt: " << selectedJets[j]->Pt() << endl;
	 
       } cout << endl;
       
       exit(1);
       }*/

     // get MC QCD yield
     
     if (eventselectedSemiEl && dataSetName.find("QCD") != string::npos) {
       yieldQCD_MC_SemiEl += (LuminosityEl*scaleFactor)/datasets[d]->EquivalentLumi();
     }

     // count the events in each runrange bin for data driven QCD
      
      if (event->runId() <= 165633) // noniso triggers
	QCDEstcounts[0]++;
      else if (event->runId() >= 165970 && event->runId() < 175860 ) // iso trigger
	QCDEstcounts[1]++;
      else if (event->runId() >= 175860) // 2011B -> high PU
	QCDEstcounts[2]++;

      //for (unsigned int j=0;j<selectedJets.size();j++)
	//cout << "after event selection jet " << j << " pt " << selectedJets[j]->Pt() << " eta " << selectedJets[j]->Eta() << endl;
      //cout << "****** end event *******"<< endl;
      //cout << endl;
      //exit(1);

     
      /*if(dataSetName.find("WJets_TuneZ2_scaledown") == 0) {
	treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class

	bool muFound=false;
	for (int o=0; o<mcParticles.size();o++) {

	  //cout << mcParticles[o]->type() << " " << mcParticles[o]->motherType() << " " << endl;
	  if (fabs( mcParticles[o]->type() ) == 13 && fabs( mcParticles[o]->motherType() ) == 24)
	    muFound=true;

	}

	// exit(1);


	if (muFound) {

	  scaleFactor = scaleFactor*0.75;
	  
	  //cout << "--- W+Jets scaledown found mu event!!!!"<< endl;

	} //else {
	  //cout << "--- W+Jets scaledown NOT a MU!!!!"<< endl;
	  //exit(1);
	//}
	}*/

      if (eventselectedSemiMu)
	nSelectedMu++;
      if (eventselectedSemiEl)
	nSelectedEl++;

      TLorentzVector* selectedLepton;

      if (eventselectedSemiMu)
	selectedLepton = (TLorentzVector*)selectedMuons[0];
      else if (eventselectedSemiEl)
	selectedLepton = (TLorentzVector*)selectedElectrons[0];

      ////////////////////////////////////
      // JET PARTON MATCHING FOR MASSES //
      ////////////////////////////////////

      // this is only necessary for the chi2 jetcomb input masses. if usemassesandresolutions == true this will not be run

      int MCPermutation[4]; 

      bool all4PartonsMatched = false; // True if the 4 ttbar semi-lep partons are matched to 4 jets (not necessarily the 4 highest pt jets)
      bool all4JetsMatched_MCdef_ = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
      bool hadronictopJetsMatched_MCdef_ = false;
      
      pair<unsigned int, unsigned int> leptonicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
      pair<unsigned int, unsigned int> hadronicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
      pair<unsigned int, unsigned int> hadronicWJet1_ = pair<unsigned int,unsigned int>(9999,9999);
      pair<unsigned int, unsigned int> hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
      
      double relDiffEJetParton_b_ = -9999;
      double relDiffEJetParton_l1_ = -9999;
      double relDiffEJetParton_l2_ = -9999;

      int pdgID_top = 6; //top quark

      if (!useMassesAndResolutions && dataSetName.find("TTbarJets") == 0 && (isSemiMu || isSemiE) ) {
	
	sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)     
		
	vector<TRootMCParticle*> mcParticlesMatching_;
	vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
	TLorentzVector topQuark, antiTopQuark;
	
	bool muPlusFromTop = false, muMinusFromTop = false, elPlusFromTop = false, elMinusFromTop = false;
	int nTTbarQuarks = 0;

	mcParticlesMatching_.clear();
	
	for(unsigned int i=0; i<mcParticles.size(); i++) {
	  //cout << i << ":  status: " << mcParticles[i]->status() << "  pdgId: " << mcParticles[i]->type()
	  //  << "  motherPdgId: " << mcParticles[i]->motherType() << "  grannyPdgId: " << mcParticles[i]->grannyType() << endl;
	  if( mcParticles[i]->status() != 3) continue;
	  if( mcParticles[i]->Pt() < 0.00000001) continue;
	  
	  if( mcParticles[i]->type() == pdgID_top )
	    topQuark = *mcParticles[i];
	  else if( mcParticles[i]->type() == -pdgID_top )
	    antiTopQuark = *mcParticles[i];
	  
	  //madgraph
	  /*if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )
	    muMinusFromTop = true;
	  if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )
	    muPlusFromTop = true;
	  if( mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )
	    elMinusFromTop = true;
	  if( mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )
	  elPlusFromTop = true;*/

	  //mc@nlo
	  if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24  )
	    muMinusFromTop = true;
	  if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24  )
	    muPlusFromTop = true;
	  if( mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24  )
	    elMinusFromTop = true;
	  if( mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24  )
	  elPlusFromTop = true;
	  
	  if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ) {  //light/b quarks, 6 should stay hardcoded
	    mcParticlesTLV.push_back(*mcParticles[i]);
	    mcParticlesMatching_.push_back(mcParticles[i]);
	    
	  }
	}
	  
	// take all the selectedJets_ to study the radiation stuff, selectedJets_ are already ordened in decreasing Pt()
	for(unsigned int i=0; i<selectedJets.size(); i++)
	  selectedJetsTLV.push_back(*selectedJets[i]);
	
	JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
	
	if(matching.getNumberOfAvailableCombinations() != 1)
	  cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
	
	vector< pair<unsigned int, unsigned int> > JetPartonPair, ISRJetPartonPair; // First one is jet number, second one is mcParticle number
	
	for(unsigned int i=0; i<mcParticlesTLV.size(); i++) {
	  int matchedJetNumber = matching.getMatchForParton(i, 0);
	  if(matchedJetNumber != -1)
	    JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
	}
	
	for(unsigned int i=0; i<JetPartonPair.size(); i++) {
	  unsigned int j = JetPartonPair[i].second;
	  
	  if( fabs(mcParticlesMatching_[j]->type()) < 6 ) {//light/b quarks, 6 should stay hardcoded
	    //if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -pdgID_top )
	    //	|| ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == pdgID_top ) ) {
	    if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == -24 )
		|| ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == 24  ) ) {
	      if(hadronicWJet1_.first == 9999) {
		hadronicWJet1_ = JetPartonPair[i];
		MCPermutation[0] = JetPartonPair[i].first;
	      } else if(hadronicWJet2_.first == 9999) {
		hadronicWJet2_ = JetPartonPair[i];
		MCPermutation[1] = JetPartonPair[i].first;
	      } else {
		cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
		cerr<<" -- isSemiMu: " << isSemiMu << " isSemiE: " << isSemiE << endl;
		cerr<<" -- muMinusFromMtop: " << muMinusFromTop << " muPlusFromMtop: " << muPlusFromTop << endl;
		cerr<<" -- pdgId: " << mcParticlesMatching_[j]->type() << " mother: " << mcParticlesMatching_[j]->motherType() << " granny: " << mcParticlesMatching_[j]->grannyType() << " Pt: " << mcParticlesMatching_[j]->Pt()<< endl;
		//exit(1);
	      }
	    }
	  }
	  if( fabs(mcParticlesMatching_[j]->type()) == 5 ) {
	    if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == -pdgID_top )
		|| ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == pdgID_top ) ) {
	      hadronicBJet_ = JetPartonPair[i];
	      MCPermutation[2] = JetPartonPair[i].first;
	    }
	    else if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticlesMatching_[j]->motherType() == pdgID_top )
		     || ( ( muMinusFromTop || elMinusFromTop ) && mcParticlesMatching_[j]->motherType() == -pdgID_top ) ) {
	      leptonicBJet_ = JetPartonPair[i];
	      MCPermutation[3] = JetPartonPair[i].first;
	    }
	  }
	  
	  // look for ISR stuff
	  //if( fabs(mcParticlesMatching_[j]->type()) != pdgID_top && fabs(mcParticlesMatching_[j]->motherType()) != 24 && fabs(mcParticlesMatching_[j]->motherType()) != pdgID_top &&
	  //    fabs(mcParticlesMatching_[j]->grannyType()) != 24 && fabs(mcParticlesMatching_[j]->grannyType()) != pdgID_top )
	  if( fabs(mcParticlesMatching_[j]->type()) != pdgID_top && fabs(mcParticlesMatching_[j]->motherType()) != 24 && fabs(mcParticlesMatching_[j]->motherType()) != pdgID_top)

	    {
	      ISRJetPartonPair.push_back(JetPartonPair[i]);
	    }
	}
	
	if(hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999 && hadronicBJet_.first != 9999 && leptonicBJet_.first != 9999) {
	  
	  all4PartonsMatched = true;
	  if(hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4 && leptonicBJet_.first < 4)
	    all4JetsMatched_MCdef_ = true;
	}
	////cout<<"   ------> according to JetCombiner: hadronicWJet1_.first = "<<hadronicWJet1_.first<<", hadronicWJet2_.first = "<<hadronicWJet2_.first<<", hadronicBJet_.first = "<<hadronicBJet_.first<<endl;
	if(hadronicWJet1_.first < 4 && hadronicWJet2_.first < 4 && hadronicBJet_.first < 4)
	  hadronictopJetsMatched_MCdef_ = true;
	
      }

      /*cout << "isSemiMu: " << isSemiMu;
      cout << " isSemiEl: " << isSemiE;
      cout << " SelSemiMu: "<<eventselectedSemiMu;
      cout << " SelSemiEl: " << eventselectedSemiEl;
      cout << " All 4 partons matched: " << all4JetsMatched_MCdef_ << endl;*/

      //if (eventselectedSemiEl) 
      //	exit(0);
    
      if (all4JetsMatched_MCdef_ && !useMassesAndResolutions && dataSetName.find("TTbarJets") == 0) {
	  	
	float WMass = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]).M();
	float TopMass = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]+*selectedJets[hadronicBJet_.first]).M();
	
	//histo1D["hadronicRecoWMass"]->Fill(WMass,lumiWeight);
	//histo1D["hadronicRecoTopMass"]->Fill(TopMass,lumiWeight);
	
	histo1D["hadronicRecoWMass"]->Fill(WMass);
	histo1D["hadronicRecoTopMass"]->Fill(TopMass);
	
	histo1D["hadronicRecoWMass-REW"]->Fill(WMass,lumiWeight);
	histo1D["hadronicRecoTopMass-REW"]->Fill(TopMass,lumiWeight);
	
	histo1D["hadronicRecoWMass-REW3D"]->Fill(WMass,lumiWeight);
	histo1D["hadronicRecoTopMass-REW3D"]->Fill(TopMass,lumiWeight);
	
	  // to be fixed

	//cout << "WMass " << WMass << " TopMass " << TopMass << endl; exit(0);
	
      }

      ///////////////////////////////////
      // !!! FILLING THE BTAGTREES !!! //
      ///////////////////////////////////

      // now we use the loaded values and do the chi2 jetcombination + fill the trees
      if (useMassesAndResolutions) {

	//cout << "Chi2Topmass -- " << Chi2Topmass << endl;
	// index convention -> i,j: jets from Hadronic W  k: Hadronic b and l: leptonic b
	float smallestChi2 = 9999999999.;
	int Permutation[4];
	for (unsigned int i=0; i<4; i++) {
  	  for (unsigned int j=0; j<4; j++) {
   	    for (unsigned int k=0; k<4; k++) {
   	      for (unsigned int l=0; l<4; l++) {
		if (i < j && i != j && i != k && i != l && j != k && j != l && k != l) {

		  //cout << i << " " << j << " " << k << " " << l << endl;

		  float Chi2WmassRec = (*selectedJets[i]+*selectedJets[j]).M();
		  float Chi2TopmassRec = (*selectedJets[i]+*selectedJets[j]+*selectedJets[k]).M();

		  float termW = pow(float(Chi2WmassRec-Chi2Wmass)/float(SigmaChi2Wmass),float(2));
		  float termTop = pow(float(Chi2TopmassRec-Chi2Topmass)/float(SigmaChi2Topmass),float(2));
		  float chi2 = termW+termTop;

		  //cout << chi2 << " " << i << " " << j << " " << k << " " << l << endl;

		  if (chi2 < smallestChi2) {
		    smallestChi2 = chi2;
		    Permutation[0]=i;
		    Permutation[1]=j;
		    Permutation[2]=k;
		    Permutation[3]=l;
		  }
		}
	      }
	    }
	  }
	}
	//cout << "smallest " << smallestChi2 << " " << Permutation[0] << " " << Permutation[1] << " " << Permutation[2] << " " << Permutation[3] << endl;

	// PERFORMANCE CHECK

	nAll++;
	if (all4JetsMatched_MCdef_ ) {
	  //cout << "first yeah " << endl;
	  nGoodCombi++;

	  if ( (Permutation[0] == hadronicWJet1_.first && Permutation[1] == hadronicWJet2_.first) || (Permutation[1] == hadronicWJet1_.first && Permutation[0] == hadronicWJet2_.first)) {
	    
	    if (Permutation[2] == hadronicBJet_.first) {

	      if (Permutation[3] == leptonicBJet_.first) {

		//cout << "O yeah" << endl;

		nGoodChi2++;

	      }

	    }

	  }
	}


	///////////////////////
	// CHECK OF BTAG EFF //
	///////////////////////

	nJets++;

	if (fabs(selectedJets[Permutation[3]]->partonFlavour()) == 5)
	  nBjets++;

	if (selectedJets[Permutation[3]]->btag_jetProbabilityBJetTags() > 0.545) {
	  nTagged++;
	  if (fabs(selectedJets[Permutation[3]]->partonFlavour()) == 5)
	    nTagged_trueb++;
	  else 
	    nTagged_nonb++;
	}

	//----------------------------//
	// PLOTS TO CHECK SYSTEMATICS //
	//----------------------------//

	double systw = (10000./datasets[d]->EquivalentLumi())*scaleFactor;

	if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // when only using one channel

	  systw=1;

	}

	if (histo1D_syst.find("weight") == histo1D_syst.end()) {
	  histo1D_syst["weight"] = new TH1D("weight_check","weight check;weight;nEvents",100,0,1);
	}
	histo1D_syst["weight"]->Fill(systw);
	
	//cout << systw << endl;

	// FOR JEC we want to look at mW, trying different definitions

	/* using the chi2 definition */

	float mW_chi2 = (*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]).M();

	if (histo1D_syst.find("JEC_Wmass_chi2def_mu") == histo1D_syst.end()) {
	  histo1D_syst["JEC_Wmass_chi2def_mu"] = new TH1D("JEC_Wmass_chi2def_mu","M_{W} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{W} (GeV/c^{2});#events",80,20,180);
	  histo1D_syst["JEC_Wmass_chi2def_el"] = new TH1D("JEC_Wmass_chi2def_el","M_{W} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{W} (GeV/c^{2});#events",80,20,180);
	}
	
	if (eventselectedSemiMu)
	  histo1D_syst["JEC_Wmass_chi2def_mu"]->Fill(mW_chi2,systw);
	if (eventselectedSemiEl)
	  histo1D_syst["JEC_Wmass_chi2def_el"]->Fill(mW_chi2,systw);

	//if (eventselectedSemiMu && event->eventId() == 7845316) {
	//  cout << event->eventId() << " " << mW_chi2 << endl;
	//  exit(1);
	//}

	/* using the non-tagged jets in a 2-btag event */

	int nBtags=0;
	vector<int> notTagged;
	for (int bt=0;bt<selectedJets.size();bt++) {
	  if (selectedJets[bt]->btag_jetProbabilityBJetTags() > 0.545)
	    nBtags++;
	  else 
	    notTagged.push_back(bt);
	}

	if (histo1D_syst.find("JEC_Wmass_doublebtagdef_mu") == histo1D_syst.end()) {
	  histo1D_syst["JEC_Wmass_doublebtagdef_mu"] = new TH1D("JEC_Wmass_doublebtagdef_mu","M_{W} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{W} (GeV/c^{2});#events",80,20,180);
	  histo1D_syst["JEC_Wmass_doublebtagdef_el"] = new TH1D("JEC_Wmass_doublebtagdef_el","M_{W} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{W} (GeV/c^{2});#events",80,20,180);

	  for (int m=5;m<45;m+=5) {

	    stringstream cutmet; cutmet << m;
	    string title_mu="JEC_Wmass_doublebtagdef_METCUT_"+cutmet.str()+"_mu";
	    string title_el="JEC_Wmass_doublebtagdef_METCUT_"+cutmet.str()+"_el";

	    histo1D_syst[title_mu.c_str()] = new TH1D(title_mu.c_str(),"M_{W} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{W} (GeV/c^{2});#events",80,20,180);
	    histo1D_syst[title_el.c_str()] = new TH1D(title_el.c_str(),"M_{W} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{W} (GeV/c^{2});#events",80,20,180);
	    
	  }
	}
	
	//cout << nBtags << " " << notTagged.size() << endl;

	if (nBtags == 2) {
	  for (int l1=0; l1<notTagged.size();l1++) {
	    for (int l2=0; l2<notTagged.size();l2++) {
	      if (l1==l2) continue;
	      if (l1>l2) continue;
	      
	      double mw_doublebtag=(*selectedJets[notTagged[l1]]+*selectedJets[notTagged[l2]]).M();
	      //cout << "     " << l1 << " " << l2 << " " << mw_doublebtag << endl;
	      if (eventselectedSemiMu)
		histo1D_syst["JEC_Wmass_doublebtagdef_mu"]->Fill(mw_doublebtag,systw);
	      if (eventselectedSemiEl)
		histo1D_syst["JEC_Wmass_doublebtagdef_el"]->Fill(mw_doublebtag,systw);
	      
	      for (int m=5;m<45;m+=5) {

		stringstream cutmet; cutmet << m;
		string title_mu="JEC_Wmass_doublebtagdef_METCUT_"+cutmet.str()+"_mu";
		string title_el="JEC_Wmass_doublebtagdef_METCUT_"+cutmet.str()+"_el";
		
		if (mets[0]->Pt() > m) {
		  if (eventselectedSemiMu)
		    histo1D_syst[title_mu.c_str()]->Fill(mw_doublebtag,systw);
		  if (eventselectedSemiEl)
		    histo1D_syst[title_el.c_str()]->Fill(mw_doublebtag,systw);
		}
	      }
	    }
	  }
	}

	notTagged.clear();
	//if (nBtags == 2 && notTagged.size() > 2) exit(1);


	//---------------//
	// CONTROL PLOTS //
	//---------------//

	double controlw = (10000./datasets[d]->EquivalentLumi())*scaleFactor;

	if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // when only using one channel

	  controlw=1;

	}

	if (histo1D_btag.find("weight") == histo1D_btag.end()) {
	  histo1D_btag["weight"] = new TH1F("weight_check","weight check;weight;nEvents",100,0,1);
	}
	histo1D_btag["weight"]->Fill(controlw);

	/*if (histo1D_btag.find("CSVTag_SS_Left_mu") == histo1D_btag.end()) {
	  histo1D_btag["CSVTag_SS_Left_mu"] = new TH1F("CSVTag_SS_Left_mu","CSVTag_SS_Left_Mu;CSV;#events",25,0,1);
	  histo1D_btag["CSVTag_SS_Right_mu"] = new TH1F("CSVTag_SS_Right_mu","CSVTag_SS_Right_Mu;CSV;#events",25,0,1);
	  histo1D_btag["CSVTag_SS_Left_b_mu"] = new TH1F("CSVTag_SS_Left_b_mu","CSVTag_SS_Left_Mu;CSV;#events",25,0,1);
	  histo1D_btag["CSVTag_SS_Right_b_mu"] = new TH1F("CSVTag_SS_Right_b_mu","CSVTag_SS_Right_Mu;CSV;#events",25,0,1);
	  histo1D_btag["CSVTag_SS_Left_l_mu"] = new TH1F("CSVTag_SS_Left_l_mu","CSVTag_SS_Left_Mu;CSV;#events",25,0,1);
	  histo1D_btag["CSVTag_SS_Right_l_mu"] = new TH1F("CSVTag_SS_Right_l_mu","CSVTag_SS_Right_Mu;CSV;#events",25,0,1);

	  histo1D_btag["CSVTag_CS_Left_mu"] = new TH1F("CSVTag_CS_Left_mu","CSVTag_CS_Left_Mu;CSV;#events",25,0,1);
	  histo1D_btag["CSVTag_CS_Right_mu"] = new TH1F("CSVTag_CS_Right_mu","CSVTag_CS_Right_Mu;CSV;#events",25,0,1);
	  histo1D_btag["CSVTag_CS_Left_b_mu"] = new TH1F("CSVTag_CS_Left_b_mu","CSVTag_CS_Left_Mu;CSV;#events",25,0,1);
	  histo1D_btag["CSVTag_CS_Right_b_mu"] = new TH1F("CSVTag_CS_Right_b_mu","CSVTag_CS_Right_Mu;CSV;#events",25,0,1);
	  histo1D_btag["CSVTag_CS_Left_l_mu"] = new TH1F("CSVTag_CS_Left_l_mu","CSVTag_CS_Left_Mu;CSV;#events",25,0,1);
	  histo1D_btag["CSVTag_CS_Right_l_mu"] = new TH1F("CSVTag_CS_Right_l_mu","CSVTag_CS_Right_Mu;CSV;#events",25,0,1);
	}
	if (eventselectedSemiMu && smallestChi2 < 125) {

	  //cout << "ok1" << endl;

	    double mlj = (*selectedJets[3]+*selectedLepton).M();
	    double mljC1 = (*(selectedJets[Permutation[0]])+*selectedLepton).M();
	    double mljC2 = (*(selectedJets[Permutation[1]])+*selectedLepton).M();

	    int pflav = selectedJets[Permutation[3]]->partonFlavour();
	    int pflavC1 = selectedJets[Permutation[0]]->partonFlavour();
	    int pflavC2 = selectedJets[Permutation[1]]->partonFlavour();

	    if (mlj > 70 && mlj < 170) {
	      //cout << "ok2" << endl;
	      histo1D_btag["CSVTag_SS_Left_mu"]->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      if (fabs(pflav)==5)
		histo1D_btag["CSVTag_SS_Left_b_mu"]->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      else
		histo1D_btag["CSVTag_SS_Left_l_mu"]->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	    }
	    else if (mlj > 170 && mlj < 300) {
	      histo1D_btag["CSVTag_SS_Right_mu"]->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      if (fabs(pflav)==5)
		histo1D_btag["CSVTag_SS_Right_b_mu"]->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      else
		histo1D_btag["CSVTag_SS_Right_l_mu"]->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	    }

	    if (mljC1 > 70 && mljC1 < 170) {
	      histo1D_btag["CSVTag_CS_Left_mu"]->Fill((*(selectedJets[Permutation[0]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      if (fabs(pflavC1)==5)
		histo1D_btag["CSVTag_CS_Left_b_mu"]->Fill((*(selectedJets[Permutation[0]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      else
		histo1D_btag["CSVTag_CS_Left_l_mu"]->Fill((*(selectedJets[Permutation[0]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	    }
	    else if (mljC1 > 170 && mljC1 < 300) {
	      histo1D_btag["CSVTag_CS_Right_mu"]->Fill((*(selectedJets[Permutation[0]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      if (fabs(pflavC1)==5)
		histo1D_btag["CSVTag_CS_Right_b_mu"]->Fill((*(selectedJets[Permutation[0]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      else
		histo1D_btag["CSVTag_CS_Right_l_mu"]->Fill((*(selectedJets[Permutation[0]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	    }

	    if (mljC2 > 70 && mljC2 < 170) {
	      histo1D_btag["CSVTag_CS_Left_mu"]->Fill((*(selectedJets[Permutation[1]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      if (fabs(pflavC2)==5)
		histo1D_btag["CSVTag_CS_Left_b_mu"]->Fill((*(selectedJets[Permutation[1]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      else
		histo1D_btag["CSVTag_CS_Left_l_mu"]->Fill((*(selectedJets[Permutation[1]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	    }
	    else if (mljC2 > 170 && mljC2 < 300) {
	      histo1D_btag["CSVTag_CS_Right_mu"]->Fill((*(selectedJets[Permutation[1]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      if (fabs(pflavC2)==5)
		histo1D_btag["CSVTag_CS_Right_b_mu"]->Fill((*(selectedJets[Permutation[1]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	      else
		histo1D_btag["CSVTag_CS_Right_l_mu"]->Fill((*(selectedJets[Permutation[1]])).btag_combinedSecondaryVertexRetrainedBJetTags());
	    }
	    
	}*/

	if (histo1D_btag.find("JetCombMinChi2_mu") == histo1D_btag.end()) {
	  histo1D_btag["JetCombMinChi2_mu"] = new TH1F("JetCombMinChi2_mu","Minimal #Chi^{2} value of the 12 combinations;#Chi^{2};#events",200,0,100);
	  histo1D_btag["JetCombMinChi2_el"] = new TH1F("JetCombMinChi2_el","Minimal #Chi^{2} value of the 12 combinations;#Chi^{2};#events",200,0,100);
	}
	if (eventselectedSemiMu)
	  histo1D_btag["JetCombMinChi2_mu"]->Fill(smallestChi2,controlw); 
	if (eventselectedSemiEl)
	  histo1D_btag["JetCombMinChi2_el"]->Fill(smallestChi2,controlw);

	if (histo1D_btag.find("BestJetCombMLB_mu") == histo1D_btag.end()) {
	  histo1D_btag["BestJetCombMLB_mu"] = new TH1F("BestJetCombMLB_mu","M_{lb} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{lb} (GeV/c^{2});#events",600,0,1200);
	  histo1D_btag["BestJetCombMLB_el"] = new TH1F("BestJetCombMLB_el","M_{lb} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{lb} (GeV/c^{2});#events",600,0,1200);
	}
	
	if (eventselectedSemiMu)
	  histo1D_btag["BestJetCombMLB_mu"]->Fill((*selectedLepton+*selectedJets[Permutation[3]]).M(),controlw);
	if (eventselectedSemiEl)
	  histo1D_btag["BestJetCombMLB_el"]->Fill((*selectedLepton+*selectedJets[Permutation[3]]).M(),controlw);

	if (histo2D_btag.find("Top_VS_W_Mass_mu") == histo2D_btag.end()) {
	  histo2D_btag["Top_VS_W_Mass_mu"] = new TH2F("Top_VS_W_Mass_mu","M_{top} VS m_{W} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{top} (GeV/c^{2});m_{W} (GeV/c^{2})",200,100,500,200,0,400);
	  histo2D_btag["Top_VS_W_Mass_el"] = new TH2F("Top_VS_W_Mass_el","M_{top} VS m_{W} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{top} (GeV/c^{2});m_{W} (GeV/c^{2})",200,100,500,200,0,400);
	}
	float mtop = (*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]+*selectedJets[Permutation[2]]).M();
	float mW = (*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]).M();
	
	if (eventselectedSemiMu)
	  histo2D_btag["Top_VS_W_Mass_mu"]->Fill(mtop,mW);
	if (eventselectedSemiEl)
	  histo2D_btag["Top_VS_W_Mass_el"]->Fill(mtop,mW);

	if (histo1D_btag.find("BestJetCombMTop_mu") == histo1D_btag.end()) {
	  histo1D_btag["BestJetCombMTop_mu"] = new TH1F("BestJetCombMTop_mu","M_{top} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{top} (GeV/c^{2});#events",80,100,500);
	  histo1D_btag["BestJetCombMTop_el"] = new TH1F("BestJetCombMTop_el","M_{top} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{top} (GeV/c^{2});#events",80,100,500);
	}
	
	if (eventselectedSemiMu)
	  histo1D_btag["BestJetCombMTop_mu"]->Fill(mtop,controlw);
	if (eventselectedSemiEl)
	  histo1D_btag["BestJetCombMTop_el"]->Fill(mtop,controlw);

	if (histo1D_btag.find("pt_jet_0_mu") == histo1D_btag.end()) {

	  // all tags

	  histo1D_btag["pt_jet_0_mu"] = new TH1F("pt_jet_0_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_0_el"] = new TH1F("pt_jet_0_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_jet_1_mu"] = new TH1F("pt_jet_1_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_1_el"] = new TH1F("pt_jet_1_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_jet_2_mu"] = new TH1F("pt_jet_2_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_2_el"] = new TH1F("pt_jet_2_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_jet_3_mu"] = new TH1F("pt_jet_3_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_3_el"] = new TH1F("pt_jet_3_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_lepton_mu"] = new TH1F("pt_lepton_mu","muon pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_lepton_el"] = new TH1F("pt_lepton_el","electron pt;p_{T} (GeV);#events",500,0,500);

	  // zero tag only

	  histo1D_btag["pt_jet_0_zerotag_mu"] = new TH1F("pt_jet_0_zerotag_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_0_zerotag_el"] = new TH1F("pt_jet_0_zerotag_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_jet_1_zerotag_mu"] = new TH1F("pt_jet_1_zerotag_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_1_zerotag_el"] = new TH1F("pt_jet_1_zerotag_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_jet_2_zerotag_mu"] = new TH1F("pt_jet_2_zerotag_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_2_zerotag_el"] = new TH1F("pt_jet_2_zerotag_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_jet_3_zerotag_mu"] = new TH1F("pt_jet_3_zerotag_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_3_zerotag_el"] = new TH1F("pt_jet_3_zerotag_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_lepton_zerotag_mu"] = new TH1F("pt_lepton_zerotag_mu","muon pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_lepton_zerotag_el"] = new TH1F("pt_lepton_zerotag_el","electron pt;p_{T} (GeV);#events",500,0,500);

	  // >= 1 tag

	  histo1D_btag["pt_jet_0_gteq1tag_mu"] = new TH1F("pt_jet_0_gteq1tag_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_0_gteq1tag_el"] = new TH1F("pt_jet_0_gteq1tag_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_jet_1_gteq1tag_mu"] = new TH1F("pt_jet_1_gteq1tag_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_1_gteq1tag_el"] = new TH1F("pt_jet_1_gteq1tag_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_jet_2_gteq1tag_mu"] = new TH1F("pt_jet_2_gteq1tag_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_2_gteq1tag_el"] = new TH1F("pt_jet_2_gteq1tag_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_jet_3_gteq1tag_mu"] = new TH1F("pt_jet_3_gteq1tag_mu","jet pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_jet_3_gteq1tag_el"] = new TH1F("pt_jet_3_gteq1tag_el","jet pt;p_{T} (GeV);#events",500,0,500);

	  histo1D_btag["pt_lepton_gteq1tag_mu"] = new TH1F("pt_lepton_gteq1tag_mu","muon pt;p_{T} (GeV);#events",500,0,500);
	  histo1D_btag["pt_lepton_gteq1tag_el"] = new TH1F("pt_lepton_gteq1tag_el","electron pt;p_{T} (GeV);#events",500,0,500);
	  
	}
	
	int nTagsControl=0;
	for (int nt=0; nt< selectedJets.size(); nt++) {

	  if (selectedJets[nt]->btag_combinedSecondaryVertexRetrainedBJetTags() > 0.679)
	    nTagsControl++;

	}

	//cout << nTagsControl << endl;

	if (eventselectedSemiMu) {
	  histo1D_btag["pt_jet_0_mu"]->Fill(selectedJets[0]->Pt(),controlw);
	  histo1D_btag["pt_jet_1_mu"]->Fill(selectedJets[1]->Pt(),controlw);
	  histo1D_btag["pt_jet_2_mu"]->Fill(selectedJets[2]->Pt(),controlw);
	  histo1D_btag["pt_jet_3_mu"]->Fill(selectedJets[3]->Pt(),controlw);
	  histo1D_btag["pt_lepton_mu"]->Fill(selectedMuons[0]->Pt(),controlw);
	  if (nTagsControl == 0) {
	    histo1D_btag["pt_jet_0_zerotag_mu"]->Fill(selectedJets[0]->Pt(),controlw);
	    histo1D_btag["pt_jet_1_zerotag_mu"]->Fill(selectedJets[1]->Pt(),controlw);
	    histo1D_btag["pt_jet_2_zerotag_mu"]->Fill(selectedJets[2]->Pt(),controlw);
	    histo1D_btag["pt_jet_3_zerotag_mu"]->Fill(selectedJets[3]->Pt(),controlw);
	    histo1D_btag["pt_lepton_zerotag_mu"]->Fill(selectedMuons[0]->Pt(),controlw);
	  } else if (nTagsControl >= 1) {
	    histo1D_btag["pt_jet_0_gteq1tag_mu"]->Fill(selectedJets[0]->Pt(),controlw);
	    histo1D_btag["pt_jet_1_gteq1tag_mu"]->Fill(selectedJets[1]->Pt(),controlw);
	    histo1D_btag["pt_jet_2_gteq1tag_mu"]->Fill(selectedJets[2]->Pt(),controlw);
	    histo1D_btag["pt_jet_3_gteq1tag_mu"]->Fill(selectedJets[3]->Pt(),controlw);
	    histo1D_btag["pt_lepton_gteq1tag_mu"]->Fill(selectedMuons[0]->Pt(),controlw);
	  }
	}
	if (eventselectedSemiEl) {
	  histo1D_btag["pt_jet_0_el"]->Fill(selectedJets[0]->Pt());
	  histo1D_btag["pt_jet_1_el"]->Fill(selectedJets[1]->Pt());
	  histo1D_btag["pt_jet_2_el"]->Fill(selectedJets[2]->Pt());
	  histo1D_btag["pt_jet_3_el"]->Fill(selectedJets[3]->Pt());
	  histo1D_btag["pt_lepton_el"]->Fill(selectedElectrons[0]->Pt());	  
	  if (nTagsControl == 0) {
	    histo1D_btag["pt_jet_0_zerotag_el"]->Fill(selectedJets[0]->Pt());
	    histo1D_btag["pt_jet_1_zerotag_el"]->Fill(selectedJets[1]->Pt());
	    histo1D_btag["pt_jet_2_zerotag_el"]->Fill(selectedJets[2]->Pt());
	    histo1D_btag["pt_jet_3_zerotag_el"]->Fill(selectedJets[3]->Pt());
	    histo1D_btag["pt_lepton_zerotag_el"]->Fill(selectedElectrons[0]->Pt());
	  } else if (nTagsControl >= 1) {
	    histo1D_btag["pt_jet_0_gteq1tag_el"]->Fill(selectedJets[0]->Pt());
	    histo1D_btag["pt_jet_1_gteq1tag_el"]->Fill(selectedJets[1]->Pt());
	    histo1D_btag["pt_jet_2_gteq1tag_el"]->Fill(selectedJets[2]->Pt());
	    histo1D_btag["pt_jet_3_gteq1tag_el"]->Fill(selectedJets[3]->Pt());
	    histo1D_btag["pt_lepton_gteq1tag_el"]->Fill(selectedElectrons[0]->Pt());
	  }
	}

	//-----------------//
	// do some data-mc //
	//-----------------//

	// when running both electron and muon data, pick the right dataset vector and lumi for the MSPlots

	if (!foundMu && !foundEl)
	  datasetsPlot = datasets;
	else if (eventselectedSemiMu) {
	  datasetsPlot = datasetsMu;
	  Luminosity = LuminosityMu;
	}
	else if (eventselectedSemiEl) {
	  datasetsPlot = datasetsEl;
	  Luminosity = LuminosityEl;
	}
	// jet plots for the BTV guys

	string leptonFlav="_other";

	if (eventselectedSemiMu)
	  leptonFlav="_mu";
	else if (eventselectedSemiEl)
	  leptonFlav="_el";

	if (MSPlot.find("Selected_Events_pT_jet1"+leptonFlav) == MSPlot.end()){
	  MSPlot["Selected_Events_pT_jet1"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet1"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	  MSPlot["Selected_Events_pT_jet2"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet2"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	  MSPlot["Selected_Events_pT_jet3"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet3"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	  MSPlot["Selected_Events_pT_jet4"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet4"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	  //MSPlot["Selected_Events_pT_jet4"] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_jet4", 30, 0, 600, "p_{T} (GeV)");
	  MSPlot["Selected_Events_pT_4leadingjets"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_4leadingjets"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	  MSPlot["Selected_Events_pT_alljets"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_pT_alljets"+leptonFlav, 30, 0, 600, "p_{T} (GeV)");
	}

	MSPlot["Selected_Events_pT_jet1"+leptonFlav]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["Selected_Events_pT_jet2"+leptonFlav]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["Selected_Events_pT_jet3"+leptonFlav]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["Selected_Events_pT_jet4"+leptonFlav]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);

	for (unsigned int q=0; q<selectedJets.size(); q++) {

	  MSPlot["Selected_Events_pT_alljets"+leptonFlav]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);

	  if (q<4)
	    MSPlot["Selected_Events_pT_4leadingjets"+leptonFlav]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);

	}

	// other data-mc

	if (MSPlot.find("Selected_Events_nPV3D"+leptonFlav) == MSPlot.end())
	  MSPlot["Selected_Events_nPV3D"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_nPV3D"+leptonFlav, 36, -0.5, 35.5, "nPV");

	MSPlot["Selected_Events_nPV3D"+leptonFlav]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("Selected_Events_nPV"+leptonFlav) == MSPlot.end())
	  MSPlot["Selected_Events_nPV"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_nPV"+leptonFlav, 36, -0.5, 35.5, "nPV");

	MSPlot["Selected_Events_nPV"+leptonFlav]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("Selected_Events_nPV_NONREW"+leptonFlav) == MSPlot.end())
	  MSPlot["Selected_Events_nPV_NONREW"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_nPV_NONREW"+leptonFlav, 36, -0.5, 35.5, "nPV");

	MSPlot["Selected_Events_nPV_NONREW"+leptonFlav]->Fill(vertex.size(), datasets[d], true, Luminosity);
	
	if (MSPlot.find("Selected_Events_JetCombMinChi2"+leptonFlav) == MSPlot.end())
	  MSPlot["Selected_Events_JetCombMinChi2"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_JetCombMinChi2"+leptonFlav, 100, 0, 100, "JetCombMinChi2");
	
	MSPlot["Selected_Events_JetCombMinChi2"+leptonFlav]->Fill(smallestChi2, datasets[d], true, Luminosity*scaleFactor);
	
	if (MSPlot.find("Selected_Events_JetCombMinChi2Exp"+leptonFlav) == MSPlot.end())
	  MSPlot["Selected_Events_JetCombMinChi2Exp"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "Selected_Events_JetCombMinChi2Exp"+leptonFlav, 100, 0, 100, "JetCombMinChi2Exp");

	MSPlot["Selected_Events_JetCombMinChi2Exp"+leptonFlav]->Fill(exp(smallestChi2), datasets[d], true, Luminosity*scaleFactor);

	// MLB
	if (MSPlot.find("BestJetCombMLB"+leptonFlav) == MSPlot.end())
	  MSPlot["BestJetCombMLB"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "BestJetCombMLB"+leptonFlav, 75, 0, 1200, "BestJetCombMLB");

	MSPlot["BestJetCombMLB"+leptonFlav]->Fill((*selectedLepton+*selectedJets[Permutation[3]]).M(), datasets[d], true, Luminosity*scaleFactor);

	// MW
	if (MSPlot.find("BestJetCombMW"+leptonFlav) == MSPlot.end())
	  MSPlot["BestJetCombMW"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "BestJetCombMW"+leptonFlav, 50, 0, 200, "BestJetCombMW");

	MSPlot["BestJetCombMW"+leptonFlav]->Fill((*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]).M(), datasets[d], true, Luminosity*scaleFactor);

	// MTOP
	if (MSPlot.find("BestJetCombMTop"+leptonFlav) == MSPlot.end())
	  MSPlot["BestJetCombMTop"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "BestJetCombMTop"+leptonFlav, 400, 0, 600, "m_{top} (GeV)","Events / 20 GeV");

	MSPlot["BestJetCombMTop"+leptonFlav]->Fill((*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]+*selectedJets[Permutation[2]]).M(), datasets[d], true, Luminosity*scaleFactor);
	
	if (MSPlot.find("BestJetCombMTop_CHI2"+leptonFlav) == MSPlot.end())
	  MSPlot["BestJetCombMTop_CHI2"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "BestJetCombMTop_CHI2"+leptonFlav, 40, 0, 600 ,"m_{top} (GeV)","Events / 20 GeV");
	
	if (smallestChi2 < 90)
	  MSPlot["BestJetCombMTop_CHI2"+leptonFlav]->Fill((*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]+*selectedJets[Permutation[2]]).M(), datasets[d], true, Luminosity*scaleFactor);

	float M3 = -1, maxPt = -1;
        for(int i=0;i<selectedJets.size();i++)
        { 
          for(int j=0;j<i;j++)
          {
            for(int k=0;k<j;k++)
            {
              float combinedPt = (*selectedJets[i]+*selectedJets[j]+*selectedJets[k]).Pt();
              if(combinedPt > maxPt)
              {
                maxPt = combinedPt;
                M3 = (*selectedJets[i]+*selectedJets[j]+*selectedJets[k]).M();
              }
            }
          }
        }
	
	if (MSPlot.find("M3"+leptonFlav) == MSPlot.end()) {
	  MSPlot["M3"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "M3"+leptonFlav,80,0,800,"M3 (GeV)","Nr. of events / 10 GeV");
	  MSPlot["M3_btag"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "M3_btag"+leptonFlav,80,0,800,"M3 (GeV)","Nr. of events / 10 GeV");
	}

	MSPlot["M3"+leptonFlav]->Fill(M3, datasets[d], true, Luminosity*scaleFactor);

	if ((*(selectedJets[Permutation[3]])).btag_trackCountingHighEffBJetTags() > 3.3)
	  MSPlot["M3_btag"+leptonFlav]->Fill(M3, datasets[d], true, Luminosity*scaleFactor);
	
	/*///////////////////////////////////////////
	// NOW WE FILL THE NTUPLE FOR THE BTAGTREE //
	///////////////////////////////////////////*/

	bool lepb_is[4]={false,false,false,false}; //order isb, is not b ishadqq, isradq //first one is decided on flavour
	bool hadb_is[4]={false,false,false,false}; //order isb, is not b ishadqq, isradq //first one is decided on flavour
	bool q1_is[4]={false,false,false,false};
	bool q2_is[4]={false,false,false,false};

	bool R_inAll=false;
	bool R_inHad=false;
	bool R_or_lepb_inHad=false;
	bool allMatched=all4PartonsMatched;

	int indexOflepbJet=-1;    
	int indexOfq1Jet=-1;
	int indexOfq2Jet=-1;
	int indexOfhadbJet=-1;

	if (dataSetName.find("TTbarJets") == 0 && (isSemiMu || isSemiE)) {
	  
	  bool GoodCombinationIsAvailable = false;
	  bool GoodCombinationIsSelected = false;
	  
	  if(all4PartonsMatched){
	    GoodCombinationIsAvailable = true;
	    GoodCombinationIsSelected = true;
	  }
      
	  if(MCPermutation[3] != Permutation[3]) GoodCombinationIsSelected=false;
	  if(MCPermutation[2] != Permutation[2]) GoodCombinationIsSelected=false;

	  //if (MCPermutation[0] != hadronicWJet1_.first) {

	  //cout << "MCPermutation[0]: " << MCPermutation[0] << " hadronicWJet1_.first: " << hadronicWJet1_.first << " hadronicWJet2_.first: " << hadronicWJet2_.first << endl;
	    
	  //}
      
	  if((MCPermutation[2] != Permutation[2] && MCPermutation[1] != Permutation[1]) && (MCPermutation[2] != Permutation[1] && MCPermutation[1] != Permutation[2])) GoodCombinationIsSelected=false;
	  
	  
	  /*if (!GoodCombinationIsAvailable) {

	    for (unsigned int i=0; i<4; i++)
	      
	         cout << "MCPermutation[" << i << "]: " << MCPermutation[i] << " ";

		 }cout << endl;*/
	  
	  for(int i=0;i<4;i++) {
	    if(MCPermutation[i] == Permutation[3]) indexOflepbJet=i;
	  }
 
	  for(int i=0;i<4;i++) {
	    if(MCPermutation[i] == Permutation[2]) indexOfhadbJet=i;
	  }
    
	  for(int i=0;i<4;i++) {
	    if(MCPermutation[i] == Permutation[0]) indexOfq1Jet=i;
	    if(MCPermutation[i] == Permutation[1]) indexOfq2Jet=i;
	  }
	} // end of dataSetName.find("TTbarJets_SemiMu") == 0;

	bool among4=false; //boolean to indicate there is radiation among all 4 selected jets (= 4 highest pt jets)
	bool among3=true; 
	bool among3_lepb=true;

	if(indexOfq1Jet==-1 || indexOfq2Jet==-1 || indexOfhadbJet==-1 || indexOflepbJet==-1){
	  among4=true;
	} else {
	  //cout << indexOfq1Jet << " " << indexOfq2Jet << " " << indexOfhadbJet << " " << indexOflepbJet << endl;
	}
	if(among4) {R_inAll=true;} else {R_inAll=false;}
    

	if(indexOfq1Jet==-1 || indexOfq2Jet==-1 || indexOfhadbJet==-1){
	  among3=false;
	} else {
	  //cout << indexOfq1Jet << " " << indexOfq2Jet << " " << indexOfhadbJet << " " << indexOflepbJet << endl;
	}
	if(among3) {R_inHad=false;} else {R_inHad=true;}
	
	if(indexOfq1Jet==-1 || indexOfq1Jet==3 || indexOfq2Jet==-1 || indexOfq2Jet==3 || indexOfhadbJet==-1 || indexOfhadbJet==3){ 
	  among3_lepb=false;
	} else {
	  //cout << indexOfq1Jet << " " << indexOfq2Jet << " " << indexOfhadbJet << " " << indexOflepbJet << endl;
	}
	
	if(among3_lepb) {R_or_lepb_inHad=false;} else {R_or_lepb_inHad=true;}
	
	//cout << among4 << endl;
	
	for(int i=0;i<4;i++) {
	  if(MCPermutation[i] == Permutation[3]) indexOflepbJet=i;
	  //cout << i << " indexOflepbJet:  " << indexOflepbJet << " MCPermutation[i]: " << MCPermutation[i] << endl;
	}

	//cout << fabs(selectedJets[Permutation[3]]->partonFlavour()) << endl;
		
	if(fabs(selectedJets[Permutation[3]]->partonFlavour())==5){
	  lepb_is[0]=true;
	} else {
	  lepb_is[1]=true;
	  if(dataSetName.find("TTbarJets_SemiMu") == 0){
	    if(indexOflepbJet == 0 || indexOflepbJet == 1)	lepb_is[2]=true;
	    if(indexOflepbJet == -1 ) lepb_is[3]=true;
	  }
	}

	//cout << lepb_is[0] << endl;

	lepb_is[0] ? NTuple->setJet_is_b(true) : NTuple->setJet_is_b(false);
	lepb_is[1] ? NTuple->setJet_is_nonb(true) : NTuple->setJet_is_nonb(false);
	lepb_is[2] ? NTuple->setJet_is_hadqq(true) : NTuple->setJet_is_hadqq(false);
	lepb_is[3] ? NTuple->setJet_is_radq(true) : NTuple->setJet_is_radq(false);
	
	if(fabs(selectedJets[Permutation[2]]->partonFlavour())==5){
	  hadb_is[0]=true;
	} else {
	  hadb_is[1]=true;
	  if(dataSetName.find("TTbarJets_SemiMu") == 0){
	    if(indexOfhadbJet == 0 || indexOfhadbJet == 1)	hadb_is[2]=true;
	    if(indexOfhadbJet == -1 ) hadb_is[3]=true;
	  }
	}
	
	hadb_is[0] ? NTuple->setJethadb_is_b(true) : NTuple->setJethadb_is_b(false);
	hadb_is[1] ? NTuple->setJethadb_is_nonb(true) : NTuple->setJethadb_is_nonb(false);
	hadb_is[2] ? NTuple->setJethadb_is_hadqq(true) : NTuple->setJethadb_is_hadqq(false);
	hadb_is[3] ? NTuple->setJethadb_is_radq(true) : NTuple->setJethadb_is_radq(false);
	
	if(fabs(selectedJets[Permutation[0]]->partonFlavour())==5){
	  q1_is[0]=true;
	} else {	
	  q1_is[1]=true;
	  if(dataSetName.find("TTbarJets_SemiMu") == 0){
	    if(indexOfq1Jet == 0 || indexOfq1Jet == 1) q1_is[2]=true;
	    if(indexOfq1Jet == -1 ) q1_is[3]=true;
	  }
	}
	
	if(fabs(selectedJets[Permutation[1]]->partonFlavour())==5){
	  q2_is[0]=true;
	} else {	 
	  q2_is[1]=true;
	  if(dataSetName.find("TTbarJets_SemiMu") == 0){
	    if(indexOfq2Jet == 0 || indexOfq2Jet == 1) q2_is[2]=true;
	    if(indexOfq2Jet == -1 ) q2_is[3]=true;
	  }
	}
	
	q1_is[3] ? NTuple->setJetControl_is_radq(true) : NTuple->setJetControl_is_radq(false);
	q1_is[2] ? NTuple->setJetControl_is_hadqq(true) : NTuple->setJetControl_is_hadqq(false); 
	q1_is[0] ? NTuple->setJetControl_is_b(true) : NTuple->setJetControl_is_b(false);
	q1_is[1] ? NTuple->setJetControl_is_nonb(true) : NTuple->setJetControl_is_nonb(false);
	
	q2_is[3] ? NTuple->setJetControl2_is_radq(true) : NTuple->setJetControl2_is_radq(false);
	q2_is[2] ? NTuple->setJetControl2_is_hadqq(true) : NTuple->setJetControl2_is_hadqq(false); 
	q2_is[0] ? NTuple->setJetControl2_is_b(true) : NTuple->setJetControl2_is_b(false);
	q2_is[1] ? NTuple->setJetControl2_is_nonb(true) : NTuple->setJetControl2_is_nonb(false);
	
	bool doRemoveMuInJets=false;
	bool doSkipMuInJets=false;
	int minId=-1;
    
	if(doRemoveMuInJets){
	  cout<<"removing muon in jets" << endl;
	  
	  //Check if there is a muon inside the lepb jet
	  //final part of the event selection
	  double deltaRjetmuons=999;
	  for(int iMuO=0; iMuO<vetoMuons.size(); iMuO++){
	    double deltaRjetmuonstmp = selectedJets[Permutation[3]]->DeltaR(*vetoMuons[iMuO]);
	    
	    //	  cout << deltaRjetmuonstmp << endl;
	    
	    if(deltaRjetmuonstmp<deltaRjetmuons) {
	      deltaRjetmuons=deltaRjetmuonstmp;
	      minId=iMuO;
	    }
	  }
	  
	  //	if(minId>-1) {cout << deltaRjetmuons << " " <<(*selectedLooseMuons[minId]).Pt() << endl;}
	  NTuple->setJetMuonOverlap(deltaRjetmuons);
	  
	  //if(deltaRjetmuons<0.5) doSkipMuInJets=true;
	  if(deltaRjetmuons>0.5) doSkipMuInJets=true;
	}

	///if(doSkipMuInJets) continue;

	NTuple->setEventId(event->eventId());
	NTuple->setRunId(event->runId());
	NTuple->setLumiBlockId(event->lumiBlockId());

	NTuple->setDecay(eventselectedSemiMu,isSemiMu,eventselectedSemiEl,isSemiE);

	NTuple->setnPV(vertex.size());

	NTuple->setMlj((*(selectedJets[Permutation[3]])+*selectedLepton).M());

	/*if ((*(selectedJets[Permutation[3]])+*selectedLepton).M() < 30) {

	  cout << endl << "event ievt" << endl;

	  cout << "Muon: (" << selectedLepton->E() << "," << selectedLepton->Px() << "," << selectedLepton->Py() << "," << selectedLepton->Pz() << ") " << "pT = " << selectedLepton->Pt() << "," << endl;
	  cout << "bCand: (" << selectedJets[Permutation[3]]->E() << "," << selectedJets[Permutation[3]]->Px() << "," << selectedJets[Permutation[3]]->Py() << "," << selectedJets[Permutation[3]]->Pz() << ") " << "pT = " << selectedJets[Permutation[3]]->Pt() << "," << endl;

	  cout << "mlj = " << (*(selectedJets[Permutation[3]])+*selectedLepton).M() << endl;

	  //cout << "bCand pflav = " << (*(selectedJets[Permutation[3]])).partonFlavour() << endl;

	  cout << "dR (mu,bcand) = " << (*selectedJets[Permutation[3]]).DeltaR(*selectedLepton) << endl;

	  }*/

	NTuple->setM3(M3);

	NTuple->setBtag(0,(*(selectedJets[Permutation[3]])).btag_trackCountingHighEffBJetTags());
	NTuple->setBtag(1,(*(selectedJets[Permutation[3]])).btag_trackCountingHighPurBJetTags());

	if (MSPlot.find("bcand_trackCountingHighEffBJetTags"+leptonFlav) == MSPlot.end())
	  MSPlot["bcand_trackCountingHighEffBJetTags"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "bcand_trackCountingHighEffBJetTags"+leptonFlav, 50, -10, 30, "Btag Discr");
	MSPlot["bcand_trackCountingHighEffBJetTags"+leptonFlav]->Fill((*(selectedJets[Permutation[3]])).btag_trackCountingHighEffBJetTags(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("bcand_trackCountingHighPurBJetTags"+leptonFlav) == MSPlot.end())
	  MSPlot["bcand_trackCountingHighPurBJetTags"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "bcand_trackCountingHighPurBJetTags"+leptonFlav, 50, -10, 30, "Btag Discr");
	MSPlot["bcand_trackCountingHighPurBJetTags"+leptonFlav]->Fill((*(selectedJets[Permutation[3]])).btag_trackCountingHighPurBJetTags(), datasets[d], true, Luminosity*scaleFactor);


	NTuple->setBtag(2,(*(selectedJets[Permutation[3]])).btag_simpleSecondaryVertexHighEffBJetTags());
	NTuple->setBtag(3,(*(selectedJets[Permutation[3]])).btag_simpleSecondaryVertexHighPurBJetTags());

	if (MSPlot.find("bcand_simpleSecondaryVertexHighEffBJetTags"+leptonFlav) == MSPlot.end())
	  MSPlot["bcand_simpleSecondaryVertexHighEffBJetTags"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "bcand_simpleSecondaryVertexHighEffBJetTags"+leptonFlav, 50, 0, 8, "Btag Discr");
	MSPlot["bcand_simpleSecondaryVertexHighEffBJetTags"+leptonFlav]->Fill((*(selectedJets[Permutation[3]])).btag_simpleSecondaryVertexHighEffBJetTags(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("bcand_simpleSecondaryVertexHighPurBJetTags"+leptonFlav) == MSPlot.end())
	  MSPlot["bcand_simpleSecondaryVertexHighPurBJetTags"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "bcand_simpleSecondaryVertexHighPurBJetTags"+leptonFlav, 50, 0, 8, "Btag Discr");
	MSPlot["bcand_simpleSecondaryVertexHighPurBJetTags"+leptonFlav]->Fill((*(selectedJets[Permutation[3]])).btag_simpleSecondaryVertexHighPurBJetTags(), datasets[d], true, Luminosity*scaleFactor);
	
	NTuple->setBtag(4,(*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexBJetTags());
	NTuple->setBtag(5,(*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexMVABJetTags());

	if (MSPlot.find("bcand_combinedSecondaryVertexBJetTags"+leptonFlav) == MSPlot.end())
	  MSPlot["bcand_combinedSecondaryVertexBJetTags"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "bcand_combinedSecondaryVertexBJetTags"+leptonFlav, 50, -1, 2, "Btag Discr");
	MSPlot["bcand_combinedSecondaryVertexBJetTags"+leptonFlav]->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexBJetTags(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("bcand_combinedSecondaryVertexMVABJetTags"+leptonFlav) == MSPlot.end())
	  MSPlot["bcand_combinedSecondaryVertexMVABJetTags"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "bcand_combinedSecondaryVertexMVABJetTags"+leptonFlav, 50, -1, 2, "Btag Discr");
	MSPlot["bcand_combinedSecondaryVertexMVABJetTags"+leptonFlav]->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexMVABJetTags(), datasets[d], true, Luminosity*scaleFactor);

	NTuple->setBtag(6,(*(selectedJets[Permutation[3]])).btag_jetBProbabilityBJetTags());
	NTuple->setBtag(7,(*(selectedJets[Permutation[3]])).btag_jetProbabilityBJetTags());

	if (MSPlot.find("bcand_jetBProbabilityBJetTags"+leptonFlav) == MSPlot.end())
	  MSPlot["bcand_jetBProbabilityBJetTags"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "bcand_jetBProbabilityBJetTags"+leptonFlav, 50, 0, 8, "Btag Discr");
	MSPlot["bcand_jetBProbabilityBJetTags"+leptonFlav]->Fill((*(selectedJets[Permutation[3]])).btag_jetBProbabilityBJetTags(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("bcand_jetProbabilityBJetTags"+leptonFlav) == MSPlot.end())
	  MSPlot["bcand_jetProbabilityBJetTags"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "bcand_jetProbabilityBJetTags"+leptonFlav, 50, 0, 3, "Btag Discr");
	MSPlot["bcand_jetProbabilityBJetTags"+leptonFlav]->Fill((*(selectedJets[Permutation[3]])).btag_jetProbabilityBJetTags(), datasets[d], true, Luminosity*scaleFactor);

	NTuple->setBtag(8,(*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexRetrainedBJetTags());

        NTuple->setBtag(9,(*(selectedJets[Permutation[3]])).btag_softPFMuonRetrainedBJetsTags());
 	NTuple->setBtag(10,(*(selectedJets[Permutation[3]])).btag_softPFElectronRetrainedBJetsTags());
        NTuple->setBtag(11,(*(selectedJets[Permutation[3]])).btag_combinedCSVSLBJetTags());
	NTuple->setBtag(12,(*(selectedJets[Permutation[3]])).btag_combinedCSVJPBJetTags());
	NTuple->setBtag(13,(*(selectedJets[Permutation[3]])).btag_combinedCSVJPSLBJetTags());
	
	// CS jet 1
	NTuple->setBtagCS1(0,(*(selectedJets[Permutation[0]])).btag_trackCountingHighEffBJetTags());
	NTuple->setBtagCS1(1,(*(selectedJets[Permutation[0]])).btag_trackCountingHighPurBJetTags());
	NTuple->setBtagCS1(2,(*(selectedJets[Permutation[0]])).btag_simpleSecondaryVertexHighEffBJetTags());
	NTuple->setBtagCS1(3,(*(selectedJets[Permutation[0]])).btag_simpleSecondaryVertexHighPurBJetTags());
	NTuple->setBtagCS1(4,(*(selectedJets[Permutation[0]])).btag_combinedSecondaryVertexBJetTags());
	NTuple->setBtagCS1(5,(*(selectedJets[Permutation[0]])).btag_combinedSecondaryVertexMVABJetTags());
	NTuple->setBtagCS1(6,(*(selectedJets[Permutation[0]])).btag_jetBProbabilityBJetTags());
	NTuple->setBtagCS1(7,(*(selectedJets[Permutation[0]])).btag_jetProbabilityBJetTags());
	NTuple->setBtagCS1(8,(*(selectedJets[Permutation[0]])).btag_combinedSecondaryVertexRetrainedBJetTags());
        NTuple->setBtagCS1(9,(*(selectedJets[Permutation[0]])).btag_softPFMuonRetrainedBJetsTags());
 	NTuple->setBtagCS1(10,(*(selectedJets[Permutation[0]])).btag_softPFElectronRetrainedBJetsTags());
        NTuple->setBtagCS1(11,(*(selectedJets[Permutation[0]])).btag_combinedCSVSLBJetTags());
	NTuple->setBtagCS1(12,(*(selectedJets[Permutation[0]])).btag_combinedCSVJPBJetTags());
	NTuple->setBtagCS1(13,(*(selectedJets[Permutation[0]])).btag_combinedCSVJPSLBJetTags());

	// CS jet 2
	NTuple->setBtagCS2(0,(*(selectedJets[Permutation[1]])).btag_trackCountingHighEffBJetTags());
	NTuple->setBtagCS2(1,(*(selectedJets[Permutation[1]])).btag_trackCountingHighPurBJetTags());
	NTuple->setBtagCS2(2,(*(selectedJets[Permutation[1]])).btag_simpleSecondaryVertexHighEffBJetTags());
	NTuple->setBtagCS2(3,(*(selectedJets[Permutation[1]])).btag_simpleSecondaryVertexHighPurBJetTags());
	NTuple->setBtagCS2(4,(*(selectedJets[Permutation[1]])).btag_combinedSecondaryVertexBJetTags());
	NTuple->setBtagCS2(5,(*(selectedJets[Permutation[1]])).btag_combinedSecondaryVertexMVABJetTags());
	NTuple->setBtagCS2(6,(*(selectedJets[Permutation[1]])).btag_jetBProbabilityBJetTags());
	NTuple->setBtagCS2(7,(*(selectedJets[Permutation[1]])).btag_jetProbabilityBJetTags());
	NTuple->setBtagCS2(8,(*(selectedJets[Permutation[1]])).btag_combinedSecondaryVertexRetrainedBJetTags());
        NTuple->setBtagCS2(9,(*(selectedJets[Permutation[1]])).btag_softPFMuonRetrainedBJetsTags());
 	NTuple->setBtagCS2(10,(*(selectedJets[Permutation[1]])).btag_softPFElectronRetrainedBJetsTags());
        NTuple->setBtagCS2(11,(*(selectedJets[Permutation[1]])).btag_combinedCSVSLBJetTags());
	NTuple->setBtagCS2(12,(*(selectedJets[Permutation[1]])).btag_combinedCSVJPBJetTags());
	NTuple->setBtagCS2(13,(*(selectedJets[Permutation[1]])).btag_combinedCSVJPSLBJetTags());

	if (MSPlot.find("bcand_combinedSecondaryVertexRetrainedBJetTags"+leptonFlav) == MSPlot.end())
	  MSPlot["bcand_combinedSecondaryVertexRetrainedBJetTags"+leptonFlav] = new MultiSamplePlot(datasetsPlot, "bcand_combinedSecondaryVertexRetrainedBJetTags"+leptonFlav, 50, -1, 2, "Btag Discr");
	MSPlot["bcand_combinedSecondaryVertexRetrainedBJetTags"+leptonFlav]->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexRetrainedBJetTags(), datasets[d], true, Luminosity*scaleFactor);

	
	//vector <TRootCaloJet*> selectedCaloJets = convertToCaloJets(selectedJets);
	
	NTuple->setPt((*(selectedJets[Permutation[3]])).Pt());
	NTuple->setPthadb((*(selectedJets[Permutation[2]])).Pt());
	/*cout << (*(selectedJets[Permutation[2]])).Pt() << " " << (*(selectedJets[Permutation[3]])).Pt() << endl;
	cout << NTuple->pthadb() << " " << NTuple->pt() << endl;
	cout << endl;
	*/
	NTuple->setE((*(selectedJets[Permutation[3]])).E());
	NTuple->setEta((*(selectedJets[Permutation[3]])).Eta());//commented for deltaR
	NTuple->setWeight(1/(datasets[d]->EquivalentLumi()));

	NTuple->setScaleFactor3D(scaleFactor);
	NTuple->setScaleFactor(scaleFactor);

	// For PDF Unc

	NTuple->setx1(event->xParton1());
	NTuple->setx2(event->xParton2());
	NTuple->setq2(event->factorizationScale());
	NTuple->setid1(event->idParton1());
	NTuple->setid2(event->idParton2());

	//cout << dataSetName << " Lumi: " <<  datasets[d]->EquivalentLumi() << endl;
	NTuple->setDataSetNumber(d);
	NTuple->setDataSetName(dataSetName);
	NTuple->setPartonFlavour(selectedJets[Permutation[3]]->partonFlavour());
	/*NTuple->setNConstituents(selectedCaloJets[Permutation[3]]->nConstituents());
	NTuple->setN90(selectedCaloJets[Permutation[3]]->n90());
	NTuple->setN60(selectedCaloJets[Permutation[3]]->n60());
	NTuple->setJetArea(selectedCaloJets[Permutation[3]]->jetArea());
	NTuple->setPileupEnergy(selectedCaloJets[Permutation[3]]->pileupEnergy());
	NTuple->setMaxDistance(selectedCaloJets[Permutation[3]]->maxDistance());
	NTuple->setMaxEInEmTowers(selectedCaloJets[Permutation[3]]->maxEInEmTowers());
	NTuple->setMaxEInHadTowers(selectedCaloJets[Permutation[3]]->maxEInHadTowers());
	NTuple->setTowersArea(selectedCaloJets[Permutation[3]]->towersArea());
	NTuple->setEcalEnergyFraction(selectedCaloJets[Permutation[3]]->ecalEnergyFraction());
	NTuple->setHcalEnergyFraction(selectedCaloJets[Permutation[3]]->hcalEnergyFraction());*/  
	NTuple->setChiSq(smallestChi2);  
	NTuple->setPtMuon((*selectedLepton).Pt());   
	NTuple->setEMuon((*selectedLepton).E());   
	NTuple->setEtaMuon((*selectedLepton).Eta());

	if (eventselectedSemiMu)
	  NTuple->setChargeMuon(selectedMuons[0]->charge());
	else if (eventselectedSemiEl)
	  NTuple->setChargeMuon(selectedElectrons[0]->charge());

	if (eventselectedSemiEl) {
	  //cout << selectedElectrons[0]->mvaTrigId() << endl;
	  NTuple->setmvaTrigID(selectedElectrons[0]->mvaTrigId());
	} else
	  NTuple->setmvaTrigID(-1);

	NTuple->setLepBlCandMass( (*selectedJets[Permutation[3]]+*selectedLepton).M() );
	NTuple->setHadTopCandMass( (*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]+*selectedJets[Permutation[2]]).M() );
	NTuple->setHadWCandMass( (*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]).M() );
	NTuple->setR_inAll(R_inAll);
	NTuple->setR_inHad(R_inHad);
	NTuple->setR_or_lepb_inHad(R_or_lepb_inHad);
	NTuple->setAllMatched(allMatched);


	//Setting some extra variables which could help reducing the amount of background:
	
	NTuple->setnJets(selectedJets.size());
	NTuple->setMET(mets[0]->Pt());
	NTuple->setBtag_trackCountingHighEffBJetTags_hadb((*(selectedJets[Permutation[2]])).btag_trackCountingHighEffBJetTags());
	
	
	TLorentzVector muon0, jet0, jet1, jet2, jet3, met0;
	TRootMET* met = (TRootMET*) mets[0]; 
	met0.SetPxPyPzE(met->Px(),met->Py(),met->Pz(),met->Energy());
	muon0.SetPxPyPzE(selectedLepton->Px(),selectedLepton->Py(),selectedLepton->Pz(),selectedLepton->Energy());
	jet0.SetPxPyPzE(selectedJets[Permutation[0]]->Px(),selectedJets[Permutation[0]]->Py(),selectedJets[Permutation[0]]->Pz(),selectedJets[Permutation[0]]->Energy());
	jet1.SetPxPyPzE(selectedJets[Permutation[1]]->Px(),selectedJets[Permutation[1]]->Py(),selectedJets[Permutation[1]]->Pz(),selectedJets[Permutation[1]]->Energy());
	jet2.SetPxPyPzE(selectedJets[Permutation[2]]->Px(),selectedJets[Permutation[2]]->Py(),selectedJets[Permutation[2]]->Pz(),selectedJets[Permutation[2]]->Energy());
	jet3.SetPxPyPzE(selectedJets[Permutation[3]]->Px(),selectedJets[Permutation[3]]->Py(),selectedJets[Permutation[3]]->Pz(),selectedJets[Permutation[3]]->Energy());
	
	double deltaR3=0;
	double deltaR0=0;
	double deltaR1=0;
	deltaR0 = ROOT::Math::VectorUtil::DeltaR(muon0,jet0);
	deltaR1 = ROOT::Math::VectorUtil::DeltaR(muon0,jet1);
	deltaR3 = ROOT::Math::VectorUtil::DeltaR(muon0,jet3);
	NTuple->setDelRlj(deltaR3);
	
	double delOmega3=0;
	delOmega3 = ROOT::Math::VectorUtil::Angle(muon0,jet3);
	NTuple->setDelOmegalj(delOmega3);
	
	double delOmega2=0;
	delOmega2 = ROOT::Math::VectorUtil::Angle(muon0,jet2);
	NTuple->setDelOmegaljhad(delOmega2);
	
	double delOmegaControl=0;
	delOmegaControl = ROOT::Math::VectorUtil::Angle(muon0,jet0);
	NTuple->setDelOmegaljControl(delOmegaControl);
	
	double delOmegaControl2=0;
	delOmegaControl2 = ROOT::Math::VectorUtil::Angle(muon0,jet1);
	NTuple->setDelOmegaljControl2(delOmegaControl2);
	
	NTuple->setPthadtop((*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]+*selectedJets[Permutation[2]]).M());
	
	double delPhitt=0;
	delPhitt = ROOT::Math::VectorUtil::DeltaPhi((muon0+met0+jet3),(jet0+jet1+jet2));
	NTuple->setDelPhitt(delPhitt);

	//Here i choose a fixed way to fill mlbcontrol and mlbcontrol2 depending on the pt of the jet
	//the first entry is the jet with highest pt, the second has lowes pt

	//if there is a fifth jet with pt higher than the 30 GeV cut also write his info
	if(selectedJets.size()>4){
	  if((*selectedJets[4]).btag_trackCountingHighEffBJetTags()<btagcut){ 
	    //if(fabs((*(ChiSqMatchSelectedJet[0])).partonFlavour())!=5){ 
	    NTuple->setMljControl3((*selectedJets[0]+*selectedLepton).M());
	    NTuple->setPtControl3((*selectedJets[0]).Pt());
	    NTuple->setEtaControl3((*selectedJets[0]).Eta());
	  }
	  
	  NTuple->setEtafifth((*selectedJets[0]).Eta());
	  NTuple->setPtfifth((*selectedJets[0]).Pt());

	}

	//instead of using eta, use the angle between the jet and the muon -> I'm abusing the eta and etacontrol, etacontrol2 to fill this (possible complication when doing the pt-eta binning, but I will not do eta binning right now)
	
	if((*(selectedJets[Permutation[0]])).Pt()>(*(selectedJets[Permutation[1]])).Pt()){
	  NTuple->setbTagControl((*(selectedJets[Permutation[0]])).btag_trackCountingHighEffBJetTags());
	  
	  if((*(selectedJets[Permutation[0]])).btag_trackCountingHighEffBJetTags()<btagcut){ 
	    //if(fabs((*(selectedJets[Permutation[0]])).partonFlavour())!=5){ 
	    NTuple->setMljControl((*(selectedJets[Permutation[0]])+*selectedLepton).M());
	    NTuple->setPtControl((*(selectedJets[Permutation[0]])).Pt());
	    NTuple->setEControl((*(selectedJets[Permutation[0]])).E());
	    NTuple->setEtaControl((*(selectedJets[Permutation[0]])).Eta());
	    NTuple->setPartonFlavourControl((*(selectedJets[Permutation[0]])).partonFlavour());
	    NTuple->setDelRljControl(deltaR0);
	    
	  } else {
	    NTuple->setMljControl(-1);
	    NTuple->setPtControl(-1);
	    NTuple->setEControl(-1);
	    NTuple->setEtaControl(999);
	  }
	  
	  NTuple->setbTagControl2((*(selectedJets[Permutation[1]])).btag_trackCountingHighEffBJetTags());
	  if((*(selectedJets[Permutation[1]])).btag_trackCountingHighEffBJetTags()<btagcut){ 
	    //if(fabs((*(selectedJets[Permutation[1]])).partonFlavour())!=5){ 
	    NTuple->setMljControl2((*(selectedJets[Permutation[1]])+*selectedLepton).M());
	    NTuple->setPtControl2((*(selectedJets[Permutation[1]])).Pt());
	    NTuple->setEControl2((*(selectedJets[Permutation[1]])).E());
	    NTuple->setEtaControl2((*(selectedJets[Permutation[1]])).Eta());
	    NTuple->setPartonFlavourControl2((*(selectedJets[Permutation[1]])).partonFlavour());
	    NTuple->setDelRljControl2(deltaR1);
	    
	  } else {
	    NTuple->setMljControl2(-1);
	    NTuple->setPtControl2(-1);
	    NTuple->setEControl2(-1);
	    NTuple->setEtaControl2(999);
	  }
	} else {
	  NTuple->setbTagControl((*(selectedJets[Permutation[1]])).btag_trackCountingHighEffBJetTags());
	  if((*(selectedJets[Permutation[1]])).btag_trackCountingHighEffBJetTags()<btagcut){ 
	    //if(fabs((*(selectedJets[Permutation[1]])).partonFlavour())!=5){ 
	    NTuple->setMljControl((*(selectedJets[Permutation[1]])+*selectedLepton).M());
	    NTuple->setPtControl((*(selectedJets[Permutation[1]])).Pt());
	    NTuple->setEControl((*(selectedJets[Permutation[1]])).E());
	    NTuple->setEtaControl((*(selectedJets[Permutation[1]])).Eta());
	    NTuple->setPartonFlavourControl((*(selectedJets[Permutation[1]])).partonFlavour());
	    NTuple->setDelRljControl(deltaR1);
	    
	  } else {
	    NTuple->setMljControl(-1);
	    NTuple->setPtControl(-1);
	    NTuple->setEControl(-1);
	    NTuple->setEtaControl(999);
	  }
	  //NTuple->setPartonFlavourControl(selectedJets[Permutation[1]]->partonFlavour());
	  //NTuple->setbTagControl(selectedJets[Permutation[1]]->btag_trackCountingHighEffBJetTags());
	  
	  NTuple->setbTagControl2((*(selectedJets[Permutation[0]])).btag_trackCountingHighEffBJetTags());
	  
	  if((*(selectedJets[Permutation[0]])).btag_trackCountingHighEffBJetTags()<btagcut){ 
	    //if(fabs((*(selectedJets[Permutation[0]])).partonFlavour())!=5){ 
	    NTuple->setMljControl2((*(selectedJets[Permutation[0]])+*selectedLepton).M());
	    NTuple->setPtControl2((*(selectedJets[Permutation[0]])).Pt());
	    NTuple->setEControl2((*(selectedJets[Permutation[0]])).E());
	    NTuple->setEtaControl2((*(selectedJets[Permutation[0]])).Eta());
	    NTuple->setPartonFlavourControl2((*(selectedJets[Permutation[0]])).partonFlavour());
	    NTuple->setDelRljControl2(deltaR0);
	    	    
	  } else {
	    NTuple->setMljControl2(-1);
	    NTuple->setPtControl2(-1);
	    NTuple->setEControl2(-1);
	    NTuple->setEtaControl2(999);
	  }
	}
	
	// Fill the TTree
	BtagTTree->Fill();
      }

      ///////////////////////////////////////
      // END OF EVENT REMOVING SOME STUFF //
      //////////////////////////////////////

    }			//loop on events

    cout<<endl;

    cout << "+> " << nSelectedMu << " mu+jets events where selected"<< endl;
    cout << "+> " << nSelectedEl << " e+jets events where selected"<< endl;

    cout << "+> mu+jets ttbar sel eff " << ((double)nSelectedMu/(double)datasets[d]->NofEvtsToRunOver())*100 << endl;

    if (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) {
      string fname = dataSetName+"_QCDEst_RunRangeCounts.txt";
      
      cout << endl << "DataDriven QCD summary: # events in each runrange (Saving results in " << fname << ")" << endl;
      
      fstream f; f.open(fname.c_str(), ios::out | ios::trunc);
      for (unsigned int r=0; r<QCDEstcounts.size(); r++) {
	cout << "++ RunRange " << r << ": " << QCDEstcounts[r] << " events." << endl;
	f << QCDEstcounts[r] << " ";
      }
      f.close();
      cout << endl;
    }

    

    ////////////////////////////////////////////////
    // STORE INPUT MASSES FOR CHI2 JETCOMBINATION //
    ////////////////////////////////////////////////

    if (!useMassesAndResolutions && dataSetName.find("TTbarJets") == 0) {
	
      string filename = "BtagMassPlots";

      //if (argc > 2)
      //filename += "_"+string(argv[2]);

      filename += ".root";
      
      cout << "INFO: Creating output file to store the masses for the Chi2 jetcombiner: " << filename << endl;    

      TFile* outfile = new TFile(filename.c_str(),"RECREATE");
      
      outfile->cd();
      
      histo1D["hadronicRecoWMass"]->Write();
      histo1D["hadronicRecoTopMass"]->Write();

      histo1D["hadronicRecoWMass-REW"]->Write();
      histo1D["hadronicRecoTopMass-REW"]->Write();
      
      histo1D["hadronicRecoWMass-REW3D"]->Write();
      histo1D["hadronicRecoTopMass-REW3D"]->Write();

      // fit the distributions
      
      for (unsigned int f=0; f<6;f++) {

	TH1F* histo;

	//if (f==0) histo=histo1D["hadronicRecoWMass"]; if (f==1) histo=histo1D["hadronicRecoTopMass"];
	//if (f==2) histo=histo1D["hadronicRecoWMass-REW"]; if (f==3) histo=histo1D["hadronicRecoTopMass-REW"];

	if (f==0) histo=histo1D["hadronicRecoWMass"]; if (f==1) histo=histo1D["hadronicRecoTopMass"];
	if (f==2) histo=histo1D["hadronicRecoWMass-REW"]; if (f==3) histo=histo1D["hadronicRecoTopMass-REW"];

	if (f==4) histo=histo1D["hadronicRecoWMass-REW3D"]; if (f==5) histo=histo1D["hadronicRecoTopMass-REW3D"];
	TF1 *fitfunc;
	
	string func_title = string(histo->GetName())+"_Fitted";
	
	double rms = histo->GetRMS();
	double maxbin =  histo->GetBinCenter(histo->GetMaximumBin());

	fitfunc = new TF1(func_title.c_str(),"gaus");
	fitfunc->SetRange(maxbin-rms,maxbin+rms);
	histo->Fit(fitfunc,"RQ");

	fitfunc->Write();

	delete fitfunc;
      
      }
            
      outfile->Close();
      
      useMassesAndResolutions=true;
	
      //exit(1);
      // go back to the ttjets to run again

      d--;

      // remove all MSPlots

      for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++) {

	delete it->second;

      } MSPlot.clear();

    } else {

      ///////////////////////////////////////////////////////////////////////
      // FILLING THE BtagTree because we have the ""GOOD"" jet combination //
      ///////////////////////////////////////////////////////////////////////
      
      // open rootfile to store the tree

      BTreeFile->cd();

      BtagTTree->Write();

      // writing histos 

      TDirectory* th1dir = BTreeFile->mkdir("ControlPlots");

      for(std::map<std::string,TH1F*>::const_iterator it = histo1D_btag.begin(); it != histo1D_btag.end(); it++) {
	  TH1F *temp = it->second;
	  int N = temp->GetNbinsX();
	  temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
	  temp->SetBinContent(N+1,0);
	  temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries..
	  th1dir->cd();
	  temp->Write();

	  delete histo1D_btag[it->first];
      }

      for(std::map<std::string,TH2F*>::const_iterator it = histo2D_btag.begin(); it != histo2D_btag.end(); it++) {
	  th1dir->cd();
	  it->second->Write();

	  delete histo2D_btag[it->first];
      }

      TDirectory* th1dirs = BTreeFile->mkdir("SystPlots");

      for(std::map<std::string,TH1D*>::const_iterator it = histo1D_syst.begin(); it != histo1D_syst.end(); it++) {
	  TH1D *temp = it->second;
	  int N = temp->GetNbinsX();
	  //temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
	  //temp->SetBinContent(N+1,0);
	  //temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries..
	  th1dirs->cd();
	  temp->Write();

	  delete histo1D_syst[it->first];
      }

      for (std::map<std::string,TGraphErrors*>::const_iterator it = graphErr_btag.begin(); it != graphErr_btag.end(); ++it) {
	th1dir->cd();
	it->second->Write();
	
	delete graphErr_btag[it->first];
      }

    }

    //////////////
    // CLEANING //
    //////////////

    if (jecUnc) delete jecUnc;
    if (jetTools) delete jetTools;

    if (NTuple) delete NTuple;
    if (BtagTTree) delete BtagTTree;
    if (BTreeFile) {
      BTreeFile->Close();
      delete BTreeFile;
    }
  
    //important: free memory
    treeLoader.UnLoadDataset();
    
    }				//loop on datasets

    //Once everything is filled ...
    if (verbose > 0)
      cout << " We ran over all the data ;-)" << endl;
    cout << "+> Chi2 performance: nAll = " << nAll << " nGoodCombi = " << nGoodCombi << " nGoodChi2 " << nGoodChi2 << endl;
   
    cout << endl;
    cout << "QCD Yield for MC in e+jets channel for L=" << LuminosityEl << endl;
    cout << "# QCD = " <<  yieldQCD_MC_SemiEl << endl;
    cout << endl;

  
    //Selection tables
    selecTableSemiEl.TableCalculator(false, true, true, true, true);
    string selectiontableEl = "SelectionTable"+postfix+"_BTAG_SEMIEL.tex";
    selecTableSemiEl.Write(selectiontableEl.c_str());
    
    if (runSpecificSample) return 0;
    
    // Do some special things with certain plots (normalize, BayesDivide, ... )
    if (verbose > 0)
      cout << "Treating the special plots." << endl;
    
    ///////////////////
    // Writing
    //////////////////
    if (verbose > 1)
      cout << " - Writing outputs on files ..." << endl;
    
    string pathPNG = "PlotsBTAG"+postfix;
    /*if (argc >= 3){
      string sample=string(argv[2]);
      pathPNG = pathPNG+"_"+sample;
      }*/
    
    pathPNG = pathPNG +"/"; 	
    
    
    
    mkdir(pathPNG.c_str(),0777);
    mkdir((pathPNG+"MSPlot/").c_str(),0777);
    
    // 1D 
    TDirectory* th1dir = fout->mkdir("1D_histograms");
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
      {
	MultiSamplePlot *temp = it->second;
	string name = it->first;
	temp->Draw(false, name, true, true, true, true, true, 1, false);

	temp->Write(fout, name, true, pathPNG+"MSPlot/");
      }
    
    //Write histograms
    fout->cd();
    th1dir->cd();
    
    fout->cd();
    
    for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
      {
	TH1F *temp = it->second;
	int N = temp->GetNbinsX();
  	temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
  	temp->SetBinContent(N+1,0);
	temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
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
    
    fout->cd();
    //add configuration info
    fout->cd();
    configTree->Fill();
    configTree->Write();
    
    //
    if (verbose > 1)
      cout << " - Done with writing the module outputs in the ouput file ..." << endl;
    cout << " - Closing the output file now..." << endl;
    //  fout->Write();
    fout->Close();
    
    delete fout;
    delete tcdatasets;
    delete tcAnaEnv;
    delete configTree;
    
    cout << "nTTbar " << nTTbar << " nTTbar_w " << nTTbar_w << endl; 
    
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "           hasn't crashed yet ;-)           " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}
