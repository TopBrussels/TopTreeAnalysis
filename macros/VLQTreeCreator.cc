///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Based on AnalysisSkeleton.cc: This macro is intended to be an example analysis macro which works out of the box. /////
/////       It should serve as the first port of call for new users of the TopBrussels framework.              /////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysis/VLQSearch/interface/VLQTree.h"

#include "Style.C"

using namespace std;
using namespace reweight;
using namespace TopTree;

int main (int argc, char *argv[])
{	
	/////////////////////////////////////
  // SWITCHES FOR SYSTAMATIC SAMPLES //
  /////////////////////////////////////
  
  int systematic = 0; // 0: off 1: JES minus, 2: JES plus, ...
  
  if (argc > 1) // systematic shift is command line parameter 1
    if (atoi(argv[1]) == 0 || atoi(argv[1]) == 1 || atoi(argv[1]) == 2 || atoi(argv[1]) == 3 || atoi(argv[1]) == 4  || atoi(argv[1]) == 5 || atoi(argv[1]) == 6  || atoi(argv[1]) == 7 || atoi(argv[1]) == 8 || atoi(argv[1]) == 9 || atoi(argv[1]) == 10 || atoi(argv[1]) == 11 || atoi(argv[1]) == 12 || atoi(argv[1]) == 13 || atoi(argv[1]) == 14 || atoi(argv[1]) == 15 || atoi(argv[1]) == 16 || atoi(argv[1]) == 17 || atoi(argv[1]) == 18 || atoi(argv[1]) == -1) // -1 == do the invertediso setup
      systematic=atoi(argv[1]);
  
  bool runSpecificSample=false;
  
  int nSpecSample=-19;

  if (argc > 2) {

    runSpecificSample=true;

    nSpecSample=atoi(argv[2]);

  }

  /*bool correctJets=false;

  //cout << systematic << " " << correctJets << endl;
  if (argc > 1) {
    if (atoi(argv[1]) >= 40) {
      systematic=atoi(argv[1])-40;
      correctJets=true;
      //cout << "oil" << endl;
    }
  }
  //cout << systematic << " " << correctJets << endl;

  //exit(1);*/


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

  cout << "systematic: " << systematic << " -> Postfix = " << postfix << endl;
	
	
	

//OLD WAY	
/*	
	//which systematic to run?
  string systematic = "Nominal";
/*  if (argc >= 2)
		systematic = string(argv[1]);
  cout << "Systematic to be used: " << systematic << endl;
  if( ! (systematic == "Nominal"  || systematic == "JESPlus" || systematic == "JESMinus" || systematic == "JERPlus" || systematic == "JERMinus") )
  {
    cout << "Unknown systematic!!!" << endl;
    cout << "Possible options are: Nominal JESPlus JESMinus JERPlus JERMinus" << endl;
    exit(-1);
  }
	string postfix = ""; // to relabel the names of the output trees  
	postfix = postfix+"_"+systematic;
*/
	//string Treespath = "VLQTrees_Summer12_PBS_v5_21Aug13"; //VLQTrees_Summer12_PBS if you use the PBB_VLQTreeCreator.py script (since it's hardcoded in there)!!
  //string Treespath = "VLQTrees_Summer12_5Mar14_FakeLepMCZ2and1Jets_rerun"; //"VLQTrees_Summer12_25Jan14"; //"VLQTrees_Summer12_PBS" //"VLQTrees_Summer12_PBS_v5_21Aug13_InvIsoRunABCD_forwardjets_PBS";
	//string Treespath = "VLQTrees_Summer12_26Jan14_FakeLepData";
	//string Treespath = "VLQTrees_Summer12_PBS";
	//string Treespath = "VLQTrees_Summer12_5Mar14_DataElD_mtop";
	string Treespath = "VLQTrees_Summer12_10Mar14_DataElD_METcut";
	Treespath = Treespath +"/";
  mkdir(Treespath.c_str(),0777);
	bool savePNG = false;



  clock_t start = clock();

  cout << "********************************************************" << endl;
  cout << " Beginning of the program for creating the VLQ Trees ! " << endl;
  cout << "********************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  string antibtagger = "CSVL";
  float antibtagWP = 9999;
	if(antibtagger == "CSVL")
     antibtagWP = 0.244 ;
  else if(antibtagger == "CSVM")
     antibtagWP = 0.679;
  else if(antibtagger == "CSVT")
     antibtagWP = 0.898;
		 
  /////////////////////
  // Configuration
  /////////////////////
  string rootFileName(Treespath+"/VLQTreeCreator_Plots_"+postfix+".root"); //some output, maybe dummy
	bool countTightLooseLeptons = true; //for data-driven fake-lepton contribution for multilepton events
	//cuts for data-driven fake-lepton contribution for multilepton events
	float LeadingJetPtCut = 200., SubLeadingJetPtCut = 100.;
	bool doZptcut = true;
	bool doMETcut = true;
	float Zptcut = 150;
	float METcut = 60;
	float massleptons = 9999.;
	float Zmass = 91.1876; //could eventually be replaced by the mean of a fit of the Z peak?
	float Zmasswindow = 7.5; //could eventually be replaced by the resolution of a fit of the Z peak (1 sigma)?
	//fake rates determined from VLQBkgEstimation.cc
	float Eff_TL_m = 0.0540;	
	float Eff_TL_m_unc = 0.0057;
	float Eff_TL_e = 0.0632;	
	float Eff_TL_e_unc = 0.0082;
	bool doAntiBtagging = true;
	
	if(countTightLooseLeptons)
	{
	   cout << "muon fake rate = " << Eff_TL_m << " +- " << Eff_TL_m_unc <<endl;
		 cout << "electron fake rate = " << Eff_TL_e << " +- " << Eff_TL_e_unc <<endl;
	}
	

	
  //xml file
  string xmlFileName ="../config/myVLQConfig.xml";	

  //if (systematic == 8)
  //  xmlFileName ="../config/myBTAGconfig_ttjetsyst.xml";
  //  //xmlFileName ="../config/myBTAGconfig_ttjetsyst_newcalib.xml";
  //else if (systematic == 9)
  //  xmlFileName ="../config/myBTAGconfig_wjetsyst.xml";
		
  if (argc > 3)
    xmlFileName = (string) argv[3];
  
  const char *xmlfile = xmlFileName.c_str();

  cout << "used config file: " << xmlfile << endl;

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
 
  cout << "analysis environment luminosity for rescaling (if not run on data) "<< oldLuminosity << endl;
  
  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  vector < Dataset* > datasetsMu;
  vector < Dataset* > datasetsEl;
  vector < Dataset* > datasetsPlot;

  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;

//  float LuminosityMu = oldLuminosity;
//  float LuminosityEl = oldLuminosity;  
	
	//hardcoded
	float LuminosityMu = 19600.84; //anaEnvLuminosity;
  float LuminosityEl = 19629.525; //anaEnvLuminosity;

  bool isSemiMu = false;
  bool isSemiE = false;

  bool foundMu = false;
  bool foundEl = false;

  for (unsigned int d = 0; d < datasets.size (); d++) {

    //if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
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

    //if( dataSetName.find("QCD") == 0 ) datasets[d]->SetColor(kYellow);
    //if( dataSetName.find("TTbar") == 0 ) datasets[d]->SetColor(8);
    
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
  
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //Global variable
  //TRootEvent* event = 0;
  
  //nof selected events
  double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
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
//  CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
//  CutsSelecTableSemiMu.push_back(string("Veto electron"));
  
  vector<string> CutsSelecTableSemiEl;
  CutsSelecTableSemiEl.push_back(string("preselected"));
  CutsSelecTableSemiEl.push_back(string("trigged"));
  CutsSelecTableSemiEl.push_back(string("Good PV"));
  CutsSelecTableSemiEl.push_back(string("1 selected electron"));
//  CutsSelecTableSemiEl.push_back(string("Veto muon"));
//  CutsSelecTableSemiEl.push_back(string("Veto 2nd electron from Z-decay"));
  CutsSelecTableSemiEl.push_back(string("Conversion veto"));
  
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

  //moved to VLQTreeAnalyzer.cc
	
	//even counts, but floats in case you run on MC and scale the events...
  float N_TmuTelLel_mu = 0, N_TmuTmuLel_mu = 0, N_TmuTelLmu_mu = 0, N_TmuTmuLmu_mu = 0;
	float N_TelTmuLmu_el = 0, N_TelTmuLel_el = 0, N_TelTelLmu_el = 0, N_TelTelLel_el = 0;
	float uncN_TmuTelLel_mu = 0, uncN_TmuTmuLel_mu = 0, uncN_TmuTelLmu_mu = 0, uncN_TmuTmuLmu_mu = 0;
	float uncN_TelTmuLmu_el = 0, uncN_TelTmuLel_el = 0, uncN_TelTelLmu_el = 0, uncN_TelTelLel_el = 0;
	
	
  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  bool changed_anaEnvCollections = false;
	
	string suffix = "";
	
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++) {

    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
		suffix = suffix + "_" +dataSetName;
		
    if (runSpecificSample && d != nSpecSample) continue;
		
    if (systematic == -1 && !(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) ) continue;
    
		if(dataSetName.find("InvIso") != string::npos && !changed_anaEnvCollections) {
			anaEnv.METCollection = anaEnv.METCollection+"NoLeptonCleaning";
			anaEnv.MuonCollection = anaEnv.MuonCollection+"NoLeptonCleaning";
			anaEnv.ElectronCollection = anaEnv.ElectronCollection+"NoLeptonCleaning";
			anaEnv.JetCollection = anaEnv.JetCollection+"NoLeptonCleaning";
			
			changed_anaEnvCollections = true;
    }
		
		int nSelectedMu=0;
    int nSelectedEl=0;
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    //open files and load
    cout<<"Load dataset"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"Dataset loaded"<<endl;
    
    
    
    /////////////////////////////////////
    /// Initialize JEC factors
    /////////////////////////////////////
   	    
    vector<JetCorrectorParameters> vCorrParam;
    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
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

    // true means redo also the L1
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
    
		//Trees
    string TreeFileName = Treespath+"VLQTree_"+dataSetName+postfix+".root";
    cout << "INFO: creating VLQTree file "+TreeFileName << endl;        
    TFile* treeFile;
		TTree* myVLQTree;
    VLQTree* myBranch_selectedEvents = 0;		
	  //if(!countTightLooseLeptons)
		//{
		   treeFile = new TFile(TreeFileName.c_str(),"RECREATE");
	     myVLQTree = new TTree("myVLQTree","Tree containing the VLQ analysis info");    
		   myVLQTree->Branch("VLQBranch_selectedEvents","VLQTree",&myBranch_selectedEvents);
		//}
		
    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    //for (unsigned int ievt = 0; ievt < 100; ievt++)
    {
      if (runSpecificSample && d != nSpecSample) continue;
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
            
      if(ievt%1000 == 0)
	         std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected: " << nSelectedMu << " (mu+jets) " << nSelectedEl << " (e+jets)" << flush<<"\r";
      //std::cout<<"Processing the "<<ievt<<"th event "<<endl;

      ////////////////
      // LOAD EVENT //
      ////////////////

      TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);  

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
      
      //////////////////////////////////////
      // Apply Jet Corrections on-the-fly //   
      //////////////////////////////////////

      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) 
			{
      	  jetTools->unCorrectMETTypeOne(init_jets_corrected, mets[0], true);
          jetTools->correctJets(init_jets_corrected, event->kt6PFJets_rho(), true); //last boolean: isData (needed for L2L3Residual...)
          jetTools->correctMETTypeOne(init_jets_corrected, mets[0], true);
      } 
			else 
			{
	        jetTools->unCorrectMETTypeOne(init_jets_corrected, mets[0], false);
          jetTools->correctJets(init_jets_corrected, event->kt6PFJets_rho(), false); //last boolean: isData (needed for L2L3Residual...)
          jetTools->correctMETTypeOne(init_jets_corrected, mets[0], false);
	    }


      ///////////////////
      // TRIGGER SETUP //
      ///////////////////

      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename)
			{
	         previousFilename = currentFilename;
	         iFile++;
	         cout<<"File changed!!! => iFile = "<<iFile<<", currentFilename = "<<currentFilename<<endl;
      }

      int currentRun = event->runId();

      if(previousRun != currentRun) {
        previousRun = currentRun;
	
	//semi-mu
	if((dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) || (dataSetName.find("InvIso_Mu") != string::npos)) {
	  //cout<<"I AM HERE (mu)"<<endl;   

	  if( event->runId() <= 190738 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);

	  else if( event->runId() >= 190782 && event->runId() <= 193621 )  //event->runId() >= 191043 && event->runId() <= 193621
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
            exit(-1); //in leptonnocleaning skims there are some files with missing runtree, i don't want the code to stop... well, commenting doesnt help
          }

	} 
	//else if (!((dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) || (dataSetName.find("InvIso_El") != string::npos) ))
  else
	{
	  itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun); // Summer12 DR53X
	}

	// semi-el
	// semi-electron
  if((dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) || (dataSetName.find("InvIso_El") != string::npos) ) {
         
	  if( event->runId() <= 190738 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v8"), currentRun, iFile);
	  
	  else if( event->runId() >= 190782 && event->runId() <= 191411 )  //event->runId() >= 191043 && event->runId() <= 191411
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
	      exit(-1); //in leptonnocleaning skims there are some files with missing runtree, i don't want the code to stop... well, commenting doesnt help
	  }
   }
	 //else if (!((dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) || (dataSetName.find("InvIso_Mu") != string::npos) ))
   else
	 {
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun); // Summer12 DR53X 
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

	         jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal",false);

           if(systematic == 6)
	           jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus",false); //false means don't use old numbers but newer ones...
	         else if(systematic == 7)
	           jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus",false);
	         else
	           jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal",false);
	         
					 // Example how to apply JES systematics	
	         if (systematic == 1)
	           jetTools->correctJetJESUnc(init_jets_corrected, mets[0], "minus",1);
	         else if (systematic == 2)
	           jetTools->correctJetJESUnc(init_jets_corrected, mets[0], "plus",1);
						 
					 if (systematic == 3)
	           jetTools->correctMETUnclusteredEnergy(mets[0],init_jets_corrected, init_muons, init_electrons, "minus");
	         else if (systematic == 4)
	           jetTools->correctMETUnclusteredEnergy(mets[0],init_jets_corrected, init_muons, init_electrons, "plus");
      }

      /////////////////////
      // EVENT SELECTION //
      /////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho());
			Selection nonstandard_selection(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho()); //mets can also be corrected... 
      Selection selectionFakeLepton(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho()); //mets can also be corrected...
			
			selection.cutHFHadronEnergyFraction();
			nonstandard_selection.cutHFHadronEnergyFraction();
			selectionFakeLepton.cutHFHadronEnergyFraction();
			
			selection.setJetCuts(30,2.4,0.01,1.,0.98,0.3,0.1); //was 2.5, but adapted to 2.4 for b-tagging... correct?
			nonstandard_selection.setJetCuts(30.,4.7,0.01,1.,0.98,0.3,0.1); //only difference: larger eta acceptance 
      selection.setMuonCuts(20,2.1,0.12,0.2,0.3,1,0.5,5,0); //Pt 20 GeV! later on, higher for leading lepton at least
      selection.setElectronCuts(20,2.5,0.1,0.02,0.5,0.3,0); //was pt 32... //Pt 20 GeV! later on, higher for leading lepton at least
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setLooseElectronCuts(15,2.5,0.15,0.); //Pt 15 GeV! first I had pt 20 -> ah damn, apparently I overwrote it again later with selection.GetSelectedLooseElectrons(20,2.5,0.15)... ?! As of 22Jan14, I will remake all trees, with putting this to 15 later on too
			float LeadingLeptonPtCut = 30.;
			
      bool triggedSemiMu = false;
      bool triggedSemiEl = false;

      //to be checked: here there is/was no mention of what to do in the InvIso case, but that was dangerous... it is still correct, but looks like a lucky accident -> update: I get HLT path 9999 not found when running on InvIso; should be related... is it a problem? Now I put the InvIso criterion...
      if( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0 || (dataSetName.find("InvIso_El") != string::npos)) )
        triggedSemiMu = treeLoader.EventTrigged (itriggerSemiMu);
      if( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0 || (dataSetName.find("InvIso_Mu") != string::npos)) )
        triggedSemiEl = treeLoader.EventTrigged (itriggerSemiEl);

      bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);


      vector<TRootJet*> selectedJets, selectedForwardJets, selectedJetsLargeEtaRange, selectedJetsNoMu, selectedJetsNoEl, selectedJetsLargeEtaRangeNoMu, selectedJetsLargeEtaRangeNoEl, selectedForwardJetsNoMu, selectedForwardJetsNoEl;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
			vector<TRootMuon*> selectedLooseMuons = selection.GetSelectedLooseMuons(10,2.5,0.2); //to veto; most probably selection.GetSelectedLooseMuons() is sufficient, the cuts were already set by setLooseMuonCuts
      //vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons(10,2.5,0.2);
			vector<TRootElectron*> selectedLooseElectrons = selection.GetSelectedLooseElectrons(15,2.5,0.15); //to veto; hmm 20 GeV is tighter than the cuts by setLooseElectronsCuts; I'm afraid this should have been GetSelectedLooseElectrons(15,2.5,0.15) ?! As of 22Jan14, I will remake all trees, with putting this to 15
      //vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(20,2.5,0.15);
      //vector<TRootElectron*> vetoElectronsSemiEl = selection.GetSelectedLooseElectrons(20,2.5,0.15);

      if( dataSetName.find("InvIso") != string::npos )  { // event selection for special Data TopTrees for ARI QCD
        vector<TRootMuon*> overlapMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0]);
        vector<TRootElectron*> overlapElectrons = selection.GetSelectedElectronsInvIso(0.2);
        selectedJetsNoMu = selection.GetSelectedJets(overlapMuons,true);
        selectedJetsNoEl = selection.GetSelectedJets(overlapElectrons,true);

	      selectedMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0], selectedJetsNoMu);
        selectedElectrons = selection.GetSelectedElectronsInvIso(0.2,selectedJetsNoEl); // Also mvaId < 0.5
				
				
				//(z)
				//something with selectedJetsLargeEtaRangeNoMu etc or not??
				selectedJetsLargeEtaRangeNoMu = nonstandard_selection.GetSelectedJets(overlapMuons,true);
        selectedJetsLargeEtaRangeNoEl = nonstandard_selection.GetSelectedJets(overlapElectrons,true);
				for(unsigned int i = 0; i<selectedJetsLargeEtaRangeNoMu.size(); i++)
      	{
	    		  if(fabs(selectedJetsLargeEtaRangeNoMu[i]->Eta()) > 2.4)
	      		  selectedForwardJetsNoMu.push_back(selectedJetsLargeEtaRangeNoMu[i]);
      	}
				for(unsigned int i = 0; i<selectedJetsLargeEtaRangeNoEl.size(); i++)
      	{
	    		  if(fabs(selectedJetsLargeEtaRangeNoEl[i]->Eta()) > 2.4)
	      		  selectedForwardJetsNoEl.push_back(selectedJetsLargeEtaRangeNoEl[i]);
      	}
      }
      else { // Normal selection
	       selectedJets = selection.GetSelectedJets(true);
				 selectedJetsLargeEtaRange = nonstandard_selection.GetSelectedJets(true);
				 for(unsigned int i = 0; i<selectedJetsLargeEtaRange.size(); i++)
      	 {
	    		  if(fabs(selectedJetsLargeEtaRange[i]->Eta()) > 2.4) //was 2.5, but adapted to 2.4 for b-tagging... correct?
	      		  selectedForwardJets.push_back(selectedJetsLargeEtaRange[i]);
      	 }

	       //selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);	
	       selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets); //no large eta range things here??
	       selectedElectrons = selection.GetSelectedElectrons(selectedJets);
      }
			
			
			//for data-driven fake-lepton contribution for multilepton events -> follows roughly InclFourthGen_TreeCreator.cc (but with modifications)
			vector<TRootMuon*> selectedMuons_FL, selectedLooseMuons_FL, selectedOnlyLooseMuons_FL;
			vector<TRootElectron*> selectedElectrons_FL, selectedLooseElectrons_FL, selectedOnlyLooseElectrons_FL;
      if(countTightLooseLeptons)
			{				
				selectionFakeLepton.setMuonCuts(20,2.1,0.12,0.2,0.3,1,0.5,5,0); //Pt 20 GeV! later on, higher for leading lepton at least
        selectionFakeLepton.setElectronCuts(20,2.5,0.1,0.02,0.5,0.3,0);						
				selectedMuons_FL = selectionFakeLepton.GetSelectedMuons(vertex[0],selectedJets);
	    	selectedElectrons_FL = selectionFakeLepton.GetSelectedElectrons(selectedJets);				
      	selectionFakeLepton.setLooseMuonCuts(10,2.5,0.2); //0.2 reliso is the maximum in our toptrees...
      	selectionFakeLepton.setLooseElectronCuts(15,2.5,0.2,0.); 								
      	selectedLooseElectrons_FL = selectionFakeLepton.GetSelectedLooseElectrons(15,2.5,0.2); //for data-driven background estimation //ATTENTION: potentially overwrites what you put in setLooseElectronCuts
      	selectedLooseMuons_FL = selectionFakeLepton.GetSelectedLooseMuons(); //for data-driven background estimation
				
				for(unsigned int i = 0; i<selectedLooseMuons_FL.size(); i++)
      	{
      		int IsTight = 0;
      		for(unsigned int j = 0; j<selectedMuons.size(); j++)
      		{
						if(fabs(selectedLooseMuons_FL[i]->Pt()-selectedMuons[j]->Pt()) < 0.00001) // the loose muon is tight
							IsTight = IsTight+1;
					}
					if(IsTight==0) selectedOnlyLooseMuons_FL.push_back(selectedLooseMuons_FL[i]); // we found no tight muon, so this loose muon is not tight!
				}
				for(unsigned int i = 0; i<selectedLooseElectrons_FL.size(); i++)
      	{
      		int IsTight = 0;
      		for(unsigned int j = 0; j<selectedElectrons.size(); j++)
      		{
						if(fabs(selectedLooseElectrons_FL[i]->Pt()-selectedElectrons[j]->Pt()) < 0.00001)
							IsTight = IsTight+1;
					}
					if(IsTight==0) selectedOnlyLooseElectrons_FL.push_back(selectedLooseElectrons_FL[i]);
				}
				//cout<<"#tight muons = "<<selectedMuons.size()<<", #tight electrons = "<<selectedElectrons.size()<<endl;
				//cout<<"#loose muons = "<<selectedLooseMuons.size()<<", #loose electrons = "<<selectedLooseElectrons.size()<<endl;
				//cout<<"#loose-not-tight muons = "<<selectedOnlyLooseMuons_FL.size()<<", #loose-not-tight electrons = "<<selectedOnlyLooseElectrons_FL.size()<<endl;
			}
			
			
			
      vector<TRootMCParticle*> mcParticles;
  
      //if(dataSetName.find("TTbarJets") == 0)
			if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) )
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }





/*
      /////////////////////////////////
      //TESTING MC particles 
			unsigned int nQD = 0, nQU = 0;
			unsigned int nQTojW = 0, nQTojZ = 0, nQTojH = 0;
      for(unsigned int u=0; u<mcParticles.size(); u++)
			{
				 unsigned int pdgidQD = 4000001;
				 unsigned int pdgidQU = 4000002;
			   if(abs(mcParticles[u]->type()) == pdgidQD)
				 {
	          nQD++;
						if(abs(mcParticles[u]->dauOneId()) == 24 || abs(mcParticles[u]->dauTwoId()) == 24)
						{
						   nQTojW++;
						}
						if(abs(mcParticles[u]->dauOneId()) == 23 || abs(mcParticles[u]->dauTwoId()) == 23)
						{
						   nQTojZ++;
						}
						if(abs(mcParticles[u]->dauOneId()) == 25 || abs(mcParticles[u]->dauTwoId()) == 25)
						{
						   nQTojH++;
						}
						//cout<<"       number of daughters = "<<mcParticles[u]->nDau()<<"; daughters pdg id: "<<mcParticles[u]->dauOneId()<<", "<<mcParticles[u]->dauTwoId()<<", "<<mcParticles[u]->dauThreeId()<<", "<<mcParticles[u]->dauFourId()<<", "<<endl;
				 }
				 if(abs(mcParticles[u]->type()) == pdgidQU)
				 {
            nQU++;
						if(abs(mcParticles[u]->dauOneId()) == 24 || abs(mcParticles[u]->dauTwoId()) == 24)
						{
						   nQTojW++;
						}
						if(abs(mcParticles[u]->dauOneId()) == 23 || abs(mcParticles[u]->dauTwoId()) == 23)
						{
						   nQTojZ++;
						}
						if(abs(mcParticles[u]->dauOneId()) == 25 || abs(mcParticles[u]->dauTwoId()) == 25)
						{
						   nQTojH++;
						}
				 }
				 
			}				
			
      if(nQD!=0 || nQU!=0) cout<<"   number of QD quarks in event = "<<nQD<<", number of QU quarks in event = "<<nQU<<endl;
			if(nQD!=0)
			{
			   if(nQD==1)
				 {
				    if(nQTojW==1)
				       cout<<"     --> QDTojW event"<<endl;
						else if(nQTojZ==1)
				       cout<<"     --> QDTojZ event"<<endl;
						else
						   cout<<"     --> WARNING: event not identified!!"<<endl;				 
				 }
				 else if(nQD==2)
				 {
				    if(nQTojW==2)
				       cout<<"     --> QDQDTojWjW event"<<endl;
						else if(nQTojZ==2)
				       cout<<"     --> QDQDTojZjZ event"<<endl;
						else if(nQTojW==1 && nQTojZ==1)
				       cout<<"     --> QDQDTojWjZ event"<<endl;
						else if(nQTojW==1 && nQTojH==1)
				       cout<<"     --> QDQDTojWjH event"<<endl;
						else if(nQTojZ==1 && nQTojH==1)
				       cout<<"     --> QDQDTojZjH event"<<endl;  
						else
						   cout<<"     --> WARNING: event not identified!!"<<endl;				 
				 }
				 else
				    cout<<"     --> WARNING: event not identified!!"<<endl;					
			}
			else if(nQU!=0)
			{
			   if(nQU==1)
				 {
				    if(nQTojW==1)
				       cout<<"     --> QUTojW event"<<endl;
						else if(nQTojZ==1)
				       cout<<"     --> QUTojZ event"<<endl;
						else
						   cout<<"     --> WARNING: event not identified!!"<<endl;				 
				 }
				 else if(nQU==2)
				 {
				    //no QU pair events generated...
				    if(nQTojW==2)
				       cout<<"     --> QUQUTojWjW event"<<endl;
						else if(nQTojZ==2)
				       cout<<"     --> QUQUTojZjZ event"<<endl;
						else if(nQTojW==1 && nQTojZ==1)
				       cout<<"     --> QUQUTojWjZ event"<<endl;
						else if(nQTojW==1 && nQTojH==1)
				       cout<<"     --> QUQUTojWjH event"<<endl;
						else if(nQTojZ==1 && nQTojH==1)
				       cout<<"     --> QUQUTojZjH event"<<endl;  
						else
						   cout<<"     --> WARNING: event not identified!!"<<endl;				 
				 }
				 else
				    cout<<"     --> WARNING: event not identified!!"<<endl;					
			}
*/





      bool eventselectedSemiMu = false;
      bool eventselectedSemiEl = false;

      unsigned int SelectednLeptons = 0, SelectednMu = 0, SelectednEl = 0;
			unsigned int SelectednLooseNotTightMu = 0, SelectednLooseNotTightEl = 0;
			
			//these booleans will never be true for the same event, unlike the trigger bit that can be true for muon trigged and electron trigged
			bool eventselected_MuChannel = false;
			bool eventselected_ElChannel = false;
			//these booleans will never be true for the same event, unlike the trigger bit that can be true for muon trigged and electron trigged
			bool eventpseudoselected_MuChannel = false; //for data-driven multilepton fake rate estimation (let's say events with 2 tight leptons and 1 loose-not-tight lepton, or 3 tight leptons and 1 loose-not-tight lepton) 
      bool eventpseudoselected_ElChannel = false;
		
      // selection with mu
      if (triggedSemiMu && (!((dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) || (dataSetName.find("InvIso_El") != string::npos) )))
			{
			  //cout<<"-> passing mu trigger"<<endl;
				if (isGoodPV) {
				  //cout<<" --> good PV"<<endl;
	  			if(dataSetName.find("InvIso") != string::npos)
					{
					   selectedJets = selectedJetsNoMu;
					   selectedJetsLargeEtaRange = selectedJetsLargeEtaRangeNoMu;
						 selectedForwardJets = selectedForwardJetsNoMu;
					}
					if(selectedJetsLargeEtaRange.size() >= 2)
					{	
					  //cout<<"  ---> >=2 jets (large eta range)"<<endl;
						/*//the electrons should not be from a conversion! -> ok, apparently now (in 53X) done in selection class...
						vector<TRootElectron*> selectedElectrons_noconversions;
						for(unsigned int iEl=0; iEl<selectedElectrons.size(); iEl++)
						{
								if(selection.passConversionRejection(selectedElectrons[iEl])
								    selectedElectrons_noconversions.push_back(selectedElectrons[iEl]);					   
						}
						selectedElectrons.clear();
						for(unsigned int iEl=0; iEl<selectedElectrons_noconversions.size(); iEl++)
						{
								    selectedElectrons.push_back(selectedElectrons_noconversions[iEl]);					   
						}*/
						
						SelectednMu = selectedMuons.size();
						SelectednEl = selectedElectrons.size();
						
						SelectednLooseNotTightMu = selectedOnlyLooseMuons_FL.size();
						SelectednLooseNotTightEl = selectedOnlyLooseElectrons_FL.size();
						
									 
						if(SelectednMu == 0)
						  continue;
						if(selectedMuons[0]->Pt() < LeadingLeptonPtCut)
						  continue;
							
						//cout<<"   ----> SelectednMu>=1 and selectedMuons[0]->Pt() >= LeadingLeptonPtCut"<<endl;
						if(SelectednMu == 1 && (selectedLooseMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 ))) // if InvertedIso, selected muon not part ofselected LooseMuons vector! (to be checked)
						{	  				  
							
							//if( vetoMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && vetoMuons.size() == 0 ) ) { // if InvertedIso, selected muon not part of vetoMuons vector!
	  				  //  if (vetoElectronsSemiMu.size() == 0) {
							//     if (selectedJets.size() >= 4) {
							//        eventselected_MuChannelSemiMu = true;
						  //    }
	            //  }
	            //}
							if(SelectednEl == 0 && (selectedLooseElectrons.size() == 0) )
							{
							     eventselected_MuChannel = true;
							}
							else if(SelectednEl == 1 && (selectedLooseElectrons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{
							     eventselected_MuChannel = true;
							}
							else if(SelectednEl == 2 && (selectedLooseElectrons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{							    
							     eventselected_MuChannel = true;
							}
							else if(SelectednEl == 3 && (selectedLooseElectrons.size() == 3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{
							     eventselected_MuChannel = true;								 
							}
							else if(SelectednEl >= 4 && (selectedLooseElectrons.size() >=  4 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )) )
							{
							     cout<<"FOUND >= 5-lepton event!"<<endl;		
							}
	         }
					 else if(SelectednMu == 2 && (selectedLooseMuons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 ))) 
					 {
					    if(SelectednEl == 0 && selectedLooseElectrons.size() == 0)
							{
							     eventselected_MuChannel = true;
							}
							else if(SelectednEl == 1 && (selectedLooseElectrons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{
							     eventselected_MuChannel = true;
							}
							else if(SelectednEl == 2 && (selectedLooseElectrons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{
							     eventselected_MuChannel = true;
							}
							else if(SelectednEl >= 3 && (selectedLooseElectrons.size() >=  3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							     cout<<"FOUND >= 5-lepton event!"<<endl;		 
					 }
					 else if(SelectednMu == 3 && (selectedLooseMuons.size() == 3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 ))) 
					 {
					    if(SelectednEl == 0 && selectedLooseElectrons.size() == 0)
							{
							     eventselected_MuChannel = true;
							}
							else if(SelectednEl == 1 && (selectedLooseElectrons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{
							     eventselected_MuChannel = true;
							}
							else if(SelectednEl >= 2 && (selectedLooseElectrons.size() >=  2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							  cout<<"FOUND >= 5-lepton event!"<<endl;
					 
					 }
					 else if(SelectednMu == 4 && (selectedLooseMuons.size() == 4 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 ))) 
					 {
					    if(SelectednEl == 0 && selectedLooseElectrons.size() == 0)
							{
							     eventselected_MuChannel = true;
							}
							else if(SelectednEl >= 1 && (selectedLooseElectrons.size() >=  1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							  cout<<"FOUND >= 5-lepton event!"<<endl;
					 
					 }
					 else if(SelectednMu >= 5 && (selectedLooseMuons.size() >=  5 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
					 {
							  cout<<"FOUND >= 5-lepton event!"<<endl;
					 }
					 //else if(countTightLooseLeptons) //&& dataSetName.find("Data")==0)   //->'else if' was bug
					 if(countTightLooseLeptons)
					 { 
					    //cout<<"    -----> counting tight and loose leptons"<<endl;
					    //for data-driven multilepton fake lepton estimation
							//for trilepton estimation
							if(SelectednMu == 1 && SelectednEl == 1)
							{
							  //cout<<"     ------> SelectednMu == 1 && SelectednEl == 1"<<endl;
							  if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1)
								{
								   //cout<<"      -------> selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1"<<endl;
								   eventpseudoselected_MuChannel = true;	//(Tmu,Tel,Lmu)							
								}
								else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1)
								{
								   //cout<<"      -------> selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1"<<endl;
								   eventpseudoselected_MuChannel = true;	//(Tmu,Tel,Lel)							
								}								
							}
							else if(SelectednMu == 2 && SelectednEl == 0)
							{
							  //cout<<"     ------> SelectednMu == 2 && SelectednEl == 0"<<endl;
							  //if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1) //->was bug
								if(selectedLooseElectrons.size() == 0 && SelectednLooseNotTightMu == 1)
								{
								   //cout<<"      -------> selectedLooseElectrons.size() == 0 && SelectednLooseNotTightMu == 1"<<endl;
								   eventpseudoselected_MuChannel = true;	//(Tmu,Tmu,Lmu)							
								}
								//else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1) //->was bug
								else if(selectedLooseMuons.size() == 2 && SelectednLooseNotTightEl == 1)
								{
								   //cout<<"      -------> selectedLooseMuons.size() == 2 && SelectednLooseNotTightEl == 1"<<endl;
								   eventpseudoselected_MuChannel = true;	//(Tmu,Tmu,Lel)							
								}								
							}
							//for fourlepton estimation -> euhm, estimating from 3-leptons? But that is a signal category itself??
							else if(SelectednMu == 1 && SelectednEl == 2)
							{
							  //if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1) //->was bug
								if(selectedLooseElectrons.size() == 2 && SelectednLooseNotTightMu == 1)
								{
								   eventpseudoselected_MuChannel = true;	//(Tmu,Tel,Tel,Lmu)							
								}
								//else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1) //->was bug
								else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1)
								{
								   eventpseudoselected_MuChannel = true;	//(Tmu,Tel,Tel,Lel)							
								}								
							}
							else if(SelectednMu == 2 && SelectednEl == 1)
							{
							  //if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1) //->was bug
								if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1)
								{
								   eventpseudoselected_MuChannel = true;	//(Tmu,Tmu,Tel,Lmu)							
								}
								//else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1) //->was bug
								if(selectedLooseMuons.size() == 2 && SelectednLooseNotTightEl == 1)
								{
								   eventpseudoselected_MuChannel = true;	//(Tmu,Tmu,Tel,Lel)							
								}								
							}
							else if(SelectednMu == 3 && SelectednEl == 0)
							{
							  //if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1) //->was bug
								if(selectedLooseElectrons.size() == 0 && SelectednLooseNotTightMu == 1) 
								{
								   eventpseudoselected_MuChannel = true;	//(Tmu,Tmu,Tmu,Lmu)							
								}
								//else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1) //->was bug
								else if(selectedLooseMuons.size() == 3 && SelectednLooseNotTightEl == 1)
								{
								   eventpseudoselected_MuChannel = true;	//(Tmu,Tmu,Tmu,Lel)							
								}								
							}
							if(eventpseudoselected_MuChannel)
							{
							  cout<<"  [mu-channel]"<<endl;
							  cout<<"   ==> eventpseudoselected_MuChannel"<<endl;
					      cout<<"   ---> SelectednMu = "<<SelectednMu<<", SelectednEl = "<<SelectednEl<<endl;
							  cout<<"        selectedLooseMuons.size() = "<<selectedLooseMuons.size()<<", selectedLooseElectrons.size() = "<<selectedLooseElectrons.size()<<endl;
							  cout<<"        SelectednLooseNotTightMu = "<<SelectednLooseNotTightMu<<", SelectednLooseNotTightEl = "<<SelectednLooseNotTightEl<<endl;
							}							 
					 }
					 
				 } //nJets
	      } //good PV
      }
			//else if or if?? double counting will be avoided anyway with the lepton pt requiements, right? 'if' might be better ?? (first I had 'else if', but I'm doubting that a bit). 
			//'Else if' in practice means that in the following, events are prohibited to have a triggered muon. So an event with one trigged muon of 28 GeV and a trigged electron of 32 GeV, will not pass our selection (because in the muon channel the muon should be above 30 GeV, and with else if it can't even get passed here. But it should be selected! I think it's no drama when you do that on data and MC, but you lose multilepton events like that...)
      if(!eventpseudoselected_MuChannel && triggedSemiEl && (!((dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) || (dataSetName.find("InvIso_Mu") != string::npos) )))
			{		
			  //cout<<"-> passing el trigger"<<endl;
        if(dataSetName.find("InvIso") != string::npos)
				{
				    selectedJets = selectedJetsNoEl;
						selectedJetsLargeEtaRange = selectedJetsLargeEtaRangeNoEl;
						selectedForwardJets = selectedForwardJetsNoEl;
				}
        if (isGoodPV ) {
					if(selectedJetsLargeEtaRange.size() >= 2)
				  {	
					  /*//the electrons should not be from a conversion! -> ok, apparently now (in 53X) done in selection class...
						vector<TRootElectron*> selectedElectrons_noconversions;
						for(unsigned int iEl=0; iEl<selectedElectrons.size(); iEl++)
						{
								if(selection.passConversionRejection(selectedElectrons[iEl])
								    selectedElectrons_noconversions.push_back(selectedElectrons[iEl]);					   
						}
						selectedElectrons.clear();
						for(unsigned int iEl=0; iEl<selectedElectrons_noconversions.size(); iEl++)
						{
								    selectedElectrons.push_back(selectedElectrons_noconversions[iEl]);					   
						}*/
						
						SelectednMu = selectedMuons.size();
						SelectednEl = selectedElectrons.size();
						
						if(SelectednEl == 0)
						  continue;
						if(selectedElectrons[0]->Pt() < LeadingLeptonPtCut)
						  continue;	
						 
	          if(SelectednEl == 1  && (selectedLooseElectrons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 ))) {
				 
	              //if( vetoMuons.size() == 0 ) {
	              //    if (vetoElectronsSemiEl.size() == 1) {
	              // }
	              //}
							
								if(SelectednMu == 0 && selectedLooseMuons.size() ==  0)
							  {
							     eventselected_ElChannel = true;
							  }
								else if(SelectednMu == 1 && (selectedLooseMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								{
								   if(selectedMuons[0]->Pt() < LeadingLeptonPtCut) //to avoid double counting with muon channel => question; was double counting still possible when requiring 'else if electron trigged' instead of 'if electron trigged'??
									 {
									   eventselected_ElChannel = true;
									 }
								}
								else if(SelectednMu == 2 && (selectedLooseMuons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								{
								   if((selectedMuons[0]->Pt() < LeadingLeptonPtCut) && (selectedMuons[1]->Pt() < LeadingLeptonPtCut)) //to avoid double counting with muon channel
									 {
									   eventselected_ElChannel = true;
									 }
								}
								else if(SelectednMu == 3 && (selectedLooseMuons.size() == 3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								{
								   if((selectedMuons[0]->Pt() < LeadingLeptonPtCut) && (selectedMuons[1]->Pt() < LeadingLeptonPtCut) && (selectedMuons[2]->Pt() < LeadingLeptonPtCut)) //to avoid double counting with muon channel
									 {
									    eventselected_ElChannel = true;
									 }
								}
								else if(SelectednMu >= 4 && (selectedLooseMuons.size() >=  4 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
							     cout<<"FOUND >= 5-lepton event!"<<endl;
							  																
	          }
						else if(SelectednEl == 2  && (selectedLooseElectrons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 ))) {
							
								if(SelectednMu == 0 && selectedLooseMuons.size() ==  0)
							  {
							     eventselected_ElChannel = true;
							  }
								else if(SelectednMu == 1 && (selectedLooseMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )) )
								{
								   if(selectedMuons[0]->Pt() < LeadingLeptonPtCut) //to avoid double counting with muon channel
									 {
									    eventselected_ElChannel = true;
									 }
								}
								else if(SelectednMu == 2 && (selectedLooseMuons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								{
								   if((selectedMuons[0]->Pt() < LeadingLeptonPtCut) && (selectedMuons[1]->Pt() < LeadingLeptonPtCut)) //to avoid double counting with muon channel
									 {
									    eventselected_ElChannel = true;
									 }
								}
								else if(SelectednMu >= 3 && (selectedLooseMuons.size() >=  3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
							     cout<<"FOUND >= 5-lepton event!"<<endl;		  																
	          }
						else if(SelectednEl == 3  && (selectedLooseElectrons.size() == 3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 ))) {
							
								if(SelectednMu == 0 && selectedLooseMuons.size() ==  0)
							  {
							     eventselected_ElChannel = true;
							  }
								else if(SelectednMu == 1 && (selectedLooseMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								{
								   if(selectedMuons[0]->Pt() < LeadingLeptonPtCut) //to avoid double counting with muon channel
									 {
									    eventselected_ElChannel = true;
									 }
								}
								else if(SelectednMu >= 2 && (selectedLooseMuons.size() >= 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
							     cout<<"FOUND >= 5-lepton event!"<<endl;						  																
	          }
						else if(SelectednEl == 4  && selectedLooseElectrons.size() == 4) {
							
								if(SelectednMu == 0 && selectedLooseMuons.size() ==  0)
							  {
							     eventselected_ElChannel = true;
							  }
								else if(SelectednMu >= 1 && (selectedLooseMuons.size() >= 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								   cout<<"FOUND >= 5-lepton event!"<<endl;							  																
	          }
						else if(SelectednEl >= 5 && (selectedLooseElectrons.size() >= 5 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
						{
							  cout<<"FOUND >= 5-lepton event!"<<endl;
						}
						//else if(countTightLooseLeptons) // && dataSetName.find("Data")==0)   //->'else if' was bug
					  if(countTightLooseLeptons)
						{ 
						  //for data-driven multilepton fake lepton estimation
							//for trilepton estimation
							if(SelectednEl == 1 && SelectednMu == 1)
							{
								if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1)
								{
								   if((selectedMuons[0]->Pt() < LeadingLeptonPtCut) && (selectedOnlyLooseMuons_FL[0]->Pt() < LeadingLeptonPtCut)) //to avoid double counting with muon channel (AND to be consistent with signal selection!)
									 {
								      eventpseudoselected_ElChannel = true;	//(Tel,Tmu,Lmu)		
									 }					
								}
								else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1)
								{
								   if(selectedMuons[0]->Pt() < LeadingLeptonPtCut) //to avoid double counting with muon channel
									 {
								      eventpseudoselected_ElChannel = true;	//(Tel,Tmu,Lel)
									 }							
								}								
							}
							else if(SelectednEl == 2 && SelectednMu == 0)
							{
							  //if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1) //->was bug
								if(selectedLooseElectrons.size() == 2 && SelectednLooseNotTightMu == 1)
								{
								   if(selectedOnlyLooseMuons_FL[0]->Pt() < LeadingLeptonPtCut) //to avoid double counting with muon channel (rather to be consistent with signal selection!)
									 {
								      eventpseudoselected_ElChannel = true;	//(Tel,Tel,Lmu)
									 }							
								}
								//else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1) //->was bug
								else if(selectedLooseMuons.size() == 0 && SelectednLooseNotTightEl == 1)
								{
								   eventpseudoselected_ElChannel = true;	//(Tel,Tel,Lel)							
								}								
							}
							//for fourlepton estimation
							else if(SelectednEl == 1 && SelectednMu == 2)
							{
							  if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1)
								{
								   if((selectedMuons[0]->Pt() < LeadingLeptonPtCut) && (selectedMuons[1]->Pt() < LeadingLeptonPtCut) && (selectedOnlyLooseMuons_FL[0]->Pt() < LeadingLeptonPtCut)) //to avoid double counting with muon channel (AND to be consistent with signal selection!)
									 {
								      eventpseudoselected_ElChannel = true;	//(Tel,Tmu,Tmu,Lmu)	
									 }						
								}
								//else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1) //->was bug
								else if(selectedLooseMuons.size() == 2 && SelectednLooseNotTightEl == 1)
								{
								   if((selectedMuons[0]->Pt() < LeadingLeptonPtCut) && (selectedMuons[1]->Pt() < LeadingLeptonPtCut)) //to avoid double counting with muon channel
									 {
								      eventpseudoselected_ElChannel = true;	//(Tel,Tmu,Tmu,Lel)		
									 }					
								}								
							}
							else if(SelectednEl == 2 && SelectednMu == 1)
							{
							  //if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1) //->was bug
								if(selectedLooseElectrons.size() == 2 && SelectednLooseNotTightMu == 1)
								{
								   if((selectedMuons[0]->Pt() < LeadingLeptonPtCut) && (selectedOnlyLooseMuons_FL[0]->Pt() < LeadingLeptonPtCut)) //to avoid double counting with muon channel (AND to be consistent with signal selection!)
									 {
								      eventpseudoselected_ElChannel = true;	//(Tel,Tel,Tmu,Lmu)				
									 }			
								}
								else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1)
								{
								   if(selectedMuons[0]->Pt() < LeadingLeptonPtCut) //to avoid double counting with muon channel
									 {
								      eventpseudoselected_ElChannel = true;	//(Tel,Tel,Tmu,Lel)		
									 }					
								}								
							}
							else if(SelectednEl == 3 && SelectednMu == 0)
							{
							  //if(selectedLooseElectrons.size() == 1 && SelectednLooseNotTightMu == 1) //->was bug
								if(selectedLooseElectrons.size() == 3 && SelectednLooseNotTightMu == 1)
								{
								   if(selectedOnlyLooseMuons_FL[0]->Pt() < LeadingLeptonPtCut) //to avoid double counting with muon channel (rather to be consistent with signal selection!)
									 {
								      eventpseudoselected_ElChannel = true;	//(Tel,Tel,Tel,Lmu)				
									 }			
								}
								//else if(selectedLooseMuons.size() == 1 && SelectednLooseNotTightEl == 1) //->was bug
								else if(selectedLooseMuons.size() == 0 && SelectednLooseNotTightEl == 1)
								{
								   eventpseudoselected_ElChannel = true;	//(Tel,Tel,Tel,Lel)							
								}								
							}	
							if(eventpseudoselected_ElChannel)
							{
							  cout<<"  [el-channel]"<<endl;
							  cout<<"   ==> eventpseudoselected_ElChannel"<<endl;
					      cout<<"   ---> SelectednMu = "<<SelectednMu<<", SelectednEl = "<<SelectednEl<<endl;
							  cout<<"        selectedLooseMuons.size() = "<<selectedLooseMuons.size()<<", selectedLooseElectrons.size() = "<<selectedLooseElectrons.size()<<endl;
							  cout<<"        SelectednLooseNotTightMu = "<<SelectednLooseNotTightMu<<", SelectednLooseNotTightEl = "<<SelectednLooseNotTightEl<<endl;
							}								 
					 }
												
				 } //nJets
       } //good PV
     }
		 
		 SelectednLeptons = SelectednMu + SelectednEl;
     
		 if(!(eventselected_MuChannel || eventselected_ElChannel) && !(eventpseudoselected_MuChannel || eventpseudoselected_ElChannel)) continue;

		if(eventselected_MuChannel || eventselected_ElChannel)
		{
     vector<float> SelectedJets_bTagCSV, SelectedForwardJets_bTagCSV;
		 vector<int> SelectedJets_partonFlavour, SelectedForwardJets_partonFlavour;
     vector<TLorentzVector> mcParticlesTLV, SelectedJetsTLV, SelectedForwardJetsTLV, SelectedMuonsTLV, SelectedElectronsTLV;
		 vector<float> SelectedMuonsRelIso, SelectedElectronsRelIso;
		 vector<int> SelectedMuonsCharge, SelectedElectronsCharge;

     /*for(unsigned int i=0; i<mcParticles.size(); i++) 
		 {
						//if( mcParticles[i]->status() != 3) continue;
						
						//if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 )
						//{
							mcParticlesTLV.push_back(*mcParticles[i]);
						//}
		 }*/
		 
		 // reweight pttop for madgraph; see slide 8 https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=254297 and https://github.com/TopBrussels/TopTreeAnalysis/blob/CMSSW_53X/macros/BTagTReeCreator.cc
     float w_pttop = 1;
     float w_pttop_top=1;
     float w_pttop_atop=1;
		 if (dataSetName.find("TTbarJets") != string::npos)
		 {        
        for (int p=0; p<mcParticles.size(); p++)
				{        
         if (mcParticles[p]->status() != 3 ) continue;
         if (mcParticles[p]->Pt() == 0 ) continue;                 
         float pt = mcParticles[p]->Pt();

         if (mcParticles[p]->type() == 6) w_pttop_top = exp(0.156-(0.00137*pt));
         if (mcParticles[p]->type() == -6) w_pttop_atop = exp(0.156-(0.00137*pt));
        }
        w_pttop = sqrt(w_pttop_top*w_pttop_atop);
		 }
		 
		 
		 
     //cout<<"event with run = "<<event->runId()<<", lumiblock = "<<event->lumiBlockId()<<", event id = "<<event->eventId()<<endl;
		 //cout<<"    => selectedForwardJets.size() = "<<selectedForwardJets.size()<<", selectedJets.size() = "<<selectedJets.size()<<endl;
		 //cout<<"    => MuTrigged = "<<triggedSemiMu<<", ElTrigged = "<<triggedSemiEl<<", SelectednEl = "<<SelectednEl<<", SelectednMu = "<<SelectednMu<<", selectedJets[0].Pt() = "<<SelectedJetsTLV[0].Pt()<<endl; 
				 
     for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++)
     {
		        //cout<<"  jet "<<iJet<<" Pt = "<<selectedJets[iJet]->Pt()<<endl;
            SelectedJetsTLV.push_back( *selectedJets[iJet] );
            SelectedJets_bTagCSV.push_back(selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags());
						SelectedJets_partonFlavour.push_back(selectedJets[iJet]->partonFlavour());
     }
     for(unsigned int iJet=0; iJet<selectedForwardJets.size(); iJet++)
     {
			      SelectedForwardJetsTLV.push_back( *selectedForwardJets[iJet] );
						SelectedForwardJets_bTagCSV.push_back(selectedForwardJets[iJet]->btag_combinedSecondaryVertexBJetTags());
						SelectedForwardJets_partonFlavour.push_back(selectedForwardJets[iJet]->partonFlavour());
		 }
		 for(unsigned int iMuon=0; iMuon<selectedMuons.size(); iMuon++)
     {
            SelectedMuonsTLV.push_back( *selectedMuons[iMuon] );
						float relIso;
						relIso = (selectedMuons[iMuon]->chargedHadronIso()+selectedMuons[iMuon]->neutralHadronIso()+selectedMuons[iMuon]->photonIso())/selectedMuons[iMuon]->Pt();      		  
						SelectedMuonsRelIso.push_back(relIso);
						SelectedMuonsCharge.push_back((int) selectedMuons[iMuon]->charge());
		 }
		 for(unsigned int iElectron=0; iElectron<selectedElectrons.size(); iElectron++)
     {
            SelectedElectronsTLV.push_back( *selectedElectrons[iElectron] );
						float relIso;
						relIso = (selectedElectrons[iElectron]->chargedHadronIso()+selectedElectrons[iElectron]->neutralHadronIso()+selectedElectrons[iElectron]->photonIso())/selectedElectrons[iElectron]->Pt();      		  
						SelectedElectronsRelIso.push_back(relIso);
						SelectedElectronsCharge.push_back((int) selectedElectrons[iElectron]->charge());	
     }
	   //cout<<"  met Et = "<<(*mets[0]).Et()<<", met Pt = "<<(*mets[0]).Pt()<<endl;
		 
     //make and fill tree
     myBranch_selectedEvents = new VLQTree();		 
     myBranch_selectedEvents->setEventID( event->eventId() );
     myBranch_selectedEvents->setRunID( event->runId() );
     myBranch_selectedEvents->setLumiBlockID( event->lumiBlockId() );
     myBranch_selectedEvents->setNPV( vertex.size() );
     myBranch_selectedEvents->setNTruePU( (unsigned int) event->nTruePU() ); 
		 myBranch_selectedEvents->setMuTrigged( triggedSemiMu );
		 myBranch_selectedEvents->setElTrigged( triggedSemiEl );
		 if(eventselected_MuChannel && eventselected_ElChannel) cout<<"WARNING: event in mu AND el channel at the same time -> check selections for bugs!"<<endl;
		 myBranch_selectedEvents->setMuChannel( eventselected_MuChannel );
		 myBranch_selectedEvents->setElChannel( eventselected_ElChannel );		 
		 myBranch_selectedEvents->setSelectednLeptons( SelectednLeptons );
     myBranch_selectedEvents->setSelectednMu( SelectednMu );
     myBranch_selectedEvents->setSelectednEl( SelectednEl );
		 myBranch_selectedEvents->setScaleFactor( scaleFactor );
		 myBranch_selectedEvents->setMET( *mets[0] );
		 myBranch_selectedEvents->setSelectedJets( SelectedJetsTLV );
		 myBranch_selectedEvents->setSelectedJets_bTagCSV( SelectedJets_bTagCSV );
		 myBranch_selectedEvents->setSelectedJets_partonFlavour( SelectedJets_partonFlavour );
		 myBranch_selectedEvents->setSelectedForwardJets( SelectedForwardJetsTLV );
		 myBranch_selectedEvents->setSelectedForwardJets_partonFlavour( SelectedForwardJets_partonFlavour );
		 myBranch_selectedEvents->setMuons( SelectedMuonsTLV );
     myBranch_selectedEvents->setElectrons( SelectedElectronsTLV );
		 myBranch_selectedEvents->setMuonsRelIso( SelectedMuonsRelIso );
     myBranch_selectedEvents->setElectronsRelIso( SelectedElectronsRelIso );
		 myBranch_selectedEvents->setMuonsCharge( SelectedMuonsCharge );
     myBranch_selectedEvents->setElectronsCharge( SelectedElectronsCharge );
		 myBranch_selectedEvents->setTopPtReWeight( w_pttop );
		 //if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) ) myBranch_selectedEvents->setmcQuarksForMatching( mcParticlesTLV );		 
		 myVLQTree->Fill();		

		 delete myBranch_selectedEvents;     
    }
		
		if(countTightLooseLeptons && (eventpseudoselected_MuChannel || eventpseudoselected_ElChannel))
		{ 
		  float eqlumi = datasets[d]->EquivalentLumi();
		  float datasetSF = 1;
		
		  //now do same 3/4-lepton signal selections as in VLQTreeAnalyzer.cc
		  //trilepton
			cout<<"          selectedJets.size() = "<<selectedJets.size()<<endl;
			if(selectedJets.size()<2) continue;
			unsigned int nAntiBtags = 0;
			for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++)
      {
			      if(doAntiBtagging) cout<<"            selectedJets["<<iJet<<"]->btag_combinedSecondaryVertexBJetTags() = "<<selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags()<<endl;
			      if((0 <= selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags()) && (selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags() < antibtagWP))
		           nAntiBtags++;
      }
			if(doAntiBtagging) cout<<"          nAntiBtags = "<<nAntiBtags<<endl;
			if((doAntiBtagging && (nAntiBtags != selectedJets.size())) ){ cout<<"... skipping"<<endl; continue;} //all jets should be anti-b-tagged...
			cout<<"          selectedJets[0]->Pt() = "<<selectedJets[0]->Pt()<<", selectedJets[1]->Pt() = "<<selectedJets[1]->Pt()<<endl;
			if((selectedJets[0]->Pt() < LeadingJetPtCut) || (selectedJets[1]->Pt() < SubLeadingJetPtCut)){ cout<<"... skipping"<<endl; continue;}
			cout<<"          MET = "<<(*mets[0]).Pt()<<endl;
			if((doMETcut && (*mets[0]).Pt()<METcut)){ cout<<"... skipping"<<endl; continue;}
			
			if(eventpseudoselected_MuChannel && eventpseudoselected_ElChannel) cout<<"WARNING: event in mu AND el channel at the same time -> check pseudoselections for bugs!"<<endl;
			
			//mu channel
			//if(triggedSemiMu)
			if(eventpseudoselected_MuChannel)
			{			
			 datasetSF = LuminosityMu / eqlumi; //needed when running on MC, scale event counts with data lumi...
			 if((dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0 || (dataSetName.find("InvIso_Mu") != string::npos)) && datasetSF != 1)
			    cout<<"WARNING: data SF not 1"<<endl;
				
			 cout<<"      ==> surviving antibtags and jet cuts!"<<endl;
			 if(SelectednMu == 1 && SelectednEl == 1 && SelectednLooseNotTightEl == 1)
			 {
			  //Tmu,Tel,Lel
				massleptons = (*selectedElectrons[0] + *selectedOnlyLooseElectrons_FL[0]).M();
				if((selectedElectrons[0]->charge() != selectedOnlyLooseElectrons_FL[0]->charge()) && (fabs(Zmass - massleptons) < Zmasswindow) && (((*selectedElectrons[0] + *selectedOnlyLooseElectrons_FL[0]).Pt() > Zptcut) || !doZptcut))
				{
						N_TmuTelLel_mu = N_TmuTelLel_mu + datasetSF;
						uncN_TmuTelLel_mu = sqrt(pow(uncN_TmuTelLel_mu,2) + pow(datasetSF,2));
						cout<<"          N_TmuTelLel_mu = "<<N_TmuTelLel_mu<<endl;
				}		
			 }
			 else if(SelectednMu == 2 && SelectednEl == 0 && SelectednLooseNotTightEl == 1)
			 {
			  //Tmu,Tmu,Lel
				massleptons = (*selectedMuons[0] + *selectedMuons[1]).M();
				if((selectedMuons[0]->charge() != selectedMuons[1]->charge()) && (fabs(Zmass - massleptons) < Zmasswindow) && (((*selectedMuons[0] + *selectedMuons[1]).Pt() > Zptcut) || !doZptcut))
				{
						N_TmuTmuLel_mu = N_TmuTmuLel_mu + datasetSF;
						uncN_TmuTmuLel_mu = sqrt(pow(uncN_TmuTmuLel_mu,2) + pow(datasetSF,2));
						cout<<"          N_TmuTmuLel_mu = "<<N_TmuTmuLel_mu<<endl;
				}		
			 }
		   else if(SelectednMu == 1 && SelectednEl == 1 && SelectednLooseNotTightMu == 1)
			 {
			  //Tmu,Tel,Lmu
				massleptons = (*selectedMuons[0] + *selectedOnlyLooseMuons_FL[0]).M();
				if((selectedMuons[0]->charge() != selectedOnlyLooseMuons_FL[0]->charge()) && (fabs(Zmass - massleptons) < Zmasswindow) && (((*selectedMuons[0] + *selectedOnlyLooseMuons_FL[0]).Pt() > Zptcut) || !doZptcut))
				{
						N_TmuTelLmu_mu = N_TmuTelLmu_mu + datasetSF;
						uncN_TmuTelLmu_mu = sqrt(pow(uncN_TmuTelLmu_mu,2) + pow(datasetSF,2));
						cout<<"         N_TmuTelLmu_mu  = "<<N_TmuTelLmu_mu<<endl;
				}		
			 }
		   else if(SelectednMu == 2 && SelectednEl == 0 && SelectednLooseNotTightMu == 1)
			 {
			            //Tmu,Tmu,Lmu
			            TLorentzVector Zcandidate; 
									if(selectedMuons[0]->charge() != selectedMuons[1]->charge())
									{
									  TLorentzVector Zcandidate_opt1 = ((TLorentzVector) *selectedMuons[0] + (TLorentzVector) *selectedMuons[1]);
										TLorentzVector Zcandidate_opt2;
										if(selectedMuons[0]->charge() != selectedOnlyLooseMuons_FL[0]->charge()) Zcandidate_opt2 = ((TLorentzVector) *selectedMuons[0] + (TLorentzVector) *selectedOnlyLooseMuons_FL[0]);
										else if(selectedMuons[1]->charge() != selectedOnlyLooseMuons_FL[0]->charge()) Zcandidate_opt2 = ((TLorentzVector) *selectedMuons[1] + (TLorentzVector) *selectedOnlyLooseMuons_FL[0]);
									  float massleptons_opt1 = Zcandidate_opt1.M(); 
									  float massleptons_opt2 = Zcandidate_opt2.M(); 
										if(massleptons_opt1 < massleptons_opt2)
										{
										   massleptons = massleptons_opt1;
											 Zcandidate = Zcandidate_opt1;
										}
										else
										{
										   massleptons = massleptons_opt2;
											 Zcandidate = Zcandidate_opt2;
										}
									}
									else if((selectedMuons[0]->charge() == selectedMuons[1]->charge()) && (selectedMuons[0]->charge() != selectedOnlyLooseMuons_FL[0]->charge()))
									{
									  TLorentzVector Zcandidate_opt1 = ((TLorentzVector) *selectedMuons[0] + (TLorentzVector) *selectedOnlyLooseMuons_FL[0]);
										TLorentzVector Zcandidate_opt2 = ((TLorentzVector) *selectedMuons[1] + (TLorentzVector) *selectedOnlyLooseMuons_FL[0]);
									  float massleptons_opt1 = Zcandidate_opt1.M(); 
									  float massleptons_opt2 = Zcandidate_opt2.M(); 
										if(massleptons_opt1 < massleptons_opt2)
										{
										   massleptons = massleptons_opt1;
											 Zcandidate = Zcandidate_opt1;
										}
										else
										{
										   massleptons = massleptons_opt2;
											 Zcandidate = Zcandidate_opt2;
										}
									}	
									if((fabs(Zmass - massleptons) < Zmasswindow) && ((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut))
									{
                    N_TmuTmuLmu_mu = N_TmuTmuLmu_mu + datasetSF;
										uncN_TmuTmuLmu_mu = sqrt(pow(uncN_TmuTmuLmu_mu,2) + pow(datasetSF,2));
										cout<<"         N_TmuTmuLmu_mu  = "<<N_TmuTmuLmu_mu<<endl;
									}									
			 }
			 
			 //TBC: 4-lepton
			 
			}
			
			//el channel
			//if(triggedSemiEl) //first I did if(triggedSemiEl) instead of else if, but that's not consistent with treeanalyzer (and rest of treecreator) -> update: changed all back to if instead of else if
			else if(eventpseudoselected_MuChannel)
			{
			 datasetSF = LuminosityEl / eqlumi; //needed when running on MC, scale event counts it with data lumi...
			 if((dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0 || (dataSetName.find("InvIso_El") != string::npos)) && datasetSF != 1)
			    cout<<"WARNING: data SF not 1"<<endl;
					
			 cout<<"      ==> surviving antibtags and jet cuts!"<<endl;
			 if(SelectednEl == 2 && SelectednMu == 0 && SelectednLooseNotTightEl == 1)
			 {
			            //Tel,Tel,Lel
			            TLorentzVector Zcandidate; 
									if(selectedElectrons[0]->charge() != selectedElectrons[1]->charge())
									{
									  TLorentzVector Zcandidate_opt1 = ((TLorentzVector) *selectedElectrons[0] + (TLorentzVector) *selectedElectrons[1]);
										TLorentzVector Zcandidate_opt2;
										if(selectedElectrons[0]->charge() != selectedOnlyLooseElectrons_FL[0]->charge()) Zcandidate_opt2 = ((TLorentzVector) *selectedElectrons[0] + (TLorentzVector) *selectedOnlyLooseElectrons_FL[0]);
										else if(selectedElectrons[1]->charge() != selectedOnlyLooseElectrons_FL[0]->charge()) Zcandidate_opt2 = ((TLorentzVector) *selectedElectrons[1] + (TLorentzVector) *selectedOnlyLooseElectrons_FL[0]);
									  float massleptons_opt1 = Zcandidate_opt1.M(); 
									  float massleptons_opt2 = Zcandidate_opt2.M(); 
										if(massleptons_opt1 < massleptons_opt2)
										{
										   massleptons = massleptons_opt1;
											 Zcandidate = Zcandidate_opt1;
										}
										else
										{
										   massleptons = massleptons_opt2;
											 Zcandidate = Zcandidate_opt2;
										}
									}
									else if((selectedElectrons[0]->charge() == selectedElectrons[1]->charge()) && (selectedElectrons[0]->charge() != selectedOnlyLooseElectrons_FL[0]->charge()))
									{
									  TLorentzVector Zcandidate_opt1 = ((TLorentzVector) *selectedElectrons[0] + (TLorentzVector) *selectedOnlyLooseElectrons_FL[0]);
										TLorentzVector Zcandidate_opt2 = ((TLorentzVector) *selectedElectrons[1] + (TLorentzVector) *selectedOnlyLooseElectrons_FL[0]);
									  float massleptons_opt1 = Zcandidate_opt1.M(); 
									  float massleptons_opt2 = Zcandidate_opt2.M(); 
										if(massleptons_opt1 < massleptons_opt2)
										{
										   massleptons = massleptons_opt1;
											 Zcandidate = Zcandidate_opt1;
										}
										else
										{
										   massleptons = massleptons_opt2;
											 Zcandidate = Zcandidate_opt2;
										}
									}	
									if((fabs(Zmass - massleptons) < Zmasswindow) && ((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut))
									{
                    N_TelTelLel_el = N_TelTelLel_el + datasetSF;
										uncN_TelTelLel_el = sqrt(pow(uncN_TelTelLel_el,2) + pow(datasetSF,2));
										cout<<"          N_TelTelLel_el = "<<N_TelTelLel_el<<endl;
									}									
			 }
			 else if(SelectednEl == 1 && SelectednMu == 1 && SelectednLooseNotTightEl == 1)
			 {
			  //Tel,Tmu,Lel
				massleptons = (*selectedElectrons[0] + *selectedOnlyLooseElectrons_FL[0]).M();
				if((selectedElectrons[0]->charge() != selectedOnlyLooseElectrons_FL[0]->charge()) && (fabs(Zmass - massleptons) < Zmasswindow) && (((*selectedElectrons[0] + *selectedOnlyLooseElectrons_FL[0]).Pt() > Zptcut) || !doZptcut))
				{
						N_TelTmuLel_el = N_TelTmuLel_el + datasetSF;
						uncN_TelTmuLel_el = sqrt(pow(uncN_TelTmuLel_el,2) + pow(datasetSF,2));
						cout<<"          N_TelTmuLel_el = "<<N_TelTmuLel_el<<endl;
				}	
			 }
			 else if(SelectednEl == 2 && SelectednMu == 0 && SelectednLooseNotTightMu == 1)
			 {
			  //Tel,Tel,Lmu
				massleptons = (*selectedElectrons[0] + *selectedElectrons[1]).M();
				if((selectedElectrons[0]->charge() != selectedElectrons[1]->charge()) && (fabs(Zmass - massleptons) < Zmasswindow) && (((*selectedElectrons[0] + *selectedElectrons[1]).Pt() > Zptcut) || !doZptcut))
				{
						N_TelTelLmu_el = N_TelTelLmu_el + datasetSF;
						uncN_TelTelLmu_el = sqrt(pow(uncN_TelTelLmu_el,2) + pow(datasetSF,2));
						cout<<"          N_TelTelLmu_el = "<<N_TelTelLmu_el<<endl;
				}
			 }	
			 else if(SelectednEl == 1 && SelectednMu == 1 && SelectednLooseNotTightMu == 1)
			 {
			  //Tel,Tmu,Lmu
				massleptons = (*selectedMuons[0] + *selectedOnlyLooseMuons_FL[0]).M();
				if((selectedMuons[0]->charge() != selectedOnlyLooseMuons_FL[0]->charge()) && (fabs(Zmass - massleptons) < Zmasswindow) && (((*selectedMuons[0] + *selectedOnlyLooseMuons_FL[0]).Pt() > Zptcut) || !doZptcut))
				{
						N_TelTmuLmu_el = N_TelTmuLmu_el + datasetSF;
						uncN_TelTmuLmu_el = sqrt(pow(uncN_TelTmuLmu_el,2) + pow(datasetSF,2));
						cout<<"          N_TelTmuLmu_el = "<<N_TelTmuLmu_el<<endl;
				}		
			 }
			 
			 //TBC: 4-lepton
			 
			}
			
		
		
		
		
		
		
		}
		
	}			//loop on events

   cout<<endl;

   //cout << "+> " << nSelectedMu << " mu+jets events where selected"<< endl;
    //if(!countTightLooseLeptons)
	  //{
     treeFile->cd();
		 TTree *configTreeFile = new TTree("configTreeFile","configuration Tree in tree File");
     TClonesArray* tcdatasettreefile = new TClonesArray("Dataset",1);
     configTreeFile->Branch("Dataset","TClonesArray",&tcdatasettreefile);
     TClonesArray* tcAnaEnvtreeFile = new TClonesArray("AnalysisEnvironment",1);
     configTreeFile->Branch("AnaEnv","TClonesArray",&tcAnaEnvtreeFile);
     new ((*tcAnaEnvtreeFile)[0]) AnalysisEnvironment(anaEnv);
     new ((*tcdatasettreefile)[0]) Dataset(*datasets[d]);    
     configTreeFile->Fill();
     configTreeFile->Write();
		 myVLQTree->Write();
		 treeFile->Close();
		 delete treeFile;
	  //}
		 
    //////////////
    // CLEANING //
    //////////////

    if (jecUnc) delete jecUnc;
    if (jetTools) delete jetTools;
  
    //important: free memory
    treeLoader.UnLoadDataset();
    
   }				//loop on datasets
	 
	 
	ofstream myfile1;
  if(countTightLooseLeptons && systematic==0)
	{		
		string myRockingFile1 = Treespath+"/VLQ_FakeLepton_Events"+suffix+".txt";
		myfile1.open(myRockingFile1.c_str());
		myfile1 << "\n";

			myfile1 << "THIS IS FOR THE MUON CHANNEL PART" << "\n";
			myfile1 << "# mu+el+el events with loose-not-tight electron: " << N_TmuTelLel_mu << " +- " << uncN_TmuTelLel_mu  << "\n";
			myfile1 << "# mu+mu+el events with loose-not-tight electron: " << N_TmuTmuLel_mu << " +- " << uncN_TmuTmuLel_mu  << "\n";
			myfile1 << "# mu+mu+el events with loose-not-tight muon: " << N_TmuTelLmu_mu << " +- " << uncN_TmuTelLmu_mu  << "\n";
			myfile1 << "# mu+mu+mu events with loose-not-tight muon: " << N_TmuTmuLmu_mu << " +- " << uncN_TmuTmuLmu_mu  << "\n";
			myfile1 << "\n";
			
			float Unc_meLe_mu = pow((Eff_TL_e/(1-Eff_TL_e))*uncN_TmuTelLel_mu,2)+pow(N_TmuTelLel_mu*(1/pow(1-Eff_TL_e,2))*Eff_TL_e_unc,2);
			myfile1 << "# predicted mu+el+Lel events: " << N_TmuTelLel_mu*Eff_TL_e/(1-Eff_TL_e) << " +- " << sqrt(Unc_meLe_mu) << "\n";
      float Unc_mmLe_mu = pow((Eff_TL_e/(1-Eff_TL_e))*uncN_TmuTmuLel_mu,2)+pow(N_TmuTmuLel_mu*(1/pow(1-Eff_TL_e,2))*Eff_TL_e_unc,2);
			myfile1 << "# predicted mu+mu+Lel events: " << N_TmuTmuLel_mu*Eff_TL_e/(1-Eff_TL_e) << " +- " << sqrt(Unc_mmLe_mu) << "\n";
      float Unc_meLm_mu = pow((Eff_TL_m/(1-Eff_TL_m))*uncN_TmuTelLmu_mu,2)+pow(N_TmuTelLmu_mu*(1/pow(1-Eff_TL_m,2))*Eff_TL_m_unc,2);
			myfile1 << "# predicted mu+el+Lmu events: " << N_TmuTelLmu_mu*Eff_TL_m/(1-Eff_TL_m) << " +- " << sqrt(Unc_meLm_mu) << "\n";
      float Unc_mmLm_mu = pow((Eff_TL_m/(1-Eff_TL_m))*uncN_TmuTmuLmu_mu,2)+pow(N_TmuTmuLmu_mu*(1/pow(1-Eff_TL_m,2))*Eff_TL_m_unc,2);
			myfile1 << "# predicted mu+mu+Lmu events: " << N_TmuTmuLmu_mu*Eff_TL_m/(1-Eff_TL_m) << " +- " << sqrt(Unc_mmLm_mu) << "\n";

      myfile1 << "\n";
			myfile1 << "THIS IS FOR THE ELECTRON CHANNEL PART" << "\n";
			myfile1 << "# el+mu+mu events with loose-not-tight muon: " << N_TelTmuLmu_el << " +- " << uncN_TelTmuLmu_el  << "\n";
			myfile1 << "# el+el+mu events with loose-not-tight muon: " << N_TelTelLmu_el << " +- " << uncN_TelTelLmu_el  << "\n";
			myfile1 << "# el+el+mu events with loose-not-tight electron: " << N_TelTmuLel_el << " +- " << uncN_TelTmuLel_el  << "\n";
			myfile1 << "# el+el+el events with loose-not-tight electron: " << N_TelTelLel_el << " +- " << uncN_TelTelLel_el  << "\n";
			myfile1 << "\n";
			
			float Unc_emLm_el = pow((Eff_TL_m/(1-Eff_TL_m))*uncN_TelTmuLmu_el,2)+pow(N_TelTmuLmu_el*(1/pow(1-Eff_TL_m,2))*Eff_TL_m_unc,2);
			myfile1 << "# predicted el+mu+Lmu events: " << N_TelTmuLmu_el*Eff_TL_m/(1-Eff_TL_m) << " +- " << sqrt(Unc_emLm_el) << "\n";
      float Unc_eeLm_el = pow((Eff_TL_m/(1-Eff_TL_m))*uncN_TelTelLmu_el,2)+pow(N_TelTelLmu_el*(1/pow(1-Eff_TL_m,2))*Eff_TL_m_unc,2);
			myfile1 << "# predicted el+el+Lmu events: " << N_TelTelLmu_el*Eff_TL_m/(1-Eff_TL_m) << " +- " << sqrt(Unc_eeLm_el) << "\n";
      float Unc_emLe_el = pow((Eff_TL_e/(1-Eff_TL_e))*uncN_TelTmuLel_el,2)+pow(N_TelTmuLel_el *(1/pow(1-Eff_TL_e,2))*Eff_TL_e_unc,2);
			myfile1 << "# predicted el+mu+Lel events: " << N_TelTmuLel_el*Eff_TL_e/(1-Eff_TL_e) << " +- " << sqrt(Unc_emLe_el) << "\n";
      float Unc_eeLe_el = pow((Eff_TL_e/(1-Eff_TL_e))*uncN_TelTelLel_el,2)+pow(N_TelTelLel_el*(1/pow(1-Eff_TL_e,2))*Eff_TL_e_unc,2);
			myfile1 << "# predicted el+el+Lel events: " << N_TelTelLel_el*Eff_TL_e/(1-Eff_TL_e) << " +- " << sqrt(Unc_eeLe_el) << "\n";

		myfile1 << "\n";
		myfile1.close();
	}
	
	
	
    //Once everything is filled ...
    if (verbose > 0)
      cout << " We ran over all the data ;-)" << endl;


  
    //Selection tables
		/*selecTableSemiMu.TableCalculator(false, true, true, true, true);
    string selectiontableMu = "SelectionTable_VLQ_SEMIMu.tex";
    selecTableSemiMu.Write(selectiontableMu.c_str());
    selecTableSemiEl.TableCalculator(false, true, true, true, true);
    string selectiontableEl = "SelectionTable_VLQ_SEMIEL.tex";
    selecTableSemiEl.Write(selectiontableEl.c_str());*/
    
		
		if (runSpecificSample) return 0;
		
		
    // Do some special things with certain plots (normalize, BayesDivide, ... )
    if (verbose > 0)
      cout << "Treating the special plots." << endl;
    
    
    delete fout;
    delete tcdatasets;
    delete tcAnaEnv;
    delete configTree;
        
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "           hasn't crashed yet ;-)           " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}
