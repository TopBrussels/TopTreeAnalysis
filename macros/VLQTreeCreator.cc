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
  string Treespath = "VLQTrees_Summer12_PBS_ExclusiveTTbarSamples";
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


  /////////////////////
  // Configuration
  /////////////////////
  string rootFileName("VLQTreeCreator_Plots_"+postfix+".root"); //some output, maybe dummy
	
  //xml file
  string xmlFileName ="../config/myVLQConfig_forexclusivettbar.xml";	

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

  float LuminosityMu = oldLuminosity;
  float LuminosityEl = oldLuminosity;

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

	
  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  bool changed_anaEnvCollections = false;
	
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++) {

    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
		
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
	  treeFile = new TFile(TreeFileName.c_str(),"RECREATE");
	  myVLQTree = new TTree("myVLQTree","Tree containing the VLQ analysis info");    
		myVLQTree->Branch("VLQBranch_selectedEvents","VLQTree",&myBranch_selectedEvents);
		
    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    //for (unsigned int ievt = 0; ievt < 20000; ievt++)
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
	         cout<<"File changed!!! => iFile = "<<iFile<<endl;
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
            exit(-1);
          }

	} else {
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
	      exit(-1);
	  }
   }
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
      selection.cutHFHadronEnergyFraction();
			nonstandard_selection.cutHFHadronEnergyFraction();
			selection.setJetCuts(30,2.4,0.01,1.,0.98,0.3,0.1); //was 2.5, but adapted to 2.4 for b-tagging... correct?
			nonstandard_selection.setJetCuts(30.,4.7,0.01,1.,0.98,0.3,0.1); //only difference: larger eta acceptance 
      selection.setMuonCuts(20,2.1,0.12,0.2,0.3,1,0.5,5,0); //Pt 20 GeV! later on, higher for leading lepton at least
      selection.setElectronCuts(20,2.5,0.1,0.02,0.5,0.3,0); //was pt 32... //Pt 20 GeV! later on, higher for leading lepton at least
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setLooseElectronCuts(15,2.5,0.15,0.); //Pt 15 GeV! first I had pt 20
			float LeadingLeptonPtCut = 30.;
			
      bool triggedSemiMu = false;
      bool triggedSemiEl = false;

      //to be checked: here there is/was no mention of what to do in the InvIso case, but that was dangerous... it is still correct, but looks like a lucky accident
      if( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) )
        triggedSemiMu = treeLoader.EventTrigged (itriggerSemiMu);
      if( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) )
        triggedSemiEl = treeLoader.EventTrigged (itriggerSemiEl);

      bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);


      vector<TRootJet*> selectedJets, selectedForwardJets, selectedJetsLargeEtaRange, selectedJetsNoMu, selectedJetsNoEl, selectedJetsLargeEtaRangeNoMu, selectedJetsLargeEtaRangeNoEl;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
			vector<TRootMuon*> selectedLooseMuons = selection.GetSelectedLooseMuons(10,2.5,0.2); //to veto
      //vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons(10,2.5,0.2);
			vector<TRootElectron*> selectedLooseElectrons = selection.GetSelectedLooseElectrons(20,2.5,0.15); //to veto
      //vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(20,2.5,0.15);
      //vector<TRootElectron*> vetoElectronsSemiEl = selection.GetSelectedLooseElectrons(20,2.5,0.15);

      if( dataSetName.find("InvIso") != string::npos )  { // event selection for special Data TopTrees for ARI QCD
        vector<TRootMuon*> overlapMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0]);
        vector<TRootElectron*> overlapElectrons = selection.GetSelectedElectronsInvIso(0.2);
        selectedJetsNoMu = selection.GetSelectedJets(overlapMuons,true);
        selectedJetsNoEl = selection.GetSelectedJets(overlapElectrons,true);

	      /*if (selectedJetsNoMu.size() >= 4) {
	         //cout << "ol" << endl;
	         if (selectedJetsNoMu[0]->Pt() < 45) selectedJetsNoMu.clear();
	         if (selectedJetsNoMu[1]->Pt() < 45) selectedJetsNoMu.clear();
	         if (selectedJetsNoMu[2]->Pt() < 40) selectedJetsNoMu.clear();
	         if (selectedJetsNoMu[3]->Pt() < 40) selectedJetsNoMu.clear();
	      }

	      if (selectedJetsNoEl.size() >= 4) {
	         //cout << "ol" << endl;
	         if (selectedJetsNoEl[0]->Pt() < 45) selectedJetsNoEl.clear();
	         if (selectedJetsNoEl[1]->Pt() < 45) selectedJetsNoEl.clear();
	         if (selectedJetsNoEl[2]->Pt() < 40) selectedJetsNoEl.clear();
	         if (selectedJetsNoEl[3]->Pt() < 40) selectedJetsNoEl.clear();
	      }*/

	      //selectedJetsNoMu = selection.GetSelectedJets(true);
        //selectedJetsNoEl = selection.GetSelectedJets(true);

	      selectedMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0], selectedJetsNoMu);
        selectedElectrons = selection.GetSelectedElectronsInvIso(0.2,selectedJetsNoEl); // Also mvaId < 0.5
				
				//(z)
				//something with selectedJetsLargeEtaRangeNoMu etc or not??
				selectedJetsLargeEtaRangeNoMu = nonstandard_selection.GetSelectedJets(overlapMuons,true);
        selectedJetsLargeEtaRangeNoEl = nonstandard_selection.GetSelectedJets(overlapElectrons,true);
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
			
			bool eventselected = false;
    
      // selection with mu
      if (triggedSemiMu) {
				if (isGoodPV) {
	  			if(dataSetName.find("InvIso") != string::npos)
					{
					   selectedJets = selectedJetsNoMu;
					   selectedJetsLargeEtaRange = selectedJetsLargeEtaRangeNoMu;
					}
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
									 
						if(SelectednMu == 0)
						  continue;
						if(selectedMuons[0]->Pt() < LeadingLeptonPtCut)
						  continue;
						
						if (SelectednMu == 1 && (selectedLooseMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 ))) // if InvertedIso, selected muon not part ofselected LooseMuons vector! (to be checked)
						{	  				  
							
							//if( vetoMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && vetoMuons.size() == 0 ) ) { // if InvertedIso, selected muon not part of vetoMuons vector!
	  				  //  if (vetoElectronsSemiMu.size() == 0) {
							//     if (selectedJets.size() >= 4) {
							//        eventselectedSemiMu = true;
						  //    }
	            //  }
	            //}
							if(SelectednEl == 0 && (selectedLooseElectrons.size() == 0) )
							{
							     eventselected = true;
							}
							else if(SelectednEl == 1 && (selectedLooseElectrons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{
							     eventselected = true;
							}
							else if(SelectednEl == 2 && (selectedLooseElectrons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{							    
							     eventselected = true;
							}
							else if(SelectednEl == 3 && (selectedLooseElectrons.size() == 3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{
							     eventselected = true;								 
							}
							else if(SelectednEl >= 4 && (selectedLooseElectrons.size() >=  4 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )) )
							{
							     cout<<"FOUND >= 5-lepton event!"<<endl;		
							}								
	         }
					 else if (SelectednMu == 2 && (selectedLooseMuons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 ))) 
					 {
					    if(SelectednEl == 0 && selectedLooseElectrons.size() == 0)
							{
							     eventselected = true;
							}
							else if(SelectednEl == 1 && (selectedLooseElectrons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{
							     eventselected = true;
							}
							else if(SelectednEl == 2 && (selectedLooseElectrons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{
							     eventselected = true;
							}
							else if(SelectednEl >= 3 && (selectedLooseElectrons.size() >=  3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							     cout<<"FOUND >= 5-lepton event!"<<endl;		 
					 }
					 else if (SelectednMu == 3 && (selectedLooseMuons.size() == 3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 ))) 
					 {
					    if(SelectednEl == 0 && selectedLooseElectrons.size() == 0)
							{
							     eventselected = true;
							}
							else if(SelectednEl == 1 && (selectedLooseElectrons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							{
							     eventselected = true;
							}
							else if(SelectednEl >= 2 && (selectedLooseElectrons.size() >=  2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							  cout<<"FOUND >= 5-lepton event!"<<endl;
					 
					 }
					 else if (SelectednMu == 4 && (selectedLooseMuons.size() == 4 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 ))) 
					 {
					    if(SelectednEl == 0 && selectedLooseElectrons.size() == 0)
							{
							     eventselected = true;
							}
							else if(SelectednEl >= 1 && (selectedLooseElectrons.size() >=  1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							  cout<<"FOUND >= 5-lepton event!"<<endl;
					 
					 }
					 else if(SelectednMu >= 5 && (selectedLooseMuons.size() >=  5 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
							  cout<<"FOUND >= 5-lepton event!"<<endl;
					 
				 } //nJets
	      } //good PV
      }
      else if( triggedSemiEl) {
        if(dataSetName.find("InvIso") != string::npos)
				{
				    selectedJets = selectedJetsNoEl;
						selectedJetsLargeEtaRange = selectedJetsLargeEtaRangeNoEl;
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
							     eventselected = true;
							  }
								else if(SelectednMu == 1 && (selectedLooseMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								{
								   if(selectedMuons[0]->Pt() < LeadingLeptonPtCut) //to avoid double counting
									 {
									   eventselected = true;
									 }
								}
								else if(SelectednMu == 2 && (selectedLooseMuons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								{
								   if((selectedMuons[0]->Pt() < LeadingLeptonPtCut) && (selectedMuons[1]->Pt() < LeadingLeptonPtCut)) //to avoid double counting
									 {
									   eventselected = true;
									 }
								}
								else if(SelectednMu == 3 && (selectedLooseMuons.size() == 3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								{
								   if((selectedMuons[0]->Pt() < LeadingLeptonPtCut) && (selectedMuons[1]->Pt() < LeadingLeptonPtCut) && (selectedMuons[2]->Pt() < LeadingLeptonPtCut)) //to avoid double counting
									 {
									    eventselected = true;
									 }
								}
								else if(SelectednMu >= 4 && (selectedLooseMuons.size() >=  4 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
							     cout<<"FOUND >= 5-lepton event!"<<endl;
							  																
	          }
						else if(SelectednEl == 2  && (selectedLooseElectrons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 ))) {
							
								if(SelectednMu == 0 && selectedLooseMuons.size() ==  0)
							  {
							     eventselected = true;
							  }
								else if(SelectednMu == 1 && (selectedLooseMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )) )
								{
								   if(selectedMuons[0]->Pt() < LeadingLeptonPtCut) //to avoid double counting
									 {
									    eventselected = true;
									 }
								}
								else if(SelectednMu == 2 && (selectedLooseMuons.size() == 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								{
								   if((selectedMuons[0]->Pt() < LeadingLeptonPtCut) && (selectedMuons[1]->Pt() < LeadingLeptonPtCut)) //to avoid double counting
									 {
									    eventselected = true;
									 }
								}
								else if(SelectednMu >= 3 && (selectedLooseMuons.size() >=  3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
							     cout<<"FOUND >= 5-lepton event!"<<endl;		  																
	          }
						else if(SelectednEl == 3  && (selectedLooseElectrons.size() == 3 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 ))) {
							
								if(SelectednMu == 0 && selectedLooseMuons.size() ==  0)
							  {
							     eventselected = true;
							  }
								else if(SelectednMu == 1 && (selectedLooseMuons.size() == 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								{
								   if(selectedMuons[0]->Pt() < LeadingLeptonPtCut) //to avoid double counting
									 {
									    eventselected = true;
									 }
								}
								else if(SelectednMu >= 2 && (selectedLooseMuons.size() >= 2 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
							     cout<<"FOUND >= 5-lepton event!"<<endl;						  																
	          }
						else if(SelectednEl == 4  && selectedLooseElectrons.size() == 4) {
							
								if(SelectednMu == 0 && selectedLooseMuons.size() ==  0)
							  {
							     eventselected = true;
							  }
								else if(SelectednMu >= 1 && (selectedLooseMuons.size() >= 1 || ( dataSetName.find("InvIso") != string::npos && selectedLooseMuons.size() == 0 )))
								   cout<<"FOUND >= 5-lepton event!"<<endl;							  																
	          }
						else if(SelectednEl >= 5 && (selectedLooseElectrons.size() >= 5 || ( dataSetName.find("InvIso") != string::npos && selectedLooseElectrons.size() == 0 )))
							  cout<<"FOUND >= 5-lepton event!"<<endl;
												
				 } //nJets
       } //good PV
     }
		 
		 SelectednLeptons = SelectednMu + SelectednEl;
     
		 if (!eventselected) continue;
		 
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
     for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++)
     {
		        //cout<<"jet "<<iJet<<" Pt = "<<selectedJets[iJet]->Pt()<<endl;
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
	
     //make and fill tree
     myBranch_selectedEvents = new VLQTree();		 
     myBranch_selectedEvents->setEventID( event->eventId() );
     myBranch_selectedEvents->setRunID( event->runId() );
     myBranch_selectedEvents->setLumiBlockID( event->lumiBlockId() );
     myBranch_selectedEvents->setNPV( vertex.size() );
     myBranch_selectedEvents->setNTruePU( (unsigned int) event->nTruePU() ); 
		 myBranch_selectedEvents->setMuTrigged( triggedSemiMu );
		 myBranch_selectedEvents->setElTrigged( triggedSemiEl );
		 myBranch_selectedEvents->setSelectednLeptons( SelectednLeptons );
     myBranch_selectedEvents->setSelectednMu( SelectednMu );
     myBranch_selectedEvents->setSelectednEl( SelectednEl );
		 myBranch_selectedEvents->setEventWeight( scaleFactor );
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
		 //if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) ) myBranch_selectedEvents->setmcQuarksForMatching( mcParticlesTLV );		 
		 myVLQTree->Fill();		 
		 		 
		 delete myBranch_selectedEvents;     
    }			//loop on events

    cout<<endl;

    //cout << "+> " << nSelectedMu << " mu+jets events where selected"<< endl;

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
		 
    //////////////
    // CLEANING //
    //////////////

    if (jecUnc) delete jecUnc;
    if (jetTools) delete jetTools;
  
    //important: free memory
    treeLoader.UnLoadDataset();
    
    }				//loop on datasets

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
