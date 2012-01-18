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
#include "../Selection/interface/SelectionTable.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../MCInformation/interface/MCWeighter.h"
#include "../Selection/interface/ElectronPlotter.h"
#include "../Selection/interface/MuonPlotter.h"
#include "../Selection/interface/JetPlotter.h"
#include "../Selection/interface/VertexPlotter.h"
#include "../Tools/interface/JetTools.h"
#include "../MCInformation/interface/ResolutionFit.h"
#include "../MCInformation/interface/JetPartonMatching.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../BtagEffAnalysis/interface/TRootNTuple.h"
#include "Style.C"
#include "../MCInformation/interface/LumiReWeighting.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"

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
    if (atoi(argv[1]) == 0 || atoi(argv[1]) == 1 || atoi(argv[1]) == 2 || atoi(argv[1]) == 3 || atoi(argv[1]) == 4  || atoi(argv[1]) == 5 || atoi(argv[1]) == 6  || atoi(argv[1]) == 7 || atoi(argv[1]) == 8 || atoi(argv[1]) == 9 || atoi(argv[1]) == 10 || atoi(argv[1]) == 11 || atoi(argv[1]) == 12)
      systematic=atoi(argv[1]);
  
  bool runSpecificSample=false;
  
  int nSpecSample=-19;

  if (argc > 2) {

    runSpecificSample=true;

    nSpecSample=atoi(argv[2]);

  }

  string postfix = ""; // to relabel the names of systematic samples

  if (systematic == 1)
    postfix="_JESMinus";
  if (systematic == 2)
    postfix="_JESPlus";
  if (systematic == 3)
    postfix="_JESMinusAfterSel";
  if (systematic == 4)
    postfix="_JESPlusAfterSel";
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

  cout << "systematic: " << systematic << " -> Postfix = " << postfix << endl;

  // return 0;

  bool doPF2PAT = false;

  bool useMassesAndResolutions = false;

  //useMassesAndResolutions = true;

  if (systematic != 0 || runSpecificSample)
    useMassesAndResolutions = true;

  double btagcut=3;

  clock_t start = clock();

  cout << "********************************************************" << endl;
  cout << " Beginning of the program for creating the BTag Trees ! " << endl;
  cout << "********************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

  //xml file
  string xmlFileName ="../config/myBTAGconfig.xml";
  //string xmlFileName ="../config/myBTAGconfig_fall11.xml";
  
  if (systematic == 8)
    xmlFileName ="../config/myBTAGconfig_ttjetsyst.xml";
  else if (systematic == 9)
    xmlFileName ="../config/myBTAGconfig_wjetsyst.xml";
  
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

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;
  bool runOnlyOverDataForJEC = false;
  bool isSemiMu = false;
  for (unsigned int d = 0; d < datasets.size (); d++)
    {
      //cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
      if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
      string dataSetName = datasets[d]->Name();
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	{
	  Luminosity = datasets[d]->EquivalentLumi();
	  if(datasets.size() > 1) runOnlyOverDataForJEC = true;
	  break;
	}
      if(dataSetName.find("TTbarJets_SemiMu") == 0) isSemiMu = true;
    }
  if(Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
	
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
  
  ////////////////////////////////////
  /// MultiSamplePlot
  ////////////////////////////////////

  map<string,MultiSamplePlot*> MSPlot;

  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTable;
  CutsSelecTable.push_back(string("initial"));
  CutsSelecTable.push_back(string("preselected"));
  CutsSelecTable.push_back(string("trigged"));
  CutsSelecTable.push_back(string("Good PV"));
  CutsSelecTable.push_back(string("1 selected muon"));
  CutsSelecTable.push_back(string("Veto 2nd muon"));
  CutsSelecTable.push_back(string("Veto electron"));
  char LabelNJets[100];
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
  CutsSelecTable.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTable.push_back(string(LabelNJets));
  CutsSelecTable.push_back(string("Events used"));
  CutsSelecTable.push_back(string("All 4 partons matched"));
 
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
      

  ////////////////////////
  // PileUp Reweighting //
  ////////////////////////

  // OLD 2D REWEIGHING

  LumiReWeighting LumiWeights;
  if (systematic == 10 || systematic == 11 || systematic == 12) { // take 1/fb histo for PU syst for now
    cout << "PU Reweighing, taking 1/fb histogram for PU Syst!!!" << endl;
    LumiWeights = LumiReWeighting("ReweightHistos/pileup/WJets-summer11.root", "ReweightHistos/pileup/data_1fb.root", "pileup2", "pileup"); // summer 11 1.1/fb
  }
  else {
    cout << "PU Reweighing, taking 2/fb histogram!!!" << endl;
    LumiWeights = LumiReWeighting("ReweightHistos/pileup/WJets-summer11.root", "ReweightHistos/pileup/data_2fb.root", "pileup2", "pileup"); // summer 11 2.14/fb
  }

  reweight::PoissonMeanShifter PShiftDown_ = reweight::PoissonMeanShifter(-0.6);
  reweight::PoissonMeanShifter PShiftUp_ = reweight::PoissonMeanShifter(0.6);

  // NEW 3D REWEIGHING
  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root", "PileUpReweighting/pileup_FineBin_2011Data_UpToRun173692.root", "pileup", "pileup"); // 2.1/fb Run2011A
  //Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root", "PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup"); // 4.6/fb Run2011A+B

  //Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root", "PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup"); // FALL114.6/fb Run2011A+B

  if (systematic == 10)
    Lumi3DWeights.weight3D_init(1.08);
  else if (systematic == 11)
    Lumi3DWeights.weight3D_init(0.92);
  else
    Lumi3DWeights.weight3D_init(1.0);

  cout << " Initialized LumiReWeighting stuff" << endl;
  
  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++) {
    //   for (unsigned int d = 2; d < 3; d++) {

    if (runSpecificSample && d != nSpecSample) continue;

    int nSelected=0;
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    //open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;
    
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
    
    /////////////////////////////////////
    /// Initialize JEC factors
    /////////////////////////////////////
   	    
    vector<JetCorrectorParameters> vCorrParam;

    // Create the JetCorrectorParameter objects, the order does not matter.
   // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L3Absolute.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L2Relative.txt");
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L1FastJet.txt");

    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);

    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
      {
	JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/START42_V17_AK5PFchs_L2L3Residual.txt");
	vCorrParam.push_back(*ResJetCorPar);
      }
    
    
    //OLD
    //JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
    //NEW
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
    
    // true means redo also the L1
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
    map<string,TH2F*> histo2D_btag;
    map<string,TGraphErrors*> graphErr_btag;
  

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
      
      if (WmassFit) {
	
	//cout << "Fitted Wmass: " << WmassFit->GetParameter(1) << "+-" << WmassFit->GetParameter(2) << endl;
	Chi2Wmass = WmassFit->GetParameter(1); 
	SigmaChi2Wmass = WmassFit->GetParameter(2); 	

      }

      TF1* TopmassFit = (TF1*)res->Get("hadronicRecoTopMass_Fitted");
      
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
    int itrigger = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
      //  for (unsigned int ievt = 0; ievt < 500; ievt++)
    {
      
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
            
      if(ievt%1000 == 0)
	std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected: " << nSelected << flush<<"\r";

      ////////////////
      // LOAD EVENT //
      ////////////////

      //TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
           
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {
        genjets = treeLoader.LoadGenJet(ievt);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }

      /////////////////////////////////
      // DETERMINE EVENT SCALEFACTOR //
      /////////////////////////////////
      
      // scale factor for the event
      float scaleFactor = 1.;
      float scaleFactorOLDREW = 1.;
      
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
      
      // Apply Jet Corrections on-the-fly
      if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) {
	//jetTools->correctJets(init_jets, vertex);
	jetTools->correctJets(init_jets_corrected,event->kt6PFJetsPF2PAT_rho());
      } else {
	jetTools->correctJets(init_jets_corrected,event->kt6PFJetsPF2PAT_rho());
      }

      // after correction
      //for (unsigned int j=0;j<init_jets_corrected.size();j++)
      //cout << "after new JES jet " << j << " pt " << init_jets_corrected[j]->Pt() << " eta " << init_jets_corrected[j]->Eta() << endl;
      //cout << "****** end event *******"<< endl;
      //cout << endl;
      //exit(1);

      // PU reweighting

      // old method

      double lumiWeightOLD = LumiWeights.ITweight( event->nPu(0) );

      if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
	lumiWeightOLD=1;

      //if (event->nPu(0) >= 25) cout << "#pu " << event->nPu(0) << " shift +0.6 " << PShiftUp_.ShiftWeight( event->nPu(0) ) << endl;

      if (systematic == 10)
	lumiWeightOLD = lumiWeightOLD*PShiftUp_.ShiftWeight( event->nPu(0) );
	
      else if (systematic == 11) 
	lumiWeightOLD = lumiWeightOLD*PShiftDown_.ShiftWeight( event->nPu(0) );

      // 3D method

      //double avPU = ((double)event->nPu(-1) + (double)event->nPu(0) + (double)event->nPu(1))/3. ; // average in 3 BX!!!, as recommended
      //double lumiWeight = LumiWeights.ITweight( (int)avPU );

      double lumiWeight = 1;//LumiWeights.ITweight( event->nPu(0) );

      if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
	lumiWeight=1;
      else if(dataSetName.find("Spring11") != string::npos) {
	lumiWeight=Lumi3DWeights.weight3D(event->nPu(0),event->nPu(0),event->nPu(0));
	//cout << "lol" << endl;
	//exit(0);
      } else
	lumiWeight=Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));

      //if (event->nPu(0) >= 25) cout << "#pu " << event->nPu(0) << " shift +0.6 " << PShiftUp_.ShiftWeight( event->nPu(0) ) << endl;

      //if (systematic == 10)
	//lumiWeight = lumiWeight*PShiftUp_.ShiftWeight( event->nPu(0) );
	//lumiWeight = lumiWeight*PShiftUp_.ShiftWeight( (int)avPU );
	//lumiWeight=Lumi3DWeights_UP.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
	
      //else if (systematic == 11) 
	//lumiWeight = lumiWeight*PShiftDown_.ShiftWeight( event->nPu(0) );
	//lumiWeight = lumiWeight*PShiftDown_.ShiftWeight( (int)avPU );
	//lumiWeight=Lumi3DWeights_DOWN.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));

      //cout << "scaleFactor before = " << scaleFactor << " " << lumiWeight <<endl;
      scaleFactor = scaleFactor*lumiWeight;
      scaleFactorOLDREW = scaleFactorOLDREW*lumiWeightOLD;
      //cout << "scaleFactor after = " << scaleFactor << " " << scaleFactorOLDREW <<  endl;

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
	
	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {

	  // the first 1.1/fb part of RUN2011A (may10->promptv4)
	  if( event->runId() <= 163261 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun, iFile);
	  else if( event->runId() >= 163270 && event->runId() <= 163869 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v6"), currentRun, iFile);
	  else if( event->runId() >= 165088 && event->runId() <= 165633 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v8"), currentRun, iFile);
	  else if( event->runId() == 166346 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v10"), currentRun, iFile);
	  else if( event->runId() >= 165970 && event->runId() <= 167043 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v9"), currentRun, iFile);
	  else if( event->runId() >= 167078 && event->runId() <= 167913 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v11"), currentRun, iFile);
	  
	  // the other part of RUN2011A (another 1/fb) (aug05,promptv6)

	  else if( event->runId() >= 170249 && event->runId() <= 172619 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_v8"), currentRun, iFile);
	  else if( event->runId() >= 172620 && event->runId() <= 173198 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_v8"), currentRun, iFile);
	  else if (event->runId() >= 173236 && event->runId() <= 173692)
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_v9"), currentRun, iFile);

	  // RUN2011B (promptv1)

	  else if( event->runId() ==  176928 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
	  else if( event->runId() == 176982 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);

	  else if( event->runId() >= 175860 && event->runId() <= 176469 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v3"), currentRun, iFile);
	  else if( event->runId() >=  176548 && event->runId() <=  176702 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v3"), currentRun, iFile);
	  else if( event->runId() >=  176797 && event->runId() <=  176889 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v3"), currentRun, iFile);
	  else if( event->runId() >=  176929 && event->runId() <=  176959 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v3"), currentRun, iFile);
	  else if( event->runId() >=  177053 && event->runId() <=  177452 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v3"), currentRun, iFile);

	  else if( event->runId() >=  176545 && event->runId() <=  176547 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
	  else if( event->runId() >=  176765 && event->runId() <=  176796 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
	  
	  else if( event->runId() >=  177718 && event->runId() <=  178380 ) // TopTree ID 804
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
	  else if( event->runId() >=  178420 && event->runId() <=  178479 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v6"), currentRun, iFile);                                
	  else if( event->runId() >=  178703 && event->runId() <=  179889 ) // TopTree ID 816
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v6"), currentRun, iFile);
	  else if( event->runId() >=  179959 && event->runId() <=  180252 )
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v7"), currentRun, iFile);

	  else
	    cout << "Unknown run for HLTpath selection: " << event->runId() << endl;
	  
	  //  HLT_IsoMu17_v5        160431-163261                          32.88(/pb)    32.88(/pb)
	  //  HLT_IsoMu17_v6        163270-163869                         163.81(/pb)   163.81(/pb)
	  //  HLT_IsoMu17_v8        165088-165633                         137.51(/pb)   137.46(/pb)
	  //  HLT_IsoMu17_v9        165970-167043, except 166346          530.06(/pb)   529.50(/pb)
	  //  HLT_IsoMu17_v10       166346                                  4.26(/pb)     4.26(/pb)
	  //  HLT_IsoMu17_v11       167078-167151                          22.07(/pb)    22.00(/pb)
	  //  HLT_IsoMu17_v11       167281-167913                         222.84(/pb)   222.84(/pb)
	} else {
	  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun);
	  if (itrigger == 9999)
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v4"), currentRun); // Spring11: HLT_Mu15_v1
	  if (itrigger == 9999)
	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu9"), currentRun); // Fall10: HLT_Mu9
	  if (itrigger == 9999)
	  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v14"), currentRun); // Fall11

	  //itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_v1"), currentRun);
	  //  Summer11 MC:        HLT_IsoMu17_v5        HLT_Mu15_v2
	  //  Spring11 MC:        HLT_IsoMu17_v4        HLT_Mu15_v1 or HLT_Mu17_v1
	  //  Fall10 MC:          HLT_IsoMu9            HLT_Mu9 or HLT_Mu11
	}
      }


      if (itrigger == 9999) {
	  
	cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT IN RUN " << event->runId() << endl;

	exit(1);
	  
      }

      /////////////////////////////////////////////////////////////////////////////
      // JES SYSTEMATICS && SMEAR JET RESOLUTION TO MIMIC THE RESOLUTION IN DATA //
      /////////////////////////////////////////////////////////////////////////////

      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) ) {
	
	// Match RecoJets with GenJets
	vector< pair<size_t, size_t> > indexVector; //first index = RecoJet, second index = GenJet
	vector<bool> mLock(genjets.size(),false);   // when locked, genJet is already matched to a recoJet
	for(size_t i=0; i<init_jets_corrected.size(); i++) {
	  pair<size_t, size_t> tmpIndex;
	  float minDR = 9999.;
	  for(size_t j=0; j<genjets.size(); j++) {
	    //cout << "Do I smell GenJets?" << endl;
	    if( ! mLock[j] ) {
	      if( init_jets_corrected[i]->DeltaR(*genjets[j]) < 0.4 && init_jets_corrected[i]->DeltaR(*genjets[j]) < minDR ) {
		minDR = init_jets_corrected[i]->DeltaR(*genjets[j]);
		tmpIndex = pair<size_t, size_t>(i,j);
	      }
	    }
	  }
	  if(minDR < 999.){
	    mLock[tmpIndex.second] = true;
	    indexVector.push_back(tmpIndex);
	  }
	}
	  
	// Apply correction for jet energy resolution on-the-fly, only for recoJets matched with a genJet
	for(size_t i=0; i<indexVector.size(); i++) {
	  if( genjets[indexVector[i].second]->Pt() < 15 ) continue;
	  float corrFactor = 0.1;
	  bool JER_minus = false, JER_plus = false; // only one of them true! never both of them!
	  
	  if (systematic == 6)
	    JER_minus=true;
	  else if (systematic == 7)
	    JER_plus=true;
	  
	  float fabsEta = fabs(init_jets_corrected[indexVector[i].first]->Eta());
	  if(JER_minus) {
	    if(fabsEta <= 1.5) corrFactor = 0.0;
	    else if(fabsEta < 2.0 && fabsEta > 1.5) corrFactor = -0.05;
	    else corrFactor = -0.1;
	  }
	  else if(JER_plus) {
	    if(fabsEta <= 1.5) corrFactor = 0.2;
	    else if(fabsEta < 2.0 && fabsEta > 1.5) corrFactor = 0.25;
	    else corrFactor = 0.3;
	  }
	  float deltapt = ( init_jets_corrected[indexVector[i].first]->Pt() - genjets[indexVector[i].second]->Pt() ) * corrFactor;
	  float ptscale = max(0.0, ( init_jets_corrected[indexVector[i].first]->Pt() + deltapt) / init_jets_corrected[indexVector[i].first]->Pt() );
	  //if(ptscale > 0.0)
	  //init_jets_corrected[indexVector[i].first]->SetPxPyPzE(init_jets_corrected[indexVector[i].first]->Px()*ptscale, init_jets_corrected[indexVector[i].first]->Py()*ptscale,init_jets_corrected[indexVector[i].first]->Pz()*ptscale, init_jets_corrected[indexVector[i].first]->E()*ptscale);
	  
	}

	// Correct jets for JES uncertainy systematics
	
	//for (unsigned int j=0;j<init_jets_corrected.size();j++)
	//cout << "jet " << j << " pt " << init_jets_corrected[j]->Pt() << " eta " << init_jets_corrected[j]->Eta() << endl;
	
	if (systematic == 1)
	  jetTools->correctJetJESUnc(init_jets_corrected, "minus",1);
	else if (systematic == 2)
	  jetTools->correctJetJESUnc(init_jets_corrected, "plus",1);
	
	//for (unsigned int j=0;j<init_jets_corrected.size();j++)
	//cout << "after jet " << j << " pt " << init_jets_corrected[j]->Pt() << " eta " << init_jets_corrected[j]->Eta() << endl;
	
	//exit(1);

	
      }

      /////////////////////
      // EVENT SELECTION //
      /////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets);

      //selection.setJetCuts();
      //selection.setLooseMuonCuts();
      //selection.setMuonCuts(); // to make shure our normal mu+jets cuts are not overwritten by dilepton cuts

      // SYNC WITH GERRIT
      selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1);

      selection.setMuonCuts(30.,2.1,0.125,10,0.02,0.3,1,1,1);
      //selection.setMuonCuts(35.,2.1,0.125,10,0.02,0.3,1,1,1);
      //selection.setMuonCuts(20.,2.1,0.125,10,0.02,0.3,1,1,1); // refSelV4 values

      
      /*if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
      {
        // Apply the scraping veto
        bool isBeamBG = true;
        if(event->nTracks() > 10)
        {
          if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
            isBeamBG = false;
        }
        if(isBeamBG) continue;*/
        
        // Apply the JSON
        // bool passedJSON = treeLoader.EventPassedJSON(datasets[d], event->runId(), event->lumiBlockId());
        // if(!passedJSON) continue;
      // }

      bool trigged;

      trigged=treeLoader.EventTrigged (itrigger);

      bool isGoodPV = false;
        
      isGoodPV = selection.isPVSelected(vertex, anaEnv.PVertexNdofCut, anaEnv.PVertexZCut, anaEnv.PVertexRhoCut);

      vector<TRootJet*> selectedJets;
      vector<TRootMuon*> selectedMuons;

      if (init_jets_corrected.size() > 0) {
	//if (init_jets_corrected[0]->jetType() == 1 || doPF2PAT) { // calojets
	  //cout << "Selecting for caloJets" << endl;
	selectedJets = selection.GetSelectedJets(true);
	selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	  //	} else {
	  //cout << "Selecting for PF/JPT jets" << endl;
	  // vector<TRootMuon*> overlapMuons = selection.GetSelectedMuons(vertex[0]);
	  //selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1); // refSelV4 values
	  //selectedJets = selection.GetSelectedJets(overlapMuons,true);
	  //selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);	
	  //	}
      }
      
      vector<TRootMuon*> selectedLooseMuons = selection.GetSelectedLooseMuons();
      vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(false);

      //vector<TRootJet*> selectedLooseJets = selection.GetSelectedJets(20., anaEnv.JetsEtaCutSR, true);
      

      vector<TRootMCParticle*> mcParticles;
      
      if(dataSetName.find("TTbarJets_SemiMu") == 0)
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }

      bool eventselected = false;
      if (trigged && isGoodPV) {
	if (selectedMuons.size() == 1) {
	  //cout << "selectedLooseMuons: " << selectedLooseMuons.size() << endl;
	  if (selectedLooseMuons.size() == 1) {
	    if (vetoElectrons.size() == 0) {
	      if (selectedJets.size() >= 4) {
	      //if (selectedJets.size() == 4) {
		eventselected = true;
	      }
	    }
	  }
	}
      }

      /*if (eventselected) {
	if (selectedJets[0]->Pt() < 70 || selectedJets[1]->Pt() < 50) 
	  eventselected=false;
	  }*/

      if (!eventselected) continue;

      //for (unsigned int j=0;j<selectedJets.size();j++)
	//cout << "after event selection jet " << j << " pt " << selectedJets[j]->Pt() << " eta " << selectedJets[j]->Eta() << endl;
      //cout << "****** end event *******"<< endl;
      //cout << endl;
      //exit(1);

      if (systematic == 3)
	jetTools->correctJetJESUnc(selectedJets, "minus");
      else if (systematic == 4)
	jetTools->correctJetJESUnc(selectedJets, "plus");

      if(dataSetName.find("WJets_TuneZ2_scaledown") == 0) {
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
      }

      nSelected++;

      ////////////////////////////////////
      // JET PARTON MATCHING FOR MASSES //
      ////////////////////////////////////

      // this is only necessary for the chi2 jetcomb input masses. if usemassesandresolutions == true this will not be run

      bool all4PartonsMatched = false; // True if the 4 ttbar semi-mu partons are matched to 4 jets (not necessarily the 4 highest pt jets)
      pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parto
      int MCPermutation[4]; 

      if (!useMassesAndResolutions && dataSetName.find("TTbarJets_SemiMu") == 0) {

	sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)     

	for (unsigned int i=0;i<4;i++) MCPermutation[i] = -1;

	vector<TRootMCParticle*> mcParticlesMatching_; 
	vector<TRootJet*> selectedJets_; 

	leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
	mcParticlesMatching_.clear();
	selectedJets_.clear();
		  
	vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
        
	bool muPlusFromTop = false, muMinusFromTop = false;
	int nTTbarQuarks = 0;
	for(unsigned int i=0; i<mcParticles.size(); i++) {
	  if( mcParticles[i]->status() != 3) continue;
	  
	  if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 ) {
	    if(muMinusFromTop) cerr<<"muMinusFromTop was already true"<<endl;
	    muMinusFromTop = true;
	  }
	  if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 ) {
	    if(muPlusFromTop) cerr<<"muPlusFromTop was already true"<<endl;
	    muPlusFromTop = true;
	  }
	  
	  if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ) {
	    mcParticlesTLV.push_back(*mcParticles[i]);
	    mcParticlesMatching_.push_back(mcParticles[i]);
	  }
	}
	
	if(muPlusFromTop && muMinusFromTop)
	  cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
	
	for(unsigned int i=0; i<4; i++)
	  selectedJetsTLV.push_back(*selectedJets[i]);
	
	JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
	
	if(matching.getNumberOfAvailableCombinations() != 1)
	  cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
	
	vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
	
	for(unsigned int i=0; i<mcParticlesTLV.size(); i++) {
	  int matchedJetNumber = matching.getMatchForParton(i, 0);
	  if(matchedJetNumber != -1)
	    JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
	}
	
	for(unsigned int i=0; i<JetPartonPair.size(); i++) {
	  unsigned int j = JetPartonPair[i].second;
	  
	  if( fabs(mcParticlesMatching_[j]->type()) < 6 )
	    {
	      if( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -6 )
		  || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == 6 ) )
		{
		  if(hadronicWJet1_.first == 9999) {
		    hadronicWJet1_ = JetPartonPair[i];
		    MCPermutation[0] = JetPartonPair[i].first;
		  } else if(hadronicWJet2_.first == 9999) {
		    hadronicWJet2_ = JetPartonPair[i];
		    MCPermutation[1] = JetPartonPair[i].first;
		  } else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
		}
	    }
	  if( fabs(mcParticlesMatching_[j]->type()) == 5 )
	    {
	      if( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -6 )
		  || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 6 ) ) {
		hadronicBJet_ = JetPartonPair[i];
		MCPermutation[2] = JetPartonPair[i].first;
	      } else if( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == 6 )
			 || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == -6 ) ) {
		leptonicBJet_ = JetPartonPair[i];
		MCPermutation[3] = JetPartonPair[i].first;
	      }
	    }
	}
	
	if(hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999 && hadronicBJet_.first != 9999 && leptonicBJet_.first != 9999) {
	  histo1D["hadronicPartonTopMass"]->Fill((*mcParticlesMatching_[hadronicWJet1_.second]+*mcParticlesMatching_[hadronicWJet2_.second]+*mcParticlesMatching_[hadronicBJet_.second]).M());
	  histo1D["hadronicPartonWMass"]->Fill((*mcParticlesMatching_[hadronicWJet1_.second]+*mcParticlesMatching_[hadronicWJet2_.second]).M());
	  
	  all4PartonsMatched = true;
	  
	  //cout << "pom" << endl;
	  
	}
	
	selecTable.Fill(d,11,scaleFactor);
	
      }

      if (all4PartonsMatched && !useMassesAndResolutions && dataSetName.find("TTbarJets_SemiMu") == 0) {

	selecTable.Fill(d,12,scaleFactor);

	float WMass = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]).M();
	float TopMass = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]+*selectedJets[hadronicBJet_.first]).M();
      
	//histo1D["hadronicRecoWMass"]->Fill(WMass,lumiWeight);
	//histo1D["hadronicRecoTopMass"]->Fill(TopMass,lumiWeight);

	histo1D["hadronicRecoWMass"]->Fill(WMass);
	histo1D["hadronicRecoTopMass"]->Fill(TopMass);

	// to be fixed

	//cout << "WMass " << WMass << " TopMass " << TopMass << endl; exit(0);
	
      }

      ///////////////////////////////////
      // !!! FILLING THE BTAGTREES !!! //
      ///////////////////////////////////

      // now we use the loaded values and do the chi2 jetcombination + fill the trees
      if (useMassesAndResolutions) {

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
	if (all4PartonsMatched) {
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

	//---------------//
	// CONTROL PLOTS //
	//---------------//

	if (histo1D_btag.find("JetCombMinChi2") == histo1D_btag.end())
	  histo1D_btag["JetCombMinChi2"] = new TH1F("JetCombMinChi2","Minimal #Chi^{2} value of the 12 combinations;#Chi^{2};#events",200,0,100);

	histo1D_btag["JetCombMinChi2"]->Fill(smallestChi2);

	if (histo1D_btag.find("BestJetCombMLB") == histo1D_btag.end())
	  histo1D_btag["BestJetCombMLB"] = new TH1F("BestJetCombMLB","M_{lb} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{lb} (GeV/c^{2});#events",600,0,1200);

	histo1D_btag["BestJetCombMLB"]->Fill((*selectedMuons[0]+*selectedJets[Permutation[3]]).M());

	if (histo2D_btag.find("Top_VS_W_Mass") == histo2D_btag.end())
	  histo2D_btag["Top_VS_W_Mass"] = new TH2F("Top_VS_W_Mass","M_{top} VS m_{W} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{top} (GeV/c^{2});m_{W} (GeV/c^{2})",200,100,500,200,0,400);
	
	float mtop = (*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]+*selectedJets[Permutation[2]]).M();
	float mW = (*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]).M();
	
	histo2D_btag["Top_VS_W_Mass"]->Fill(mtop,mW);

	if (histo1D_btag.find("BestJetCombMTop") == histo1D_btag.end())
	  histo1D_btag["BestJetCombMTop"] = new TH1F("BestJetCombMTop","M_{top} for jetcomb with smallest #Chi^{2} value of the 12 combinations;m_{top} (GeV/c^{2});#events",80,100,500);
	
	histo1D_btag["BestJetCombMTop"]->Fill(mtop);

	//-----------------//
	// do some data-mc //
	//-----------------//

	// jet plots for the BTV guys
	
	if (MSPlot.find("Selected_Events_pT_jet1") == MSPlot.end()){
	  MSPlot["Selected_Events_pT_jet1"] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet1", 30, 0, 600, "p_{T} (GeV)");
	  MSPlot["Selected_Events_pT_jet2"] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet2", 30, 0, 600, "p_{T} (GeV)");
	  MSPlot["Selected_Events_pT_jet3"] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet3", 30, 0, 600, "p_{T} (GeV)");
	  MSPlot["Selected_Events_pT_jet4"] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet4", 30, 0, 600, "p_{T} (GeV)");
	  //MSPlot["Selected_Events_pT_jet4"] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet4", 30, 0, 600, "p_{T} (GeV)");
	  MSPlot["Selected_Events_pT_4leadingjets"] = new MultiSamplePlot(datasets, "Selected_Events_pT_4leadingjets", 30, 0, 600, "p_{T} (GeV)");
	  MSPlot["Selected_Events_pT_alljets"] = new MultiSamplePlot(datasets, "Selected_Events_pT_alljets", 30, 0, 600, "p_{T} (GeV)");
	}

	MSPlot["Selected_Events_pT_jet1"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["Selected_Events_pT_jet2"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["Selected_Events_pT_jet3"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["Selected_Events_pT_jet4"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);

	for (unsigned int q=0; q<selectedJets.size(); q++) {

	  MSPlot["Selected_Events_pT_alljets"]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);

	  if (q<4)
	    MSPlot["Selected_Events_pT_4leadingjets"]->Fill(selectedJets[q]->Pt(), datasets[d], true, Luminosity*scaleFactor);

	}

	// other data-mc

	if (MSPlot.find("Selected_Events_nPV3D") == MSPlot.end())
	  MSPlot["Selected_Events_nPV3D"] = new MultiSamplePlot(datasets, "Selected_Events_nPV3D", 36, -0.5, 35.5, "nPV");

	MSPlot["Selected_Events_nPV3D"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("Selected_Events_nPV") == MSPlot.end())
	  MSPlot["Selected_Events_nPV"] = new MultiSamplePlot(datasets, "Selected_Events_nPV", 36, -0.5, 35.5, "nPV");

	MSPlot["Selected_Events_nPV"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactorOLDREW);

	if (MSPlot.find("Selected_Events_nPV_NONREW") == MSPlot.end())
	  MSPlot["Selected_Events_nPV_NONREW"] = new MultiSamplePlot(datasets, "Selected_Events_nPV_NONREW", 36, -0.5, 35.5, "nPV");

	MSPlot["Selected_Events_nPV_NONREW"]->Fill(vertex.size(), datasets[d], true, Luminosity);
	
	if (MSPlot.find("Selected_Events_JetCombMinChi2") == MSPlot.end())
	  MSPlot["Selected_Events_JetCombMinChi2"] = new MultiSamplePlot(datasets, "Selected_Events_JetCombMinChi2", 100, 0, 100, "JetCombMinChi2");
	
	MSPlot["Selected_Events_JetCombMinChi2"] ->Fill(smallestChi2, datasets[d], true, Luminosity*scaleFactor);
	
	if (MSPlot.find("Selected_Events_JetCombMinChi2Exp") == MSPlot.end())
	  MSPlot["Selected_Events_JetCombMinChi2Exp"] = new MultiSamplePlot(datasets, "Selected_Events_JetCombMinChi2Exp", 100, 0, 100, "JetCombMinChi2Exp");

	MSPlot["Selected_Events_JetCombMinChi2Exp"] ->Fill(exp(smallestChi2), datasets[d], true, Luminosity*scaleFactor);

	// MLB
	if (MSPlot.find("BestJetCombMLB") == MSPlot.end())
	  MSPlot["BestJetCombMLB"] = new MultiSamplePlot(datasets, "BestJetCombMLB", 75, 0, 1200, "BestJetCombMLB");

	MSPlot["BestJetCombMLB"] ->Fill((*selectedMuons[0]+*selectedJets[Permutation[3]]).M(), datasets[d], true, Luminosity*scaleFactor);

	// MW
	if (MSPlot.find("BestJetCombMW") == MSPlot.end())
	  MSPlot["BestJetCombMW"] = new MultiSamplePlot(datasets, "BestJetCombMW", 50, 0, 200, "BestJetCombMW");

	MSPlot["BestJetCombMW"] ->Fill((*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]).M(), datasets[d], true, Luminosity*scaleFactor);

	// MTOP
	if (MSPlot.find("BestJetCombMTop") == MSPlot.end())
	  MSPlot["BestJetCombMTop"] = new MultiSamplePlot(datasets, "BestJetCombMTop", 400, 0, 600, "m_{top} (GeV)","Events / 20 GeV");

	MSPlot["BestJetCombMTop"] ->Fill((*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]+*selectedJets[Permutation[2]]).M(), datasets[d], true, Luminosity*scaleFactor);
	
	if (MSPlot.find("BestJetCombMTop_CHI2") == MSPlot.end())
	  MSPlot["BestJetCombMTop_CHI2"] = new MultiSamplePlot(datasets, "BestJetCombMTop_CHI2", 40, 0, 600 ,"m_{top} (GeV)","Events / 20 GeV");
	
	if (smallestChi2 < 90)
	  MSPlot["BestJetCombMTop_CHI2"] ->Fill((*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]+*selectedJets[Permutation[2]]).M(), datasets[d], true, Luminosity*scaleFactor);
	
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

	if (dataSetName.find("TTbarJets_SemiMu") == 0) {
	  
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
	  for(int iMuO=0; iMuO<selectedLooseMuons.size(); iMuO++){
	    double deltaRjetmuonstmp = selectedJets[Permutation[3]]->DeltaR(*selectedMuons[iMuO]);
	    
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
	
	NTuple->setMlj((*(selectedJets[Permutation[3]])+*selectedMuons[0]).M());

	NTuple->setBtag(0,(*(selectedJets[Permutation[3]])).btag_trackCountingHighEffBJetTags());
	NTuple->setBtag(1,(*(selectedJets[Permutation[3]])).btag_trackCountingHighPurBJetTags());

	if (MSPlot.find("bcand_trackCountingHighEffBJetTags") == MSPlot.end())
	  MSPlot["bcand_trackCountingHighEffBJetTags"] = new MultiSamplePlot(datasets, "bcand_trackCountingHighEffBJetTags", 50, -10, 30, "Btag Discr");
	MSPlot["bcand_trackCountingHighEffBJetTags"] ->Fill((*(selectedJets[Permutation[3]])).btag_trackCountingHighEffBJetTags(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("bcand_trackCountingHighPurBJetTags") == MSPlot.end())
	  MSPlot["bcand_trackCountingHighPurBJetTags"] = new MultiSamplePlot(datasets, "bcand_trackCountingHighPurBJetTags", 50, -10, 30, "Btag Discr");
	MSPlot["bcand_trackCountingHighPurBJetTags"] ->Fill((*(selectedJets[Permutation[3]])).btag_trackCountingHighPurBJetTags(), datasets[d], true, Luminosity*scaleFactor);


	NTuple->setBtag(2,(*(selectedJets[Permutation[3]])).btag_simpleSecondaryVertexHighEffBJetTags());
	NTuple->setBtag(3,(*(selectedJets[Permutation[3]])).btag_simpleSecondaryVertexHighPurBJetTags());

	if (MSPlot.find("bcand_simpleSecondaryVertexHighEffBJetTags") == MSPlot.end())
	  MSPlot["bcand_simpleSecondaryVertexHighEffBJetTags"] = new MultiSamplePlot(datasets, "bcand_simpleSecondaryVertexHighEffBJetTags", 50, 0, 8, "Btag Discr");
	MSPlot["bcand_simpleSecondaryVertexHighEffBJetTags"] ->Fill((*(selectedJets[Permutation[3]])).btag_simpleSecondaryVertexHighEffBJetTags(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("bcand_simpleSecondaryVertexHighPurBJetTags") == MSPlot.end())
	  MSPlot["bcand_simpleSecondaryVertexHighPurBJetTags"] = new MultiSamplePlot(datasets, "bcand_simpleSecondaryVertexHighPurBJetTags", 50, 0, 8, "Btag Discr");
	MSPlot["bcand_simpleSecondaryVertexHighPurBJetTags"] ->Fill((*(selectedJets[Permutation[3]])).btag_simpleSecondaryVertexHighPurBJetTags(), datasets[d], true, Luminosity*scaleFactor);
	
	NTuple->setBtag(4,(*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexBJetTags());
	NTuple->setBtag(5,(*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexMVABJetTags());

	if (MSPlot.find("bcand_combinedSecondaryVertexBJetTags") == MSPlot.end())
	  MSPlot["bcand_combinedSecondaryVertexBJetTags"] = new MultiSamplePlot(datasets, "bcand_combinedSecondaryVertexBJetTags", 50, -1, 2, "Btag Discr");
	MSPlot["bcand_combinedSecondaryVertexBJetTags"] ->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexBJetTags(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("bcand_combinedSecondaryVertexMVABJetTags") == MSPlot.end())
	  MSPlot["bcand_combinedSecondaryVertexMVABJetTags"] = new MultiSamplePlot(datasets, "bcand_combinedSecondaryVertexMVABJetTags", 50, -1, 2, "Btag Discr");
	MSPlot["bcand_combinedSecondaryVertexMVABJetTags"] ->Fill((*(selectedJets[Permutation[3]])).btag_combinedSecondaryVertexMVABJetTags(), datasets[d], true, Luminosity*scaleFactor);

	NTuple->setBtag(6,(*(selectedJets[Permutation[3]])).btag_jetBProbabilityBJetTags());
	NTuple->setBtag(7,(*(selectedJets[Permutation[3]])).btag_jetProbabilityBJetTags());

	if (MSPlot.find("bcand_jetBProbabilityBJetTags") == MSPlot.end())
	  MSPlot["bcand_jetBProbabilityBJetTags"] = new MultiSamplePlot(datasets, "bcand_jetBProbabilityBJetTags", 50, 0, 8, "Btag Discr");
	MSPlot["bcand_jetBProbabilityBJetTags"] ->Fill((*(selectedJets[Permutation[3]])).btag_jetBProbabilityBJetTags(), datasets[d], true, Luminosity*scaleFactor);

	if (MSPlot.find("bcand_jetProbabilityBJetTags") == MSPlot.end())
	  MSPlot["bcand_jetProbabilityBJetTags"] = new MultiSamplePlot(datasets, "bcand_jetProbabilityBJetTags", 50, 0, 3, "Btag Discr");
	MSPlot["bcand_jetProbabilityBJetTags"] ->Fill((*(selectedJets[Permutation[3]])).btag_jetProbabilityBJetTags(), datasets[d], true, Luminosity*scaleFactor);
	
	//vector <TRootCaloJet*> selectedCaloJets = convertToCaloJets(selectedJets);
	
	NTuple->setPt((*(selectedJets[Permutation[3]])).Pt());
	NTuple->setPthadb((*(selectedJets[Permutation[3]])).Pt());
	NTuple->setE((*(selectedJets[Permutation[3]])).E());
	NTuple->setEta((*(selectedJets[Permutation[3]])).Eta());//commented for deltaR
	NTuple->setWeight(1/(datasets[d]->EquivalentLumi()));
	NTuple->setScaleFactor3D(scaleFactor);
	NTuple->setScaleFactor(scaleFactorOLDREW);

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
	NTuple->setPtMuon((*selectedMuons[0]).Pt());   
	NTuple->setEMuon((*selectedMuons[0]).E());   
	NTuple->setEtaMuon((*selectedMuons[0]).Eta());  
	
	NTuple->setLepBlCandMass( (*selectedJets[Permutation[3]]+*selectedMuons[0]).M() );
	NTuple->setHadTopCandMass( (*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]+*selectedJets[Permutation[2]]).M() );
	NTuple->setHadWCandMass( (*selectedJets[Permutation[0]]+*selectedJets[Permutation[1]]).M() );
	NTuple->setR_inAll(R_inAll);
	NTuple->setR_inHad(R_inHad);
	NTuple->setR_or_lepb_inHad(R_or_lepb_inHad);
	NTuple->setAllMatched(allMatched);


	//Setting some extra variables which could help reducing the amount of background:
	
	NTuple->setnJets(init_jets_corrected.size());
	NTuple->setBtag_trackCountingHighEffBJetTags_hadb((*(selectedJets[Permutation[2]])).btag_trackCountingHighEffBJetTags());
	
	
	TLorentzVector muon0, jet0, jet1, jet2, jet3, met0;
	TRootMET* met = (TRootMET*) mets[0]; 
	met0.SetPxPyPzE(met->Px(),met->Py(),met->Pz(),met->Energy());
	muon0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
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
	    NTuple->setMljControl3((*selectedJets[0]+*selectedMuons[0]).M());
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
	    NTuple->setMljControl((*(selectedJets[Permutation[0]])+*selectedMuons[0]).M());
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
	    NTuple->setMljControl2((*(selectedJets[Permutation[1]])+*selectedMuons[0]).M());
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
	    NTuple->setMljControl((*(selectedJets[Permutation[1]])+*selectedMuons[0]).M());
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
	    NTuple->setMljControl2((*(selectedJets[Permutation[0]])+*selectedMuons[0]).M());
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

    cout << "+> " << nSelected << " events where selected"<< endl;

    ////////////////////////////////////////////////
    // STORE INPUT MASSES FOR CHI2 JETCOMBINATION //
    ////////////////////////////////////////////////

    if (!useMassesAndResolutions && dataSetName.find("TTbarJets_SemiMu") == 0) {
	
      string filename = "BtagMassPlots";

      //if (argc > 2)
      //filename += "_"+string(argv[2]);

      filename += ".root";
      
      cout << "INFO: Creating output file to store the masses for the Chi2 jetcombiner: " << filename << endl;    

      TFile* outfile = new TFile(filename.c_str(),"RECREATE");
      
      outfile->cd();
      
      histo1D["hadronicRecoWMass"]->Write();
      histo1D["hadronicRecoTopMass"]->Write();

      // fit the distributions
      
      for (unsigned int f=0; f<2;f++) {

	TH1F* histo;

	if (f==0) histo=histo1D["hadronicRecoWMass"]; if (f==1) histo=histo1D["hadronicRecoTopMass"];

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
   
  
    //Selection tables
    selecTable.TableCalculator(false, true, true, true);
    string selectiontable = "SelectionTable"+postfix+"_BTAG";
    /*if (argc >= 3){
      string sample=string(argv[2]);
      selectiontable = selectiontable +"_"+sample;
      }*/
    selectiontable = selectiontable +".tex"; 	
    selecTable.Write(selectiontable.c_str());
    
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
	temp->Draw(false, name, false, true, true, true, true);
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
