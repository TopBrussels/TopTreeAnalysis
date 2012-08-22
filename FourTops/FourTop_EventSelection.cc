/////////////
////////////
// TODO + COMMENTS
//1.  Validate  all those electron ID variables...
//2a. Add JER systematic calculation, using GEN Jets.
//2b. Add btag scale factor systematic calculation, need to update number for appropriate taggers.
//3. Change to 52X JES correction
//4. Need  to do jet-parton matching to plot hadronic W mass and  estimate sigma(M_{w}), example in Petra's code.]
//5. Add e + jets final state, prob leave this till last.


#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../Tools/interface/JetTools.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../Reconstruction/interface/MakeBinning.h"
#include "../MCInformation/interface/LumiReWeighting.h"
#include "../Reconstruction/interface/MEzCalculator.h"

//#include "FourTopTree/interface/InclFourthGenTree.h"


#include "../macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
map<string,MultiSamplePlot*> MultiPadPlot;

struct HighestTCHEBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    	return j1->btag_trackCountingHighEffBJetTags() > j2->btag_trackCountingHighEffBJetTags();
    }
};
struct HighestCVSBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    	return j1->btag_combinedSecondaryVertexBJetTags() > j2->btag_combinedSecondaryVertexBJetTags();
    }
};

bool match;

int main (int argc, char *argv[])
{

  int doJESShift = 0; // 0: off 1: minus 2: plus
  cout << "doJESShift: " << doJESShift << endl;

  int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
  cout << "doJERShift: " << doJERShift << endl;

  int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2: plus
  cout << "dobTagEffShift: " << dobTagEffShift << endl;

  int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
  cout << "domisTagEffShift: " << domisTagEffShift << endl;

  int doPUShift = 0; //0: off (except nominal PU reweighting) 1: minus 2: plus
  cout << "doPUShift: " << doPUShift << endl;

  string btagger = "CSVM";
// b-tag scalefactor => TCHEL: data/MC scalefactor = 0.95 +- 0.10,    TCHEM: data/MC scalefactor = 0.94 +- 0.09
// mistag scalefactor => TCHEL: data/MC scalefactor = 1.11 +- 0.12,    TCHEM: data/MC scalefactor = 1.21 +- 0.17
  float scalefactorbtageff, mistagfactor;
  if(btagger == "TCHPM"  || btagger == "TCHET"  ||  btagger == "SSV" ){
    cout<<"This tagger ("<< btagger <<")is not commisioned in 2012, please use CSV, TCHP or JetProb"<<endl;
    exit(1);
}
  else if(btagger == "TCHEM") //redundant for now, but will use as skeleton for CSVM
  {
  	  if(dobTagEffShift == 0)
		scalefactorbtageff = 0.94;
	  if(dobTagEffShift == 1)
		scalefactorbtageff = 0.85;
	  if(dobTagEffShift == 2)
		scalefactorbtageff = 1.03;
		
	  if(domisTagEffShift == 0)
		mistagfactor = 1.21;
	  if(domisTagEffShift == 1)
		mistagfactor = 1.04;
	  if(domisTagEffShift == 2)
		mistagfactor = 1.38;
  }
  float workingpointvalue = 9999; //working points updated to 2012 BTV-POG recommendations.
 
  if(btagger == "TCHPM"  || btagger == "TCHET"  ||  btagger == "SSV" ){
    cout<<"This tagger ("<< btagger <<")is not commisioned in 2012, please use CSV, TCHP or JetProb"<<endl;
    exit(1); 
  }
  else if(btagger == "TCHPL")
     workingpointvalue = 1.470;
  else if(btagger == "TCHPT")
    workingpointvalue = 3.42;
  else if(btagger == "CSVL")
     workingpointvalue = .244;	
  else if(btagger == "CSVM")
    workingpointvalue = .679;
  else if(btagger == "CSVT")
    workingpointvalue = .898;

  clock_t start = clock();

 cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FourTop search ! "           << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  setTDRStyle();
  //setGregStyle();
  //setMyStyle();

  string postfix = "_EventSelection"; // to relabel the names of the output file

  if (doJESShift == 1)
    postfix= postfix+"_JESMinus";
  if (doJESShift == 2)
    postfix= postfix+"_JESPlus";
  if (doJERShift == 1)
    postfix= postfix+"_JERMinus";
  if (doJERShift == 2)
    postfix= postfix+"_JERPlus";
  if (dobTagEffShift == 1)
    postfix= postfix+"_bTagMinus";
  if (dobTagEffShift == 2)
    postfix= postfix+"_bTagPlus";
  if(domisTagEffShift == 1)
    postfix= postfix+"_misTagMinus";
  if(domisTagEffShift == 2)
    postfix= postfix+"_misTagPlus";

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////// Configuration ///////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 

  string channelpostfix = "";
  string xmlFileName = "";

  bool Electron = false; // use Electron channel?
  bool Muon = true; // use Muon channel?
  if(Electron && Muon){
	cout << "  --> Using both Muon and Electron channel? Choose only one ( since different samples/skims are required)!" << endl;
	exit(1);
  }

  if(Muon){
	cout << " --> Using the Muon channel..." << endl;
	channelpostfix = "_Mu";
	xmlFileName = "../config/myTopFCNCconfig_Muon.xml";
  }
  else if(Electron){
	cout << " --> Using the Electron channel..." << endl;
	channelpostfix = "_El";
	xmlFileName = "../config/myTopFCNCconfig_Electron.xml";
  }

      xmlFileName = "config/test_fullsamples.xml";
 //     xmlFileName = "config/test_2.xml";
      // xmlFileName = "config/refsel.xml";

  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////// AnalysisEnvironment /////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = 2;//anaEnv.Verbose;


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  float Luminosity = 5343.64; //pb^-1??
  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
     cout <<"found sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
     string dataSetName = datasets[d]->Name();
   if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
	{
		  Luminosity = datasets[d]->EquivalentLumi();
		  cout <<"found DATA sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
		  break;
	 }
  }
  cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;

  //Output ROOT file
  string rootFileName ("FourTop"+postfix+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* >   vertex;
  vector < TRootMuon* >     init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* >      init_jets;
  vector < TRootMET* >      mets;

  //Global variable
  TRootEvent* event = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ////////////////// MultiSample plots  //////////////////////////////
  ////////////////////////////////////////////////////////////////////

  MSPlot["RhoCorrection"]              = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
  MSPlot["NbOfVertices"]               = new MultiSamplePlot(datasets, "NbOfVertices", 40, 0, 40, "Nb. of vertices");

  //Muons
  MSPlot["NbOfIsolatedMuons"]          = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfIsolatedElectrons"]      = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
  MSPlot["NbOfExtraIsolatedMuons"]     = new MultiSamplePlot(datasets, "NbOfExtraIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfExtraIsolatedElectrons"] = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
  MSPlot["MuonRelIsolation"] = new MultiSamplePlot(datasets, "MuonRelIsolation", 50, 0, .15, "RelIso");
  MSPlot["MuonPt"]              = new MultiSamplePlot(datasets, "MuonPt", 75, 0, 300, "PT_{#mu}");
  MSPlot["MuonEta"]              = new MultiSamplePlot(datasets, "MuonEta", 25, -2.4, 2.4, "#eta_{#mu}");
  MSPlot["MuonPhi"]              = new MultiSamplePlot(datasets, "MuonPhi", 50, -4, 4, "#phi_{#mu}");
  MSPlot["MuonNValidHits"]              = new MultiSamplePlot(datasets, "MuonNValidHits", 30, 0, 30, "NValidHits_{#mu}");
  MSPlot["Muond0"]              = new MultiSamplePlot(datasets, "Muond0", 50, 0, .05, "d0_{#mu}");
  MSPlot["MuondRJets"]              = new MultiSamplePlot(datasets, "MuondRJets", 50, 0, 10, "dRJets_{#mu}");
  MSPlot["MuonNMatchedStations"]              = new MultiSamplePlot(datasets, "MuonNMatchedStations", 10, 0, 10, "NMatchedStations_{#mu}");
  MSPlot["MuonDistVzPVz"]              = new MultiSamplePlot(datasets, "MuonDistVzPVz", 75, -.3,.3, "DistVzPVz_{#mu}");
  MSPlot["MuonTrackerLayersWithMeasurement"]    = new MultiSamplePlot(datasets, "MuonTrackerLayersWithMeasurement", 25, 0, 25, "nLayers");
  MSPlot["DiMuon_InvMass"]     = new MultiSamplePlot(datasets, "DiMuon_InvMass", 60, 0, 120, "DiMuon_InvMass");
  MSPlot["NbOfLooseMuon"]     = new MultiSamplePlot(datasets, "NbOfLooseMuon", 10, 0, 10, "Nb. of loose muons");
    
  //B-tagging discriminators
  MSPlot["BdiscBJetCand_CSV"]  = new MultiSamplePlot(datasets, "HighestBdisc_CSV", 75, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_m_ch_CSV"]            = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_e_ch_CSV"]            = new MultiSamplePlot(datasets, "HighestBdisc_ee_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_m_ch_TCHE"]           = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_TCHE",100, 0, 50, "TCHE b-disc.");
  MSPlot["HighestBdisc_e_ch_TCHE"]           = new MultiSamplePlot(datasets, "HighestBdisc_ee_ch_TCHE",100, 0, 50, "TCHE b-disc.");
  
  //Jets
  MSPlot["NbOfSelectedJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedLightJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
  MSPlot["NbOfSelectedBJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0, 8, "Nb. of jets");
  MSPlot["JetEta"]                  = new MultiSamplePlot(datasets, "JetEta", 30,-3, 3, "Jet #eta");
  MSPlot["JetPhi"]                  = new MultiSamplePlot(datasets, "JetPhi", 50, -4,4 , "Jet #phi");
  MSPlot["FirstTriJetMass_1BJet"] = new MultiSamplePlot(datasets, "FirstTriJetMass_1BJet", 50, 0, 1000, "m_{bjj}");
  MSPlot["FirstTriJetPt_1BJet"] = new MultiSamplePlot(datasets, "FirstTriJetPt_1BJet", 50, 0, 1000, "Pt_{bjj}");
  MSPlot["FirstDiJetMass"] = new MultiSamplePlot(datasets, "FirstDiJetMass", 50, 0, 200, "m_{jj}");
  MSPlot["SecondTriJetMass_1BJet"] = new MultiSamplePlot(datasets, "SecondTriJetMass_1BJet", 50, 0, 1000, "m_{bjj}");
  MSPlot["SecondTriJetPt_1BJet"] = new MultiSamplePlot(datasets, "SecondTriJetPt_1BJet", 50, 0, 1000, "Pt_{bjj}");
  MSPlot["SecondDiJetMass"] = new MultiSamplePlot(datasets, "SecondTriJetMass", 50, 0, 200, "m_{jj}");
  MSPlot["WMt"] = new MultiSamplePlot(datasets, "WMt", 50, 0, 250, "W Transverse Mass");
  MSPlot["LepWMass"] = new MultiSamplePlot(datasets, "LepWMass", 50, 0, 200, "MuMET");
  MSPlot["MuMetBMasses"] = new MultiSamplePlot(datasets, "MuMetBMasses", 50, 0, 1000, "m_{muMETb}");
  MSPlot["MuMetBPt"] = new MultiSamplePlot(datasets, "MuMetBPt", 50, 0, 1000, "Pt_{muMETb}");
  MSPlot["MuMetBMasses_chi2"] = new MultiSamplePlot(datasets, "MuMetBMasses_chi2", 50, 0, 1000, "\\chi^{2}");
  MSPlot["SelectedJetPt"] = new MultiSamplePlot(datasets, "JetPt", 50, 0, 1000, "PT_{jet}");
  MSPlot["4thJetPt"] = new MultiSamplePlot(datasets, "4thJetPt", 60, 0, 400, "PT_{jet}");
  MSPlot["5thJetPt"] = new MultiSamplePlot(datasets, "5thJetPt", 60, 0, 400, "PT_{jet}");
  MSPlot["6thJetPt"] = new MultiSamplePlot(datasets, "6thJetPt", 60, 0, 400, "PT_{jet}");
  MSPlot["7thJetPt"] = new MultiSamplePlot(datasets, "7thJetPt", 60, 0, 400, "PT_{jet}");
  MSPlot["8thJetPt"] = new MultiSamplePlot(datasets, "8thJetPt", 60, 0, 400, "PT_{jet}");
  MSPlot["SelectedJetPt_light"] = new MultiSamplePlot(datasets, "JetPt_light", 50, 0, 1000, "PT_{lightjet}");
  MSPlot["SelectedJetPt_b"] = new MultiSamplePlot(datasets, "JetPt_b", 50, 0, 1000, "PT_{bjet}");
  MSPlot["HT_SelectedJets"] = new MultiSamplePlot(datasets, "HT_SelectedJets", 50, 0, 1500, "HT");
  MSPlot["MHT_SelectedJets"] = new MultiSamplePlot(datasets, "MHT_SelectedJets", 75, 0, 1000, "MHT");
  MSPlot["MHTSig_SelectedJets"] = new MultiSamplePlot(datasets, "MHTSig_SelectedJets", 75, 0, 30, "MHTSig");
  MSPlot["MET"] = new MultiSamplePlot(datasets, "MET", 75, 0, 700, "MET");
  MSPlot["MET_MHT"]= new MultiSamplePlot(datasets, "MET_MHT", 75, 0, 200, "MET_MHT");
  MSPlot["STLep"] = new MultiSamplePlot(datasets, "STLep", 75, 0, 1000, "STLep");
  MSPlot["STJet"] = new MultiSamplePlot(datasets, "STJet", 75, 0, 1000, "STJet");

  //Electrons
  MSPlot["ElectronPt"]              = new MultiSamplePlot(datasets, "ElectronPt", 50, 0, 100, "PT_{e}");
  MSPlot["NbOfLooseElectron"] = new MultiSamplePlot(datasets, "NbOfLooseElectron", 10, 0, 10, "Nb. of loose electrons");

  //Declare arrays of MSPlots
  Int_t minNJets=4, maxNJets=6, minNBJets=2, maxNBJets=4;

  for (Int_t q = minNJets; q <= maxNJets; q++){
    for (Int_t p = minNBJets; p<= maxNBJets; p++){

    string NJets_str = static_cast<ostringstream*>( &(ostringstream() << q) )->str();
    string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << p) )->str();
    string HT_Name = "HT_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string HTX_Name = "HTX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string EventMass_Name = "EventMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string EventMassX_Name = "EventMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

    MSPlot[HT_Name.c_str() ] = new MultiSamplePlot(datasets, HT_Name.c_str() , 35, 0, 1400, "HT");
    MSPlot[HTX_Name.c_str() ] = new MultiSamplePlot(datasets, HTX_Name.c_str() , 35, 0, 1400, "HTX");
    MSPlot[EventMass_Name.c_str() ] = new MultiSamplePlot(datasets, EventMass_Name.c_str() , 35, 0, 1400, "EventMass");
    MSPlot[EventMassX_Name.c_str() ] = new MultiSamplePlot(datasets, EventMassX_Name.c_str() , 35, 0, 1400, "EventMassX");
   
}
}

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
  for (unsigned int d = 0; d < datasets.size(); d++){
    histo2D[("RelIso_vs_MET_"+datasets[d]->Name()).c_str()] = new TH2F(("RelIso_vs_MET_"+datasets[d]->Name()).c_str(),"RelIso:MET",100,0,1000, 100, 0,1);
	
  }
  //  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "FourTop"+postfix+channelpostfix;
  pathPNG += "_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Selection Tables ///////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µ + jets  //////////////////////////
  ////////////////////////////////////////////////////////////////////
  vector<string> CutsSelecTableMu;
  CutsSelecTableMu.push_back(string("initial"));
  //  CutsSelecTableMu.push_back(string("PU reweighting"));
  CutsSelecTableMu.push_back(string("Event cleaning and Trigger"));
  CutsSelecTableMu.push_back(string("Exactly 1 isolated muon"));
  CutsSelecTableMu.push_back(string("Loose muon veto"));
  CutsSelecTableMu.push_back(string("Electron veto"));
  CutsSelecTableMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableMu.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableMu.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableMu.push_back(string("$b\\texttt{-}disc \\geq 0.679$ (CSVM)"));

  SelectionTable selecTableMu(CutsSelecTableMu, datasets);
  selecTableMu.SetLuminosity(Luminosity);
  selecTableMu.SetPrecision(1);


  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : e + jets    ////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableEl;
  CutsSelecTableEl.push_back(string("initial"));
  //  CutsSelecTableEl.push_back(string("PU reweighting"));
  CutsSelecTableEl.push_back(string("Event cleaning and Trigger"));
  CutsSelecTableEl.push_back(string("Exactly 1 isolated electron"));
  CutsSelecTableEl.push_back(string("Loose muon veto"));
  CutsSelecTableEl.push_back(string("Loose electron veto"));
  CutsSelecTableEl.push_back(string("Conversion veto"));
  CutsSelecTableEl.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableEl.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableEl.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableEl.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableEl.push_back(string("$b\\texttt{-}disc \\geq 0.679$ (CSVM)"));
  //CutsSelecTableEl.push_back(string("MET > 40 GeV"));


  SelectionTable selecTableEl(CutsSelecTableEl, datasets);
  selecTableEl.SetLuminosity(Luminosity);
  selecTableEl.SetPrecision(1);
  
  //cout << " - SelectionTable instantiated ..." << endl;


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// PileUp Reweighting - N True interactions, recommended method for 2012 ///////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  //NEW METHOD (TRUE INTERACTIONS)


  LumiReWeighting LumiWeights;
  
  LumiWeights = LumiReWeighting("../macros/PileUpReweighting/pileup_MC_Summer12.root", "../macros/PileUpReweighting/pileup_2012Data_UpToRun191810.root", "pileup", "pileup");

  reweight::PoissonMeanShifter PShiftDown_ = reweight::PoissonMeanShifter(-0.6);
  reweight::PoissonMeanShifter PShiftUp_ = reweight::PoissonMeanShifter(0.6);

  cout << " Initialized LumiReWeighting stuff" << endl;
  
  //exit(1);


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Loop on datasets ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

  for (unsigned int d = 0; d < datasets.size(); d++) //d < datasets.size()
  {
    if (verbose > 1){
      cout << "   Dataset " << d << " name : " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;
      cout << " - Cross section = " << datasets[d]->Xsection() << endl;
      cout << " - IntLumi = " << datasets[d]->EquivalentLumi() << "  NormFactor = " << datasets[d]->NormFactor() << endl;
      cout << " - Nb of events : " << datasets[d]->NofEvtsToRunOver() << endl;
    }
    //open files and load
    cout<<"Load Dataset"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
		
    string previousFilename = "";
    int iFile = -1;
    
    string dataSetName = datasets[d]->Name();	
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////// Initialize JEC factors /////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
    vector<JetCorrectorParameters> vCorrParam;

    // Create the JetCorrectorParameter objects, the order does not matter.
    // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L1FastJet.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L2Relative.txt");
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L3Absolute.txt");

    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);

    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
      {
	JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L2L3Residual.txt");
	vCorrParam.push_back(*ResJetCorPar);
      }
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../macros/JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); // last boolean ('startFromRaw') = false!    


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////                      Loop on events                                                    ///////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    int itrigger = -1, previousRun = -1;
   
    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();

    cout <<"Number of events = "<<  end  <<endl;

    bool debug = false;

    if (verbose > 1) cout << " - Loop over events " << endl;     
    for (unsigned int ievt = start; ievt < end; ievt++)
    {  


	if(ievt%1000 == 0)
		std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

     selecTableMu.Fill(d,0,1.);
     selecTableEl.Fill(d,0,1.);

	//load event
	event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

	vector<TRootGenJet*> genjets;
	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
	  // loading GenJets as I need them for JER
	  		genjets = treeLoader.LoadGenJet(ievt);
	}

	// check which file in the dataset it is to have the HLTInfo right
	string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
	if(previousFilename != currentFilename)
	{
		previousFilename = currentFilename;
        	iFile++;
		cout<<"File changed!!! => iFile = "<<iFile<<endl;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////// trigger /////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool trigged = false;
	int currentRun = event->runId();
	if(previousRun != currentRun)
	{
	  cout <<"What run? "<< currentRun<<endl;
      		previousRun = currentRun;
		if(Muon)
		{
			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
			{
		      
			  if(event->runId() >= 190456 && event->runId() <= 190761)
				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
			  else if ( event->runId() >= 190762 && event->runId() <= 191511 )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v3"), currentRun, iFile);
			  else if ( event->runId() >= 191512  && event->runId() <= 193805  )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v4"), currentRun, iFile);
			  //for  top triggers change PD at this point.
			  else if ( event->runId() >= 193806  && event->runId() <= 194269  )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v1"), currentRun, iFile);
			  else if ( event->runId() >= 194269   )
			    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);

			
  		  if(itrigger == 9999)
				{
    		  cout << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
		    exit(1);
 	 			}
	   		}
	   		else 
	   		{
			  //for refsel, HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2 recommended not HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2
			  if(dataSetName == "TTJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);

			  else  if(dataSetName == "TTTT") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v2"),currentRun, iFile);

			  // 	    if(dataSetName == "TTJets") itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8"), currentRun, iFile);
      			else if (dataSetName == "WJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);

                        else if (dataSetName == "ZJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v3"), currentRun, iFile);

			else if (dataSetName == "SingleTop_t") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);

       			else if (dataSetName == "SingleToptW_T") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
       			else if (dataSetName == "SingleToptW_TBar") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);

				else if (dataSetName == "MultiJet") {

	                 	  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v4"), currentRun, iFile); 
				
				  // if(event->runId() >= 193575 && event->runId() <= 193621)			
				  // if(event->runId() >= 190645 && event->runId() <= 190738)
				  //   itrigger = treeLoader.iTrigger (string ("HLT_Mu20_eta2p1_TriCentralPFJet30_v3"), currentRun, iFile);  
				  //                                  if(event->runId() >= 191057 && event->runId() <= 191411)
				  //  itrigger = treeLoader.iTrigger (string ("HLT_Mu20_eta2p1_TriCentralPFJet30_v4"), currentRun, iFile);

				  //         			  if(event->runId() >= 191695 && event->runId() <= 191810)
				  //  itrigger = treeLoader.iTrigger (string ("HLT_Mu20_eta2p1_TriCentralPFJet30_v5"), currentRun, iFile);    
				  	
				  // if(event->runId() >= 193093 && event->runId() <= 191859)
				  //    itrigger = treeLoader.iTrigger (string ("HLT_Mu20_eta2p1_TriCentralPFJet30_v5"), currentRun, iFile);    
				  }
		    
				       
  				if(itrigger == 9999)
				{
    			  		cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
					//		exit(1);
				}
			}

		} //end if Muon
		else if(Electron)
		{
			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
			{
			
			  if(event->runId() >= 190456 && event->runId() <= 190738)
				 itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8"), currentRun, iFile);
			  else if ( event->runId() >= 190762 && event->runId() <= 191511 )
			         itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v9"), currentRun, iFile);
		  else if ( event->runId() >= 191512  )
			         itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v10"), currentRun, iFile);
	  if(itrigger == 9999)
				{
    		  cout << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
		    exit(1);
 	 			}		
			}
	   		else 
			  {
				if(dataSetName == "TTJets") itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8"), currentRun, iFile);
				else if (dataSetName == "WJets") itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8"), currentRun, iFile);
    
  				if(itrigger == 9999)
				{
    			  		cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
						exit(1);
				}
			}
		} //end if Electron
	} //end previousRun != currentRun

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Jet energy scale corrections     /////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Apply Jet Corrections on-the-fly
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
	if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
	else
		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Type I MET corrections: (Only for |eta| <=4.7 ) //////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");      
	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
		jetTools->correctMETTypeOne(init_jets,mets[0],true);
	else
		jetTools->correctMETTypeOne(init_jets,mets[0],false);
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After MET type I correction:");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Jet energy smearing and systematic uncertainty ///////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
	  //JER 
	  doJERShift == 0;
		if(doJERShift == 1)
			jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
		else if(doJERShift == 2)
			jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
		else
			jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
	  
		//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");	       

		// JES systematic! 
		if (doJESShift == 1)
			jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
		else if (doJESShift == 2)
			jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
	
		//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");
	
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////// Beam scrapping veto and PU reweighting ///////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// scale factor for the event
	float scaleFactor = 1.;

	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	{
		// Apply the scraping veto. (Is it still needed?)
        	bool isBeamBG = true;
        	if(event->nTracks() > 10)
        	{
			if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
			isBeamBG = false;
		}
      		if(isBeamBG) continue;
	}
	else{
	double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
	double lumiWeightOLD=lumiWeight;
	if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
	  lumiWeight=1;
	scaleFactor = scaleFactor*lumiWeight;

	}

	histo1D["lumiWeights"]->Fill(scaleFactor);	
			
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////// Event selection ////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//	selecTableMu.Fill(d,0, 1.);//datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
	//	selecTableEl.Fill(d,1,scaleFactor);
	       
	// Apply trigger selection
	trigged = treeLoader.EventTrigged (itrigger);
	if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;


	//Applying trigger selection again with 2012 Muon+Jets triggers.
	if(!trigged)		   continue;

	// Declare selection instance    
	Selection selection(init_jets, init_muons, init_electrons, mets);

	// Define object selection cuts
	selection.setJetCuts(20.,2.5,0.01,1.,0.98,0,0);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets

	selection.setElectronCuts(30.,2.5,0.1,0.02,0.,999.,0); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets
	selection.setLooseElectronCuts(20,2.5,0.2,0.);

	selection.setMuonCuts(26.,2.1,.12,0,0.2,0,1,0.5,5 ); //Pt,Eta,RelIso,NValidMuHits,d0, dRJets, NMatchedStations,DistVzPVz,NTrackerLayersWithMeas 
	selection.setLooseMuonCuts(10,2.5,0.2);
	  
	//Select objects 
	//	vector<TRootElectron*> selectedElectrons_NoIso = selection.GetSelectedElectrons(20,2.4,999.,vertex[0]);
	vector<TRootElectron*> selectedElectrons       = selection.GetSelectedElectrons(vertex[0]);
	vector<TRootElectron*> selectedExtraElectrons;
	vector<TRootMuon*>     selectedMuons_NoIso = selection.GetSelectedMuons(26,2.4,999.); 
	vector<TRootMuon*>     selectedMuons       = selection.GetSelectedMuons(vertex[0]);
	vector<TRootMuon*>     selectedExtraMuons;
	vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId
    vector<TRootJet*>      selectedSoftJets        = selection.GetSelectedJets(20.,2.5, selectedMuons, 0., true); // ApplyJetId
	vector<TRootMuon*>     selectedLooseMuons     = selection.GetSelectedLooseMuons();
    vector<TRootElectron*> selectedLooseElectrons = selection.GetSelectedLooseElectrons(); // VBTF ID
    vector<TRootJet*>      selectedBJets; // B-Jets
    vector<TRootJet*>      selectedLightJets; // light-Jets

	//order jets wrt to Pt, then set bool corresponding to RefSel cuts.                                                                                 
    sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt Pt.                                                                    

    int JetCut =0;
    int nMu = selectedMuons.size();
    int nEl = selectedElectrons.size();


	// Apply primary vertex selection
	bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
        if(!isGoodPV) continue;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////// Sync'ing cutflow ///////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
        
        
	selecTableMu.Fill(d,1,1.) ;
	selecTableEl.Fill(d,1,1.) ;
	int nTags = 0;
        
        
	  if (nMu ==1){selecTableMu.Fill(d,2,1.) ;
	    if (selectedLooseMuons.size()==1){  selecTableMu.Fill(d,3,1.);
	         if (selectedLooseElectrons.size()==0)  {selecTableMu.Fill(d,4,1.);     
		   if (selectedJets.size() >= 1 && selectedJets[0]->Pt() >45.) {  selecTableMu.Fill(d,5,1.) ;
		     if (selectedJets.size() >= 2 && selectedJets[1]->Pt() >45.) { selecTableMu.Fill(d,6,1.) ;
		       if (selectedJets.size() >= 3 && selectedJets[2]->Pt() >45.){  selecTableMu.Fill(d,7,1.) ;
			 if (selectedJets.size() >= 4 && selectedJets[3]->Pt() >20.) { selecTableMu.Fill(d,8,1.) ;

			   for (int testjet=0; testjet<selectedJets.size(); testjet++) {
			     if (selectedJets[testjet]->btag_combinedSecondaryVertexBJetTags() > 0.679)  nTags++;
			   }
        if (nTags >= 1) selecTableMu.Fill(d,9,1.);
	}
		}
		}
	      }
	    }
	  }
	  }

	  //Sync'ing cutflow -electrons
	   nTags = 0;

	   if (nEl ==1){selecTableEl.Fill(d,2,1.) ;//one isolated electron
	     if (selectedLooseMuons.size()==0){  selecTableEl.Fill(d,3,1.);// loose muon veto
	       if (selectedLooseElectrons.size()==1)  {selecTableEl.Fill(d,4,1.);    // di-lepton veto 
		if( selectedElectrons[0]->passConversion() == true  ){selecTableEl.Fill(d,5,1.); //put conversion rejection cut here
		  if (selectedJets.size() >= 1 && selectedJets[0]->Pt() >45.) {  selecTableEl.Fill(d,6,1.) ; //1 jet 45
		    if (selectedJets.size() >= 2 && selectedJets[1]->Pt() >45.) { selecTableEl.Fill(d,7,1.) ; //2 jet 45
		      if (selectedJets.size() >= 3 && selectedJets[2]->Pt() >45.){  selecTableEl.Fill(d,8,1.) ; //3 jet 45
			if (selectedJets.size() >= 4 && selectedJets[3]->Pt() >20.) { selecTableEl.Fill(d,9,1.) ; //4 jet 20
		      for (int testjet=0; testjet<selectedJets.size(); testjet++) {
			if (selectedJets[testjet]->btag_combinedSecondaryVertexBJetTags() > 0.679)  nTags++;
        }
                      if (nTags >= 1) selecTableEl.Fill(d,10,1.);  //1 CSVM Btag
	}
	}
	}
		}
	      }
	    }
	  }
	  }
        
        
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Applying baseline offline event selection here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (debug)	cout <<" applying baseline event selection..."<<endl;
        
    if (!(selectedJets.size() >= 4)) continue;

    //Apply the staggered Jet Pt cut
	if  ( selectedJets[0]->Pt() >45. &&  selectedJets[1]->Pt() >45. && selectedJets[2]->Pt() >45. && selectedJets[3]->Pt() >20.   ) JetCut = 1;
	if (debug) cout<<" jet1 pt =  "<<selectedJets[0]->Pt() << "   "<< " jet2 pt =  "<<selectedJets[1]->Pt() << "   "<< " jet2 pt =  "<<selectedJets[2]->Pt() << "   "<< " jet3 pt =  "<<selectedJets[3]->Pt() << "  JetCut?"  << JetCut  <<endl;
 
	if (debug)	cout <<" filling bjet vec "<<endl;
	//filling vector of b-jets
	for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ ){
	  if( selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue) selectedBJets.push_back(selectedJets[seljet]);
	  else selectedLightJets.push_back(selectedJets[seljet]);
	}


	if(Muon){

    if (debug) cout <<"Number of Muons, Jets, BJets, JetCut  ===>  "<< selectedMuons.size() <<"  "  << selectedJets.size()   <<"  " <<  selectedBJets.size()   <<"  "<<JetCut  <<endl;


	  //Apply the selection
	  if  (  !( JetCut==1 && selectedBJets.size() >= 2 &&  nMu == 1  )) continue; 

	  if (debug) cout<< "Event passed..."<<endl;

        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////Filling histograms / plotting
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    
      //Order of plotting
	  // 0. Vertices
	  // 1. Muons
	  // 2. Jets: per jet plots, event level variables, jet combinations,discriminants.


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Vertices
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Muons
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        
    //loop to get inv. masses of muon-loosemuon combinations
    for (Int_t mu1 = 0; mu1 < selectedMuons.size(); mu1++){	 
	  for (Int_t mu2 = 0; mu2 < selectedLooseMuons.size(); mu2++){
     	    TRootMuon DiMu (TLorentzVector ( selectedLooseMuons[mu1]->Px() + selectedLooseMuons[mu2]->Px() ,selectedLooseMuons[mu1]->Py() + selectedLooseMuons[mu2]->Py()    ,selectedLooseMuons[mu1]->Pz() + selectedLooseMuons[mu2]->Pz() , selectedLooseMuons[mu1]->E() + selectedLooseMuons[mu2]->E() ));
            MSPlot["DiMuon_InvMass"]->Fill( DiMu.M(), datasets[d], true, Luminosity*scaleFactor );  
	  }
	  }
        
	  MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);

	  //Fill Muon ID plots
      double muzPVz = selectedMuons[0]->vz() - vertex[0]->Z();
      double Mtrans =  (*selectedMuons[0] + *mets[0]).Mt();
	  MSPlot["MuonRelIsolation"]->Fill(selectedMuons[0]->relativePfIso03(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonPt"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonEta"]->Fill(selectedMuons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonPhi"]->Fill(selectedMuons[0]->Phi(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonNValidHits"]->Fill(selectedMuons[0]->nofValidHits(), datasets[d], true, Luminosity*scaleFactor);	 
	  MSPlot["Muond0"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["MuonDistVzPVz"]->Fill(muzPVz, datasets[d], true, Luminosity*scaleFactor );
      MSPlot["MuonNMatchedStations"]->Fill(selectedMuons[0]->nofMatchedStations(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["MuonTrackerLayersWithMeasurement"]->Fill(selectedMuons[0]->nofTrackerLayersWithMeasurement(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["WMt"]->Fill(Mtrans,datasets[d], true, Luminosity*scaleFactor);
      histo2D[("RelIso_vs_MET_"+datasets[d]->Name()).c_str()]->Fill(mets[0]->Et(),selectedMuons_NoIso[0]->relativePfIso03(), Luminosity*scaleFactor );

      if (debug) cout <<"filled all muID plots  .."<<endl;

   
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Jets
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfSelectedBJets"]->Fill(selectedBJets.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["RhoCorrection"]->Fill(event->kt6PFJetsPF2PAT_rho(), datasets[d], true, Luminosity*scaleFactor);
	  if (debug) cout <<"per jet plots.."<<endl;

   	  double ljetpt;
	  double bjetpt;
      double jetpt;
      double HT =0;
	  double MHT =0;
	  double sumpx =0, sumpy=0, sumpz=0, sume=0; 

	//plots to to inspire staggered Jet Pt selection
      if (selectedJets.size()>=4) MSPlot["4thJetPt"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=5) MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=6) MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=7) MSPlot["7thJetPt"]->Fill(selectedJets[6]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	
	for (Int_t seljet1 =0; seljet1 < selectedLightJets.size(); seljet1++ ){
		  MSPlot["SelectedJetPt_light"]->Fill(  selectedLightJets[seljet1]->Pt()  , datasets[d], true, Luminosity*scaleFactor);
		  MSPlot["BdiscBJetCand_CSV"]->Fill(selectedJets[seljet1]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
		  MSPlot["SelectedJetPt"]->Fill(jetpt, datasets[d], true, Luminosity*scaleFactor);
		  MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
		  MSPlot["JetPhi"]->Fill(selectedJets[seljet1]->Phi() , datasets[d], true, Luminosity*scaleFactor);

		  //Event-level variables
          jetpt = selectedJets[seljet1]->Pt();
		  HT = HT + jetpt;
		  sumpx = sumpx + selectedJets[seljet1]->Px();
		  sumpy = sumpy + selectedJets[seljet1]->Py();
		  sumpz = sumpz + selectedJets[seljet1]->Pz();
		  sume = sume + selectedJets[seljet1]->E();
	}

	    TRootJet sumjet (TLorentzVector (sumpx, sumpy, sumpz,sume ));
        MHT =  sumjet.Pt();
     	double MHTSig = MHT/sqrt(HT);
	    double STJet = HT + mets[0]->Et() + selectedMuons[0]->Pt();
        double EventMass = sumjet.M();
        
        MSPlot["MHT_SelectedJets"]->Fill(MHT, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
	    MSPlot["MHTSig_SelectedJets"]->Fill(MHTSig, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MET_MHT"]->Fill(mets[0]->Et()/MHT, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["STJet"]->Fill( STJet , datasets[d], true, Luminosity*scaleFactor);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Jet combinations, Mass Reconstructions...
	//// Logic: 1) Reconstruct best hadronic top candidate (simulataneously minimize Mass_{Topcandidate} - Mass_{Top} && Mass_{Wcandidate} - Mass_{W})
	////        2) Reconstruct best leptonic top candidate without considering b-jet used in hadronic candidate.
	////        3) Reconstruct best hadronic top candidate as in 1) but only use jets not already used.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (debug) cout <<"jet comb."<<endl;
	int firstDiMass = 0;
	int secondDiMass = 0;

	int index_j1_fir, index_j2_fir, index_b_fir, index_j1_sec, index_j2_sec, index_b_sec;
	double Wmass = 80.4;
	double topmass = 176;
	double candchi2 ;
	double firstchi2 =1000;
	double secondchi2 =1000;
	double lepchi2 =1000;

	//These dummy sigmas need to replaced with the fitted sigma of the
	// appropriate dijet and trijet mass distributions.
	double sigW = 25;
	double sigTop = 50;

//////////////////////////////////////////////////////
// First hadronic top candidate
//////////////////////////////////////////////////////
    double fir_top_Pt_chosen;
    double fir_top_mass_chosen;
    double fir_w_mass_chosen;
	for (Int_t jet1 =0; jet1 < selectedLightJets.size(); jet1++ ){
	  for (Int_t jet2 =0; jet2 < selectedLightJets.size(); jet2++ ){
          TRootJet* j1 = (TRootJet*) selectedLightJets[jet1];
          TRootJet* j2 = (TRootJet*) selectedLightJets[jet2];
          double  fir_w_mass = (*j1 + *j2).M();

            for (Int_t bjet1 =0; bjet1 < selectedBJets.size(); bjet1++ ){

                TRootJet* bj1 = (TRootJet*) selectedBJets[bjet1];
                double fir_top_mass  = (*j1 + *j2 + *bj1).M(); 

                double candchi2 =  ( fabs(fir_w_mass  - Wmass) ) + (fabs(fir_top_mass - topmass));
                double fir_top_Pt = (*j1 + *j2 + *bj1).Pt();

                    if ( candchi2 < firstchi2   ){
                        firstchi2 = candchi2;
                        fir_top_mass_chosen = fir_top_mass;
                        fir_w_mass_chosen = fir_w_mass ;
                        fir_top_Pt_chosen = fir_top_Pt;

                        //set indices of best comb
                        index_j1_fir = jet1;
                        index_j2_fir = jet2;
                        index_b_fir = bjet1;	       
	      } 

	    }
	}
	}//end outer dijet loop

//////////////////////////////////////////////////////
// Leptonic top candidate
//////////////////////////////////////////////////////
    double lep_top_Pt_chosen;
    double lep_w_mass_chosen;
    double lep_top_mass_chosen;
        
    MEzCalculator NuPz;
    NuPz.SetMET(*mets[0]);
    NuPz.SetMuon(*selectedMuons[0]);          

    int index_b_lep =0;//index of B-jet in lep top cand.
    double bestlepmass =1000;
       
    TLorentzVector mumet;
    TLorentzVector mumb;

    //loop on Bjets to to make best leptonic top candidate.
	for (Int_t bj1 =0; bj1 < selectedBJets.size(); bj1++ ){

	  //don't consider b-jet from hadronic candidate
	  if (bj1==index_b_fir) continue;

        if (debug) cout <<"in leptonic top bjet loop "<<endl;
	   //Form Lep W
        double   mumpx =   selectedMuons[0]->Px() + mets[0]->Px();
        double   mumpy =   selectedMuons[0]->Py() + mets[0]->Py();
        double   mumpz =   selectedMuons[0]->Pz() + NuPz.Calculate();
        double   mumpe =   selectedMuons[0]->E()  + mets[0]->E();
        mumet.SetPx(mumpx  );
        mumet.SetPy(mumpy  );
        mumet.SetPz(mumpz  );
        mumet.SetE(mumpe  );

        ///Can now fill MuMETB/LepTop mass
        double   mumbpx =   selectedBJets[bj1]->Px() + selectedMuons[0]->Px() + mets[0]->Px();
        double   mumbpy =   selectedJets[bj1]->Py() + selectedMuons[0]->Py() + mets[0]->Py();
        double   mumbpz =   selectedJets[bj1]->Pz() + selectedMuons[0]->Pz() + NuPz.Calculate();
        double   mumbpe =   selectedJets[bj1]->E()  + selectedMuons[0]->E()  + mets[0]->E();
        mumb.SetPx(mumbpx );
        mumb.SetPy(mumbpy );
        mumb.SetPz(mumbpz );
        mumb.SetE(mumbpe );

        double lep_top_mass =mumb.M();
        double lep_top_Pt = mumb.Pt();
        double lep_w_mass =mumet.M();
        double candchi = fabs(lep_top_mass - topmass)  + fabs( lep_w_mass - Wmass);
	   
        if (candchi  < lepchi2 ){
            candchi = lepchi2 ;
            lep_top_mass_chosen = lep_top_mass ;
            lep_w_mass_chosen = lep_w_mass;
            lep_top_Pt_chosen = lep_top_Pt ;
            index_b_lep = bj1;
	   }
	}

	  //////////////////////////////////////////////////////
	  // Second hadronic top candidate
	  //////////////////////////////////////////////////////
	double sec_top_mass_chosen =0 ;
	double sec_top_Pt_chosen = 0 ;
	double sec_w_mass_chosen = 0;

	for (Int_t jet1 =0; jet1 < selectedLightJets.size(); jet1++ ){
	  for (Int_t jet2 =0; jet2 < selectedLightJets.size(); jet2++ ){
          if (debug) cout <<"in second had top  loop "<<endl;

          //disregard any previously used light jets
          if ( jet1 == index_j1_fir ||  jet2 == index_j1_fir ) continue;
          TRootJet* j1 = (TRootJet*) selectedLightJets[jet1];
          TRootJet* j2 = (TRootJet*) selectedLightJets[jet1];
	    
          double sec_w_mass = (*j1 + *j2).M();
	    
          for (Int_t bjet1 =0; bjet1 < selectedBJets.size(); bjet1++ ){

              //disregard any previously used  jets
              if ( bjet1 == index_b_lep ||  bjet1 == index_b_fir ) continue;
              
                TRootJet* bj1 = (TRootJet*) selectedBJets[bjet1];

                double sec_top_mass = (*j1 + *j2 + *bj1).M();	     
                double candchi2 =  ( fabs(sec_w_mass - Wmass)   ) + (fabs(sec_top_mass - topmass));
                double sec_top_Pt = (*j1 + *j2 + *bj1).Pt(); 

                    if ( candchi2 < secondchi2   ){
                        secondchi2 = candchi2;
                        sec_top_mass_chosen = sec_top_mass ;
                        sec_top_Pt_chosen = sec_top_Pt;
                        sec_w_mass_chosen = sec_w_mass;
	
	      } 

	    }
	}
	}//end outer dijet loop

	//leptonic top/W
    MSPlot["MuMetBMasses"]->Fill(lep_top_mass_chosen, datasets[d], true, Luminosity*scaleFactor );
    MSPlot["LepWMass"]->Fill(lep_w_mass_chosen, datasets[d], true, Luminosity*scaleFactor );
	MSPlot["MuMetBPt"]->Fill(lep_top_Pt_chosen,  datasets[d], true, Luminosity*scaleFactor );

	//1st hadronic top/W
	MSPlot["FirstDiJetMass"]->Fill(fir_w_mass_chosen,  datasets[d], true, Luminosity*scaleFactor );
	MSPlot["FirstTriJetMass_1BJet"]->Fill(fir_top_mass_chosen,  datasets[d], true, Luminosity*scaleFactor );
	MSPlot["FirstTriJetPt_1BJet"]->Fill(fir_top_Pt_chosen,  datasets[d], true, Luminosity*scaleFactor );

	//2nd hadronic top/W
	MSPlot["SecondTriJetPt_1BJet"]->Fill(sec_top_Pt_chosen,  datasets[d], true, Luminosity*scaleFactor );
	MSPlot["SecondTriJetMass_1BJet"]->Fill(sec_top_mass_chosen,  datasets[d], true, Luminosity*scaleFactor );
	MSPlot["SecondDiJetMass"]->Fill(sec_w_mass_chosen,  datasets[d], true, Luminosity*scaleFactor );

	if((dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)&& (fabs(sec_w_mass_chosen - Wmass) < 5. )){
	  cout<<"  "<<endl;
          cout <<"Run: " <<event->runId() << "  Lumi block: " << event->lumiBlockId() <<"  Event: " << event->eventId()   <<endl;
	  cout <<"Interesting event!  1st  DiJet Mass = "<< fir_w_mass_chosen << "  1st trijet mass =  "  <<  fir_top_mass_chosen  <<endl;
	  cout <<"Interesting event!  MuMet Mass = "<< lep_w_mass_chosen  << "  MuMetB mass =  "  << lep_top_mass_chosen  <<endl;
	  cout <<"Interesting event!  2nd  DiJet Mass = "<< sec_w_mass_chosen << "  2nd trijet mass =  "  << sec_top_mass_chosen <<endl;
	  cout<<" N Jets: " << selectedJets.size()  <<"  N Tags: " << selectedBJets.size() <<endl;
	  cout<<"  "<<endl;
	}
        
        
//////////////////////////////////////////////////////////////////////////////////
/////          Calculate derived variables for rejecting ttbar + X          //////
/////                                                                       //////
///// 1. HT, 2. HTX, 3.EventMass, 4.EventMassX                              //////
//////////////////////////////////////////////////////////////////////////////////
        
    double HTX = 0, sumpx_X = 0, sumpy_X= 0, sumpz_X =0, sume_X= 0;
    //Calculate HTX: first caculate HT_light
    for (Int_t seljet1 =0; seljet1 < selectedLightJets.size(); seljet1++ ){
            
        //disregard any previously used light jets
        if ( seljet1 == index_j1_fir ||  seljet1 == index_j2_fir  ) continue;
            //Event-level variables
            double ljetpt = selectedLightJets[seljet1]->Pt();
            HTX = HTX + ljetpt;
            sumpx_X = sumpx_X + selectedLightJets[seljet1]->Px();
            sumpy_X = sumpy_X + selectedLightJets[seljet1]->Py();
            sumpz_X = sumpz_X + selectedLightJets[seljet1]->Pz();
            sume_X = sume_X + selectedLightJets[seljet1]->E();
        }
        
        
    //Calculate HTX: second caculate HT_b
    for (Int_t seljet1 =0; seljet1 < selectedBJets.size(); seljet1++ ){
            
        //disregard any previously used b jets
        if ( seljet1 == index_b_fir ||  seljet1 == index_b_lep   ) continue;
            
            double  bjetpt = selectedBJets[seljet1]->Pt();
            HTX = HTX + bjetpt;
            sumpx_X = sumpx_X + selectedBJets[seljet1]->Px();
            sumpy_X = sumpy_X + selectedBJets[seljet1]->Py();
            sumpz_X = sumpz_X + selectedBJets[seljet1]->Pz();
            sume_X =  sume_X  + selectedBJets[seljet1]->E();
        
        }
        
        TRootJet sumjet_X (TLorentzVector (sumpx_X, sumpy_X, sumpz_X,sume_X ));        
        double EventMassX = sumjet_X.M();
    
//////////////////////////////////////////////////////////////////////////////////
/////          Filling NJet Vs NBJet arrays for                             //////
/////                                                                       //////
///// 1. HT, 2. STLep, 3.STJet, 4.HT(w/o ttbar), 5.InvMass (w/o) ttbar      //////
//////////////////////////////////////////////////////////////////////////////////

  for (Int_t b = minNJets; b <= maxNJets; b++){
      for (Int_t c = minNBJets; c<= maxNBJets; c++){
          string NJets_str = static_cast<ostringstream*>( &(ostringstream() << b) )->str();
          string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << c) )->str();
          string HT_Name = "HT_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string HTX_Name = "HTX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string EventMass_Name = "EventMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string EventMassX_Name = "EventMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

        if(b<6 && c<3){
            if(selectedJets.size() == b && selectedBJets.size() == c  ) {
                MSPlot[HT_Name.c_str() ]->Fill(HT,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[HTX_Name.c_str() ]->Fill(HTX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[EventMass_Name.c_str() ]->Fill(EventMass,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[EventMassX_Name.c_str() ]->Fill(EventMassX,datasets[d], true, Luminosity*scaleFactor);
                      }
                        }
        else{
            if(selectedJets.size() >= b && selectedBJets.size() >= c  ) {
                MSPlot[HT_Name.c_str() ]->Fill(HT,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[HTX_Name.c_str() ]->Fill(HTX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[EventMass_Name.c_str() ]->Fill(EventMass,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[EventMassX_Name.c_str() ]->Fill(EventMassX,datasets[d], true, Luminosity*scaleFactor);
            }
        }
}
}

	}
	else if(Electron){
	  //	MSPlot["1stLeadingElectronRelIsolation"]->Fill(selectedElectrons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);

	}

    }//loop on events
    
    //important: free memory
    treeLoader.UnLoadDataset();

    if(jetTools) delete jetTools;
    
  } //loop on datasets
    
  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;

  //Selection tables
  if(Muon){ 
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	selecTableMu.TableCalculator(  false, true, true, true, true);

        cout << "sel tables" << endl;

  //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
  //	selecTableMu.Write(  "FourTop"+postfix+"Table_Mu.tex",    true,true,true,true,false,false,true);
  }
    else if(Electron){
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	selecTableEl.TableCalculator(  false, true, true, true, true);
   //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
    selecTableEl.Write(  "FourTop"+postfix+"Table_El.tex",  true,true,true,true,false,false,true);

  }

 
  fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
	MultiSamplePlot *temp = it->second;
	TH1F *tempHisto_data;
	TH1F *tempHisto_TTTT;
	//	temp->addText("CMS preliminary");
	string name = it->first;
	temp->Draw(false, name, true, true, true, true, true,1,false); // merge TT/QCD/W/Z/ST/
	//Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
	temp->Write(fout, name, true, pathPNG, "pdf");
      
  }

  cout <<"1D  "<< histo1D.size()  <<"2D   "  <<  histo2D.size() <<endl;

  TDirectory* th1dir = fout->mkdir("Histos1D");
  th1dir->cd();
  for(map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {

    
	TH1F *temp = it->second;
	temp->Write();
	//TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	//tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
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

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

