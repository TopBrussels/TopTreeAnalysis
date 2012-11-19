#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
//#include "TKey.h"
//#include "TRandom3.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../../Selection/interface/SelectionTable.h"
#include "../../Content/interface/AnalysisEnvironment.h"
#include "../../Content/interface/Dataset.h"
#include "../../Tools/interface/JetTools.h"
//#include "../../Tools/interface/PlottingTools.h"
#include "../../Tools/interface/MultiSamplePlot.h"
#include "../../Tools/interface/TTreeLoader.h"
#include "../../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../../Reconstruction/interface/JetCorrectorParameters.h"
//#include "../../Reconstruction/interface/JetCorrectionUncertainty.h"
//#include "../../Reconstruction/interface/MakeBinning.h"
#include "../../MCInformation/interface/LumiReWeighting.h"
#include "../interface/TopFCNC_Evt.h"

#include "../../macros/Style.C"

using namespace std;
using namespace reweight;
using namespace TopTree;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

struct HighestJPBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    	return j1->btag_jetProbabilityBJetTags() > j2->btag_jetProbabilityBJetTags();
    }
};
struct HighestCVSBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    	return j1->btag_combinedSecondaryVertexBJetTags() > j2->btag_combinedSecondaryVertexBJetTags();
    }
};

double MuEffSF_Id_Run2012(double eta, double pt);
double MuEffSF_Iso04_Run2012(double eta, double pt);
double MuEffSF_TrgMu8_Run2012(double eta, double pt);
double MuEffSF_TrgMu17_Run2012(double eta, double pt);

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

  /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
  Tagger name  	            WP name WP Discr cut
  TrackCountingHighPur 	    TCHPT 	3.41
  JetProbability 	          JPL 	  0.275
  JetProbability 	          JPM 	  0.545
  JetProbability 	          JPT 	  0.790
  CombinedSecondaryVertex 	CSVL 	  0.244
  CombinedSecondaryVertex 	CSVM 	  0.679
  CombinedSecondaryVertex 	CSVT 	  0.898
  */
  float btagcut     = 0.679;
	float Zmass       = 91.2;
	float Zwindowsize = 30.;
	bool applyAsymmJetPtCut = true;
	float JetPtCuts[4]={50.,30.,20.,20.};

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FCNC search ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
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

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string channelpostfix = "";
  string comments = "_Run2012B";
  string xmlFileName = "";

  bool diElectron = false; // use diElectron channel?
  bool diMuon = true; // use diMuon channel?
  if(diElectron && diMuon){
	cout << "  --> Using both diMuon and diElectron channel? Choose only one (for the moment, since this requires running on different samples/skims)!" << endl;
	exit(1);
  }

  if(diMuon){
	cout << " --> Using the diMuon channel..." << endl;
	channelpostfix = "_DiMuTrigger";
	xmlFileName = "../config/myTopFCNCconfig_Muon.xml";
	//xmlFileName = "../config/TopFCNCconfig_EventSelection_Muon.xml";
  }
  else if(diElectron){
	cout << " --> Using the diElectron channel..." << endl;
	channelpostfix = "_DiElTrigger";
	xmlFileName = "../config/myTopFCNCconfig_Electron.xml";
  }

  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////// AnalysisEnvironment /////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = 2;//anaEnv.Verbose;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  float Luminosity = 5000;
  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    	string dataSetName = datasets[d]->Name();
	if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
	{
		  Luminosity = datasets[d]->EquivalentLumi();
		  break;
	 }
  }
  cout << "Rescaled to an integrated luminosity of "<< Luminosity << endl;

  //Output ROOT file
  string rootFileName ("TopFCNC"+postfix+channelpostfix+comments+"_Plots.root");
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
  Float_t     LeptIdSF = 0;
  Float_t     LeptIsoSF= 0;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ////////////////// MultiSample plots  //////////////////////////////
  ////////////////////////////////////////////////////////////////////

  MSPlot["RhoCorrection"]                     = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
  MSPlot["NbOfVertices"]                      = new MultiSamplePlot(datasets, "NbOfVertices", 30, 0, 30, "Nb. of vertices");

  MSPlot["LeptonDzero"]                       = new MultiSamplePlot(datasets, "LeptonDzero", 500, -0.02, 0.02, "d_{0} [cm]");

  MSPlot["1stLeadingLeptonPt"]                = new MultiSamplePlot(datasets, "1stLeadingLeptonPt", 300, 0, 150, "p_{T} [GeV/c]");
  MSPlot["2ndLeadingLeptonPt"]                = new MultiSamplePlot(datasets, "2ndLeadingLeptonPt", 300, 0, 150, "p_{T} [GeV/c]");

  MSPlot["3rdLeadingMuonPt"]                  = new MultiSamplePlot(datasets, "3rdLeadingMuonPt", 300, 0, 150, "p_{T} [GeV/c]");
  MSPlot["3rdLeadingElectronPt"]              = new MultiSamplePlot(datasets, "3rdLeadingElectronPt", 300, 0, 150, "p_{T} [GeV/c]");

  MSPlot["1stLeadingLeptonRelIsolation"]      = new MultiSamplePlot(datasets, "1stLeadingLeptonRelIsolation", 100, 0, 0.5, "RelIso");
  MSPlot["2ndLeadingLeptonRelIsolation"]      = new MultiSamplePlot(datasets, "2ndLeadingLeptonRelIsolation", 100, 0, 0.5, "RelIso");

  MSPlot["3rdLeadingMuonRelIsolation"]        = new MultiSamplePlot(datasets, "3rdLeadingMuonRelIsolation", 100, 0, 0.5, "RelIso");
  MSPlot["3rdLeadingElectronRelIsolation"]    = new MultiSamplePlot(datasets, "3rdLeadingElectronRelIsolation", 100, 0, 0.5, "RelIso");

  MSPlot["NbOfIsolatedMuons"]                 = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfIsolatedElectrons"]             = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");

  MSPlot["DiLeptonInvMass"]                   = new MultiSamplePlot(datasets, "DiLeptonInvMass", 400, 50, 130, "m_{ll}");
  MSPlot["DiLeptonSystPt"]                    = new MultiSamplePlot(datasets, "DiLeptonSystPt", 400, 0, 400, "p_{T}^{ll} [GeV/c]");
  MSPlot["DiLeptonDR"]                        = new MultiSamplePlot(datasets, "DiLeptonDR", 400, 0, 4, "#Delta R(l^+l^-)");

  MSPlot["NbOfExtraIsolatedMuons"]            = new MultiSamplePlot(datasets, "NbOfExtraIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfExtraIsolatedElectrons"]        = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");

  MSPlot["NbOfSelectedJets_Before3rdLeptCut"] = new MultiSamplePlot(datasets, "NbOfSelectedJets_Before3rdLeptCut", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_mm_ch"]            = new MultiSamplePlot(datasets, "NbOfSelectedJets_mm_ch",  15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_mme_ch"]           = new MultiSamplePlot(datasets, "NbOfSelectedJets_mme_ch", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_mmm_ch"]           = new MultiSamplePlot(datasets, "NbOfSelectedJets_mmm_ch", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_ee_ch"]            = new MultiSamplePlot(datasets, "NbOfSelectedJets_ee_ch",  15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_eee_ch"]           = new MultiSamplePlot(datasets, "NbOfSelectedJets_eee_ch", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedJets_eem_ch"]           = new MultiSamplePlot(datasets, "NbOfSelectedJets_eem_ch", 15, 0, 15, "Nb. of jets");

  MSPlot["FirstLeadingJetPt"]                 = new MultiSamplePlot(datasets, "FirstLeadingJetPt", 100, 0, 500, "Jet p_{T} [GeV/c]");
  MSPlot["SecondLeadingJetPt"]                = new MultiSamplePlot(datasets, "SecondLeadingJetPt", 80, 0, 400, "Jet p_{T} [GeV/c]");
  MSPlot["ThirdLeadingJetPt_mm_ch"]           = new MultiSamplePlot(datasets, "ThirdLeadingJetPt_mm_ch",  60, 0, 300, "Jet p_{T} [GeV/c]");
  MSPlot["ThirdLeadingJetPt_ee_ch"]           = new MultiSamplePlot(datasets, "ThirdLeadingJetPt_ee_ch",  60, 0, 300, "Jet p_{T} [GeV/c]");
  MSPlot["FourthLeadingJetPt_mm_ch"]          = new MultiSamplePlot(datasets, "FourthLeadingJetPt_mm_ch", 50, 0, 250, "Jet p_{T} [GeV/c]");
  MSPlot["FourthLeadingJetPt_ee_ch"]          = new MultiSamplePlot(datasets, "FourthLeadingJetPt_ee_ch", 50, 0, 250, "Jet p_{T} [GeV/c]");

  MSPlot["NbOfVertices_AtLeastFourJets_mm_ch"]= new MultiSamplePlot(datasets, "NbOfVertices_AtLeastFourJets_mm_ch", 30, 0, 30, "Nb. of vertices");
  MSPlot["NbOfVertices_AtLeastTwoJets_mmm_ch"]= new MultiSamplePlot(datasets, "NbOfVertices_AtLeastTwoJets_mmm_ch", 30, 0, 30, "Nb. of vertices");
  MSPlot["NbOfVertices_AtLeastTwoJets_mme_ch"]= new MultiSamplePlot(datasets, "NbOfVertices_AtLeastTwoJets_mme_ch", 30, 0, 30, "Nb. of vertices");
  MSPlot["NbOfVertices_AtLeastFourJets_ee_ch"]= new MultiSamplePlot(datasets, "NbOfVertices_AtLeastFourJets_ee_ch", 30, 0, 30, "Nb. of vertices");
  MSPlot["NbOfVertices_AtLeastTwoJets_eem_ch"]= new MultiSamplePlot(datasets, "NbOfVertices_AtLeastTwoJets_eem_ch", 30, 0, 30, "Nb. of vertices");
  MSPlot["NbOfVertices_AtLeastTwoJets_eee_ch"]= new MultiSamplePlot(datasets, "NbOfVertices_AtLeastTwoJets_eee_ch", 30, 0, 30, "Nb. of vertices");
/*
  MSPlot["BdiscBJetCand_mm_ch_CVS"]           = new MultiSamplePlot(datasets, "BdiscBJetCand_mm_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_ee_ch_CVS"]           = new MultiSamplePlot(datasets, "BdiscBJetCand_ee_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_mm_ch_JP"]            = new MultiSamplePlot(datasets, "BdiscBJetCand_mm_ch_JP",  50, 0, 1, "JP b-disc.");
  MSPlot["BdiscBJetCand_ee_ch_JP"]            = new MultiSamplePlot(datasets, "BdiscBJetCand_ee_ch_JP",  50, 0, 1, "JP b-disc.");

  MSPlot["BdiscBJetCand_mme_ch_CVS"]          = new MultiSamplePlot(datasets, "BdiscBJetCand_mme_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_eee_ch_CVS"]          = new MultiSamplePlot(datasets, "BdiscBJetCand_eee_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_mme_ch_JP"]           = new MultiSamplePlot(datasets, "BdiscBJetCand_mme_ch_JP",  50, 0, 1, "JP b-disc.");
  MSPlot["BdiscBJetCand_eee_ch_JP"]           = new MultiSamplePlot(datasets, "BdiscBJetCand_eee_ch_JP",  50, 0, 1, "JP b-disc.");

  MSPlot["BdiscBJetCand_mmm_ch_CVS"]          = new MultiSamplePlot(datasets, "BdiscBJetCand_mmm_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_eem_ch_CVS"]          = new MultiSamplePlot(datasets, "BdiscBJetCand_eem_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["BdiscBJetCand_mmm_ch_JP"]           = new MultiSamplePlot(datasets, "BdiscBJetCand_mmm_ch_JP",  50, 0, 1, "JP b-disc.");
  MSPlot["BdiscBJetCand_eem_ch_JP"]           = new MultiSamplePlot(datasets, "BdiscBJetCand_eem_ch_JP",  50, 0, 1, "JP b-disc.");
*/
  MSPlot["HighestBdisc_mm_ch_CVS"]            = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_ee_ch_CVS"]            = new MultiSamplePlot(datasets, "HighestBdisc_ee_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_mm_ch_JP"]             = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_JP",  50, 0, 1, "JP b-disc.");
  MSPlot["HighestBdisc_ee_ch_JP"]             = new MultiSamplePlot(datasets, "HighestBdisc_ee_ch_JP",  50, 0, 1, "JP b-disc.");

  MSPlot["HighestBdisc_mmm_ch_CVS"]           = new MultiSamplePlot(datasets, "HighestBdisc_mmm_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_eem_ch_CVS"]           = new MultiSamplePlot(datasets, "HighestBdisc_eem_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_mmm_ch_JP"]            = new MultiSamplePlot(datasets, "HighestBdisc_mmm_ch_JP",  50, 0, 1, "JP b-disc.");
  MSPlot["HighestBdisc_eem_ch_JP"]            = new MultiSamplePlot(datasets, "HighestBdisc_eem_ch_JP",  50, 0, 1, "JP b-disc.");

  MSPlot["HighestBdisc_mme_ch_CVS"]           = new MultiSamplePlot(datasets, "HighestBdisc_mme_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_eee_ch_CVS"]           = new MultiSamplePlot(datasets, "HighestBdisc_eee_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_mme_ch_JP"]            = new MultiSamplePlot(datasets, "HighestBdisc_mme_ch_JP",  50, 0, 1, "JP b-disc.");
  MSPlot["HighestBdisc_eee_ch_JP"]            = new MultiSamplePlot(datasets, "HighestBdisc_eee_ch_JP",  50, 0, 1, "JP b-disc.");
  
  MSPlot["MET_mm_ch"]                         = new MultiSamplePlot(datasets, "MET_mm_ch",  50, 0, 200, "\\slashE_{T} [GeV]");
  MSPlot["MET_mme_ch"]                        = new MultiSamplePlot(datasets, "MET_mme_ch", 50, 0, 200, "\\slashE_{T} [GeV]");
  MSPlot["MET_mmm_ch"]                        = new MultiSamplePlot(datasets, "MET_mmm_ch", 50, 0, 200, "\\slashE_{T} [GeV]");
  MSPlot["MET_ee_ch"]                         = new MultiSamplePlot(datasets, "MET_ee_ch",  50, 0, 200, "\\slashE_{T} [GeV]");
  MSPlot["MET_eee_ch"]                        = new MultiSamplePlot(datasets, "MET_eee_ch", 50, 0, 200, "\\slashE_{T} [GeV]");
  MSPlot["MET_eem_ch"]                        = new MultiSamplePlot(datasets, "MET_eem_ch", 50, 0, 200, "\\slashE_{T} [GeV]");
/*
  MSPlot["Mtt_mm_ch"]                         = new MultiSamplePlot(datasets, "Mtt_mm_ch",  50, 0, 1000, "m_{t#bar{t}} [GeV/c^{2}]");
  MSPlot["Mtt_mme_ch"]                        = new MultiSamplePlot(datasets, "Mtt_mme_ch", 50, 0, 1000, "m_{t#bar{t}} [GeV/c^{2}]");
  MSPlot["Mtt_mmm_ch"]                        = new MultiSamplePlot(datasets, "Mtt_mmm_ch", 50, 0, 1000, "m_{t#bar{t}} [GeV/c^{2}]");
  MSPlot["Mtt_ee_ch"]                         = new MultiSamplePlot(datasets, "Mtt_ee_ch",  50, 0, 1000, "m_{t#bar{t}} [GeV/c^{2}]");
  MSPlot["Mtt_eee_ch"]                        = new MultiSamplePlot(datasets, "Mtt_eee_ch", 50, 0, 1000, "m_{t#bar{t}} [GeV/c^{2}]");
  MSPlot["Mtt_eem_ch"]                        = new MultiSamplePlot(datasets, "Mtt_eem_ch", 50, 0, 1000, "m_{t#bar{t}} [GeV/c^{2}]");

  MSPlot["Mzq_mm_ch"]                         = new MultiSamplePlot(datasets, "Mzq_mm_ch",  50, 80, 300, "m_{Zq} [GeV/c^{2}]");
  MSPlot["Mzq_mme_ch"]                        = new MultiSamplePlot(datasets, "Mzq_mme_ch", 50, 80, 300, "m_{Zq} [GeV/c^{2}]");
  MSPlot["Mzq_mmm_ch"]                        = new MultiSamplePlot(datasets, "Mzq_mmm_ch", 50, 80, 300, "m_{Zq} [GeV/c^{2}]");
  MSPlot["Mzq_ee_ch"]                         = new MultiSamplePlot(datasets, "Mzq_ee_ch",  50, 80, 300, "m_{Zq} [GeV/c^{2}]");
  MSPlot["Mzq_eee_ch"]                        = new MultiSamplePlot(datasets, "Mzq_eee_ch", 50, 80, 300, "m_{Zq} [GeV/c^{2}]");
  MSPlot["Mzq_eem_ch"]                        = new MultiSamplePlot(datasets, "Mzq_eem_ch", 50, 80, 300, "m_{Zq} [GeV/c^{2}]");

  MSPlot["MassChi2_mm_ch"]                    = new MultiSamplePlot(datasets, "MassChi2_mm_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_mmm_ch"]                   = new MultiSamplePlot(datasets, "MassChi2_mmm_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_mme_ch"]                   = new MultiSamplePlot(datasets, "MassChi2_mme_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_ee_ch"]                    = new MultiSamplePlot(datasets, "MassChi2_ee_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_eem_ch"]                   = new MultiSamplePlot(datasets, "MassChi2_eem_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_eee_ch"]                   = new MultiSamplePlot(datasets, "MassChi2_eee_ch", 20, 0, 20, "\\chi^{2}");

  MSPlot["MassChi2_SM_mm_ch"]                 = new MultiSamplePlot(datasets, "MassChi2_SM_mm_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_SM_mmm_ch"]                = new MultiSamplePlot(datasets, "MassChi2_SM_mmm_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_SM_mme_ch"]                = new MultiSamplePlot(datasets, "MassChi2_SM_mme_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_SM_ee_ch"]                 = new MultiSamplePlot(datasets, "MassChi2_SM_ee_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_SM_eem_ch"]                = new MultiSamplePlot(datasets, "MassChi2_SM_eem_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_SM_eee_ch"]                = new MultiSamplePlot(datasets, "MassChi2_SM_eee_ch", 20, 0, 20, "\\chi^{2}");

  MSPlot["MassChi2_FCNC_mm_ch"]               = new MultiSamplePlot(datasets, "MassChi2_FCNC_mm_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_FCNC_mmm_ch"]              = new MultiSamplePlot(datasets, "MassChi2_FCNC_mmm_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_FCNC_mme_ch"]              = new MultiSamplePlot(datasets, "MassChi2_FCNC_mme_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_FCNC_ee_ch"]               = new MultiSamplePlot(datasets, "MassChi2_FCNC_ee_ch",  20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_FCNC_eem_ch"]              = new MultiSamplePlot(datasets, "MassChi2_FCNC_eem_ch", 20, 0, 20, "\\chi^{2}");
  MSPlot["MassChi2_FCNC_eee_ch"]              = new MultiSamplePlot(datasets, "MassChi2_FCNC_eee_ch", 20, 0, 20, "\\chi^{2}");
*/
  MSPlot["TriLeptonInvMass_mmm_ch"]           = new MultiSamplePlot(datasets, "TriLeptonInvMass_mmm_ch", 160, 50, 130, "m_{lll} [GeV/c^{2}]");
  MSPlot["TriLeptonInvMass_mme_ch"]           = new MultiSamplePlot(datasets, "TriLeptonInvMass_mme_ch", 160, 50, 130, "m_{lll} [GeV/c^{2}]");
  MSPlot["TriLeptonInvMass_eem_ch"]           = new MultiSamplePlot(datasets, "TriLeptonInvMass_eem_ch", 160, 50, 130, "m_{lll} [GeV/c^{2}]");
  MSPlot["TriLeptonInvMass_eee_ch"]           = new MultiSamplePlot(datasets, "TriLeptonInvMass_eee_ch", 160, 50, 130, "m_{lll} [GeV/c^{2}]");
//  MSPlot["NbOfLooseMuon"]     = new MultiSamplePlot(datasets, "NbOfLooseMuon", 10, 0, 10, "Nb. of loose muons");
//  MSPlot["NbOfLooseElectron"] = new MultiSamplePlot(datasets, "NbOfLooseElectron", 10, 0, 10, "Nb. of loose electrons");

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",50,0,4);
  for (unsigned int d = 0; d < datasets.size(); d++){
  	histo2D[("d0_vs_phi_1stleadinglepton_"+datasets[d]->Name()).c_str()] = new TH2F(("d0_vs_phi_1stleadinglepton_"+datasets[d]->Name()).c_str(),"d_{0}:#phi",500,-0.02,0.02,500,0,4);
  	histo2D[("d0_vs_phi_2ndleadinglepton_"+datasets[d]->Name()).c_str()] = new TH2F(("d0_vs_phi_2ndleadinglepton_"+datasets[d]->Name()).c_str(),"d_{0}:#phi",500,-0.02,0.02,500,0,4);
  }
/*
  for (unsigned int d = 0; d < datasets.size(); d++){
	histo2D[("MET_vs_Mzq_mm_ch_"+datasets[d]->Name()).c_str()] = new TH2F(("MET_vs_Mzq_mm_ch_"+datasets[d]->Name()).c_str(),"MET:m_{zq}",50,0,200,50,80,300);
	histo2D[("MET_vs_Mzq_mmm_ch_"+datasets[d]->Name()).c_str()] = new TH2F(("MET_vs_Mzq_mmm_ch_"+datasets[d]->Name()).c_str(),"MET:m_{zq}",50,0,200,50,80,300);
	histo2D[("MET_vs_Mzq_mme_ch_"+datasets[d]->Name()).c_str()] = new TH2F(("MET_vs_Mzq_mme_ch_"+datasets[d]->Name()).c_str(),"MET:m_{zq}",50,0,200,50,80,300);
  }
*/
  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "TopFCNC"+postfix+channelpostfix+comments+"_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Selection Tables ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  char zmasscutname[100];
  sprintf(zmasscutname,"$|m_{ll}-m_Z|<%f$ GeV",Zwindowsize);
  char btagcutname[100];
  sprintf(btagcutname,"$b\\texttt{-}disc \\geq %f$ (CSV)",btagcut);

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µµ  ////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  vector<string> CutsSelecTableDiMu;
  CutsSelecTableDiMu.push_back(string("initial"));
  CutsSelecTableDiMu.push_back(string("PU reweighting"));
  CutsSelecTableDiMu.push_back(string("Trigger"));
  CutsSelecTableDiMu.push_back(string("Good PV"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 2 isolated muons"));
  CutsSelecTableDiMu.push_back(string(zmasscutname));//"$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableDiMu.push_back(string("Veto on 3rd iso. lept."));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableDiMu.push_back(string(btagcutname));//"$b\\texttt{-}disc \\geq 0.7$ (CSV)"));
//  CutsSelecTableDiMu.push_back(string("$MET \\leq 30$ GeV"));  

  SelectionTable selecTableDiMu(CutsSelecTableDiMu, datasets);
  selecTableDiMu.SetLuminosity(Luminosity);
  selecTableDiMu.SetPrecision(1);

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µµµ ////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableTriMu;
  CutsSelecTableTriMu.push_back(string("initial"));
  CutsSelecTableTriMu.push_back(string("PU reweighting"));
  CutsSelecTableTriMu.push_back(string("Trigger"));
  CutsSelecTableTriMu.push_back(string("Good PV"));
  CutsSelecTableTriMu.push_back(string("$\\geq$ 2 isolated muons"));
  CutsSelecTableTriMu.push_back(string(zmasscutname));//"$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableTriMu.push_back(string("3rd iso. mu."));
  CutsSelecTableTriMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableTriMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableTriMu.push_back(string("$\\geq$ 3 jet"));
//  CutsSelecTableTriMu.push_back(string("$MET \\leq 30$ GeV"));  

  SelectionTable selecTableTriMu(CutsSelecTableTriMu, datasets);
  selecTableTriMu.SetLuminosity(Luminosity);
  selecTableTriMu.SetPrecision(1);

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µµe ////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsselecTableDiMuEl;
  CutsselecTableDiMuEl.push_back(string("initial"));
  CutsselecTableDiMuEl.push_back(string("PU reweighting"));
  CutsselecTableDiMuEl.push_back(string("Trigger"));
  CutsselecTableDiMuEl.push_back(string("Good PV"));
  CutsselecTableDiMuEl.push_back(string("$\\geq$ 2 isolated muons"));
  CutsselecTableDiMuEl.push_back(string(zmasscutname));//"$|m_{ll}-m_Z|<30$ GeV"));
  CutsselecTableDiMuEl.push_back(string("3rd iso. elec."));
  CutsselecTableDiMuEl.push_back(string("$\\geq$ 1 jet"));
  CutsselecTableDiMuEl.push_back(string("$\\geq$ 2 jet"));
  CutsselecTableDiMuEl.push_back(string("$\\geq$ 3 jet"));
//  CutsselecTableDiMuEl.push_back(string("$MET \\leq 30$ GeV"));  

  SelectionTable selecTableDiMuEl(CutsselecTableDiMuEl, datasets);
  selecTableDiMuEl.SetLuminosity(Luminosity);
  selecTableDiMuEl.SetPrecision(1);

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : ee /////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableDiEl;
  CutsSelecTableDiEl.push_back(string("initial"));
  CutsSelecTableDiEl.push_back(string("PU reweighting"));
  CutsSelecTableDiEl.push_back(string("Trigger"));
  CutsSelecTableDiEl.push_back(string("Good PV"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 2 isolated electrons"));
  CutsSelecTableDiEl.push_back(string(zmasscutname));//"$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableDiEl.push_back(string("Veto on 3rd iso. lept."));
  //CutsSelecTableDiEl.push_back(string("Conversion veto"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableDiEl.push_back(string("$b\\texttt{-}disc \\geq 0.7$ (CSV)"));
  //CutsSelecTableDiEl.push_back(string("MET > 40 GeV"));


  SelectionTable selecTableDiEl(CutsSelecTableDiEl, datasets);
  selecTableDiEl.SetLuminosity(Luminosity);
  selecTableDiEl.SetPrecision(1);
  
  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : eee ////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableTriEl;
  CutsSelecTableTriEl.push_back(string("initial"));
  CutsSelecTableTriEl.push_back(string("PU reweighting"));
  CutsSelecTableTriEl.push_back(string("Trigger"));
  CutsSelecTableTriEl.push_back(string("Good PV"));
  CutsSelecTableTriEl.push_back(string("$\\geq$ 2 isolated electrons"));
  CutsSelecTableTriEl.push_back(string(zmasscutname));//"$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableTriEl.push_back(string("3rd iso. elec."));
  CutsSelecTableTriEl.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableTriEl.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableTriEl.push_back(string("$\\geq$ 3 jet"));
//  CutsSelecTableTriEl.push_back(string("$MET \\leq 30$ GeV"));  

  SelectionTable selecTableTriEl(CutsSelecTableTriEl, datasets);
  selecTableTriEl.SetLuminosity(Luminosity);
  selecTableTriEl.SetPrecision(1);
  

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : eeµ ////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableDiElMu;
  CutsSelecTableDiElMu.push_back(string("initial"));
  CutsSelecTableDiElMu.push_back(string("PU reweighting"));
  CutsSelecTableDiElMu.push_back(string("Trigger"));
  CutsSelecTableDiElMu.push_back(string("Good PV"));
  CutsSelecTableDiElMu.push_back(string("$\\geq$ 2 isolated electrons"));
  CutsSelecTableDiElMu.push_back(string(zmasscutname));//"$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableDiElMu.push_back(string("3rd iso. muon"));
  CutsSelecTableDiElMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableDiElMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableDiElMu.push_back(string("$\\geq$ 3 jet"));
  //CutsSelecTableDiElMu.push_back(string("Conversion veto"));


  SelectionTable selecTableDiElMu(CutsSelecTableDiElMu, datasets);
  selecTableDiElMu.SetLuminosity(Luminosity);
  selecTableDiElMu.SetPrecision(1);

  cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////// PileUp Reweighting ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LumiReWeighting LumiWeights = LumiReWeighting("../pileup/pileup_MC_S10.root", "../pileup/DoubleMu_Run2012A_13Jul_TopTreeID_2013_PileupHistogram.root", "pileup", "pileup");

  cout << " - Initialized LumiReWeighting stuff" << endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Loop on datasets ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
    cout<<"Load Dataset"<<endl;
		
    string previousFilename = "";
    int iFile = -1;
    
    string dataSetName = datasets[d]->Name();	
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////// Initialize JEC factors /////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
    vector<JetCorrectorParameters> vCorrParam;

    // Create the JetCorrectorParameter objects, the order does not matter.
    // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../../macros/JECFiles/Jec11V2_db_AK5PFchs_L1FastJet.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../../macros/JECFiles/Jec11V2_db_AK5PFchs_L2Relative.txt");
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../../macros/JECFiles/Jec11V2_db_AK5PFchs_L3Absolute.txt");

    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);
    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("../../macros/JECFiles/Jec11V2_db_AK5PFchs_L2L3Residual.txt");
      vCorrParam.push_back(*ResJetCorPar);
    }

    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../../macros/JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); // last boolean ('startFromRaw') = false!    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////// Create TTree ///////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    string TTreeFileName ("TopFCNC"+postfix+channelpostfix+comments+"_TTree_"+dataSetName+".root");
    
    cout << "INFO: creating file : "+TTreeFileName << endl;
    
    TFile* TTreeFile = new TFile(TTreeFileName.c_str(),"RECREATE");
        
    TTree* Tree = new TTree("Tree","Tree containing the TopFCNC event candidate");

    TopFCNC_Evt* MyTopFCNC_EvtCand = 0;

    Tree->Branch("TheTopFCNC_Evt","TopFCNC_Evt",&MyTopFCNC_EvtCand);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////// Loop on events //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int itrigger1 = -1, itrigger2 = -1, previousRun = -1;
    int fourIsoLeptCounter = 0;
      
    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();
     
    if (verbose > 1) cout << " - Loop over events " << endl;      
    
    for (unsigned int ievt = start; ievt < end; ievt++)
    {        

      if(ievt%1000 == 0)
		  std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

      //load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

      vector<TRootGenJet*> genjets;
      if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
	    {
		    genjets = treeLoader.LoadGenJet(ievt,false);
	    }
      //cout << "run: " << event->runId() << "  lumi: " << event->lumiBlockId() << "  event: " << event->eventId() << endl;

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
        previousRun = currentRun;
        if(diMuon)
        {
          if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
          {
				    /*------------------------------------------------------------------
				    Dataset : DoubleMu/Run2012A-13Jul2012-v1
				    --------------------------------------------------------------------
            Trigger HLT_Mu13_Mu8_v16 available for runs 190645-193621
            Trigger HLT_Mu17_Mu8_v16 available for runs 190645-193621
				    ------------------------------------------------------------------*/
				    if(currentRun >= 190645 && currentRun <= 193621){
					    itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Mu8_v16"), currentRun, iFile);
					    itrigger2 = treeLoader.iTrigger (string ("HLT_Mu17_TkMu8_v9"), currentRun, iFile);
					  }
				    /*--------------------------------------------------------------------
				    Sub-Total integrated luminosity = 778,2(/pb)
				        Total integrated luminosity = 778,2(/pb)
				    ------------------------------------------------------------------*/

				    /*------------------------------------------------------------------
				    Dataset : DoubleMu/Run2012B-13Jul2012-v4
				    --------------------------------------------------------------------
				    ------------------------------------------------------------------*/
            else if (currentRun >= 193806 && currentRun <= 196027){
              itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Mu8_v17"), currentRun, iFile);
              itrigger2 = treeLoader.iTrigger (string ("HLT_Mu17_TkMu8_v10"), currentRun, iFile);
              // int. lumi = 3401/pb
            }
            else if (currentRun >= 196046 && currentRun <= 196531){
              itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Mu8_v18"), currentRun, iFile);
              itrigger2 = treeLoader.iTrigger (string ("HLT_Mu17_TkMu8_v11"), currentRun, iFile);
              // int. lumi = 939.288/pb
            }
				    /*--------------------------------------------------------------------
				    Sub-Total integrated luminosity = 4340.3(/pb)
				        Total integrated luminosity = 5118.5(/pb)
				    ------------------------------------------------------------------*/
									   
  		      if(itrigger1 == 9999 && itrigger2 == 9999)
				    {
    		      cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
    		      exit(1);
 	 			    }
				    //trigged = treeLoader.EventTrigged (itrigger);
	   		  }
	   		  else{
            itrigger1 = -1;
            itrigger2 = -1;
          }
	   		  /*{
				    if(dataSetName != "ttbar_fcnc"){
				      itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Mu8_v17"), currentRun, iFile);
				      itrigger2 = treeLoader.iTrigger (string ("HLT_Mu17_TkMu8_v10"), currentRun, iFile);
				    }
				    else{
				      itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Mu8_v18"), currentRun, iFile);
				      itrigger2 = treeLoader.iTrigger (string ("HLT_Mu17_TkMu8_v11"), currentRun, iFile);
            }
  				  if(itrigger1 == 9999 || itrigger2 == 9999)
				    {
    			    cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
              exit(1);
            }
            cout<<"Trigger1 bit nr : "<<itrigger1<<endl;
            cout<<"Trigger2 bit nr : "<<itrigger2<<endl;
          }*/
		    } //end if diMuon
		    else if(diElectron)
		    {
		      itrigger2 = -1;
			    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
			    {
				    /*------------------------------------------------------------------
				    Dataset : DoubleElectron/Run2012A-13Jul2012-v1
				    --------------------------------------------------------------------
    				------------------------------------------------------------------*/
				    if(currentRun >=190645  && currentRun <= 190738)
					    itrigger1 = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15", currentRun);
					    // Int. lumi. = 91.063/pb
				    else if(currentRun >= 191043 && currentRun <= 191411)
					    itrigger1 = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16", currentRun);
					    // Int. lumi. = 266.980/pb
				    else if(currentRun >= 191695 && currentRun <= 193621)
					    itrigger1 = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17", currentRun);
					    // Int. lumi. = 347.429/pb
				    /*--------------------------------------------------------------------
				    Sub-Total integrated luminosity = XXX,XXX(/pb)
				        Total integrated luminosity = XXX,XXX(/pb)
				    ------------------------------------------------------------------*/

				    /*--------------------------------------------------------------------
				    Dataset : DoubleElectron/Run2012B-13Jul-v1
				    --------------------------------------------------------------------
    				--------------------------------------------------------------------*/
		    		/*--------------------------------------------------------------------
		    		Sub-Total integrated luminosity =  907,012(/pb)
		    		    Total integrated luminosity = 1108,132(/pb)
		    		------------------------------------------------------------------*/

			    }
          else{
            itrigger1 = -1;
            itrigger2 = -1;
          }
	   		  /*{
				    if(dataSetName != "ttbar_fcnc") itrigger1 = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17"), currentRun, iFile);
				    else itrigger1 = treeLoader.iTrigger (string ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18"), currentRun, iFile);
    
            if(itrigger1 == 9999)
				    {
    			    cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
    			    exit(1);
				    }
				    cout<<"Trigger bit nr : "<<itrigger1<<endl;
			    }*/
		    } //end if diElectron
	    } //end previousRun != currentRun

	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    //////////////////////////////////////////// Jet energy scale corrections     /////////////////////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	    // Apply Jet Corrections on-the-fly
	    //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
/*
	    if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
		    jetTools->correctJets(init_jets,event->kt6PFJets_rho(),true); //last boolean: isData (needed for L2L3Residual...)
	    else
		    jetTools->correctJets(init_jets,event->kt6PFJets_rho(),false); //last boolean: isData (needed for L2L3Residual...)
	    //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");

*/
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    //////////////////////////////////////////// Type I MET corrections: (Only for |eta| <=4.7 ) //////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	    //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");      
	    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
		    jetTools->correctMETTypeOne(init_jets,mets[0],true);
	    else
		    jetTools->correctMETTypeOne(init_jets,mets[0],false);
	    //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After MET type I correction:");
*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Jet energy smearing and systematic uncertainty ///////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
  if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
	{
		// JER systematic! 
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
/*
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
*/
	if(dataSetName != "Data" && dataSetName != "data" && dataSetName != "DATA")
	{
		// Apply pile-up reweighting
    double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
   	scaleFactor *= lumiWeight;
	}

	histo1D["lumiWeights"]->Fill(scaleFactor);	
	MSPlot["RhoCorrection"]->Fill(event->kt6PFJets_rho(), datasets[d], true, Luminosity*scaleFactor);
			
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////// Event selection ////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	selecTableDiMu.Fill(d,0, 1.);//datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
	selecTableTriMu.Fill(d,0, 1.);
	selecTableDiMuEl.Fill(d,0, 1.);
	selecTableDiEl.Fill(d,0, 1.);
	selecTableTriEl.Fill(d,0, 1.);
	selecTableDiElMu.Fill(d,0, 1.);

	selecTableDiMu.Fill(d,1,scaleFactor);
	selecTableTriMu.Fill(d,1,scaleFactor);
	selecTableDiMuEl.Fill(d,1,scaleFactor);
	selecTableDiEl.Fill(d,1,scaleFactor);
	selecTableTriEl.Fill(d,1,scaleFactor);
	selecTableDiElMu.Fill(d,1,scaleFactor);
		
	// Apply trigger selection
      if (itrigger1 == -1 && itrigger2 == -1) {
        trigged = true;
      }
      else {
        trigged = (treeLoader.EventTrigged (itrigger1) || treeLoader.EventTrigged (itrigger2));
      }

	if(!trigged) continue;
	selecTableDiMu.Fill(d,2,scaleFactor);
	selecTableTriMu.Fill(d,2,scaleFactor);
	selecTableDiMuEl.Fill(d,2,scaleFactor);
	selecTableDiEl.Fill(d,2,scaleFactor);
	selecTableTriEl.Fill(d,2,scaleFactor);
	selecTableDiElMu.Fill(d,2,scaleFactor);

	// Declare selection instance    
	Selection selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...
/*            
	vector<TRootMCParticle*> mcParticles;
	if(dataSetName.find("ttbar_fcnc") == 0)
	{
		treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
		sort(mcParticles.begin(),mcParticles.end(),HighestPt());
	}
*/
	// Define object selection cuts
	selection.setJetCuts(20.,2.5,0.01,1.,0.98,0.3,0.1);

	selection.setDiElectronCuts(20,2.5,0.15,0.04,0.,1,0.3,1); //Et,Eta,RelIso,d0,MVAId,DistVzPVz,DRJets,MaxMissHits
	//selection.setLooseElectronCuts(15,2.5,0.2);

	selection.setDiMuonCuts(20,2.4,0.20,999.); //Et,Eta,RelIso,d0
	//selection.setLooseMuonCuts(15,2.4,0.2);
	  
	//Select objects 
	vector<TRootElectron*> selectedElectrons_NoIso = selection.GetSelectedDiElectrons(20,2.5,999.); //Et,Eta,RelIso
	vector<TRootElectron*> selectedElectrons       = selection.GetSelectedDiElectrons();
	vector<TRootElectron*> selectedExtraElectrons;

	vector<TRootMuon*>     selectedMuons_NoIso = selection.GetSelectedDiMuons(20,2.4,999.); //Et,Eta,RelIso
	vector<TRootMuon*>     selectedMuons       = selection.GetSelectedDiMuons();
	vector<TRootMuon*>     selectedExtraMuons;

	vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId

	//vector<TRootMuon*>     looseMuons     = selection.GetSelectedLooseMuons();
	//vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseElectrons(true); // VBTF Id

	// Apply primary vertex selection
	bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
  if(!isGoodPV) continue;
	selecTableDiMu.Fill(d,3,scaleFactor);
	selecTableTriMu.Fill(d,3,scaleFactor);
	selecTableDiMuEl.Fill(d,3,scaleFactor);
	selecTableDiEl.Fill(d,3,scaleFactor);
	selecTableTriEl.Fill(d,3,scaleFactor);
	selecTableDiElMu.Fill(d,3,scaleFactor);

	// Select events with at least two leptons and apply muon Id/Iso/Trg SF
	if(diMuon){
		if(selectedMuons_NoIso.size()<2) continue;
    scaleFactor *= MuEffSF_Id_Run2012(selectedMuons[0]->Eta(), selectedMuons[0]->Pt());// Id SF
    scaleFactor *= MuEffSF_Id_Run2012(selectedMuons[1]->Eta(), selectedMuons[1]->Pt());// Id SF

    scaleFactor *= MuEffSF_Iso04_Run2012(selectedMuons[0]->Eta(), selectedMuons[0]->Pt());// Iso SF
    scaleFactor *= MuEffSF_Iso04_Run2012(selectedMuons[1]->Eta(), selectedMuons[1]->Pt());// Iso SF

    scaleFactor *= MuEffSF_TrgMu17_Run2012(selectedMuons[0]->Eta(), selectedMuons[0]->Pt());// Trigger SF
    scaleFactor *= MuEffSF_TrgMu8_Run2012( selectedMuons[1]->Eta(), selectedMuons[1]->Pt());// Trigger SF
	}
	else if(diElectron){
		if(selectedElectrons_NoIso.size()<2) continue;
	}

	MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

	if(diMuon){
	  histo2D[("d0_vs_phi_1stleadinglepton_"+datasets[d]->Name()).c_str()]->Fill(selectedMuons_NoIso[0]->d0(),selectedMuons_NoIso[0]->Phi());
	  histo2D[("d0_vs_phi_2ndleadinglepton_"+datasets[d]->Name()).c_str()]->Fill(selectedMuons_NoIso[1]->d0(),selectedMuons_NoIso[1]->Phi());
	  for(unsigned int i=0;i<selectedMuons_NoIso.size();i++) MSPlot["LeptonDzero"]->Fill(selectedMuons_NoIso[i]->d0(), datasets[d], true, Luminosity*scaleFactor);
		MSPlot["1stLeadingLeptonPt"]->Fill(selectedMuons_NoIso[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		MSPlot["2ndLeadingLeptonPt"]->Fill(selectedMuons_NoIso[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		MSPlot["1stLeadingLeptonRelIsolation"]->Fill(selectedMuons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
		MSPlot["2ndLeadingLeptonRelIsolation"]->Fill(selectedMuons_NoIso[1]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
	}
	else if(diElectron){
	  histo2D[("d0_vs_phi_1stleadinglepton_"+datasets[d]->Name()).c_str()]->Fill(selectedElectrons_NoIso[0]->d0(),selectedElectrons_NoIso[0]->Phi());
	  histo2D[("d0_vs_phi_2ndleadinglepton_"+datasets[d]->Name()).c_str()]->Fill(selectedElectrons_NoIso[1]->d0(),selectedElectrons_NoIso[1]->Phi());
		MSPlot["1stLeadingLeptonPt"]->Fill(selectedElectrons_NoIso[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		MSPlot["2ndLeadingLeptonPt"]->Fill(selectedElectrons_NoIso[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		MSPlot["1stLeadingLeptonRelIsolation"]->Fill(selectedElectrons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
		MSPlot["2ndLeadingLeptonRelIsolation"]->Fill(selectedElectrons_NoIso[1]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
	}

	MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);

	// Select events with at least two isolated muons
	if(diMuon){
		if(selectedMuons.size()<2) continue;
	}
	else if(diElectron){
		if(selectedElectrons.size()<2) continue;
	}

	selecTableDiMu.Fill(d,4,scaleFactor);
	selecTableTriMu.Fill(d,4, scaleFactor);
	selecTableDiMuEl.Fill(d,4, scaleFactor);
	selecTableDiEl.Fill(d,4,scaleFactor);
	selecTableTriEl.Fill(d,4,scaleFactor);
	selecTableDiElMu.Fill(d,4,scaleFactor);

  bool foundZ = false;
	int idx_Z_1 = -1, idx_Z_2 = -1;
	// Calculate the invariant mass for each isolated lepton pairs
	// - return true if the mass is the Z boson mass window 
	// - return the indices of the lepton candidates
	if(diMuon){
    for(unsigned int i=0;i<selectedMuons.size()-1;i++)
    {
      for(unsigned int j=i+1;j<selectedMuons.size();j++)
      {
        TRootMuon* mu1 = (TRootMuon*) selectedMuons[i];
        TRootMuon* mu2 = (TRootMuon*) selectedMuons[j];
        if(mu1->charge() == mu2->charge()) continue;
        double invMass = (*mu1 + *mu2).M();
        MSPlot["DiLeptonInvMass"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
        if( invMass >= (Zmass-Zwindowsize) && invMass <= (Zmass+Zwindowsize) )
        {
          idx_Z_1 = i;
          idx_Z_2 = j;
          foundZ = true;
        }
      }
    }
	}
	else if(diElectron){
	  	for(unsigned int i=0;i<selectedElectrons.size()-1;i++)
	  	{
	  		for(unsigned int j=i+1;j<selectedElectrons.size();j++)
	  		{
	   			TRootElectron* el1 = (TRootElectron*) selectedElectrons[i];
	   			TRootElectron* el2 = (TRootElectron*) selectedElectrons[j];
				if(el1->charge() == el2->charge()) continue;
				double invMass = (*el1 + *el2).M();
				MSPlot["DiLeptonInvMass"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
				if( invMass >= (Zmass-Zwindowsize) && invMass <= (Zmass+Zwindowsize) )
				{
					idx_Z_1 = i;
					idx_Z_2 = j;
					foundZ = true;
				}
	    		}
	  	}
	}
      // Select events with at least one pair of opposite charge leptons with |mll-mz|<windowsize
      if(!foundZ)
        continue;
      
      float DiLeptonSystPt = 0;
      float DiLeptonDR = 4;
      if(diMuon){
        DiLeptonSystPt = (*selectedMuons[idx_Z_1] + *selectedMuons[idx_Z_2]).Pt();
        DiLeptonDR     =   selectedMuons[idx_Z_1]->DeltaR(*selectedMuons[idx_Z_2])
      }
      else if(diElectron){
        DiLeptonSystPt = (*selectedElectrons[idx_Z_1] + *selectedElectrons[idx_Z_2]).Pt();
        DiLeptonDR     =   selectedElectrons[idx_Z_1]->DeltaR(*selectedElectrons[idx_Z_2])
      }
      MSPlot["DiLeptonSystPt"]->Fill(DiLeptonSystPt, datasets[d], true, Luminosity*scaleFactor);
      MSPlot["DiLeptonDR"]    ->Fill(DiLeptonDR,     datasets[d], true, Luminosity*scaleFactor);

      selecTableDiMu.Fill(d,5,scaleFactor);
      selecTableTriMu.Fill(d,5, scaleFactor);
      selecTableDiMuEl.Fill(d,5, scaleFactor);
      selecTableDiEl.Fill(d,5,scaleFactor);
      selecTableTriEl.Fill(d,5,scaleFactor);
      selecTableDiElMu.Fill(d,5,scaleFactor);

      selectedExtraMuons = selectedMuons;
      selectedExtraElectrons = selectedElectrons;
      // Erase Z boson lepton candidates
	if(diMuon){
		if(selectedMuons_NoIso.size()>2){
		  MSPlot["3rdLeadingMuonPt"]->Fill(selectedMuons_NoIso[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["3rdLeadingMuonRelIsolation"]->Fill(selectedMuons_NoIso[2]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
    }
		if(selectedElectrons_NoIso.size()>0){
		  MSPlot["3rdLeadingElectronPt"]->Fill(selectedElectrons_NoIso[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["3rdLeadingElectronRelIsolation"]->Fill(selectedElectrons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
    }
		selectedExtraMuons.erase(selectedExtraMuons.begin()+idx_Z_2);
		selectedExtraMuons.erase(selectedExtraMuons.begin()+idx_Z_1);
	}
	else if(diElectron){
		if(selectedMuons_NoIso.size()>0){
		  MSPlot["3rdLeadingMuonPt"]->Fill(selectedMuons_NoIso[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["3rdLeadingMuonRelIsolation"]->Fill(selectedMuons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
    }
		if(selectedElectrons_NoIso.size()>2){
		  MSPlot["3rdLeadingElectronPt"]->Fill(selectedElectrons_NoIso[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["3rdLeadingElectronRelIsolation"]->Fill(selectedElectrons_NoIso[2]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
    }
		selectedExtraElectrons.erase(selectedExtraElectrons.begin()+idx_Z_2);
		selectedExtraElectrons.erase(selectedExtraElectrons.begin()+idx_Z_1);

	}

	MSPlot["NbOfExtraIsolatedMuons"]->Fill(selectedExtraMuons.size(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["NbOfExtraIsolatedElectrons"]->Fill(selectedExtraElectrons.size(), datasets[d], true, Luminosity*scaleFactor);

	MSPlot["NbOfSelectedJets_Before3rdLeptCut"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);

	MyTopFCNC_EvtCand = 0;
	double invMass = 0;
	double highestbtagdisc = 0;
	// Select events based on the presence of *exactly one* extra isolated lepton
	if(diMuon && selectedExtraMuons.size()==0 && selectedElectrons.size()==0){
		selecTableDiMu.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_mm_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
		  MSPlot["FirstLeadingJetPt"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		  if(applyAsymmJetPtCut && selectedJets[0]->Pt()<JetPtCuts[0]) continue;
			selecTableDiMu.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
		    MSPlot["SecondLeadingJetPt"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		    if(applyAsymmJetPtCut && selectedJets[1]->Pt()<JetPtCuts[1]) continue;
				selecTableDiMu.Fill(d,8,scaleFactor);
				if(selectedJets.size()>2){ //at least 3 jets
   		    MSPlot["ThirdLeadingJetPt_mm_ch"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
					selecTableDiMu.Fill(d,9,scaleFactor);
					if(selectedJets.size()>3){ //at least 4 jets
     		    MSPlot["FourthLeadingJetPt_mm_ch"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
     		    MSPlot["NbOfVertices_AtLeastFourJets_mm_ch"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
						selecTableDiMu.Fill(d,10,scaleFactor);

						sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
						MSPlot["HighestBdisc_mm_ch_CVS"]->Fill(selectedJets[0]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
						highestbtagdisc = selectedJets[0]->btag_combinedSecondaryVertexBJetTags();
						sort(selectedJets.begin(),selectedJets.end(),HighestJPBtag());
						MSPlot["HighestBdisc_mm_ch_JP"]->Fill(selectedJets[0]->btag_jetProbabilityBJetTags(),datasets[d], true, Luminosity*scaleFactor);

            // Create TopFCNC_Evt object
						MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kMuon);
						MyTopFCNC_EvtCand->SetLepton1FromZ(*selectedMuons[idx_Z_1]);
						MyTopFCNC_EvtCand->SetLepton2FromZ(*selectedMuons[idx_Z_2]);
						MyTopFCNC_EvtCand->SetSelectedJets(selectedJets);
						MyTopFCNC_EvtCand->SetMET(*mets[0]);

						if(highestbtagdisc>btagcut){
  						selecTableDiMu.Fill(d,11,scaleFactor);
  						MSPlot["MET_mm_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
  				  }
					}
				}
			}
		}
	}
	else if(diMuon && selectedExtraMuons.size()==0 && selectedElectrons.size()==1){
		selecTableDiMuEl.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_mme_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
		  MSPlot["FirstLeadingJetPt"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		  if(applyAsymmJetPtCut && selectedJets[0]->Pt()<JetPtCuts[0]) continue;
			selecTableDiMuEl.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
		    MSPlot["SecondLeadingJetPt"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		    if(applyAsymmJetPtCut && selectedJets[1]->Pt()<JetPtCuts[1]) continue;
        MSPlot["NbOfVertices_AtLeastTwoJets_mme_ch"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
				selecTableDiMuEl.Fill(d,8,scaleFactor);
  			if(selectedJets.size()>2) selecTableDiMuEl.Fill(d,9,scaleFactor);//at least 3 jets

				sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
				MSPlot["HighestBdisc_mme_ch_CVS"]->Fill(selectedJets[0]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
				sort(selectedJets.begin(),selectedJets.end(),HighestJPBtag());
				MSPlot["HighestBdisc_mme_ch_JP"]->Fill(selectedJets[0]->btag_jetProbabilityBJetTags(),datasets[d], true, Luminosity*scaleFactor);

        // Create TopFCNC_Evt object
  			MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kMuon,TopFCNC_Evt::kElec);
				MyTopFCNC_EvtCand->SetLepton1FromZ(*selectedMuons[idx_Z_1]); 
				MyTopFCNC_EvtCand->SetLepton2FromZ(*selectedMuons[idx_Z_2]);
				MyTopFCNC_EvtCand->SetLeptonFromW(*selectedElectrons[0]);
				MyTopFCNC_EvtCand->SetSelectedJets(selectedJets), 
				MyTopFCNC_EvtCand->SetNeutrino(*mets[0]);
				MyTopFCNC_EvtCand->SetMET(*mets[0]);

				MSPlot["MET_mme_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);

				invMass = (*selectedMuons[idx_Z_1]+*selectedMuons[idx_Z_2]+*selectedElectrons[0]).M();
				MSPlot["TriLeptonInvMass_mme_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
			}
		}
	}
	else if(diMuon && selectedExtraMuons.size()==1 && selectedElectrons.size()==0){
		selecTableTriMu.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_mmm_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
		  MSPlot["FirstLeadingJetPt"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		  if(applyAsymmJetPtCut && selectedJets[0]->Pt()<JetPtCuts[0]) continue;
			selecTableTriMu.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
		    MSPlot["SecondLeadingJetPt"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		    if(applyAsymmJetPtCut && selectedJets[1]->Pt()<JetPtCuts[1]) continue;
		    MSPlot["NbOfVertices_AtLeastTwoJets_mmm_ch"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
				selecTableTriMu.Fill(d,8,scaleFactor);
  			if(selectedJets.size()>2) selecTableTriMu.Fill(d,9,scaleFactor); //at least 3 jets

				sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
				MSPlot["HighestBdisc_mmm_ch_CVS"]->Fill(selectedJets[0]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
				sort(selectedJets.begin(),selectedJets.end(),HighestJPBtag());
				MSPlot["HighestBdisc_mmm_ch_JP"]->Fill(selectedJets[0]->btag_jetProbabilityBJetTags(),datasets[d], true, Luminosity*scaleFactor);

        // Create TopFCNC_Evt object
				MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kMuon,TopFCNC_Evt::kMuon);
				MyTopFCNC_EvtCand->SetLepton1FromZ(*selectedMuons[idx_Z_1]); 
				MyTopFCNC_EvtCand->SetLepton2FromZ(*selectedMuons[idx_Z_2]);
				MyTopFCNC_EvtCand->SetLeptonFromW(*selectedExtraMuons[0]);
				MyTopFCNC_EvtCand->SetSelectedJets(selectedJets), 
				MyTopFCNC_EvtCand->SetNeutrino(*mets[0]);
				MyTopFCNC_EvtCand->SetMET(*mets[0]);

				MSPlot["MET_mmm_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);

				invMass = (*selectedMuons[idx_Z_1]+*selectedMuons[idx_Z_2]+*selectedExtraMuons[0]).M();
				MSPlot["TriLeptonInvMass_mmm_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
			}
		}
	}
	else if(diElectron && selectedExtraElectrons.size()==0 && selectedMuons.size()==0){
		selecTableDiEl.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_ee_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
		  MSPlot["FirstLeadingJetPt"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		  if(applyAsymmJetPtCut && selectedJets[0]->Pt()<JetPtCuts[0]) continue;
			selecTableDiEl.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
		    MSPlot["SecondLeadingJetPt"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		    if(applyAsymmJetPtCut && selectedJets[1]->Pt()<JetPtCuts[1]) continue;
				selecTableDiEl.Fill(d,8,scaleFactor);
				if(selectedJets.size()>2){ //at least 3 jets
		      MSPlot["ThirdLeadingJetPt_ee_ch"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
					selecTableDiEl.Fill(d,9,scaleFactor);
					if(selectedJets.size()>3){ //at least 4 jets
		        MSPlot["FourthLeadingJetPt_ee_ch"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
     		    MSPlot["NbOfVertices_AtLeastFourJets_ee_ch"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
						selecTableDiEl.Fill(d,10,scaleFactor);

						sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
						MSPlot["HighestBdisc_ee_ch_CVS"]->Fill(selectedJets[0]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
						highestbtagdisc = selectedJets[0]->btag_combinedSecondaryVertexBJetTags();
						sort(selectedJets.begin(),selectedJets.end(),HighestJPBtag());
						MSPlot["HighestBdisc_ee_ch_JP"]->Fill(selectedJets[0]->btag_jetProbabilityBJetTags(),datasets[d], true, Luminosity*scaleFactor);

            // Create TopFCNC_Evt object
						MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kElec);
						MyTopFCNC_EvtCand->SetLepton1FromZ(*selectedElectrons[idx_Z_1]);
						MyTopFCNC_EvtCand->SetLepton2FromZ(*selectedElectrons[idx_Z_2]);
						MyTopFCNC_EvtCand->SetSelectedJets(selectedJets);
						MyTopFCNC_EvtCand->SetMET(*mets[0]);

						MSPlot["MET_ee_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
					}
				}
			}
		}
	}
	else if(diElectron && selectedExtraElectrons.size()==0 && selectedMuons.size()==1){
		selecTableDiElMu.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_eem_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
		  if(applyAsymmJetPtCut && selectedJets[0]->Pt()<JetPtCuts[0]) continue;
		  MSPlot["FirstLeadingJetPt"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			selecTableDiElMu.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
		    MSPlot["SecondLeadingJetPt"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		    if(applyAsymmJetPtCut && selectedJets[1]->Pt()<JetPtCuts[1]) continue;
		    MSPlot["NbOfVertices_AtLeastTwoJets_eem_ch"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
				selecTableDiElMu.Fill(d,8,scaleFactor);
				if(selectedJets.size()>2) selecTableDiElMu.Fill(d,9,scaleFactor); //at least 3 jets

				sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
				MSPlot["HighestBdisc_eem_ch_CVS"]->Fill(selectedJets[0]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
				sort(selectedJets.begin(),selectedJets.end(),HighestJPBtag());
				MSPlot["HighestBdisc_eem_ch_JP"]->Fill(selectedJets[0]->btag_jetProbabilityBJetTags(),datasets[d], true, Luminosity*scaleFactor);

        // Create TopFCNC_Evt object
  			MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kElec,TopFCNC_Evt::kMuon);
				MyTopFCNC_EvtCand->SetLepton1FromZ(*selectedElectrons[idx_Z_1]); 
				MyTopFCNC_EvtCand->SetLepton2FromZ(*selectedElectrons[idx_Z_2]);
				MyTopFCNC_EvtCand->SetLeptonFromW(*selectedMuons[0]);
				MyTopFCNC_EvtCand->SetSelectedJets(selectedJets), 
				MyTopFCNC_EvtCand->SetNeutrino(*mets[0]);
				MyTopFCNC_EvtCand->SetMET(*mets[0]);

				MSPlot["MET_eem_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);

				invMass = (*selectedElectrons[idx_Z_1]+*selectedElectrons[idx_Z_2]+*selectedMuons[0]).M();
				MSPlot["TriLeptonInvMass_eem_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
			}
		}
	}
	else if(diElectron && selectedExtraElectrons.size()==1 && selectedMuons.size()==0){
		selecTableTriEl.Fill(d,6,scaleFactor);
		MSPlot["NbOfSelectedJets_eee_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		if(selectedJets.size()>0){ //at least 1 jet
		  MSPlot["FirstLeadingJetPt"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		  if(applyAsymmJetPtCut && selectedJets[0]->Pt()<JetPtCuts[0]) continue;
			selecTableTriEl.Fill(d,7,scaleFactor); 
			if(selectedJets.size()>1){ //at least 2 jets
		    MSPlot["SecondLeadingJetPt"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
		    if(applyAsymmJetPtCut && selectedJets[1]->Pt()<JetPtCuts[1]) continue;
		    MSPlot["NbOfVertices_AtLeastTwoJets_eee_ch"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
				selecTableTriEl.Fill(d,8,scaleFactor);
				if(selectedJets.size()>2) selecTableTriEl.Fill(d,9,scaleFactor); //at least 3 jets

				sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
				MSPlot["HighestBdisc_eee_ch_CVS"]->Fill(selectedJets[0]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
				sort(selectedJets.begin(),selectedJets.end(),HighestJPBtag());
				MSPlot["HighestBdisc_eee_ch_JP"]->Fill(selectedJets[0]->btag_jetProbabilityBJetTags(),datasets[d], true, Luminosity*scaleFactor);

        // Create TopFCNC_Evt object
				MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kElec,TopFCNC_Evt::kElec);
				MyTopFCNC_EvtCand->SetLepton1FromZ(*selectedElectrons[idx_Z_1]); 
				MyTopFCNC_EvtCand->SetLepton2FromZ(*selectedElectrons[idx_Z_2]);
				MyTopFCNC_EvtCand->SetLeptonFromW(*selectedExtraElectrons[0]);
				MyTopFCNC_EvtCand->SetSelectedJets(selectedJets), 
				MyTopFCNC_EvtCand->SetNeutrino(*mets[0]);
				MyTopFCNC_EvtCand->SetMET(*mets[0]);

				MSPlot["MET_eee_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);

				invMass = (*selectedElectrons[idx_Z_1]+*selectedElectrons[idx_Z_2]+*selectedExtraElectrons[0]).M();
				MSPlot["TriLeptonInvMass_eee_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
			}
		}
	}
	else fourIsoLeptCounter++;
  if(MyTopFCNC_EvtCand){
    MyTopFCNC_EvtCand->SetEventID( event->eventId() );
    MyTopFCNC_EvtCand->SetRunID( event->runId() );
    MyTopFCNC_EvtCand->SetLumiBlockID( event->lumiBlockId() );
    MyTopFCNC_EvtCand->SetIdParton1( event->idParton1() );
    MyTopFCNC_EvtCand->SetIdParton2( event->idParton2() );
    MyTopFCNC_EvtCand->SetxParton1( event->xParton1() );
    MyTopFCNC_EvtCand->SetxParton2( event->xParton2() );
    MyTopFCNC_EvtCand->SetFactorizationScale( event->factorizationScale() );
    MyTopFCNC_EvtCand->SetnPV(vertex.size());
    MyTopFCNC_EvtCand->SetnTruePU(event->nTruePU());
    MyTopFCNC_EvtCand->SetEventWeight(scaleFactor);
  }
  Tree->Fill();

	//delete selection;
	if(MyTopFCNC_EvtCand) delete MyTopFCNC_EvtCand;
  }//loop on events
    
    cout<<endl;
    cout<<"FYI ; nb of events with at least four isolated leptons = "<<fourIsoLeptCounter<<endl;

    TTreeFile->cd();
    
    TTree *configTree = new TTree("configTree","configuration Tree");
    TClonesArray* tcdataset = new TClonesArray("Dataset",1);
    configTree->Branch("Dataset","TClonesArray",&tcdataset);
    new ((*tcdataset)[0]) Dataset(*datasets[d]);
    
    configTree->Fill();
    configTree->Write();
    Tree->Write();
    TTreeFile->Close();
    delete TTreeFile;
    
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
  if(diMuon){ 
	  //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	  selecTableDiMu.TableCalculator(  false, true, true, true, true);
	  selecTableDiMuEl.TableCalculator(false, true, true, true, true);
	  selecTableTriMu.TableCalculator( false, true, true, true, true);
    //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
	  selecTableDiMu.Write(  "TopFCNC"+postfix+channelpostfix+comments+"_SelectionTable_DiMu.tex",    true,true,true,true,false,false,true);
	  selecTableDiMuEl.Write("TopFCNC"+postfix+channelpostfix+comments+"_SelectionTable_DiMuElec.tex",true,true,true,true,false,false,true);
	  selecTableTriMu.Write( "TopFCNC"+postfix+channelpostfix+comments+"_SelectionTable_TriMu.tex",   true,true,true,true,false,false,true);
  }
  else if(diElectron){
	  //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	  selecTableDiEl.TableCalculator(  false, true, true, true, true);
	  selecTableDiElMu.TableCalculator(false, true, true, true, true);
	  selecTableTriEl.TableCalculator( false, true, true, true, true);
    //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
	  selecTableDiEl.Write(  "TopFCNC"+postfix+channelpostfix+comments+"_SelectionTable_DiEl.tex",  true,true,true,true,false,false,true);
	  selecTableDiElMu.Write("TopFCNC"+postfix+channelpostfix+comments+"_SelectionTable_DiElMu.tex",true,true,true,true,false,false,true);
  	selecTableTriEl.Write( "TopFCNC"+postfix+channelpostfix+comments+"_SelectionTable_TriEl.tex", true,true,true,true,false,false,true);
  }
  fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
	  MultiSamplePlot *temp = it->second;
	  //temp->addText("CMS preliminary");
	  string name = it->first;
	  name += comments;
	  temp->Draw(false, name, true, true, true, true, true,1,false); // merge TT/QCD/W/Z/ST/
	  //Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
	  temp->Write(fout, name, true, pathPNG, "pdf");
  }
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
double MuEffSF_Id_Run2012(double eta, double pt){
  if( fabs(eta) >= 0 && fabs(eta) < 0.9 ){
    if( pt >= 10 && pt < 20)
      return 1.00747;
    else if( pt >= 20 && pt < 25)
      return 1.0102;
    else if( pt >= 25 && pt < 30)
      return 0.999714;
    else if( pt >= 30 && pt < 35)
      return 0.991798;
    else if( pt >= 35 && pt < 40)
      return 0.989783;
    else if( pt >= 40 && pt < 50)
      return 0.987796;
    else if( pt >= 50 && pt < 60)
      return 0.995079;
    else if( pt >= 60 && pt < 90)
      return 0.990791;
    else if( pt >= 90 && pt < 140)
      return 0.99111;
    else if( pt >= 140 && pt < 500)
      return 1.01688;
    else
      return 1;
  }
  if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ){
    if( pt >= 10 && pt < 20)
      return 1.02897;
    else if( pt >= 20 && pt < 25)
      return 1.00923;
    else if( pt >= 25 && pt < 30)
      return 0.993101;
    else if( pt >= 30 && pt < 35)
      return 0.988653;
    else if( pt >= 35 && pt < 40)
      return 0.989246;
    else if( pt >= 40 && pt < 50)
      return 0.986642;
    else if( pt >= 50 && pt < 60)
      return 0.991663;
    else if( pt >= 60 && pt < 90)
      return 0.994594;
    else if( pt >= 90 && pt < 140)
      return 0.957334;
    else if( pt >= 140 && pt < 500)
      return 0.876914;
    else
      return 1;
  }
  if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ){
    if( pt >= 10 && pt < 20)
      return 1.04513;
    else if( pt >= 20 && pt < 25)
      return 1.00195;
    else if( pt >= 25 && pt < 30)
      return 0.99624;
    else if( pt >= 30 && pt < 35)
      return 0.998091;
    else if( pt >= 35 && pt < 40)
      return 0.995041;
    else if( pt >= 40 && pt < 50)
      return 0.995362;
    else if( pt >= 50 && pt < 60)
      return 0.995984;
    else if( pt >= 60 && pt < 90)
      return 0.974529;
    else if( pt >= 90 && pt < 140)
      return 0.959219;
    else if( pt >= 140 && pt < 500)
      return 1.01609;
    else
      return 1;
  }
}
double MuEffSF_Iso04_Run2012(double eta, double pt){
  if( fabs(eta) >= 0 && fabs(eta) < 0.9 ){
    if( pt >= 10 && pt < 20)
      return 0.967563;
    else if( pt >= 20 && pt < 25)
      return 0.990809;
    else if( pt >= 25 && pt < 30)
      return 1.0029;
    else if( pt >= 30 && pt < 35)
      return 1.00061;
    else if( pt >= 35 && pt < 40)
      return 0.999133;
    else if( pt >= 40 && pt < 50)
      return 0.99914;
    else if( pt >= 50 && pt < 60)
      return 0.999507;
    else if( pt >= 60 && pt < 90)
      return 1.00458;
    else if( pt >= 90 && pt < 140)
      return 0.999191;
    else if( pt >= 140 && pt < 500)
      return 1.00421;
    else
      return 1;
  }
  if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ){
    if( pt >= 10 && pt < 20)
      return 0.973839;
    else if( pt >= 20 && pt < 25)
      return 0.990743;
    else if( pt >= 25 && pt < 30)
      return 1.00316;
    else if( pt >= 30 && pt < 35)
      return 1.00558;
    else if( pt >= 35 && pt < 40)
      return 1.00145;
    else if( pt >= 40 && pt < 50)
      return 0.999834;
    else if( pt >= 50 && pt < 60)
      return 1.00084;
    else if( pt >= 60 && pt < 90)
      return 1.00294;
    else if( pt >= 90 && pt < 140)
      return 0.995513;
    else if( pt >= 140 && pt < 500)
      return 0.985999;
    else
      return 1;
  }
  if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ){
    if( pt >= 10 && pt < 20)
      return 0.983218;
    else if( pt >= 20 && pt < 25)
      return 0.996806;
    else if( pt >= 25 && pt < 30)
      return 1.00645;
    else if( pt >= 30 && pt < 35)
      return 1.00497;
    else if( pt >= 35 && pt < 40)
      return 1.0019;
    else if( pt >= 40 && pt < 50)
      return 1.00078;
    else if( pt >= 50 && pt < 60)
      return 1.00091;
    else if( pt >= 60 && pt < 90)
      return 1.00086;
    else if( pt >= 90 && pt < 140)
      return 0.995858;
    else if( pt >= 140 && pt < 500)
      return 1.0017;
    else
      return 1;
  }
}
double MuEffSF_TrgMu8_Run2012(double eta, double pt){
  if( fabs(eta) >= 0 && fabs(eta) < 0.9 ){
    if( pt >= 10 && pt < 20)
      return 0.997657;
    else if( pt >= 20 && pt < 30)
      return 0.990708;
    else if( pt >= 30 && pt < 35)
      return 0.989016;
    else if( pt >= 35 && pt < 40)
      return 0.989852;
    else if( pt >= 40 && pt < 50)
      return 0.989148;
    else if( pt >= 50 && pt < 60)
      return 0.988999;
    else if( pt >= 60 && pt < 90)
      return 0.98839;
    else if( pt >= 90 && pt < 140)
      return 0.98725;
    else if( pt >= 140 && pt < 500)
      return 0.989903;
    else
      return 1;
  }
  if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ){
    if( pt >= 10 && pt < 20)
      return 1.00496;
    else if( pt >= 20 && pt < 30)
      return 0.990941;
    else if( pt >= 30 && pt < 35)
      return 0.986394;
    else if( pt >= 35 && pt < 40)
      return 0.986699;
    else if( pt >= 40 && pt < 50)
      return 0.98872;
    else if( pt >= 50 && pt < 60)
      return 0.985978;
    else if( pt >= 60 && pt < 90)
      return 0.980763;
    else if( pt >= 90 && pt < 140)
      return 0.987128;
    else if( pt >= 140 && pt < 500)
      return 0.954924;
    else
      return 1;
  }
  if( fabs(eta) >= 1.2 && fabs(eta) < 2.4 ){
    if( pt >= 10 && pt < 20)
      return 0.994495;
    else if( pt >= 20 && pt < 30)
      return 0.992625;
    else if( pt >= 30 && pt < 35)
      return 0.989671;
    else if( pt >= 35 && pt < 40)
      return 0.988672;
    else if( pt >= 40 && pt < 50)
      return 0.989098;
    else if( pt >= 50 && pt < 60)
      return 0.984751;
    else if( pt >= 60 && pt < 90)
      return 0.990953;
    else if( pt >= 90 && pt < 140)
      return 0.975961;
    else if( pt >= 140 && pt < 500)
      return 0.99754;
    else
      return 1;
  }
}
double MuEffSF_TrgMu17_Run2012(double eta, double pt){
  if( fabs(eta) >= 0 && fabs(eta) < 0.9 ){
    if( pt >= 10 && pt < 20)
      return 0.987431;
    else if( pt >= 20 && pt < 30)
      return 0.990293;
    else if( pt >= 30 && pt < 35)
      return 0.988535;
    else if( pt >= 35 && pt < 40)
      return 0.98961;
    else if( pt >= 40 && pt < 50)
      return 0.988796;
    else if( pt >= 50 && pt < 60)
      return 0.988964;
    else if( pt >= 60 && pt < 90)
      return 0.988338;
    else if( pt >= 90 && pt < 140)
      return 0.986894;
    else if( pt >= 140 && pt < 500)
      return 0.98723;
    else
      return 1;
  }
  if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ){
    if( pt >= 10 && pt < 20)
      return 1.0148;
    else if( pt >= 20 && pt < 30)
      return 0.987789;
    else if( pt >= 30 && pt < 35)
      return 0.984216;
    else if( pt >= 35 && pt < 40)
      return 0.984829;
    else if( pt >= 40 && pt < 50)
      return 0.987931;
    else if( pt >= 50 && pt < 60)
      return 0.985025;
    else if( pt >= 60 && pt < 90)
      return 0.977625;
    else if( pt >= 90 && pt < 140)
      return 0.98885;
    else if( pt >= 140 && pt < 500)
      return 0.954924;
    else
      return 1;
  }
  if( fabs(eta) >= 1.2 && fabs(eta) < 2.4 ){
    if( pt >= 10 && pt < 20)
      return 1.00149;
    else if( pt >= 20 && pt < 30)
      return 0.998639;
    else if( pt >= 30 && pt < 35)
      return 0.992417;
    else if( pt >= 35 && pt < 40)
      return 0.990839;
    else if( pt >= 40 && pt < 50)
      return 0.989054;
    else if( pt >= 50 && pt < 60)
      return 0.98503;
    else if( pt >= 60 && pt < 90)
      return 0.990908;
    else if( pt >= 90 && pt < 140)
      return 0.973882;
    else if( pt >= 140 && pt < 500)
      return 1.00649;
    else
      return 1;
  }
}
