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
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
//#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
//#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
//#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysis/TopFCNC/interface/TopFCNC_Evt.h"

//#include "../../macros/Style.C"

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

double MuEffSF_Id_Run2012(string period, double eta, double pt);
double MuEffSF_Iso04_Run2012(string period, double eta, double pt);

double MuEffSF_Id_Run2012ABCD(double eta, double pt);
double MuEffSF_Iso04_Run2012ABCD(double eta, double pt);

double ElEffSF_Id_Run2012(string period, double eta, double pt);
double ElEffSF_Iso04_Run2012(string period, double eta, double pt);

double ElEffSF_Id_Run2012ABCD(double eta, double pt);
double ElEffSF_Iso04_Run2012ABCD(double eta, double pt);

double MuonRelIso(TRootMuon* muon);
double ElecRelIso(TRootElectron* elec, TRootEvent* event);

bool ZCandInvMass(vector<TRootMuon*>& muons, vector<TRootElectron*>& electrons, int& idx_Z_1, int& idx_Z_2, double& invMass, double Zmass, double Zwindowsize);

int main (int argc, char *argv[])
{
  clock_t start = clock();

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
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
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
  int   btagAlgo     = 6; //CSV
  string btagName[8] = {"TCHE","TCHP","JP","JBP","SVHE","SVHP","CSV","CSV_MVA"};
  float btagCut      = 0.679;
	const double Zmass        = 91.2;
	const double Zwindowsize  = 30.;
  const double MuPtCut      = 20.;
  const double MuEtaCut     = 2.4;
  const double MuRelIsoCut  = 0.20;
  const double ElPtCut      = 20.;
  const double ElEtaCut     = 2.5;
  const double ElRelIsoCut  = 0.15;
  const double JetPtCut     = 20;
  const double JetEtaCut    = 2.4;
	bool applyAsymmJetPtCut   = true;
	const double JetPtCuts[4] = {50.,40.,30.,20.};
  bool applyLeptonSF        = true;
  
  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FCNC search ! " << endl;
  cout << "*************************************************************" << endl;
  cout << " - Selection cuts " << endl;
  cout << " -- Lepton ID&ISO SF: " << (applyLeptonSF ? "TRUE" : "FALSE" ) << endl;
  cout << " -- Muons: " << endl;
  cout << " --- pt: "  << MuPtCut << endl;
  cout << " --- eta: " << MuEtaCut << endl;
  cout << " --- iso: " << MuRelIsoCut << endl;
  cout << " -- Electrons: " << endl;
  cout << " --- pt: "  << ElPtCut << endl;
  cout << " --- eta: " << ElEtaCut << endl;
  cout << " --- iso: " << ElRelIsoCut << endl;
  cout << " -- Jets: " << endl;
  cout << " --- lower pt cut : " << JetPtCut << endl;
  cout << " --- eta: " << JetEtaCut << endl;
  cout << " --- asymmetric pt cuts ? " << (applyAsymmJetPtCut ? "TRUE" : "FALSE" ) << endl;
  cout << " --- pt: " << JetPtCuts[0] << " / " << JetPtCuts[1] << " / " << JetPtCuts[2] << " / " << JetPtCuts[3] << endl;
  cout << " -- B-tagging: " << endl;
  cout << " --- algorithm: " << btagName[btagAlgo] << endl;
  cout << " --- threshold: " << btagCut << endl;
  cout << " -- Z boson mass = " << Zmass<< endl;
  cout << " -- Z boson mass window = " << Zwindowsize<< endl;
    
  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FCNC search ! " << endl;
  cout << "*************************************************************" << endl;
  
  //SetStyle if needed
  //setTDRStyle();
  //setGregStyle();
  //setMyStyle();
  
  string postfix = "_BckgdEstimation"; // to relabel the names of the output file
  
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
  
  string channelpostfix = "_MuEGTrigger";
  string comments = "_Run2012ABCD_SkimMuEG_3Jets";
  string xmlFileName = "../config/TopFCNCconfig_EventSelection_MuEG_SkimMuEG3Jets.xml";

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
  
  MSPlot["MuonDzero"]                         = new MultiSamplePlot(datasets, "MuonDzero", 500, -0.02, 0.02, "d_{0} [cm]");
  MSPlot["ElecDzero"]                         = new MultiSamplePlot(datasets, "ElecDzero", 500, -0.02, 0.02, "d_{0} [cm]");
  
  MSPlot["1stLeadingMuonPt"]                  = new MultiSamplePlot(datasets, "1stLeadingMuonPt", 300, 0, 150, "p_{T} [GeV/c]");
  MSPlot["1stLeadingElecPt"]                  = new MultiSamplePlot(datasets, "1stLeadingElecPt", 300, 0, 150, "p_{T} [GeV/c]");
  
  MSPlot["1stLeadingMuonRelIso"]              = new MultiSamplePlot(datasets, "1stLeadingMuonRelIso", 100, 0, 0.5, "RelIso");
  MSPlot["1stLeadingElecRelIso"]              = new MultiSamplePlot(datasets, "1stLeadingElecRelIso", 100, 0, 0.5, "RelIso");
  
  MSPlot["NbOfIsolatedMuons"]                 = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfIsolatedElectrons"]             = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
  
  MSPlot["DiLeptonInvMass_me_ch"]            = new MultiSamplePlot(datasets, "DiLeptonInvMass_me_ch", 400, 50, 130, "m_{ll}");
  MSPlot["DiLeptonSystPt_me_ch"]             = new MultiSamplePlot(datasets, "DiLeptonSystPt_me_ch", 400, 0, 400, "p_{T}^{ll} [GeV/c]");
  MSPlot["DiLeptonDR_me_ch"]                 = new MultiSamplePlot(datasets, "DiLeptonDR_me_ch", 400, 0, 10, "#Delta R(l^{+}l^{-})");

  MSPlot["DiLeptonInvMass_AtLeastFourJets_me_ch"]= new MultiSamplePlot(datasets, "DiLeptonInvMass_AtLeastFourJets_me_ch", 400, 50, 130, "m_{ll}");
  MSPlot["DiLeptonDR_AtLeastFourJets_me_ch"]     = new MultiSamplePlot(datasets, "DiLeptonDR_AtLeastFourJets_me_ch", 100, 0, 5, "#Delta R(l^{+}l^{-})");
  MSPlot["DiLeptonDPhi_AtLeastFourJets_me_ch"]   = new MultiSamplePlot(datasets, "DiLeptonDPhi_AtLeastFourJets_me_ch", 70, 0, 3.5, "#Delta #Phi(l^{+}l^{-})");
  MSPlot["DiLeptonSystPt_AtLeastFourJets_me_ch"] = new MultiSamplePlot(datasets, "DiLeptonSystPt_AtLeastFourJets_me_ch", 100, 0, 400, "p_{T}^{ll} [GeV/c]");
  
  MSPlot["DiLeptonSystDPhi_LeadingJet_AtLeastFourJets_me_ch"]   = new MultiSamplePlot(datasets, "DiLeptonSystDPhi_LeadingJet_AtLeastFourJets_me_ch", 70, 0, 3.5, "#Delta #Phi(l^{+}+l^{-},#Sigma jet))");
  MSPlot["DiLeptonSystDPhi_JetSyst_AtLeastFourJets_me_ch"]      = new MultiSamplePlot(datasets, "DiLeptonSystDPhi_JetSyst_AtLeastFourJets_me_ch", 70, 0, 3.5, "#Delta #Phi(l^{+}+l^{-},#Sigma jet))");
  MSPlot["DiLeptonSystPt_Over_JetSystPt_AtLeastFourJets_me_ch"] = new MultiSamplePlot(datasets, "DiLeptonSystPt_Over_JetSystPt_AtLeastFourJets_me_ch", 500, 0, 10, "p_{T}^{ll}/p_{T}^{#Sigma jets}");
  
  MSPlot["NbOfExtraIsolatedMuons"]            = new MultiSamplePlot(datasets, "NbOfExtraIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfExtraIsolatedElectrons"]        = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
  
  MSPlot["NbOfSelectedJets_me_ch"]            = new MultiSamplePlot(datasets, "NbOfSelectedJets_me_ch",  15, 0, 15, "Nb. of jets");
  
  MSPlot["FirstLeadingJetPt_me_ch"]           = new MultiSamplePlot(datasets, "FirstLeadingJetPt_me_ch", 100, 0, 500, "Jet p_{T} [GeV/c]");
  MSPlot["SecondLeadingJetPt_me_ch"]          = new MultiSamplePlot(datasets, "SecondLeadingJetPt_me_ch", 80, 0, 400, "Jet p_{T} [GeV/c]");
  MSPlot["ThirdLeadingJetPt_me_ch"]           = new MultiSamplePlot(datasets, "ThirdLeadingJetPt_me_ch",  60, 0, 300, "Jet p_{T} [GeV/c]");
  MSPlot["FourthLeadingJetPt_me_ch"]          = new MultiSamplePlot(datasets, "FourthLeadingJetPt_me_ch", 50, 0, 250, "Jet p_{T} [GeV/c]");
  
  MSPlot["NbOfVertices_AtLeastFourJets_me_ch"]= new MultiSamplePlot(datasets, "NbOfVertices_AtLeastFourJets_me_ch", 30, 0, 30, "Nb. of vertices");
  
  MSPlot["HighestBdisc_me_ch_CVS"]            = new MultiSamplePlot(datasets, "HighestBdisc_me_ch_CVS", 50, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_me_ch_JP"]             = new MultiSamplePlot(datasets, "HighestBdisc_me_ch_JP" , 50, 0, 1, "JP b-disc.");

  MSPlot["MET_AtLeastFourJets_me_ch"]         = new MultiSamplePlot(datasets, "MET_AtLeastFourJets_me_ch",  50, 0, 200, "\\slashE_{T} [GeV]");
  
  MSPlot["NbOfSelectedBJets_AtLeastFourJets_me_ch_CSV"]= new MultiSamplePlot(datasets, "NbOfSelectedBJets_AtLeastFourJets_me_ch_CSV",  10, 0, 10, "Nb. of b-tagged jets (CSVM)");
  
  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",50,0,4);
  for (unsigned int d = 0; d < datasets.size(); d++){
    histo1D[("scaleFactors_"+datasets[d]->Name()).c_str()] = new TH1F(("ScaleFactors_"+datasets[d]->Name()).c_str(),";Total SF;#events",500,0,2);
  	histo2D[("d0_vs_phi_1stleadingmuon_"+datasets[d]->Name()).c_str()] = new TH2F(("d0_vs_phi_1stleadingmuon_"+datasets[d]->Name()).c_str(),"d_{0}:#phi",500,-0.02,0.02,500,0,4);
  	histo2D[("d0_vs_phi_1stleadingelec_"+datasets[d]->Name()).c_str()] = new TH2F(("d0_vs_phi_1stleadingelec_"+datasets[d]->Name()).c_str(),"d_{0}:#phi",500,-0.02,0.02,500,0,4);
  }
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
  sprintf(zmasscutname,"$|m_{ll}-m_Z|<%.0f$ GeV",Zwindowsize);
  char btagcutname[100];
  sprintf(btagcutname,"$b\\texttt{-}disc \\geq %1.3f$ (CSV)",btagCut);
  
  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µe  ////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  vector<string> CutsSelecTableMuEl;
  CutsSelecTableMuEl.push_back(string("initial"));
  CutsSelecTableMuEl.push_back(string("PU reweighting"));
  CutsSelecTableMuEl.push_back(string("Trigger"));
  CutsSelecTableMuEl.push_back(string("Good PV"));
  CutsSelecTableMuEl.push_back(string("$\\geq$ 2 isolated leptons"));
  CutsSelecTableMuEl.push_back(string(zmasscutname));//"$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableMuEl.push_back(string("Veto on 3rd iso. lept."));
  CutsSelecTableMuEl.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableMuEl.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableMuEl.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableMuEl.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableMuEl.push_back(string(btagcutname));//"$b\\texttt{-}disc \\geq 0.7$ (CSV)"));
  //  CutsSelecTableMuEl.push_back(string("$MET \\leq 30$ GeV"));
  
  SelectionTable selecTableMuEl(CutsSelecTableMuEl, datasets);
  selecTableMuEl.SetLuminosity(Luminosity);
  selecTableMuEl.SetPrecision(1);
    
  cout << " - SelectionTable instantiated ..." << endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////// PileUp Reweighting ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  string PileUpFile = "";
  string PileUpSyst = ( doPUShift==0 ? "" : ( doPUShift==1 ? "_Systematic_Down_5perc" : "_Systematic_Up_5perc" ) );
  if( comments.find("2012ABCD") != string::npos ) PileUpFile = "../pileup/MuEG_Run2012ABCD_22Jan2013_TopTreeID_3496to3499_PileupHistogram"+PileUpSyst+".root";
  /*
   else if( comments.find("2012AB")   != string::npos ) PileUpFile = (diMuon ? "../pileup/DoubleMu_Run2012AB_22Jan2013_TopTreeID_3494_3495_PileupHistogram.root" : "../pileup/XXX.root");
   else if( comments.find("2012C")    != string::npos ) PileUpFile = (diMuon ? "../pileup/DoubleMu_Run2012C_Promptv2_TopTreeID_2226_PileupHistogram.root" : "../pileup/XXX.root");
   else if( comments.find("2012D")    != string::npos ) PileUpFile = (diMuon ? "../pileup/DoubleMu_Run2012D_Promptv1_TopTreeID_2218_PileupHistogram.root" : "../pileup/XXX.root");
   */
  else{
    cerr << "Cannot determine which pile up root should be used. Make sure the string comment is filled correctly." << endl;
    exit(1);
  }
  LumiReWeighting LumiWeights = LumiReWeighting("../pileup/pileup_MC_S10.root", PileUpFile, "pileup", "pileup");
  
  cout << " - Initialized LumiReWeighting with file : " << PileUpFile << endl;
  
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
    
    bool isData = false;
    string dataSetName = datasets[d]->Name();
    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
      isData = true;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////// Initialize JEC factors /////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
    vector<JetCorrectorParameters> vCorrParam;
    
    // Create the JetCorrectorParameter objects, the order does not matter.
    // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    
    string CalibPath = "../../../TopTreeAnalysisBase/Calibrations/JECFiles/";
    
    string L1Corr = (isData?"FT_53_V21_AN4_Summer13_Data_L1FastJet_AK5PFchs.txt":"START53_V23_Summer13_L1FastJet_AK5PFchs.txt");
    string L2Corr = (isData?"FT_53_V21_AN4_Summer13_Data_L2Relative_AK5PFchs.txt":"START53_V23_Summer13_L2Relative_AK5PFchs.txt");
    string L3Corr = (isData?"FT_53_V21_AN4_Summer13_Data_L3Absolute_AK5PFchs.txt":"START53_V23_Summer13_L3Absolute_AK5PFchs.txt");
    string L2L3Rs = (isData?"FT_53_V21_AN4_Summer13_Data_L2L3Residual_AK5PFchs.txt":"START53_V23_Summer13_L2L3Residual_AK5PFchs.txt");
    string JECUnc = (isData?"FT_53_V21_AN4_Summer13_Data_Uncertainty_AK5PFchs.txt":"START53_V23_Summer13_Uncertainty_AK5PFchs.txt");
    
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(CalibPath+L1Corr);
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(CalibPath+L2Corr);
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(CalibPath+L3Corr);
    
    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);
    if(isData) // Data!
    {
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters(CalibPath+L2L3Rs);
      vCorrParam.push_back(*ResJetCorPar);
    }
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(CalibPath+JECUnc);
    
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
    
    // Set cache size.
    Int_t cachesize=10000000; // 10 MB = 10000000
    
    datasets[d]->eventTree()->SetCacheSize(cachesize);
    datasets[d]->eventTree()->SetCacheEntryRange(start,end);
    datasets[d]->eventTree()->AddBranchToCache("*",kTRUE);
    //  datasets[d]->eventTree()->SetBranchStatus("*",0);
    datasets[d]->eventTree()->StopCacheLearningPhase();
    
    datasets[d]->runTree()->SetCacheSize(cachesize);
    datasets[d]->runTree()->SetCacheEntryRange(start,end);
    datasets[d]->runTree()->AddBranchToCache("*",kTRUE);
    //  datasets[d]->runTree()->SetBranchStatus("*",0);
    datasets[d]->runTree()->StopCacheLearningPhase();
    
    //ps[d] = new TTreePerfStats(("ioperf_"+dataSetName).c_str(),datasets[d]->eventTree());
    
    for (unsigned int ievt = start; ievt < end; ievt++)
    {
      
      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";
      
      //load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      
      vector<TRootGenJet*> genjets;
      if(!isData)
		    genjets = treeLoader.LoadGenJet(ievt,false);

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
        if(isData)
        {
          /*------------------------------------------------------------------
           Dataset : MuEG/Run2012A-22Jan2013-v1
           ------------------------------------------------------------------*/
          if(currentRun >= 190645 && currentRun <= 190738){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4"), currentRun, iFile);
					  itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4"), currentRun, iFile);
            // int. lumi = 94.005/pb
					}
          else if (currentRun >= 190782 && currentRun <= 191367){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5"), currentRun, iFile);
            itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5"), currentRun, iFile);
            // int. lumi = 385.304/pb
          }
          else if (currentRun >= 191695 && currentRun <= 193621){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6"), currentRun, iFile);
            itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6"), currentRun, iFile);
            // int. lumi = 384.511/pb
          }
				  /*------------------------------------------------------------------
           Sub-Total integrated luminosity = 863.8(/pb)
           Total integrated luminosity = 863.8(/pb)
           ------------------------------------------------------------------*/
          
				  /*------------------------------------------------------------------
           Dataset : MuEG/Run2012B-22Jan2013-v1
           ------------------------------------------------------------------*/
          else if (currentRun >= 193806 && currentRun <= 196531){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"), currentRun, iFile);
            itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"), currentRun, iFile);
          }
          /*------------------------------------------------------------------
           Sub-Total integrated luminosity = 4375(/pb)
           Total integrated luminosity = 5238.8(/pb)
           -----------------------------------------------------------------*/

				  /*------------------------------------------------------------------
           Dataset : MuEG/Run2012C-22Jan2013-v1
           -------------------------------------------------------------------
           Will search for runs with trigger HLT_Mu*_Ele*_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v* available
           Trigger HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8 available for runs 198049-199608
           Trigger HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9 available for runs 199698-203742
           Trigger HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8 available for runs 198049-199608
           Trigger HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9 available for runs 199698-203742

           ------------------------------------------------------------------*/
          else if (currentRun >= 198049 && currentRun <= 199608){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8"), currentRun, iFile);
            itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8"), currentRun, iFile);
            // int. lumi = 1783/pb
          }
          else if (currentRun >= 199698 && currentRun <= 203742){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9"), currentRun, iFile);
            itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9"), currentRun, iFile);
            // int. lumi = 5363/pb
          }
          /*------------------------------------------------------------------
           Sub-Total integrated luminosity = 7146(/pb)
           Total integrated luminosity = 12384.8(/pb)
           -----------------------------------------------------------------*/

				  /*------------------------------------------------------------------
           Dataset : MuEG/Run2012D-22Jan2013-v1
           --------------------------------------------------------------------
           Will search for runs with trigger HLT_Mu*_Ele*_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v* available
           Trigger HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9 available for runs 203777-208686
           Trigger HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9 available for runs 203777-208686
           ------------------------------------------------------------------*/
          else if (currentRun >= 203777 && currentRun <= 208686){
            itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9"), currentRun, iFile);
            itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9"), currentRun, iFile);
          }
          /*------------------------------------------------------------------
           Sub-Total integrated luminosity = 7292(/pb)
           Total integrated luminosity = 19676.8(/pb)
           -----------------------------------------------------------------*/

          if(itrigger1 == 9999 && itrigger2 == 9999)
          {
            cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
            exit(1);
          }
          //trigged = treeLoader.EventTrigged (itrigger);
        }
        else{
          itrigger1 = treeLoader.iTrigger (string ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"), currentRun, iFile);
          itrigger2 = treeLoader.iTrigger (string ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"), currentRun, iFile);
          if(itrigger1 == 9999 || itrigger2 == 9999)
          {
            cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
            exit(1);
          }
          cout<<"Trigger1 bit nr : "<<itrigger1<<endl;
          cout<<"Trigger2 bit nr : "<<itrigger2<<endl;
        }
	    } //end previousRun != currentRun
      
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    //////////////////////////////////////////// Jet energy scale corrections /////////////////////////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
	    // Apply Jet Corrections on-the-fly
	    //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
      jetTools->correctJets(init_jets,event->kt6PFJets_rho(),isData); //last boolean: isData (needed for L2L3Residual...)
      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");
      
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    //////////////////////////////////////////// Type I MET corrections: (Only for |eta| <=4.7 ) //////////////////////////////////////////////
	    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");
      jetTools->correctMETTypeOne(init_jets,mets[0],isData);
      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After MET type I correction:");
      
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////// Jet energy smearing and systematic uncertainty ///////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      if(!isData)
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
      
      if(isData)
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
        // Apply pile-up reweighting
        double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
        scaleFactor *= lumiWeight;
      }
      
      histo1D["lumiWeights"]->Fill(scaleFactor);
      MSPlot["RhoCorrection"]->Fill(event->kt6PFJets_rho(), datasets[d], true, Luminosity*scaleFactor);
			
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////// Event selection ////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      selecTableMuEl.Fill(d,0, 1.);
      selecTableMuEl.Fill(d,1,scaleFactor);
      
      //__Trigger selection_________________________________________________________________
      //      if (itrigger1 == -1 && itrigger2 == -1)
      //        trigged = true;
      //      else
      trigged = (treeLoader.EventTrigged (itrigger1) || treeLoader.EventTrigged (itrigger2));
      if(!trigged) continue;
      
      //__Trigger SF_______________________________________________________________________
      if(!isData){
        if(comments.find("2012ABCD") != string::npos){
          scaleFactor *= 1.; // To be modified
        }
        else if(comments.find("2012AB") != string::npos){
          scaleFactor *= 1.; // To be modified
        }
        else if(comments.find("2012C") != string::npos){
          scaleFactor *= 1.; // To be modified
        }
        else{
          cerr << "__Lepton trigger SF: wrong period." << endl;
          exit(1);
        }
      }
      
      selecTableMuEl.Fill(d,2,scaleFactor);
      
      //__Selection instances_________________________________________________________________
      Selection selection(init_jets, init_muons, init_electrons, mets, event->kt6PFJets_rho()); //mets can also be corrected...
      //Selection selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...
      
      // Specify the correction to be applied to the electron isolation
      selection.setElectronIsoCorrType(1); // 0: no corr, 1: EA corr, 2: dB corr
      
      //__Object election cuts_________________________________________________________________
      selection.setJetCuts(JetPtCut,JetEtaCut,0.01,1.,0.98,0.3,0.1);
      
      selection.setDiElectronCuts(ElPtCut,ElEtaCut,ElRelIsoCut,0.04,0.,1,0.3,1); //Et,Eta,RelIso,d0,MVAId,DistVzPVz,DRJets,MaxMissHits
      //selection.setLooseElectronCuts(15,2.5,0.2);
      
      selection.setDiMuonCuts(MuPtCut,MuEtaCut,MuRelIsoCut,999.); //Et,Eta,RelIso,d0
      //selection.setLooseMuonCuts(15,2.4,0.2);
      
      //__Object selection_________________________________________________________________
      //vector<TRootElectron*> selectedElectrons_NoIso = selection.GetSelectedDiElectrons(ElPtCut,ElEtaCut,999.); //Et,Eta,RelIso
      vector<TRootElectron*> selectedElectrons       = selection.GetSelectedDiElectrons();
      vector<TRootElectron*> selectedExtraElectrons;
      
      //vector<TRootMuon*>     selectedMuons_NoIso = selection.GetSelectedDiMuons(MuPtCut,MuEtaCut,999.); //Et,Eta,RelIso
      vector<TRootMuon*>     selectedMuons       = selection.GetSelectedDiMuons();
      vector<TRootMuon*>     selectedExtraMuons;
      
      vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId
      
      //vector<TRootMuon*>     looseMuons     = selection.GetSelectedLooseMuons();
      //vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseElectrons(true); // VBTF Id
      
      //__Primary vertex selection_________________________________________________________________
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
      if(!isGoodPV) continue;
      selecTableMuEl.Fill(d,3,scaleFactor);
      
      //__Lepton ID&&ISO SF_____________________________________________________________________________
      if(!isData && applyLeptonSF){
        for(unsigned int i=0;i<selectedMuons.size();i++){
          scaleFactor *= MuEffSF_Id_Run2012(comments, selectedMuons[i]->Eta(), selectedMuons[i]->Pt());// Id SF
          scaleFactor *= MuEffSF_Iso04_Run2012(comments, selectedMuons[i]->Eta(), selectedMuons[i]->Pt());// Iso SF
        }
        for(unsigned int i=0;i<selectedElectrons.size();i++){
          scaleFactor *= ElEffSF_Id_Run2012(comments, selectedElectrons[i]->Eta(), selectedElectrons[i]->Pt());// Id SF
          scaleFactor *= ElEffSF_Iso04_Run2012(comments, selectedElectrons[i]->Eta(), selectedElectrons[i]->Pt());// Iso SF
        }
      }

      MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
      histo1D[("scaleFactors_"+datasets[d]->Name()).c_str()]->Fill(scaleFactor);

      //__Dilepton selection_________________________________________________________________
      if( selectedMuons.size()<1 || selectedElectrons.size()<1 ) continue;

      histo2D[("d0_vs_phi_1stleadingmuon_"+datasets[d]->Name()).c_str()]->Fill(selectedMuons[0]->d0(),selectedMuons[0]->Phi());
      histo2D[("d0_vs_phi_1stleadingelec_"+datasets[d]->Name()).c_str()]->Fill(selectedElectrons[0]->d0(),selectedElectrons[0]->Phi());
      for(unsigned int i=0;i<selectedMuons.size();i++)
        MSPlot["MuonDzero"]->Fill(selectedMuons[i]->d0(), datasets[d], true, Luminosity*scaleFactor);
      for(unsigned int i=0;i<selectedElectrons.size();i++)
        MSPlot["ElecDzero"]->Fill(selectedElectrons[i]->d0(), datasets[d], true, Luminosity*scaleFactor);

      MSPlot["1stLeadingMuonPt"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["1stLeadingElecPt"]->Fill(selectedElectrons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);

      float reliso_mu = MuonRelIso(selectedMuons[0]);
      float reliso_el = ElecRelIso(selectedElectrons[0], event);
      MSPlot["1stLeadingMuonRelIso"]->Fill(reliso_mu, datasets[d], true, Luminosity*scaleFactor);
      MSPlot["1stLeadingElecRelIso"]->Fill(reliso_el, datasets[d], true, Luminosity*scaleFactor);
      
      selecTableMuEl.Fill(d,4,scaleFactor);
      
      //__Z candidate selection_________________________________________________________________
      bool foundZ = false;
      int idx_Z_mu = -1, idx_Z_el = -1;
      double invMass = 0;
      // Calculate the invariant mass for each isolated lepton pairs
      // - return true if the mass is the Z boson mass window
      // - return the indices of the lepton candidates
      foundZ = ZCandInvMass(selectedMuons, selectedElectrons, idx_Z_mu, idx_Z_el, invMass, Zmass, Zwindowsize);
      MSPlot["DiLeptonInvMass_me_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
      // Select events with at least one pair of opposite charge leptons with |mll-mz|<windowsize
      if(!foundZ)
        continue;
      selecTableMuEl.Fill(d,5,scaleFactor);

      float DiLeptonSystPt = (*selectedMuons[idx_Z_mu] + *selectedElectrons[idx_Z_el]).Pt();
      float DiLeptonDR     =   selectedMuons[idx_Z_mu]->DeltaR(*selectedElectrons[idx_Z_el]);
      float DiLeptonDPhi   =   selectedMuons[idx_Z_mu]->DeltaPhi(*selectedElectrons[idx_Z_el]);
      TLorentzVector JetSyst; // initialized by (0.,0.,0.,0.)
      for(unsigned int i=0; i<selectedJets.size();i++){
        JetSyst += *selectedJets[i];
      }
      float JetSystPt = JetSyst.Pt();
      float DiLeptonSystDPhi_JetSyst = (*selectedMuons[idx_Z_mu] + *selectedElectrons[idx_Z_el]).DeltaPhi(JetSyst);
      float DiLeptonSystDPhi_LeadingJet = 0;
      if(selectedJets.size()>0) DiLeptonSystDPhi_LeadingJet = (*selectedMuons[idx_Z_mu] + *selectedElectrons[idx_Z_el]).DeltaPhi(*selectedJets[0]);

      MSPlot["DiLeptonSystPt_me_ch"]->Fill(DiLeptonSystPt, datasets[d], true, Luminosity*scaleFactor);
      MSPlot["DiLeptonDR_me_ch"]    ->Fill(DiLeptonDR,     datasets[d], true, Luminosity*scaleFactor);
      MSPlot["DiLeptonDPhi_me_ch"]  ->Fill(DiLeptonDPhi,   datasets[d], true, Luminosity*scaleFactor);
      MSPlot["DiLeptonSystDPhi_JetSyst_me_ch"]     ->Fill(DiLeptonSystDPhi_JetSyst, datasets[d], true, Luminosity*scaleFactor);
      MSPlot["DiLeptonSystPt_Over_JetSystPt_me_ch"]->Fill(DiLeptonSystPt/JetSystPt, datasets[d], true, Luminosity*scaleFactor);
      
      selectedExtraMuons = selectedMuons;
      selectedExtraElectrons = selectedElectrons;
      // Erase Z boson lepton candidates
      selectedExtraMuons.erase(selectedExtraMuons.begin()+idx_Z_mu);
      selectedExtraElectrons.erase(selectedExtraElectrons.begin()+idx_Z_el);
      
      MSPlot["NbOfExtraIsolatedMuons"]->Fill(selectedExtraMuons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfExtraIsolatedElectrons"]->Fill(selectedExtraElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
      
      MyTopFCNC_EvtCand = 0;
      double highestbtagdisc = 0;
      int NbOfBtagged = 0;

      // Select events based on the presence of *exactly one* extra isolated lepton
      if(selectedExtraMuons.size()==0 && selectedElectrons.size()==0){
        selecTableMuEl.Fill(d,6,scaleFactor);
        MSPlot["NbOfSelectedJets_me_ch"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
        if(selectedJets.size()>0){ //at least 1 jet
          MSPlot["FirstLeadingJetPt_me_ch"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
          if(applyAsymmJetPtCut && selectedJets[0]->Pt()<JetPtCuts[0]) continue;
          selecTableMuEl.Fill(d,7,scaleFactor);
          if(selectedJets.size()>1){ //at least 2 jets
            MSPlot["SecondLeadingJetPt_me_ch"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            if(applyAsymmJetPtCut && selectedJets[1]->Pt()<JetPtCuts[1]) continue;
            selecTableMuEl.Fill(d,8,scaleFactor);
            if(selectedJets.size()>2){ //at least 3 jets
              MSPlot["ThirdLeadingJetPt_me_ch"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
              if(applyAsymmJetPtCut && selectedJets[2]->Pt()<JetPtCuts[2]) continue;
              selecTableMuEl.Fill(d,9,scaleFactor);
              if(selectedJets.size()>3){ //at least 4 jets
                MSPlot["FourthLeadingJetPt_me_ch"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                if(applyAsymmJetPtCut && selectedJets[3]->Pt()<JetPtCuts[3]) continue;
                selecTableMuEl.Fill(d,10,scaleFactor);

                MSPlot["NbOfVertices_AtLeastFourJets_me_ch"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["DiLeptonInvMass_AtLeastFourJets_me_ch"]->Fill(invMass, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["DiLeptonDR_AtLeastFourJets_me_ch"]->Fill(DiLeptonDR, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["DiLeptonDPhi_AtLeastFourJets_me_ch"]->Fill(DiLeptonDPhi, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["DiLeptonSystPt_AtLeastFourJets_me_ch"]->Fill(DiLeptonSystPt, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["DiLeptonSystDPhi_LeadingJet_AtLeastFourJets_me_ch"]->Fill(DiLeptonSystDPhi_LeadingJet, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["DiLeptonSystDPhi_JetSyst_AtLeastFourJets_me_ch"]->Fill(DiLeptonSystDPhi_JetSyst, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["DiLeptonSystPt_Over_JetSystPt_AtLeastFourJets_me_ch"]->Fill(DiLeptonSystPt/JetSystPt, datasets[d], true, Luminosity*scaleFactor);

                sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
                MSPlot["HighestBdisc_me_ch_CVS"]->Fill(selectedJets[0]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
                highestbtagdisc = selectedJets[0]->btag_combinedSecondaryVertexBJetTags();
                sort(selectedJets.begin(),selectedJets.end(),HighestJPBtag());
                MSPlot["HighestBdisc_me_ch_JP"]->Fill(selectedJets[0]->btag_jetProbabilityBJetTags(),datasets[d], true, Luminosity*scaleFactor);
 
                MSPlot["MET_AtLeastFourJets_me_ch"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
                
                NbOfBtagged = selection.GetSelectedBJets(selectedJets, btagAlgo, btagCut).size();
                MSPlot["NbOfSelectedBJets_AtLeastFourJets_me_ch_CSV"]->Fill(NbOfBtagged, datasets[d], true, Luminosity*scaleFactor);
                if(NbOfBtagged>=1)
                  selecTableMuEl.Fill(d,11,scaleFactor);
                
                // Create TopFCNC_Evt object
                MyTopFCNC_EvtCand = new TopFCNC_Evt(TopFCNC_Evt::kMuon);
                MyTopFCNC_EvtCand->SetLepton1FromZ(*selectedMuons[idx_Z_mu]);
                MyTopFCNC_EvtCand->SetLepton2FromZ(*selectedElectrons[idx_Z_el]);
                MyTopFCNC_EvtCand->SetSelectedJets(selectedJets);
                MyTopFCNC_EvtCand->SetMET(*mets[0]);
              }
            }
          }
        }
      }

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
        Tree->Fill();
        delete MyTopFCNC_EvtCand;
      }
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
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////// Writing histograms //////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " - Writing outputs to the files ..." << endl;
  
  //Selection tables
  //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)	
  selecTableMuEl.TableCalculator(false, true, true, true, true);
  //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
  selecTableMuEl.Write("TopFCNC"+postfix+channelpostfix+comments+"_SelectionTable_MuEl.tex",    true,true,true,true,false,false,true);

  //MultiSample plots
  fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
	  MultiSamplePlot *temp = it->second;
	  //temp->addText("CMS preliminary");
	  string name = it->first;
	  name += comments;
	  temp->Draw(false, name, true, true, true, true, true,1,true); // merge TT/QCD/W/Z/ST/
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

double MuEffSF_Id_Run2012(string period, double eta, double pt){
  if(period.find("2012ABCD") != string::npos)
    return MuEffSF_Id_Run2012ABCD(eta,pt);
  else
    return 1.;
}
double MuEffSF_Iso04_Run2012(string period, double eta, double pt){
  if(period.find("2012ABCD") != string::npos)
    return MuEffSF_Iso04_Run2012ABCD(eta,pt);
  else
    return 1.;
}

double MuEffSF_Id_Run2012ABCD(double eta, double pt){
  if( fabs(eta) >= 0 && fabs(eta) < 0.9 ){
    if( pt >= 10 && pt < 20)
      return 0.975506;
    else if( pt >= 20 && pt < 25)
      return 0.986202;
    else if( pt >= 25 && pt < 30)
      return 0.99337;
    else if( pt >= 30 && pt < 35)
      return 0.996769;
    else if( pt >= 35 && pt < 40)
      return 0.997442;
    else if( pt >= 40 && pt < 50)
      return 0.998198;
    else if( pt >= 50 && pt < 60)
      return 0.995042;
    else if( pt >= 60 && pt < 90)
      return 0.987628;
    else if( pt >= 90 && pt < 140)
      return 0.998926;
    else if( pt >= 140 && pt < 300)
      return 1.00231;
    else
      return 1.;
  }
  else if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ){
    if( pt >= 10 && pt < 20)
      return 1.00453;
    else if( pt >= 20 && pt < 25)
      return 0.999235;
    else if( pt >= 25 && pt < 30)
      return 0.998903;
    else if( pt >= 30 && pt < 35)
      return 0.997493;
    else if( pt >= 35 && pt < 40)
      return 0.998633;
    else if( pt >= 40 && pt < 50)
      return 0.998107;
    else if( pt >= 50 && pt < 60)
      return 0.997018;
    else if( pt >= 60 && pt < 90)
      return 0.989447;
    else if( pt >= 90 && pt < 140)
      return 1.00105;
    else if( pt >= 140 && pt < 300)
      return 1.00431;
    else
      return 1.;
  }
  else if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ){
    if( pt >= 10 && pt < 20)
      return 1.00221;
    else if( pt >= 20 && pt < 25)
      return 1.00319;
    else if( pt >= 25 && pt < 30)
      return 0.998141;
    else if( pt >= 30 && pt < 35)
      return 0.998023;
    else if( pt >= 35 && pt < 40)
      return 0.99822;
    else if( pt >= 40 && pt < 50)
      return 0.99835;
    else if( pt >= 50 && pt < 60)
      return 0.994846;
    else if( pt >= 60 && pt < 90)
      return 0.984053;
    else if( pt >= 90 && pt < 140)
      return 1.0043;
    else if( pt >= 140 && pt < 300)
      return 0.985374;
    else
      return 1.;
  }
  else return 1.;
}
double MuEffSF_Iso04_Run2012ABCD(double eta, double pt){
  if( fabs(eta) >= 0 && fabs(eta) < 0.9 ){
    if( pt >= 10 && pt < 20)
      return 0.963644;
    else if( pt >= 20 && pt < 25)
      return 0.988643;
    else if( pt >= 25 && pt < 30)
      return 0.999383;
    else if( pt >= 30 && pt < 35)
      return 0.998656;
    else if( pt >= 35 && pt < 40)
      return 0.998344;
    else if( pt >= 40 && pt < 50)
      return 0.99817;
    else if( pt >= 50 && pt < 60)
      return 0.999064;
    else if( pt >= 60 && pt < 90)
      return 1.00044;
    else if( pt >= 90 && pt < 140)
      return 1.0008;
    else if( pt >= 140 && pt < 300)
      return 1.00188;
    else
      return 1.;
  }
  else if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ){
    if( pt >= 10 && pt < 20)
      return 0.963851;
    else if( pt >= 20 && pt < 25)
      return 0.9879;
    else if( pt >= 25 && pt < 30)
      return 1.00148;
    else if( pt >= 30 && pt < 35)
      return 1.00213;
    else if( pt >= 35 && pt < 40)
      return 1.00105;
    else if( pt >= 40 && pt < 50)
      return 1.00005;
    else if( pt >= 50 && pt < 60)
      return 0.999964;
    else if( pt >= 60 && pt < 90)
      return 1.00101;
    else if( pt >= 90 && pt < 140)
      return 1.00153;
    else if( pt >= 140 && pt < 300)
      return 0.998248;
    else
      return 1.;
  }
  else if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ){
    if( pt >= 10 && pt < 20)
      return 0.977632;
    else if( pt >= 20 && pt < 25)
      return 0.994509;
    else if( pt >= 25 && pt < 30)
      return 1.00207;
    else if( pt >= 30 && pt < 35)
      return 1.00289;
    else if( pt >= 35 && pt < 40)
      return 1.00185;
    else if( pt >= 40 && pt < 50)
      return 1.;
    else if( pt >= 50 && pt < 60)
      return 1.00011;
    else if( pt >= 60 && pt < 90)
      return 1.00037;
    else if( pt >= 90 && pt < 140)
      return 0.99977;
    else if( pt >= 140 && pt < 300)
      return 0.997083;
    else
      return 1.;
  }
  else return 1.;
}

double ElEffSF_Id_Run2012(string period, double eta, double pt){
  if(period.find("2012ABCD") != string::npos)
    return ElEffSF_Id_Run2012ABCD(eta,pt);
  else
    return 1.;
}
double ElEffSF_Iso04_Run2012(string period, double eta, double pt){
  if(period.find("2012ABCD") != string::npos)
    return ElEffSF_Iso04_Run2012ABCD(eta,pt);
  else
    return 1.;
}

double ElEffSF_Id_Run2012ABCD(double eta, double pt){
  if( fabs(eta) >= 0.000 && fabs(eta) < 0.800 ){
    if( pt >= 10 && pt < 15)
      return 0.818;
    else if( pt >= 15 && pt < 20)
      return 0.928;
    else if( pt >= 20 && pt < 30)
      return 0.973;
    else if( pt >= 30 && pt < 40)
      return 0.979;
    else if( pt >= 40 && pt < 50)
      return 0.984;
    else if( pt >= 50)// && pt < 200)
      return 0.983;
    else
      return 1.;
  }
  else if( fabs(eta) >= 0.800 && fabs(eta) < 1.442 ){
    if( pt >= 10 && pt < 15)
      return 0.84;
    else if( pt >= 15 && pt < 20)
      return 0.914;
    else if( pt >= 20 && pt < 30)
      return 0.948;
    else if( pt >= 30 && pt < 40)
      return 0.961;
    else if( pt >= 40 && pt < 50)
      return 0.972;
    else if( pt >= 50)// && pt < 200)
      return 0.977;
    else
      return 1.;
  }
  else if( fabs(eta) >= 1.442 && fabs(eta) < 1.556 ){
    if( pt >= 10 && pt < 15)
      return 1.008;
    else if( pt >= 15 && pt < 20)
      return 0.877;
    else if( pt >= 20 && pt < 30)
      return 0.983;
    else if( pt >= 30 && pt < 40)
      return 0.983;
    else if( pt >= 40 && pt < 50)
      return 0.957;
    else if( pt >= 50)// && pt < 200)
      return 0.978;
    else
      return 1.;
  }
  else if( fabs(eta) >= 1.556 && fabs(eta) < 2.000 ){
    if( pt >= 10 && pt < 15)
      return 0.906;
    else if( pt >= 15 && pt < 20)
      return 0.907;
    else if( pt >= 20 && pt < 30)
      return 0.957;
    else if( pt >= 30 && pt < 40)
      return 0.962;
    else if( pt >= 40 && pt < 50)
      return 0.985;
    else if( pt >= 50)// && pt < 200)
      return 0.986;
    else
      return 1.;
  }
  else if( fabs(eta) >= 2.000 && fabs(eta) < 2.500 ){
    if( pt >= 10 && pt < 15)
      return 0.991;
    else if( pt >= 15 && pt < 20)
      return 0.939;
    else if( pt >= 20 && pt < 30)
      return 1.001;
    else if( pt >= 30 && pt < 40)
      return 1.002;
    else if( pt >= 40 && pt < 50)
      return 0.999;
    else if( pt >= 50)// && pt < 200)
      return 0.995;
    else
      return 1.;
  }
  else return 1.;
}
double ElEffSF_Iso04_Run2012ABCD(double eta, double pt){
  return 0.9912;
}

double MuonRelIso(TRootMuon* muon){
  double reliso = ( muon->chargedHadronIso() + max( 0.0, muon->neutralHadronIso() + muon->photonIso() - 0.5*muon->puChargedHadronIso() ) ) / muon->Pt();
  return reliso;
}

double ElecRelIso(TRootElectron* elec, TRootEvent* event){
  if(!elec){
    cerr << "Pointer to electron is NULL" << endl;
    exit(1);
  }
  // HCP 2012 updated for electron conesize = 0.3
  // Taken from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h?revision=1.4&view=markup
  double EffectiveArea = 0.;
  if (fabs(elec->superClusterEta()) >= 0.0   && fabs(elec->superClusterEta()) < 1.0   ) EffectiveArea = 0.130;
  if (fabs(elec->superClusterEta()) >= 1.0   && fabs(elec->superClusterEta()) < 1.479 ) EffectiveArea = 0.137;
  if (fabs(elec->superClusterEta()) >= 1.479 && fabs(elec->superClusterEta()) < 2.0   ) EffectiveArea = 0.067;
  if (fabs(elec->superClusterEta()) >= 2.0   && fabs(elec->superClusterEta()) < 2.2   ) EffectiveArea = 0.089;
  if (fabs(elec->superClusterEta()) >= 2.2   && fabs(elec->superClusterEta()) < 2.3   ) EffectiveArea = 0.107;
  if (fabs(elec->superClusterEta()) >= 2.3   && fabs(elec->superClusterEta()) < 2.4   ) EffectiveArea = 0.110;
  if (fabs(elec->superClusterEta()) >= 2.4) EffectiveArea = 0.138;
  
  double isocorr = EffectiveArea*event->kt6PFJets_rho();
  
  double reliso = (elec->chargedHadronIso() + max(0.0 , elec->neutralHadronIso() + elec->photonIso() - isocorr) )/ elec->Pt();
  
  return reliso;
}

bool ZCandInvMass(vector<TRootMuon*>& muons, vector<TRootElectron*>& electrons, int& idx_Z_mu, int& idx_Z_el, double& invMass, double Zmass, double Zwindowsize){
  TRootMuon* mu = 0;
  TRootElectron* el = 0;
  double invMass_tmp = 0;
  for(unsigned int i=0;i<muons.size();i++)
  {
    for(unsigned int j=0;j<electrons.size();j++)
    {
      mu = (TRootMuon*) muons[i];
      el = (TRootElectron*) electrons[j];
      if(mu->charge() == el->charge()) continue;
      invMass_tmp = (*mu + *el).M();
      if(fabs(invMass_tmp-Zmass)<fabs(invMass-Zmass)){
        idx_Z_mu = i;
        idx_Z_el = j;
        invMass = invMass_tmp;
      }
    }
  }
  if(fabs(invMass-Zmass)<Zwindowsize)
    return true;
  else
    return false;
}