/////////////
////////////
// TODO + COMMENTS
//1.  Validate  all those electron ID variables...
//2a. Add JER systematic calculation, using GEN Jets.
//2b. Add btag scale factor systematic calculation, need to update number for appropriate taggers.
//2c. Add MET corrections
//3. Data-driven estimation of QCD component
//4. Change to 52X JES correction
//5. Object selection cuts to be reviewed
//6. Need  to do jet-parton matching to plot hadronic W mass and  estimate sigma(M_{w}), example in Petra's code.]
//7. Dcriminator values of 3/4 most 'b-like' jets, how do I measure efficiency?? Template fit?
//8. Test inputs for MVA/BDT/NN
//9. Limit setting, CLS interval, possibly likelihood ratio as test stat.
//10. Add e + jets final state, prob leave this till last.
///11.

#include "TStyle.h"
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
//#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../MCInformation/interface/LumiReWeighting.h"
//#include "../TopFCNC/interface/TopFCNC_Evt.h"

#include "../macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

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
  if(btagger == "TCHEL") //track counting high eff loose working point, need to replace numbers with the correct ones for
    // the CSV and TCHP taggers for 2012
  {
	  if(dobTagEffShift == 0)
		scalefactorbtageff = 0.95;
	  if(dobTagEffShift == 1)
		scalefactorbtageff = 0.85;
	  if(dobTagEffShift == 2)
		scalefactorbtageff = 1.05;
	  
	  if(domisTagEffShift == 0)
		mistagfactor = 1.11;
	  if(domisTagEffShift == 1)
		mistagfactor = 0.99;
	  if(domisTagEffShift == 2)
		mistagfactor = 1.23;
		
  }
  else if(btagger == "TCHEM") //track counting high eff medium working point
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
  float workingpointvalue = 9999; //{1.7,3.3,10.2}; trackcountinghighefficiency working points: loose, medium, tight
  if(btagger == "TCHPL")
     workingpointvalue = 1.470;
  else if(btagger == "TCHPM")
     workingpointvalue = 2.360;	
  else if(btagger == "TCHPT")
    workingpointvalue = 5.360;
  else if(btagger == "CSVL")
     workingpointvalue = .387;	
  else if(btagger == "CSVM")
    workingpointvalue = .838;
  else if(btagger == "CSVT")
    workingpointvalue = .940;

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

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string channelpostfix = "";
  string xmlFileName = "";

  //different final states, for now, just mu+jets and e+jets
  // possible add multileptons later ...maybe.

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

  //  xmlFileName = "../config/myFourthGenconfig_Muon_Fall11.xml";

  //   xmlFileName = "test_fullsamples.xml";
   xmlFileName = "test_2.xml";

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
  float Luminosity = 1000; //pb^-1??
  
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
  //  TopFCNC_Evt* MyTopFCNC_EvtCand = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ////////////////// MultiSample plots  //////////////////////////////
  ////////////////////////////////////////////////////////////////////

  MSPlot["RhoCorrection"]              = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
  MSPlot["NbOfVertices"]               = new MultiSamplePlot(datasets, "NbOfVertices", 50, 0, 50, "Nb. of vertices");

  MSPlot["NbOfIsolatedMuons"]          = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfIsolatedElectrons"]      = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");

  MSPlot["NbOfExtraIsolatedMuons"]     = new MultiSamplePlot(datasets, "NbOfExtraIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
  MSPlot["NbOfExtraIsolatedElectrons"] = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");

  //Muon ID variables
  MSPlot["MuonRelIsolation"] = new MultiSamplePlot(datasets, "MuonRelIsolation", 150, 0, .7, "RelIso");
  MSPlot["MuonPt"]              = new MultiSamplePlot(datasets, "MuonPt", 100, 0, 200, "PT_{#mu}");
  MSPlot["MuonEta"]              = new MultiSamplePlot(datasets, "MuonEta", 50, 0, 6, "#eta_{#mu}");
  MSPlot["MuonNValidHits"]              = new MultiSamplePlot(datasets, "MuonNValidHits", 30, 0, 30, "NValidHits_{#mu}");
  MSPlot["Muond0"]              = new MultiSamplePlot(datasets, "Muond0", 50, 0, .05, "d0_{#mu}");
  MSPlot["MuondRJets"]              = new MultiSamplePlot(datasets, "MuondRJets", 50, 0, 10, "dRJets_{#mu}");
  MSPlot["MuonNMatchedStations"]              = new MultiSamplePlot(datasets, "MuonNMatchedStations", 10, 0, 10, "NMatchedStations_{#mu}");
  MSPlot["MuonDistVzPVz"]              = new MultiSamplePlot(datasets, "MuonDistVzPVz", 75, -.3,.3, "DistVzPVz_{#mu}");
  MSPlot["MuonTrackerLayersWithMeasurement"]    = new MultiSamplePlot(datasets, "MuonTrackerLayersWithMeasurement", 25, 0, 25, "nLayers");

  MSPlot["BdiscBJetCand_CSV"]  = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_m_ch_CSV"]            = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_e_ch_CSV"]            = new MultiSamplePlot(datasets, "HighestBdisc_ee_ch_CVS", 100, 0, 1, "CSV b-disc.");
  MSPlot["HighestBdisc_m_ch_TCHE"]           = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_TCHE",100, 0, 50, "TCHE b-disc.");
  MSPlot["HighestBdisc_e_ch_TCHE"]           = new MultiSamplePlot(datasets, "HighestBdisc_ee_ch_TCHE",100, 0, 50, "TCHE b-disc.");
   
  MSPlot["NbOfSelectedJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedJets", 10, 0, 10, "Nb. of jets");
  MSPlot["NbOfSelectedLightJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
  MSPlot["NbOfSelectedBJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 10, 0, 10, "Nb. of jets");
  MSPlot["JetEta"]                  = new MultiSamplePlot(datasets, "JetEta", 50,-5, 5, "Jet #eta");
  MSPlot["JetPhi"]                  = new MultiSamplePlot(datasets, "JetPhi", 50, -4,4 , "Jet #phi");
  MSPlot["DiJetMasses_noBJets"] = new MultiSamplePlot(datasets, "DiJetMasses_noBJets", 50, 0, 1000, "m_{jj}");
  MSPlot["DiJetMass_chi2"] = new MultiSamplePlot(datasets, "DiJetMass_chi2", 20, 0, 20, "\\chi^{2}");
  MSPlot["TriJetMasses_1BJet"] = new MultiSamplePlot(datasets, "TriJetMasses_1BJet", 50, 0, 1000, "m_{jjj}");
  MSPlot["BestTriJetMass_1BJet"] = new MultiSamplePlot(datasets, "BestTriJetMass_1BJet", 50, 0, 1000, "m_{jjj}");
  MSPlot["TriJetMass_chi2"] = new MultiSamplePlot(datasets, "TriJetMass_chi2", 50, 0, 1000, "\\chi^{2}");
  MSPlot["WMt"] = new MultiSamplePlot(datasets, "WMt", 50, 0, 250, "W Transverse Mass");
  MSPlot["MuMetBMasses"] = new MultiSamplePlot(datasets, "MuMetBMasses", 50, 0, 1000, "m_{muMETj}");
  MSPlot["MuMetBMasses_chi2"] = new MultiSamplePlot(datasets, "MuMetBMasses_chi2", 50, 0, 1000, "\\chi^{2}");
  MSPlot["NbOf_WCands"]                  = new MultiSamplePlot(datasets, "NbOf_WCands", 5, 0, 5, "Nb. of W candidates");
  MSPlot["NbOf_HadTops"]                  = new MultiSamplePlot(datasets, "NbOf_HadTops", 5, 0, 5, "Nb. of Hadronic Tops");
  MSPlot["SummedDiscr_QuadJet_CSV"]   = new MultiSamplePlot(datasets, "SummedDiscr_QuadJet_CSV", 50, 0, 5, "SummedDiscr_QuadJet_CSV");
  MSPlot["SummedDiscr_TriJet_CSV"]    = new MultiSamplePlot(datasets, "SummedDiscr_TriJet_CSV", 50, 0, 5, "SummedDiscr_TriJet_CSV");
  MSPlot["SelectedJetPt"] = new MultiSamplePlot(datasets, "JetPt", 50, 0, 1000, "PT_{jet}");
  MSPlot["SelectedJetPt_light"] = new MultiSamplePlot(datasets, "JetPt_light", 50, 0, 1000, "PT_{lightjet}");
  MSPlot["SelectedJetPt_b"] = new MultiSamplePlot(datasets, "JetPt_b", 50, 0, 1000, "PT_{bjet}");
  MSPlot["HT_SelectedJets"] = new MultiSamplePlot(datasets, "HT_SelectedJets", 50, 0, 1000, "HT");
  MSPlot["MHT_SelectedJets"] = new MultiSamplePlot(datasets, "MHT_SelectedJets", 75, 0, 1000, "MHT");
  MSPlot["MHTSig_SelectedJets"] = new MultiSamplePlot(datasets, "MHTSig_SelectedJets", 75, 0, 30, "MHTSig");
  MSPlot["MET"] = new MultiSamplePlot(datasets, "MET", 75, 0, 1000, "MET");
  MSPlot["MET_MHT"]= new MultiSamplePlot(datasets, "MET_MHT", 75, 0, 200, "MET_MHT");
  MSPlot["STLep"] = new MultiSamplePlot(datasets, "STLep", 75, 0, 1000, "STLep");
  MSPlot["ElectronPt"]              = new MultiSamplePlot(datasets, "ElectronPt", 50, 0, 100, "PT_{e}");
 
  MSPlot["DiMuon_InvMass"]     = new MultiSamplePlot(datasets, "DiMuon_InvMass", 60, 0, 120, "DiMuon_InvMass");
  MSPlot["NbOfLooseMuon"]     = new MultiSamplePlot(datasets, "NbOfLooseMuon", 10, 0, 10, "Nb. of loose muons");
  MSPlot["NbOfLooseElectron"] = new MultiSamplePlot(datasets, "NbOfLooseElectron", 10, 0, 10, "Nb. of loose electrons");

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
  for (unsigned int d = 0; d < datasets.size(); d++){
    histo2D[("RelIso_vs_MET_"+datasets[d]->Name()).c_str()] = new TH2F(("RelIso_vs_MET_"+datasets[d]->Name()).c_str(),"RelIso:MET",100,0,1000, 100, 0,1);

	//	histo2D[("MET_vs_Mzq_mm_ch_"+datasets[d]->Name()).c_str()] = new TH2F(("MET_vs_Mzq_mm_ch_"+datasets[d]->Name()).c_str(),"MET:m_{zq}",100,0,200,100,80,300);
	//histo2D[("MET_vs_Mzq_mmm_ch_"+datasets[d]->Name()).c_str()] = new TH2F(("MET_vs_Mzq_mmm_ch_"+datasets[d]->Name()).c_str(),"MET:m_{zq}",100,0,200,100,80,300);
	//	histo2D[("MET_vs_Mzq_mme_ch_"+datasets[d]->Name()).c_str()] = new TH2F(("MET_vs_Mzq_mme_ch_"+datasets[d]->Name()).c_str(),"MET:m_{zq}",100,0,200,100,80,300);
  }
  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "FourTop"+postfix+channelpostfix;
  pathPNG += "_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Selection Tables ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : µ + jets  //////////////////////////
  ////////////////////////////////////////////////////////////////////
  vector<string> CutsSelecTableMu;
  CutsSelecTableMu.push_back(string("initial"));
  CutsSelecTableMu.push_back(string("PU reweighting"));
  CutsSelecTableMu.push_back(string("Trigger"));
  CutsSelecTableMu.push_back(string("Good PV"));
  CutsSelecTableMu.push_back(string("$\\geq$ 2 isolated muons"));
  CutsSelecTableMu.push_back(string("$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableMu.push_back(string("Veto on 3rd iso. lept."));
  CutsSelecTableMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableMu.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableMu.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableMu.push_back(string("$b\\texttt{-}disc \\geq 0.7$ (CSV)"));
//  CutsSelecTableMu.push_back(string("$MET \\leq 30$ GeV"));  

  SelectionTable selecTableMu(CutsSelecTableMu, datasets);
  selecTableMu.SetLuminosity(Luminosity);
  selecTableMu.SetPrecision(1);


  ////////////////////////////////////////////////////////////////////
  ///////////////////// Channel : e + jets    ////////////////////////
  ////////////////////////////////////////////////////////////////////

  vector<string> CutsSelecTableEl;
  CutsSelecTableEl.push_back(string("initial"));
  CutsSelecTableEl.push_back(string("PU reweighting"));
  CutsSelecTableEl.push_back(string("Trigger"));
  CutsSelecTableEl.push_back(string("Good PV"));
  CutsSelecTableEl.push_back(string("$\\geq$ 2 isolated electrons"));
  CutsSelecTableEl.push_back(string("$|m_{ll}-m_Z|<30$ GeV"));
  CutsSelecTableEl.push_back(string("Veto on 3rd iso. lept."));
  //CutsSelecTableEl.push_back(string("Conversion veto"));
  CutsSelecTableEl.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableEl.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableEl.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableEl.push_back(string("$\\geq$ 4 jet"));
  CutsSelecTableEl.push_back(string("$b\\texttt{-}disc \\geq 0.7$ (CSV)"));
  //CutsSelecTableEl.push_back(string("MET > 40 GeV"));


  SelectionTable selecTableEl(CutsSelecTableEl, datasets);
  selecTableEl.SetLuminosity(Luminosity);
  selecTableEl.SetPrecision(1);
  
  cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// PileUp Reweighting - N True interactions, recommended method for 2012 ///////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  //OLD METHOD (3D)
  //  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  //  Lumi3DWeights.weight3D_init(1.0);
  //
  // if(doPUShift == 1)
  //	Lumi3DWeights.weight3D_init(0.92);
  //else if(doPUShift == 2)
  //	Lumi3DWeights.weight3D_init(1.08);
  // else
  // 	Lumi3DWeights.weight3D_init(1.0);

  //  cout << " - Initialized LumiReWeighting stuff" << endl;
  

  //NEW METHOD (TRUE INTERACTIONS)


  LumiReWeighting LumiWeights;
  
  LumiWeights = LumiReWeighting("../macros/PileUpReweighting/pileup_MC_Summer12.root", "../macros/PileUpReweighting/pileup_2012Data_UpToRun191810.root", "pileup", "pileup");

  reweight::PoissonMeanShifter PShiftDown_ = reweight::PoissonMeanShifter(-0.6);
  reweight::PoissonMeanShifter PShiftUp_ = reweight::PoissonMeanShifter(0.6);

  cout << " Initialized LumiReWeighting stuff" << endl;
  
  //exit(1);

  //Get fitted widths of hadronic W and top distributions
  //
  //
  //These dummy sigmas need to replaced with the fitted sigma of the
  // appropriate truth-matched dijet and trijet mass distributions. Maybe the sigma for
  // lepton + bJet + MET objects is needed also.
	double sigW = 10;
	double sigTop = 30;



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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////// Loop on events //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int itrigger = -1, previousRun = -1;
    int fourIsoLeptCounter = 0;
      
    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();

    bool debug = true;

    if (verbose > 1) cout << " - Loop over events " << endl;      
    
    for (unsigned int ievt = start; ievt < end; ievt++)
    {        

	if(ievt%1000 == 0)
		std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

	//load event
	event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

	vector<TRootGenJet*> genjets;
	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
	  // loading GenJets until I need them for JER
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

				     cout <<" checking trig data "<<endl;


			  if(event->runId() >= 190456 && event->runId() <= 190761)
				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
			  else if ( event->runId() >= 190762 && event->runId() <= 191511 )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v3"), currentRun, iFile);
			  else if ( event->runId() >= 191512  && event->runId() <= 193805  )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v4"), currentRun, iFile);
			  //for  top triggers change PD at this point, seems no trigger for run 194305??

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

				if(dataSetName == "TTJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
				else if (dataSetName == "WJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);

                                else if (dataSetName == "ZJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v3"), currentRun, iFile);

				else if (dataSetName == "SingleTop") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);

				else if (dataSetName == "SingleToptW_T") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
				else if (dataSetName == "SingleToptW_TBar") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);

				else if (dataSetName == "MultiJet") {

	                 	  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v4"), currentRun, iFile); 
				  cout <<" found trig "<<  itrigger   <<endl;


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
					//					if (dataSetName != "MultiJet")
					exit(1);
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
	  //JK
	  doJERShift == 0;
	  //JK
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
		// Apply pile-up reweighting
	  //	double lumiWeight3D = Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
	  //	scaleFactor *= lumiWeight3D;
	  //replacing 3D reweighting with N true interactions...


	double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
	double lumiWeightOLD=lumiWeight;
	if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
	  lumiWeight=1;

	//if (event->nPu(0) >= 25) cout << "#pu " << event->nPu(0) << " shift +0.6 " << PShiftUp_.ShiftWeight( event->nPu(0) ) << endl;

	//will enable flags here for PU systematic calculation
	//	if (systematic == 10) //less PU
	//	  lumiWeight = lumiWeight*PShiftUp_.ShiftWeight( event->nPu(0) );
	
	  //	  else if (systematic == 11) //more PU
	// lumiWeight = lumiWeight*PShiftDown_.ShiftWeight( event->nPu(0) );

	//cout << "scaleFactor before = " << scaleFactor << " " << lumiWeight <<endl;
	scaleFactor = scaleFactor*lumiWeight;

	}

	histo1D["lumiWeights"]->Fill(scaleFactor);	
	MSPlot["RhoCorrection"]->Fill(event->kt6PFJetsPF2PAT_rho(), datasets[d], true, Luminosity*scaleFactor);
			
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////// Event selection ////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	selecTableMu.Fill(d,0, 1.);//datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
	selecTableEl.Fill(d,1,scaleFactor);
	       
	// Apply trigger selection
	trigged = treeLoader.EventTrigged (itrigger);
	cout<<"triggered? Y/N?  "<< trigged  <<endl;

	//Applying trigger selection again with 2012 Muon+Jets triggers.
	if(!trigged)		   continue;
	//Turning off until I find right triggers for secondary bbkgs

	selecTableMu.Fill(d,2,scaleFactor);
	selecTableEl.Fill(d,2,scaleFactor);

	// Declare selection instance    
	Selection selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...

	//Object selection cuts to be revised

	// Define object selection cuts
	selection.setJetCuts(45.,2.5,0.01,1.,0.98,0.3,0.1);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets

	selection.setElectronCuts(20.,2.5,0.17,0.02,1.,1,0.3); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets
	selection.setLooseElectronCuts(15.,2.5,0.2,0.2);

	selection.setMuonCuts(20.,2.4,0.125,10,0.02,0.3,1,1.,1 ); //Pt,Eta,RelIso,NValidHits,d0, dRJets, NMatchedStations,DistVzPVz,NTrackerLayersWithMeas 
	selection.setLooseMuonCuts(15,2.4,0.2);
	  
	//Select objects 
	vector<TRootElectron*> selectedElectrons_NoIso = selection.GetSelectedElectrons(20,2.4,999.,vertex[0]);
	vector<TRootElectron*> selectedElectrons       = selection.GetSelectedElectrons(vertex[0]);
	vector<TRootElectron*> selectedExtraElectrons;
	vector<TRootMuon*>     selectedMuons_NoIso = selection.GetSelectedMuons(20,2.4,999.);
	vector<TRootMuon*>     selectedMuons       = selection.GetSelectedMuons();
	vector<TRootMuon*>     selectedExtraMuons;
	vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId
        vector<TRootJet*>      selectedSoftJets        = selection.GetSelectedJets(20.,2.5, selectedMuons, 0., true); // ApplyJetId

	vector<TRootMuon*>     selectedLooseMuons     = selection.GetSelectedLooseMuons();
	vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseElectrons(true); // VBTF Id

	sort(selectedMuons.begin(),selectedMuons.end(),HighestPt()); //order muons wrt Pt.

	//Will impement more elegant retrieval of container of b-jets, for now, will take from selected jets
	//        vector<TRootJet*>      selectedBJets   =  selection.GetSelectedBJets(selectedJets,CSV, ); // B-Jets

        vector<TRootJet*>      selectedBJets; // B-Jets
        vector<TRootJet*>      selectedLightJets; // light-Jets

	for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ ){
	  if( selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue) selectedBJets.push_back(selectedJets[seljet]);
	  else selectedLightJets.push_back(selectedJets[seljet]);
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Applying baseline offline event selection here
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (debug)	cout <<" applying baseline event selection..."<<endl;

	// Apply primary vertex selection
	bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
       if(!isGoodPV) continue;
	selecTableMu.Fill(d,3,scaleFactor);
	selecTableEl.Fill(d,3,scaleFactor);

	if (debug) cout <<"good vertex.."<<endl;

	if(Muon){

	  int nMu = selectedMuons.size();

	  if (debug) cout <<"Number of Muons, Jets, BJets  ===>  "<< selectedMuons.size() <<"  "  << selectedJets.size()   <<"  " <<  selectedBJets.size()   <<endl;

	  if  (  !( selectedJets.size() >= 3 && selectedBJets.size() >= 1 && selectedSoftJets.size() >= 1 &&     nMu == 1   )) continue; 
	
	  //add cut on dimu inv mass to supress Z+jets...maybe not neccessary now


	  //Order of plotting
	  // 0. Vertices
	  // 1. Muons
	  // 2. Jets: per jet plots, event level variables, jet combinations,summed discriminants.


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Vertices
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Muons
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  if (debug) cout <<"muons.."<<  selectedMuons.size()  <<"   "<< selectedLooseMuons.size()   <<endl;
	  for (Int_t mu1 = 0; mu1 < selectedMuons.size(); mu1++){
	if (debug) cout <<"first muon loop.."<<endl;
	  for (Int_t mu2 = 0; mu2 < selectedLooseMuons.size(); mu2++){
	if (debug) cout <<"second muon loop.."<<endl;
	    TRootMuon DiMu (TLorentzVector ( selectedLooseMuons[mu1]->Px() + selectedLooseMuons[mu2]->Px() ,selectedLooseMuons[mu1]->Py() + selectedLooseMuons[mu2]->Py()    ,selectedLooseMuons[mu1]->Pz() + selectedLooseMuons[mu2]->Pz() , selectedLooseMuons[mu1]->E() + selectedLooseMuons[mu2]->E() ));

	if (debug) cout <<"made di-muon .."<<endl;

	    //fill plots for ABCD here...
	    histo2D[("RelIso_vs_MET_"+datasets[d]->Name()).c_str()]->Fill(mets[0]->Et(),selectedMuons_NoIso[0]->relativePfIso03(), Luminosity*scaleFactor );

	if (debug) cout <<"filled ABCD plot  .."<<endl;


	    MSPlot["DiMuon_InvMass"]->Fill( DiMu.M(), datasets[d], true, Luminosity*scaleFactor );  
	  }
	  }
	  MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);

	  //Fill Muon ID plots
	  MSPlot["MuonRelIsolation"]->Fill(selectedMuons_NoIso[0]->relativePfIso03(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonPt"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonEta"]->Fill(selectedMuons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonNValidHits"]->Fill(selectedMuons[0]->nofValidHits(), datasets[d], true, Luminosity*scaleFactor);	 
	  MSPlot["Muond0"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);

	  MSPlot["MuonNMatchedStations"]->Fill(selectedMuons[0]->nofMatchedStations(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["MuonTrackerLayersWithMeasurement"]->Fill(selectedMuons[0]->nofTrackerLayersWithMeasurement(), datasets[d], true, Luminosity*scaleFactor);
	  double muzPVz = selectedMuons[0]->vz() - vertex[0]->Z();
          MSPlot["MuonDistVzPVz"]->Fill(muzPVz, datasets[d], true, Luminosity*scaleFactor );

	  //          double   mumpx =   selectedMuons[0]->Px() + mets[0]->Px();
	  //double   mumpy =   selectedMuons[0]->Py() + mets[0]->Py();
	  // double   mumpz =   selectedMuons[0]->Pz() + mets[0]->Pz();
	  //double   mumpe =   selectedMuons[0]->E()  + mets[0]->E();
	  //          TRootParticle mum(TLorentzVector (mumpx, mumpy, mumpz,mumpe ));

	  double Mtrans =  (*selectedMuons[0] + *mets[0]).Mt();
	  MSPlot["WMt"]->Fill(Mtrans,datasets[d], true, Luminosity*scaleFactor);

	if (debug) cout <<"filled all muID plots  .."<<endl;
   
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Per jet plots
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
          MSPlot["NbOfSelectedBJets"]->Fill(selectedBJets.size(), datasets[d], true, Luminosity*scaleFactor);

	if (debug) cout <<"per jet plots.."<<endl;

	double ljetpt;
	double bjetpt;
	double jetpt;
	double HT =0;
	double MHT =0;
	double sumpx =0, sumpy=0, sumpz=0, sume=0; 
	

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
        MSPlot["MHT_SelectedJets"]->Fill(MHT, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
	MSPlot["MHTSig_SelectedJets"]->Fill(MHTSig, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MET_MHT"]->Fill(mets[0]->Et()/MHT, datasets[d], true, Luminosity*scaleFactor);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Jet combinations, keeping it simple for now.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (debug) cout <<"jet comb."<<endl;
	int nHadTops =0;
	int nWCand   =0;
	int bestMass = 0;
	double topmass = 176;
	//These dummy sigmas need to replaced with the fitted sigma of the
	// appropriate dijet and trijet mass distributions. Maybe the sigma for
	// lepton + bJet + MET objects is needed also.
	double sigW = 25;
	double sigTop = 50;
	for (Int_t jet1 =0; jet1 < selectedLightJets.size(); jet1++ ){
	  for (Int_t jet2 =0; jet2 < selectedLightJets.size(); jet2++ ){
	    TRootJet* j1 = (TRootJet*) selectedLightJets[jet1];
	    TRootJet* j2 = (TRootJet*) selectedLightJets[jet2];
	    double dijetMass = (*j1 + *j2).M();
	    
	    	    MSPlot["DiJetMasses_noBJets"]->Fill(dijetMass, datasets[d], true, Luminosity*scaleFactor);
	   	    if (fabs(( 80 - dijetMass )/ sigW) <  1 ) nWCand++;

	    //calc trijet mass here, restrict to one candidate per event (closest to top mass), then use other B-Jet for MuMetB mass
	    for (Int_t bjet1 =0; bjet1 < selectedBJets.size(); bjet1++ ){
	      TRootJet* bj1 = (TRootJet*) selectedBJets[bjet1];
	      double trijetMass = (*j1 + *j2 + *bj1).M();
	      double cand_chi2 = fabs(( 176 - trijetMass )/ sigTop) +  fabs(( 80 - dijetMass )/ sigW);
	      MSPlot["DiJetMass_chi2"]->Fill(cand_chi2, datasets[d], true, Luminosity*scaleFactor );
	      MSPlot["TriJetMasses_1BJet"]->Fill(trijetMass, datasets[d], true, Luminosity*scaleFactor);
	      if (fabs(( 176 - trijetMass )/ sigTop) <  1 && fabs(( 80 - dijetMass )/ sigW) <  1   ) nHadTops++;
	      if ( fabs(trijetMass - topmass) < fabs( bestMass - topmass ) ) bestMass = trijetMass;
 
	    }
	}
	}//end outer dijet loop
	MSPlot["BestTriJetMass_1BJet"]->Fill(bestMass,  datasets[d], true, Luminosity*scaleFactor );
	MSPlot["NbOf_WCands"]->Fill(nWCand, datasets[d], true, Luminosity*scaleFactor);
	MSPlot["NbOf_HadTops"]->Fill(nHadTops, datasets[d], true, Luminosity*scaleFactor);   

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Summed discriminant calculations
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 	float sumdiscr_3 =0;
        float sumdiscr_4 =0;

	if (debug) cout <<"summed discr.."<<endl;

	   if (selectedJets.size() >= 4 ){


     	   sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag()); //order jets wrt CSV discr value...
	   //           sumdiscr_3 = selectedJets[0]->btag_combinedSecondaryVertexBJetTags();

	   //make sure not to include default vals of discr.
	   for (Int_t ordJet =0; ordJet < 3; ordJet++ ){
	     if(selectedJets[ordJet]->btag_combinedSecondaryVertexBJetTags() <0 ) continue;
	     sumdiscr_3 = sumdiscr_3 + selectedJets[ordJet]->btag_combinedSecondaryVertexBJetTags();
	   }

	   for (Int_t ordJet2 =0; ordJet2 < 4; ordJet2++ ){
	     if(selectedJets[ordJet2]->btag_combinedSecondaryVertexBJetTags() <0 ) continue;
	     sumdiscr_4 = sumdiscr_4 + selectedJets[ordJet2]->btag_combinedSecondaryVertexBJetTags();
	   }

           MSPlot["SummedDiscr_QuadJet_CSV"]->Fill(sumdiscr_4,datasets[d], true, Luminosity*scaleFactor);
           MSPlot["SummedDiscr_TriJet_CSV"]->Fill(sumdiscr_3,datasets[d], true, Luminosity*scaleFactor);

	   if( sumdiscr_3 >  sumdiscr_4){
           cout <<"Sorted CSV values:   "<< selectedJets[0]->btag_combinedSecondaryVertexBJetTags()     <<  "  " << selectedJets[1]->btag_combinedSecondaryVertexBJetTags()  <<"  "<< selectedJets[2]->btag_combinedSecondaryVertexBJetTags()  <<"  "<< selectedJets[3]->btag_combinedSecondaryVertexBJetTags()<<endl;

	   cout <<"Summed discr. values:  TriJet =   "<<   sumdiscr_3  <<  "  QuadJet = "<<  sumdiscr_4 <<endl;
	   }
	   }                          


	   ///Can now fill MuMETB mass, could do this earlier, for tidier code, need to get some method
	   //of calculating mets[0].Pz()

	   double   mumbpx =   selectedJets[0]->Px() + selectedMuons[0]->Px() + mets[0]->Px();
	   double   mumbpy =   selectedJets[0]->Py() + selectedMuons[0]->Py() + mets[0]->Py();
	   double   mumbpz =   selectedJets[0]->Pz() + selectedMuons[0]->Pz() + mets[0]->Pz();
	   double   mumbpe =   selectedJets[0]->E()  + selectedMuons[0]->E()  + mets[0]->E();

           TRootParticle mumb(TLorentzVector (mumbpx, mumbpy, mumbpz,mumbpe ));

	   MSPlot["MuMetBMasses"]->Fill(mumb.M(), datasets[d], true, Luminosity*scaleFactor );

	}
	else if(Electron){
		MSPlot["1stLeadingElectronRelIsolation"]->Fill(selectedElectrons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);
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

  //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
	selecTableMu.Write(  "FourTop"+postfix+"Table_Mu.tex",    true,true,true,true,false,false,true);

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
	//	temp->addText("CMS preliminary");
	string name = it->first;
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

