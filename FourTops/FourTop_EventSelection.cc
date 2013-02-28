/////////////
////////////
// TODO + COMMENTS
//1.  Validate  all those electron ID variables...
//2a. Add JER systematic calculation, using GEN Jets.
//2b. Add btag scale factor systematic calculation, need to update number for appropriate taggers.
//3. Change to 52X JES correction
//4. Need  to do jet-parton matching to plot hadronic W mass and  estimate sigma(M_{w}), example in Petra's code.]
//5. Add e + jets final state, prob leave this till last.
//6. The setElectroncuts, setmuoncuts and getselectedelectrons methods need to be updated with the correct arguements for 53X_v4

#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"
#include "TRandom.h"

#include <iostream>
#include <map>

#include <cstdlib> 

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"

#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"
//#include "TopTreeAnalysis/InclFourthGenSearch/interface/InclFourthGenTree.h"
//#include "TopTreeAnalysis/InclFourthGenSearch/interface/InclFourthGenSearchTools.h"


#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"


//#include "interface/FourTopTree.h"

#include "TopTreeAnalysis/macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

bool split_ttbar = false;

pair<float, vector<unsigned int> > MVAvals1;	
pair<float, vector<unsigned int> > MVAvals2;	
pair<float, vector<unsigned int> > MVAvals2ndPass;
int nMVASuccesses=0;
int nMatchedEvents=0;


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


int Factorial(int N);


int main (int argc, char *argv[])
{

    BTagWeightTools * bTool = new BTagWeightTools("SFb-pt_payload_Moriond13.txt", "CSVM") ;
    
    TRandom *rand = new TRandom();

  
    
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
//  setTDRStyle();
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
     // xmlFileName = "config/test_dimuon.xml";
     //xmlFileName = "config/test_2.xml";
      // xmlFileName = "config/refsel.xml";

  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////// AnalysisEnvironment////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = 2;//anaEnv.Verbose;


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  float Luminosity = 5343.64; //pb^-1??
  vector<string> MVAvars;

    
  //A few bools to steer the MassReco and Event MVAs
  bool Tprime = false; // If false, regular variables are used in MVA
  string MVAmethod = "BDT"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)
  bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff
  bool trainEventMVA = false; // If false, the previously trained MVA will be used to calculate stuff
  bool computeEventMVA = false;
    
  JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, Tprime);
    
    double bestTopMass1 =0.;
    double bestTopMass2 = 0.;
    double bestTopMass2ndPass = 0.;
    double bestTopPt =0.;
  //  MVATrainer* Eventtrainer_ =0;
    MVAComputer* Eventcomputer_ =0;
    
    MVATrainer* Eventtrainer_ = new MVATrainer("BDT","EventMVA", "EventMVA.root");

    
    if(trainEventMVA){
        cout<<"instantiating trainer..."<<endl;
    
        cout<<"instantiated trainer..."<<endl;
        
    Eventtrainer_->bookInputVar("H");
    Eventtrainer_->bookInputVar("HT");
    Eventtrainer_->bookInputVar("HTX");
    Eventtrainer_->bookInputVar("HTH");
    Eventtrainer_->bookInputVar("HTXHX");
    Eventtrainer_->bookInputVar("EventMass");
    Eventtrainer_->bookInputVar("EventMassX");
    Eventtrainer_->bookInputVar("SumJetMass");
    Eventtrainer_->bookInputVar("SumJetMassX");
    Eventtrainer_->bookInputVar("PTBalTopEventX");
    Eventtrainer_->bookInputVar("PTBalTopSumJetX");
    Eventtrainer_->bookInputVar("bestTopMass2");
    
    Eventtrainer_->bookInputVar("nJets");
    Eventtrainer_->bookInputVar("bTagRatio13");
    Eventtrainer_->bookInputVar("bTagRatio14");
    
    Eventtrainer_->bookInputVar("Jet1Pt");
    Eventtrainer_->bookInputVar("Jet2Pt");
    Eventtrainer_->bookInputVar("Jet3Pt");
    Eventtrainer_->bookInputVar("Jet4Pt");
    Eventtrainer_->bookInputVar("Jet5Pt");
    Eventtrainer_->bookInputVar("Jet6Pt");
    }else if (computeEventMVA){
    
    MVAvars.push_back("H");
    MVAvars.push_back("HT");
    MVAvars.push_back("HTX");
    MVAvars.push_back("HTH");
    MVAvars.push_back("HTXHX");
    MVAvars.push_back("EventMass");
    MVAvars.push_back("EventMassX");
    MVAvars.push_back("SumJetMass");
    MVAvars.push_back("SumJetMassX");
    MVAvars.push_back("PTBalTopEventX");
    MVAvars.push_back("PTBalTopSumJetX");
    MVAvars.push_back("bestTopMass2");
    MVAvars.push_back("nJets");
    MVAvars.push_back("bTagRatio13");
    MVAvars.push_back("bTagRatio14");
    MVAvars.push_back("Jet1Pt");
    MVAvars.push_back("Jet2Pt");
    MVAvars.push_back("Jet3Pt");
    MVAvars.push_back("Jet4Pt");
    MVAvars.push_back("Jet5Pt");
    MVAvars.push_back("Jet6Pt");
    cout << " Initialized Eventcomputer_" << endl;

    MVAComputer* Eventcomputer_ = new MVAComputer("BDT","EventMVA.root","EventMVA",MVAvars, "test");

    }

        
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

    // for splitting the ttbar sample, it is essential to have the ttjets sample as the last
    //dataset loaded
    if (split_ttbar){
        int ndatasets = datasets.size() - 1 ;
        cout << " - splitting TTBar dataset ..." << ndatasets   << endl;
        vector<string> ttbar_filenames = datasets[ndatasets]->Filenames();
        cout <<"ttbar filenames =  "<< ttbar_filenames[0] <<endl;
        
        Dataset* ttbar_ll = new Dataset("TTJets_ll","tt + ll" , true, 633, 2, 2, 1, 213.4,ttbar_filenames );
        Dataset* ttbar_cc = new Dataset("TTJets_cc","tt + cc" , true, 633, 2, 2, 1, 6.9, ttbar_filenames );
        Dataset* ttbar_bb = new Dataset("TTJets_bb","tt + bb" , true, 633, 2, 2, 1, 4.8, ttbar_filenames );
        
        ttbar_ll->SetEquivalentLuminosity(223000.0);
        ttbar_cc->SetEquivalentLuminosity(223000.0);
        ttbar_bb->SetEquivalentLuminosity(223000.0);
        
        ttbar_ll->SetColor(kRed);
        ttbar_cc->SetColor(kRed-3);
        ttbar_bb->SetColor(kRed+2);
        
        
        datasets.pop_back();
        datasets.push_back(ttbar_bb);
        datasets.push_back(ttbar_cc);
        datasets.push_back(ttbar_ll);     
    }     
    
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

    
    ///Histos
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////
    /////////////////////////////////
    ////////////////// MultiSample plots  //////////////////////////////
    ////////////////////////////////////////////////////////////////////
    
    MSPlot["RhoCorrection"]              = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
    MSPlot["NbOfVertices"]               = new MultiSamplePlot(datasets, "NbOfVertices", 40, 0, 40, "Nb. of vertices");
    
    //Muons
    MSPlot["NbOfIsolatedMuons"]          = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
    MSPlot["NbOfIsolatedElectrons"]      = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["NbOfExtraIsolatedMuons"]     = new MultiSamplePlot(datasets, "NbOfExtraIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
    MSPlot["NbOfExtraIsolatedElectrons"] = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["MuonRelIsolation"] = new MultiSamplePlot(datasets, "MuonRelIsolation", 50, 0, .25, "RelIso");
    MSPlot["MuonRelIsolation_PreTrig"] = new MultiSamplePlot(datasets, "MuonRelIsolation_PreTrig", 15, 0, .25, "RelIso");
    MSPlot["MuonRelIsolation_PostTrig"] = new MultiSamplePlot(datasets, "MuonRelIsolation_PreTrig", 15, 0, .25, "RelIso");

    
    MSPlot["MuonPt"]              = new MultiSamplePlot(datasets, "MuonPt", 75, 0, 300, "PT_{#mu}");

    MSPlot["TriggerMuonPt"]              = new MultiSamplePlot(datasets, "MuonPt", 75, 0, 300, "PT_{#mu}");

    
    MSPlot["MuonPt_PreTrig"]              = new MultiSamplePlot(datasets, "MuonPt_PreTrig", 75, 0, 300, "PT_{#mu}");
    MSPlot["MuonPt_PostTrig"]              = new MultiSamplePlot(datasets, "MuonPt_PostTrig", 75, 0, 300, "PT_{#mu}");

    
    MSPlot["Muond0_PreTrig"]              = new MultiSamplePlot(datasets, "Muond0_PreTrig", 50, 0, .05, "d0_{#mu}");
    MSPlot["Muond0_PostTrig"]              = new MultiSamplePlot(datasets, "Muond0_PostTrig", 50, 0, .05, "d0_{#mu}");
    
    MSPlot["MuonEta"]              = new MultiSamplePlot(datasets, "MuonEta", 25, -2.4, 2.4, "#eta_{#mu}");
    MSPlot["MuonPhi"]              = new MultiSamplePlot(datasets, "MuonPhi", 50, -4, 4, "#phi_{#mu}");
    MSPlot["MuonNValidHits"]              = new MultiSamplePlot(datasets, "MuonNValidHits", 30, 0, 30, "NValidHits_{#mu}");
    MSPlot["Muond0"]              = new MultiSamplePlot(datasets, "Muond0", 50, 0, .05, "d0_{#mu}");
    MSPlot["MuondZPVz"]              = new MultiSamplePlot(datasets, "MuondZPVz", 50, 0, .5, "dZPVZ_{#mu}");

    MSPlot["MuondRJets"]              = new MultiSamplePlot(datasets, "MuondRJets", 50, 0, 10, "dRJets_{#mu}");
    MSPlot["MuonNMatchedStations"]              = new MultiSamplePlot(datasets, "MuonNMatchedStations", 10, 0, 10, "NMatchedStations_{#mu}");
    MSPlot["MuonDistVzPVz"]              = new MultiSamplePlot(datasets, "MuonDistVzPVz", 50, 0 ,.3, "DistVzPVz_{#mu}");
    MSPlot["MuonDz"]              = new MultiSamplePlot(datasets, "MuonDz", 25, -.6 ,.6, "Dz_{#mu}");

    MSPlot["MuonTrackerLayersWithMeasurement"]    = new MultiSamplePlot(datasets, "MuonTrackerLayersWithMeasurement", 25, 0, 25, "nLayers");
    MSPlot["DiMuon_InvMass"]     = new MultiSamplePlot(datasets, "DiMuon_InvMass", 60, 0, 120, "DiMuon_InvMass");
    MSPlot["NbOfLooseMuon"]     = new MultiSamplePlot(datasets, "NbOfLooseMuon", 10, 0, 10, "Nb. of loose muons");
    
    //B-tagging discriminators
    MSPlot["BdiscBJetCand_CSV"]  = new MultiSamplePlot(datasets, "HighestBdisc_CSV", 75, 0, 1, "CSV b-disc.");
    MSPlot["HighestBdisc_m_ch_CSV"]            = new MultiSamplePlot(datasets, "HighestBdisc_mm_ch_CVS", 100, 0, 1, "CSV b-disc.");
    MSPlot["HighestBdisc_e_ch_CSV"]            = new MultiSamplePlot(datasets, "HighestBdisc_ee_ch_CVS", 100, 0, 1, "CSV b-disc.");
    
    //Jets
    MSPlot["NbOfSelectedJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0, 15, "Nb. of jets");
    MSPlot["NbOfSelectedLightJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
    MSPlot["NbOfSelectedBJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0, 8, "Nb. of jets");
    MSPlot["JetEta"]                  = new MultiSamplePlot(datasets, "JetEta", 30,-3, 3, "Jet #eta");
    MSPlot["JetPhi"]                  = new MultiSamplePlot(datasets, "JetPhi", 50, -4,4 , "Jet #phi");
    
    MSPlot["NbOfBadTrijets"]                  = new MultiSamplePlot(datasets, "NbOfBadTriJets", 150, 0, 150, "Nb. of Bad Combs");
    
    MSPlot["FirstTriJetMass_1BJet_g"] = new MultiSamplePlot(datasets, "FirstTriJetMass_1BJet_g", 15, 0, 1000, "m_{bjj}");
    MSPlot["FirstTriJetMass_1BJet_All"] = new MultiSamplePlot(datasets, "FirstTriJetMass_1BJet_All", 50, 0, 1000, "m_{bjj}");
    MSPlot["FirstTriJetPt_1BJet_g"] = new MultiSamplePlot(datasets, "FirstTriJetPt_1BJet_g", 50, 0, 1000, "Pt_{bjj}");
    MSPlot["FirstDiJetMass"] = new MultiSamplePlot(datasets, "FirstDiJetMass", 50, 0, 200, "m_{jj}");
    MSPlot["FirstDiJetMass_g"] = new MultiSamplePlot(datasets, "FirstDiJetMass_g", 50, 0, 200, "m_{jj}");
    MSPlot["SecondTriJetMass_1BJet"] = new MultiSamplePlot(datasets, "SecondTriJetMass_1BJet", 15, 0, 1000, "m_{bjj}");
    MSPlot["SecondTriJetMass_1BJet_g"] = new MultiSamplePlot(datasets, "SecondTriJetMass_1BJet_g", 10, 0, 1000, "m_{bjj}");
    MSPlot["SecondTriJetMass_1BJet_g_chi2"] = new MultiSamplePlot(datasets, "SecondTriJetMass_1BJet_g_chi2", 30, 0, 800, "chi^{2}");
    MSPlot["SecondTriJetMass_1BJet_All"] = new MultiSamplePlot(datasets, "SecondTriJetMass_1BJet_All", 10, 0, 1000, "m_{bjj}");
    MSPlot["SecondTriJetPt_1BJet"] = new MultiSamplePlot(datasets, "SecondTriJetPt_1BJet", 50, 0, 1000, "Pt_{bjj}");
    MSPlot["SecondDiJetMass"] = new MultiSamplePlot(datasets, "SecondTriJetMass", 20, 0, 200, "m_{jj}");
    MSPlot["SecondDiJetMass_g"] = new MultiSamplePlot(datasets, "SecondTriJetMass_g", 20, 0, 200, "m_{jj}");
    
    MSPlot["TriJetMass_Matched"] = new MultiSamplePlot(datasets, "TriJetMassMatched", 100, 0, 1000, "m_{bjj}");
    MSPlot["TriJetMass_UnMatched"] = new MultiSamplePlot(datasets, "TriJetMassUnMatched", 100, 0, 1000, "m_{bjj}");

    MSPlot["MVA1TriJetMass"] = new MultiSamplePlot(datasets, "MVA1TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA2TriJetMass"] = new MultiSamplePlot(datasets, "MVA2TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA2ndPassTriJetMass"] = new MultiSamplePlot(datasets, "MVA2ndPassTriJetMass", 30, 0, 1000, "m_{bjj}");
    
    MSPlot["MVA1TriJetMassMatched"] = new MultiSamplePlot(datasets, "MVA1TriJetMassMatched", 75, 0, 500, "m_{bjj}");

    
    MSPlot["BDisc_Asym"] = new MultiSamplePlot(datasets, "BDisc_Asym", 100, -5, 5, "BDiscAsym");
    MSPlot["RPer"] = new MultiSamplePlot(datasets, "RPer", 100, 0, 18, "Rper");
    MSPlot["HTComb"] = new MultiSamplePlot(datasets, "HTComb", 100, 0, 800, "HTComb");
    MSPlot["WMassComb"] = new MultiSamplePlot(datasets, "WMassComb", 75, 0, 650, "WMassComb");

    MSPlot["WMt"] = new MultiSamplePlot(datasets, "WMt", 50, 0, 250, "W Transverse Mass");
    MSPlot["LepWMass"] = new MultiSamplePlot(datasets, "LepWMass", 50, 0, 200, "MuMET");
    MSPlot["LepWMass_g"] = new MultiSamplePlot(datasets, "LepWMass_g", 50, 0, 200, "MuMET");
    MSPlot["MuMetBMasses"] = new MultiSamplePlot(datasets, "MuMetBMasses", 50, 0, 600, "m_{muMETb}");
    MSPlot["MuMetBMasses_g"] = new MultiSamplePlot(datasets, "MuMetBMasses_g", 50, 0, 1000, "m_{muMETb}");
    MSPlot["MuMetBPt"] = new MultiSamplePlot(datasets, "MuMetBPt", 50, 0, 1000, "Pt_{muMETb}");
    MSPlot["MuMetBPt_g"] = new MultiSamplePlot(datasets, "MuMetBPt_g", 50, 0, 1000, "Pt_{muMETb}"); 
    MSPlot["MuMetBMasses_chi2"] = new MultiSamplePlot(datasets, "MuMetBMasses_chi2", 50, 0, 1000, "\\chi^{2}");
    MSPlot["SelectedJetPt"] = new MultiSamplePlot(datasets, "JetPt", 50, 0, 100, "PT_{jet}");
    MSPlot["4thJetPt"] = new MultiSamplePlot(datasets, "4thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["5thJetPt"] = new MultiSamplePlot(datasets, "5thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["6thJetPt"] = new MultiSamplePlot(datasets, "6thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["7thJetPt"] = new MultiSamplePlot(datasets, "7thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["SelectedJetPt_light"] = new MultiSamplePlot(datasets, "JetPt_light", 50, 0, 1000, "PT_{lightjet}");
    MSPlot["SelectedJetPt_b"] = new MultiSamplePlot(datasets, "JetPt_b", 50, 0, 1000, "PT_{bjet}");
    MSPlot["HT_SelectedJets"] = new MultiSamplePlot(datasets, "HT_SelectedJets", 50, 0, 1500, "HT");
    
    MSPlot["H"] = new MultiSamplePlot(datasets, "H", 50, 0, 3000, "H");
    MSPlot["HX"] = new MultiSamplePlot(datasets, "HX", 50, 0, 1500, "HX");
    
    MSPlot["HTH"] = new MultiSamplePlot(datasets, "HT/H", 50, 0, 1, "HT/H");
    MSPlot["HTXHX"] = new MultiSamplePlot(datasets, "HTX/HX", 50, 0, 1, "HTX/HX");

    MSPlot["MHT_SelectedJets"] = new MultiSamplePlot(datasets, "MHT_SelectedJets", 75, 0, 1000, "MHT");
    MSPlot["MHTSig_SelectedJets"] = new MultiSamplePlot(datasets, "MHTSig_SelectedJets", 75, 0, 30, "MHTSig");
    MSPlot["MET"] = new MultiSamplePlot(datasets, "MET", 75, 0, 700, "MET");
    MSPlot["MET_MHT"]= new MultiSamplePlot(datasets, "MET_MHT", 75, 0, 200, "MET_MHT");
    MSPlot["STLep"] = new MultiSamplePlot(datasets, "STLep", 40, 0, 2500, "STLep");
    MSPlot["STJet"] = new MultiSamplePlot(datasets, "STJet", 40, 0, 2500, "STJet");
    
    MSPlot["EventMass"] = new MultiSamplePlot(datasets, "EventMass", 40, 0, 3000, "EventMass");
    MSPlot["EventMassX"] = new MultiSamplePlot(datasets, "EventMassX", 30, 0, 2000, "EventMassX");
    MSPlot["SumJetMass"] = new MultiSamplePlot(datasets, "SumJetMass", 40, 0, 3000, "SumJetMass");
    MSPlot["SumJetMassX"] = new MultiSamplePlot(datasets, "SumJetMassX", 30, 0, 2000, "SumJetMassX");
    MSPlot["HTX"] = new MultiSamplePlot(datasets, "HTX", 20, 0, 1000, "HTX");
    
    MSPlot["PTBalTopEventX"] = new MultiSamplePlot(datasets, "PTBal_TopX", 35, 0, 500, "PTBal_TopX");
    MSPlot["PTBalTopSumJetX"] = new MultiSamplePlot(datasets, "PTBal_TopSumJetX", 35, 0, 500, "PTBal_TopSumJetX");
    
    MSPlot["PTBalTopMuMet"] = new MultiSamplePlot(datasets, "PTBal_TopMuMet", 35, 0, 600, "PTBal_TopMuMet");
    MSPlot["PTBalTopMuMetB"] = new MultiSamplePlot(datasets, "PTBal_TopMuMetB", 35, 0, 600, "PTBal_TopMuMetB");

    
    MSPlot["deltaMTopMuMet"] = new MultiSamplePlot(datasets, "deltaMTopMuMet", 75, -500, 500, "deltaMTopMuMet");

    MSPlot["deltaMTopMuMetB"] = new MultiSamplePlot(datasets, "deltaMTopMuMetB", 75, -500, 500, "deltaMTopMuMetB");

    
    
    MSPlot["HT_CombinationRatio"] = new MultiSamplePlot(datasets, "HT_CombinationRatio", 50, 0, 1, "HT_Ratio");
    
    MSPlot["MVA"] = new MultiSamplePlot(datasets, "MVA", 50, -1, 1, "BDT Disciminator");

    
    //Electrons
    MSPlot["ElectronPt"]              = new MultiSamplePlot(datasets, "ElectronPt", 50, 0, 100, "PT_{e}");
    MSPlot["NbOfLooseElectron"] = new MultiSamplePlot(datasets, "NbOfLooseElectron", 10, 0, 10, "Nb. of loose electrons");

    //Declare arrays of MSPlots
    Int_t minNJets=6, maxNJets=9, minNBJets=1, maxNBJets=4;

  for (Int_t q = minNJets; q <= maxNJets; q++){
    for (Int_t p = minNBJets; p<= maxNBJets; p++){

    string NJets_str = static_cast<ostringstream*>( &(ostringstream() << q) )->str();
    string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << p) )->str();
    string HT_Name = "HT_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string HTX_Name = "HTX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string H_Name = "H_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string HX_Name = "HX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string HTH_Name = "HTH_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string HTXHX_Name = "HTXHX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;


    string EventMass_Name = "EventMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string EventMassX_Name = "EventMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string SumJetMass_Name = "SumJetMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string SumJetMassX_Name = "SumJetMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string PTBalTopEventX_Name = "PTBalTopEventX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string PTBalTopSumJetX_Name = "PTBalTopSumJetX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
        
    string FirstTriJetMass_1BJet_g_Name = "FirstTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string FirstDiJetMass_1BJet_g_Name = "FirstDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string SecondDiJetMass_1BJet_g_Name = "SecondDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

    string SecondTriJetMass_1BJet_g_Name = "SecondTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string SecondTriJetMass_1BJet_g_chi2_Name = "SecondTriJetMass_1BJet_g_chi2_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MuMetBMasses_g_Name = "MuMetBMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MuMetMasses_g_Name = "MuMetMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MET_Name = "MET"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MVA1TriJetMass_Name = "MVA1TriJetMass"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MVA2TriJetMass_Name = "MVA2TriJetMass"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
        
    MSPlot[MET_Name.c_str() ] = new MultiSamplePlot(datasets, MET_Name.c_str() , 50, 0, 700, "MET");
    MSPlot[MuMetMasses_g_Name.c_str() ] = new MultiSamplePlot(datasets, MuMetMasses_g_Name.c_str() , 50, 0, 700, "M_{muMET}");
    MSPlot[MuMetBMasses_g_Name.c_str() ] = new MultiSamplePlot(datasets, MuMetBMasses_g_Name.c_str() , 50, 0, 700, "M_{muMETb}");
    MSPlot[HT_Name.c_str() ] = new MultiSamplePlot(datasets, HT_Name.c_str() , 20, 0, 1700, "HT");
    MSPlot[HTX_Name.c_str() ] = new MultiSamplePlot(datasets, HTX_Name.c_str() , 20, 0, 1200, "HTX");
    MSPlot[H_Name.c_str() ] = new MultiSamplePlot(datasets, H_Name.c_str() , 20, 0, 1700, "H");
    MSPlot[HX_Name.c_str() ] = new MultiSamplePlot(datasets, HX_Name.c_str() , 20, 0, 1700, "HX");
    MSPlot[HTH_Name.c_str() ] = new MultiSamplePlot(datasets, HTH_Name.c_str() , 20, 0, 1, "HT/H");
    MSPlot[HTXHX_Name.c_str() ] = new MultiSamplePlot(datasets, HTXHX_Name.c_str() , 20, 0, 1, "HTX/HX");
        
    MSPlot[EventMass_Name.c_str() ] = new MultiSamplePlot(datasets, EventMass_Name.c_str() , 15, 0, 2500, "EventMass");
    MSPlot[EventMassX_Name.c_str() ] = new MultiSamplePlot(datasets, EventMassX_Name.c_str() , 15, 0, 1700, "EventMassX");
    
    MSPlot[SumJetMass_Name.c_str() ] = new MultiSamplePlot(datasets, SumJetMass_Name.c_str() , 15, 0, 2500, "SumJetMass");
    MSPlot[SumJetMassX_Name.c_str() ] = new MultiSamplePlot(datasets, SumJetMassX_Name.c_str() , 15, 0, 1700, "SumJetMassX");
    MSPlot[PTBalTopEventX_Name.c_str() ] = new MultiSamplePlot(datasets, PTBalTopEventX_Name.c_str() , 15, 0, 500, "PTBalTopEvent_Name");
    MSPlot[PTBalTopSumJetX_Name.c_str() ] = new MultiSamplePlot(datasets, PTBalTopSumJetX_Name.c_str() , 15, 0, 500, "PTBalTopSumJetX_Name");
        
    MSPlot[FirstTriJetMass_1BJet_g_Name.c_str() ] = new MultiSamplePlot(datasets, FirstTriJetMass_1BJet_g_Name.c_str() , 60, 0 ,350 , "M_{bjj}");
    MSPlot[FirstDiJetMass_1BJet_g_Name.c_str() ] = new MultiSamplePlot(datasets, FirstDiJetMass_1BJet_g_Name.c_str() , 60, 0 ,350 , "M_{jj}");
    MSPlot[SecondDiJetMass_1BJet_g_Name.c_str() ] = new MultiSamplePlot(datasets, SecondDiJetMass_1BJet_g_Name.c_str() , 60, 0 ,350 , "M_{jj}");
    MSPlot[SecondTriJetMass_1BJet_g_Name.c_str() ] = new MultiSamplePlot(datasets, SecondTriJetMass_1BJet_g_Name.c_str() , 60, 0 ,350 , "M_{bjj}");
    MSPlot[SecondTriJetMass_1BJet_g_chi2_Name.c_str() ] = new MultiSamplePlot(datasets, SecondTriJetMass_1BJet_g_Name.c_str() , 25, 0 ,250 , "#chi^{2}");
    MSPlot[MVA1TriJetMass_Name.c_str() ] = new MultiSamplePlot(datasets, MVA1TriJetMass_Name.c_str() , 40, 20 ,500 , "M_{bjj}");
    MSPlot[MVA2TriJetMass_Name.c_str() ] = new MultiSamplePlot(datasets, MVA2TriJetMass_Name.c_str() , 40, 20 ,500 , "M_{bjj}");
}
}

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////

    histo1D["RelIso_PreTrig"] = new TH1F("RelIso_PreTrig","RelIso_PreTrig;RelIso_PreTrig;#events",40,0,0.75);
    histo1D["RelIso_PostTrig"] = new TH1F("RelIso_PostTrig","RelIso_PostTrig;RelIso_PostTrig;#events",40,0,0.75);

    histo1D["Pt_PreTrig"] = new TH1F("Pt_PreTrig","Pt_PreTrig;Pt_PreTrig;#events",25,0,500);
    histo1D["Pt_PostTrig"] = new TH1F("Pt_PostTrig","Pt_PostTrig;Pt_PostTrig;#events",25,0,500);
    
    histo1D["d0_PreTrig"] = new TH1F("d0_PreTrig","d0_PreTrig;d0_PreTrig;#events",20,0,0.02);
    histo1D["d0_PostTrig"] = new TH1F("d0_PostTrig","d0_PostTrig;d0_PostTrig;#events",20,0,0.02);
    
    histo1D["Eta_PreTrig"] = new TH1F("Eta_PreTrig","Eta_PreTrig;Eta_PreTrig;#events",20,-2.5,2.5);
    histo1D["Eta_PostTrig"] = new TH1F("Eta_PostTrig","Eta_PostTrig;Eta_PostTrig;#events",20,-2.5,2.5);
    
    histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
    
    
  for (unsigned int d = 0; d < datasets.size(); d++){
    histo2D[("RelIso_vs_MET_"+datasets[d]->Name()).c_str()] = new TH2F(("RelIso_vs_MET_"+datasets[d]->Name()).c_str(),"RelIso:MET",100,0,1000, 100, 0,1);
    histo2D[("ThirdTopMass_vs_HT_"+datasets[d]->Name()).c_str()] = new TH2F(("ThirdTopMass_vs_HT_"+datasets[d]->Name()).c_str(),"ThirdTopMass_vs_HT",30,0,1000, 30, 0,1500);
    histo2D[("ThirdTopMass_vs_SecondTopMass"+datasets[d]->Name()).c_str()] = new TH2F(("ThirdTopMass_vs_SecondTopMass"+datasets[d]->Name()).c_str(),"ThirdTopMass_vs_SecondTopMass",30,0,1000, 30, 0,1000);
      
    histo2D[("MassChi2_vs_HT"+datasets[d]->Name()).c_str()] = new TH2F(("MassChi2_vs_HT"+datasets[d]->Name()).c_str(),"MassChi2_vs_HT",15,0,400, 15, 0,1400);
    histo2D[("EventMassX_vs_HT"+datasets[d]->Name()).c_str()] = new TH2F(("EventMassX_vs_HT"+datasets[d]->Name()).c_str(),"EventMassX_vs_HT",20,0,1500, 15, 0,1200);
    histo2D[("EventMassX_vs_HTX"+datasets[d]->Name()).c_str()] = new TH2F(("EventMassX_vs_HTX"+datasets[d]->Name()).c_str(),"EventMassX_vs_HTX",20,0,1500, 15, 0,1200);

      histo2D[("OfflineMuonPt_vs_TriggerMuonPt"+datasets[d]->Name()).c_str()] = new TH2F(("OfflineMuonPt_vs_TriggerMuonPt"+datasets[d]->Name()).c_str(),"OfflineMuonPt_vs_TriggerMuonPt",75,0,300, 75, 0,300);
      
     histo2D[("OfflineMuonEta_vs_TriggerMuonEta"+datasets[d]->Name()).c_str()] = new TH2F(("OfflineMuonEta_vs_TriggerMuonEta"+datasets[d]->Name()).c_str(),"OfflineMuonEta_vs_TriggerMuonEta",15,-2.5,2.5, 15, -2.5,2.5);

  }
    
    
  //  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

    string pathPNG_MVA = "MVAPlots_"+postfix+channelpostfix;

    
  string pathPNG = "FourTop"+postfix+channelpostfix;
  pathPNG += "_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////// Selection Table/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

    ////////////////////////////////////////////////////////////////////
    ///////////////////// Tree for MVA    ////////////////////////
    ////////////////////////////////////////////////////////////////////
    
   // string TreeFileName = "FourTopTree.root";
    
   // TFile* treeFile;
    TTree* myFourTopTree;
   // FourTopTree* myBranch_selectedEvents = 0;	

    //cout << "INFO: creating FourTopTree file "+TreeFileName << endl;
    
    bool datadriven = false;
    
    if(!datadriven){
     //treeFile = new TFile(TreeFileName.c_str(),"RECREATE");
      //  myFourTopTree = new TTree("FourTopTree;","Tree containing the Four top information");
      //  myFourTopTree->Branch("FourTopBranch_selectedEvents","FourTopTree",&myBranch_selectedEvents);
    }
    
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// PileUp Reweighting - N True interactions, recommended method for 2012 //////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  //NEW METHOD (TRUE INTERACTIONS)
    LumiReWeighting LumiWeights;
 
    LumiWeights = LumiReWeighting("PUReweighting/pileup_MC_S10.root", "PUReweighting/PURewtoptree_id_2014_1350565897_PileupHistogram.root", "pileup", "pileup");
    
//  reweight::PoissonMeanShifter PShiftDown_ = reweight::PoissonMeanShifter(-0.6);
  //reweight::PoissonMeanShifter PShiftUp_ = reweight::PoissonMeanShifter(0.6);

  cout << " Initialized LumiReWeighting stuff" << endl;
  
  //exit(1);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Loop on datasets
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
      int event_start;
      if (dataSetName == "Data")TrainMVA=false, trainEventMVA=false;
      
      if (TrainMVA){
        event_start = 0;
        end = 200000;
      }
 //     else if (trainEventMVA){
          if(  (dataSetName == "NP_overlay_TTTT" ) ){
        event_start = 0;
     //   end = 9000;
        end = datasets[d]->NofEvtsToRunOver();

          }
          else if (dataSetName == "Data"){
              event_start = 0;
              
             // end = 100;
              end = datasets[d]->NofEvtsToRunOver();
          }
          
          
          else{
              event_start = 0;
          //    end = 100000;
            end = datasets[d]->NofEvtsToRunOver();
          }

     // }
      
      
      cout <<"Starting event = = = = "<< event_start  << endl;

    if (verbose > 1) cout << " - Loop over events " << endl;
      
      int nBBBar, nCCBar, nLLBar;
      nBBBar=  nCCBar = nLLBar = 0;
      
      
      double MHT = 0.,MHTSig = 0., STJet = 0., EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0.  ,H = 0., HX =0., HT = 0., HTX = 0.,HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0., PTBalTopEventX = 0., PTBalTopSumJetX =0., PTBalTopMuMet=0;
      

      
      double end_d = end;
      double endfrac = end_d/1.;
      
    for (unsigned int ievt = event_start; ievt < end ; ievt++)
    {
        
MHT = 0.,MHTSig = 0., STJet = 0., EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0.  ,H = 0., HX =0., HT = 0., HTX = 0.,HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0., PTBalTopEventX = 0., PTBalTopSumJetX =0.;
        
	if(ievt%1000 == 0)
		std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

     selecTableMu.Fill(d,0,1.);
  //   selecTableEl.Fill(d,0,1.);

	//load event
	event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

  //  cout <<"test rho :  " <<event->kt6PFJets_rho() <<endl;   
        float rho = event->kt6PFJets_rho();
        
    //load mcparticles to check jet flavour for ttjets events
        
    vector<TRootMCParticle*> mcParticles_flav;
    Int_t ttbar_flav = -1;

        
        double nExB,nExC,nExL;
        nExB = nExC = nExL = 0.;
        
    TRootGenEvent* genEvt_flav = 0;
    genEvt_flav = treeLoader.LoadGenEvent(ievt,false);
    treeLoader.LoadMCEvent(ievt, genEvt_flav, 0, mcParticles_flav,false); 
        
       // cout <<" mc parts "<< mcParticles_flav.size()  <<endl;

        if(  (dataSetName == "TTJets_ll" || dataSetName == "TTJets_cc" || dataSetName == "TTJets_bb" ) )
        {
        for(unsigned int p=0; p<mcParticles_flav.size(); p++) {
          
            if(mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())==5 && abs(mcParticles_flav[p]->motherType())!=6) {
               // ttbar_flav=2;
                nExB++;  
            }
            
            else if (mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())==4 && abs(mcParticles_flav[p]->motherType())!=6
                && abs(mcParticles_flav[p]->motherType())!=5 && abs(mcParticles_flav[p]->motherType())!=24
                ){
           // ttbar_flav=1;
                 nExC++; 
            }
            
            else if (mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())<4 && abs(mcParticles_flav[p]->motherType())!=6){
                // ttbar_flav=1;
                nExL++; 
            }

        }
            
     //   cout <<"TTBar flav composition : " << nExL  <<"  light, " << nExC <<"  C, " << nExB<< "  B" <<  endl;
            
     //   if (ttbar_flav != 1 && ttbar_flav != 2 ) ttbar_flav = 0;
       
            if (nExB >= 2.){
            ttbar_flav =2; 
            nBBBar++ ; //  bbbar
            }
            else if ( nExC >=2.) {
            ttbar_flav =1; 
            nCCBar++ ; //ccbar
            }
            else{
            ttbar_flav =0.; 
                nLLBar++;  //llbar   
            }
            
            if (ttbar_flav ==0 && (dataSetName == "TTJets_cc"  || dataSetName == "TTJets_bb"))  continue;
            if (ttbar_flav ==1 && (dataSetName == "TTJets_ll"  || dataSetName == "TTJets_bb" ))  continue;
            if (ttbar_flav ==2 && (dataSetName == "TTJets_ll"  || dataSetName == "TTJets_cc" ))  continue;
        
        }
        
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
    std::string filterName = "";
	int currentRun = event->runId();
	if(previousRun != currentRun)
	{
	 // cout <<"What run? "<< currentRun<<endl;
      		previousRun = currentRun;
		if(Muon)
		{
			//if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
		//	{
                
                
                // semi-muon
                if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
                {
                    if( event->runId() <= 190738 ){
                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";
                    }
                    else if( event->runId() >= 191043 && event->runId() <= 193621 ){
                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v12"), currentRun, iFile);
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";
                    }
                    else if( event->runId() >= 193834 && event->runId() <= 196531 ){
                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                    }
                    else if( event->runId() >= 198049 && event->runId() <= 199608){
                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v14"), currentRun, iFile);
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                    }
                    else if( event->runId() >= 199698 && event->runId() <= 208357){
                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v15"), currentRun, iFile);
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                    }
                    else{
                        cout << "Unknown run for SemiMu HLTpath selection: " << event->runId() << endl;
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";

                    }
                    if( itrigger == 9999 )
                    {
                        cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
                        exit(-1);
                    }
                }
                
                
                
		      /*
			  if(event->runId() >= 190456 && event->runId() <= 190761)
				 itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
			  else if ( event->runId() >= 190762 && event->runId() <= 191511 )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v3"), currentRun, iFile);
			  else if ( event->runId() >= 191512  && event->runId() <= 193805  )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v4"), currentRun, iFile);
			  //for  top triggers change PD at this point.
			  else if ( event->runId() >= 193806  && event->runId() <= 194269  )
			         itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v1"), currentRun, iFile);
                //extendidng the 52X 1 fb
			  else if ( event->runId() >= 194269 && event->runId() <= 196531  )
			    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);
                
              else if ( event->runId() >= 196532 && event->runId() <= 198522  )
                  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v2"), currentRun, iFile);

              else if ( event->runId() >= 198523 && event->runId() <= 199608  )
                  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v2"), currentRun, iFile);
                
              else if ( event->runId() >= 199609 && event->runId() <= 201678  )
                  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25_v1"), currentRun, iFile);
                
              else if ( event->runId() >= 201679 && event->runId() <= 202016  )
                  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25_v1"), currentRun, iFile);
			
  		  if(itrigger == 9999)
				{
    		  cout << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
		    exit(1);
 	 			}
	   		}
               */
	   		else 
	   		{
			  //for refsel, HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2 recommended not HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2
			 // if(dataSetName == "TTJets" || dataSetName == "TTJets_ll" || dataSetName == "TTJets_cc"  || dataSetName == "TTJets_bb" ) itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
                
                  if(dataSetName == "TTJets" || dataSetName == "TTJets_ll" || dataSetName == "TTJets_cc"  || dataSetName == "TTJets_bb" ) itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);
                

			 // else  if(dataSetName == "TTTT") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v2"),currentRun, iFile);
                //else  if(dataSetName == "TTTT") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_eta2p1_TriCentralPFNoPUJet30_v2"),currentRun, iFile);
                else  if(dataSetName == "NP_overlay_TTTT") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"),currentRun, iFile);
                
                else  if(dataSetName == "T1TTTT") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2"),currentRun, iFile);

			  // 	    if(dataSetName == "TTJets") itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v8"), currentRun, iFile);
      			else if (dataSetName == "WJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);

                else if (dataSetName == "ZJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);

			    else if (dataSetName == "SingleTop_t") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);

       			else if (dataSetName == "SingleTop_tW_T") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);
       			else if (dataSetName == "SingleTop_tW_TBar") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile);

				else if (dataSetName == "MultiJet") {
                                                                   
	                 	  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v1"), currentRun, iFile); 
				
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
//	if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
//		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
//	else
//		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// Type I MET corrections: (Only for |eta| <=4.7 ) //////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");      
	//if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
	//	jetTools->correctMETTypeOne(init_jets,mets[0],true);
	//else
	//	jetTools->correctMETTypeOne(init_jets,mets[0],false);
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
//////////////////////////////////////////////// Beam scrapping veto and PU reweighting ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
			
////////////////////////////////////////////////////////////////////////////////////////////// Event selection///////////////////////////////////////////////////////////////////////////////////////////////////////

	//	selecTableMu.Fill(d,0, 1.);//datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
	//	selecTableEl.Fill(d,1,scaleFactor);
	       
	// Apply trigger selection
	trigged = treeLoader.EventTrigged (itrigger);
	if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;


	//Applying trigger selection again with 2012 Muon+Jets triggers.
//	if(!trigged)		   continue;

	// Declare selection instance    
	Selection selection(init_jets, init_muons, init_electrons, mets);

	// Define object selection cuts
	selection.setJetCuts(40.,2.5,0.01,1.,0.98,0,0);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets

        
        selection.setElectronCuts();
//	selection.setElectronCuts(30.,2.5,0.1,0.02,0.,999.,0,1); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets, MaxMissingHits
	selection.setLooseElectronCuts(20,2.5,0.2,0.);

        //for muon + jets channel
	//selection.setMuonCuts(26.,2.1,.12,0,0.2,0,1,0.5,5,1 ); //Pt,Eta,RelIso,NValidMuHits,d0, dRJets, NMatchedStations,DistVzPVz,NTrackerLayersWithMeas 
        
        //for dimuon + jets channel
   // selection.setMuonCuts(20.,2.4,.2,0,0.2,0,1,0.5,5,1 ); //Pt,Eta,RelIso,NValidMuHits,d0, dRJets, NMatchedStations,DistVzPVz,
	//selection.setMuonCuts();
        
    selection.setMuonCuts(26,2.1,0.12,0.2,0.3,1,0.5,5,0); // refSel 2012 values: Pt, Eta, RelIso, d0, dRJets, NMatchedStations, Dz, NTrackerLayerswithMeas, nValidPixelHits

        
    selection.setLooseMuonCuts(10,2.5,0.2);
	  
	//Select objects 
	//	vector<TRootElectron*> selectedElectrons_NoIso = selection.GetSelectedElectrons(20,2.4,999.,vertex[0]);
	//vector<TRootElectron*> selectedElectrons        = selection.GetSelectedElectrons(vertex[0]);
//    vector<TRootElectron*> selectedElectrons        = selection.GetSelectedElectrons(vertex[0], rho);
   
        vector<TRootElectron*> selectedElectrons        = selection.GetSelectedElectrons();

        
	vector<TRootElectron*> selectedExtraElectrons;
	vector<TRootMuon*>     selectedMuons_NoIso      = selection.GetSelectedMuons(26,2.4,999.); 
	vector<TRootMuon*>     selectedMuons            = selection.GetSelectedMuons();
	vector<TRootMuon*>     selectedExtraMuons;
	vector<TRootJet*>      selectedJets             = selection.GetSelectedJets(true); // ApplyJetId
    vector<TRootJet*>      selectedJets2ndPass;
    vector<TRootJet*>      selectedSoftJets         = selection.GetSelectedJets(20.,2.5, selectedMuons, 0., true); // ApplyJetId
	vector<TRootMuon*>     selectedLooseMuons       = selection.GetSelectedLooseMuons();
    vector<TRootElectron*> selectedLooseElectrons   = selection.GetSelectedElectrons(); // VBTF ID
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
////////////////////////////////////////////////////////// Sync'ing cutflow//////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
                
	selecTableMu.Fill(d,1,1.) ;
	//selecTableEl.Fill(d,1,1.) ;
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
/*
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
        */
        
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Applying baseline offline event selection here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (debug)	cout <<" applying baseline event selection..."<<endl;
        
    if (!(selectedJets.size() >= 6)) continue;
        
        int nbadcombs = Factorial(selectedJets.size())    /  (  (Factorial(3)) * (Factorial(  selectedJets.size() - 3     ) )       ) ;
        
    //Apply the staggered Jet Pt cut (lep + jets)
//	if  ( selectedJets[0]->Pt() >45. &&  selectedJets[1]->Pt() >45. && selectedJets[2]->Pt() >45. && selectedJets[3]->Pt() >45. && selectedJets[4]->Pt() >45. ) JetCut = 1;
        
        
        //Apply the staggered Jet Pt cut (dilep + jets)
        //if  ( selectedJets[0]->Pt() >45. &&  selectedJets[1]->Pt() >45. && selectedJets[2]->Pt() >45. && selectedJets[3]->Pt() >45. && selectedJets[4]->Pt() >45. && selectedJets[5]->Pt() >45. ) JetCut = 1;
        JetCut=1;
        
	if (debug) cout<<" jet1 pt =  "<<selectedJets[0]->Pt() << "   "<< " jet2 pt =  "<<selectedJets[1]->Pt() << "   "<< " jet2 pt =  "<<selectedJets[2]->Pt() << "   "<< " jet3 pt =  "<<selectedJets[3]->Pt() << "  JetCut?"  << JetCut  <<endl;
 
	if (debug)	cout <<" filling bjet vec "<<endl;
	//filling vector of b-jets
	for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ ){
	  if( selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue) selectedBJets.push_back(selectedJets[seljet]);
	  else selectedLightJets.push_back(selectedJets[seljet]);
	}
                
        
        double temp_HT = 0;
        double temp_jetpt = 0;
        for (Int_t seljet0 =0; seljet0 < selectedJets.size(); seljet0++ ){
          
            //Event-level variables
            temp_jetpt = selectedJets[seljet0]->Pt();
            temp_HT = temp_HT + temp_jetpt;
        }
        
        
	if(Muon){

    if (debug) cout <<"Number of Muons, Jets, BJets, JetCut  ===>  "<< selectedMuons.size() <<"  "  << selectedJets.size()   <<"  " <<  selectedBJets.size()   <<"  "<<JetCut  <<endl;


	  //Apply the selection
	  if  (  !( JetCut==1 && selectedBJets.size() >= 2 && temp_HT > 400. &&  nMu == 1 )) continue;
        double trigmuonPT=500.;
        double trigmuonEta=500.;
        double trigmuonPhi=500.;
        bool foundfilter = false;

        double offlinemuonPT=500;;
        int matchmu_index;
       // if(dataSetName!="Data"){
     //   double bTagSF = (  bTool->getWeight(selectedBJets[0]->Pt(), selectedBJets[0]->Eta(), 5,0 )*  bTool->getWeight(selectedBJets[0]->Pt(), selectedBJets[0]->Eta(), 5,0 ));

      //  double bTagSF_b =   bTool->getWeight(selectedBJets[0]->Pt(), selectedBJets[0]->Eta(), 5,1 );
      //  double bTagSF_c =   bTool->getWeight(selectedBJets[0]->Pt(), selectedBJets[0]->Eta(), 4,1 );
       // double bTagSF_l =   bTool->getWeight(selectedBJets[0]->Pt(), selectedBJets[0]->Eta(), 1,1 );
        
       // cout <<" SF = >"<< bTagSF_b <<"  "<< bTagSF_c<< "  "<< bTagSF_l <<endl;
       // scaleFactor = scaleFactor*bTagSF;
        //}
        
        std::map<std::string,std::vector<TopTree::triggeredObject> > filters= event->getTriggerFilters(); //previously: double instead of TopTree::triggeredObject

        double delr = 999.;
       // std::map<std::string,std::vector<TopTree::triggeredObject> >::iterator itr = filters.begin();
        
        for(std::map<std::string, std::vector<TopTree::triggeredObject> > ::iterator itr = filters.begin(), itr_end = filters.end(); itr != itr_end; ++itr) {
            //itr->second.M;
            string sub_string = itr->first.substr(0,6);
            string sub_filterName = filterName.substr(0,6);
            if(sub_string == sub_filterName ){
                foundfilter = true;
    
      if (itr->second.size() ==2 ) cout <<"n trigger objects "<< itr->second[0].pt <<"   "<< itr->second[1].pt  <<" "<< itr->second[0].id   << "  " << itr->second[1].id<<endl;
                
            //cout<<"filter name "<< itr->first <<"  ID = " << itr->second[0].id  << "  object pt = "<< itr->second[0].pt  <<endl;
            trigmuonPT = itr->second[0].pt;
            trigmuonEta = itr->second[0].eta;
            trigmuonPhi = itr->second[0].phi;
                
                
                double delr_temp =   sqrt(   pow(selectedMuons[0]->Eta() -  trigmuonEta ,2) + pow(selectedMuons[0]->Phi() - trigmuonPhi  ,2 )   ) ;
                if (delr_temp<delr) {
                    delr = delr_temp;
                }
                
            }
        }
        
        
        float reliso = (selectedMuons[0]->chargedHadronIso()+selectedMuons[0]->neutralHadronIso()+selectedMuons[0]->photonIso())/selectedMuons[0]->Pt();
        if (foundfilter){
        histo1D["RelIso_PreTrig"]->Fill(reliso);
        histo1D["Pt_PreTrig"]->Fill(selectedMuons[0]->Pt());
        histo1D["d0_PreTrig"]->Fill(selectedMuons[0]->d0());
        histo1D["Eta_PreTrig"]->Fill(selectedMuons[0]->Eta());
        }
        
        if ((foundfilter) && (delr < 0.25) &&(trigged) ){
       // cout <<"matched offline muon to hlt muon..."<<endl;
        histo1D["RelIso_PostTrig"]->Fill(reliso);
        histo1D["Pt_PostTrig"]->Fill(selectedMuons[0]->Pt());
        histo1D["d0_PostTrig"]->Fill(selectedMuons[0]->d0());
        histo1D["Eta_PostTrig"]->Fill(selectedMuons[0]->Eta());
        }
        
        offlinemuonPT=selectedMuons[0]->Pt();
        

        //best offline muon which matches best to trigger muon

        
        histo2D[("OfflineMuonPt_vs_TriggerMuonPt"+datasets[d]->Name()).c_str()]->Fill( trigmuonPT,selectedMuons[matchmu_index]->Pt()  );        
        histo2D[("OfflineMuonEta_vs_TriggerMuonEta"+datasets[d]->Name()).c_str()]->Fill( trigmuonEta,selectedMuons[matchmu_index]->Eta()  );

        MSPlot["TriggerMuonPt"]->Fill(trigmuonPT , datasets[d], true, Luminosity*scaleFactor);
        
        
//From ConfDB:   path HLT_IsoMu24_eta2p1_v15 = HLTBeginSequence + hltL1sMu16Eta2p1 + hltPreIsoMu24eta2p1 + hltL1fL1sMu16Eta2p1L1Filtered0 + HLTL2muonrecoSequence + hltL2fL1sMu16Eta2p1L1f0L2Filtered16Q + HLTL3muonrecoSequence + hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered24Q + HLTL3muoncaloisorecoSequenceNoBools + HLTL3muonisorecoSequence + hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15 + HLTEndSequence
        
        
        //for (std::map<std::string,std::vector<TopTree::triggeredObject>::iterator it=filters.begin(); it != filters.end(); ++it){
       //     std::cout <<"trig filter  = = >   "<< it->first << " => " << it->second << '\n';
        
       // return 0;
   // }
      //  cout <<<" trig filter  "<< filters.first.c_str()    << endl;
        
        
	  if(debug) cout<< "Event passed..."<<endl;
        
        vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
        vector<TRootMCParticle*> mcParticlesMatching_;
        bool muPlusFromTop = false, muMinusFromTop = false;
        bool elPlusFromTop = false, elMinusFromTop = false;
        pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
        leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);


        
////////////////////////////////////////////////////////////////////////////////////
//// Getting Gen Event
////////////////////////////////////////////////////////////////////////////////////
         if(dataSetName != "data" && dataSetName != "Data" && dataSetName != "Data"){
        vector<TRootMCParticle*> mcParticles;
        vector<TRootMCParticle*> mcTops;
             
        mcParticlesMatching_.clear();
        mcParticlesTLV.clear();
        selectedJetsTLV.clear();
                     
             
        int leptonPDG, muonPDG = 13, electronPDG = 11;
        leptonPDG = muonPDG;
             
             
        TRootGenEvent* genEvt = 0;
        genEvt = treeLoader.LoadGenEvent(ievt,false);
        sort(selectedJets.begin(),selectedJets.end(),HighestPt()); 
        treeLoader.LoadMCEvent(ievt, genEvt, 0, mcParticles,false);  
    
        if (debug) cout <<"size   "<< mcParticles.size()<<endl;

         
         
             //Pick out MCParticles of interest
             for (Int_t i=0; i<mcParticles.size(); i++){
                 if( mcParticles[i]->status() != 3) continue;
                 
                 if( (abs(mcParticles[i]->type())< 6||mcParticles[i]->type()==21)&&mcParticles[i]->status() ==3){
                     mcParticlesTLV.push_back(*mcParticles[i]);
                     mcParticlesMatching_.push_back(mcParticles[i]);
                 }
                 
                 // check if there is a mu +/- or a el +/-
                 if( mcParticles[i]->type() == leptonPDG && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 )
                 {
                     //if(muMinusFromTop) cerr<<"muMinusFromTop was already true"<<endl;
                     if(leptonPDG==muonPDG) muMinusFromTop = true;
                     else if(leptonPDG==electronPDG) elMinusFromTop = true;
                 }
                 if( mcParticles[i]->type() == -leptonPDG && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 )
                 {
                     //if(muPlusFromTop) cerr<<"muPlusFromTop was already true"<<endl;
                     if(leptonPDG==muonPDG) muPlusFromTop = true;
                     else if(leptonPDG==electronPDG) elPlusFromTop = true;
                 }
                 
             }
         
         
         
         }
        
        if(dataSetName == "Data"){

      
        MSPlot["MuonPt_PreTrig"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["Muond0_PreTrig"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
        
        MSPlot["MuonRelIsolation_PreTrig"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
    
        if (trigged) {
            
            MSPlot["MuonRelIsolation_PostTrig"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["MuonPt_PostTrig"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["Muond0_PostTrig"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
            
            
         //   histo1D["RelIso_PostTrig"]->Fill(reliso);
           // histo1D["Pt_PostTrig"]->Fill(selectedMuons[0]->Pt());
           // histo1D["d0_PostTrig"]->Fill(selectedMuons[0]->d0());
        }
        }
    
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

        if (debug) cout <<"  nMUons  "<< selectedMuons.size()<< "  n loose muons " << selectedLooseMuons.size()    <<endl;
    //loop to get inv. masses of muon-loosemuon combinations
    for (Int_t mu1 = 0; mu1 < selectedMuons.size(); mu1++){	 
	  for (Int_t mu2 = 0; mu2 < selectedMuons.size(); mu2++){
          
          if (debug) cout <<"making dimuosss   "<< endl;
          
     	    TRootMuon DiMu (TLorentzVector ( selectedMuons[mu1]->Px() + selectedMuons[mu2]->Px() ,selectedMuons[mu1]->Py() + selectedMuons[mu2]->Py()    ,selectedMuons[mu1]->Pz() + selectedMuons[mu2]->Pz() , selectedMuons[mu1]->E() + selectedMuons[mu2]->E() ));
            MSPlot["DiMuon_InvMass"]->Fill( DiMu.M(), datasets[d], true, Luminosity*scaleFactor );  
	  }
	  }
        
        
        if (debug) cout <<"made dimuons...   "<< endl;
        
	  MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);

        if (debug) cout <<"filling muid   "<< selectedMuons[0]->Pt()<<endl;
        
	  //Fill Muon ID plots
        double muzPVz = fabs(selectedMuons[0]->vz() - vertex[0]->Z());
        MSPlot["MuonPt"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["Muond0"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MuonDistVzPVz"]->Fill(muzPVz, datasets[d], true, Luminosity*scaleFactor );
        MSPlot["MuonDz"]->Fill( selectedMuons[0]->dz() , datasets[d], true, Luminosity*scaleFactor );

        
        /*
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

         */
      if (debug) cout <<"filled all muID plots  .."<<endl;

   
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////Jets
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfSelectedBJets"]->Fill(selectedBJets.size(), datasets[d], true, Luminosity*scaleFactor);
	//  MSPlot["RhoCorrection"]->Fill(event->kt6PFJetsPF2PAT_rho(), datasets[d], true, Luminosity*scaleFactor);
	  if (debug) cout <<"per jet plots.."<<endl;

	//plots to to inspire staggered Jet Pt selection
      if (selectedJets.size()>=4) MSPlot["4thJetPt"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=5) MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=6) MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=7) MSPlot["7thJetPt"]->Fill(selectedJets[6]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	
    if (debug) cout <<"got muons and mets"<<endl;
        
        
        /////////////////////////////////
        /// Find indices of jets from Tops
        ////////////////////////////////
            
        for(unsigned int i=0; i<selectedJets.size(); i++) selectedJetsTLV.push_back(*selectedJets[i]);
        JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
        
        vector< pair<unsigned int, unsigned int> > JetPartonPair;
        for(unsigned int i=0; i<mcParticlesTLV.size(); i++) //loop through mc particles and find matched jets
        {
            int matchedJetNumber = matching.getMatchForParton(i, 0);
            if(matchedJetNumber != -1)
                JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );// Jet index, MC Particle index
        }
 
        
       if (debug) cout <<"n sel jets  "<<selectedJets.size()  << "   n mc particles tlv : "<< mcParticlesTLV.size() << " jet parton pari size :   "<< JetPartonPair.size()<<"  "<< muPlusFromTop<<muMinusFromTop<<endl;
        
        for(unsigned int i=0; i<JetPartonPair.size(); i++)//looping through matched jet-parton pairs
        {
            unsigned int j = JetPartonPair[i].second;	  //get index of matched mc particle
        
            
            if( fabs(mcParticlesMatching_[j]->type()) < 5 )
            {

                if( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -6 )
                   || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == 6 ) )
                {
                    if(hadronicWJet1_.first == 9999)
                    {
                        hadronicWJet1_ = JetPartonPair[i];
                       // MCPermutation[0] = JetPartonPair[i].first;
                    }
                    else if(hadronicWJet2_.first == 9999)
                    {
                        hadronicWJet2_ = JetPartonPair[i];
                        //MCPermutation[1] = JetPartonPair[i].first;
                    }
                    //else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
                }
            }
            
            else if( fabs(mcParticlesMatching_[j]->type()) == 5 )
            {
               

                if(  ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -6) || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 6 ) )
                {
                    hadronicBJet_ = JetPartonPair[i];
                    //MCPermutation[2] = JetPartonPair[i].first;
                }
                else if((muPlusFromTop && mcParticlesMatching_[j]->motherType() == 6) ||  ( muMinusFromTop &&mcParticlesMatching_[j]->motherType() == -6) )
                {
                    leptonicBJet_ = JetPartonPair[i];
                    //MCPermutation[3] = JetPartonPair[i].first;
                }
            }
            
        }
 
//        cout <<"setting w1  "<<endl;
   
  //  cout <<"  "<<endl;
   if (debug) cout <<"Indices of matched jets are :  "<< hadronicBJet_.first<<"  "<< hadronicWJet1_.first  <<" " << hadronicWJet2_.first <<endl;
   
        
        /////////////////////////////////
        /// TMVA for mass reconstruction
        ////////////////////////////////

        
        jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticles_flav, selectedJets,  genEvt_flav, scaleFactor);
        if (debug) cout <<"Processing event with jetcombiner :  "<< endl;
      
        if(!TrainMVA){
        MVAvals1 = jetCombiner->getMVAValue(MVAmethod, 1); // 1 means the highest MVA value
        MVAvals2 = jetCombiner->getMVAValue(MVAmethod, 8); // 2 means the highest MVA value
        selectedJets2ndPass.clear();
            
            

          //  jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticles_flav, selectedJets2ndPass,  genEvt_flav, scaleFactor);

            MVAvals2ndPass = jetCombiner->getMVAValue(MVAmethod, 1); // 2 means the highest MVA value

          //  cout <<"sel jets 1st pass = "<< selectedJets.size() << "sel jets 2nd pass =  " << selectedJets2ndPass.size() << endl;
            
         bestTopMass1 =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]).M();
         bestTopMass2 =( *selectedJets[MVAvals2.second[0]] + *selectedJets[MVAvals2.second[1]] + *selectedJets[MVAvals2.second[2]]).M();
         bestTopMass2ndPass =( *selectedJets[MVAvals2ndPass.second[0]] + *selectedJets[MVAvals2ndPass.second[1]] + *selectedJets[MVAvals2ndPass.second[2]]).M();

         bestTopPt =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]).Pt();

            
            
     if (debug)   cout <<"Indices of best MVA jets are :  "<< MVAvals1.second[0] <<"  "<< MVAvals1.second[1]  <<" " << MVAvals1.second[2]<<endl;

            
          //  cout <<"MVA Mass 1 = "<< bestTopMass1 << " MVA Mass 2 = "<< bestTopMass2 << endl;
            
        MSPlot["MVA1TriJetMass"]->Fill(bestTopMass1,  datasets[d], true, Luminosity*scaleFactor );
        MSPlot["MVA2TriJetMass"]->Fill(bestTopMass2,  datasets[d], true, Luminosity*scaleFactor );
        MSPlot["MVA2ndPassTriJetMass"]->Fill(bestTopMass2ndPass,  datasets[d], true, Luminosity*scaleFactor );

       if (debug)  cout <<"MVA Mass 1 = "<< bestTopMass1 << " MVA Mass 2 = "<< bestTopMass2 << endl;

            
        //////////////////////////////////////////////////////////////////////////////////
        /////          Calculating how well the MVA jet selection is doing         //////
        /////                                                                       //////
        //////////////////////////////////////////////////////////////////////////////////
        
        if(   ( hadronicBJet_.first == MVAvals1.second[0] || hadronicBJet_.first == MVAvals1.second[1] || hadronicBJet_.first == MVAvals1.second[2]   )  && ( hadronicWJet1_.first == MVAvals1.second[0] || hadronicWJet1_.first == MVAvals1.second[1] || hadronicWJet1_.first == MVAvals1.second[2]   )    && ( hadronicWJet2_.first == MVAvals1.second[0] || hadronicWJet2_.first == MVAvals1.second[1] || hadronicWJet2_.first == MVAvals1.second[2]   )      ){
            nMVASuccesses++;
        }
        
            if (debug)  cout <<"checked mva success..."<< endl;

            double matchedMass, unmatchedMass;
            
            // Pick some random jets to investigate wrong combinations...
            int badjet1 =-1 , badjet2 =-1 , badjet3=-1;
            
           while((badjet1==-1)||(badjet2==-1)||(badjet3==-1)) {
                Int_t rand_jet_index  = rand->Integer(selectedJets.size());
                if ((rand_jet_index != hadronicBJet_.first) ||(rand_jet_index != hadronicWJet1_.first ) || (rand_jet_index != hadronicWJet2_.first)){

                    if (badjet1== -1) {
                        badjet1 = rand_jet_index;
                        continue;
                    }
                    if ((badjet2== -1) && (rand_jet_index!= badjet1)) {
                        badjet2 = rand_jet_index;
                        continue;
    
                    }
                    if ((badjet3== -1)  && (rand_jet_index!= badjet1) &&( rand_jet_index!= badjet2)) {
                        badjet3 = rand_jet_index;
                        continue;
   
                    }
                }
                
            } 
            unmatchedMass =( *selectedJets[badjet1] + *selectedJets[badjet2] + *selectedJets[badjet3]).M();

            MSPlot["TriJetMass_UnMatched"]->Fill(unmatchedMass,  datasets[d], true, Luminosity*scaleFactor );

           // cout <<"Indices of matched jets are :  "<< hadronicBJet_.first<<"  "<< hadronicWJet1_.first  <<" " << hadronicWJet2_.first <<endl;
           // cout <<"Indices of random  jets are :  "<<badjet1 <<"  "<< badjet2  <<" " << badjet3  <<endl;
            //cout<<" "<<endl;
        

        if(   ( hadronicBJet_.first != 9999 )  && ( hadronicWJet1_.first != 9999   )    && ( hadronicWJet2_.first != 9999  )      ){
            nMatchedEvents++;
            if (debug) cout <<"matched..." <<endl;
            
            matchedMass =( *selectedJets[hadronicBJet_.first] + *selectedJets[hadronicWJet1_.first] + *selectedJets[hadronicWJet2_.first]).M();


            MSPlot["MVA1TriJetMassMatched"]->Fill(bestTopMass1,  datasets[d], true, Luminosity*scaleFactor );

            MSPlot["TriJetMass_Matched"]->Fill(matchedMass,  datasets[d], true, Luminosity*scaleFactor );
        }
        }
        
	       
//////////////////////////////////////////////////////////////////////////////////
/////          Calculate derived variables for rejecting ttbar + X          //////
/////                                                                       //////
///// 1. HT, 2. HTX, 3.EventMass, 4.EventMassX                              //////
//////////////////////////////////////////////////////////////////////////////////
        
        if (debug) cout <<"Caculating event level variables...  "<< endl;

        MEzCalculator NuPz;
        NuPz.SetMET(*mets[0]);
        NuPz.SetMuon(*selectedMuons[0]);
        if (debug) cout <<"Created NuPz object "<< endl;
        
        for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ){
            MSPlot["BdiscBJetCand_CSV"]->Fill(selectedJets[seljet1]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
            MSPlot["SelectedJetPt"]->Fill(selectedJets[seljet1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
            MSPlot["JetPhi"]->Fill(selectedJets[seljet1]->Phi() , datasets[d], true, Luminosity*scaleFactor);
            
            //Event-level variables
            jetpt = selectedJets[seljet1]->Pt();
            HT = HT + jetpt;
            H = H + selectedJets[seljet1]->P();
            sumpx = sumpx + selectedJets[seljet1]->Px();
            sumpy = sumpy + selectedJets[seljet1]->Py();
            sumpz = sumpz + selectedJets[seljet1]->Pz();
            sume = sume + selectedJets[seljet1]->E();  
            
            if(!TrainMVA){
            if(seljet1 == MVAvals1.second[0]  || seljet1 == MVAvals1.second[1] || seljet1 == MVAvals1.second[2]  ) continue;
            selectedJets2ndPass.push_back(selectedJets[seljet1]);
            HTX = HTX + jetpt;
            HX = HX + selectedJets[seljet1]->P();
            sumpx_X = sumpx_X + selectedJets[seljet1]->Px();
            sumpy_X = sumpy_X + selectedJets[seljet1]->Py();
            sumpz_X = sumpz_X + selectedJets[seljet1]->Pz();
            sume_X = sume_X + selectedJets[seljet1]->E();
            }
        }
        
        
        sort(selectedJets2ndPass.begin(),selectedJets2ndPass.end(),HighestCVSBtag()); //order muons wrt dsicriminator

        if (debug) cout <<"Creating sumjets "<< endl;

        
        TRootJet sumjet (TLorentzVector (sumpx, sumpy, sumpz,sume )); //Object representing all the jets summed
        TRootJet sumjet_X (TLorentzVector (sumpx_X, sumpy_X, sumpz_X,sume_X )); //Object representing all the jets summed minus the hadronic system

        //now extend sumjet object to form event and reduced event
        //1. add muon
        sumpx_X = sumpx_X + selectedMuons[0]->Px();
        sumpy_X = sumpy_X + selectedMuons[0]->Py();
        sumpz_X = sumpz_X + selectedMuons[0]->Pz();
        sume_X  = sume_X  + selectedMuons[0]->E();

        
        //2. add MET
        sumpx_X = sumpx_X + mets[0]->Px();
        sumpy_X = sumpy_X + mets[0]->Py();
        sumpz_X = sumpz_X + NuPz.Calculate();
        sume_X = sume_X + mets[0]->E();
        
        
        //now extend sumjet object to form event
        //1. add muon
        sumpx = sumpx + selectedMuons[0]->Px();
        sumpy = sumpy + selectedMuons[0]->Py();
        sumpz = sumpz + selectedMuons[0]->Pz();
        sume  = sume  + selectedMuons[0]->E();
        
        //2. add MET
        sumpx = sumpx + mets[0]->Px();
        sumpy = sumpy + mets[0]->Py();
        sumpz = sumpz + NuPz.Calculate();
        sume = sume + mets[0]->E();
    
        TLorentzVector bestMVATop;
        if(!TrainMVA) bestMVATop =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]);

        TRootJet Event (TLorentzVector (sumpx, sumpy, sumpz,sume )); //Object representing the whole event

        TRootJet Event_X (TLorentzVector (sumpx_X, sumpy_X, sumpz_X,sume_X )); //Object representing the event minus the hadronic top

    
        EventMass = Event.M();
        EventMassX = Event_X.M();
        SumJetMass = sumjet.M();
        SumJetMassX = sumjet_X.M();
        
        PTBalTopEventX = bestMVATop.Pt() - Event_X.Pt()  ;
        PTBalTopSumJetX = bestMVATop.Pt() - sumjet_X.Pt()  ;
        
        HTH = HT/H;
        HTXHX = HTX/HX;
        

        
         TLorentzVector mumet;

        
        //Form Lep W
         double   mumpx =   selectedMuons[0]->Px() + mets[0]->Px();
         double   mumpy =   selectedMuons[0]->Py() + mets[0]->Py();
         double   mumpz =   selectedMuons[0]->Pz() + NuPz.Calculate();
         double   mumpe =   selectedMuons[0]->E()  + mets[0]->E();
         mumet.SetPx(mumpx);
         mumet.SetPy(mumpy);
         mumet.SetPz(mumpz);
         mumet.SetE(mumpe);
        
        TLorentzVector muMETHadTop;
        TLorentzVector muMETB;
        TLorentzVector muMETBHadTop;
        
        if(!TrainMVA){
        muMETHadTop =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]  + mumet);
        muMETB =( *selectedJets2ndPass[0]  + mumet);
        muMETBHadTop =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]  + mumet + *selectedJets2ndPass[0] );
        }
        
        double deltaMTopMuMetB = muMETB.M() - bestMVATop.M() ;

        double deltaMTopMuMet = mumet.M() - bestMVATop.M() ;

        //cout <<"  mumetB mass  "<<  muMETB.M()   << "  MVA Tops Mass "<< bestMVATop.M()   <<  "  difference : "<< deltaMTopMuMetB <<endl;
        
        //cout<<" H  = "<< H  << " HT  = "<< HT <<" EventMass =   "<< EventMass<< " Event mass X  "<< EventMassX <<  "SumJet mass " <<  SumJetMass<< "SumJet mass X " << SumJetMassX  <<endl;
        
        if (debug) cout <<"arrays  "<<endl;
        
        MHT = sumjet.Pt();
        MHTSig = MHT/sqrt(HT);
        STJet = HT + mets[0]->Et() + selectedMuons[0]->Pt();
        EventMass = sumjet.M();
    
        sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag()); //order muons wrt dsicriminator
        
        double bTagRatio13 = selectedJets[0]->btag_combinedSecondaryVertexBJetTags()/selectedJets[2]->btag_combinedSecondaryVertexBJetTags();
        double bTagRatio14 = selectedJets[0]->btag_combinedSecondaryVertexBJetTags()/selectedJets[3]->btag_combinedSecondaryVertexBJetTags();
        
        sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt dsicriminator
        
        if (debug) cout <<"jets sorted...  "<<endl;

        
        if(trainEventMVA){
        
        if(dataSetName=="NP_overlay_TTTT"){

        Eventtrainer_->Fill("S","H", H);
        Eventtrainer_->Fill("S","HT", HT);
        Eventtrainer_->Fill("S","HTX", HTX);
        Eventtrainer_->Fill("S","HTH", HTH);
        Eventtrainer_->Fill("S","HTXHX", HTXHX);
        Eventtrainer_->Fill("S","EventMass", EventMass);
        Eventtrainer_->Fill("S","EventMassX", EventMassX);
        Eventtrainer_->Fill("S","SumJetMass", SumJetMass);
        Eventtrainer_->Fill("S","SumJetMassX", SumJetMassX);
        Eventtrainer_->Fill("S","PTBalTopEventX", PTBalTopEventX);
        Eventtrainer_->Fill("S","PTBalTopSumJetX", PTBalTopSumJetX);
        Eventtrainer_->Fill("S","bestTopMass2", bestTopMass2);
            
        Eventtrainer_->Fill("S","bTagRatio13", bTagRatio13);
        Eventtrainer_->Fill("S","bTagRatio14", bTagRatio14);
        Eventtrainer_->Fill("S","nJets",selectedJets.size()  );
        Eventtrainer_->Fill("S","Jet1Pt",selectedJets[0]->Pt());
        Eventtrainer_->Fill("S","Jet2Pt",selectedJets[1]->Pt());
        Eventtrainer_->Fill("S","Jet3Pt",selectedJets[2]->Pt());
        Eventtrainer_->Fill("S","Jet4Pt",selectedJets[3]->Pt());
        Eventtrainer_->Fill("S","Jet5Pt",selectedJets[4]->Pt());
        Eventtrainer_->Fill("S","Jet6Pt",selectedJets[5]->Pt());
        }
        
        if(dataSetName=="TTJets"){
            cout <<"jets sorted.  ttjets..  "<<endl;

            
            Eventtrainer_->Fill("B","H", H);
            cout <<"jets sorted.  ttjets.. cc0 "<<endl;

            Eventtrainer_->Fill("B","HT", HT);
            Eventtrainer_->Fill("B","HTX", HTX);
            Eventtrainer_->Fill("B","HTH", HTH);
            Eventtrainer_->Fill("B","HTXHX", HTXHX);
            Eventtrainer_->Fill("B","EventMass", EventMass);
            Eventtrainer_->Fill("B","EventMassX", EventMassX);
            Eventtrainer_->Fill("B","SumJetMass", SumJetMass);
            Eventtrainer_->Fill("B","SumJetMassX", SumJetMassX);
            Eventtrainer_->Fill("B","PTBalTopEventX", PTBalTopEventX);
            Eventtrainer_->Fill("B","PTBalTopSumJetX", PTBalTopSumJetX);
            Eventtrainer_->Fill("B","bestTopMass2", bestTopMass2);
            
            cout <<"jets sorted.  ttjets.. cc1 "<<endl;

            Eventtrainer_->Fill("B","bTagRatio13", bTagRatio13);
            Eventtrainer_->Fill("B","bTagRatio14", bTagRatio14);
            Eventtrainer_->Fill("B","nJets", selectedJets.size() );
            cout <<"jets sorted.  ttjets.. cc2 "<<endl;

            Eventtrainer_->Fill("B","Jet1Pt", selectedJets[0]->Pt() );
            cout <<"jets sorted.  ttjets.. cc 3"<<endl;

            Eventtrainer_->Fill("B","Jet2Pt", selectedJets[1]->Pt() );
            cout <<"jets sorted.  ttjets.. cc4 "<<endl;

            Eventtrainer_->Fill("B","Jet3Pt", selectedJets[2]->Pt() );
            cout <<"jets sorted.  ttjets.. cc 5"<<endl;

            Eventtrainer_->Fill("B","Jet4Pt", selectedJets[3]->Pt() );
            cout <<"jets sorted.  ttjets.. cc 6"<<endl;

            Eventtrainer_->Fill("B","Jet5Pt", selectedJets[4]->Pt() );
            cout <<"jets sorted.  ttjets.. cc 7"<<endl;

            Eventtrainer_->Fill("B","Jet6Pt", selectedJets[5]->Pt() );
            cout <<"jets sorted.  ttjets.. cc "<<endl;

            
        }
            
        }
        
        else if (computeEventMVA){
            if (debug) cout <<"filling computer...."<<H <<endl;

            if (Eventcomputer_ == 0) cout <<"null computer...."<<H <<endl;

            
        Eventcomputer_->FillVar("H", H);
            
        if (debug) cout <<"filling computer...."<<endl;

        Eventcomputer_->FillVar("HT", HT);
        Eventcomputer_->FillVar("HTX", HTX);
        Eventcomputer_->FillVar("HTH", HTH);
        Eventcomputer_->FillVar("HTXHX", HTXHX);
        Eventcomputer_->FillVar("EventMass", EventMass);
        Eventcomputer_->FillVar("EventMassX", EventMassX);
        Eventcomputer_->FillVar("SumJetMass", SumJetMass);
        Eventcomputer_->FillVar("SumJetMassX", SumJetMassX);
        Eventcomputer_->FillVar("PTBalTopEventX", PTBalTopEventX);
        Eventcomputer_->FillVar("PTBalTopSumJetX", PTBalTopSumJetX);
        Eventcomputer_->FillVar("bestTopMass2", bestTopMass2);
        Eventcomputer_->FillVar("bTagRatio13", bTagRatio13);
        Eventcomputer_->FillVar("bTagRatio14", bTagRatio14);
        Eventcomputer_->FillVar("nJets", selectedJets.size() );
            
            if (debug) cout <<"pre jets file///."<<endl;

        Eventcomputer_->FillVar("Jet1Pt", selectedJets[0]->Pt() );
        Eventcomputer_->FillVar("Jet2Pt", selectedJets[1]->Pt() );
        Eventcomputer_->FillVar("Jet3Pt", selectedJets[2]->Pt() );
        Eventcomputer_->FillVar("Jet4Pt", selectedJets[3]->Pt() );
        Eventcomputer_->FillVar("Jet5Pt", selectedJets[4]->Pt() );
        Eventcomputer_->FillVar("Jet6Pt", selectedJets[5]->Pt() );
            if (debug) cout <<"filled inside loop...."<<endl;

            
            }
        if (debug) cout <<"computer filled...."<<endl;
        
        double BDTscore;
        
        if(computeEventMVA){  std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
        
        for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it){
        
          //  cout <<"MVA Method : "<< it->first    <<" Score : "<< it->second <<endl;
            BDTscore = it->second;

        }
        
    }
       // cout <<"  "<<endl;
        
        
        MSPlot["MVA"]->Fill(BDTscore, datasets[d], true, Luminosity*scaleFactor);
        
        if (debug) cout <<"filling event level histos  "<<endl;

        
        MSPlot["EventMass"]->Fill(EventMass, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["EventMassX"]->Fill(EventMassX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["SumJetMass"]->Fill(SumJetMass, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["SumJetMassX"]->Fill(SumJetMassX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTX"]->Fill(HTX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["H"]->Fill(H, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HX"]->Fill(HX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTXHX"]->Fill(HTXHX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTH"]->Fill(HTH, datasets[d], true, Luminosity*scaleFactor);
        
        MSPlot["deltaMTopMuMet"]->Fill(deltaMTopMuMet, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["deltaMTopMuMetB"]->Fill(deltaMTopMuMetB, datasets[d], true, Luminosity*scaleFactor);

        MSPlot["MuMetBMasses"]->Fill(muMETB.M(), datasets[d], true, Luminosity*scaleFactor);
        
        MSPlot["PTBalTopEventX"]->Fill(PTBalTopEventX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["PTBalTopSumJetX"]->Fill(PTBalTopSumJetX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["PTBalTopMuMet"]->Fill(muMETHadTop.Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["PTBalTopMuMetB"]->Fill(muMETBHadTop.Pt(), datasets[d], true, Luminosity*scaleFactor);

        
        
        MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
	    MSPlot["MHTSig_SelectedJets"]->Fill(MHTSig, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MET_MHT"]->Fill(mets[0]->Et()/MHT, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["STJet"]->Fill( STJet , datasets[d], true, Luminosity*scaleFactor);
        
        
        if (debug) cout <<"filled event level histos  "<<endl;
        

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
          string H_Name = "H_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string HX_Name = "HX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string HTH_Name = "HTH_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string HTXHX_Name = "HTXHX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

          
          string EventMass_Name = "EventMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string EventMassX_Name = "EventMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string FirstTriJetMass_1BJet_g_Name = "FirstTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string FirstDiJetMass_1BJet_g_Name = "FirstDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string SecondDiJetMass_1BJet_g_Name = "SecondDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

          string SumJetMass_Name = "SumJetMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string SumJetMassX_Name = "SumJetMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string PTBalTopEventX_Name = "PTBalTopEventX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string PTBalTopSumJetX_Name = "PTBalTopSumJetX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

          string MuMetBMasses_g_Name = "MuMetBMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string MuMetMasses_g_Name = "MuMetMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string MET_Name = "MET"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

          string SecondTriJetMass_1BJet_g_Name = "SecondTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string SecondTriJetMass_1BJet_g_chi2_Name = "SecondTriJetMass_1BJet_g_chi2_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string MVA1TriJetMass_Name = "MVA1TriJetMass"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          string MVA2TriJetMass_Name = "MVA2TriJetMass"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

          
          
        if(b<=8 && c<=2){
            if(selectedJets.size() == b && selectedBJets.size() == c  ) {
                MSPlot[HT_Name.c_str() ]->Fill(HT,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[HTX_Name.c_str() ]->Fill(HTX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[H_Name.c_str() ]->Fill(H,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[HX_Name.c_str() ]->Fill(HX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[HTH_Name.c_str() ]->Fill(HTH,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[HTXHX_Name.c_str() ]->Fill(HTXHX,datasets[d], true, Luminosity*scaleFactor);

                
                MSPlot[EventMass_Name.c_str() ]->Fill(EventMass,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[EventMassX_Name.c_str() ]->Fill(EventMassX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[SumJetMass_Name.c_str() ]->Fill(SumJetMass,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[SumJetMassX_Name.c_str() ]->Fill(SumJetMassX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[PTBalTopEventX_Name.c_str() ]->Fill(PTBalTopEventX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[PTBalTopSumJetX_Name.c_str() ]->Fill(PTBalTopSumJetX,datasets[d], true, Luminosity*scaleFactor);
            
             //   MSPlot[FirstTriJetMass_1BJet_g_Name.c_str()]->Fill(fir_top_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[FirstDiJetMass_1BJet_g_Name.c_str()]->Fill(fir_w_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);
               // MSPlot[SecondDiJetMass_1BJet_g_Name.c_str()]->Fill(sec_w_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[SecondTriJetMass_1BJet_g_Name.c_str()]->Fill(sec_top_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor );
              //  MSPlot[SecondTriJetMass_1BJet_g_chi2_Name.c_str()]->Fill(firstchi2_g ,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot[MET_Name.c_str()]->Fill(mets[0]->Et() ,  datasets[d], true, Luminosity*scaleFactor );
               // MSPlot[MuMetMasses_g_Name.c_str()]->Fill(lep_w_mass_chosen_g ,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot[MVA1TriJetMass_Name.c_str()]->Fill(bestTopMass1,  datasets[d], true, Luminosity*scaleFactor);
                MSPlot[MVA2TriJetMass_Name.c_str()]->Fill(bestTopMass2,  datasets[d], true, Luminosity*scaleFactor);
                      }
                        }
        else if (b>8 && c > 2) {
            if(selectedJets.size() >= b && selectedBJets.size() >= c  ) {
                MSPlot[HT_Name.c_str() ]->Fill(HT,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[HTX_Name.c_str() ]->Fill(HTX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[H_Name.c_str() ]->Fill(H,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[HX_Name.c_str() ]->Fill(HX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[HTH_Name.c_str() ]->Fill(HTH,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[HTXHX_Name.c_str() ]->Fill(HTXHX,datasets[d], true, Luminosity*scaleFactor);
                
                MSPlot[EventMass_Name.c_str() ]->Fill(EventMass,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[EventMassX_Name.c_str() ]->Fill(EventMassX,datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[FirstTriJetMass_1BJet_g_Name.c_str()]->Fill(fir_top_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[FirstDiJetMass_1BJet_g_Name.c_str()]->Fill(fir_w_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[SecondDiJetMass_1BJet_g_Name.c_str()]->Fill(sec_w_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor);
                MSPlot[MET_Name.c_str()]->Fill(mets[0]->Et() ,  datasets[d], true, Luminosity*scaleFactor );
                
                MSPlot[SumJetMass_Name.c_str() ]->Fill(SumJetMass,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[SumJetMassX_Name.c_str() ]->Fill(SumJetMassX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[PTBalTopEventX_Name.c_str() ]->Fill(PTBalTopEventX,datasets[d], true, Luminosity*scaleFactor);
                MSPlot[PTBalTopSumJetX_Name.c_str() ]->Fill(PTBalTopSumJetX,datasets[d], true, Luminosity*scaleFactor);
                
                //MSPlot[SecondTriJetMass_1BJet_g_Name.c_str()]->Fill(sec_top_mass_chosen_g,  datasets[d], true, Luminosity*scaleFactor );
                //MSPlot[SecondTriJetMass_1BJet_g_chi2_Name.c_str()]->Fill(firstchi2_g ,  datasets[d], true, Luminosity*scaleFactor );
               // MSPlot[MuMetBMasses_g_Name.c_str()]->Fill(lep_top_mass_chosen_g ,  datasets[d], true, Luminosity*scaleFactor );
               // MSPlot[MuMetMasses_g_Name.c_str()]->Fill(lep_w_mass_chosen_g ,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot[MVA1TriJetMass_Name.c_str()]->Fill(bestTopMass1,  datasets[d], true, Luminosity*scaleFactor);
                MSPlot[MVA2TriJetMass_Name.c_str()]->Fill(bestTopMass2,  datasets[d], true, Luminosity*scaleFactor);


            }
        }
}
}
 
        //histo2D[("ThirdTopMass_vs_HT_"+datasets[d]->Name()).c_str()]->Fill(sec_top_mass_chosen, HT);
      //  histo2D[("ThirdTopMass_vs_SecondTopMass"+datasets[d]->Name()).c_str()]->Fill(sec_top_mass_chosen,fir_top_mass_chosen_g );
        //histo2D[("MassChi2_vs_HT"+datasets[d]->Name()).c_str()]->Fill(firstchi2,HT );
        histo2D[("EventMassX_vs_HT"+datasets[d]->Name()).c_str()]->Fill(EventMassX, HT);
        histo2D[("EventMassX_vs_HTX"+datasets[d]->Name()).c_str()]->Fill(EventMassX, HTX);
        
    
        
        

        if(!datadriven)
        {
           
    //       myBranch_selectedEvents = new FourTopTree();
            //myBranch_selectedEvents->setEventID( 100 );
            //myBranch_selectedEvents->setRunID( 100 );
            //yBranch_selectedEvents->setLumiBlockID( 100 );
           // myBranch_selectedEvents->setNPV( vertex.size() );
            //myBranch_selectedEvents->setNPUBXm1( event->nPu(-1) );
        }
        
        
        if (debug) cout <<"filling tree  "<<endl;        
       // myInclFourthGenTree->Fill(); 
        
        if (debug) cout <<"filling tree  "<<endl;        
    //    delete myBranch_selectedEvents;
        
	}
	else if(Electron){
	  //	MSPlot["1stLeadingElectronRelIsolation"]->Fill(selectedElectrons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);

	}
        
        
    //myTree->Fill();
        
    }//loop on events
    
    if (debug)cout <<"N BBar = = " << nBBBar <<"  N CCBar = = " << nCCBar <<"  N LLBar = =  " << nLLBar << endl;
      
     cout <<"setting tree"<<endl;
      string crannfilename = dataSetName + "_Crann.root";
      
     // TFile *crann = new TFile (crannfilename.c_str(),"RECREATE");
      //crann->cd();
      //myTree->Write();

      //crann->Write();
      //crann->Close();
      
    //  treeFile->Write();
      //treeFile->Close();
    //  delete treeFile;
      
    if(jetTools) delete jetTools;
    
    
      //important: free memory
      treeLoader.UnLoadDataset();
      
      double nMVASuccessesd = nMVASuccesses;
      double nMatchedEventsd = nMatchedEvents;
      
      
      cout <<"Efficiency of MVA jet selection = = "<<  nMVASuccessesd/nMatchedEventsd   <<endl;
      
  } //loop on datasets
    

    if(trainEventMVA) Eventtrainer_->TrainMVA("Block","",0,0,"",0,0,"test");
    
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
  	selecTableMu.Write(  "FourTop"+postfix+"Table_Mu.tex",    true,false,true,true,false,false,false);
  }
    else if(Electron){
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	selecTableEl.TableCalculator(  false, true, true, true, true);
   //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
    selecTableEl.Write(  "FourTop"+postfix+"Table_El.tex",  true,true,true,true,false,false,true);

  }
    
    fout->cd();
   TFile *foutmva = new TFile ("foutMVA.root","RECREATE");
   //foutmva->cd();
    cout <<" after cd .."<<endl;
    
    //Output ROOT file
    string mvarootFileName ("MVA"+postfix+channelpostfix+".root");
  //  TFile * foutmva;
  //  TFile* foutmva = new TFile(mvarootFileName.c_str(), "RECREATE");
    //foutmva->cd();
    
   // TH1F * h1 = new TH1F();
   // h1->Write();
    //foutmva->Write();
  //  foutmva->Close();
    
    string pathPNGJetCombi = pathPNG + "JetCombination/";
    mkdir(pathPNGJetCombi.c_str(),0777);
    
   // if (foutmva) cout <<"fout "<<endl;
    //cout <<"pointer to file -->"<< foutmva << endl;
    //foutmva->Print();
    
    if(TrainMVA)jetCombiner->Write(foutmva, true, pathPNGJetCombi.c_str());
   
//   jetCombiner->Write(fout, true, pathPNGJetCombi.c_str());
    
    
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        
        MultiSamplePlot *temp = it->second;
        TH1F *tempHisto_data;
        TH1F *tempHisto_TTTT;

      //Option to write ROOT files containing histograms with systematics varied +/-
      //TFile * tempErrorFile;
        
        if (doJESShift == 1){
            string filename = "ErrorBands/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
            TH1F * tempHisto = temp->getTH1F("TTJets");
            tempHisto->Write("Minus");
           // tempErrorFile->Write();
            tempErrorFile->Close();
        }
        else if  (doJESShift ==2){
            
            string filename = "ErrorBands/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
            TH1F * tempHisto = temp->getTH1F("TTJets");
            tempHisto->Write("Plus");
            tempErrorFile->Write();
            tempErrorFile->Close();

            cout <<"JES sys down"<<endl;
            
        }
        else if  (doJESShift ==0){
            cout <<"JES sys off "<<endl;
      

	//	temp->addText("CMS preliminary");
	
	temp->Draw(false, name, true, true, true, true, true,100.,true); // merge TT/QCD/W/Z/ST/
	//Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
      
      cout <<" looping plots..., name ... "<< name<<endl;
        
        temp->Write(fout, name, true, pathPNG, "pdf");
    }
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




int Factorial(int N = 1)
{
    int fact = 1;
    for( int i=1; i<=N; i++ )
        fact = fact * i;  // OR fact *= i;
    
    return fact;
}


