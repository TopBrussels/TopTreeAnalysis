#include "TStyle.h"
#include "TF2.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"

// RooFit librairies

#include "RooArgSet.h"
#include "RooAddition.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooFitResult.h"

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
#include "../Tools/interface/MVATrainer.h"
#include "../Tools/interface/MVAComputer.h"
#include "../Tools/interface/JetTools.h"
#include "../JESMeasurement/interface/JetCombiner.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../Reconstruction/interface/MakeBinning.h"
#include "../Reconstruction/interface/TTreeObservables.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../InclFourthGenSearch/interface/InclFourthGenTree.h"
#include "../InclFourthGenSearch/interface/InclFourthGenSearchTools.h"
#include "../InclFourthGenSearch/interface/TwoDimTemplateTools.h"
//#include "../MCInformation/interface/LumiReWeighting.h" 
//for Kinematic Fit
#include "../MCInformation/interface/ResolutionFit.h"
#include "../KinFitter/interface/TKinFitter.h"
#include "../KinFitter/interface/TFitConstraintM.h"
#include "../KinFitter/interface/TFitParticleEtThetaPhi.h"

#include "Style.C"

using namespace std;
using namespace TopTree;
using namespace RooFit;
//using namespace reweight;

//defined in InclFourthGenSearchTools.h
//std::string IntToStr( int n ){	std::ostringstream result;	result << n;	return result.str();}

struct sort_pair_decreasing
{
    bool operator()(const std::pair<int,float> &left, const std::pair<int,float> &right)
    {
        return left.second > right.second;
    }
};

//To cout the Px, Py, Pz, E and Pt of objects
void coutObjectsFourVector(vector < TRootMuon* > init_muons, vector < TRootElectron* > init_electrons, vector < TRootJet* > init_jets, vector < TRootMET* > mets, string Comment);

float jetprob(float jetpt, float btagvalue);

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

int main (int argc, char *argv[])
{
  //which systematic to run?
  string systematic = "Nominal";
  if (argc >= 2)
		systematic = string(argv[1]);
  cout << "Systematic to be used: " << systematic << endl;
  if( ! (systematic == "Nominal"  || systematic == "JESPlus" || systematic == "JESMinus" || systematic == "JERPlus" || systematic == "JERMinus" || systematic == "bTagPlus" || systematic == "bTagMinus" || systematic == "misTagPlus" || systematic == "misTagMinus" || systematic == "PUPlus" || systematic == "PUMinus"))
  {
    cout << "Unknown systematic!!!" << endl;
    cout << "Possible options are: " << endl;
    exit(-1);
  }	

  string btagger = "TCHPM";
//from BTV-11-001; WARNING: btag SFs are for jets of 20-240 GeV, averaged over eta; and mistag SFs are for jets of 50-80 GeV
// b-tag scalefactor => TCHEL: data/MC scalefactor = 0.95 +- 0.10,    TCHEM: data/MC scalefactor = 0.94 +- 0.09,	TCHPM: data/MC scalefactor =  0.91 +- 0.09
// mistag scalefactor => TCHEL: data/MC scalefactor = 1.11 +- 0.12,    TCHEM: data/MC scalefactor = 1.21 +- 0.17,	TCHPM: data/MC scalefactor =  1.27 +- 0.15
  float scalefactorbtageff = 1, mistagfactor = 1;
  if(btagger == "TCHEL") //track counting high eff loose working point
  {	  
		scalefactorbtageff = 0.95;
	  if(systematic == "bTagMinus")
			scalefactorbtageff = 0.85;
	  if(systematic == "bTagPlus")
			scalefactorbtageff = 1.05;
	  
		mistagfactor = 1.11;
	  if(systematic == "misTagMinus")
			mistagfactor = 0.99;
	  if(systematic == "misTagPlus")
			mistagfactor = 1.23;
		
  }
  else if(btagger == "TCHEM") //track counting high eff medium working point
  {
		scalefactorbtageff = 0.94;
	  if(systematic == "bTagMinus")
			scalefactorbtageff = 0.85;
	  if(systematic == "bTagPlus")
			scalefactorbtageff = 1.03;
		
		mistagfactor = 1.21;
	  if(systematic == "misTagMinus")
			mistagfactor = 1.04;
	  if(systematic == "misTagPlus")
			mistagfactor = 1.38;
  }
	else if(btagger == "TCHPM") //track counting high pur medium working point
  {
		scalefactorbtageff = 0.91;
	  if(systematic == "bTagMinus")
			scalefactorbtageff = 0.8194;
	  if(systematic == "bTagPlus")
			scalefactorbtageff = 1.0006;
		
		mistagfactor = 1.27;
	  if(systematic == "misTagMinus")
			mistagfactor = 1.1187;
	  if(systematic == "misTagPlus")
			mistagfactor = 1.4213;
  }
	else if(btagger == "TCHPT") //track counting high pur tight working point
  {
	  cout<<"WARNING: look up SFs for TCHPT"<<endl;
  }
  float workingpointvalue = 9999; //{1.7,3.3,10.2}; trackcountinghighefficiency working points: loose, medium, tight
  if(btagger == "TCHEL")
     workingpointvalue = 1.7;
  else if(btagger == "TCHEM")
     workingpointvalue = 3.3;
  else if(btagger == "TCHPM")
     workingpointvalue = 1.93;
	else if(btagger == "TCHPT")
     workingpointvalue = 3.41;

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the fourth generation search ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  setTDRStyle();
  //setMyStyle();

  string postfixOld = "_Fall11_Round4"; // should be same as postifix in TreeCreator of the trees
	string postfix= postfixOld+"_"+systematic;

  string Treespath = "InclFourthGenTrees_Fall11_Round4";
  Treespath = Treespath +"/"; 		
  //mkdir(TreespathPNG.c_str(),0777);
	bool savePNG = false;

  /////////////////////
  // Configuration
  /////////////////////
  bool useMassesAndResolutions = true;
  bool doMVAjetcombination = true; //when false, the jet combination and the top mass will not be reconstructed, and nothing will be trained
  bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff. Note: there is an MVA output file with the training, but also some files in the ./weights directory!!
  if (systematic != "Nominal")
  {
		useMassesAndResolutions = true;
		doMVAjetcombination = true;
		TrainMVA = false;
  }
  bool doPUreweighting = true;
	bool doMCtriggerEffreweighting = false;
	bool doBtagSFreweighting = true;
	bool doRun2011AOnly = false;
	bool doRun2011BOnly = false;
	if(doRun2011AOnly) cout<<" RUNNING ON RUN2011A ONLY: MAKE SURE TO PUT IN CORRECT LUMINOSITY IN XML CONFIG!!"<<endl;
	else if(doRun2011BOnly) cout<<" RUNNING ON RUN2011B ONLY: MAKE SURE TO PUT IN CORRECT LUMINOSITY IN XML CONFIG!!"<<endl;
	
  //bool TrainwithTprime = false; //temporarily not supported
  string MVAmethod = "Likelihood"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)
  string channelpostfix = "";
  bool semiElectron = false; // use semiElectron channel?
  bool semiMuon = true; // use semiMuon channel?
  if(semiElectron && semiMuon)
  {
     cout << "  --> Using both semiMuon and semiElectron channel? Choose only one (for the moment, since this requires running on different samples/skims)!" << endl;
     exit(1);
  }
  else
  {
    if(semiMuon){
       cout << " --> Using the semiMuon channel..." << endl;
       channelpostfix = "_semiMu";
    }
    else if(semiElectron){
       cout << " --> Using the semiElectron channel..." << endl;
       channelpostfix = "_semiEl";
    }
  }
  
  //xml file
  string xmlFileName = "";
  //if(semiElectron) xmlFileName = "../config/myFourthGenconfig_Electron.xml";
  //else if(semiMuon) xmlFileName = "../config/myFourthGenconfig.xml";
	if(semiElectron) xmlFileName = "../config/myFourthGenconfig_Electron_Fall11.xml";
  else if(semiMuon) xmlFileName = "../config/myFourthGenconfig_Muon_Fall11.xml";	
  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    
  
  
  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = anaEnv.Verbose;
  float anaEnvLuminosity = anaEnv.Luminosity;	// in 1/pb 
  cout << "analysis environment luminosity for rescaling "<< anaEnvLuminosity << endl;

  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  
  //is this block needed?
  float Luminosity = anaEnvLuminosity;
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    	string dataSetName = datasets[d]->Name();
	if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
	{
		  Luminosity = datasets[d]->EquivalentLumi();
		  break;
	 }
  }
  if(Luminosity != anaEnvLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  
	
	vector<string> inputTrees; //fill with the tree files you want to read!!!
  //inputTrees.push_back();
	for (unsigned int d = 0; d < datasets.size(); d++) //d < datasets.size()
  {
    cout << "   Dataset " << d << " name : " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;    
    string dataSetName = datasets[d]->Name();
		string inputTreeFileName; //should follow convention of TreeFileName in InclFourthGen_TreeCreator.cc
		if(systematic == "JESPlus" || systematic == "JESMinus" || systematic == "JERPlus" || systematic == "JERMinus") inputTreeFileName = Treespath+"InclFourthGenTree_"+dataSetName+postfix+channelpostfix+".root";
		else inputTreeFileName = Treespath+"InclFourthGenTree_"+dataSetName+postfixOld+"_Nominal"+channelpostfix+".root";
		inputTrees.push_back(inputTreeFileName);
	}
	
  //Output ROOT file
  string rootFileName ("InclFourthGenSearch_TreeAnalyzer"+postfix+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

 
  //vector of objects
  cout << " - Variable declaration ..." << endl;

  //Global variable

  string pathPNG = "InclFourthGenSearchPlots_TreeAnalyzer"+postfix+channelpostfix;
  pathPNG = pathPNG +"/"; 	
  pathPNG = pathPNG +"/"; 	
  if(savePNG) mkdir(pathPNG.c_str(),0777);

  //Most 1D and MS plots are declared inside this class
  InclFourthGenSearchTools myInclFourthGenSearchTools(semiMuon, semiElectron, datasets, Luminosity,false); //last argument is doKinematicFit
  
  MSPlot["allDiJetMasses"] = new MultiSamplePlot(datasets, "allDiJetMasses", 50, 0, 1000, "m_{jj}");
  MSPlot["hadronicRecoWMass_chosenWjets"] = new MultiSamplePlot(datasets, "hadronicRecoWMass_chosenWjets", 50, 0, 500, "m_{W}"); 
  histo1D["hadronicPartonWMass"] = new TH1F("hadronicPartonWMass","Hadronic W Mass, using the Partons",100,0,200);
  histo1D["hadronicRecoWMass"] = new TH1F("hadronicRecoWMass","Hadronic W Mass, using the RecoJets",100,0,200);
  
  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
  histo1D["LeptonPt_TTbar"] = new TH1F("leptonspt ttbar","leptonspt ttbar;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_Tprime500"] = new TH1F("leptonspt tprime500","leptonspt tprime500;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_Bprime500"] = new TH1F("leptonspt bprime500","leptonspt bprime500;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_SBprime500"] = new TH1F("leptonspt sbprime500","leptonspt sbprime500;pt leptons;#events",250,0,500);
  
  string multileptons[7] = {"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"};
  string histoName,histo_dataset;
  for(int i = 0; i<7; i++)
  {
		histoName = "NbEvents_"+multileptons[i];
		for(unsigned int d = 0; d < datasets.size (); d++){
			histo_dataset = histoName+(datasets[d]->Name()).c_str(); 
			histo1D[histo_dataset.c_str()] = new TH1F(histo_dataset.c_str(),histo_dataset.c_str(), 1, 0.5, 1.5);
		}
  }
	
  MSPlot["MS_NbSSevents"] = new MultiSamplePlot(datasets,"# events with SS leptons", 1, 0.5, 1.5, "");
  MSPlot["MS_NbSSElElevents"] = new MultiSamplePlot(datasets,"# events with SS electrons", 1, 0.5, 1.5, "");
  MSPlot["MS_NbSSElMuevents"] = new MultiSamplePlot(datasets,"# events with SS electron+muon", 1, 0.5, 1.5, "");
  MSPlot["MS_NbSSMuMuevents"] = new MultiSamplePlot(datasets,"# events with SS muons", 1, 0.5, 1.5, "");
  MSPlot["MS_NbTrievents"] = new MultiSamplePlot(datasets,"# events with 3 leptons", 1, 0.5, 1.5, "");
  MSPlot["MS_MET"] = new MultiSamplePlot(datasets,"MET", 75, 0, 200, "");
	MSPlot["MS_METafterbtagSF"] = new MultiSamplePlot(datasets,"MET", 75, 0, 200, "");
  MSPlot["MS_LeptonPt"] = new MultiSamplePlot(datasets,"lepton pt", 150, 0, 300, "");
  MSPlot["MS_nPV"] = new MultiSamplePlot(datasets, "nPrimaryVertices", 21, -0.5, 20.5, "Nr. of primary vertices");
  MSPlot["MS_JetMultiplicity_SingleLepton"] = new MultiSamplePlot(datasets, "JetMultiplicity", 10, -0.5, 9.5, "Jet Multiplicity");
  MSPlot["MS_BtaggedJetMultiplicity_SingleLepton"] = new MultiSamplePlot(datasets, "BtaggedJetMultiplicity", 7, -0.5, 6.5, "b-tagged jet multiplicity");
	MSPlot["MS_JetMultiplicityAtleast1Btag_SingleLepton"] = new MultiSamplePlot(datasets, "JetMultiplicityAtleast1Btag", 10, -0.5, 9.5, "Jet multiplicity (>=1 b-tag)");

  MSPlot["MS_JetPt_all_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_all", 50, 0, 300, "Pt of all jets (GeV)");
	MSPlot["MS_JetPt_btagged_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_btagged", 50, 0, 300, "Pt of b-tagged jets (GeV)");
	MSPlot["MS_JetPt_nonbtagged_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_nonbtagged", 50, 0, 300, "Pt of non b-tagged jets (GeV)");
	
	
	cout << " - Declared histograms ..." <<  endl;

  float NbSSevents = 0;   
  float NbTrievents = 0;  

  /////////////////////////////////////////////////////////
  //Configuration and variables for 2D HT-Mtop distribution. There is a 2D distribution for 2 boxes seperately: 1B_2W and 2B_2W
  /////////////////////////////////////////////////////////
  string xvariable = "HT", yvariable = "MTop"; //these are the two variables for which the 2D plane is made
  int nbinsxvariable = 6, nbinsyvariable = 12; //if the binning is already created, make sure these are the same as before!
  string binningFileName_HTvsMTop_1B_2W = "Binning_InclFourthGenSearch_1B_2W_TTbarJetsFlat"+channelpostfix+".root";
  string binningFileName_HTvsMTop_2B_2W = "Binning_InclFourthGenSearch_2B_2W_TTbarJetsFlat"+channelpostfix+".root";
  TwoDimTemplateTools HTvsMTop_1B_2W("1B_2W",xvariable,nbinsxvariable,yvariable,nbinsyvariable);
  TwoDimTemplateTools HTvsMTop_2B_2W("2B_2W",xvariable,nbinsxvariable,yvariable,nbinsyvariable); 
  if(doMVAjetcombination && !TrainMVA && useMassesAndResolutions) //Note: the boolean 'useMassesAndResolutions'is not needed when a Wmass file exists
  { 
    HTvsMTop_1B_2W.SetDatasets(datasets);
    HTvsMTop_2B_2W.SetDatasets(datasets);
    if(anaEnv.MCRound!=0)
    {
       HTvsMTop_1B_2W.LoadTwoDimBinning(binningFileName_HTvsMTop_1B_2W);
       HTvsMTop_2B_2W.LoadTwoDimBinning(binningFileName_HTvsMTop_2B_2W);
    }
    else
    {
       cout<<" ... The 2D binnings will be created!"<<endl;
    }
  }	
	
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTableSemiLep;
  CutsSelecTableSemiLep.push_back(string("selected")); //0
  CutsSelecTableSemiLep.push_back(string("box 1B 1W"));
  CutsSelecTableSemiLep.push_back(string("box 1B 2W"));
  CutsSelecTableSemiLep.push_back(string("box 1B 3W"));
  CutsSelecTableSemiLep.push_back(string("box 1B $\\geq$ 4W")); //4
  CutsSelecTableSemiLep.push_back(string("box 2B 1W"));
  CutsSelecTableSemiLep.push_back(string("box 2B 2W"));
  CutsSelecTableSemiLep.push_back(string("box 2B 3W"));
  CutsSelecTableSemiLep.push_back(string("box 2B $\\geq$ 4W")); //8

  vector<string> CutsSelecTableMultiLep;
  CutsSelecTableMultiLep.push_back(string("SS leptons in all boxes combined")); //0
  CutsSelecTableMultiLep.push_back(string("trileptons in all boxes combined")); //1
  
	SelectionTable selecTableSemiLep(CutsSelecTableSemiLep, datasets);
  selecTableSemiLep.SetLuminosity(Luminosity);
  selecTableSemiLep.SetPrecision(1);

	SelectionTable selecTableMultiLep(CutsSelecTableMultiLep, datasets);
  selecTableMultiLep.SetLuminosity(Luminosity);
  selecTableMultiLep.SetPrecision(2);
  
  cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////////////////////
  // PileUp Reweighting - 3D//
  ////////////////////////////////////////////////////  
	Lumi3DReWeighting Lumi3DWeights;
	//if(!(doRun2011AOnly || doRun2011BOnly)) Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  //else if(doRun2011AOnly) Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun173692.root", "pileup", "pileup");
	if(!(doRun2011AOnly || doRun2011BOnly)) Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  else if(doRun2011AOnly) Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun173692.root", "pileup", "pileup");	
	else if(doRun2011BOnly)
	{
	  cout<<"WARNING: currently no histogram for 2011B only available: using Run2011A+B pilup reweighting!!"<<endl;
		//Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root","pileup","pileup");
		Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root","pileup","pileup");
	}
	
	Lumi3DWeights.weight3D_init(1.0);
	if(systematic == "PUMinus")
  	Lumi3DWeights.weight3D_init(0.92);
	else if(systematic == "PUPlus")
  	Lumi3DWeights.weight3D_init(1.08);
	else
  	Lumi3DWeights.weight3D_init(1.0);

//  LumiReWeighting LumiWeights = LumiReWeighting("PileUpReweighting/pileup_WJets_36bins.root", "PileUpReweighting/pileup_2011Data_UpToRun180252.root", "pileup2", "pileup");
//  PoissonMeanShifter PShiftUp = PoissonMeanShifter(0.6); // PU-systematic
//  PoissonMeanShifter PShiftDown = PoissonMeanShifter(-0.6); // PU-systematic
  cout << " - Initialized LumiReWeighting stuff" << endl;
  
  
  ///////////////
  // JetCombiner
  ///////////////
  JetCombiner* jetCombiner;
  if(!doMVAjetcombination) TrainMVA = false;
  else if(doMVAjetcombination)
  {
    jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, true, "",channelpostfix); //boolean is basically to use also the W mass as constraint    
  }
  

  if(doMVAjetcombination && TrainMVA) useMassesAndResolutions = true; //just to make sure the W mass plot is not produced (is not used anyway)
  
  cout << " - JetCombiner instantiated ..." << endl;

    
  ////////////////////////////////////
  ////////////////////////////////////
  ///////// Loop on datasets
  ////////////////////////////////////
  ////////////////////////////////////
  cout << " - Loop over datasets ... " << inputTrees.size() << " datasets !" << endl;

  for (unsigned int d = 0; d < inputTrees.size(); d++) //d < datasets.size()
  {
    if (verbose > 1)
      cout << "file: " << inputTrees[d] << endl;
	 
    TFile* inFile = new TFile(inputTrees[d].c_str(),"READ");
    TTree* inConfigTree = (TTree*) inFile->Get("configTreeFile");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);	//will not be used (problem with MVA 1/3 and 2/3 splitting of sample)	
		
		string dataSetName = datasets[d]->Name();
		
		TTree* inInclFourthGenTree = (TTree*) inFile->Get("myInclFourthGenTree");
    TBranch* m_br = (TBranch*) inInclFourthGenTree->GetBranch("InclFourthGenBranch_selectedEvents");
    //TBranch* m_br = (TBranch*) inInclFourthGenTree->GetBranch("TheInclFourthGenTree");
    InclFourthGenTree* myBranch_selectedEvents = 0;
    m_br->SetAddress(&myBranch_selectedEvents);
    int nEvents = inInclFourthGenTree->GetEntries();

		cout << " Processing DataSet: " << dataSetName << "  containing " << nEvents << " events" << endl;
    cout << " Cross section = " << datasets[d]->Xsection() << "  intLumi = " << datasets[d]->EquivalentLumi() << "  NormFactor = " << datasets[d]->NormFactor() << endl;

	 	if(!useMassesAndResolutions && dataSetName != "TTbarJets_SemiMuon" && dataSetName != "TTbarJets_SemiElectron") 
			continue;
		
    ///////////////////////////////////////////////////
    // calculate W reconstructed mass and resolution: load the mass values and resolutions for chi2 jetcombination.
		// Can this be put outside dataset loop?
    ///////////////////////////////////////////////////
    
    float Wmass = -9999;
    float SigmaWmass = -9999;
    
    if(useMassesAndResolutions)
    {
      string resfilename = "MassPlot.root";
      
      cout << "INFO: Using mass and width for reconstructed W mass from: " << resfilename << endl;
      TFile* res = new TFile(resfilename.c_str(),"READ");
      res->cd();
	
      TF1* WmassFit = (TF1*)res->Get("hadronicRecoWMass_Fitted");
      if (WmassFit)
      {
	 			Wmass = WmassFit->GetParameter(1); 
	 			SigmaWmass = WmassFit->GetParameter(2); 	
      }

      res->Close();
      delete res;
      //delete WmassFit;
    }


    ////////////////////////////////////
    ////////////////////////////////////
    ///////// Loop on events
    ////////////////////////////////////
    ////////////////////////////////////
      
    //block of code to arrange that only on 1/3th of the ttbar semimu or semiel sample is run for the training, and on 2/3th of the sample for the evaluation. 
    //Don't forget to change the eqLumi in the config depending on training or evaluation (not safe, to be changed? Not so trivial...)
    int start = 0;
    int end = nEvents;
    if(TrainMVA && (dataSetName == "TTbarJets_SemiMuon" || dataSetName == "TTbarJets_SemiElectron"))
    { 
        start = 0;
        end = int(nEvents/3);
    }
    else if (!TrainMVA && (dataSetName == "TTbarJets_SemiMuon" || dataSetName == "TTbarJets_SemiElectron"))
    {    
        start = int(nEvents/3);
        end = nEvents;
    }
    else{
        start = 0;
        end = nEvents;
    }
     
    cout << " - Loop over events " << endl;      
    //cout<<"start = "<<start<<"end = "<<end<<endl;
		
    for (int ievt = start; ievt < end; ievt++)
    {        

      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";
        
      //load event
			inInclFourthGenTree->GetEvent(ievt);
   
	 
	    vector<TLorentzVector> InitJets;
			vector<float> InitJetsbTagValues;
			if(btagger.find("TCHE")<=0)
			  InitJetsbTagValues = myBranch_selectedEvents->InitJetsbTagTCHE();
			else if(btagger.find("TCHP")<=0)
			  InitJetsbTagValues = myBranch_selectedEvents->InitJetsbTagTCHP(); 
			InitJets = myBranch_selectedEvents->InitJets();
	 
      // scale factor for the event
      float scaleFactor = 1.;

			double lumiWeight3D = 1.0;
			if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")){
				////////////////////////////
      	// apply trigger Reweighting
      	////////////////////////////
				float mceventtriggerweight = 1;
				float NEWmceventtriggerweight = 1;
				float prob=0;
      	if(semiElectron && doMCtriggerEffreweighting){
					std::vector<float> probabilities;
					for(size_t i=0; i<InitJets.size(); ++i){
    				if(fabs(InitJets[i].Eta())>2.6) continue;
    				probabilities.push_back(jetprob(InitJets[i].Pt(),InitJetsbTagValues[i]));
					}
				
					//use binary code for objects to be triggered or not triggered
					for(int i=0; i<pow(2.,(double)probabilities.size());++i){
    				int ntrigobj=0;
    				for(unsigned int j=0; j<probabilities.size();++j){
							if((int)(i/pow(2.,(double)j))%2) ntrigobj++;
						}
						if(ntrigobj<1) continue;  
						float newprob=1;
						for(unsigned int j=0; j<probabilities.size();++j){
							if((int)(i/pow(2.,(double)j))%2) newprob*=probabilities[j];
							else newprob*=1-probabilities[j];
						}
						prob+=newprob;
					}
					mceventtriggerweight*=prob;

					//stupid workaround, because part (single electron path corresponding to 200/pb) of the MC needs SF of 1 for trigger and part (btag trigger) needs the procedure above
					NEWmceventtriggerweight = (218./Luminosity +(mceventtriggerweight*(1-(218./Luminosity))));

					scaleFactor = scaleFactor*NEWmceventtriggerweight;
 					//cout << "mcevent triggerweight " << mceventtriggerweight << endl;
 					//cout << "scalefactor (only triggerweight) " << scaleFactor << endl;
      	}
        
      	////////////////////////////
      	// apply PU Reweighting
      	////////////////////////////
				lumiWeight3D = Lumi3DWeights.weight3D(myBranch_selectedEvents->nPUBXm1(),myBranch_selectedEvents->nPU(),myBranch_selectedEvents->nPUBXp1());
	 			if(doPUreweighting)
				{
					scaleFactor = scaleFactor*lumiWeight3D;
      		histo1D["lumiWeights"]->Fill(scaleFactor);
				}	
			}
			
			//reading variables from the tree
      bool isSingleLepton = myBranch_selectedEvents->SelectedSingleLepton();
      bool isSSLepton = myBranch_selectedEvents->SelectedSSLepton();
      bool isTriLepton = myBranch_selectedEvents->SelectedTriLepton();
      bool isSingleMuon = myBranch_selectedEvents->SelectedSingleMu();
      bool isSingleElectron = myBranch_selectedEvents->SelectedSingleEl();
      bool isSSMuon = myBranch_selectedEvents->SelectedSSMu();
      bool isSSElectron = myBranch_selectedEvents->SelectedSSEl();
      bool isSSMuEl = myBranch_selectedEvents->SelectedSSMuEl();
      bool isTriMuon = myBranch_selectedEvents->SelectedMuMuMu();
      bool isTriElectron = myBranch_selectedEvents->SelectedElElEl();
      bool isTriElMuMu = myBranch_selectedEvents->SelectedMuMuEl();
      bool isTriElElMu = myBranch_selectedEvents->SelectedMuElEl();
		  
			bool isSemiLep_MC = false;
			if(myBranch_selectedEvents->semiMuDecay() || myBranch_selectedEvents->semiElDecay())
			  isSemiLep_MC = true;
			
			float met = (myBranch_selectedEvents->met()).Et();
      vector<TLorentzVector> selectedJets = myBranch_selectedEvents->selectedJets();
			vector<TLorentzVector> selectedForwardJets = myBranch_selectedEvents->selectedForwardJets();
			vector<TLorentzVector> selectedMuons = myBranch_selectedEvents->selectedMuons();
			vector<TLorentzVector> selectedElectrons = myBranch_selectedEvents->selectedElectrons();
			vector<float> bTagValues, bTagValuesForMVA;
			bTagValuesForMVA = myBranch_selectedEvents->bTagTCHP(); //choice...
			if(btagger.find("TCHE")<=0)
			  bTagValues = myBranch_selectedEvents->bTagTCHE();
			else if(btagger.find("TCHP")<=0)
			  bTagValues = myBranch_selectedEvents->bTagTCHP();

						
			MSPlot["MS_nPV"]->Fill(myBranch_selectedEvents->nPV(),datasets[d], true, Luminosity*scaleFactor);				

			if(isSingleLepton){
				MSPlot["MS_MET"]->Fill(met,datasets[d], true, Luminosity*scaleFactor);				
				if(semiElectron) MSPlot["MS_LeptonPt"]->Fill(selectedElectrons[0].Pt(),datasets[d], true, Luminosity*scaleFactor);				
				if(semiMuon) MSPlot["MS_LeptonPt"]->Fill(selectedMuons[0].Pt(),datasets[d], true, Luminosity*scaleFactor);				
			
				for(unsigned int j=0;j<selectedJets.size();j++)
				{
				  MSPlot["MS_JetPt_all_SingleLepton"]->Fill(selectedJets[j].Pt(),datasets[d], true, Luminosity*scaleFactor);
					if(bTagValues[j] > workingpointvalue)
					  MSPlot["MS_JetPt_btagged_SingleLepton"]->Fill(selectedJets[j].Pt(),datasets[d], true, Luminosity*scaleFactor);
					else
					  MSPlot["MS_JetPt_nonbtagged_SingleLepton"]->Fill(selectedJets[j].Pt(),datasets[d], true, Luminosity*scaleFactor);
				}			
			}

			if(isSSLepton)
			{
				//cout << "IS SAME-SIGN LEPTON EVENT" << endl;
				NbSSevents = NbSSevents + datasets[d]->NormFactor()*Luminosity*scaleFactor;
			}
			if(isTriLepton)
			{
				//cout << "IS TRI-LEPTON EVENT" << endl;
				NbTrievents = NbTrievents + datasets[d]->NormFactor()*Luminosity*scaleFactor;
			}

			bool TprimeEvaluation = false; //temporarily not supported
			//////////////////////////////////////////////////////////////////////////
      // MVA training
      //////////////////////////////////////////////////////////////////////////
      if(doMVAjetcombination && selectedJets.size()>=4 && TrainMVA && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron))) //otherwise, if the jets vector only has 2 jets selected the jetcombiner crashes... 
      {
				if(!isSingleLepton) continue;				
        //the jets are sorted according to Pt in the TreeCreator, with everything conistently sorted along (should be), because the sorting was done before matching and pushing back b-tag values in vectors... Check and be careful!!
				//sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class 				
				if(semiMuon) 
					jetCombiner->ProcessEvent(datasets[d], myBranch_selectedEvents->mcQuarksForMatching(), selectedJets, bTagValuesForMVA, selectedMuons[0], isSemiLep_MC, scaleFactor, TprimeEvaluation);	        
				else if(semiElectron)
					jetCombiner->ProcessEvent(datasets[d], myBranch_selectedEvents->mcQuarksForMatching(), selectedJets, bTagValuesForMVA, selectedElectrons[0], isSemiLep_MC, scaleFactor, TprimeEvaluation);	        
      }

       
      if(TrainMVA) continue; //for the training, only the jetcombiner is relevant, so the following can be skipped (to the next event in the event loop)			
			//////////////////////////////////////////////////////////////////////////
      // the Wmassplot for the W counting
      //////////////////////////////////////////////////////////////////////////
      if(selectedJets.size()>=4 && !useMassesAndResolutions && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron)))
			{
      	if(myBranch_selectedEvents->Wbosonpartonsmatched())
	  		{ 					
					float WMassmatched_ = myBranch_selectedEvents->WMassmatched();
      		histo1D["hadronicRecoWMass"]->Fill(WMassmatched_,scaleFactor);
	  		}
			}
      
			

			/////////////////////////////////////////////////////////////////////////////
			//// core of the analysis, need to check b-tagging and reconstruct W bosons!	
			////////////////////////////////////////////////////////////////////////////
      if(useMassesAndResolutions)
      {
				vector<TLorentzVector> selectedJets_MVAinput; //for the MVA jet combination for the mass reconstruction in the 1B_2W and 2B_2W boxes
			  vector<float> bTagValuesForMVA_1B_2W, bTagValuesForMVA_2B_2W, bTagValuesForMVA_FromW, bTagValuesForMVA_FromW_DropUsedJets, bTagValuesForMVA_FromW_DropUsedJets_tmp;// the idea is to follow the exact same tricks as done with the jets, to keep track of the ordering
				//////////////////////////
				//b-tagging stuff here
				//////////////////////////
				vector<TLorentzVector> selectedJetsForBtagging; //need the jets within the tracker acceptance
				for(unsigned int i = 0; i<selectedJets.size(); i++)
				{
					selectedJetsForBtagging.push_back(selectedJets[i]);//this is the same collection as the selected jets, right?
				}
			
				vector< pair< int, float > > jetindex_btagvalue;
				vector< pair< int, bool > > jetindex_isb;
				for(unsigned int i = 0; i<selectedJetsForBtagging.size(); i++)
				{
					pair<int,float> dummy (i,bTagValues[i]);
					jetindex_btagvalue.push_back(dummy);
				}
				/*cout<<"BEFORE SORTING"<<endl;
				for(int k=0;k<jetindex_btagvalue.size();k++)
				{
					cout<<" jetindex_btagvalue["<<k<<"].first = "<<jetindex_btagvalue[k].first<<", jetindex_btagvalue["<<k<<"].second = "<<jetindex_btagvalue[k].second<<endl;
				}*/
 				
				// this is sorted according to the second value of the pair and not the first
 				std::sort(jetindex_btagvalue.begin(), jetindex_btagvalue.end(), sort_pair_decreasing());

				// check if jet is truly from a b-quark or not and make a pair with the result
				if(! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ))
				{
					for(unsigned int i = 0; i<jetindex_btagvalue.size(); i++)
					{
						if(fabs((myBranch_selectedEvents->partonFlavourJet())[jetindex_btagvalue[i].first]) == 5 )
						{// is a true b-jet
							pair<int,float> dummy (jetindex_btagvalue[i].first,true); //not: pair<int,float> dummy (i,true)
							jetindex_isb.push_back(dummy);
						}
						else
						{// is not a b-jet
							pair<int,float> dummy (jetindex_btagvalue[i].first,false); //not: pair<int,float> dummy (i,false)
							jetindex_isb.push_back(dummy);							
						}
					}
				}
							

				//////////////////
				// rescale events according to b-tag/mistag eff scalefactors
				//////////////////
      	double HT = 0.;
				HT = HT + met;

				int nbOfBtags = 0;
				int bJet1 = 999;
				int bJet2 = 999;
				for(unsigned int i = 0; i<jetindex_btagvalue.size(); i++)
				{
			  	//cout<<" jetindex_btagvalue["<<i<<"].first = "<<jetindex_btagvalue[i].first<<", jetindex_btagvalue["<<i<<"].second = "<<jetindex_btagvalue[i].second<<", selectedJetsForBtagging[jetindex_btagvalue["<<i<<"].first]->partonFlavour() = "<<selectedJetsForBtagging[jetindex_btagvalue[i].first]->partonFlavour()<<endl;
					if(jetindex_btagvalue[i].second > workingpointvalue)
					{
						nbOfBtags++;
						if(nbOfBtags<=2)
						{
					  	HT = HT + selectedJetsForBtagging[jetindex_btagvalue[i].first].Pt(); //only add pt of b-tagged jet for the 2 jets with the highest b-tag value (if more than 2 b-tagged jets, only the 2 highest are taken into account)
					  	selectedJets_MVAinput.push_back(selectedJetsForBtagging[jetindex_btagvalue[i].first]);
						}
						if(nbOfBtags == 1)
						{
							bJet1 = jetindex_btagvalue[i].first;
							bTagValuesForMVA_1B_2W.push_back(bTagValuesForMVA[bJet1]);
							bTagValuesForMVA_2B_2W.push_back(bTagValuesForMVA[bJet1]);
							for(unsigned int j = 0; j < jetindex_isb.size(); j++)
							{
								if(doBtagSFreweighting && bJet1 == jetindex_isb[j].first)
								{
									if(!jetindex_isb[j].second)
									{ // jet is b-tagged but NOT a true b
										scaleFactor = scaleFactor * mistagfactor;
									}
									else
									{ // jet is b-tagged and a true b
										scaleFactor = scaleFactor * scalefactorbtageff;
									}										
								}
							}
						}
						if(nbOfBtags == 2)
						{
							bJet2 = jetindex_btagvalue[i].first;
							bTagValuesForMVA_2B_2W.push_back(bTagValuesForMVA[bJet2]);
							for(unsigned int j = 0; j < jetindex_isb.size(); j++)
							{								
								if(doBtagSFreweighting &&bJet2 == jetindex_isb[j].first)
								{       
							  	//cout<<"    jetindex_isb["<<j<<"].second = "<<jetindex_isb[j].second<<endl;
									if(!jetindex_isb[j].second)
									{ // jet is b-tagged but NOT a true b
										//cout<<"jet is b-tagged but NOT a true b"<<endl;
										scaleFactor = scaleFactor * mistagfactor;
									} 
									else
									{ // jet is b-tagged and a true b
								  	//cout<<"jet is b-tagged and a true b"<<endl;
										scaleFactor = scaleFactor * scalefactorbtageff;
									}									
								}
							}
						}
					}
				} //end loop over jets (~ b-tagging)
  		
			if(isSingleLepton)
			{
				MSPlot["MS_METafterbtagSF"]->Fill(met,datasets[d], true, Luminosity*scaleFactor);
			}					
				
		
		  vector<TLorentzVector> selectedJetsFromW,selectedJetsFromW_DropUsedJets,selectedJetsFromW_DropUsedJets_tmp;
			///////////////////////
			//copy the collection of selected jets except for the one or two b-tagged jets (used in the HT calculation)
			//this collection is used for W counting				
			//////////////////////
			for(unsigned int i = 0; i<selectedJetsForBtagging.size(); i++)
			{
				if((int)i != bJet1 && (int)i != bJet2)
				{
					selectedJetsFromW.push_back(selectedJetsForBtagging[i]);
					bTagValuesForMVA_FromW.push_back(bTagValuesForMVA[i]);
				}
			}

	
			///////////////////////////////
			// start the W counting procedure here
			///////////////////////////////
			int nbOfWs = 0; 
			if(isSingleLepton) nbOfWs = 1;
			if(isSSLepton) nbOfWs = 2;
			if(isTriLepton) nbOfWs = 3;
			
			if(nbOfWs == 0){ cout << "WARNING - no W counted yet, this can NOT be!!!" << endl; break;}
			
			float recoWmass = 0.; float newrecoWmass = 0.;
			int indexWjet1 = 999;
			int indexWjet2 = 999;
			float massdifference = 100;
			int previoussize = 0;
			
			for(unsigned int i = 0; i<selectedJetsFromW.size(); i++)
			{
				selectedJetsFromW_DropUsedJets.push_back(selectedJetsFromW[i]);
				bTagValuesForMVA_FromW_DropUsedJets.push_back(bTagValuesForMVA_FromW[i]);
			}
			
			int firstWboson_jet1 = 9999; int firstWboson_jet2 = 9999;
			bool firsttime = true;
			do
			{
				previoussize = selectedJetsFromW_DropUsedJets.size();
			  indexWjet1 = 999; indexWjet2 = 999;
				massdifference = 100;
				for(unsigned int i = 0; i<selectedJetsFromW_DropUsedJets.size(); i++)
				{
					for(unsigned int j = 0; j<i; j++)
					{
						recoWmass = (selectedJetsFromW_DropUsedJets[i] + selectedJetsFromW_DropUsedJets[j]).M();
						if(!(dataSetName.find("NP_")<=dataSetName.size()))
 							MSPlot["allDiJetMasses"]->Fill(recoWmass, datasets[d], true, Luminosity*scaleFactor);
						if( fabs(recoWmass-Wmass)<massdifference)
						{
							massdifference = fabs(recoWmass-Wmass);
							newrecoWmass = recoWmass;
							indexWjet1 = i;
							indexWjet2 = j;
						}
					}
				}
				
				// fill a plot with the W mass from the chosen jets (only for the best jet couple, even if there are more W's reconstructed)
				if(firsttime && indexWjet1<999 && indexWjet2<999)
				{
					firsttime = false;
					firstWboson_jet1 = indexWjet1;
					firstWboson_jet2 = indexWjet2;
					float chosenWMass = (selectedJetsFromW_DropUsedJets[firstWboson_jet1] + selectedJetsFromW_DropUsedJets[firstWboson_jet2]).M();
					MSPlot["hadronicRecoWMass_chosenWjets"]->Fill(chosenWMass, datasets[d], true, Luminosity*scaleFactor);
				}

				if(massdifference < SigmaWmass && indexWjet1<999 && indexWjet2<999)
				{
					nbOfWs++;

					HT = HT + (selectedJetsFromW_DropUsedJets[indexWjet1] + selectedJetsFromW_DropUsedJets[indexWjet2]).Pt();
					
					if(selectedJets_MVAinput.size()<=2) //you're only using it for the nbOfWs=2 boxes anyway...
					{
					        selectedJets_MVAinput.push_back(selectedJetsFromW_DropUsedJets[indexWjet1]);
					        selectedJets_MVAinput.push_back(selectedJetsFromW_DropUsedJets[indexWjet2]);
									bTagValuesForMVA_1B_2W.push_back(bTagValuesForMVA_FromW_DropUsedJets[indexWjet1]);
									bTagValuesForMVA_1B_2W.push_back(bTagValuesForMVA_FromW_DropUsedJets[indexWjet2]);
									bTagValuesForMVA_2B_2W.push_back(bTagValuesForMVA_FromW_DropUsedJets[indexWjet1]);
									bTagValuesForMVA_2B_2W.push_back(bTagValuesForMVA_FromW_DropUsedJets[indexWjet2]);
					}
					
					for(unsigned int i = 0; i<selectedJetsFromW_DropUsedJets.size(); i++)
					{
						if((unsigned int)indexWjet1!=i && (unsigned int)indexWjet2!=i)
						{
							selectedJetsFromW_DropUsedJets_tmp.push_back(selectedJetsFromW_DropUsedJets[i]);
							bTagValuesForMVA_FromW_DropUsedJets_tmp.push_back(bTagValuesForMVA_FromW_DropUsedJets[i]);
						}
					}
					selectedJetsFromW_DropUsedJets.clear();
					bTagValuesForMVA_FromW_DropUsedJets.clear();
					for(unsigned int i = 0; i<selectedJetsFromW_DropUsedJets_tmp.size(); i++)
					{
						selectedJetsFromW_DropUsedJets.push_back(selectedJetsFromW_DropUsedJets_tmp[i]);
						bTagValuesForMVA_FromW_DropUsedJets.push_back(bTagValuesForMVA_FromW_DropUsedJets_tmp[i]);
					}
					selectedJetsFromW_DropUsedJets_tmp.clear();
					bTagValuesForMVA_FromW_DropUsedJets_tmp.clear();
				}

			}
			while(selectedJetsFromW_DropUsedJets.size()>=2 && selectedJetsFromW_DropUsedJets.size()!= (unsigned int)previoussize);

				
			if(nbOfBtags>=2) nbOfBtags = 2; //to make sure that events with more than 2 b-jets end up in the 2 b-jet bin			
			
			////////////////////////////////////
			//// here we put the events into boxes of #b's and #W's
			////////////////////////////////////
			if(nbOfWs>=4) nbOfWs = 4; // to make sure that events with more than 4 W bosons end up in the 4W bin

			selecTableSemiLep.Fill(d,0,scaleFactor);
	
			
			//cout<<"nbOfBtags = "<<nbOfBtags<<endl;						
			if(nbOfBtags==1)
			{				
				if(nbOfWs==1)
				{
					//cout << "in 1B 1W box" << endl;
					//Requirements: a leptonically decaying W, and a forward jet and a 'central' one
					if(selectedForwardJets.size() != 1) continue;
					if(selectedJets.size() != 1) continue;
				
					float DeltaPhi = fabs(selectedJets[0].DeltaPhi(selectedForwardJets[0]));
					float RelPt = fabs(selectedJets[0].Pt()-selectedForwardJets[0].Pt())/(selectedJets[0].Pt()+selectedForwardJets[0].Pt());
				
					double const Pi=4*atan(1);
					if(DeltaPhi > ((Pi/2)+RelPt*Pi)) continue;
					
					if(isSingleLepton) 
					{						
						HT = HT + selectedForwardJets[0].Pt();

						myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,met,selectedJets,scaleFactor);
						selecTableSemiLep.Fill(d,1,scaleFactor);
					}
					//cout << "done in 1B 1W box" << endl;
				}
				else if(nbOfWs==2)
				{
					//cout << "in 1B 2W box" << endl;
					if(isSingleLepton && selectedJets.size()>=4)
					{
						HT = HT + selectedJetsFromW_DropUsedJets[0].Pt();

						myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,met,selectedJets,scaleFactor);
						selecTableSemiLep.Fill(d,2,scaleFactor);
			
						//////////////
						// MVA stuff, only to be done when this is decided in the configuration of this macro (the doMVAjetcombination boolean)
						//////////////
						pair<float, vector<unsigned int> > MVAvals;	
						if(doMVAjetcombination)
						{	
					   //you should put the highest Pt jet in here THAT IS NOT YET IN THE VECTOR!
       			 if(selectedJets_MVAinput.size()==3)
					   {
					     selectedJets_MVAinput.push_back(selectedJetsFromW_DropUsedJets[0]);
							 bTagValuesForMVA_1B_2W.push_back(bTagValuesForMVA_FromW_DropUsedJets[0]);
					   }
					   else
					     cout<<"Something is not right, the vector of selected jets for MVA input in this stage should be 3"<<endl;
					   
					   if(selectedJets_MVAinput.size()==4)
					   {	
					
					     //sort(selectedJets_MVAinput.begin(),selectedJets_MVAinput.end(),HighestPt()); // HighestPt() is included from the Selection class
 	                          
			   			 if(semiMuon) jetCombiner->ProcessEvent(datasets[d], myBranch_selectedEvents->mcQuarksForMatching(), selectedJets_MVAinput, bTagValuesForMVA_1B_2W, selectedMuons[0], isSemiLep_MC, scaleFactor, TprimeEvaluation); //datasets[d],mcParticles,selectedJets_MVAinput,selectedMuons[0],init_electrons,init_muons,genEvt,scaleFactor);	//OLD WAY (class has changed since then): jetCombiner->ProcessEvent(datasets[d], mcParticles, selectedJets, selectedMuons[0], vertex[0], eventSelected, init_electrons, init_muons, scaleFactor);
               else if(semiElectron) jetCombiner->ProcessEvent(datasets[d], myBranch_selectedEvents->mcQuarksForMatching(), selectedJets_MVAinput, bTagValuesForMVA_1B_2W, selectedElectrons[0], isSemiLep_MC, scaleFactor, TprimeEvaluation);
	  
			   			 //vector<unsigned int> goodCombi = jetCombiner_1B_2W->GetGoodJetCombination(); //get the MC matched jet combination, not the MVA best matched		   	
			   			 MVAvals = jetCombiner->getMVAValue(MVAmethod, 1); // 1 means the highest MVA value
					   }
					   else
					     cout<<"WARNING: vector of selected jets for MVA input is not equal to 4 (but to "<<selectedJets_MVAinput.size()<<"); fix this!!"<<endl;					   
					
					   //reconstructing the mass of the top/t', and fill vector for binning
					   myInclFourthGenSearchTools.CalculateTopMass(selectedJets_MVAinput[MVAvals.second[0]],selectedJets_MVAinput[MVAvals.second[1]],selectedJets_MVAinput[MVAvals.second[2]]); //(light 1, light 2, hadronic b-jet);
					   myInclFourthGenSearchTools.FillMassPlots(d,nbOfBtags,nbOfWs,scaleFactor);
					   	
					   //for the 2D binning (HT vs Mtop)
					   float fillweight = datasets[d]->NormFactor () * Luminosity * scaleFactor;
					   if(anaEnv.MCRound==0)
					   {
					        HTvsMTop_1B_2W.Fill_for2DBinning(myInclFourthGenSearchTools.GetHT(),myInclFourthGenSearchTools.GetMtop(),fillweight);
					   }
					   else if(anaEnv.MCRound!=0)
					   {
					        HTvsMTop_1B_2W.Fill(myInclFourthGenSearchTools.GetHT(),myInclFourthGenSearchTools.GetMtop(),fillweight,d);
					   }
					 }
					
				    } 
				    else if(isSSLepton && selectedJets.size() >=2) 
				    {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							selecTableMultiLep.Fill(d,0,scaleFactor);
							
							if(isSSMuon){
								histoName = "NbEvents_"+multileptons[2];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSMuMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSElectron){
								histoName = "NbEvents_"+multileptons[0];
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElElevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSMuEl){
								histoName = "NbEvents_"+multileptons[1];
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);							
							}	
				    }
 	
					//cout << "done in 1B 2W box" << endl;
				}
				else if(nbOfWs==3)
				{
				    //cout << "in 1B 3W box" << endl;
				    if(isSingleLepton && selectedJets.size() >=6) 
				    {
				      myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,met,selectedJets,scaleFactor);				
							selecTableSemiLep.Fill(d,3,scaleFactor);
				    }
				    else if(isSSLepton && selectedJets.size() >=4) 
			      {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
							selecTableSemiLep.Fill(d,0,scaleFactor);
							if(isSSMuon){
								histoName = "NbEvents_"+multileptons[2];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSMuMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSElectron){
								histoName = "NbEvents_"+multileptons[0];
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElElevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSMuEl){
								histoName = "NbEvents_"+multileptons[1];
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);							
							}	
						}
				    else if(isTriLepton && selectedJets.size() >=2) 
				    {
							//cout << "IS TRI-LEPTON" << endl;
							//count number of events with three leptons
							MSPlot["MS_NbTrievents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
				    	selecTableMultiLep.Fill(d,1,scaleFactor);
							if(isTriElectron) histoName = "NbEvents_"+multileptons[3];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriMuon) histoName = "NbEvents_"+multileptons[6];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriElElMu) histoName = "NbEvents_"+multileptons[4];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriElMuMu) histoName = "NbEvents_"+multileptons[5];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							histo_dataset = histoName+dataSetName.c_str(); 
							histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
						}
				    //cout << "done in 1B 3W box" << endl;
				}
				else
				{
					//cout << "in 1B 4W box" << endl;
				   if(isSingleLepton && selectedJets.size() >=8) 
				   {
							myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,met,selectedJets,scaleFactor);
							selecTableSemiLep.Fill(d,4,scaleFactor);
					
				   }
				   else if(isSSLepton && selectedJets.size() >=6) 
				   {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
			     		selecTableMultiLep.Fill(d,0,scaleFactor);
							if(isSSMuon){
								histoName = "NbEvents_"+multileptons[2];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSMuMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSElectron){
								histoName = "NbEvents_"+multileptons[0];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElElevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSMuEl){
								histoName = "NbEvents_"+multileptons[1];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);							
							}	
					 }
				   else if(isTriLepton && selectedJets.size() >=4) 
				   {
							//cout << "IS TRI-LEPTON" << endl;
							//count number of events with three leptons
							MSPlot["MS_NbTrievents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
				      selecTableMultiLep.Fill(d,1,scaleFactor);
							if(isTriElectron) histoName = "NbEvents_"+multileptons[3];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriMuon) histoName = "NbEvents_"+multileptons[6];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriElElMu) histoName = "NbEvents_"+multileptons[4];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriElMuMu) histoName = "NbEvents_"+multileptons[5];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							histo_dataset = histoName+dataSetName.c_str(); 
							histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
					  }
				   //cout << "done in 1B 4W box" << endl;
				}
				
			} //end number of btags == 1
			else if(nbOfBtags == 2)
			{
				if(nbOfWs>=4) nbOfWs = 4;// to make sure that events with more than 4 W bosons end up in the 4W bin
				
				if(nbOfWs==1)
				{
				   //cout << "in 2B 1W box" << endl;
				   if(selectedForwardJets.size() != 0) continue;
				   if(selectedJets.size() != 2) continue;

				   float DeltaPhi = fabs(selectedJets[0].DeltaPhi(selectedJets[1]));
				   float RelPt = fabs(selectedJets[0].Pt()-selectedJets[1].Pt())/(selectedJets[0].Pt()+selectedJets[1].Pt());

				   double const Pi=4*atan(1);
				   if(DeltaPhi > ((Pi/2)+RelPt*Pi)) continue;
				   
				   //cout << "in 2B 1W box" << endl;

				   if(isSingleLepton) 
				   {
				      myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,met,selectedJets,scaleFactor);
							selecTableSemiLep.Fill(d,5,scaleFactor);
				   }
					//cout << "done in 2B 1W box" << endl;
				}
				else if(nbOfWs==2)
				{
					//cout << "in 2B 2W box" << endl;
				   if(isSingleLepton && selectedJets.size() >=4) 
				   {
					
							myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,met,selectedJets,scaleFactor);
							selecTableSemiLep.Fill(d,6,scaleFactor);
					
						//////////////
						// MVA stuff, only to be done when this is decided in the configuration of this macro (the doMVAjetcombination boolean)
						//////////////
					pair<float, vector<unsigned int> > MVAvals;
					//bool TprimeEvaluation = false;
					//if(dataSetName.find("NP_Tprime")==0 || dataSetName.find("NP_overlay_Tprime")==0)
					//  TprimeEvaluation = true;
					if(doMVAjetcombination)
					{	
					   if(selectedJets_MVAinput.size()==4)
					   {
       			   			//sort(selectedJets_MVAinput.begin(),selectedJets_MVAinput.end(),HighestPt()); // HighestPt() is included from the Selection class
	                          			
			  						if(semiMuon) jetCombiner->ProcessEvent(datasets[d], myBranch_selectedEvents->mcQuarksForMatching(), selectedJets_MVAinput, bTagValuesForMVA_2B_2W, selectedMuons[0], isSemiLep_MC, scaleFactor, TprimeEvaluation); //datasets[d],mcParticles,selectedJets_MVAinput,selectedMuons[0],init_electrons,init_muons,genEvt,scaleFactor);	//OLD WAY (class has changed since then): jetCombiner->ProcessEvent(datasets[d], mcParticles, selectedJets, selectedMuons[0], vertex[0], eventSelected, init_electrons, init_muons, scaleFactor);
              			else if(semiElectron) jetCombiner->ProcessEvent(datasets[d], myBranch_selectedEvents->mcQuarksForMatching(), selectedJets_MVAinput, bTagValuesForMVA_2B_2W, selectedElectrons[0], isSemiLep_MC, scaleFactor, TprimeEvaluation);
		  
			   					//vector<unsigned int> goodCombi = jetCombiner_2B_2W->GetGoodJetCombination(); //get the MC matched jet combination, not the MVA best matched		   	
			   					MVAvals = jetCombiner->getMVAValue(MVAmethod, 1); // 1 means the highest MVA value

	   					//coutObjectsFourVector(init_muons,init_electrons,selectedJets_MVAinput,mets,"*** Before kinematic fit ***");
						
						////to test 'purity' of good jet combinations with a 'simple' method (not MVA)... for the moment only testable in 2B_2W box
					        //myInclFourthGenSearchTools.TestPurityGoodCombinations(d,nbOfBtags,nbOfWs,genEvt,mcParticles,selectedJets_MVAinput,TprimeEvaluation,scaleFactor);
					        
						
					        //reconstructing the mass of the top/t', and fill vector for binning
					   			myInclFourthGenSearchTools.CalculateTopMass(selectedJets_MVAinput[MVAvals.second[0]],selectedJets_MVAinput[MVAvals.second[1]],selectedJets_MVAinput[MVAvals.second[2]]); //(light 1, light 2, hadronic b-jet);
					        myInclFourthGenSearchTools.FillMassPlots(d,nbOfBtags,nbOfWs,scaleFactor);
					   	
						//for the 2D binning (HT vs Mtop)
					   	float fillweight = datasets[d]->NormFactor () * Luminosity * scaleFactor;
					   	if(anaEnv.MCRound==0)
					   	{
					   	  HTvsMTop_2B_2W.Fill_for2DBinning(myInclFourthGenSearchTools.GetHT(),myInclFourthGenSearchTools.GetMtop(),fillweight);
					   	}
					   	else if(anaEnv.MCRound!=0)
					   	{
					   	  HTvsMTop_2B_2W.Fill(myInclFourthGenSearchTools.GetHT(),myInclFourthGenSearchTools.GetMtop(),fillweight,d);
					   	}
					   }
					   else
					     cout<<"WARNING: vector of selected jets for MVA input is not equal to 4; fix this!!"<<endl;	
					}
			             } 
				     else if(isSSLepton && selectedJets.size() >=2) 
				     {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
						  selecTableMultiLep.Fill(d,0,scaleFactor);
							if(isSSMuon){
								histoName = "NbEvents_"+multileptons[2];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSMuMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSElectron){
								histoName = "NbEvents_"+multileptons[0];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElElevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSMuEl){
								histoName = "NbEvents_"+multileptons[1];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);							
							}	
					  } 
				    // cout << "done in 2B 2W box" << endl;
				}
				else if(nbOfWs==3)
				{    
				    // cout << "in 2B 3W box" << endl;
				     if(isSingleLepton && selectedJets.size() >=6) 
				     {
					myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,met,selectedJets,scaleFactor);
							selecTableSemiLep.Fill(d,7,scaleFactor);
				     }
				     else if(isSSLepton && selectedJets.size() >=4) 
				     {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
						  selecTableMultiLep.Fill(d,0,scaleFactor);
							if(isSSMuon){
								histoName = "NbEvents_"+multileptons[2];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSMuMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSElectron){
								histoName = "NbEvents_"+multileptons[0];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElElevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSMuEl){
								histoName = "NbEvents_"+multileptons[1];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);							
							}	
						 }
				     else if(isTriLepton && selectedJets.size() >=2) 
				     {
							//cout << "IS TRI-LEPTON" << endl;
							//count number of events with three leptons
							MSPlot["MS_NbTrievents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
						  selecTableMultiLep.Fill(d,1,scaleFactor);
							if(isTriElectron) histoName = "NbEvents_"+multileptons[3];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriMuon) histoName = "NbEvents_"+multileptons[6];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriElElMu) histoName = "NbEvents_"+multileptons[4];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriElMuMu) histoName = "NbEvents_"+multileptons[5];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							histo_dataset = histoName+dataSetName.c_str(); 
							histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
						 }
					//cout << "done in 2B 3W box" << endl;
				}
				else
				{
				   //  cout << "in 2B 4W box" << endl;
				     if(isSingleLepton && selectedJets.size() >=8 )
				     {
							myInclFourthGenSearchTools.FillPlots(d,nbOfBtags,nbOfWs,HT,selectedMuons,selectedElectrons,met,selectedJets,scaleFactor);
							selecTableSemiLep.Fill(d,8,scaleFactor);
				     }
				     else if(isSSLepton && selectedJets.size() >=6) 
				     {
							//cout << "IS SAME-SIGN LEPTON" << endl;
							//count number of events with SS leptons
							MSPlot["MS_NbSSevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
						  selecTableMultiLep.Fill(d,0,scaleFactor);
							if(isSSMuon){
								histoName = "NbEvents_"+multileptons[2];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSMuMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSElectron){
								histoName = "NbEvents_"+multileptons[0];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElElevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);	
							}else if(isSSMuEl){
								histoName = "NbEvents_"+multileptons[1];//({"SSElEl","SSElMu","SSMuMu","TriLeptons"})
								histo_dataset = histoName+dataSetName.c_str(); 
								histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
								MSPlot["MS_NbSSElMuevents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);							
							}	
						 }
				     else if(isTriLepton && selectedJets.size() >=4) 
				     {
						//cout << "IS TRI-LEPTON" << endl;
								//count number of events with three leptons
							MSPlot["MS_NbTrievents"]->Fill(1.0,datasets[d], true, Luminosity*scaleFactor);
						  selecTableMultiLep.Fill(d,1,scaleFactor);
							if(isTriElectron) histoName = "NbEvents_"+multileptons[3];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriMuon) histoName = "NbEvents_"+multileptons[6];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriElElMu) histoName = "NbEvents_"+multileptons[4];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							if(isTriElMuMu) histoName = "NbEvents_"+multileptons[5];//({"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"})
							histo_dataset = histoName+dataSetName.c_str(); 
							histo1D[histo_dataset.c_str()]-> Fill(1.0,datasets[d]->NormFactor()*Luminosity*scaleFactor);
						}
				  //   cout << "done in 2B 4W box" << endl;
				}
			} //end number of btags == 2
		
		 }//useMassesAndResolutions = true
	//delete selection;
   // cout << "done processing event" << endl;
		}//loop on events

    ////to test 'purity' of good jet combinations with a 'simple' method (not MVA). For the moment not supported
    //myInclFourthGenSearchTools.PrintPurityGoodCombinations();
    
    cout<<endl;
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
		if (!useMassesAndResolutions && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) ||  (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron)))
    {
	
      string resfilename = "MassPlot.root";
     
      cout << "INFO: Creating output file to store the mass of the W boson: " << resfilename << endl;    

      TFile* outfile = new TFile(resfilename.c_str(),"RECREATE");
      outfile->cd();
      
      histo1D["hadronicRecoWMass"]->Write();

      // fit the distributions
      TH1F* histo = histo1D["hadronicRecoWMass"];
      TF1 *fitfunc;
      string func_title = string(histo->GetName())+"_Fitted";

      double rms = histo->GetRMS();
      double maxbin =  histo->GetBinCenter(histo->GetMaximumBin());

      fitfunc = new TF1(func_title.c_str(),"gaus");
      fitfunc->SetRange(maxbin-rms,maxbin+rms);
      histo->Fit(fitfunc,"RQ");

      fitfunc->Write();
      delete fitfunc;
      
      outfile->Close();
     //delete outfile; 
    }
		
		inFile->Close();
		delete inFile;
		
  } //loop on 'datasets'

  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;
  fout->cd();
	
  if(doMVAjetcombination)
  {
	  string pathPNGJetCombi = pathPNG+"JetCombination/";
    if(savePNG)
		{
      mkdir(pathPNGJetCombi.c_str(),0777);
		}
		jetCombiner->Write(fout, savePNG, pathPNGJetCombi);
  }  

  
  if(!TrainMVA)
  {
    fout->cd();
    //Write histograms: MSPlots
    if(savePNG) mkdir((pathPNG+"MSPlot/").c_str(),0777);
		cout << "Running over all MS plots" << endl;
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        temp->Draw(false, name, true, true, true, true, true,5,false, true, true);//(bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal, bool addRatio, bool mergeVV, bool mergeTTV)
        temp->Write(fout, name, savePNG, pathPNG+"MSPlot/");//bool savePNG
    }
    cout << "MultiSamplePlots written" << endl;
	
    //Write histograms: 1D 
    TDirectory* th1dir = fout->mkdir("1D_histograms");
    fout->cd();
    th1dir->cd();
    for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
			TH1F *temp = it->second;
 	//		int N = temp->GetNbinsX();
 	//  	temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
 	//  	temp->SetBinContent(N+1,0);
 	//		temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
			temp->Write();
			TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
			if(savePNG) tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
    }    
    cout << "1D plots written" << endl;
    
    myInclFourthGenSearchTools.WritePlots(fout,th1dir,savePNG,pathPNG);
    
    //delete th1dir;
    // 2D
    TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
    th2dir->cd();
    for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {
			TH2F *temp = it->second;
			temp->Write();
			TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
			if(savePNG) tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
    }
    //delete th2dir;
    cout << "2D plots written" << endl;
    fout->cd();
    
    //Selection tables
    selecTableSemiLep.TableCalculator(true, true, true, true, true, true, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
    string selectiontableSemiLep = "InclFourthGenSearch_SelectionTable_TreeAnalyzer_"+postfix+channelpostfix;
    selectiontableSemiLep = selectiontableSemiLep +".tex"; 	
    selecTableSemiLep.Write(selectiontableSemiLep.c_str(),false, true, false, false, false, false, false); //(filename, error, merged, lines, unscaled, eff, totaleff, landscape)
    
		selecTableMultiLep.TableCalculator(true, true, true, true, true, true, true, true);
    string selectiontableMultiLep = "InclFourthGenSearch_SelectionTable_TreeAnalyzer_MultiLepton_"+postfix+channelpostfix;
    selectiontableMultiLep = selectiontableMultiLep +".tex"; 	
    selecTableMultiLep.Write(selectiontableMultiLep.c_str(),false, true, false, false, false, false, false);
		
     //regarding binning, which is only relevant when you do the MVA jet combination to reconstruct the top mass   
     if (doMVAjetcombination && anaEnv.MCRound==0)
     {	
        HTvsMTop_1B_2W.Write_for2DBinning(binningFileName_HTvsMTop_1B_2W);
				HTvsMTop_2B_2W.Write_for2DBinning(binningFileName_HTvsMTop_2B_2W);
     }    
     else if(doMVAjetcombination && anaEnv.MCRound==1)
     {  
				HTvsMTop_1B_2W.Convert2Dto1D(postfix);
        HTvsMTop_2B_2W.Convert2Dto1D(postfix);
        HTvsMTop_1B_2W.Write(fout,th1dir);
				HTvsMTop_2B_2W.Write(fout,th1dir);
     }
     
    fout->cd();
    
    cout << " - Closing the output file now..." << endl;
    fout->Close();
  } //end !trainMVA
  
  //delete
  if(jetCombiner) delete jetCombiner; //IMPORTANT!! file for training otherwise not filled... (?) //crashes when calculating resolutions for kinfit
  
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

//To cout the Px, Py, Pz, E and Pt of objects
void coutObjectsFourVector(vector < TRootMuon* > init_muons, vector < TRootElectron* > init_electrons, vector < TRootJet* > init_jets, vector < TRootMET* > mets, string Comment)
{
     cout<<Comment<<endl;
     
     for(unsigned int k=0; k<init_muons.size(); k++)
     {
	   cout<<" init_muons["<<k<<"] -> Px() = "<<init_muons[k]->Px()<<endl;
	   cout<<"              -> Py() = "<<init_muons[k]->Py()<<endl;
	   cout<<"              -> Pz() = "<<init_muons[k]->Pz()<<endl;
	   cout<<"                -> Pt() = "<<init_muons[k]->Pt()<<endl;
	   cout<<"              -> E() = "<<init_muons[k]->E()<<endl;   
     }
     for(unsigned int k=0; k<init_electrons.size(); k++)
     {
	   cout<<" init_electrons["<<k<<"] -> Px() = "<<init_electrons[k]->Px()<<endl;
	   cout<<"              -> Py() = "<<init_electrons[k]->Py()<<endl;
	   cout<<"              -> Pz() = "<<init_electrons[k]->Pz()<<endl;
	   cout<<"                -> Pt() = "<<init_electrons[k]->Pt()<<endl;
	   cout<<"              -> E() = "<<init_electrons[k]->E()<<endl;   
     }         
     for(unsigned int k=0; k<init_jets.size(); k++)
     {
	   cout<<" init_jets["<<k<<"] -> Px() = "<<init_jets[k]->Px()<<endl;
	   cout<<"              -> Py() = "<<init_jets[k]->Py()<<endl;
	   cout<<"              -> Pz() = "<<init_jets[k]->Pz()<<endl;
	   cout<<"                -> Pt() = "<<init_jets[k]->Pt()<<endl;
	   cout<<"              -> E() = "<<init_jets[k]->E()<<endl;	   
     }
     for(unsigned int k=0; k<mets.size(); k++)
     {
           cout<<" mets["<<k<<"] -> Px() = "<<mets[k]->Px()<<endl;
           cout<<"         ->  Py() = "<<mets[k]->Py()<<endl;
	   cout<<"         ->  Pz() = "<<mets[k]->Pz()<<endl;
	   cout<<"              -> Pt() = "<<mets[k]->Pt()<<endl;
	   cout<<"         ->  E() = "<<mets[k]->E()<<endl;
	   cout<<"              -> Et() = "<<mets[k]->Et()<<endl;
     }
};

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopTurnOnCurves
float jetprob(float jetpt, float btagvalue){
	float prob=0.982*exp(-30.6*exp(-0.151*jetpt));
  prob*=0.844*exp((-6.72*exp(-0.720*btagvalue))); //"for the offline TCHP tagger"
	//prob*=0.736*exp((-8.01*exp(-0.540*btagvalue))); //"for the offline TCHE tagger"
	return prob;
};
