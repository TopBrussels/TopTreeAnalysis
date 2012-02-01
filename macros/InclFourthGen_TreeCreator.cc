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
  if( ! (systematic == "Nominal"  || systematic == "JESPlus" || systematic == "JESMinus" || systematic == "JERPlus" || systematic == "JERMinus") )
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
		mistagfactor = 1.11;		
  }
  else if(btagger == "TCHEM") //track counting high eff medium working point
  {
		scalefactorbtageff = 0.94;
		mistagfactor = 1.21;
  }
	else if(btagger == "TCHPM") //track counting high pur medium working point
  {
		scalefactorbtageff = 0.91;		
		mistagfactor = 1.27;
  }
	else if(btagger == "TCHPT") //track counting high pur tight working point
  {
	  cout<<"WARNING: look up SFs for TCHPT"<<endl;
  	 /* if(dobTagEffShift == 0)
		scalefactorbtageff = 0.91;
	  if(dobTagEffShift == 1)
		scalefactorbtageff = 0.8194;
	  if(dobTagEffShift == 2)
		scalefactorbtageff = 1.0006;
		
	  if(domisTagEffShift == 0)
		mistagfactor = 1.27;
	  if(domisTagEffShift == 1)
		mistagfactor = 1.1187;
	  if(domisTagEffShift == 2)
		mistagfactor = 1.4213;*/
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

  string postfix = "_TEST"; // to relabel the names of the output file  
	postfix= postfix+"_"+systematic;

  string Treespath = "InclFourthGenTrees";
  Treespath = Treespath +"/";
  mkdir(Treespath.c_str(),0777);
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
  
  //bool TrainwithTprime = false;
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
  if(semiElectron) xmlFileName = "../config/myFourthGenconfig_Electron.xml";
  else if(semiMuon) xmlFileName = "../config/myFourthGenconfig.xml";
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
  
  //Output ROOT file
  string rootFileName ("InclFourthGenSearch_TreeCreator"+postfix+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

 
  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootMET* > mets;

  //Global variable
  TRootEvent* event = 0;

  string pathPNG = "InclFourthGenSearchPlots_TreeCreator"+postfix+channelpostfix;
  pathPNG = pathPNG +"/"; 	
  pathPNG = pathPNG +"/"; 	
  if(savePNG) mkdir(pathPNG.c_str(),0777);


  MSPlot["allDiJetMasses"] = new MultiSamplePlot(datasets, "allDiJetMasses", 50, 0, 1000, "m_{jj}");
  MSPlot["hadronicRecoWMass_chosenWjets"] = new MultiSamplePlot(datasets, "hadronicRecoWMass_chosenWjets", 50, 0, 500, "m_{W}"); 
  histo1D["hadronicPartonWMass"] = new TH1F("hadronicPartonWMass","Hadronic W Mass, using the Partons",100,0,200);
  histo1D["hadronicRecoWMass"] = new TH1F("hadronicRecoWMass","Hadronic W Mass, using the RecoJets",100,0,200);
  
  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
  histo1D["LeptonPt_TTbar"] = new TH1F("leptonspt ttbar","leptonspt ttbar;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_Tprime500"] = new TH1F("leptonspt tprime500","leptonspt tprime500;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_Bprime500"] = new TH1F("leptonspt bprime500","leptonspt bprime500;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_SBprime500"] = new TH1F("leptonspt sbprime500","leptonspt sbprime500;pt leptons;#events",250,0,500);
  
  string multileptons[2] = {"SSLeptons","TriLeptons"};
  string histoName,histo_dataset;
  for(int i = 0; i<2; i++)
  {
		histoName = "NbEvents_"+multileptons[i];
		for(unsigned int d = 0; d < datasets.size (); d++){
			histo_dataset = histoName+(datasets[d]->Name()).c_str(); 
			histo1D[histo_dataset.c_str()] = new TH1F(histo_dataset.c_str(),histo_dataset.c_str(), 1, 0.5, 1.5);
		}
  }
	
  MSPlot["MS_NbSSevents"] = new MultiSamplePlot(datasets,"# events with SS leptons", 1, 0.5, 1.5, "");
  MSPlot["MS_NbTrievents"] = new MultiSamplePlot(datasets,"# events with 3 leptons", 1, 0.5, 1.5, "");
  MSPlot["MS_MET"] = new MultiSamplePlot(datasets,"MET", 75, 0, 150, "");
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
  
	
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("initial")); //0
  CutsSelecTableSemiMu.push_back(string("preselected"));
  CutsSelecTableSemiMu.push_back(string("trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 1 muon"));
  CutsSelecTableSemiMu.push_back(string("$\\geq$ 1 b-tagged jet"));
  CutsSelecTableSemiMu.push_back(string("MET $>$ 40 GeV"));  
  CutsSelecTableSemiMu.push_back(string("single muon"));

  vector<string> CutsSelecTableSemiEl;
  CutsSelecTableSemiEl.push_back(string("initial")); //0
  CutsSelecTableSemiEl.push_back(string("preselected"));
  CutsSelecTableSemiEl.push_back(string("trigged"));
  CutsSelecTableSemiEl.push_back(string("Good PV"));
  CutsSelecTableSemiEl.push_back(string("$\\geq$ 1 selected electron"));
  CutsSelecTableSemiEl.push_back(string("Conversion veto"));
  CutsSelecTableSemiEl.push_back(string("$\\geq$ 1 b-tagged jet"));
  CutsSelecTableSemiEl.push_back(string("MET $>$ 40 GeV"));
  CutsSelecTableSemiEl.push_back(string("single electron"));

  vector<string> CutsSelecTableMultiLepton;
  CutsSelecTableMultiLepton.push_back(string("SS leptons"));
  CutsSelecTableMultiLepton.push_back(string("trileptons"));

  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  selecTableSemiMu.SetLuminosity(Luminosity);
  selecTableSemiMu.SetPrecision(1);
  SelectionTable selecTableSemiEl(CutsSelecTableSemiEl, datasets);
  selecTableSemiEl.SetLuminosity(Luminosity);
  selecTableSemiEl.SetPrecision(1);
  SelectionTable selecTableMultiLepton(CutsSelecTableMultiLepton, datasets);
  selecTableMultiLepton.SetLuminosity(Luminosity);
  selecTableMultiLepton.SetPrecision(1);
  
  cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////////////////////
  // PileUp Reweighting - 3D//
  ////////////////////////////////////////////////////
  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DWeights.weight3D_init(1.0);

//  LumiReWeighting LumiWeights = LumiReWeighting("PileUpReweighting/pileup_WJets_36bins.root", "PileUpReweighting/pileup_2011Data_UpToRun180252.root", "pileup2", "pileup");
//  PoissonMeanShifter PShiftUp = PoissonMeanShifter(0.6); // PU-systematic
//  PoissonMeanShifter PShiftDown = PoissonMeanShifter(-0.6); // PU-systematic
  cout << " - Initialized LumiReWeighting stuff" << endl;
  
  

/*	
  ///////////////
  // JetCombiner, to test
  ///////////////
  JetCombiner* jetCombiner;
  if(!doMVAjetcombination) TrainMVA = false;
  else if(doMVAjetcombination)
  {
    jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, true, "",channelpostfix); //last bool is basically to use also the W mass as constraint    
  }
    
  cout << " - JetCombiner instantiated ..." << endl;
*/  
    
	
		
  ////////////////////////////////////
  ////////////////////////////////////
  ///////// Loop on datasets
  ////////////////////////////////////
  ////////////////////////////////////
  cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

  for (unsigned int d = 0; d < datasets.size(); d++) //d < datasets.size()
  {
    if (verbose > 1)
      cout << "   Dataset " << d << " name : " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
		
    //open files and load
    cout<<"Load Dataset"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
		
    string previousFilename = "";
    int iFile = -1;
    
    string dataSetName = datasets[d]->Name();	
    selecTableSemiMu.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
    selecTableSemiEl.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
    
    /////////////////////////////////////
    /// Initialize JEC factors 
    /////////////////////////////////////
		
    //L2L3 residual corrections already in data Toptrees now! (because a global tag is used where these corrections are included)
    vector<JetCorrectorParameters> vCorrParam;
    /*if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
	   JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Jec11V2_db_AK5PFchs_L2L3Residual.txt");
	   vCorrParam.push_back(*ResJetCorPar);
    }*/
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false); // last boolean ('startFromRaw') = false!

    /*if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
     	JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L2Relative.txt");
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L3Absolute.txt");
      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/GR_R_38X_V15_AK5PF_L2L3Residual.txt");
      vCorrParam.push_back(*L2JetCorPar);
      vCorrParam.push_back(*L3JetCorPar);
      vCorrParam.push_back(*ResJetCorPar);
    }
    else // MC!
    {
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L2Relative.txt");
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("JECFiles/START38_V14_AK5PF_L3Absolute.txt");
      vCorrParam.push_back(*L2JetCorPar);
      vCorrParam.push_back(*L3JetCorPar);
    }*/     


    string TreeFileName;
		TreeFileName = Treespath+"InclFourthGenTree_"+dataSetName+postfix+channelpostfix+".root";
    cout << "INFO: creating InclFourthGenTree file "+TreeFileName << endl;        
    TFile* treeFile = new TFile(TreeFileName.c_str(),"RECREATE");      
		
    InclFourthGenTree* myTree = 0;		
    TTree* myInclFourthGenTree = new TTree("myInclFourthGenTree","Tree containing the FourthGen information");
    myInclFourthGenTree->Branch("TheInclFourthGenTree","InclFourthGenTree",&myTree);



    ////////////////////////////////////
    ////////////////////////////////////
    ///////// Loop on events
    ////////////////////////////////////
    ////////////////////////////////////
    int itrigger = -1, previousRun = -1;
      
    //block of code to arrange that only on 1/3th of the ttbar semimu or semiel sample is run for the training, and on 2/3th of the sample for the evaluation. 
    //Don't forget to change the eqLumi in the config depending on training or evaluation (not safe, to be changed? Not so trivial...)
    int start = 0;
    int end = datasets[d]->NofEvtsToRunOver();
    
		/*	
		//to test
		if(TrainMVA && (datasets[d]->Name () == "TTbarJets_SemiMuon" || datasets[d]->Name () == "TTbarJets_SemiElectron"))
    { 
        start = 0;
        end = int(datasets[d]->NofEvtsToRunOver()/3);
    }
    else if (!TrainMVA && (datasets[d]->Name () == "TTbarJets_SemiMuon" || datasets[d]->Name () == "TTbarJets_SemiElectron"))
    {    
        start = int(datasets[d]->NofEvtsToRunOver()/3);
        end = datasets[d]->NofEvtsToRunOver();
    }
    else{
        start = 0;
        end = datasets[d]->NofEvtsToRunOver();
    }
*/			
		
     
    if (verbose > 1)
      cout << " - Loop over events " << endl;      
    
    for (int ievt = start; ievt < end; ievt++)
    {        

      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";
      
			//load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      vector<TRootGenJet*> genjets;
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {
        genjets = treeLoader.LoadGenJet(ievt,false);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      // scale factor for the event
      float scaleFactor = 1.;
                
      // check which file in the dataset it is to have the HLTInfo right
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename)
      {
      	previousFilename = currentFilename;
        iFile++;
	 			cout<<"File changed!!! => iFile = "<<iFile<<endl;
      }
      
      ///////////////////////////////
      // trigger
      ///////////////////////////////
			int currentRun = event->runId();
			if(previousRun != currentRun)
			{
			previousRun = currentRun;
				if(semiMuon)
				{
					if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
					{
						if (event->runId() >= 160431 && event->runId() <= 163261)//May10ReReco
							itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun, iFile);
  					else if (event->runId() >= 163270 && event->runId() <= 163869)
    				  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v6"), currentRun, iFile);
  					else if (event->runId() >= 165088 && event->runId() <= 165633)//PromptReco_v4; splitted over 2 toptrees: 565 and 641
    					itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v8"), currentRun, iFile);
  					else if (event->runId() >= 165970 && event->runId() <= 167043 && event->runId() != 166346)
    					itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v9"), currentRun, iFile);
  					else if (event->runId() == 166346)
    				  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v10"), currentRun, iFile);
  					else if (event->runId() >= 167078 && event->runId() <= 167913)
    				  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v11"), currentRun, iFile);
						else if (event->runId() >= 170249 && event->runId() <= 172619) //Aug05ReReco: equivalent to the run range of PromptReco_v5 normally, but Aug05 replaces this. Warning: somewhere we last about 5/pb in this data?
				  		itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_v8"), currentRun, iFile);
						else if (event->runId() >= 172620 && event->runId() <= 173198) //first part of PromptReco_v6, same as previous trigger
            	itrigger = treeLoader.iTrigger (string ("HLT_IsoMu20_v8"), currentRun, iFile);
						else if (event->runId() >= 173236 && event->runId() <= 173692) //second part of PromptReco_v6
				  		itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_v9"), currentRun, iFile);
				
        			// RUN2011B (promptv1)
   					else if( event->runId() >= 175860 && event->runId() <= 177452 )// TopTree ID 722
   				  	itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
   					else if( event->runId() >=  177718 && event->runId() <=  178380 ) // TopTree ID 804
   				  	itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);
   					else if( event->runId() >=  178420 && event->runId() <=  178479 )
   				  	itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v6"), currentRun, iFile);								
						else if( event->runId() >=  178703 && event->runId() <=  179889 ) // TopTree ID 816
							itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v6"), currentRun, iFile);
						else if( event->runId() >=  179959 && event->runId() <=  180252 )
							itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v7"), currentRun, iFile); 
									   
  					if(itrigger == 9999)
						{
    				  cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
    				  exit(1);
  					}
	   			}
	   			else 
	   			{  
   					itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun, iFile);//Summer11 MC has other triggers!	
    
  					if(itrigger == 9999)
						{
    			  	cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
    			  	exit(1);
						}
					}
				} //end if semiMuon
	 			else if(semiElectron)
				{
	  			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	   			{      		
						// /SingleElectron/Run2011A-May10ReReco-v1/AOD 
						if (event->runId() >= 160404 && event->runId() < 161217)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1"), currentRun, iFile);
						else if (event->runId() >= 161217 && event->runId() < 163270)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2"), currentRun, iFile);
     				else if (event->runId() >= 163270 && event->runId() <= 163869)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3"), currentRun, iFile);				 
						// /ElectronHad/Run2011A-PromptReco-v4/AOD
						else if (event->runId() >= 165088 && event->runId() < 165970)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_BTagIP_v4"), currentRun, iFile);
						else if (event->runId() >= 165970 && event->runId() < 167038)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v1"), currentRun, iFile);
						else if (event->runId() >= 167038 && event->runId() <= 167913)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v2"), currentRun, iFile);			  
						// /ElectronHad/Run2011A-05Aug2011-v1/AOD
						else if (event->runId() >= 170249 && event->runId() <= 172619)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v4"), currentRun, iFile);  
						// /ElectronHad/Run2011A-PromptReco-v6/AOD 
						else if (event->runId() >= 172620 && event->runId() < 173212)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v4"), currentRun, iFile);  
						else if (event->runId() >= 173212 && event->runId() <= 173692)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v5"), currentRun, iFile);  				   				   	
						// RUN2011B (promptv1)
						else if (event->runId() >= 175832 && event->runId() < 178411)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v5"), currentRun, iFile);  				   					else if (event->runId() >= 178411 && event->runId() < 179942)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v8"), currentRun, iFile);  				   					else if (event->runId() >= 179942 && event->runId() <= 180296)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v9"), currentRun, iFile);  				   											   
  					if(itrigger == 9999)
						{
    					cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
    					exit(1);
  					}// semi-electron
 	   			}
	   			else 
	   			{
   					itrigger = treeLoader.iTrigger (string ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2"), currentRun, iFile);//Summer11 MC has other triggers!	
						if(itrigger == 9999)
						{
							cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;	
							exit(1);
						}
					}	 
				} //end if semiElectron	
			} //end previousRun != currentRun


      //ordering is relevant; most probably 1) Type I MET correction, 2) JER where jet corrections are propagated to MET, 3) JES systematics where jet corrections are propagated to MET
      //----------------------------------------------------------
      // Apply type I MET corrections:  (Only for |eta| <= 4.7 )
      //---------------------------------------------------------
      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");      
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
        jetTools->correctMETTypeOne(init_jets,mets[0],true);
      else
        jetTools->correctMETTypeOne(init_jets,mets[0],false);
      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After MET type I correction:");
      
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
      {	
      	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JER correction:");
				if(systematic == "JERMinus")
					jetTools->correctJetJER(init_jets, genjets, mets[0], "minus",false); //false means don't use old numbers but newer ones...
				else if(systematic == "JERPlus")
					jetTools->correctJetJER(init_jets, genjets, mets[0], "plus",false);
				else
					jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal",false);
				//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");	       
		
				// JES systematic! 
				if (systematic == "JESMinus")
					jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
				else if (systematic == "JESPlus")
					jetTools->correctJetJESUnc(init_jets, mets[0], "plus");	       
      }

			double lumiWeight3D = 1.0;
			if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")){
				////////////////////////////
      	// apply trigger Reweighting
      	////////////////////////////
				float mceventtriggerweight = 1;
				float NEWmceventtriggerweight = 1;
				float prob=0;
      	if(semiElectron){
					std::vector<float> probabilities;
					for(size_t i=0; i<init_jets.size(); ++i){
    				if(fabs(init_jets[i]->Eta())>2.6) continue;
    				probabilities.push_back(jetprob(init_jets[i]->Pt(),init_jets[i]->btag_trackCountingHighPurBJetTags())); //highEff when this is used offline
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
				lumiWeight3D = Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
	 			scaleFactor = scaleFactor*lumiWeight3D;
      	histo1D["lumiWeights"]->Fill(scaleFactor);	
			}
						
								
      /////////////////////////////
      // Selection
      /////////////////////////////

     //Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...
      Selection nonstandard_selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected... 
      
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
      {
        // Apply the scraping veto. Note: should be checked if still necessary, maybe already done in toptree production
        bool isBeamBG = true;
        if(event->nTracks() > 10)
        {
          if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
            isBeamBG = false;
      	}
      	if(isBeamBG) continue;
      }	

      bool trigged, isGoodPV;
      trigged = treeLoader.EventTrigged (itrigger);			
      isGoodPV = selection.isPVSelected(vertex, 4, 24., 2); //in the past this was put in the config, but this is not very useful, since the PV cuts are quite standard
			
      bool eventSelected = false;
      bool isSingleLepton = false;
      bool isSSLepton = false;
      bool isTriLepton = false;
      bool isSingleMuon = false;
      bool isSingleElectron = false;
      bool isSSMuon = false;
      bool isSSElectron = false;
      bool isSSMuEl = false;
      bool isTriMuon = false;
      bool isTriElectron = false;
      bool isTriMu2El1 = false;
      bool isTriMu1El2 = false;
      
      vector<TRootJet*> selectedJets;//selectedJetsFromW,selectedJetsFromW_DropUsedJets,selectedJetsFromW_DropUsedJets_tmp;
      vector<TRootJet*> selectedForwardJets, selectedJetsLargeEtaRange;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootElectron*> selectedLooseElectronsNoVBTFid;
      vector<TRootElectron*> selectedLooseElectronsVBTFid;
      vector<TRootMuon*> selectedLooseMuons;
      vector<TRootMCParticle*> mcParticles;
      
      if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_SemiElectron") == 0 || dataSetName.find("NP_Tprime")==0 || dataSetName.find("NP_overlay_Tprime")==0)
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      float METCut = 40;
      selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1);
      nonstandard_selection.setJetCuts(30.,4.7,0.01,1.,0.98,0.3,0.1); //only difference: larger eta acceptance 
      selection.setMuonCuts(20,2.1,0.125,10,0.02,0.3,1,1,1);
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setElectronCuts(20,2.5,0.1,0.02,1,0.3);
      selection.setLooseElectronCuts(15,2.5,0.2);	 				 
      
      if (init_jets.size() > 0)
      {
	    	selectedJets = selection.GetSelectedJets(true);				
	    	selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	    	selectedElectrons = selection.GetSelectedElectrons(vertex[0],selectedJets);
	    	selectedJetsLargeEtaRange = nonstandard_selection.GetSelectedJets(true);	    	    
      }
      selectedLooseElectronsNoVBTFid = selection.GetSelectedLooseElectrons(false); //no vbtfid is required
      selectedLooseElectronsVBTFid = selection.GetSelectedLooseElectrons(true); //vbtfid is required 
      selectedLooseMuons = selection.GetSelectedLooseMuons(); //veto muons	
      
      for(unsigned int i = 0; i<selectedJetsLargeEtaRange.size(); i++)
      {
	    	if(selectedJetsLargeEtaRange[i]->Eta() > 2.4)
	      	selectedForwardJets.push_back(selectedJetsLargeEtaRange[i]);
      }
     
      selecTableSemiMu.Fill(d,1,scaleFactor);
      selecTableSemiEl.Fill(d,1,scaleFactor);
		
			
			
			//// EVENTS TRIGGERED BY MUON TRIGGER			
      if(trigged && semiMuon)
      { 
				selecTableSemiMu.Fill(d,2,scaleFactor);
        if(isGoodPV)
				{
					selecTableSemiMu.Fill(d,3,scaleFactor);
					if(selectedMuons.size()>=1 && selectedMuons[0]->Pt()>40)
					{
						selecTableSemiMu.Fill(d,4,scaleFactor);
						sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class
						
						if(selectedJets.size()>=(unsigned int)anaEnv.NofJets)
						{  //at least 1 jet!
						
						  //block for the jet multiplicity plot
							if(mets[0]->Et()> METCut)
							{							
							      if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsVBTFid.size() == 0)
										{
										  int nBtags = 0;
											for(unsigned int j=0;j<selectedJets.size();j++)
											{
											   if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
												 {
												    nBtags++;
												 }
											}
											MSPlot["MS_JetMultiplicity_SingleLepton"]->Fill(selectedJets.size(),datasets[d], true, Luminosity*scaleFactor);
											MSPlot["MS_BtaggedJetMultiplicity_SingleLepton"]->Fill(nBtags,datasets[d], true, Luminosity*scaleFactor);
										  if(nBtags>0)
											  MSPlot["MS_JetMultiplicityAtleast1Btag_SingleLepton"]->Fill(selectedJets.size(),datasets[d], true, Luminosity*scaleFactor);
										}
							}
								
							//continuing for the selection 	
							for(unsigned int j=0;j<selectedJets.size();j++)
							{
								//now require at least a b-tagged jet larger than a certain pre-defined cut
								if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue && !eventSelected)
								{
									selecTableSemiMu.Fill(d,5,scaleFactor); 
									if(mets[0]->Et()> METCut)
									{
										selecTableSemiMu.Fill(d,6,scaleFactor);
										eventSelected = true; 

										//cout << "event is selected according to the baseline selection!" << endl;
										
										////for single muon require exactly 1 muon, veto for other loose muons and veto for very loose electrons
										if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsVBTFid.size() == 0)
										{
											isSingleLepton = true; // we have a single muon event
											isSingleMuon = true;
											//std::cout<<"Processing the "<<ievt<<"th event" << endl;
											//cout << "is single muon!" << endl;
											//cout << "-> muon pt: " << selectedMuons[0]->Pt() << endl;
			     					}
			     					
										
										////for same-sign muons require exactly 2 muons, veto for very loose electrons
										else if(selectedMuons.size() == 2 && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsVBTFid.size() == 0)
										{ 
											//require that there are the two muons do not form the Z mass
											if( !selection.foundZCandidate(selectedMuons, selectedMuons, 10.) )
											{
												//require the same charge
												if(selectedMuons[0]->charge()== selectedMuons[1]->charge())
												{
													isSSLepton = true; // we have two same-sign muons
													isSSMuon = true;
													//std::cout<<"Processing the "<<ievt<<"th event" << endl;
													//cout << "is same-sign muon!" << endl;
													//cout << "-> muon 1 pt: " << selectedMuons[0]->Pt() << endl;
													//cout << "-> muon 2 pt: " << selectedMuons[1]->Pt() << endl;
												}
											}
										}
										
										////for same-sign muon and electron require exactly 1 muon, veto for other loose muons
										else if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size())
										{
											//require exactly 1 electron
											if(selectedElectrons.size() == 1 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
											{
				    						//require that there are no two electrons forming the Z mass
												if(!selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid,10.))
												{
													//it should not be an electron from a conversion!
													if( selection.passConversionRejection(selectedElectrons[0]) )
													{
														//require the same charge for muon and electron
														if(selectedElectrons[0]->charge()== selectedMuons[0]->charge())
														{
															isSSLepton = true; // we have a same-sign electron and muon
															isSSMuEl = true;
															//std::cout<<"Processing the "<<ievt<<"th event" << endl;
															//cout << "is same-sign muon+electron!" << endl;
															//cout << "-> muon pt: " << selectedMuons[0]->Pt() << endl;
															//cout << "-> electron pt: " << selectedElectrons[0]->Pt() << endl;
														}
													}
												}
											}
										}
										
										////three leptons
										
										//require at least 3 muons, veto for very loose electrons
										else if(selectedMuons.size() == 3  && selectedLooseMuons.size() == selectedMuons.size() && selectedLooseElectronsVBTFid.size() == 0)
										{
											//require that there are no two muons forming the Z mass
											if( !selection.foundZCandidate(selectedMuons, selectedMuons, 10.) )
											{
												isTriLepton = true; // at least three muons
												isTriMuon = true;
												std::cout<<"Processing the "<<ievt<<"th event" << endl;
												//cout << "is trilepton: 3 muons!" << endl;
												//cout << "-> muon 1 pt: " << selectedMuons[0]->Pt() << endl;
												//cout << "-> muon 2 pt: " << selectedMuons[1]->Pt() << endl;
												//cout << "-> muon 3 pt: " << selectedMuons[2]->Pt() << endl;
											}
					 					}
										
										//require 2 muons and an electron
										else if(selectedMuons.size() == 2  && selectedLooseMuons.size() == selectedMuons.size() && selectedElectrons.size() == 1 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
										{
											//require that there are no two muons forming the Z mass
											if( !selection.foundZCandidate(selectedMuons, selectedMuons, 10.) )
											{
												//require exactly 1 electron
													//require that there are no two electrons forming the Z mass
													if(!selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid,10.))
													{
														//it should not be an electron from a conversion!
														if( selection.passConversionRejection(selectedElectrons[0]) )
														{
															isTriLepton = true; //at least two muons and one electron
															isTriMu2El1 = true;
															//std::cout<<"Processing the "<<ievt<<"th event" << endl;
															//cout << "is trilepton: 2 muons + 1 electron!" << endl;
															//cout << "-> muon 1 pt: " << selectedMuons[0]->Pt() << endl;
															//cout << "-> muon 2 pt: " << selectedMuons[1]->Pt() << endl;
															//cout << "-> electron pt: " << selectedElectrons[0]->Pt() << endl;
														}
													}
											}
										}
										
										//require exactly 1 muon
										else if(selectedMuons.size() == 1 && selectedLooseMuons.size() == selectedMuons.size() )
										{
											//require at least 2 electrons
											if(selectedElectrons.size() == 2 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
											{
												//require that there are no two electrons forming the Z mass
												if(!selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.))
												{
													if(selection.passConversionRejection(selectedElectrons[0]) && selection.passConversionRejection(selectedElectrons[1]))
													{
														isTriLepton = true; //one muon and at least two electrons
														isTriMu1El2 = true;
														//std::cout<<"Processing the "<<ievt<<"th event" << endl;
														//cout << "is trilepton: 1 muon + 2 electrons!" << endl;
														//cout << "-> muon pt: " << selectedMuons[0]->Pt() << endl;
														//cout << "-> electron 1 pt: " << selectedElectrons[0]->Pt() << endl;
														//cout << "-> electron 2 pt: " << selectedElectrons[1]->Pt() << endl;
													}
												}
											}
										}
									} // end MET cut
								} // end requirement of at least a b-tagged jet larger than a certain pre-defined cut																		
							} // end 'loop' on jets
						} //end 'at least one jet'  
          } // end if selectedMuons.size()>=1
        } // end good PV
      }// end trigged & semiMuon
      
			///// EVENTS TRIGGERED BY ELECTRON TRIGGER
			else if(trigged && semiElectron)
      {
        selecTableSemiEl.Fill(d,2,scaleFactor);
    	 	
        if( isGoodPV )
        {
          selecTableSemiEl.Fill(d,3,scaleFactor);
          if( selectedElectrons.size() >= 1 && selectedElectrons[0]->Pt()>40)
          {
            selecTableSemiEl.Fill(d,4,scaleFactor);
              if( selection.passConversionRejection(selectedElectrons[0]) )
              {
								selecTableSemiEl.Fill(d,5,scaleFactor);
								sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class

								if( selectedJets.size()>=(unsigned int)anaEnv.NofJets)
								{
									//block for the jet multiplicity plot
									if(mets[0]->Et()> METCut)
									{							
										if(selectedElectrons.size() == 1 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseMuons.size() == 0 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
										{
										  int nBtags = 0;
											for(unsigned int j=0;j<selectedJets.size();j++)
											{
											   if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
												 {
												    nBtags++;
												 }
											}
											MSPlot["MS_JetMultiplicity_SingleLepton"]->Fill(selectedJets.size(),datasets[d], true, Luminosity*scaleFactor);
											MSPlot["MS_BtaggedJetMultiplicity_SingleLepton"]->Fill(nBtags,datasets[d], true, Luminosity*scaleFactor);
											if(nBtags>0)
											  MSPlot["MS_JetMultiplicityAtleast1Btag_SingleLepton"]->Fill(selectedJets.size(),datasets[d], true, Luminosity*scaleFactor);
										}
									}
									
									for(unsigned int j=0;j<selectedJets.size();j++)
									{
										//now require at least a b-tagged jet larger than a certain pre-defined cut
										if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue && !eventSelected)
										{
		             			selecTableSemiEl.Fill(d,6,scaleFactor);
											if(mets[0]->Et()> METCut)
											{
			        					selecTableSemiEl.Fill(d,7,scaleFactor);
                        eventSelected = true;
		
												//cout << "event is selected according to the baseline selection!" << endl;
														
												//// single electron
												if(selectedElectrons.size() == 1 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseMuons.size() == 0 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													isSingleLepton = true;
													isSingleElectron = true;
													//std::cout<<"Processing the "<<ievt<<"th event" << endl;
													//cout << "is single electron!" << endl;
													//cout << "-> electron pt: " << selectedElectrons[0]->Pt() << endl;
												}
												
												//// two same-sign electrons
												else if(selectedElectrons.size() == 2 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													if(selection.passConversionRejection(selectedElectrons[1]))
													{
														if(selectedElectrons[0]->charge()== selectedElectrons[1]->charge())
														{
															isSSLepton = true;
															isSSElectron = true;
															//std::cout<<"Processing the "<<ievt<<"th event" << endl;
															//cout << "is same-sign electron!" << endl;
															//cout << "-> electron 1 pt: " << selectedElectrons[0]->Pt() << endl;
															//cout << "-> electron 2 pt: " << selectedElectrons[1]->Pt() << endl;
														}
													}
												}
												
												//// a same-sign electron and muon
												else if(selectedElectrons.size() == 1 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													if(selectedMuons.size() == 1 && selectedMuons[0]->Pt()<40 && selectedLooseMuons.size() == selectedMuons.size())
													{
														if(selectedElectrons[0]->charge()== selectedMuons[1]->charge())
														{
															isSSLepton = true;
															isSSMuEl = true;
															//std::cout<<"Processing the "<<ievt<<"th event" << endl;
															//cout << "is same-sign electron + muon!" << endl;
															//cout << "-> electron pt: " << selectedElectrons[0]->Pt() << endl;
															//cout << "-> muon pt: " << selectedMuons[0]->Pt() << endl;
														}
													}
												}
											
												//// 1 electron and 2 muons
												else if(selectedElectrons.size() == 1 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													if(selectedMuons.size() == 2 && selectedMuons[1]->Pt()<40 && selectedLooseMuons.size() == selectedMuons.size())
													{
														isTriLepton = true;
														isTriMu2El1 = true;
														//std::cout<<"Processing the "<<ievt<<"th event" << endl;
														//cout << "is tri-lepton: 1 electron and 2 muons!" << endl;
														//cout << "-> electron pt: " << selectedElectrons[0]->Pt() << endl;
														//cout << "-> muon 1 pt: " << selectedMuons[0]->Pt() << endl;
														//cout << "-> muon 2 pt: " << selectedMuons[1]->Pt() << endl;
													}
												}
												
												//// 2 electrons and 1 muon
												else if(selectedElectrons.size() == 2 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													if(selection.passConversionRejection(selectedElectrons[1])&& selectedElectrons[1]->Pt()<40)
													{
														if(selectedMuons.size() == 1 && selectedMuons[0]->Pt()<40 && selectedLooseMuons.size() == selectedMuons.size())
														{
															isTriLepton = true;
															isTriMu1El2 = true;
															//std::cout<<"Processing the "<<ievt<<"th event" << endl;
															//cout << "is tri-lepton: 2 electrons and 1 muon!" << endl;
															//cout << "-> electron 1 pt: " << selectedElectrons[0]->Pt() << endl;
															//cout << "-> electron 2 pt: " << selectedElectrons[1]->Pt() << endl;
															//cout << "-> muon 1 pt: " << selectedMuons[0]->Pt() << endl;
														}
													}
												}
												
												//// three electrons
												else if(selectedElectrons.size() == 3 && !selection.foundZCandidate(selectedElectrons, selectedLooseElectronsNoVBTFid, 10.) && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
												{
													if(selection.passConversionRejection(selectedElectrons[1]) && selection.passConversionRejection(selectedElectrons[2]))
													{
														isTriLepton = true; // at least three electrons
														isTriElectron = true;
														//std::cout<<"Processing the "<<ievt<<"th event" << endl;
														//cout << "is trilepton: 3 electrons!" << endl;
														//cout << "-> electron 1 pt: " << selectedElectrons[0]->Pt() << endl;
														//cout << "-> electron 2 pt: " << selectedElectrons[1]->Pt() << endl;
														//cout << "-> electron 3 pt: " << selectedElectrons[2]->Pt() << endl;
													}
												}
											} // end MET cut
										} // end requirement of at least a b-tagged jet larger than a certain pre-defined cut
									} // end 'loop' on jets
								} // end 'at least one jet'
							} // end conversion rejection for leading electron
          } // end if selectedElectrons.size()>=1
        } // end good PV
      } // end trigged & semiElectron

						
      if(!isSingleLepton && !isSSLepton && !isTriLepton) continue; //same as all cuts just above (baseline selection is there) 

			MSPlot["MS_nPV"]->Fill(vertex.size(),datasets[d], true, Luminosity*scaleFactor);				

			if(isSingleLepton){
				MSPlot["MS_MET"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);				
				if(semiElectron) MSPlot["MS_LeptonPt"]->Fill(selectedElectrons[0]->Pt(),datasets[d], true, Luminosity*scaleFactor);				
				if(semiMuon) MSPlot["MS_LeptonPt"]->Fill(selectedMuons[0]->Pt(),datasets[d], true, Luminosity*scaleFactor);				
			
				for(unsigned int j=0;j<selectedJets.size();j++)
				{
				  MSPlot["MS_JetPt_all_SingleLepton"]->Fill(selectedJets[j]->Pt(),datasets[d], true, Luminosity*scaleFactor);
					if(selectedJets[j]->btag_trackCountingHighPurBJetTags() > workingpointvalue)
					  MSPlot["MS_JetPt_btagged_SingleLepton"]->Fill(selectedJets[j]->Pt(),datasets[d], true, Luminosity*scaleFactor);
					else
					  MSPlot["MS_JetPt_nonbtagged_SingleLepton"]->Fill(selectedJets[j]->Pt(),datasets[d], true, Luminosity*scaleFactor);
				}			
				if(semiElectron) selecTableSemiEl.Fill(d,8,scaleFactor);
				if(semiMuon) selecTableSemiMu.Fill(d,7,scaleFactor);				
			}

			if(isSSLepton)
			{
				//cout << "IS SAME-SIGN LEPTON EVENT" << endl;
				NbSSevents = NbSSevents + datasets[d]->NormFactor()*Luminosity*scaleFactor;
				selecTableMultiLepton.Fill(d,0,scaleFactor);
			}
			if(isTriLepton)
			{
				//cout << "IS TRI-LEPTON EVENT" << endl;
				NbTrievents = NbTrievents + datasets[d]->NormFactor()*Luminosity*scaleFactor;
				selecTableMultiLepton.Fill(d,1,scaleFactor);				
			}			
			
			
			TRootGenEvent* genEvt = 0;
			if((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron))
			  genEvt = treeLoader.LoadGenEvent(ievt,false);

						
			bool Wbosonpartonsmatched = false; // True if the Wboson ttbar semi-mu or semi-el partons are matched to their 2 jets (not necessarily the 4 highest pt jets)
			float WMassmatched_ = -9999;
			//////////////////////////////////////////////////////////////////////////
      // jet-parton matching needed to make the Wmassplot for the W counting
      //////////////////////////////////////////////////////////////////////////
      sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)   
			
			if(selectedJets.size()>=4 && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron)))
			{
			
			  pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
			  leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
			
				//sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class)     
				int MCPermutation[4]; for (unsigned int i=0;i<4;i++) MCPermutation[i] = -1;

				vector<TRootMCParticle*> mcParticlesMatching_; 
				vector<TRootJet*> selectedJets_; 

				mcParticlesMatching_.clear();
				selectedJets_.clear();
      
					vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
            
					bool muPlusFromTop = false, muMinusFromTop = false;
					bool elPlusFromTop = false, elMinusFromTop = false;
					int leptonPDG, muonPDG = 13, electronPDG = 11;
					if(semiMuon)
					  leptonPDG = muonPDG;
					else if(semiElectron)
					  leptonPDG = electronPDG; 
					for(unsigned int i=0; i<mcParticles.size(); i++) {
						if( mcParticles[i]->status() != 3) continue;
						
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
						
						if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 )
						{
							mcParticlesTLV.push_back(*mcParticles[i]);
							mcParticlesMatching_.push_back(mcParticles[i]);
						}
					}
	
					//if(muPlusFromTop && muMinusFromTop)
	  				//cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
	
					for(unsigned int i=0; i<4; i++)
	  				  selectedJetsTLV.push_back(*selectedJets[i]);
	
					JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
	
					if(matching.getNumberOfAvailableCombinations() != 1)
	  				  cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
	
					vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
	
					for(unsigned int i=0; i<mcParticlesTLV.size(); i++) 
					{
	  				  int matchedJetNumber = matching.getMatchForParton(i, 0);
	  				  if(matchedJetNumber != -1)
	    				    JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
					}

					for(unsigned int i=0; i<JetPartonPair.size(); i++)
					{
	  				  unsigned int j = JetPartonPair[i].second;	  
	  				  if( fabs(mcParticlesMatching_[j]->type()) < 6 )
	    			          {
	      			            if( ( ((muPlusFromTop && semiMuon) || (elPlusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -6 )
		  				|| ( ((muMinusFromTop && semiMuon) || (elMinusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == 6 ) )
					    {
		  				if(hadronicWJet1_.first == 9999)
						{
		    					hadronicWJet1_ = JetPartonPair[i];
		    					MCPermutation[0] = JetPartonPair[i].first;
		  				}
						else if(hadronicWJet2_.first == 9999)
						{
		   	 				hadronicWJet2_ = JetPartonPair[i];
		    					MCPermutation[1] = JetPartonPair[i].first;
		  				} 
						else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
					    }
	    			          }
	  				  if( fabs(mcParticlesMatching_[j]->type()) == 5 )
	    			          {
	      			            if( ( ((muPlusFromTop && semiMuon) || (elPlusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == -6 )
		    			|| ( ((muMinusFromTop && semiMuon) || (elMinusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == 6 ) )
					    {
								hadronicBJet_ = JetPartonPair[i];
								MCPermutation[2] = JetPartonPair[i].first;
	      				     }
					     else if( ( ((muPlusFromTop && semiMuon) || (elPlusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == 6 )
			 				|| ( ((muMinusFromTop && semiMuon) || (elMinusFromTop && semiElectron)) && mcParticlesMatching_[j]->motherType() == -6 ) )
					     {
			 					leptonicBJet_ = JetPartonPair[i];
			 					MCPermutation[3] = JetPartonPair[i].first;
	      				     }
	    				  }
					}					
	
					if(hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999)
					{
	  				  histo1D["hadronicPartonWMass"]->Fill((*mcParticlesMatching_[hadronicWJet1_.second]+*mcParticlesMatching_[hadronicWJet2_.second]).M(),scaleFactor);	  
	  				  Wbosonpartonsmatched = true;	  
					}					

      	if(Wbosonpartonsmatched)
	  		{
					WMassmatched_ = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]).M();
      		histo1D["hadronicRecoWMass"]->Fill(WMassmatched_,scaleFactor);
	  		}	
			
			} //end selectedJets.size()>=4 !useMassesAndResolutions && dataSet semiMu or semiEl

     
		  //now do some stuff for later use in the MVA (obtain the quarks that should be used for matching later on)
			//bool all4JetsMatched_MCdef = false; // True if the 4 highest pt jets are matched to the 4 ttbar semi-lep partons
  		//bool hadronictopJetsMatched_MCdef = false;
			vector<TLorentzVector> mcQuarksForMatching; //ordering: hadronicWQuark1,hadronicWQuark2,hadronicbQuark,leptonicbQuark
			int pdgID_top = 6; //top quark
  		//if(TprimeEvaluation)
    	//	pdgID_top = 8; //4th generation t' quark
			if(selectedJets.size()>=4 && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron)))
			{
    		TLorentzVector hadronicWQuark1,hadronicWQuark2,hadronicbQuark,leptonicbQuark;
    		bool muPlusFromTop = false, muMinusFromTop = false, elPlusFromTop = false, elMinusFromTop = false;
    
    		for(unsigned int i=0; i<mcParticles.size(); i++)
    		{
    		  //cout << i << ":  status: " << mcParticles[i]->status() << "  pdgId: " << mcParticles[i]->type()
    		  //  << "  motherPdgId: " << mcParticles[i]->motherType() << "  grannyPdgId: " << mcParticles[i]->grannyType() << endl;
    		  if( mcParticles[i]->status() != 3) continue;
          
      		if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )
      		  muMinusFromTop = true;
      		if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )
      		  muPlusFromTop = true;
      		if( mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )
      		  elMinusFromTop = true;
      		if( mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )
      		  elPlusFromTop = true;
				}
				
				bool WQuark1Found = false;
        for(unsigned int i=0; i<mcParticles.size(); i++)
    		{
			  	if( mcParticles[i]->status() != 3) continue;
					
      		if( abs(mcParticles[i]->type()) < 6) //light/b quarks, 6 should stay hardcoded
      		{
					  if( ((( muPlusFromTop || elPlusFromTop ) && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top) 
						  || (( muMinusFromTop || elMinusFromTop ) && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top)) 
						&& !WQuark1Found)
        		{							  											
							hadronicWQuark1 = (TLorentzVector) *mcParticles[i];
							//cout<<"hadronicWQuark1.Pt() = "<<hadronicWQuark1.Pt()<<endl;
							WQuark1Found = true;
						}
						else if ( (( muPlusFromTop || elPlusFromTop ) && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top) 
						  || (( muMinusFromTop || elMinusFromTop ) && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top ))
					  {
							hadronicWQuark2 = (TLorentzVector) *mcParticles[i];
							//cout<<"hadronicWQuark2.Pt() = "<<hadronicWQuark2.Pt()<<endl;
						}
						else if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticles[i]->motherType() == -pdgID_top )
          		|| ( ( muMinusFromTop || elMinusFromTop ) && mcParticles[i]->motherType() == pdgID_top ) )
						{
						  hadronicbQuark = (TLorentzVector) *mcParticles[i];
							//cout<<"hadronicbQuark.Pt() = "<<hadronicbQuark.Pt()<<endl;
						}
						else if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticles[i]->motherType() == pdgID_top )
          		|| ( ( muMinusFromTop || elMinusFromTop ) && mcParticles[i]->motherType() == -pdgID_top ) )
						{
							leptonicbQuark = (TLorentzVector) *mcParticles[i];
							//cout<<"leptonicbQuark.Pt() = "<<leptonicbQuark.Pt()<<endl;
						}						
					}
					if( abs(mcParticles[i]->type()) == 5) //light/b quarks, 6 should stay hardcoded
      		{
					  if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticles[i]->motherType() == -pdgID_top )
          		|| ( ( muMinusFromTop || elMinusFromTop ) && mcParticles[i]->motherType() == pdgID_top ) )
						{
						  hadronicbQuark = (TLorentzVector) *mcParticles[i];
							//cout<<"hadronicbQuark.Pt() = "<<hadronicbQuark.Pt()<<endl;
						}
						else if( ( ( muPlusFromTop || elPlusFromTop ) && mcParticles[i]->motherType() == pdgID_top )
          		|| ( ( muMinusFromTop || elMinusFromTop ) && mcParticles[i]->motherType() == -pdgID_top ) )
						{
							leptonicbQuark = (TLorentzVector) *mcParticles[i];
							//cout<<"leptonicbQuark.Pt() = "<<leptonicbQuark.Pt()<<endl;
						}
				  }
				}
				//respect this ordering
				mcQuarksForMatching.push_back(hadronicWQuark1);
				mcQuarksForMatching.push_back(hadronicWQuark2);
				mcQuarksForMatching.push_back(hadronicbQuark);
				mcQuarksForMatching.push_back(leptonicbQuark);
			} //end selectedJets.size()>=4 && dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_SemiElectron") == 0

/*
	    //////////////////////////////////////////////////////////////////////////
      // MVA training, to test, OLD way
      //////////////////////////////////////////////////////////////////////////
      if(doMVAjetcombination && selectedJets.size()>=4 && TrainMVA && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron))) //otherwise, if the jets vector only has 2 jets selected the jetcombiner crashes... 
      {
				if(!isSingleLepton) continue;
	
        			TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
        			//sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class 
				if(semiMuon) 
					jetCombiner->ProcessEvent(datasets[d],mcParticles,selectedJets,selectedMuons[0],init_electrons,init_muons,genEvt,scaleFactor,false);	        
				else if(semiElectron)
					jetCombiner->ProcessEvent(datasets[d],mcParticles,selectedJets,selectedElectrons[0],init_electrons,init_muons,genEvt,scaleFactor,false);
      }

      if(TrainMVA) continue; //for the training, only the jetcombiner is relevant, so the following can be skipped (to the next event in the event loop)
*/

      vector<float> bTagTCHE, bTagTCHP, InitJetsbTagTCHE, InitJetsbTagTCHP;
			vector<int> partonFlavourJet;
      vector<TLorentzVector> SelectedJetsTLV, SelectedForwardJetsTLV, SelectedMuonsTLV, SelectedElectronsTLV, InitJets;
      for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++)
      {
            SelectedJetsTLV.push_back( *selectedJets[iJet] );
            bTagTCHE.push_back(selectedJets[iJet]->btag_trackCountingHighEffBJetTags());
            bTagTCHP.push_back(selectedJets[iJet]->btag_trackCountingHighPurBJetTags());
						partonFlavourJet.push_back(selectedJets[iJet]->partonFlavour());
      }
			for(unsigned int iJet=0; iJet<selectedForwardJets.size(); iJet++)
      {
			      SelectedForwardJetsTLV.push_back( *selectedForwardJets[iJet] );			  	
			}
			for(unsigned int iJet=0; iJet<init_jets.size(); iJet++)
			{
						InitJets.push_back( *init_jets[iJet] ); //needed for MC trigger reweighting in e-channel...
						InitJetsbTagTCHE.push_back(init_jets[iJet]->btag_trackCountingHighEffBJetTags());
            InitJetsbTagTCHP.push_back(init_jets[iJet]->btag_trackCountingHighPurBJetTags());			
			}
			for(unsigned int iMuon=0; iMuon<selectedMuons.size(); iMuon++)
      {
            SelectedMuonsTLV.push_back( *selectedMuons[iMuon] );
      }
			for(unsigned int iElectron=0; iElectron<selectedElectrons.size(); iElectron++)
      {
            SelectedElectronsTLV.push_back( *selectedElectrons[iElectron] );
      }
						
				
/*					
			//////////////////////////////////////////////////////////////////////////
      // MVA training , to test, NEW way
      //////////////////////////////////////////////////////////////////////////
      if(doMVAjetcombination && selectedJets.size()>=4 && TrainMVA && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) || (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron))) //otherwise, if the jets vector only has 2 jets selected the jetcombiner crashes... 
      {
				if(!isSingleLepton) continue;
				
				TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);	//well, needed now but not for the new way eventually		
        //the jets are sorted according to Pt in the TreeCreator, with everything conistently sorted along (should be), because the sorting was done before matching and pushing back b-tag values in vectors... Check and be careful!!
				//sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class 				
				if(semiMuon) 
					jetCombiner->ProcessEvent(datasets[d], mcQuarksForMatching, SelectedJetsTLV, bTagTCHE, SelectedMuonsTLV[0], genEvt->isSemiLeptonic( TRootGenEvent::kMuon ), scaleFactor, false);	        
				else if(semiElectron)
					jetCombiner->ProcessEvent(datasets[d], mcQuarksForMatching, SelectedJetsTLV, bTagTCHE, SelectedElectronsTLV[0], genEvt->isSemiLeptonic( TRootGenEvent::kElec ), scaleFactor, false);	        
      }

       
      if(TrainMVA) continue; //for the training, only the jetcombiner is relevant, so the following can be skipped (to the next event in the event loop)			
*/


			myTree = new InclFourthGenTree();
      myTree->setEventID( event->eventId() );
      myTree->setRunID( event->runId() );
      myTree->setLumiBlockID( event->lumiBlockId() );
      myTree->setNPV( vertex.size() );
      myTree->setNPUBXm1( event->nPu(-1) );
      myTree->setNPU( event->nPu(0) );
      myTree->setNPUBXp1( event->nPu(1) );
			
      myTree->setSelectedSingleLepton( isSingleLepton );
			myTree->setSelectedSingleMu( isSingleMuon );
			myTree->setSelectedSingleEl( isSingleElectron );
			myTree->setSelectedSSLepton( isSSLepton );
			myTree->setSelectedSSMu( isSSMuon );
			myTree->setSelectedSSEl( isSSElectron );
			myTree->setSelectedSSMuEl( isSSMuEl );
			myTree->setSelectedTriLepton( isTriLepton );
			myTree->setSelectedMuMuMu( isTriMuon );
			myTree->setSelectedMuMuEl( isTriMu2El1 );
			myTree->setSelectedMuElEl( isTriMu1El2 );
			myTree->setSelectedElElEl( isTriElectron );
			
      if(dataSetName.find("TTbarJets_Semi") == 0)
      {
      	myTree->setSemiMuDecay(genEvt->isSemiLeptonic( TRootGenEvent::kMuon ));
      	myTree->setSemiElDecay(genEvt->isSemiLeptonic( TRootGenEvent::kElec ));
				myTree->setWbosonpartonsmatched(Wbosonpartonsmatched);
				if(Wbosonpartonsmatched)
				  myTree->setWMassmatched(WMassmatched_);
				myTree->setmcQuarksForMatching( mcQuarksForMatching );
				//myTree->setAll4JetsMCMatched( all4JetsMatched_MCdef );
        //myTree->setAllHadronicJetsMCMatched( hadronictopJetsMatched_MCdef );
				//vector<unsigned int> MatchedJetsIndices;
				//MatchedJetsIndices.push_back(hadronicWJet1_.first);
				//MatchedJetsIndices.push_back(hadronicWJet2_.first);
				//MatchedJetsIndices.push_back(hadronicBJet_.first);
				//MatchedJetsIndices.push_back(leptonicBJet_.first);
				//myTree->setMatchedJetsIndices( MatchedJetsIndices );
      }
			myTree->setEventWeight( scaleFactor );
			myTree->setMET( *mets[0] );
      myTree->setSelectedJets( SelectedJetsTLV );
			myTree->setBTagTCHE( bTagTCHE );
      myTree->setBTagTCHP( bTagTCHP );
			myTree->setpartonFlavourJet( partonFlavourJet );
			myTree->setSelectedForwardJets( SelectedForwardJetsTLV );
			myTree->setInitJets( InitJets );
			myTree->setInitJetsBTagTCHE( InitJetsbTagTCHE );
      myTree->setInitJetsBTagTCHP( InitJetsbTagTCHP );				     						
      myTree->setMuons( SelectedMuonsTLV );
      myTree->setElectrons( SelectedElectronsTLV );
			
			
			//myTree->setTopDecayedLept( topDecayedLept );
      //myTree->setAll4JetsMCMatched( jetCombiner->All4JetsMatched_MCdef() );
      //myTree->setAllHadronicJetsMCMatched( jetCombiner->HadronicTopJetsMatched_MCdef() );
      //myTree->setMvaVals(mvaValsVector);
      //myTree->setMvaResults(mvaResultsVector);
      //myTree->setHadrBJet( hadrBJetIndex );
      //myTree->setHadrLJet1( lightJet1Index );
      //myTree->setHadrLJet2( lightJet2Index );
      //myTree->setLeptBJet( leptBJetIndex );      
      //myTree->setHadrBQuark( jetCombiner->GetHadrBQuark() );
      //myTree->setHadrLQuark1( jetCombiner->GetLightQuark1() );
      //myTree->setHadrLQuark2( jetCombiner->GetLightQuark2() );
      //myTree->setLeptBQuark( jetCombiner->GetLeptBQuark() );
			
			
            
      myInclFourthGenTree->Fill();
	    delete myTree;

	
    }//loop on events

    
    cout<<endl;
		
		
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

		myInclFourthGenTree->Write();
		treeFile->Close();
		delete treeFile;
		
    if(jetTools) delete jetTools;
		
		//important: free memory
    treeLoader.UnLoadDataset();
	
    
  } //loop on datasets
  
	
  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;
 
/* 
  //to test
  if(doMVAjetcombination)
  {
    string pathPNGJetCombi = pathPNG+"JetCombination/";
    if(savePNG) mkdir(pathPNGJetCombi.c_str(),0777);
    jetCombiner->Write(fout, savePNG, pathPNGJetCombi);
  }
*/	
	
	
  if(!TrainMVA)
  {
    fout->cd();
    //Write histograms: MSPlots
    if(savePNG) mkdir((pathPNG+"MSPlot/").c_str(),0777);
    //cout << "mkdir " << (pathPNG+"MSPlot/").c_str()<< endl;
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
			if(savePNG) tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() ); //well, is actually not png but pdf...
    }    
    cout << "1D plots written" << endl;
    
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
    selecTableSemiMu.TableCalculator(true, true, true, true, true, true, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
    string selectiontableSemiMu = "InclFourthGenSearch_SelectionTable_SemiMu"+postfix;
    selectiontableSemiMu = selectiontableSemiMu +".tex"; 	
    if(semiMuon) selecTableSemiMu.Write(selectiontableSemiMu.c_str(),false, true, false, false, false, false, false); //(filename, error, merged, lines, unscaled, eff, totaleff, landscape)
	
    selecTableSemiEl.TableCalculator(true, true, true, true, true, true, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
    string selectiontableSemiEl = "InclFourthGenSearch_SelectionTable_SemiEl"+postfix;		
    selectiontableSemiEl = selectiontableSemiEl +".tex"; 	
    if(semiElectron) selecTableSemiEl.Write(selectiontableSemiEl.c_str(),false, true, false, false, false, false, false);
    
    selecTableMultiLepton.TableCalculator(true, true, true, true, true, true, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
    string selectiontableMultiLepton = "InclFourthGenSearch_SelectionTable_MultiLepton"+postfix+channelpostfix;
    selectiontableMultiLepton = selectiontableMultiLepton +".tex"; 	
    selecTableMultiLepton.Write(selectiontableMultiLepton.c_str(),false, true, false, false, false, false, false);

    cout << " - Closing the output file now..." << endl;
    fout->Close();
  } //end !trainMVA
 
/*  
  //delete, to test
  if(jetCombiner) delete jetCombiner; //IMPORTANT!! file for training otherwise not filled... (?) //crashes when calculating resolutions for kinfit
*/ 
 
 
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
