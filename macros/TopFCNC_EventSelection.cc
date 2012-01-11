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
/*
#include "../Selection/interface/ElectronPlotter.h"
#include "../Selection/interface/MuonPlotter.h"
#include "../Selection/interface/JetPlotter.h"
#include "../Selection/interface/VertexPlotter.h"
#include "../Tools/interface/MVATrainer.h"
#include "../Tools/interface/MVAComputer.h"
*/
#include "../Tools/interface/JetTools.h"
//#include "../Tools/interface/TwoDimTemplateTools.h"
//#include "../Tools/interface/InclFourthGenSearchTools.h"
#include "../JESMeasurement/interface/JetCombiner.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../Reconstruction/interface/MakeBinning.h"
#include "../Reconstruction/interface/TTreeObservables.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
//#include "../MCInformation/interface/LumiReWeighting.h" 
//for Kinematic Fit
//#include "../MCInformation/interface/ResolutionFit.h"
//#include "../KinFitter/interface/TKinFitter.h"
//#include "../KinFitter/interface/TFitConstraintM.h"
//#include "../KinFitter/interface/TFitParticleEtThetaPhi.h"

#include "Style.C"

using namespace std;
using namespace TopTree;
using namespace RooFit;
//using namespace reweight;


struct sort_pair_decreasing
{
    bool operator()(const std::pair<int,float> &left, const std::pair<int,float> &right)
    {
        return left.second > right.second;
    }
};


/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

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


  string btagger = "TCHEM";
// b-tag scalefactor => TCHEL: data/MC scalefactor = 0.95 +- 0.10,    TCHEM: data/MC scalefactor = 0.94 +- 0.09
// mistag scalefactor => TCHEL: data/MC scalefactor = 1.11 +- 0.12,    TCHEM: data/MC scalefactor = 1.21 +- 0.17
  float scalefactorbtageff, mistagfactor;
  if(btagger == "TCHEL") //track counting high eff loose working point
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
  if(btagger == "TCHEL")
     workingpointvalue = 1.7;
  else if(btagger == "TCHEM")
     workingpointvalue = 3.3;	

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the fourth generation search ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  setTDRStyle();
  //setMyStyle();

  string postfix = "TEST"; // to relabel the names of the output file

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


  /////////////////////
  // Configuration
  /////////////////////

  string channelpostfix = "";
  bool diElectron = false; // use diElectron channel?
  bool diMuon = true; // use diMuon channel?
  if(diElectron && diMuon)
  {
     cout << "  --> Using both diMuon and diElectron channel? Choose only one (for the moment, since this requires running on different samples/skims)!" << endl;
     exit(1);
  }
  else
  {
    if(diMuon){
       cout << " --> Using the diMuon channel..." << endl;
       channelpostfix = "_diMu";
    }
    else if(diElectron){
       cout << " --> Using the diElectron channel..." << endl;
       channelpostfix = "_diEl";
    }
  }

  //xml file
  string xmlFileName = "";
  float Luminosity = 1;
  if(diElectron){
	xmlFileName = "../config/myTopFCNCconfig_Electron.xml";
	Luminosity  = 4593.348;
  }
  else if(diMuon){
	xmlFileName = "../config/myTopFCNCconfig_Muon.xml";
	Luminosity  = 4534.871;
  }
  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    
  
  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = anaEnv.Verbose;

  //float anaEnvLuminosity = anaEnv.Luminosity;	// in 1/pb 
  //cout << "analysis environment luminosity for rescaling "<< anaEnvLuminosity << endl;

  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  
  //is this block needed?
/*
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
*/  
  //Output ROOT file
  string rootFileName ("TopFCNCSearch"+postfix+channelpostfix+".root");
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

  string pathPNG = "TopFCNCSearchPlots"+postfix+channelpostfix;
  pathPNG = pathPNG +"/"; 	
  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  //Most 1D and MS plots are declared inside this class
  
  MSPlot["MET"]               = new MultiSamplePlot(datasets, "MET", 500, 0, 500, "\\slashed{E_T}");
  MSPlot["NbOfLooseMuon"]     = new MultiSamplePlot(datasets, "NbOfLooseMuon", 10, 0, 10, "Nb. of loose muons");
  MSPlot["NbOfLooseElectron"] = new MultiSamplePlot(datasets, "NbOfLooseElectron", 10, 0, 10, "Nb. of loose electrons");

  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
  
  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTableDiMu;
  CutsSelecTableDiMu.push_back(string("initial"));
  CutsSelecTableDiMu.push_back(string("preselected"));
  CutsSelecTableDiMu.push_back(string("trigged"));
  CutsSelecTableDiMu.push_back(string("Good PV"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 2 isolated muon"));
  CutsSelecTableDiMu.push_back(string("opposite charge and $|m_{ll}-m_Z|<20$"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 1 jet"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 2 jet"));
  CutsSelecTableDiMu.push_back(string("$\\geq$ 3 jet"));
  CutsSelecTableDiMu.push_back(string("MET < 30 GeV"));  

  vector<string> CutsSelecTableDiEl;
/*
  CutsSelecTableDiEl.push_back(string("initial"));
//CutsSelecTableDiEl.push_back(string("preselected"));
  CutsSelecTableDiEl.push_back(string("trigged"));
  CutsSelecTableDiEl.push_back(string("Good PV"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 1 selected electron"));
  CutsSelecTableDiEl.push_back(string("Veto muon"));
  CutsSelecTableDiEl.push_back(string("Conversion veto"));
  CutsSelecTableDiEl.push_back(string("$\\geq$ 1 b-tagged jet with pt > 50"));
  CutsSelecTableDiEl.push_back(string("MET > 40 GeV"));
*/
  SelectionTable selecTableDiMu(CutsSelecTableDiMu, datasets);
  selecTableDiMu.SetLuminosity(Luminosity);
  SelectionTable selecTableSemiEl(CutsSelecTableDiEl, datasets);
  selecTableSemiEl.SetLuminosity(Luminosity);
  
  cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////////////////////
  // PileUp Reweighting - 3D//
  ////////////////////////////////////////////////////
	Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun173692.root", "pileup", "pileup");
	Lumi3DWeights.weight3D_init(1.0);

  cout << " - Initialized LumiReWeighting stuff" << endl;
  
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
    //selecTableDiMu.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
    //selecTableSemiEl.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );
    
    // scale factor for the event
    float scaleFactor = 1.;
    
    // cross sections and weights
    if (dataSetName == "data")			scaleFactor = 1;
    else if (dataSetName == "tt")		scaleFactor = Luminosity*163/3160707;
    else if (dataSetName == "tt2l")		scaleFactor = Luminosity*17.10/8576584;
    else if (dataSetName == "twdr")		scaleFactor = Luminosity*7.87/813743;
    else if (dataSetName == "atwdr")		scaleFactor = Luminosity*7.87/689462;
    else if (dataSetName == "twds")		scaleFactor = Luminosity*7.87/794802;
    else if (dataSetName == "atwds")		scaleFactor = Luminosity*7.87/784764;
    else if (dataSetName == "t")		scaleFactor = Luminosity*41.92/3337875;
    else if (dataSetName == "at")		scaleFactor = Luminosity*22.65/1943627;
    else if (dataSetName == "ts")		scaleFactor = Luminosity*3.19/259777;
    else if (dataSetName == "ats")		scaleFactor = Luminosity*1.44/137889;
    else if (dataSetName == "wjets")		scaleFactor = Luminosity*31314/49708092;
    else if (dataSetName == "zjets")		scaleFactor = Luminosity*3048/26523984;
 
    else if (dataSetName == "dymumu")		scaleFactor = Luminosity*1666/1;
    else if (dataSetName == "dyee")		scaleFactor = Luminosity*1666/1;
    else if (dataSetName == "dytautau")		scaleFactor = Luminosity*1666/1;
    else if (dataSetName == "ww")		scaleFactor = Luminosity*42.9/4223785;
    else if (dataSetName == "wz")		scaleFactor = Luminosity*18.3/3863081;
    else if (dataSetName == "zz")		scaleFactor = Luminosity*7.67/4188624;
    else if (dataSetName == "qcd_mu")		scaleFactor = Luminosity*84679.3/25079892;

    /////////////////////////////////////
    /// Initialize JEC factors 
    /////////////////////////////////////
		
    //L2L3 residual corrections already in data Toptrees now! (because a global tag is used where these corrections are included)
    vector<JetCorrectorParameters> vCorrParam;
    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false); // last boolean ('startFromRaw') = false!    

    ////////////////////////////////////
    //////////////////////////////////// 
    ///////// Loop on events
    ////////////////////////////////////
    ////////////////////////////////////

    int itrigger = -1, previousRun = -1;
      
    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();
     
    if (verbose > 1) cout << " - Loop over events " << endl;      
    
    for (unsigned int ievt = start; ievt < end; ievt++)
    {        

	//if(ievt%1000 == 0)
	//std::cout<<"Processing the "<<ievt<<"th event ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

	//load event
	event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
	vector<TRootGenJet*> genjets;
	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
		genjets = treeLoader.LoadGenJet(ievt,false);
		sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
	}

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
		if(diMuon)
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
		} //end if diMuon
	} //end previousRun != currentRun

	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
		// JES systematic! 
		if (doJESShift == 1)
			jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
		else if (doJESShift == 2)
			jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
	
		//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JER correction:");
	
		if(doJERShift == 1)
			jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
		else if(doJERShift == 2)
			jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
		else
			jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
	  
		//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");	       
	}

	////////////////////////////
	// if data, apply beam scrapping veto and PU reweighting
	////////////////////////////

	if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"))
	{
		// Apply the scraping veto. Note: should be checked if still necessary, maybe already done in toptree production
        	bool isBeamBG = true;
        	if(event->nTracks() > 10)
        	{
			if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
			isBeamBG = false;
		}
      		if(isBeamBG) continue;

		double lumiWeight3D = Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
	 	scaleFactor *= lumiWeight3D;
	}
	histo1D["lumiWeights"]->Fill(scaleFactor);	
			
	/////////////////////////////
	// Selection
	/////////////////////////////

	//Declare selection instance    
	Selection selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...
      
	bool trigged = treeLoader.EventTrigged (itrigger);			
	bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);

      
/*      
      if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("XXX") == 0)
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
*/
	//Define object selection cuts
	selection.setJetCuts(20.,2.4,0.01,1.,0.98,0.3,0.1);
	selection.setDiElectronCuts(20,2.5,0.15,0.02,1);
	selection.setLooseElectronCuts(15,2.5,0.2);
	selection.setDiMuonCuts(20,2.4,0.15,10,0.02);
	selection.setLooseMuonCuts(10,2.5,0.2);
	  
	//Select objects 
	vector<TRootElectron*> selectedElectrons = selection.GetSelectedDiElectrons(vertex[0]);
	vector<TRootMuon*>     selectedMuons     = selection.GetSelectedDiMuons();
	vector<TRootJet*>      selectedJets      = selection.GetSelectedJets(20,2.4,true); // pt,eta,ApplyJetId
	sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class

	vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseElectrons(true); // VBTF Id
	vector<TRootMuon*>     looseMuons     = selection.GetSelectedLooseMuons();
	//vector<TRootMCParticle*> mcParticles;

	selecTableDiMu.Fill(d,1,scaleFactor);
	selecTableSemiEl.Fill(d,1,scaleFactor);
		
	//// EVENTS TRIGGERED BY MUON TRIGGER			

	if(!trigged)		   continue;
		selecTableDiMu.Fill(d,2,scaleFactor);

        if(!isGoodPV)   	   continue;
		selecTableDiMu.Fill(d,3,scaleFactor);

	if(selectedMuons.size()<2) continue;
		selecTableDiMu.Fill(d,4,scaleFactor);

  	bool foundZ = false;
	float windowsize = 20.;
  	for(unsigned int j=0;j<selectedMuons.size();j++)
  	{
  		for(unsigned int i=0;i<selectedMuons.size();i++)
  		{
   			TRootMuon* mu1 = (TRootMuon*) selectedMuons[j];
   			TRootMuon* mu2 = (TRootMuon*) selectedMuons[i];
    		if( fabs(mu2->Pt() - mu1->Pt()) < 0.001 || fabs(mu2->Eta() - mu1->Eta()) < 0.001 ) continue;
		if(mu1->charge() == mu2->charge()) continue;
      		double zMass = (*mu1 + *mu2).M();
      		if( zMass >= (91.-windowsize) && zMass <= (91.+windowsize) ) foundZ = true;
    		}
  	}

	if(!foundZ) continue; // opposite charge leptons with |mll-mz|<windowsize
		selecTableDiMu.Fill(d,5,scaleFactor);

	if(selectedJets.size()<1) continue; //at least 1 jet
		selecTableDiMu.Fill(d,6,scaleFactor); 

	if(selectedJets.size()<2) continue; //at least 2 jets
		selecTableDiMu.Fill(d,7,scaleFactor); 

	if(selectedJets.size()<3) continue; //at least 3 jets
		selecTableDiMu.Fill(d,8,scaleFactor); 

	MSPlot["MET"]->Fill(mets[0]->Et(),datasets[d], true, Luminosity*scaleFactor);
	if(mets[0]->Et()> 30.)    continue;
		selecTableDiMu.Fill(d,9,scaleFactor);

	MSPlot["NbOfLooseMuon"]->Fill(looseMuons.size(),datasets[d], true, Luminosity*scaleFactor);

//	MSPlot["NbOfLooseElectron"]->Fill(,datasets[d], true, Luminosity*scaleFactor);

// opposite charge leptons
//if(selectedMuons[0]->charge()== selectedMuons[1]->charge())
//require that there are no two electrons forming the Z mass
//if( selection.foundZCandidate(selectedMuons, selectedMuons, 20.) )
//it should not be an electron from a conversion!


      //////////////////////////////////////////////////////////////////////////
      // jet-parton matching
      //////////////////////////////////////////////////////////////////////////

	//delete selection;
    }//loop on events
    
    cout<<endl;
    
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

  
  // Fill the histograms
  
  //delete
  
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

