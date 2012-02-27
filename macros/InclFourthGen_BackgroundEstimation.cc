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

  string postfix = "_21Feb2012"; // to relabel the names of the output file  
	postfix= postfix+"_"+systematic;

  /////////////////////
  // Configuration
  /////////////////////
	string channelpostfix = "";
  bool semiElectron = false; // use semiElectron channel?
  bool semiMuon = true; // use semiMuon channel?
	if (argc >= 3)
	{	
	  semiMuon = atoi(argv[2]);
		semiElectron = !semiMuon;
	}
	
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

  //Output ROOT file
  string rootFileName ("InclFourthGenSearch_BackgroundEstimation"+postfix+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
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
  
 
  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootMET* > mets;

  //Global variable
  TRootEvent* event = 0;
	
  histo1D["LeptonPt_loose"] = new TH1F("leptonpt loose","leptonpt;pt leptons;#events",250,0,500);
  histo1D["LeptonPt_tight"] = new TH1F("leptonpt tight","leptonpt;pt leptons;#events",250,0,500);
	histo1D["LeptonReliso_loose"] = new TH1F("Lepton reliso loose", "Lepton reliso; reliso; # events", 50, 0, 1);
	histo1D["LeptonReliso_tight"] = new TH1F("Lepton reliso tight", "Lepton reliso; reliso; # events", 50, 0, 1);

	cout << " - Declared histograms ..." <<  endl;

  float NbSSevents = 0;   
  float NbTrievents = 0;  
	float Nb_Zpeak_EB_SS_MC = 0; float Nb_Zpeak_EE_SS_MC = 0; float Nb_Zpeak_EB_OS_MC = 0; float Nb_Zpeak_EE_OS_MC = 0;
	float Nb_Zpeak_EB_SS_data = 0; float Nb_Zpeak_EE_SS_data = 0; float Nb_Zpeak_EB_OS_data = 0; float Nb_Zpeak_EE_OS_data = 0;
	int NbOfLooseMuons = 0; int NbOfTightMuons = 0; int NbOfLooseElectrons = 0; int NbOfTightElectrons = 0; 
	

  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
  vector<string> CutsSelecTableChargeMisId_2El;
  CutsSelecTableChargeMisId_2El.push_back(string("Zpeak: SS el EB"));
  CutsSelecTableChargeMisId_2El.push_back(string("Zpeak: SS el EE"));
  CutsSelecTableChargeMisId_2El.push_back(string("Zpeak: OS el EB"));
  CutsSelecTableChargeMisId_2El.push_back(string("Zpeak: OS el EE"));

  SelectionTable selecTableChargeMisId_2El(CutsSelecTableChargeMisId_2El, datasets);
  selecTableChargeMisId_2El.SetLuminosity(Luminosity);
  selecTableChargeMisId_2El.SetPrecision(2);

  cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////////////////////////////////
  // PileUp Reweighting - 3D//
  ////////////////////////////////////////////////////
//  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
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
		
    /////////////////////////////////////
    /// Initialize JEC factors 
    /////////////////////////////////////
		
		//"OLD way, no JEC corrections on the fly"
    /*//L2L3 residual corrections already in data Toptrees now! (because a global tag is used where these corrections are included)
    vector<JetCorrectorParameters> vCorrParam;    
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false);
    */
   	    
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
 
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
    
    // true means redo also the L1
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);// last boolean ('startFromRaw') = true!
		


    ////////////////////////////////////
    ////////////////////////////////////
    ///////// Loop on events
    ////////////////////////////////////
    ////////////////////////////////////
    int itrigger = -1, previousRun = -1;
     
    if (verbose > 1)
      cout << " - Loop over events " << endl;      
    
    for (int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    {        

      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event ("<<100*ievt/datasets[d]->NofEvtsToRunOver()<<"%)"<<flush<<"\r";
      
			//load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"First cout after loading event:");
			
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
			//cout<<" currentFilename = "<<currentFilename<<", previousFilename = "<<previousFilename<<endl;
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
			//cout<<" currentRun = "<<currentRun<<", previousRun = "<<previousRun<<endl;
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
   					if(dataSetName == "ttW" || dataSetName == "ttZ" || dataSetName == "samesignWWjj" || dataSetName == "TTbarJets_scaleup" || dataSetName == "TTbarJets_scaledown" || dataSetName == "TTbarJets_matchingup" || dataSetName == "TTbarJets_matchingdown")
						  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_v5"), currentRun, iFile);//Summer11 MC! also the TTJets systematic samples...!
						else
						  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu30_eta2p1_v3"), currentRun, iFile);//Fall11 MC!
						
    
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
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v5"), currentRun, iFile);  				   					
						else if (event->runId() >= 178411 && event->runId() < 179942)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v8"), currentRun, iFile);  				   					
						else if (event->runId() >= 179942 && event->runId() <= 180296)
    					itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v9"), currentRun, iFile);  				   											   
  					if(itrigger == 9999)
						{
    					cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
    					exit(1);
  					}// semi-electron
 	   			}
	   			else 
	   			{
					  //Problem: a trigger reweighting procedure for MC should be done when using the summer11 electron trigger...
					  if(dataSetName == "ttW" || dataSetName == "ttZ" || dataSetName == "samesignWWjj" || dataSetName == "TTbarJets_scaleup" || dataSetName == "TTbarJets_scaledown" || dataSetName == "TTbarJets_matchingup" || dataSetName == "TTbarJets_matchingdown")
   						itrigger = treeLoader.iTrigger (string ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2"), currentRun, iFile);//Summer11 MC has other triggers!	
						else
						  itrigger = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v5"), currentRun, iFile);//Fall11 MC!
						
						if(itrigger == 9999)
						{
							cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;	
							exit(1);
						}
					}	 
				} //end if semiElectron	
			} //end previousRun != currentRun


//		cout << "bla 1" << endl;

			// JES CORRECTION   
      // Apply Jet Corrections on-the-fly: not if already in toptrees! (our first Fall11 round)
			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
//			if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
//				jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
//			else
//				jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
		  //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");

      //ordering is relevant; most probably 1) Type I MET correction, 2) JER where jet corrections are propagated to MET, 3) JES systematics where jet corrections are propagated to MET
      //----------------------------------------------------------
      // Apply type I MET corrections:  (Only for |eta| <= 4.7 )
      //---------------------------------------------------------
      
			//not if already in toptrees! (our first Fall11 round)
			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");      
//      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
//        jetTools->correctMETTypeOne(init_jets,mets[0],true);
//      else
//        jetTools->correctMETTypeOne(init_jets,mets[0],false);
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

//		cout << "bla 2" << endl;


			double lumiWeight3D = 1.0;
			if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA"))
			{				
      	////////////////////////////
      	// apply PU Reweighting
      	////////////////////////////
				lumiWeight3D = Lumi3DWeights.weight3D(event->nPu(-1),event->nPu(0),event->nPu(+1));
	 			scaleFactor = scaleFactor*lumiWeight3D;
      	//histo1D["lumiWeights"]->Fill(scaleFactor);	
			}
						
//		cout << "bla 3" << endl;
								
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
			
//		cout << "bla 4" << endl;
      bool eventSelected = false;
      
      vector<TRootJet*> selectedJets;//selectedJetsFromW,selectedJetsFromW_DropUsedJets,selectedJetsFromW_DropUsedJets_tmp;
      vector<TRootJet*> selectedForwardJets, selectedJetsLargeEtaRange;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootElectron*> selectedLooseElectronsNoVBTFid;
      vector<TRootElectron*> selectedLooseElectronsVBTFid;
      vector<TRootMuon*> selectedLooseMuons;
      vector<TRootMCParticle*> mcParticles;
      vector<TRootMuon*> selectedFakeMuons;
      vector<TRootElectron*> selectedFakeElectrons;
//		cout << "bla 5" << endl;


      //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before treeLoader.LoadMCEvent:");      
      if(dataSetName.find("TTbarJets_SemiMu") == 0 || dataSetName.find("TTbarJets_SemiElectron") == 0 || dataSetName.find("NP_Tprime")==0 || dataSetName.find("NP_overlay_Tprime")==0)
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
			//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After treeLoader.LoadMCEvent:");
//		cout << "bla 6" << endl;

      
      float METCut = 40;
      selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1);
      selection.setMuonCuts(20,2.1,0.125,10,0.02,0.3,1,1,1);
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setElectronCuts(20,2.5,0.1,0.02,1,0.3);
      selection.setLooseElectronCuts(15,2.5,0.2);	 				
      
//		cout << "bla 7" << endl;
      if (init_jets.size() > 0)
      {
	    	selectedJets = selection.GetSelectedJets(true);				
	    	selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	    	selectedElectrons = selection.GetSelectedElectrons(vertex[0],selectedJets);
      }
      selectedLooseElectronsNoVBTFid = selection.GetSelectedLooseElectrons(false); //no vbtfid is required
      selectedLooseElectronsVBTFid = selection.GetSelectedLooseElectrons(true); //loose vbtfid is required 
      selectedLooseMuons = selection.GetSelectedLooseMuons(); //veto muons	
//		cout << "bla 8" << endl;
      
      selection.setMuonCuts(10,2.5,0.2,0,99999.,1.0,0,1,0);
      selectedFakeMuons = selection.GetSelectedMuons(); //muons that could fake tight muons
      selection.setLooseElectronCuts(15,2.5,0.2);
      selectedFakeElectrons = selection.GetSelectedLooseElectrons(false); //electrons that could fake tight electrons

			vector<TRootMuon*> selectedMuons_semiMuon, selectedFakeMuons_semiMuon;
			vector<TRootElectron*> selectedElectrons_semiElectron, selectedFakeElectrons_semiElectron;
			
//		cout << "test1 " << endl;
			if(semiMuon){
				for(unsigned int j=0;j<selectedMuons.size();j++)
					if(j>0) selectedMuons_semiMuon.push_back(selectedMuons[j]);
 				for(unsigned int i=0;i<selectedFakeMuons.size();i++)
				{
					if(selectedMuons.size()>0)
					{
						TRootMuon* mu1 = (TRootMuon*) selectedMuons[0];
						TRootMuon* mu2 = (TRootMuon*) selectedFakeMuons[i];
						if( fabs(mu2->Pt() - mu1->Pt()) > 0.001 && fabs(mu2->Eta() - mu1->Eta()) > 0.001){
							selectedFakeMuons_semiMuon.push_back(selectedFakeMuons[i]);
						}
					}else
						selectedFakeMuons_semiMuon.push_back(selectedFakeMuons[i]);
				}
			}
			
//		cout << "test2 " << endl;
			if(semiElectron){
				for(unsigned int j=0;j<selectedElectrons.size();j++)
					if(j>0) selectedElectrons_semiElectron.push_back(selectedElectrons[j]);			
 				for(unsigned int i=0;i<selectedFakeElectrons.size();i++)
				{
					if(selectedElectrons.size()>0)
					{
						TRootElectron* el1 = (TRootElectron*) selectedElectrons[0];
						TRootElectron* el2 = (TRootElectron*) selectedFakeElectrons[i];
						if( fabs(el2->Pt() - el1->Pt()) > 0.001 && fabs(el2->Eta() - el1->Eta()) > 0.001){
							selectedFakeElectrons_semiElectron.push_back(selectedFakeElectrons[i]);
						}
					}else
						selectedFakeElectrons_semiElectron.push_back(selectedFakeElectrons[i]);
					
				}
			}
			
//		cout << "test3 " << endl;
			//// EVENTS TRIGGERED BY MUON TRIGGER			
      if(trigged && semiMuon)
      { 
        if(isGoodPV)
				{
					if(selectedMuons.size()>=1 && selectedMuons[0]->Pt()>40)
					{
						sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class
						
						if(selectedJets.size()>=(unsigned int)anaEnv.NofJets)
						{  //at least 1 jet!
						
							////////////////////// FAKE LEPTON RATE ESTIMATION ///////////////////
									//if(dataSetName == "Data")
									//{
  									bool foundZ = false;
  									if(selectedFakeMuons.size()>=1){
											for(unsigned int j=0;j<selectedFakeMuons.size();j++)
  										{
  											for(unsigned int i=0;i<selectedFakeMuons.size();i++)
										 		{
   												TRootMuon* mu1 = (TRootMuon*) selectedFakeMuons[j];
   												TRootMuon* mu2 = (TRootMuon*) selectedFakeMuons[i];
    											if( fabs(mu2->Pt() - mu1->Pt()) > 0.001 && fabs(mu2->Eta() - mu1->Eta()) > 0.001 && mu1->charge() != mu2->charge())
    											{	
      											if( (*mu1 + *mu2).M() >= (91.-20) && (*mu1 + *mu2).M() <= (91.+20) )
        											foundZ = true;
													}
										    }
										  }
										}
 										if(selectedFakeElectrons.size()>=1)
										{
  										for(unsigned int j=0;j<selectedFakeElectrons.size();j++)
  										{
  											for(unsigned int i=0;i<selectedFakeElectrons.size();i++)
										 		{
   												TRootElectron* el1 = (TRootElectron*) selectedFakeElectrons[j];
   												TRootElectron* el2 = (TRootElectron*) selectedFakeElectrons[i];
    											if( fabs(el2->Pt() - el1->Pt()) > 0.001 && fabs(el2->Eta() - el1->Eta()) > 0.001 && el1->charge() != el2->charge())
    											{	
      											if( (*el1 + *el2).M() >= (91.-20) && (*el1 + *el2).M() <= (91.+20) )
        											foundZ = true;
													}
										    }
										  }
										}
//		cout << "test6 " << endl;
										
										if(mets[0]->Et() < 20 && !foundZ) // reject events with MET>20 and with a Z boson
										{											
											float MT = 999999.;
											if(selectedFakeMuons_semiMuon.size()>=1)
											{
												//cout << "selectedFakeMuons.size() "<< selectedFakeMuons.size() << endl;
												//cout << "selectedMuons.size() "<< selectedMuons.size() << endl;
												for(unsigned int j=0;j<selectedFakeMuons_semiMuon.size();j++){
													MT = sqrt(pow(mets[0]->Et()+selectedFakeMuons_semiMuon[j]->Et(),2)-pow(mets[0]->Pt()+selectedFakeMuons_semiMuon[j]->Pt(),2));
													//if(selectedFakeMuons[j]->Pt()<35 && 
													//if(MT<25) 
														NbOfLooseMuons++;
												}
												for(unsigned int j=0;j<selectedMuons_semiMuon.size();j++){
													MT = sqrt(pow(mets[0]->Et()+selectedMuons_semiMuon[j]->Et(),2)-pow(mets[0]->Pt()+selectedMuons_semiMuon[j]->Pt(),2));
													//if(selectedMuons[j]->Pt()<35 && 	
													//if(MT<25) 
														NbOfTightMuons++;
												}
											}
//		cout << "test7 " << endl;
											if(selectedFakeElectrons_semiElectron.size()>=1)
											{
												//cout << "selectedFakeElectrons_semiElectron.size() "<< selectedFakeElectrons_semiElectron.size() << endl;
												//cout << "selectedElectrons_semiElectron.size() "<< selectedElectrons_semiElectron.size() << endl;
												float MT = 999999.;
												for(unsigned int j=0;j<selectedFakeElectrons_semiElectron.size();j++){
													MT = sqrt(pow(mets[0]->Et()+selectedFakeElectrons[j]->Et(),2)-pow(mets[0]->Pt()+selectedFakeElectrons[j]->Pt(),2));
													//if(selectedFakeElectrons[j]->Pt()<55 && 
													//if(MT<25)
														NbOfLooseElectrons++;      
												}
												float relIso = (selectedFakeElectrons_semiElectron[0]->chargedHadronIso()+selectedFakeElectrons_semiElectron[0]->neutralHadronIso()+selectedFakeElectrons_semiElectron[0]->photonIso())/selectedFakeElectrons_semiElectron[0]->Pt();
												for(unsigned int j=0;j<selectedElectrons_semiElectron.size();j++){
													MT = sqrt(pow(mets[0]->Et()+selectedElectrons_semiElectron[j]->Et(),2)-pow(mets[0]->Pt()+selectedElectrons_semiElectron[j]->Pt(),2));
													//if(selectedElectrons_semiElectron[j]->Pt()<55 && 
													//if(MT<25)
														NbOfTightElectrons++;
												}
											}
//		cout << "test8 " << endl;
										}
									//}
									////////////////////// FAKE LEPTON RATE ESTIMATION ///////////////////
						} //end 'at least one jet'  
          } // end if selectedMuons.size()>=1
        } // end good PV
      }// end trigged & semiMuon
      
			///// EVENTS TRIGGERED BY ELECTRON TRIGGER
			else if(trigged && semiElectron)
      {
//		cout << "test4 " << endl;
    	 	
        if( isGoodPV )
        {
          if( selectedElectrons.size() >= 1 && selectedElectrons[0]->Pt()>40)
          {
              if( selection.passConversionRejection(selectedElectrons[0]) )
              {
								sort(selectedJets.begin(),selectedJets.end(),HighestPt()); // HighestPt() is included from the Selection class

								if( selectedJets.size()>=(unsigned int)anaEnv.NofJets)
								{

									////////////////////// CHARGE MIS-ID RATE (NO B-TAG CUT & MET CUT!!!) ///////////////////////
 									if(selectedElectrons.size() == 2 && selectedLooseElectronsVBTFid.size() == selectedElectrons.size())
									{
										if(selection.passConversionRejection(selectedElectrons[1]))
										{
											if(selection.foundZCandidate(selectedElectrons, selectedElectrons, 10.))
											{								
												if(selectedElectrons[0]->charge()== selectedElectrons[1]->charge())
												{ 
													if(fabs(selectedElectrons[0]->superClusterEta())<1.4442 && fabs(selectedElectrons[1]->superClusterEta())<1.4442){
														selecTableChargeMisId_2El.Fill(d,0,scaleFactor);
														if(dataSetName.find("Data") == 0) Nb_Zpeak_EB_SS_data+=scaleFactor;
														else Nb_Zpeak_EB_SS_MC+=scaleFactor;
													}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660 && fabs(selectedElectrons[1]->superClusterEta())>1.5660){
														selecTableChargeMisId_2El.Fill(d,1,scaleFactor);
														if(dataSetName.find("Data") == 0) Nb_Zpeak_EE_SS_data+=scaleFactor;
														else Nb_Zpeak_EE_SS_MC+=scaleFactor;
													}
												}else{
													if(fabs(selectedElectrons[0]->superClusterEta())<1.4442 && fabs(selectedElectrons[1]->superClusterEta())<1.4442){
														selecTableChargeMisId_2El.Fill(d,2,scaleFactor);
														if(dataSetName.find("Data") == 0) Nb_Zpeak_EB_OS_data+=scaleFactor;
														else Nb_Zpeak_EB_OS_MC+=scaleFactor;
													}else if(fabs(selectedElectrons[0]->superClusterEta())>1.5660 && fabs(selectedElectrons[1]->superClusterEta())>1.5660){
														selecTableChargeMisId_2El.Fill(d,3,scaleFactor);
														if(dataSetName.find("Data") == 0) Nb_Zpeak_EE_OS_data+=scaleFactor;
														else Nb_Zpeak_EE_OS_MC+=scaleFactor;
													}															
												}
											}	
										}
									}
									////////////////////// CHARGE MIS-ID RATE ----- END //////////////////
								
//		cout << "test5 " << endl;
								
									////////////////////// FAKE LEPTON RATE ESTIMATION ///////////////////
									//if(dataSetName == "Data")
									//{
  									bool foundZ = false;
  									if(selectedFakeMuons.size()>=1){
											for(unsigned int j=0;j<selectedFakeMuons.size();j++)
  										{
  											for(unsigned int i=0;i<selectedFakeMuons.size();i++)
										 		{
   												TRootMuon* mu1 = (TRootMuon*) selectedFakeMuons[j];
   												TRootMuon* mu2 = (TRootMuon*) selectedFakeMuons[i];
    											if( fabs(mu2->Pt() - mu1->Pt()) > 0.001 && fabs(mu2->Eta() - mu1->Eta()) > 0.001 && mu1->charge() != mu2->charge())
    											{	
      											if( (*mu1 + *mu2).M() >= (91.-20) && (*mu1 + *mu2).M() <= (91.+20) )
        											foundZ = true;
													}
										    }
										  }
										}
 										if(selectedFakeElectrons.size()>=1)
										{
  										for(unsigned int j=0;j<selectedFakeElectrons.size();j++)
  										{
  											for(unsigned int i=0;i<selectedFakeElectrons.size();i++)
										 		{
   												TRootElectron* el1 = (TRootElectron*) selectedFakeElectrons[j];
   												TRootElectron* el2 = (TRootElectron*) selectedFakeElectrons[i];
    											if( fabs(el2->Pt() - el1->Pt()) > 0.001 && fabs(el2->Eta() - el1->Eta()) > 0.001 && el1->charge() != el2->charge())
    											{	
      											if( (*el1 + *el2).M() >= (91.-20) && (*el1 + *el2).M() <= (91.+20) )
        											foundZ = true;
													}
										    }
										  }
										}
//		cout << "test6 " << endl;
										
										if(mets[0]->Et() < 20 && !foundZ) // reject events with MET>20 and with a Z boson
										{											
											float MT = 999999.;
											if(selectedFakeMuons.size()>=1)
											{
												cout << "selectedFakeMuons.size() "<< selectedFakeMuons.size() << endl;
												cout << "selectedMuons.size() "<< selectedMuons.size() << endl;
												for(unsigned int j=0;j<selectedFakeMuons.size();j++){
													MT = sqrt(pow(mets[0]->Et()+selectedFakeMuons[j]->Et(),2)-pow(mets[0]->Pt()+selectedFakeMuons[j]->Pt(),2));
													//if(selectedFakeMuons[j]->Pt()<35 && 
													//if(MT<25) 
														NbOfLooseMuons++;
												}
												for(unsigned int j=0;j<selectedMuons.size();j++){
													MT = sqrt(pow(mets[0]->Et()+selectedMuons[j]->Et(),2)-pow(mets[0]->Pt()+selectedMuons[j]->Pt(),2));
													//if(selectedMuons[j]->Pt()<35 && 	
													//if(MT<25) 
														NbOfTightMuons++;
												}
											}
//		cout << "test7 " << endl;
											if(selectedFakeElectrons_semiElectron.size()>=1)
											{
												cout << "selectedFakeElectrons_semiElectron.size() "<< selectedFakeElectrons_semiElectron.size() << endl;
												cout << "selectedElectrons_semiElectron.size() "<< selectedElectrons_semiElectron.size() << endl;
												float MT = 999999.;
												for(unsigned int j=0;j<selectedFakeElectrons_semiElectron.size();j++){
													MT = sqrt(pow(mets[0]->Et()+selectedFakeElectrons[j]->Et(),2)-pow(mets[0]->Pt()+selectedFakeElectrons[j]->Pt(),2));
													//if(selectedFakeElectrons[j]->Pt()<55 && 
													//if(MT<25)
														NbOfLooseElectrons++;      
												}
												histo1D["LeptonPt_loose"]->Fill(selectedFakeElectrons_semiElectron[0]->Pt(),scaleFactor);	
												float relIso = (selectedFakeElectrons_semiElectron[0]->chargedHadronIso()+selectedFakeElectrons_semiElectron[0]->neutralHadronIso()+selectedFakeElectrons_semiElectron[0]->photonIso())/selectedFakeElectrons_semiElectron[0]->Pt();
												histo1D["LeptonReliso_loose"]->Fill(relIso,scaleFactor);	
												for(unsigned int j=0;j<selectedElectrons_semiElectron.size();j++){
													MT = sqrt(pow(mets[0]->Et()+selectedElectrons_semiElectron[j]->Et(),2)-pow(mets[0]->Pt()+selectedElectrons_semiElectron[j]->Pt(),2));
													//if(selectedElectrons_semiElectron[j]->Pt()<55 && 
													//if(MT<25)
														NbOfTightElectrons++;
												}
												if(selectedElectrons_semiElectron.size()>=1){
													histo1D["LeptonPt_tight"]->Fill(selectedElectrons_semiElectron[0]->Pt(),scaleFactor);	
													relIso = (selectedElectrons_semiElectron[0]->chargedHadronIso()+selectedElectrons_semiElectron[0]->neutralHadronIso()+selectedElectrons_semiElectron[0]->photonIso())/selectedElectrons_semiElectron[0]->Pt();
												histo1D["LeptonReliso_tight"]->Fill(relIso,scaleFactor);	
												}
											}
//		cout << "test8 " << endl;
										}
									//}
									////////////////////// FAKE LEPTON RATE ESTIMATION ///////////////////
								} // end 'at least one jet'
							} // end conversion rejection for leading electron
          } // end if selectedElectrons.size()>=1
        } // end good PV
      } // end trigged & semiElectron
						
    
	
    }//loop on events
    
    cout<<endl;
		
    if(jetTools) delete jetTools;
		
		//important: free memory
    treeLoader.UnLoadDataset();
	
    
  } //loop on datasets

	if(semiElectron){
		float chargeMisId_Barrel_MC = (float)Nb_Zpeak_EB_SS_MC/(2*(float)Nb_Zpeak_EB_OS_MC);
		float chargeMisId_Barrel_data = (float)Nb_Zpeak_EB_SS_data/(2*(float)Nb_Zpeak_EB_OS_data);
		float chargeMisId_Endcap_MC = (float)Nb_Zpeak_EE_SS_MC/(2*(float)Nb_Zpeak_EE_OS_MC);
		float chargeMisId_Endcap_data = (float)Nb_Zpeak_EE_SS_data/(2*(float)Nb_Zpeak_EE_OS_data);
  
		cout << "chargeMisId_Barrel_MC " << chargeMisId_Barrel_MC << endl;
		cout << "chargeMisId_Barrel_data " << chargeMisId_Barrel_data << endl;
		cout << "chargeMisId_Endcap_MC " << chargeMisId_Endcap_MC << endl;
		cout << "chargeMisId_Endcap_data " << chargeMisId_Endcap_data << endl;
	
		float MuonFakeRate = (float)NbOfTightMuons/(float)(NbOfLooseMuons);
		float ElectronFakeRate = (float)NbOfTightElectrons/(float)(NbOfLooseElectrons);
				
		cout << "Number of tight muons " << NbOfTightMuons << endl;
		cout << "Number of loose muons " << NbOfLooseMuons << endl;
		cout << "Fake rate for muons " << MuonFakeRate << endl;
		cout << "Number of tight electrons " << NbOfTightElectrons << endl;
		cout << "Number of loose electrons " << NbOfLooseElectrons << endl;
		cout << "Fake rate for electrons " << ElectronFakeRate << endl;
	}

  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;	
  selecTableChargeMisId_2El.TableCalculator(true, true, true, true, true, false, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
  string selectiontableChargeMisId_2El = "InclFourthGenSearch_SelectionTable_ChargeMisId2El"+postfix;
  selectiontableChargeMisId_2El = selectiontableChargeMisId_2El +".tex"; 	
	if(semiElectron) selecTableChargeMisId_2El.Write(selectiontableChargeMisId_2El.c_str(),false, true, false, false, false, false, false);
 
    fout->cd();
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
			tempCanvas->SaveAs( (it->first+".pdf").c_str() ); //well, is actually not png but pdf...
    }    
    cout << "1D plots written" << endl;
    fout->Close();
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
     for(unsigned int k=0; k<init_jets.size(); k++) //init_jets.size()
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
