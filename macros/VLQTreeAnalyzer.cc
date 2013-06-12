//////////////////////////////////////////////////
///// Based on InclFourthGen_TreeAnalyzer.cc /////
//////////////////////////////////////////////////


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
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysis/VLQSearch/interface/VLQTree.h"

#include "Style.C"

using namespace std;
using namespace TopTree;
using namespace RooFit;
using namespace reweight;

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

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

unsigned int nBjetsPresent(vector < float > btagvalues, float btagWP);

unsigned int nNonBjetsPresent(vector < float > btagvalues, float btagWP);

int main (int argc, char *argv[])
{	
	cout<<"The arguments passed to the executable are: "<<endl;
  for(unsigned int i=1;i<argc;i++)
	{
		cout<<argv[i]<<endl;
	}
	
  //which systematic to run?
  string systematic = "Nominal";
	string systematicnew = "";
  if (argc >= 2)
		systematic = string(argv[1]);
  cout << "Systematic to be used: " << systematic << endl;
  if( ! (systematic == "Nominal"  || systematic == "JESPlus" || systematic == "JESMinus" || systematic == "JERPlus" || systematic == "JERMinus" || systematic == "bTagPlus" || systematic == "bTagMinus" || systematic == "misTagPlus" || systematic == "misTagMinus" || systematic == "PUPlus" || systematic == "PUMinus"))
  {
    cout << "Unknown systematic!!!" << endl;
    cout << "Possible options are: Nominal JESPlus JESMinus JERPlus JERMinus bTagPlus bTagMinus misTagPlus misTagMinus PUPlus PUMinus" << endl;
    exit(-1);
  }	  

  //maybe not used in analysis at all
  string btagger = "CSVM", antibtagger = "CSVM";
  float btagWP = 9999, antibtagWP = 9999; //{1.7,3.3,10.2}; trackcountinghighefficiency working points: loose, medium, tight
  if(btagger == "CSVL")
     btagWP = 0.244 ;
  else if(btagger == "CSVM")
     btagWP = 0.679;
  else if(btagger == "CSVT")
     btagWP = 0.898;
	if(antibtagger == "CSVL")
     antibtagWP = 0.244 ;
  else if(antibtagger == "CSVM")
     antibtagWP = 0.679;
  else if(antibtagger == "CSVT")
     antibtagWP = 0.898;

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the VLQ search ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  setTDRStyle();
  //setMyStyle();

  string inputpostfixOld = ""; // should be same as postfix in TreeCreator of the trees
	string inputpostfix= inputpostfixOld+"_"+systematic;		

  string Treespath = "VLQTrees_Summer12_PBS";
  Treespath = Treespath + "/"; 		
	bool savePNG = false;
	string outputpostfix = "";
	outputpostfix = outputpostfix+"_"+systematic;
	string Outputpath = "OutputFiles_VLQAnalyzer_29May2013_withoutQCDestimation_leadingjetcut_mutabletest";
	Outputpath = Outputpath + "/";
	mkdir(Outputpath.c_str(),0777);
		
	
  /////////////////////
  // Configuration
  /////////////////////

  bool doPUreweighting = true;
	
	//bool make2Dbinning = false;
	//if (argc >= 6)
  //	make2Dbinning = atoi(argv[5]);
	//bool SplitTTbarSample_for2Dbinning = false; //split ttbar sample in two pieces to make the 2D binning and to do the evalution
	//bool UseBinningTTbarSample_forEvaluation = false; //i.e. use for the evaluation the piece of the ttbar sample that is used to create the 2D binning 
  
  //xml file
  string xmlFileName = "";
	xmlFileName = "../config/myVLQConfig.xml";	
  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;   
	
	float LeadingJetPtCut = 200., SubLeadingJetPtCut = 100.; //subleading jet pt cut only for the pair production signal boxes...
	bool doBtagging = true; //Higgs->bbar decay
	bool doAntiBtagging = true; //all other jets
  
  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = anaEnv.Verbose;
  float anaEnvLuminosity = anaEnv.Luminosity;	// in 1/pb 
  cout << "analysis environment luminosity for rescaling (if not run on data) "<< anaEnvLuminosity << endl;

  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
	vector < Dataset* > datasetsMu;
  vector < Dataset* > datasetsEl;
  vector < Dataset* > datasetsPlot;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  
  float Luminosity = anaEnvLuminosity;

  float LuminosityMu = anaEnvLuminosity;
  float LuminosityEl = anaEnvLuminosity;

  bool foundMu = false;
  bool foundEl = false;
	
	vector<string> inputTrees; //fill with the tree files you want to read!!!
	for (unsigned int d = 0; d < datasets.size(); d++) //d < datasets.size()
  {				
    cout << "   Dataset " << d << " name : " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;    
    string dataSetName = datasets[d]->Name();
		
		if(dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
      LuminosityMu = datasets[d]->EquivalentLumi();
      foundMu=true;
    }  
    if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) {
      LuminosityEl = datasets[d]->EquivalentLumi();
      foundEl=true;  
    } 
		
		if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // when only using one channel
      Luminosity = datasets[d]->EquivalentLumi();
    } 
		
		string inputTreeFileName; //should follow convention of TreeFileName in VLQTreeCreator.cc
		if(systematic == "JESPlus" || systematic == "JESMinus" || systematic == "JERPlus" || systematic == "JERMinus")
		{
				if ( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) )
				{
					inputTreeFileName = Treespath+"VLQTree_"+dataSetName+inputpostfixOld+"_Nominal"+".root"; //is actually dummy
					cout<<"  Running systematics: data will be skipped later on"<<endl;
				}
				else
					inputTreeFileName = Treespath+"VLQTree_"+dataSetName+inputpostfix+".root";		
		}
		else
		{
		    ////TEMPORARY 
				//if( (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso_Mu") != string::npos))
				//{
				   inputTreeFileName = Treespath+"VLQTree_"+dataSetName+inputpostfixOld+".root";
				//}
				//else inputTreeFileName = Treespath+"VLQTree_"+dataSetName+inputpostfixOld+"_Nominal"+".root";
		}
		inputTrees.push_back(inputTreeFileName);
	}
	
	if(!foundMu && !foundEl && Luminosity != anaEnvLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  else {
    if(LuminosityMu != anaEnvLuminosity) cout << "Muon PD: changed analysis environment luminosity to "<< LuminosityMu << endl;
    if(LuminosityEl != anaEnvLuminosity) cout << "Electron PD: changed analysis environment luminosity to "<< LuminosityEl << endl;
  }
	
	// make a datasets vector only for 
  if (foundMu) {
      for (unsigned int d = 0; d < datasets.size (); d++) {
	string dataSetName = datasets[d]->Name();
	if ( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) && !(dataSetName.find("InvIso_El") != string::npos)) 
	  datasetsMu.push_back(datasets[d]);
      }
  }
      
  if (foundEl) {
    for (unsigned int d = 0; d < datasets.size (); d++) {
      string dataSetName = datasets[d]->Name();
      if ( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) && !(dataSetName.find("InvIso_Mu") != string::npos)) 
	datasetsEl.push_back(datasets[d]);
    }
  }
	
  //Output ROOT file
  string rootFileName (Outputpath+"VLQTreeAnalyzer"+outputpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //Global variable

  string pathPNG = "VLQSearchPlots_TreeAnalyzer"+inputpostfix;
  pathPNG = pathPNG +"/"; 		
  if(savePNG) mkdir(pathPNG.c_str(),0777);
  
  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);

  
  /*string multileptons[7] = {"SSElEl","SSElMu","SSMuMu","ElElEl","ElElMu","ElMuMu","MuMuMu"};
  string histoName,histo_dataset;
  for(int i = 0; i<7; i++)
  {
		histoName = "NbEvents_"+multileptons[i];
		for(unsigned int d = 0; d < datasets.size (); d++){
			histo_dataset = histoName+(datasets[d]->Name()).c_str(); 
			histo1D[histo_dataset.c_str()] = new TH1F(histo_dataset.c_str(),histo_dataset.c_str(), 1, 0.5, 1.5);
		}
  }*/
	
  /*MSPlot["MS_NbSSevents"] = new MultiSamplePlot(datasets,"# events with SS leptons", 1, 0.5, 1.5, "");
  MSPlot["MS_NbSSElElevents"] = new MultiSamplePlot(datasets,"# events with SS electrons", 1, 0.5, 1.5, "");
  MSPlot["MS_NbSSElMuevents"] = new MultiSamplePlot(datasets,"# events with SS electron+muon", 1, 0.5, 1.5, "");
  MSPlot["MS_NbSSMuMuevents"] = new MultiSamplePlot(datasets,"# events with SS muons", 1, 0.5, 1.5, "");
  MSPlot["MS_NbTrievents"] = new MultiSamplePlot(datasets,"# events with 3 leptons", 1, 0.5, 1.5, "");
  MSPlot["MS_MET"] = new MultiSamplePlot(datasets,"MET", 75, 0, 200, "Missing transverse energy (GeV)");
  MSPlot["MS_LeptonPt"] = new MultiSamplePlot(datasets,"lepton pt", 75, 0, 250, "Pt lepton (GeV)");
  MSPlot["MS_nPV_beforePUreweighting"] = new MultiSamplePlot(datasets, "nPrimaryVertices_beforePUreweighting", 21, -0.5, 20.5, "Nr. of primary vertices");
	MSPlot["MS_nPV"] = new MultiSamplePlot(datasets, "nPrimaryVertices", 21, -0.5, 20.5, "Nr. of primary vertices");
  MSPlot["MS_JetMultiplicity_SingleLepton"] = new MultiSamplePlot(datasets, "JetMultiplicity", 10, -0.5, 9.5, "Jet Multiplicity");
  MSPlot["MS_BtaggedJetMultiplicity_SingleLepton"] = new MultiSamplePlot(datasets, "BtaggedJetMultiplicity", 7, -0.5, 6.5, "b-tagged jet multiplicity");
	MSPlot["MS_JetMultiplicityAtleast1Btag_SingleLepton"] = new MultiSamplePlot(datasets, "JetMultiplicityAtleast1Btag", 10, -0.5, 9.5, "Jet multiplicity (>=1 b-tag)");
  MSPlot["MS_Reliso_Lepton"] = new MultiSamplePlot(datasets, "Lepton reliso", 40, 0, 0.25, "Lepton reliso");

  MSPlot["MS_JetPt_all_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_all", 50, 0, 300, "Pt of all jets (GeV)");
	MSPlot["MS_JetPt_btagged_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_btagged", 50, 0, 300, "Pt of b-tagged jets (GeV)");
	MSPlot["MS_JetPt_nonbtagged_SingleLepton"] = new MultiSamplePlot(datasets,"JetPt_nonbtagged", 50, 0, 300, "Pt of non b-tagged jets (GeV)");
  
	MSPlot["MS_MET_nobtag"] = new MultiSamplePlot(datasets,"MET", 75, 0, 200, "Missing transverse energy (GeV)");
  MSPlot["MS_LeptonPt_nobtag"] = new MultiSamplePlot(datasets,"lepton pt", 75, 0, 250, "Pt lepton (GeV)");*/
	
	
	
	
	//'mu channel' = when at least one muon is trigged 
	MSPlot["MS_nPV_noPUreweighting_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices before PU reweighting", 31, -0.5, 30.5, "Nr. of primary vertices");
  MSPlot["MS_nPV_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices", 31, -0.5, 30.5, "Nr. of primary vertices");  //after PU reweighting
	MSPlot["MS_nPV_noPUreweighting_4jets_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices before PU reweighting", 31, -0.5, 30.5, "Nr. of primary vertices");
  MSPlot["MS_nPV_4jets_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices", 31, -0.5, 30.5, "Nr. of primary vertices");  //after PU reweighting
	MSPlot["MS_nLeptons_mu"] = new MultiSamplePlot(datasetsMu,"number of leptons", 6, -0.5, 5.5, "number of leptons");
	MSPlot["MS_nMuons_mu"] = new MultiSamplePlot(datasetsMu,"number of muons", 6, -0.5, 5.5, "number of muons");
	MSPlot["MS_nElectrons_mu"] = new MultiSamplePlot(datasetsMu,"number of electrons", 6, -0.5, 5.5, "number of electrons");
	MSPlot["MS_MET_mu"] = new MultiSamplePlot(datasetsMu,"MET", 200, 0, 1000, "Missing transverse energy (GeV)");
  MSPlot["MS_MET_singleMu_mu"] = new MultiSamplePlot(datasetsMu,"MET", 200, 0, 1000, "Missing transverse energy (GeV)");
	//jets (not forward)
	MSPlot["MS_JetMultiplicity_mu"] = new MultiSamplePlot(datasetsMu, "JetMultiplicity", 10, -0.5, 9.5, "jet multiplicity");
	MSPlot["MS_Pt_alljets_mu"] = new MultiSamplePlot(datasetsMu,"Pt_alljets", 50, 0, 300, "Pt of all jets (GeV)");
	MSPlot["MS_Pt_jet1_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet1", 50, 0, 300, "Pt of jet 1 (GeV)");
	MSPlot["MS_Pt_jet2_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet2", 50, 0, 300, "Pt of jet 2 (GeV)");
	MSPlot["MS_Pt_jet3_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet3", 50, 0, 300, "Pt of jet 3 (GeV)");
	MSPlot["MS_Pt_jet4_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet4", 50, 0, 300, "Pt of jet 4 (GeV)");
	MSPlot["MS_Eta_alljets_mu"] = new MultiSamplePlot(datasetsMu,"Eta_alljets", 50, -2.5 , 2.5, "Eta of all jets");
	MSPlot["MS_Eta_jet1_mu"] = new MultiSamplePlot(datasetsMu,"Eta_jet1", 50, -2.5 , 2.5, "Eta of jet 1");
	MSPlot["MS_Eta_jet2_mu"] = new MultiSamplePlot(datasetsMu,"Eta_jet2", 50, -2.5 , 2.5, "Eta of jet 2");
	MSPlot["MS_Eta_jet3_mu"] = new MultiSamplePlot(datasetsMu,"Eta_jet3", 50, -2.5 , 2.5, "Eta of jet 3");
	MSPlot["MS_Eta_jet4_mu"] = new MultiSamplePlot(datasetsMu,"Eta_jet4", 50, -2.5 , 2.5, "Eta of jet 4");
	//forward jets
	MSPlot["MS_ForwardJetMultiplicity_mu"] = new MultiSamplePlot(datasetsMu, "JetMultiplicity", 10, -0.5, 9.5, "Forward jet multiplicity");
	MSPlot["MS_Pt_allForwardjets_mu"] = new MultiSamplePlot(datasetsMu,"Pt_alljets", 50, 0, 300, "Pt of all forward jets (GeV)");
	MSPlot["MS_Pt_Forwardjet1_mu"] = new MultiSamplePlot(datasetsMu,"Pt_Forwardjet1", 50, 0, 300, "Pt of forward jet 1 (GeV)");
	MSPlot["MS_Pt_Forwardjet2_mu"] = new MultiSamplePlot(datasetsMu,"Pt_Forwardjet2", 50, 0, 300, "Pt of forward jet 2 (GeV)");
	MSPlot["MS_Eta_allForwardjets_mu"] = new MultiSamplePlot(datasetsMu,"Eta_allForwardjets", 100, -5 , 5, "Eta of all forward jets");
	MSPlot["MS_Eta_Forwardjet1_mu"] = new MultiSamplePlot(datasetsMu,"Eta_Forwardjet1", 100, -5 , 5, "Eta of forward jet 1");
	MSPlot["MS_Eta_Forwardjet2_mu"] = new MultiSamplePlot(datasetsMu,"Eta_Forwardjet2", 100, -5 , 5, "Eta of forward jet 2");
	//leptons
  MSPlot["MS_Pt_allMu_mu"] = new MultiSamplePlot(datasetsMu,"muon pt", 75, 0, 250, "Pt of all muons (GeV)");
	MSPlot["MS_Pt_mu1_mu"] = new MultiSamplePlot(datasetsMu,"Pt_mu1", 50, 0, 300, "Pt of muon 1 (GeV)");
	MSPlot["MS_Pt_mu2_mu"] = new MultiSamplePlot(datasetsMu,"Pt_mu2", 50, 0, 300, "Pt of muon 2 (GeV)");
	MSPlot["MS_Pt_mu3_mu"] = new MultiSamplePlot(datasetsMu,"Pt_mu3", 50, 0, 300, "Pt of muon 3 (GeV)");
	MSPlot["MS_Pt_mu4_mu"] = new MultiSamplePlot(datasetsMu,"Pt_mu4", 50, 0, 300, "Pt of muon 4 (GeV)");
	MSPlot["MS_Pt_allEl_mu"] = new MultiSamplePlot(datasetsMu,"electron pt", 75, 0, 250, "Pt of all electrons (GeV)");
	MSPlot["MS_Pt_el1_mu"] = new MultiSamplePlot(datasetsMu,"Pt_el1", 50, 0, 300, "Pt of electron 1 (GeV)");
	MSPlot["MS_Pt_el2_mu"] = new MultiSamplePlot(datasetsMu,"Pt_el2", 50, 0, 300, "Pt of electron 2 (GeV)");
	MSPlot["MS_Pt_el3_mu"] = new MultiSamplePlot(datasetsMu,"Pt_el3", 50, 0, 300, "Pt of electron 3 (GeV)");
	//MSPlot["MS_Pt_el4_mu"] = new MultiSamplePlot(datasetsMu,"Pt_el4", 50, 0, 300, "Pt of electron 4 (GeV)");
	MSPlot["MS_Eta_allMu_mu"] = new MultiSamplePlot(datasetsMu,"muon eta", 50, -2.5 , 2.5, "Eta of all muons");
	MSPlot["MS_Eta_mu1_mu"] = new MultiSamplePlot(datasetsMu,"Eta_mu1", 50, -2.5 , 2.5, "Eta of muon 1");
	MSPlot["MS_Eta_mu2_mu"] = new MultiSamplePlot(datasetsMu,"Eta_mu2", 50, -2.5 , 2.5, "Eta of muon 2");
	MSPlot["MS_Eta_mu3_mu"] = new MultiSamplePlot(datasetsMu,"Eta_mu3", 50, -2.5 , 2.5, "Eta of muon 3");
	MSPlot["MS_Eta_mu4_mu"] = new MultiSamplePlot(datasetsMu,"Eta_mu4", 50, -2.5 , 2.5, "Eta of muon 4");
	MSPlot["MS_Eta_allEl_mu"] = new MultiSamplePlot(datasetsMu,"electron eta", 50, -2.5 , 2.5, "Eta of all electrons");
	MSPlot["MS_Eta_el1_mu"] = new MultiSamplePlot(datasetsMu,"Eta_el1", 50, -2.5 , 2.5, "Eta of electron 1");
	MSPlot["MS_Eta_el2_mu"] = new MultiSamplePlot(datasetsMu,"Eta_el2", 50, -2.5 , 2.5, "Eta of electron 2");
	MSPlot["MS_Eta_el3_mu"] = new MultiSamplePlot(datasetsMu,"Eta_el3", 50, -2.5 , 2.5, "Eta of electron 3");
	//MSPlot["MS_Eta_el4_mu"] = new MultiSamplePlot(datasetsMu,"Eta_el4", 50, -2.5 , 2.5, "Eta of muon 4");
	//'boxes'
	MSPlot["MS_St_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box Wqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_MQ_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Mass Q (box Wqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxWHqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box WHqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxZqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box Zqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_MQ_boxZqq_mu"] = new MultiSamplePlot(datasetsMu,"Mass Q (box Zqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxZHqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box ZHqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxWWqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box WWqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_nEvts_boxWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxZZqq_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box ZZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqqZZqq_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq+ZZqq)", 1, 0.5, 1.5, "");
	
	//'el channel' = when NO muon is trigged, and at least one electron
	MSPlot["MS_nPV_noPUreweighting_el"] = new MultiSamplePlot(datasetsEl, "nPrimaryVertices before PU reweighting", 31, -0.5, 30.5, "Nr. of primary vertices");
  MSPlot["MS_nPV_el"] = new MultiSamplePlot(datasetsEl, "nPrimaryVertices", 31, -0.5, 30.5, "Nr. of primary vertices");  //after PU reweighting
	MSPlot["MS_nPV_noPUreweighting_4jets_el"] = new MultiSamplePlot(datasetsEl, "nPrimaryVertices before PU reweighting", 31, -0.5, 30.5, "Nr. of primary vertices");
  MSPlot["MS_nPV_4jets_el"] = new MultiSamplePlot(datasetsEl, "nPrimaryVertices", 31, -0.5, 30.5, "Nr. of primary vertices");  //after PU reweighting
	MSPlot["MS_nLeptons_el"] = new MultiSamplePlot(datasetsEl,"number of leptons", 6, -0.5, 5.5, "number of leptons");
	MSPlot["MS_nMuons_el"] = new MultiSamplePlot(datasetsEl,"number of muons", 6, -0.5, 5.5, "number of muons");
	MSPlot["MS_nElectrons_el"] = new MultiSamplePlot(datasetsEl,"number of electrons", 6, -0.5, 5.5, "number of electrons");
	MSPlot["MS_MET_el"] = new MultiSamplePlot(datasetsEl,"MET", 200, 0, 1000, "Missing transverse energy (GeV)");  
	MSPlot["MS_MET_singleEl_el"] = new MultiSamplePlot(datasetsEl,"MET", 200, 0, 1000, "Missing transverse energy (GeV)");
	//jets (not forward)
	MSPlot["MS_JetMultiplicity_el"] = new MultiSamplePlot(datasetsEl, "JetMultiplicity", 10, -0.5, 9.5, "jet multiplicity");
	MSPlot["MS_Pt_alljets_el"] = new MultiSamplePlot(datasetsEl,"Pt_alljets", 50, 0, 300, "Pt of all jets (GeV)");
	MSPlot["MS_Pt_jet1_el"] = new MultiSamplePlot(datasetsEl,"Pt_jet1", 50, 0, 300, "Pt of jet 1 (GeV)");
	MSPlot["MS_Pt_jet2_el"] = new MultiSamplePlot(datasetsEl,"Pt_jet2", 50, 0, 300, "Pt of jet 2 (GeV)");
	MSPlot["MS_Pt_jet3_el"] = new MultiSamplePlot(datasetsEl,"Pt_jet3", 50, 0, 300, "Pt of jet 3 (GeV)");
	MSPlot["MS_Pt_jet4_el"] = new MultiSamplePlot(datasetsEl,"Pt_jet4", 50, 0, 300, "Pt of jet 4 (GeV)");
	MSPlot["MS_Eta_alljets_el"] = new MultiSamplePlot(datasetsEl,"Eta_alljets", 50, -2.5 , 2.5, "Eta of all jets");
	MSPlot["MS_Eta_jet1_el"] = new MultiSamplePlot(datasetsEl,"Eta_jet1", 50, -2.5 , 2.5, "Eta of jet 1");
	MSPlot["MS_Eta_jet2_el"] = new MultiSamplePlot(datasetsEl,"Eta_jet2", 50, -2.5 , 2.5, "Eta of jet 2");
	MSPlot["MS_Eta_jet3_el"] = new MultiSamplePlot(datasetsEl,"Eta_jet3", 50, -2.5 , 2.5, "Eta of jet 3");
	MSPlot["MS_Eta_jet4_el"] = new MultiSamplePlot(datasetsEl,"Eta_jet4", 50, -2.5 , 2.5, "Eta of jet 4");
	//forward jets
	MSPlot["MS_ForwardJetMultiplicity_el"] = new MultiSamplePlot(datasetsEl, "JetMultiplicity", 10, -0.5, 9.5, "Forward jet multiplicity");
	MSPlot["MS_Pt_allForwardjets_el"] = new MultiSamplePlot(datasetsEl,"Pt_alljets", 50, 0, 300, "Pt of all forward jets (GeV)");
	MSPlot["MS_Pt_Forwardjet1_el"] = new MultiSamplePlot(datasetsEl,"Pt_Forwardjet1", 50, 0, 300, "Pt of forward jet 1 (GeV)");
	MSPlot["MS_Pt_Forwardjet2_el"] = new MultiSamplePlot(datasetsEl,"Pt_Forwardjet2", 50, 0, 300, "Pt of forward jet 2 (GeV)");
	MSPlot["MS_Eta_allForwardjets_el"] = new MultiSamplePlot(datasetsEl,"Eta_allForwardjets", 100, -5 , 5, "Eta of all forward jets");
	MSPlot["MS_Eta_Forwardjet1_el"] = new MultiSamplePlot(datasetsEl,"Eta_Forwardjet1", 100, -5 , 5, "Eta of forward jet 1");
	MSPlot["MS_Eta_Forwardjet2_el"] = new MultiSamplePlot(datasetsEl,"Eta_Forwardjet2", 100, -5 , 5, "Eta of forward jet 2");
	//leptons
  MSPlot["MS_Pt_allMu_el"] = new MultiSamplePlot(datasetsEl,"muon pt", 75, 0, 250, "Pt of all muons (GeV)");
	MSPlot["MS_Pt_el1_el"] = new MultiSamplePlot(datasetsEl,"Pt_el1", 50, 0, 300, "Pt of muon 1 (GeV)");
	MSPlot["MS_Pt_el2_el"] = new MultiSamplePlot(datasetsEl,"Pt_el2", 50, 0, 300, "Pt of muon 2 (GeV)");
	MSPlot["MS_Pt_el3_el"] = new MultiSamplePlot(datasetsEl,"Pt_el3", 50, 0, 300, "Pt of muon 3 (GeV)");
	MSPlot["MS_Pt_el4_el"] = new MultiSamplePlot(datasetsEl,"Pt_el4", 50, 0, 300, "Pt of muon 4 (GeV)");
	MSPlot["MS_Pt_allEl_el"] = new MultiSamplePlot(datasetsEl,"electron pt", 75, 0, 250, "Pt of all electrons (GeV)");
	MSPlot["MS_Pt_el1_el"] = new MultiSamplePlot(datasetsEl,"Pt_el1", 50, 0, 300, "Pt of electron 1 (GeV)");
	MSPlot["MS_Pt_el2_el"] = new MultiSamplePlot(datasetsEl,"Pt_el2", 50, 0, 300, "Pt of electron 2 (GeV)");
	MSPlot["MS_Pt_el3_el"] = new MultiSamplePlot(datasetsEl,"Pt_el3", 50, 0, 300, "Pt of electron 3 (GeV)");
	//MSPlot["MS_Pt_el4_el"] = new MultiSamplePlot(datasetsEl,"Pt_el4", 50, 0, 300, "Pt of electron 4 (GeV)");
	MSPlot["MS_Eta_allMu_el"] = new MultiSamplePlot(datasetsEl,"muon eta", 50, -2.5 , 2.5, "Eta of all muons");
	MSPlot["MS_Eta_el1_el"] = new MultiSamplePlot(datasetsEl,"Eta_el1", 50, -2.5 , 2.5, "Eta of muon 1");
	MSPlot["MS_Eta_el2_el"] = new MultiSamplePlot(datasetsEl,"Eta_el2", 50, -2.5 , 2.5, "Eta of muon 2");
	MSPlot["MS_Eta_el3_el"] = new MultiSamplePlot(datasetsEl,"Eta_el3", 50, -2.5 , 2.5, "Eta of muon 3");
	MSPlot["MS_Eta_el4_el"] = new MultiSamplePlot(datasetsEl,"Eta_el4", 50, -2.5 , 2.5, "Eta of muon 4");
	MSPlot["MS_Eta_allEl_el"] = new MultiSamplePlot(datasetsEl,"electron eta", 50, -2.5 , 2.5, "Eta of all electrons");
	MSPlot["MS_Eta_el1_el"] = new MultiSamplePlot(datasetsEl,"Eta_el1", 50, -2.5 , 2.5, "Eta of electron 1");
	MSPlot["MS_Eta_el2_el"] = new MultiSamplePlot(datasetsEl,"Eta_el2", 50, -2.5 , 2.5, "Eta of electron 2");
	MSPlot["MS_Eta_el3_el"] = new MultiSamplePlot(datasetsEl,"Eta_el3", 50, -2.5 , 2.5, "Eta of electron 3");
	//MSPlot["MS_Eta_el4_el"] = new MultiSamplePlot(datasetsEl,"Eta_el4", 50, -2.5 , 2.5, "Eta of muon 4");
	//'boxes'
	MSPlot["MS_St_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box Wqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_MQ_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Mass Q (box Wqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxWHqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box WHqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxZqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box Zqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_MQ_boxZqq_el"] = new MultiSamplePlot(datasetsEl,"Mass Q (box Zqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxZHqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box ZHqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxWWqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box WWqq)", 75, 0, 1500, "S_{T} (GeV)");
	MSPlot["MS_nEvts_boxWZqq_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxZZqq_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box ZZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqqZZqq_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq+ZZqq)", 1, 0.5, 1.5, "");
	
	
	
	
	
	
	
	
	/*MSPlot["MS_JetPt_all_SingleLepton_nobtag"] = new MultiSamplePlot(datasets,"JetPt_all", 50, 0, 300, "Pt of all jets (GeV)");
	MSPlot["MS_JetPt_btagged_SingleLepton_nobtag"] = new MultiSamplePlot(datasets,"JetPt_btagged", 50, 0, 300, "Pt of b-tagged jets (GeV)");
	MSPlot["MS_JetPt_nonbtagged_SingleLepton_nobtag"] = new MultiSamplePlot(datasets,"JetPt_nonbtagged", 50, 0, 300, "Pt of non b-tagged jets (GeV)");
  MSPlot["MS_Reliso_Lepton_nobtag"] = new MultiSamplePlot(datasets, "Lepton reliso", 40, 0, 0.25, "Lepton reliso");

  //4 jets, no b-tag requirement yet
  MSPlot["MS_MET_nobtag_4jets"] = new MultiSamplePlot(datasets,"MET", 75, 0, 200, "Missing transverse energy (GeV)");
  MSPlot["MS_LeptonPt_nobtag_4jets"] = new MultiSamplePlot(datasets,"lepton pt", 75, 0, 250, "Pt lepton (GeV)");
  MSPlot["MS_JetPt_all_SingleLepton_nobtag_4jets"] = new MultiSamplePlot(datasets,"JetPt_all", 50, 0, 300, "Pt of all jets (GeV)");
	//4 jets, with b-tag requirement
  MSPlot["MS_MET_4jets"] = new MultiSamplePlot(datasets,"MET", 75, 0, 200, "Missing transverse energy (GeV)");
  MSPlot["MS_LeptonPt_4jets"] = new MultiSamplePlot(datasets,"lepton pt", 75, 0, 250, "Pt lepton (GeV)");
  MSPlot["MS_JetPt_all_SingleLepton_4jets"] = new MultiSamplePlot(datasets,"JetPt_all", 50, 0, 300, "Pt of all jets (GeV)");	

  //for the events that are actually used
	MSPlot["MS_MET_used"] = new MultiSamplePlot(datasets,"MET", 75, 0, 200, "Missing transverse energy (GeV)");
  MSPlot["MS_LeptonPt_used"] = new MultiSamplePlot(datasets,"lepton pt", 75, 0, 250, "Pt lepton (GeV)");
  MSPlot["MS_JetPt_all_used"] = new MultiSamplePlot(datasets,"JetPt_all", 50, 0, 300, "Pt of all jets (GeV)");	
	MSPlot["MS_JetHt_all_used"] = new MultiSamplePlot(datasets,"JetHt_all", 80, 0, 800, "Ht of all jets (GeV)");
	MSPlot["MS_nPV_used"] = new MultiSamplePlot(datasets, "nPrimaryVertices", 21, -0.5, 20.5, "Nr. of primary vertices");
	
	histo2D["HTvsMTop_1B_2W_TTbarJets"] = new TH2F("HTvsMTop_1B_2W_TTbarJets","HTvsMTop_1B_2W_TTbarJets",400,0,1500,400,0,1500);
	histo2D["HTvsMTop_2B_2W_TTbarJets"] = new TH2F("HTvsMTop_2B_2W_TTbarJets","HTvsMTop_2B_2W_TTbarJets",400,0,1500,400,0,1500);
	histo2D["HTvsMTop_1B_2W_TprimeTprime"] = new TH2F("HTvsMTop_1B_2W_TprimeTprime","HTvsMTop_1B_2W_TprimeTprime",400,0,1500,400,0,1500);
	histo2D["HTvsMTop_2B_2W_TprimeTprime"] = new TH2F("HTvsMTop_2B_2W_TprimeTprime","HTvsMTop_2B_2W_TprimeTprime",400,0,1500,400,0,1500);*/
	
	cout << " - Declared histograms ..." <<  endl;

  //Configuration and variables for 2D distribution.

	
  ////////////////////////////////////
  /// Selection table   --> to be modified
  ////////////////////////////////////

  vector<string> CutsSelecTableLep; 
  CutsSelecTableLep.push_back(string("baseline selected")); //0
  CutsSelecTableLep.push_back(string("leading jet Pt > 200 GeV")); // semi-leptonic events in all the boxes together
	
	vector<string> CutsSelecTableMu; 
  CutsSelecTableMu.push_back(string("baseline selected (mu trigged)")); //0
  CutsSelecTableMu.push_back(string("leading jet Pt > 200 GeV")); // semi-leptonic events in all the boxes together
	
	
	

  vector<string> CutsSelecTableMultiLep;
  CutsSelecTableMultiLep.push_back(string("SS leptons: all")); //0
  CutsSelecTableMultiLep.push_back(string("SS leptons: 2 muons")); //1
  CutsSelecTableMultiLep.push_back(string("SS leptons: 2 electrons")); //2
  CutsSelecTableMultiLep.push_back(string("SS leptons: electron+muon")); //3
  CutsSelecTableMultiLep.push_back(string("trileptons in all boxes combined")); //4
  
	SelectionTable selecTableLep(CutsSelecTableLep, datasets);
  selecTableLep.SetLuminosity(Luminosity);
  selecTableLep.SetPrecision(1);
	
	SelectionTable selecTableMu(CutsSelecTableMu, datasets);
  selecTableMu.SetLuminosity(Luminosity);
  selecTableMu.SetPrecision(1);
	

	SelectionTable selecTableMultiLep(CutsSelecTableMultiLep, datasets);
  selecTableMultiLep.SetLuminosity(Luminosity);
  selecTableMultiLep.SetPrecision(2);
  
  cout << " - SelectionTable instantiated ..." << endl;

  
	////////////////////////
  // PileUp Reweighting //
  ////////////////////////

  //cout << Luminosity << endl;

  LumiReWeighting LumiWeights, LumiWeightsUp, LumiWeightsDown;

  //For Run2012A+B only! (temporary)
  LumiWeights = LumiReWeighting("../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../PileupHistos/Run2012AB_SingleMu/pileup_2012Data53X_UpToRun196531_RunAB/nominal.root", "pileup", "pileup");
  LumiWeightsUp = LumiReWeighting("../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../PileupHistos/Run2012AB_SingleMu/pileup_2012Data53X_UpToRun196531_RunAB/sys_up.root", "pileup", "pileup");
  LumiWeightsDown = LumiReWeighting("../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../PileupHistos/Run2012AB_SingleMu/pileup_2012Data53X_UpToRun196531_RunAB/sys_down.root", "pileup", "pileup");

  cout << " Initialized LumiReWeighting stuff" << endl;
	  
  // initialize lepton SF
	//WARNING: should be checked if correct if only run on RunA+B (temporary)
  LeptonTools* leptonTools = new LeptonTools(false);
  leptonTools->readMuonSF("../../TopTreeAnalysisBase/Calibrations/LeptonSF/Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root", "../../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonEfficiencies_Run_2012A_2012B_53X.root", "../../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonEfficiencies_Run_2012C_53X.root", "../../TopTreeAnalysisBase/Calibrations/LeptonSF/TriggerMuonEfficiencies_Run_2012D_53X.root");
  leptonTools->readElectronSF();
    
		
	////////////////////////
  // HCAL laser filter; because not yet applied on PAT level //
  ////////////////////////
  string fname = "AllBadHCALLaser.txt";
  //string fnameS = "AllBadHCALLaser"+data_postfix+".txt";	
	string fnameS = "AllBadHCALLaser_smaller.txt";	
	ifstream testF(fnameS.c_str());
	if (!testF) 
	{
        cout << "BadHCALLaserEvents:: using large file "+fname+" to remove bad events, creating subset "+fnameS << endl;
        testF.close();            
  } 
  else 
	{     
        fname = fnameS;        
        cout << "BadHCALLaserEvents:: using smaller file "+fname+" to remove bad events" << endl;        
  }
	
	ifstream badev(fname.c_str());
  vector<string> bad_event;
  while (!badev.eof()) {        
        string ev_;        
        badev >> ev_;        
        bad_event.push_back(ev_);        
        //cout << ev_ << endl;        
  }
  badev.close();	
		
		
  ////////////////////////////////////
  ////////////////////////////////////
  ///////// Loop on datasets
  ////////////////////////////////////
  ////////////////////////////////////
  cout << " - Loop over datasets ... " << inputTrees.size() << " datasets !" << endl;
  ofstream myfile1, myfile2;
	string myRockingFile1 = Outputpath+"InterestingEvents"+".txt";
	if(systematic=="Nominal") myfile1.open(myRockingFile1.c_str());

  
	for (unsigned int d = 0; d < inputTrees.size(); d++)
  {
		string dataSetName = datasets[d]->Name();
			
		//if(make2Dbinning && dataSetName != "TTbarJets_SemiMuon" && dataSetName != "TTbarJets_SemiElectron" && dataSetName != "TTbarJets_Other")
		//	continue;
			
		if(systematic != "Nominal" && !(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0))
		{
		  cout<<"  Running systematics: skipping data"<<endl;
		  continue;
		}	
			
    if (verbose > 1)
      cout << "file: " << inputTrees[d] << endl;
	 
    TFile* inFile = new TFile(inputTrees[d].c_str(),"READ");
   // TTree* inConfigTree = (TTree*) inFile->Get("configTreeFile");
   // TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
   // TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
   // d_br->SetAddress(&tc_dataset);
   // inConfigTree->GetEvent(0);
   // Dataset* dataSet = (Dataset*) tc_dataset->At(0);	//will not be used (problem with MVA 1/3 and 2/3 splitting of sample)	
		
		TTree* inVLQTree = (TTree*) inFile->Get("myVLQTree");
    TBranch* m_br = (TBranch*) inVLQTree->GetBranch("VLQBranch_selectedEvents");
    VLQTree* myBranch_selectedEvents = 0;
    m_br->SetAddress(&myBranch_selectedEvents);
    int nEvents = inVLQTree->GetEntries();			
		
		cout << " Processing DataSet: " << dataSetName << "  containing " << nEvents << " events" << endl;
    cout << " Cross section = " << datasets[d]->Xsection() << "  intLumi = " << datasets[d]->EquivalentLumi() << "  NormFactor = " << datasets[d]->NormFactor() << endl;

		
    ///////////////////////////////////////////////////
    // calculate W reconstructed mass and resolution: load the mass values and resolutions for chi2 jetcombination.
		// Can this be put outside dataset loop?
    ///////////////////////////////////////////////////
    
    /*float Wmass = -9999;
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
    }*/


    ////////////////////////////////////
    ////////////////////////////////////
    ///////// Loop on events
    ////////////////////////////////////
    ////////////////////////////////////

    int start = 0;
    int end = nEvents;
     
    cout << " - Loop over events " << endl;      
    //cout<<"start = "<<start<<"end = "<<end<<endl;
		
    for (int ievt = start; ievt < end; ievt++)
    {        

      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";
      systematicnew = "";
			//load event
			inVLQTree->GetEvent(ievt);
	 
	 
	    //reading variables from the tree
			unsigned int nPV = myBranch_selectedEvents->nPV();
			unsigned int nTruePU = myBranch_selectedEvents->nTruePU();
			
			bool MuTrigged = myBranch_selectedEvents->MuTrigged();
			bool ElTrigged = myBranch_selectedEvents->ElTrigged();
      unsigned int SelectednLeptons = myBranch_selectedEvents->SelectednLeptons();
      unsigned int SelectednMu = myBranch_selectedEvents->SelectednMu();
      unsigned int SelectednEl = myBranch_selectedEvents->SelectednEl();			
						
			float met = (myBranch_selectedEvents->met()).Et();
      vector<TLorentzVector> selectedJets = myBranch_selectedEvents->selectedJets();
			vector<TLorentzVector> selectedForwardJets = myBranch_selectedEvents->selectedForwardJets();
			vector<TLorentzVector> selectedMuons = myBranch_selectedEvents->selectedMuons();
			vector<TLorentzVector> selectedElectrons = myBranch_selectedEvents->selectedElectrons();
			vector<float> selectedJets_bTagCSV = myBranch_selectedEvents->selectedJets_bTagCSV();
			vector<float> selectedForwardJets_bTagCSV = myBranch_selectedEvents->selectedForwardJets_bTagCSV();
			
	 
	 
      // scale factor for the event
      float scaleFactor = 1.;

      double lumiWeight = LumiWeights.ITweight( (int) nTruePU );

      if((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos))
			{
	       lumiWeight=1.;
				 
				 stringstream r; r<<myBranch_selectedEvents->runID();
         stringstream t; t<<myBranch_selectedEvents->lumiBlockID();
         stringstream u; u<<myBranch_selectedEvents->eventID();
				 string ev = r.str()+":"+t.str()+":"+u.str();
				 int pos = std::find(bad_event.begin(),bad_event.end(),ev)-bad_event.begin();                
         if (pos < bad_event.size()) 
				 {
               cout << "BadHCALLaserEvents:: Skipping " << ev << " " << pos << endl;
               if (fnameS != fname) 
							 {
                        ofstream fout;
                        fout.open(fnameS.c_str(),ios::out | ios::app);
                        fout << ev << "\n";
                        fout.close();
                }
                    
               continue;
         }
			}
      
      // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int) myBranch_selectedEvents->nTruePU() );
      // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int) myBranch_selectedEvents->nTruePU() );

      scaleFactor = scaleFactor*lumiWeight;
      histo1D["lumiWeights"]->Fill(scaleFactor);
			   

			//jet Pt cut
			/*vector<TLorentzVector> selectedJets_cuts, selectedForwardJets_cuts;
			float JetPtCut = 30;
			for(unsigned int i=0;i<selectedJets.size();i++)
			{
				if(selectedJets[i].Pt() > JetPtCut)
					selectedJets_cuts.push_back(selectedJets[i]);													
			}
			selectedJets.clear();
			selectedJets = selectedJets_cuts;
			for(unsigned int i=0;i<selectedForwardJets.size();i++)
			{
				if(selectedForwardJets[i].Pt() > JetPtCut)
					selectedForwardJets_cuts.push_back(selectedForwardJets[i]);													
			}
			selectedForwardJets.clear();
			selectedForwardJets = selectedForwardJets_cuts;			
			*/
			
			/*// when running both electron and muon data, pick the right dataset vector and lumi for the MSPlots
     
      if (!foundMu && !foundEl)
        datasetsPlot = datasets;
      else if (eventselectedSemiMu) {
        datasetsPlot = datasetsMu;
        Luminosity = LuminosityMu;
      }
      else if (eventselectedSemiEl) {
        datasetsPlot = datasetsEl;
        Luminosity = LuminosityEl;
      } */

      //'muon channel'
      //if(SelectednMu >= 1) //not anymore
			if(MuTrigged)
			{
			
			  //lepton scale factor; apply on MC
				if(!((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			  {
				   scaleFactor = scaleFactor*leptonTools->getMuonSF(selectedMuons[0].Eta(), selectedMuons[0].Pt());
				}
				 
			  //MSPlot["MS_nPV_noPUreweighting_mu"]->Fill(nPV,datasets[d], true, LuminosityMu*scaleFactor / lumiWeight);	
			  MSPlot["MS_nPV_noPUreweighting_mu"]->Fill(nPV,datasets[d], true, LuminosityMu);	
			  MSPlot["MS_nPV_mu"]->Fill(nPV,datasets[d], true, LuminosityMu*scaleFactor);
				MSPlot["MS_nLeptons_mu"]->Fill(SelectednLeptons,datasets[d], true, LuminosityMu*scaleFactor);
	      MSPlot["MS_nMuons_mu"]->Fill(SelectednMu,datasets[d], true, LuminosityMu*scaleFactor);
	      MSPlot["MS_nElectrons_mu"]->Fill(SelectednEl,datasets[d], true, LuminosityMu*scaleFactor);
	      MSPlot["MS_MET_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
				
				if(SelectednMu == 1 && SelectednEl == 0) 
				{
				   if(selectedJets[0].Pt() >= LeadingJetPtCut*0.5)
					    MSPlot["MS_MET_singleMu_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
				}
				
				if(selectedJets.size() == 4)
			  {
			    MSPlot["MS_nPV_noPUreweighting_4jets_mu"]->Fill(nPV,datasets[d], true, LuminosityMu);	
			    MSPlot["MS_nPV_4jets_mu"]->Fill(nPV,datasets[d], true, LuminosityMu*scaleFactor);
			  }
				
				MSPlot["MS_JetMultiplicity_mu"]->Fill(selectedJets.size(),datasets[d], true, LuminosityMu*scaleFactor);
	      for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
			  {
				   MSPlot["MS_Pt_alljets_mu"]->Fill(selectedJets[iJet].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
				}
	      if( selectedJets.size() >= 1 ) MSPlot["MS_Pt_jet1_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);				
	      if( selectedJets.size() >= 2 ) MSPlot["MS_Pt_jet2_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
	      if( selectedJets.size() >= 3 ) MSPlot["MS_Pt_jet3_mu"]->Fill(selectedJets[2].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
	      if( selectedJets.size() >= 4 ) MSPlot["MS_Pt_jet4_mu"]->Fill(selectedJets[3].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
	      for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
			  {
				   MSPlot["MS_Eta_alljets_mu"]->Fill(selectedJets[iJet].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
				}
	      if( selectedJets.size() >= 1 ) MSPlot["MS_Eta_jet1_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);				
	      if( selectedJets.size() >= 2 ) MSPlot["MS_Eta_jet2_mu"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
	      if( selectedJets.size() >= 3 ) MSPlot["MS_Eta_jet3_mu"]->Fill(selectedJets[2].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
	      if( selectedJets.size() >= 4 ) MSPlot["MS_Eta_jet4_mu"]->Fill(selectedJets[3].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
	      MSPlot["MS_ForwardJetMultiplicity_mu"]->Fill(selectedForwardJets.size(),datasets[d], true, LuminosityMu*scaleFactor);
				for(unsigned int iJet = 0; iJet < selectedForwardJets.size(); iJet++)
			  {
				   MSPlot["MS_Pt_allForwardjets_mu"]->Fill(selectedForwardJets[iJet].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
				}
	      if( selectedForwardJets.size() >= 1 ) MSPlot["MS_Pt_Forwardjet1_mu"]->Fill(selectedForwardJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
	      if( selectedForwardJets.size() >= 2 ) MSPlot["MS_Pt_Forwardjet2_mu"]->Fill(selectedForwardJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
	      for(unsigned int iJet = 0; iJet < selectedForwardJets.size(); iJet++)
			  {
				   MSPlot["MS_Eta_allForwardjets_mu"]->Fill(selectedForwardJets[iJet].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
				}
	      if( selectedForwardJets.size() >= 1 ) MSPlot["MS_Eta_Forwardjet1_mu"]->Fill(selectedForwardJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
	      if( selectedForwardJets.size() >= 2 ) MSPlot["MS_Eta_Forwardjet2_mu"]->Fill(selectedForwardJets[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
	      //leptons
				for(unsigned int iMu = 0; iMu < selectedMuons.size(); iMu++)
			  {
				   MSPlot["MS_Pt_allMu_mu"]->Fill(selectedMuons[iMu].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
					 MSPlot["MS_Eta_allMu_mu"]->Fill(selectedMuons[iMu].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
			  }
	      if( selectedMuons.size() >= 1 )
				{   
				    MSPlot["MS_Pt_mu1_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
						MSPlot["MS_Eta_mu1_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
				}
				if( selectedMuons.size() >= 2 )
				{   
				    MSPlot["MS_Pt_mu2_mu"]->Fill(selectedMuons[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
						MSPlot["MS_Eta_mu2_mu"]->Fill(selectedMuons[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
				}
				if( selectedMuons.size() >= 3 )
				{   
				    MSPlot["MS_Pt_mu3_mu"]->Fill(selectedMuons[2].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
						MSPlot["MS_Eta_mu3_mu"]->Fill(selectedMuons[2].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
				}
				if( selectedMuons.size() >= 4 )
				{   
				    MSPlot["MS_Pt_mu3_mu"]->Fill(selectedMuons[2].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
						MSPlot["MS_Eta_mu3_mu"]->Fill(selectedMuons[2].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
				}
				
				for(unsigned int iEl = 0; iEl < selectedElectrons.size(); iEl++)
			  {
				   MSPlot["MS_Pt_allEl_mu"]->Fill(selectedElectrons[iEl].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
					 MSPlot["MS_Eta_allEl_mu"]->Fill(selectedElectrons[iEl].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
			  }
	      if( selectedElectrons.size() >= 1 )
				{   
				    MSPlot["MS_Pt_el1_mu"]->Fill(selectedElectrons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
						MSPlot["MS_Eta_el1_mu"]->Fill(selectedElectrons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
				}
				if( selectedElectrons.size() >= 2 )
				{   
				    MSPlot["MS_Pt_el2_mu"]->Fill(selectedElectrons[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
						MSPlot["MS_Eta_el2_mu"]->Fill(selectedElectrons[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
				}
				if( selectedElectrons.size() >= 3 )
				{   
				    MSPlot["MS_Pt_el3_mu"]->Fill(selectedElectrons[2].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
						MSPlot["MS_Eta_el3_mu"]->Fill(selectedElectrons[2].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
				}
				//if( selectedElectrons.size() >= 4 )
				//{   
				//    MSPlot["MS_Pt_el3_mu"]->Fill(selectedElectrons[2].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
				//		MSPlot["MS_Eta_el3_mu"]->Fill(selectedElectrons[2].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
				//}
 
				selecTableMu.Fill(d,0,scaleFactor);
				if(selectedJets[0].Pt() >= LeadingJetPtCut) selecTableMu.Fill(d,1,scaleFactor);
			}
			
			//'electron channel'
			//if(SelectednMu == 0 && SelectednEl >= 1) //not anymore
			else if(ElTrigged)  //or an 'if'??
			{
			  //lepton scale factor; apply on MC
				if(!((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			  {
				   scaleFactor = scaleFactor*leptonTools->getElectronSF(selectedElectrons[0].Eta(), selectedElectrons[0].Pt());
				}
				
			  //cout<<"d = "<<d<<endl;
			  MSPlot["MS_nPV_noPUreweighting_el"]->Fill(nPV,datasets[d], true, LuminosityEl);	
			  MSPlot["MS_nPV_el"]->Fill(nPV,datasets[d], true, LuminosityEl*scaleFactor);
				MSPlot["MS_nLeptons_el"]->Fill(SelectednLeptons,datasets[d], true, LuminosityEl*scaleFactor);
	      MSPlot["MS_nMuons_el"]->Fill(SelectednMu,datasets[d], true, LuminosityEl*scaleFactor);
	      MSPlot["MS_nElectrons_el"]->Fill(SelectednEl,datasets[d], true, LuminosityEl*scaleFactor);
	      MSPlot["MS_MET_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
				
				if(SelectednEl == 1 && SelectednMu == 0)
				{
				   if(selectedJets[0].Pt() >= LeadingJetPtCut*0.5)
				        MSPlot["MS_MET_singleEl_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
				}
				
				if(selectedJets.size() == 4)
			  {
			    MSPlot["MS_nPV_noPUreweighting_4jets_el"]->Fill(nPV,datasets[d], true, LuminosityEl);	
			    MSPlot["MS_nPV_4jets_el"]->Fill(nPV,datasets[d], true, LuminosityEl*scaleFactor);
			  }
				MSPlot["MS_JetMultiplicity_el"]->Fill(selectedJets.size(),datasets[d], true, LuminosityEl*scaleFactor);
	      for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
			  {
				   MSPlot["MS_Pt_alljets_el"]->Fill(selectedJets[iJet].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
				}
	      if( selectedJets.size() >= 1 ) MSPlot["MS_Pt_jet1_el"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);				
	      if( selectedJets.size() >= 2 ) MSPlot["MS_Pt_jet2_el"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
	      if( selectedJets.size() >= 3 ) MSPlot["MS_Pt_jet3_el"]->Fill(selectedJets[2].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
	      if( selectedJets.size() >= 4 ) MSPlot["MS_Pt_jet4_el"]->Fill(selectedJets[3].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
	      for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
			  {
				   MSPlot["MS_Eta_alljets_el"]->Fill(selectedJets[iJet].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
				}
	      if( selectedJets.size() >= 1 ) MSPlot["MS_Eta_jet1_el"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);				
	      if( selectedJets.size() >= 2 ) MSPlot["MS_Eta_jet2_el"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
	      if( selectedJets.size() >= 3 ) MSPlot["MS_Eta_jet3_el"]->Fill(selectedJets[2].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
	      if( selectedJets.size() >= 4 ) MSPlot["MS_Eta_jet4_el"]->Fill(selectedJets[3].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
	      MSPlot["MS_ForwardJetMultiplicity_el"]->Fill(selectedForwardJets.size(),datasets[d], true, LuminosityEl*scaleFactor);
				for(unsigned int iJet = 0; iJet < selectedForwardJets.size(); iJet++)
			  {
				   MSPlot["MS_Pt_allForwardjets_el"]->Fill(selectedForwardJets[iJet].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
				}
	      if( selectedForwardJets.size() >= 1 ) MSPlot["MS_Pt_Forwardjet1_el"]->Fill(selectedForwardJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
	      if( selectedForwardJets.size() >= 2 ) MSPlot["MS_Pt_Forwardjet2_el"]->Fill(selectedForwardJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
	      for(unsigned int iJet = 0; iJet < selectedForwardJets.size(); iJet++)
			  {
				   MSPlot["MS_Eta_allForwardjets_el"]->Fill(selectedForwardJets[iJet].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
				}
	      if( selectedForwardJets.size() >= 1 ) MSPlot["MS_Eta_Forwardjet1_el"]->Fill(selectedForwardJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
	      if( selectedForwardJets.size() >= 2 ) MSPlot["MS_Eta_Forwardjet2_el"]->Fill(selectedForwardJets[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
	      //leptons
				for(unsigned int iMu = 0; iMu < selectedMuons.size(); iMu++)
			  {
				   MSPlot["MS_Pt_allMu_el"]->Fill(selectedMuons[iMu].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
					 MSPlot["MS_Eta_allMu_el"]->Fill(selectedMuons[iMu].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
			  }
	      if( selectedMuons.size() >= 1 )
				{   
				    MSPlot["MS_Pt_el1_el"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
						MSPlot["MS_Eta_el1_el"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
				}
				if( selectedMuons.size() >= 2 )
				{   
				    MSPlot["MS_Pt_el2_el"]->Fill(selectedMuons[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
						MSPlot["MS_Eta_el2_el"]->Fill(selectedMuons[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
				}
				if( selectedMuons.size() >= 3 )
				{   
				    MSPlot["MS_Pt_el3_el"]->Fill(selectedMuons[2].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
						MSPlot["MS_Eta_el3_el"]->Fill(selectedMuons[2].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
				}
				if( selectedMuons.size() >= 4 )
				{   
				    MSPlot["MS_Pt_el3_el"]->Fill(selectedMuons[2].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
						MSPlot["MS_Eta_el3_el"]->Fill(selectedMuons[2].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
				}
				
				for(unsigned int iEl = 0; iEl < selectedElectrons.size(); iEl++)
			  {
				   MSPlot["MS_Pt_allEl_el"]->Fill(selectedElectrons[iEl].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
					 MSPlot["MS_Eta_allEl_el"]->Fill(selectedElectrons[iEl].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
			  }
	      if( selectedElectrons.size() >= 1 )
				{   
				    MSPlot["MS_Pt_el1_el"]->Fill(selectedElectrons[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
						MSPlot["MS_Eta_el1_el"]->Fill(selectedElectrons[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
				}
				if( selectedElectrons.size() >= 2 )
				{   
				    MSPlot["MS_Pt_el2_el"]->Fill(selectedElectrons[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
						MSPlot["MS_Eta_el2_el"]->Fill(selectedElectrons[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
				}
				if( selectedElectrons.size() >= 3 )
				{   
				    MSPlot["MS_Pt_el3_el"]->Fill(selectedElectrons[2].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
						MSPlot["MS_Eta_el3_el"]->Fill(selectedElectrons[2].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
				}
				//if( selectedElectrons.size() >= 4 )
				//{   
				//    MSPlot["MS_Pt_el3_el"]->Fill(selectedElectrons[2].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
				//		MSPlot["MS_Eta_el3_el"]->Fill(selectedElectrons[2].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
				//}
				
					
				
				
			}
			
			
			float massleptons = 9999.;
			float Zmass = 91.1876; //could eventually be replaced by the mean of a fit of the Z peak?
			float Zmasswindow = 15; //could eventually be replaced by the resolution of a fit of the Z peak (1 sigma)?
					
			float ST = met; //test search variable			
			for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
			{
						    ST = ST + selectedJets[iJet].Pt();
			}
			for(unsigned int iMu = 0; iMu < selectedMuons.size(); iMu++)
			{
						    ST = ST + selectedMuons[iMu].Pt();
			}
			for(unsigned int iEl = 0; iEl < selectedElectrons.size(); iEl++)
			{
						    ST = ST + selectedElectrons[iEl].Pt();
			}
					
			//define the 'boxes'
			//'muon channel'
			if(MuTrigged)
			{
			 if(SelectednMu == 1)
			 {
			  if(SelectednEl == 0)
				{
			    if( (selectedJets.size() == 1 || selectedJets.size() == 2) && selectedForwardJets.size() == 1 )
			    {
					  if(selectedJets[0].Pt() >= LeadingJetPtCut)
						{
						  if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0 && nBjetsPresent(selectedForwardJets_bTagCSV,antibtagWP)==0))
							{
					  	   //Wqq category						
							   ST = ST + selectedForwardJets[0].Pt();
							   MSPlot["MS_St_boxWqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
							   //MQ = 
							   //MSPlot["MS_MQ_boxWqq_mu"]->Fill(MQ,datasets[d], true, LuminosityMu*scaleFactor);
							}
						}				
				  }
				  else if( selectedJets.size() >= 3 )
			    {
					  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
						{
						  //b jet(s) from Higgs decay... requiring at least 1; not yet a mass constraint...
						  if(!doBtagging || (doBtagging && nBjetsPresent(selectedJets_bTagCSV,btagWP)>=1))
							{		
							   //WARNING: maybe better to make sure the high-Pt jets are not b-tagged (however, those from a Higgs decay are also high in Pt...)
							   if(!doAntiBtagging || (doAntiBtagging && nNonBjetsPresent(selectedJets_bTagCSV,antibtagWP)==2))
							   {					  
				    	      //WHqq category
							      MSPlot["MS_St_boxWHqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
								 }
							}	
						}		
				  }
					else
					{
						    //cout<<"do nothing...!?"<<endl;
					}
				}
				else if(SelectednEl == 1)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
					   if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
						 {
			          //WWqq category (no need to check the invariant mass of the leptons, they should not come from a Z...)
					      MSPlot["MS_St_boxWWqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);	
						 }
					}
				}
				else if(SelectednEl == 2)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
					  if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
						{
			    	   //WZqq category
						   MSPlot["MS_nEvts_boxWZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						   MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						}
					}
				}
				/*else if(SelectednEl == 3)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
					   if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
						 {
			          //ZZqq category
					      MSPlot["MS_nEvts_boxZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						    MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						 }
					}
				}*/
				else
				{
				   cout<<"FOUND >= 5-lepton event!"<<endl;		
				}
									
			 }
			 else if(SelectednMu == 2)
			 {
			  if(SelectednEl == 0)
				{
			    massleptons = (selectedMuons[0] + selectedMuons[1]).M();
					
					if( fabs(Zmass - massleptons) < Zmasswindow)
					{					
				     if( (selectedJets.size() == 1 || selectedJets.size() == 2) && selectedForwardJets.size() == 1 )
			       {
						   if(selectedJets[0].Pt() >= LeadingJetPtCut)
				       {
							    if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0 && nBjetsPresent(selectedForwardJets_bTagCSV,antibtagWP)==0))
						      {
				             //Zqq category
									   ST = ST + selectedForwardJets[0].Pt();
							       MSPlot["MS_St_boxZqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
									   //calculate mass of (two-lepton and leading jet)-system
									   float massZq = (selectedMuons[0] + selectedMuons[1] + selectedJets[0]).M();
									   MSPlot["MS_MQ_boxZqq_mu"]->Fill(massZq,datasets[d], true, LuminosityMu*scaleFactor);
								  }
							 }
				     }
				     else if( selectedJets.size() >= 3 )
			       {
						   if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				       {
							    if(!doBtagging || (doBtagging && nBjetsPresent(selectedJets_bTagCSV,btagWP)>=1))
							    {
									    //WARNING: maybe better to make sure the high-Pt jets are not b-tagged (however, those from a Higgs decay are also high in Pt...)
							        if(!doAntiBtagging || (doAntiBtagging && nNonBjetsPresent(selectedJets_bTagCSV,antibtagWP)==2))
							        {
				                 //ZHqq category
							           MSPlot["MS_St_boxZHqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
											}
									}
							 }
				     }
						 else
					   {
						    //cout<<"do nothing...!?"<<endl;
					   }
					}
					else
					{
					   if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				     {
						    if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
						    {
					         //WWqq category
						       MSPlot["MS_St_boxWWqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);			
								}
						 }		
					}
				}
				else if(SelectednEl == 1)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
					   if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
						 {
			          //WZqq category
					      MSPlot["MS_nEvts_boxWZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						    MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						 }
					}					
				}
				else if(SelectednEl == 2)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
					   if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
						 {
			          //ZZqq category
					      MSPlot["MS_nEvts_boxZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						    MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						 }
					}
				}
				else 
				{
				   cout<<"FOUND >= 5-lepton event!"<<endl;		
				}
			 }
			 else if(SelectednMu == 3)
			 {
			  if(SelectednEl == 0)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
					   if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
						 {
			          //WZqq category
					      MSPlot["MS_nEvts_boxWZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						    MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						 }
					}
				}
				/*else if(SelectednEl == 1)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
					   if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
						 {
			          //ZZqq category
					      MSPlot["MS_nEvts_boxZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						    MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						 }
					}
				}*/
				else 
				{
				   cout<<"FOUND >= 5-lepton event!"<<endl;		
				}
			 }
			 else if(SelectednMu == 4)
			 {
			  if(SelectednEl == 0)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
			       //ZZqq category
					   MSPlot["MS_nEvts_boxZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						 MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
					}
				}
				else 
				{
				   cout<<"FOUND >= 5-lepton event!"<<endl;		
				}
			 }
			}			
			
			//'electron channel'
			else if(ElTrigged) //or not 'else if' but 'if'???? (then it should be consistent with the TreeCreator...)
			{
			 if(SelectednMu == 0)
			 {
			  if(SelectednEl == 1)
				{	
				  if( (selectedJets.size() == 1 || selectedJets.size() == 2) && selectedForwardJets.size() == 1 )
			    { 
					   if(selectedJets[0].Pt() >= LeadingJetPtCut)
				     {
				        //Wqq category
						    ST = ST + selectedForwardJets[0].Pt();
						    MSPlot["MS_St_boxWqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
						 }
				  }
				  else if( selectedJets.size() >= 3 )
			    {
					   if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				     {
						    //b jet(s) from Higgs decay... requiring at least 1; not yet a mass constraint...
						    if(!doBtagging || (doBtagging && nBjetsPresent(selectedJets_bTagCSV,btagWP)>=1))
							  {		
							     //WARNING: maybe better to make sure the high-Pt jets are not b-tagged (however, those from a Higgs decay are also high in Pt...)
							     if(!doAntiBtagging || (doAntiBtagging && nNonBjetsPresent(selectedJets_bTagCSV,antibtagWP)==2))
							     {	
							        //WHqq category
						          MSPlot["MS_St_boxWHqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
							     }
								}
						 }
				  }
					else
					{
						    //cout<<"do nothing...!?"<<endl;
					}
				}
				else if(SelectednEl == 2)
				{									
					massleptons = (selectedElectrons[0] + selectedElectrons[1]).M();
					
					if( fabs(Zmass - massleptons) < Zmasswindow)
					{					
				     if( (selectedJets.size() == 1 || selectedJets.size() == 2) && selectedForwardJets.size() == 1 )
			       {
						 
						    if(selectedJets[0].Pt() >= LeadingJetPtCut)
				        {
				          //Zqq category
							    ST = ST + selectedForwardJets[0].Pt();
							    MSPlot["MS_St_boxZqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
									//calculate mass of (two-lepton and leading jet)-system
									float massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedJets[0]).M();
									MSPlot["MS_MQ_boxZqq_mu"]->Fill(massZq,datasets[d], true, LuminosityMu*scaleFactor);
								}
				     }
				     else if( selectedJets.size() >= 3 )
			       {
						    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				        {								
								   //b jet(s) from Higgs decay... requiring at least 1; not yet a mass constraint...
						       if(!doBtagging || (doBtagging && nBjetsPresent(selectedJets_bTagCSV,btagWP)>=1))
							     {	
									    //WARNING: maybe better to make sure the high-Pt jets are not b-tagged (however, those from a Higgs decay are also high in Pt...)
							        if(!doAntiBtagging || (doAntiBtagging && nNonBjetsPresent(selectedJets_bTagCSV,antibtagWP)==2))
							        {	
				                 //ZHqq category
							           //calculate mass of (two-lepton and leading jet)-system (CHECK THE JETS ARE SORTED!!!)
							           //float massZq1 = (selectedElectrons[0] + selectedElectrons[1] + selectedJets[0]).M();
							 
							           MSPlot["MS_St_boxZHqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
											}
									 }
								}
					 
				     }
						 else
						 {
						    //cout<<"do nothing...!?"<<endl;
						 }
					}
					else
					{
					   if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				     {
					      //WWqq category
						    MSPlot["MS_St_boxWWqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);		
						 }							
					}
				}
				else if(SelectednEl == 3)
				{
				   if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				   {
				      //WZqq category
					    MSPlot["MS_nEvts_boxWZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
							MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);	
					 }							
				}
				else if(SelectednEl == 4)
				{
				   if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				   {
				      //ZZqq category
					    MSPlot["MS_nEvts_boxZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
							MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
					 }	
				}
				else
				{
				   cout<<"FOUND >= 5-lepton event!"<<endl;		
				}				
			 }
			 else if(SelectednMu == 1)
			 {
			  if(SelectednEl == 1)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
			       //WWqq category (no need to check the invariant mass of the leptons, they should not come from a Z...)
					   MSPlot["MS_St_boxWWqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
					}	
				}
				else if(SelectednEl == 2)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
			       //WZqq category
					   MSPlot["MS_nEvts_boxWZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						 MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
					}
				}
				/*else if(SelectednEl == 3)
				{
				  if(selectedJets[0].Pt() >= LeadingJetPtCut)
				  {
			       //ZZqq category
					   MSPlot["MS_nEvts_boxZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						 MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
					}
				}*/
				else
				{
				   cout<<"FOUND >= 5-lepton event!"<<endl;		
				}									
			 }
			 else if(SelectednMu == 2)
			 {			  
				if(SelectednEl == 1)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
			      //WZqq category
					  MSPlot["MS_nEvts_boxWZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);	
						MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);	
				  }			
				}
				else if(SelectednEl == 2)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
			      //ZZqq category
					  MSPlot["MS_nEvts_boxZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
					}
				}
				else 
				{
				   cout<<"FOUND >= 5-lepton event!"<<endl;		
				}
			 }
			 else if(SelectednMu == 3)
			 {
				/*if(SelectednEl == 1)
				{
				  if(selectedJets[0].Pt() >= LeadingJetPtCut)
				  {
			      //ZZqq category
					  MSPlot["MS_nEvts_boxZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
					}
				}
				else 
				{
				   cout<<"FOUND >= 5-lepton event!"<<endl;		
				}*/
			 }
		 
			  
			}	
			

		 
	//delete selection;
   // cout << "done processing event" << endl;
		}//loop on events

    ////to test 'purity' of good jet combinations with a 'simple' method (not MVA). For the moment not supported
    //myInclFourthGenSearchTools.PrintPurityGoodCombinations();
    
    cout<<endl;
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
		/*if (!useMassesAndResolutions && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon) ||  (dataSetName.find("TTbarJets_SemiElectron") == 0 && semiElectron)))
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
    }*/
		
		inFile->Close();
		delete inFile;

  } //loop on 'datasets'
	
	if(systematic=="Nominal") myfile1.close();
	if(systematic=="Nominal") myfile2.close();

  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;
  fout->cd();
  


    fout->cd();
    //Write histograms: MSPlots
    if(savePNG) mkdir((pathPNG+"MSPlot/").c_str(),0777);
		cout << "Running over all MS plots" << endl;
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        temp->Draw(false, name, true, true, true, true, true,5,false, true, true, true, true);//(bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal, bool addRatio, bool mergeVV, bool mergeTTV, bool mergeVVV, bool mergeSameSignWW)
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
    selecTableMu.TableCalculator(true, true, true, true, true, true, true, true);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
    string selectiontableMu = "VLQSearch_SelectionTable_TreeAnalyzer_Mu";
    selectiontableMu = Outputpath+selectiontableMu+outputpostfix+".tex"; 	
    selecTableMu.Write(selectiontableMu.c_str(),false, true, false, false, false, false, false); //(filename, error, merged, lines, unscaled, eff, totaleff, landscape)
   /* 
		selecTableMultiLep.TableCalculator(true, true, true, true, true, false, true, true);
    string selectiontableMultiLep = "InclFourthGenSearch_SelectionTable_TreeAnalyzer_MultiLepton";
    selectiontableMultiLep = Outputpath+selectiontableMultiLep+outputpostfix+".tex"; 	
    selecTableMultiLep.Write(selectiontableMultiLep.c_str(),false, true, false, false, false, false, false);
		*/
     
    fout->cd();
  
	delete fout;
	
	delete leptonTools;

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

unsigned int nBjetsPresent(vector < float > btagvalues, float btagWP)
{
   unsigned int nBjets = 0;
   for(unsigned int iJet = 0; iJet < btagvalues.size(); iJet++)
	 {
	    if(btagvalues[iJet] >= btagWP)
			{
					 nBjets++; //because the current jet would be b tagged
			}
	 }
	 return nBjets;
};

unsigned int nNonBjetsPresent(vector < float > btagvalues, float btagWP)
{
   unsigned int nNonBjets = 0;
   for(unsigned int iJet = 0; iJet < btagvalues.size(); iJet++)
	 {
	    if(btagvalues[iJet] < btagWP)
			{
					 nNonBjets++; //because the current jet would not be b tagged
			}
	 }
	 return nNonBjets;
};
