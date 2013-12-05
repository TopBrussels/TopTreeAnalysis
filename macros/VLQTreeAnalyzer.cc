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
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
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
} pair_decrease;

struct sort_decreasing {
  bool operator() (float a,float b) { return (a>b);}
} decrease;

//To cout the Px, Py, Pz, E and Pt of objects
void coutObjectsFourVector(vector < TRootMuon* > init_muons, vector < TRootElectron* > init_electrons, vector < TRootJet* > init_jets, vector < TRootMET* > mets, string Comment);

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/*unsigned int nBjetsPresent(vector < float > btagvalues, float btagWP);

unsigned int nNonBjetsPresent(vector < float > btagvalues, float btagWP);
*/
//float calculateMEz(TLorentzVector lepton, TLorentzvector metvector);

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
  string btagger = "CSVM", antibtagger = "CSVL"; //string btagger = "CSVM", antibtagger = "CSVL"
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

  //string Treespath = "VLQTrees_Summer12_PBS_24June2013"; //v4
	string Treespath = "VLQTrees_Summer12_PBS_v5_21Aug13";
  Treespath = Treespath + "/"; 		
	bool savePNG = false;
	string outputpostfix = "";
	outputpostfix = outputpostfix+"_"+systematic;
	string Outputpath = "OutputFiles_VLQAnalyzer_withoutQCDestimation_5Dec13_LeptonSF_signal700";
	Outputpath = Outputpath + "/";
	mkdir(Outputpath.c_str(),0777);
		
	
  /////////////////////
  // Configuration
  /////////////////////

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
	
	bool doPUreweighting = true;
	bool verbosityBox = false; //to cout in which signal box you end up
	float LeadingJetPtCut = 200., SubLeadingJetPtCut = 100.; //subleading jet pt cut only for the pair production signal boxes...
	bool doBtagging = true; //Higgs->bbar decay
	bool doAntiBtagging = true; //all other jets
	bool applyBtagSF = true;
	bool useBtagInfoFromFile = true;
	int syst_btag = 0; //-1 for SF up, +1 for SF down
	bool doMETcutboxWqq = false; //would only be in Wqq categories (i.e. single lepton with only 1 high-Pt jet)
	float METcutboxWqq = 60.; //problem: QCD gets mostly rejected with this, but it is also cutting away quite a lot of signal!!
	bool doMTcut = true;//would only be in Wqq categories (i.e. single lepton with only 1 high-Pt jet), instead of MET cut...
	float MTcut = 40.;
	bool doMETcutboxWWqq = true; //would only be in WWqq categories (i.e. double lepton)
	float METcutboxWWqq = 60.;
	bool applyLeptonSF = true;
	
	//particular signal model (not only for the overlay mass point, but for all mass points, for consistency. Make sure you don't forget this in the scanning later on!!)
	bool doModelScaling = true; //to reweight according to given Branching Fractions and kappatilde values or not
	float BF_W = 0.5, BF_Z = 0.25;
	float BF_H = 1.0 - BF_W - BF_Z;
	float kappatilde_W = 1.;
	float kappatilde_Z = kappatilde_W * sqrt( 2 * BF_Z / BF_W);  
	//normalization factors for the different signal samples
	//pair production
	float F_QQTojWjW = BF_W*BF_W;
	float F_QQTojWjZ = 2*BF_W*BF_Z;
	float F_QQTojWjH = 2*BF_W*BF_H;
	float F_QQTojZjZ = BF_Z*BF_Z;
	float F_QQTojZjH = 2*BF_Z*BF_H;
	//single EW production
	float F_QTojW_CC = kappatilde_W*kappatilde_W*BF_W;
	float F_QTojW_NC = kappatilde_Z*kappatilde_Z*BF_W;
	float F_QTojZ_CC = kappatilde_W*kappatilde_W*BF_Z;
	float F_QTojZ_NC = kappatilde_Z*kappatilde_Z*BF_Z;
	
	
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
	
//  if (foundMu)
//	{
      for (unsigned int d = 0; d < datasets.size (); d++)
			{
	       string dataSetName = datasets[d]->Name();
	       if ( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) && !(dataSetName.find("InvIso_El") != string::npos)) 
	          datasetsMu.push_back(datasets[d]);
      }
//  }
      
//  if (foundEl)
//	{
    for (unsigned int d = 0; d < datasets.size (); d++)
		{
        string dataSetName = datasets[d]->Name();
        if ( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) && !(dataSetName.find("InvIso_Mu") != string::npos)) 
	          datasetsEl.push_back(datasets[d]);
    }
//  }
	
	
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
	MSPlot["MS_nPV_noPUreweighting_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices before PU reweighting", 36, -0.5, 35.5, "Nr. of primary vertices");
  MSPlot["MS_nPV_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices", 36, -0.5, 35.5, "Nr. of primary vertices");  //after PU reweighting
	MSPlot["MS_nPV_noPUreweighting_4jets_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices before PU reweighting", 36, -0.5, 35.5, "Nr. of primary vertices");
  MSPlot["MS_nPV_4jets_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices", 36, -0.5, 35.5, "Nr. of primary vertices");  //after PU reweighting
	MSPlot["MS_nLeptons_mu"] = new MultiSamplePlot(datasetsMu,"number of leptons", 6, -0.5, 5.5, "number of leptons");
	MSPlot["MS_nMuons_mu"] = new MultiSamplePlot(datasetsMu,"number of muons", 6, -0.5, 5.5, "number of muons");
	MSPlot["MS_nElectrons_mu"] = new MultiSamplePlot(datasetsMu,"number of electrons", 6, -0.5, 5.5, "number of electrons");
	MSPlot["MS_MET_mu"] = new MultiSamplePlot(datasetsMu,"MET", 200, 0, 1000, "Missing transverse energy (GeV)");
  MSPlot["MS_MET_singleMu_mu"] = new MultiSamplePlot(datasetsMu,"MET", 200, 0, 1000, "Missing transverse energy (GeV)");
	MSPlot["MS_METoverST_mu"] = new MultiSamplePlot(datasetsMu,"METoverST", 100, 0, 1, "Missing transverse energy / S_{T}");
	MSPlot["MS_MTleptonmet_preboxWqq_mu"] = new MultiSamplePlot(datasetsMu,"MT(lepton,MET)", 40, 0, 300, "Transverse mass(lepton,MET) (GeV)");
	MSPlot["MS_HT_mu"] = new MultiSamplePlot(datasetsMu,"HT", 200, 0, 2000, "HT (GeV)");
	//jets (not forward)
	MSPlot["MS_JetMultiplicity_mu"] = new MultiSamplePlot(datasetsMu, "JetMultiplicity", 10, -0.5, 9.5, "jet multiplicity");
	MSPlot["MS_JetMultiplicity_singleMu_mu"] = new MultiSamplePlot(datasetsMu, "JetMultiplicity", 10, -0.5, 9.5, "jet multiplicity");
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
	MSPlot["MS_Eta_allForwardjets_mu"] = new MultiSamplePlot(datasetsMu,"Eta_allForwardjets", 50, -5 , 5, "Eta of all forward jets");
	MSPlot["MS_Eta_Forwardjet1_mu"] = new MultiSamplePlot(datasetsMu,"Eta_Forwardjet1", 50, -5 , 5, "Eta of forward jet 1");
	MSPlot["MS_Eta_Forwardjet2_mu"] = new MultiSamplePlot(datasetsMu,"Eta_Forwardjet2", 50, -5 , 5, "Eta of forward jet 2");
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
	MSPlot["MS_St_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box Wqq)", 50, 200, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxWplusqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box W^{+}qq)", 40, 200, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxWminusqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box W^{-}qq)", 40, 200, 1500, "S_{T} (GeV)");
	MSPlot["MS_MQ_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Mass Q (box Wqq)", 40, 200, 1500, "Mass(lepton,neutrino,leadingjet) (GeV)");
	MSPlot["MS_MQ_boxWplusqq_mu"] = new MultiSamplePlot(datasetsMu,"Mass Q (box W^{+}qq)", 40, 200, 1500, "Mass(lepton,neutrino,leadingjet) (GeV)");
	MSPlot["MS_MQ_boxWminusqq_mu"] = new MultiSamplePlot(datasetsMu,"Mass Q (box W^{-}qq)", 40, 200, 1500, "Mass(lepton,neutrino,leadingjet) (GeV)");
	MSPlot["MS_MTQ_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Transverse mass Q (box Wqq)", 40, 200, 1500, "Transverse mass (GeV)");
	MSPlot["MS_MTQ_boxWplusqq_mu"] = new MultiSamplePlot(datasetsMu,"Transverse mass Q (box W^{+}qq)", 40, 200, 1500, "Transverse mass (GeV)");
	MSPlot["MS_MTQ_boxWminusqq_mu"] = new MultiSamplePlot(datasetsMu,"Transverse mass Q (box W^{-}qq)", 40, 200, 1500, "Transverse mass (GeV)");
	MSPlot["MS_St_boxWHqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box WHqq)", 50, 200, 1700, "S_{T} (GeV)");
	MSPlot["MS_MHq_boxWHqq_mu"] = new MultiSamplePlot(datasetsMu,"Mass Q (box WHqq)", 50, 0, 1500, "Mass(H,highptjet) (GeV)");
	MSPlot["MS_St_boxZqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box Zqq)", 40, 200, 1500, "S_{T} (GeV)");
	MSPlot["MS_MQ_boxZqq_mu"] = new MultiSamplePlot(datasetsMu,"Mass Q (box Zqq)", 40, 200, 1500, "Mass(lepton,lepton,leadingjet) (GeV)");
	MSPlot["MS_St_boxZHqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box ZHqq)", 15, 200, 1700, "S_{T} (GeV)");
	MSPlot["MS_MQ_boxZHqq_mu"] = new MultiSamplePlot(datasetsMu,"Mass Q (box ZHqq)", 15, 0, 1500, "Mass(lepton,lepton,highptjet) (GeV)");
	////MSPlot["MS_MHq_boxZHqq_mu"] = new MultiSamplePlot(datasetsMu,"Mass Hq (box ZHqq)", 15, 0, 1500, "Mass(H,highptjet) (GeV)");
	MSPlot["MS_St_boxWWqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box WWqq)", 15, 200, 1700, "S_{T} (GeV)");
//	MSPlot["MS_AveragePseudoMQ_boxWWqq_mu"] = new MultiSamplePlot(datasetsMu,"Average Pseudomass Q (box WWqq)", 40, 200, 1500, "Average Pseudomass(lepton,neutrino,jet) (GeV)");
	MSPlot["MS_nEvts_boxWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxZZqq_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box ZZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqqZZqq_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq+ZZqq)", 1, 0.5, 1.5, "");
	//btag or antibtag MC event weights in 'boxes'
	MSPlot["MS_BtagWeight_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"BtagWeight (box Wqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxWHqq_mu"] = new MultiSamplePlot(datasetsMu,"BtagWeight (box WHqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxZqq_mu"] = new MultiSamplePlot(datasetsMu,"BtagWeight (box Zqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxZHqq_mu"] = new MultiSamplePlot(datasetsMu,"BtagWeight (box ZHqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxWWqq_mu"] = new MultiSamplePlot(datasetsMu,"BtagWeight (box WWqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxWZqq_mu"] = new MultiSamplePlot(datasetsMu,"BtagWeight (box WZqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxZZqq_mu"] = new MultiSamplePlot(datasetsMu,"BtagWeight (box ZZqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxWZqqZZqq_mu"] = new MultiSamplePlot(datasetsMu,"BtagWeight (box WZqq+ZZqq)", 80, 0.8, 1.2, "weight");
	//plots in boxes
	MSPlot["MS_MHcandidate_boxWHqq_mu"] = new MultiSamplePlot(datasetsMu,"MH (box WHqq)", 50, 0, 500, "reconstructed Higgs mass (GeV)");
	////MSPlot["MS_Dijetmass_boxWHqq_mu"] = new MultiSamplePlot(datasetsMu,"Dijetmass (box WHqq)", 50, 0, 500, "Dijet mass (GeV)");
	////MSPlot["MS_Dibjetmass_boxWHqq_mu"] = new MultiSamplePlot(datasetsMu,"Dibjetmass (box WHqq)", 50, 0, 500, "Dibjet mass (GeV)");
	////MSPlot["MS_Dijet34mass_boxWHqq_mu"] = new MultiSamplePlot(datasetsMu,"Dijet34mass (box WHqq)", 50, 0, 500, "Dijet34 mass (GeV)"); //jet 3 and 4 meaning the third and fourth jet (sorted in pt)
	////MSPlot["MS_nJets_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"nJets (box Wqq)", 9, -0.5, 8.5, "Number of jets (GeV)");
	////MSPlot["MS_nJets_preboxWqq_noBveto_mu"] = new MultiSamplePlot(datasetsMu,"nJets (prebox Wqq)", 9, -0.5, 8.5, "Number of Jets (GeV)");
	MSPlot["MS_MET_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"MET (box Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	MSPlot["MS_MET_preboxWqq_noleadingjetcut_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	MSPlot["MS_MET_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	MSPlot["MS_MET_preboxWqq_noBveto_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	
	//'el channel' = when NO muon is trigged, and at least one electron
	MSPlot["MS_nPV_noPUreweighting_el"] = new MultiSamplePlot(datasetsEl, "nPrimaryVertices before PU reweighting", 36, -0.5, 35.5, "Nr. of primary vertices");
  MSPlot["MS_nPV_el"] = new MultiSamplePlot(datasetsEl, "nPrimaryVertices", 36, -0.5, 35.5, "Nr. of primary vertices");  //after PU reweighting
	MSPlot["MS_nPV_noPUreweighting_4jets_el"] = new MultiSamplePlot(datasetsEl, "nPrimaryVertices before PU reweighting", 36, -0.5, 35.5, "Nr. of primary vertices");
  MSPlot["MS_nPV_4jets_el"] = new MultiSamplePlot(datasetsEl, "nPrimaryVertices", 36, -0.5, 35.5, "Nr. of primary vertices");  //after PU reweighting
	MSPlot["MS_nLeptons_el"] = new MultiSamplePlot(datasetsEl,"number of leptons", 6, -0.5, 5.5, "number of leptons");
	MSPlot["MS_nMuons_el"] = new MultiSamplePlot(datasetsEl,"number of muons", 6, -0.5, 5.5, "number of muons");
	MSPlot["MS_nElectrons_el"] = new MultiSamplePlot(datasetsEl,"number of electrons", 6, -0.5, 5.5, "number of electrons");
	MSPlot["MS_MET_el"] = new MultiSamplePlot(datasetsEl,"MET", 200, 0, 1000, "Missing transverse energy (GeV)");  
	MSPlot["MS_MET_singleEl_el"] = new MultiSamplePlot(datasetsEl,"MET", 200, 0, 1000, "Missing transverse energy (GeV)");
	MSPlot["MS_METoverST_el"] = new MultiSamplePlot(datasetsEl,"METoverST", 100, 0, 1, "Missing transverse energy / S_{T}");
	MSPlot["MS_MTleptonmet_preboxWqq_el"] = new MultiSamplePlot(datasetsEl,"MT(lepton,MET)", 40, 0, 300, "Transverse mass(lepton,MET) (GeV)");
	MSPlot["MS_HT_el"] = new MultiSamplePlot(datasetsEl,"HT", 200, 0, 2000, "HT (GeV)");
	//jets (not forward)
	MSPlot["MS_JetMultiplicity_el"] = new MultiSamplePlot(datasetsEl, "JetMultiplicity", 10, -0.5, 9.5, "jet multiplicity");
	MSPlot["MS_JetMultiplicity_singleEl_el"] = new MultiSamplePlot(datasetsEl, "JetMultiplicity", 10, -0.5, 9.5, "jet multiplicity");
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
	MSPlot["MS_Eta_allForwardjets_el"] = new MultiSamplePlot(datasetsEl,"Eta_allForwardjets", 50, -5 , 5, "Eta of all forward jets");
	MSPlot["MS_Eta_Forwardjet1_el"] = new MultiSamplePlot(datasetsEl,"Eta_Forwardjet1", 50, -5 , 5, "Eta of forward jet 1");
	MSPlot["MS_Eta_Forwardjet2_el"] = new MultiSamplePlot(datasetsEl,"Eta_Forwardjet2", 50, -5 , 5, "Eta of forward jet 2");
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
	MSPlot["MS_St_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box Wqq)", 40, 200, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxWplusqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box W^{+}qq)", 40, 200, 1500, "S_{T} (GeV)");
	MSPlot["MS_St_boxWminusqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box W^{-}qq)", 40, 200, 1500, "S_{T} (GeV)");
	MSPlot["MS_MQ_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Mass Q (box Wqq)", 40, 200, 1500, "Mass(lepton,neutrino,leadingjet) (GeV)");
	MSPlot["MS_MQ_boxWplusqq_el"] = new MultiSamplePlot(datasetsEl,"Mass Q (box W^{+}qq)", 40, 200, 1500, "Mass(lepton,neutrino,leadingjet) (GeV)");
	MSPlot["MS_MQ_boxWminusqq_el"] = new MultiSamplePlot(datasetsEl,"Mass Q (box W^{-}qq)", 40, 200, 1500, "Mass(lepton,neutrino,leadingjet) (GeV)");
	MSPlot["MS_MTQ_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Transverse mass Q (box Wqq)", 40, 200, 1500, "Transverse mass (GeV)");
	MSPlot["MS_MTQ_boxWplusqq_el"] = new MultiSamplePlot(datasetsEl,"Transverse mass Q (box W^{+}qq)", 40, 200, 1500, "Transverse mass (GeV)");
	MSPlot["MS_MTQ_boxWminusqq_el"] = new MultiSamplePlot(datasetsEl,"Transverse mass Q (box W^{-}qq)", 40, 200, 1500, "Transverse mass (GeV)");
	MSPlot["MS_St_boxWHqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box WHqq)", 50, 200, 1700, "S_{T} (GeV)");
	MSPlot["MS_MHq_boxWHqq_el"] = new MultiSamplePlot(datasetsEl,"Mass Q (box WHqq)", 50, 0, 1500, "Mass(H,highptjet) (GeV)");
	MSPlot["MS_St_boxZqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box Zqq)", 40, 200, 1500, "S_{T} (GeV)");
	MSPlot["MS_MQ_boxZqq_el"] = new MultiSamplePlot(datasetsEl,"Mass Q (box Zqq)", 40, 200, 1500, "Mass(lepton,lepton,leadingjet) (GeV)");
	MSPlot["MS_St_boxZHqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box ZHqq)", 15, 200, 1700, "S_{T} (GeV)");
	MSPlot["MS_MQ_boxZHqq_el"] = new MultiSamplePlot(datasetsEl,"Mass Q (box ZHqq)", 15, 0, 1500, "Mass(lepton,lepton,highptjet) (GeV)");
	////MSPlot["MS_MHq_boxZHqq_el"] = new MultiSamplePlot(datasetsEl,"Mass Hq (box ZHqq)", 15, 0, 1500, "Mass(H,highptjet) (GeV)");
	MSPlot["MS_St_boxWWqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box WWqq)", 15, 200, 1700, "S_{T} (GeV)");
	MSPlot["MS_nEvts_boxWZqq_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxZZqq_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box ZZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqqZZqq_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq+ZZqq)", 1, 0.5, 1.5, "");
	//btag or antibtag MC event weights in 'boxes'
	MSPlot["MS_BtagWeight_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"BtagWeight (box Wqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxWHqq_el"] = new MultiSamplePlot(datasetsEl,"BtagWeight (box WHqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxZqq_el"] = new MultiSamplePlot(datasetsEl,"BtagWeight (box Zqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxZHqq_el"] = new MultiSamplePlot(datasetsEl,"BtagWeight (box ZHqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxWWqq_el"] = new MultiSamplePlot(datasetsEl,"BtagWeight (box WWqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxWZqq_el"] = new MultiSamplePlot(datasetsEl,"BtagWeight (box WZqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxZZqq_el"] = new MultiSamplePlot(datasetsEl,"BtagWeight (box ZZqq)", 80, 0.8, 1.2, "weight");
	MSPlot["MS_BtagWeight_boxWZqqZZqq_el"] = new MultiSamplePlot(datasetsEl,"BtagWeight (box WZqq+ZZqq)", 80, 0.8, 1.2, "weight");
  //plots in boxes
	MSPlot["MS_MHcandidate_boxWHqq_el"] = new MultiSamplePlot(datasetsEl,"MH (box WHqq)", 50, 0, 500, "reconstructed Higgs mass (GeV)");
	////MSPlot["MS_Dijetmass_boxWHqq_el"] = new MultiSamplePlot(datasetsEl,"Dijetmass (box WHqq)", 50, 0, 1700, "Dijet mass (GeV)");
	////MSPlot["MS_Dibjetmass_boxWHqq_el"] = new MultiSamplePlot(datasetsEl,"Dibjetmass (box WHqq)", 50, 0, 1700, "Dibjet mass (GeV)");
	////MSPlot["MS_Dijet34mass_boxWHqq_el"] = new MultiSamplePlot(datasetsEl,"Dijet34mass (box WHqq)", 50, 0, 1700, "Dijet34 mass (GeV)"); //jet 3 and 4 meaning the third and fourth jet (sorted in pt)
	////MSPlot["MS_nJets_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"nJets (box Wqq)", 9, -0.5, 8.5, "Number of jets (GeV)");
	////MSPlot["MS_nJets_preboxWqq_noBveto_el"] = new MultiSamplePlot(datasetsEl,"nJets (prebox Wqq)", 9, -0.5, 8.5, "Number of Jets (GeV)");
	MSPlot["MS_MET_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"MET (box Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	MSPlot["MS_MET_preboxWqq_noleadingjetcut_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	MSPlot["MS_MET_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	MSPlot["MS_MET_preboxWqq_noBveto_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");

	
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
	CutsSelecTableMu.push_back(string("two leptons"));
	CutsSelecTableMu.push_back(string("#jets>=2"));
	CutsSelecTableMu.push_back(string("leading and subleading jet Pt cuts"));
	CutsSelecTableMu.push_back(string("b-veto"));

	
	
	

  vector<string> CutsSelecTableMultiLep;
  CutsSelecTableMultiLep.push_back(string("SS leptons: all")); //0
  CutsSelecTableMultiLep.push_back(string("SS leptons: 2 muons")); //1
  CutsSelecTableMultiLep.push_back(string("SS leptons: 2 electrons")); //2
  CutsSelecTableMultiLep.push_back(string("SS leptons: electron+muon")); //3
  CutsSelecTableMultiLep.push_back(string("trileptons in all boxes combined")); //4
  
	SelectionTable selecTableLep(CutsSelecTableLep, datasets);
  selecTableLep.SetLuminosity(Luminosity);
  selecTableLep.SetPrecision(1);
	
	SelectionTable selecTableMu(CutsSelecTableMu, datasets); //or datasets?? to be checked!!
  selecTableMu.SetLuminosity(LuminosityMu);
  selecTableMu.SetPrecision(1);
	

	SelectionTable selecTableMultiLep(CutsSelecTableMultiLep, datasets);
  selecTableMultiLep.SetLuminosity(Luminosity);
  selecTableMultiLep.SetPrecision(2);
  
  cout << " - SelectionTable instantiated ..." << endl;

  
	////////////////////////
  // PileUp Reweighting //
  ////////////////////////

  //cout << Luminosity << endl;

	string MCpileuphistofile = "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root";
//  if ( dataSetName.find("NP_") == 0 && dataSetName.find("fullsim") != 0)
	string MCpileuphistofile_S7 = "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Spring12_S7.root";

  //For Run2012A+B_C_D (rereco) SingleMu
  LumiReWeighting LumiWeights_mu = LumiReWeighting(MCpileuphistofile, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleMu/nominal.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsUp_mu = LumiReWeighting(MCpileuphistofile, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleMu/sys_up.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsDown_mu = LumiReWeighting(MCpileuphistofile, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleMu/sys_down.root", "pileup", "pileup");
  //temporary fix for fastsim samples produced with the Spring12 PU, which is an S7 scenario, not S10 like in official Summer12
	LumiReWeighting LumiWeights_S7_mu = LumiReWeighting(MCpileuphistofile_S7, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleMu/nominal.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsUp_S7_mu = LumiReWeighting(MCpileuphistofile_S7, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleMu/sys_up.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsDown_S7_mu = LumiReWeighting(MCpileuphistofile_S7, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleMu/sys_down.root", "pileup", "pileup");
 
	//For Run2012A+B_C_D (rereco) SingleElectron
  LumiReWeighting LumiWeights_el = LumiReWeighting(MCpileuphistofile, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleElectron/nominal.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsUp_el = LumiReWeighting(MCpileuphistofile, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleElectron/sys_up.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsDown_el = LumiReWeighting(MCpileuphistofile, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleElectron/sys_down.root", "pileup", "pileup");
  //temporary fix for fastsim samples produced with the Spring12 PU, which is an S7 scenario, not S10 like in official Summer12
	LumiReWeighting LumiWeights_S7_el = LumiReWeighting(MCpileuphistofile_S7, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleElectron/nominal.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsUp_S7_el = LumiReWeighting(MCpileuphistofile_S7, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleElectron/sys_up.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsDown_S7_el = LumiReWeighting(MCpileuphistofile_S7, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleElectron/sys_down.root", "pileup", "pileup");
 

  cout << " Initialized LumiReWeighting stuff " << endl;
	cout << "     (used " << MCpileuphistofile << " in MC PU reweighting)" << endl;
	  
		
		
		
		
		
		
		
  ////////////////////////////////////
	//  Initialise BTag ScaleFactors  //
	////////////////////////////////////
	//will use method 1a) of https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
	
	BTagWeightTools *bTool_L; //for creation of MC eff histos and application of SF for anti-b-tagging...
	BTagWeightTools *bTool_T; //for creation of MC eff histos and application of SF for b-tagging for H->bbar...
	BTagWeightTools *bTool_LT; //for application of SF (using formerly created MC eff histos) for anti-b-tagging, and b-tagging for H->bbar...
	
	if(!useBtagInfoFromFile)
	{
	   bTool_L = new BTagWeightTools("../../TopTreeAnalysisBase/Calibrations/BTagging/SFb-pt_WITHttbar_payload_EPS13_noCSVSL.txt",antibtagger,"EPS2013"); //for creation of MC eff histos for anti-b-tagging...
	   bTool_T = new BTagWeightTools("../../TopTreeAnalysisBase/Calibrations/BTagging/SFb-pt_WITHttbar_payload_EPS13_noCSVSL.txt",btagger,"EPS2013"); //for creation of MC eff histos for b-tagging for H->bbar...
		 bTool_L->InitializeMCEfficiencyHistos(15,30.,340.,2); //NofPtBins,PtMin,PtMax,NofEtaBins
		 bTool_T->InitializeMCEfficiencyHistos(15,30.,340.,2); //NofPtBins,PtMin,PtMax,NofEtaBins
	}
	else if(useBtagInfoFromFile)
	{
	   bTool_L = new BTagWeightTools("../../TopTreeAnalysisBase/Calibrations/BTagging/SFb-pt_WITHttbar_payload_EPS13_noCSVSL.txt",antibtagger,"EPS2013"); //for anti-b-tagging, with application of Sfs and created MC eff histos...
	   bTool_L->ReadMCEfficiencyHistos("PlotsForBtagWeights_CSVL_ttbar.root");
		 //bTool_T = new BTagWeightTools("../../TopTreeAnalysisBase/Calibrations/BTagging/SFb-pt_WITHttbar_payload_EPS13_noCSVSL.txt",btagger,"EPS2013");  //for b-tagging for H->bbar with application of Sfs and created MC eff histos... never needed without b-tagging?
	   //bTool_T->ReadMCEfficiencyHistos("PlotsForBtagWeights_btag.root");
		 bTool_LT = new BTagWeightTools("../../TopTreeAnalysisBase/Calibrations/BTagging/SFb-pt_WITHttbar_payload_EPS13_noCSVSL.txt",antibtagger,btagger,"EPS2013"); //for anti-b-tagging, and b-tagging for H->bbar with application of Sfs and created MC eff histos...
	   bTool_LT->ReadMCEfficiencyHistos("PlotsForBtagWeights_CSVL_ttbar.root","PlotsForBtagWeights_CSVM_ttbar.root");
	}
		
		
		
  // initialize lepton SF
  LeptonTools* leptonTools = new LeptonTools(false); //boolean indicates verbosity
  //leptonTools->readMuonSF("../../TopTreeAnalysisBase/Calibrations/LeptonSF/Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root", "../../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonEfficiencies_Run_2012A_2012B_53X.root", "../../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonEfficiencies_Run_2012C_53X.root", "../../TopTreeAnalysisBase/Calibrations/LeptonSF/TriggerMuonEfficiencies_Run_2012D_53X.root"); //old way (promptreco)
  leptonTools->readMuonSF("../../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonEfficiencies_Run2012ReReco_53X.root","../../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonEfficiencies_ISO_Run_2012ReReco_53X.root","../../TopTreeAnalysisBase/Calibrations/LeptonSF/SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root");
	leptonTools->readElectronSF();
    
/*		
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
*/	
		
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
		//if(inFile->IsZombie() && )
		
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
		
		
		
		
		/*
		//////////////////WARNING: JUST  FOR Vjets CHECKS!!!!!
		if((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)) end = 2000;
		///////////////////////////////////////////////////////
		*/
		
		
		
		
    for (int ievt = start; ievt < end; ievt++)
		//for (int ievt = start; ievt < 5000; ievt++)
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
						
			TLorentzVector metvector = myBranch_selectedEvents->met();
			float met = (myBranch_selectedEvents->met()).Pt();  //Pt instead of Et, because the Et is not corrected when applying JEC or systematic shifts
      vector<TLorentzVector> selectedJets = myBranch_selectedEvents->selectedJets();
			vector<int> selectedJets_partonFlavour = myBranch_selectedEvents->selectedJets_partonFlavour();
			vector<TLorentzVector> selectedForwardJets = myBranch_selectedEvents->selectedForwardJets();
			vector<TLorentzVector> selectedMuons = myBranch_selectedEvents->selectedMuons();
			vector<TLorentzVector> selectedElectrons = myBranch_selectedEvents->selectedElectrons();
			vector<float> selectedJets_bTagCSV = myBranch_selectedEvents->selectedJets_bTagCSV();
			vector<int> selectedMuonsCharge = myBranch_selectedEvents->selectedMuonsCharge();	
			vector<int> selectedElectronsCharge = myBranch_selectedEvents->selectedElectronsCharge();	
				   
					 
			//cout<<"event "<<ievt<<" with id = "<<myBranch_selectedEvents->eventID()<<", met = "<<met<<endl;		 
					 
					 
			//temporary (dirty) fix to make sure selectedJets only has jets with |eta| < 2.4, best way to fix is in treecreator, but then trees need to be rerun
			vector<TLorentzVector> selectedJets_tmp;
			vector<float> selectedJets_bTagCSV_tmp;
			vector<int> selectedJets_partonFlavour_tmp;
			if(selectedJets_bTagCSV.size()!=selectedJets.size() || selectedJets_partonFlavour.size()!=selectedJets.size())
			{
			   cout<<"ERROR: selectedJets_bTagCSV.size() != selectedJets.size()"<<endl;
				 exit(0);
		  }
			for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
			{
				   if(fabs(selectedJets[iJet].Eta())<2.4)
					 {
					    selectedJets_tmp.push_back(selectedJets[iJet]);
							selectedJets_bTagCSV_tmp.push_back(selectedJets_bTagCSV[iJet]);	
						  selectedJets_partonFlavour_tmp.push_back(selectedJets_partonFlavour[iJet]);	
					 }
					 else
					 {
					    selectedForwardJets.push_back(selectedJets[iJet]); //WARNING: pt ordering might be screwed up!! (but maybe not important)
							//btag value or flavor not applicable or needed...						
					 }
			}
			selectedJets.clear();
			selectedJets_bTagCSV.clear();
			selectedJets_partonFlavour.clear();			
			selectedJets = selectedJets_tmp;
			selectedJets_bTagCSV = selectedJets_bTagCSV_tmp;
			selectedJets_partonFlavour = selectedJets_partonFlavour_tmp;
			//end temporary (dirty) fix
					 
	 
	
	    vector<TLorentzVector> selectedBtaggedJets, selectedAntiBtaggedJets, selectedOthertaggedJets;
			vector<int> selectedAntiBtaggedJets_indices;
			for(int iJet = 0; iJet < selectedJets.size(); iJet++)
			{
			   if(selectedJets_bTagCSV[iJet] > btagWP)
				    selectedBtaggedJets.push_back(selectedJets[iJet]);
				 else if((0 <= selectedJets_bTagCSV[iJet]) && (selectedJets_bTagCSV[iJet] < antibtagWP))
				 {
				    selectedAntiBtaggedJets.push_back(selectedJets[iJet]);
						//cout<<" found antitagged jet: pushing back "<<iJet<<endl;
						selectedAntiBtaggedJets_indices.push_back(iJet);
				 }
				 else  
				    selectedOthertaggedJets.push_back(selectedJets[iJet]); //also non-taggable jets? (discr value -1 or -10)				 
			}
			

	 
	 
	 
      // scale factor for the event
      float scaleFactor = 1.;
      float lumiWeight = 1., lumiWeightUp = 1., lumiWeightDown = 1.;			
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("InvIso") != string::npos) )
			{
			  if(MuTrigged)
				{	
				  if ( dataSetName.find("NP_") == 0 && dataSetName.find("fullsim") == string::npos && dataSetName.find("NP_TprimeTprime") == string::npos)
					{
					   lumiWeight = LumiWeights_S7_mu.ITweight( (int) nTruePU );
             lumiWeightUp = LumiWeightsUp_S7_mu.ITweight( (int) nTruePU );
             lumiWeightDown = LumiWeightsDown_S7_mu.ITweight( (int) nTruePU );
					}
					else
					{			
				     lumiWeight = LumiWeights_mu.ITweight( (int) nTruePU );
             lumiWeightUp = LumiWeightsUp_mu.ITweight( (int) nTruePU );
             lumiWeightDown = LumiWeightsDown_mu.ITweight( (int) nTruePU );
					}
				}
				else if(ElTrigged)
				{
				  if ( dataSetName.find("NP_") == 0 && dataSetName.find("fullsim") == string::npos && dataSetName.find("NP_TprimeTprime") == string::npos)
					{
					   lumiWeight = LumiWeights_S7_el.ITweight( (int) nTruePU );
             lumiWeightUp = LumiWeightsUp_S7_el.ITweight( (int) nTruePU );
             lumiWeightDown = LumiWeightsDown_S7_el.ITweight( (int) nTruePU );
					}
					else
					{
				     lumiWeight = LumiWeights_el.ITweight( (int) nTruePU );
             lumiWeightUp = LumiWeightsUp_el.ITweight( (int) nTruePU );
             lumiWeightDown = LumiWeightsDown_el.ITweight( (int) nTruePU );
					}
				}
			
			}
			
			
			
      if((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos))
			{
	       lumiWeight=1.; //safe, not needed
				 
				/* stringstream r; r<<myBranch_selectedEvents->runID();
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
         }*/
			}
      
      // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int) myBranch_selectedEvents->nTruePU() );
      // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int) myBranch_selectedEvents->nTruePU() );

      if(doPUreweighting) scaleFactor = scaleFactor*lumiWeight;
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



      //signal sample scale factors to scale with particular point in the BF triangle (and kappa parameter)
			if(doModelScaling && (dataSetName.find("NP") == 0))
			{			
        if(dataSetName.find("QDQDTojWjW") != string::npos || dataSetName.find("NP_TprimeTprimeToBWBWinc") != string::npos)
			  {
			     scaleFactor = scaleFactor*F_QQTojWjW;
			  }
			  else if(dataSetName.find("QDQDTojWjZ") != string::npos)
			  {
			     scaleFactor = scaleFactor*F_QQTojWjZ;
			  }
			  else if(dataSetName.find("QDQDTojWjH") != string::npos)
			  {
			     scaleFactor = scaleFactor*F_QQTojWjH;
			  }
			  else if(dataSetName.find("QDQDTojZjZ") != string::npos)
			  {
			     scaleFactor = scaleFactor*F_QQTojZjZ;
			  }
			  else if(dataSetName.find("QDQDTojZjH") != string::npos)
			  {
			     scaleFactor = scaleFactor*F_QQTojZjH;
			  }
			  else if((dataSetName.find("QDTojW") != string::npos) && (dataSetName.find("CC") != string::npos))
			  {
			     scaleFactor = scaleFactor*F_QTojW_CC;
			  }
			  else if((dataSetName.find("QDTojW") != string::npos) && (dataSetName.find("NC") != string::npos))
			  {
			     scaleFactor = scaleFactor*F_QTojW_NC;
			  }
				else if((dataSetName.find("QDTojZ") != string::npos) && (dataSetName.find("CC") != string::npos))
			  {
			     scaleFactor = scaleFactor*F_QTojZ_CC;
			  }
			  else if((dataSetName.find("QDTojZ") != string::npos) && (dataSetName.find("NC") != string::npos))
			  {
			     scaleFactor = scaleFactor*F_QTojZ_NC;
			  }
			}

     
		  float HT = 0; //scalar sum of jet pt
      float ST = met; //test search variable			
			for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
			{
						    ST = ST + selectedJets[iJet].Pt();
								HT = ST + selectedJets[iJet].Pt();
			}			
			for(unsigned int iMu = 0; iMu < selectedMuons.size(); iMu++)
			{
						    ST = ST + selectedMuons[iMu].Pt();
			}
			for(unsigned int iEl = 0; iEl < selectedElectrons.size(); iEl++)
			{
						    ST = ST + selectedElectrons[iEl].Pt();
			}


      //'muon channel'
      //if(SelectednMu >= 1) //not anymore
			if(MuTrigged)
			{
			  //cout<<"**** Muon channel **** "<<endl;
			  if(!useBtagInfoFromFile)
				{
				   //to determine MC b-tag efficiency. Double-check if right place to put it
				   bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
        	 bTool_T->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);  
				}
				
									
			  //lepton scale factor; apply on MC
				if(!((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			  {								
				   if(applyLeptonSF)
					 {
					   //cout<<"(selectedMuons[0].Eta(), selectedMuons[0].Pt()) = ("<<selectedMuons[0].Eta()<<","<< selectedMuons[0].Pt()<<")"<<endl;
					   scaleFactor = scaleFactor*leptonTools->getMuonSF(selectedMuons[0].Eta(), selectedMuons[0].Pt()); //note: this is assuming the leading muon is the trigged one...
					   //cout<<"leptonTools->getMuonSF(selectedMuons[0].Eta(), selectedMuons[0].Pt()) = "<<leptonTools->getMuonSF(selectedMuons[0].Eta(), selectedMuons[0].Pt())<<endl;
						 for(unsigned int iMu = 0; iMu < selectedMuons.size(); iMu++)
			       {
					     if(iMu!=0 && (selectedMuons[iMu].Pt()>=20))
							 {
							    scaleFactor = scaleFactor*leptonTools->getMuonIdIsoSF(selectedMuons[iMu].Eta(), selectedMuons[iMu].Pt()); //apply id*iso scale factor for non-leading muons.
							    //cout<<"   leptonTools->getMuonIdIsoSF(selectedMuons["<<iMu<<"].Eta(), selectedMuons["<<iMu<<"].Pt()) = "<<leptonTools->getMuonIdIsoSF(selectedMuons[iMu].Eta(), selectedMuons[iMu].Pt())<<endl;
							 }
					   }
					   for(unsigned int iEl = 0; iEl < selectedElectrons.size(); iEl++)
			       {
					     if(selectedElectrons[iEl].Pt()>=20) 
							 {
							   scaleFactor = scaleFactor*leptonTools->getElectronIdIsoSF(selectedElectrons[iEl].Eta(), selectedElectrons[iEl].Pt()); //apply idiso scale factor for electrons
					       //cout<<"   leptonTools->getElectronIdIsoSF(selectedElectrons["<<iEl<<"].Eta(), selectedElectrons["<<iEl<<"].Pt()) = "<<leptonTools->getElectronIdIsoSF(selectedElectrons[iEl].Eta(), selectedElectrons[iEl].Pt())<<endl;
						   }
						 }
					 }					 
				}
				
				 
			  //MSPlot["MS_nPV_noPUreweighting_mu"]->Fill(nPV,datasets[d], true, LuminosityMu*scaleFactor / lumiWeight);	
			  MSPlot["MS_nPV_noPUreweighting_mu"]->Fill(nPV,datasets[d], true, LuminosityMu);	
			  MSPlot["MS_nPV_mu"]->Fill(nPV,datasets[d], true, LuminosityMu*scaleFactor);
				MSPlot["MS_nLeptons_mu"]->Fill(SelectednLeptons,datasets[d], true, LuminosityMu*scaleFactor);
	      MSPlot["MS_nMuons_mu"]->Fill(SelectednMu,datasets[d], true, LuminosityMu*scaleFactor);
	      MSPlot["MS_nElectrons_mu"]->Fill(SelectednEl,datasets[d], true, LuminosityMu*scaleFactor);
	      MSPlot["MS_MET_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);		
				MSPlot["MS_METoverST_mu"]->Fill(met/ST,datasets[d], true, LuminosityMu*scaleFactor);
				MSPlot["MS_HT_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);
				
				if(SelectednMu == 1 && SelectednEl == 0) 
				{
				   if(selectedJets[0].Pt() >= LeadingJetPtCut*0.5)
					 {
					    MSPlot["MS_MET_singleMu_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
					    MSPlot["MS_JetMultiplicity_singleMu_mu"]->Fill(selectedJets.size(),datasets[d], true, LuminosityMu*scaleFactor);
					 }
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
 
				//selecTableMu.Fill(d,0,scaleFactor);
				//if(selectedJets[0].Pt() >= LeadingJetPtCut) selecTableMu.Fill(d,1,scaleFactor);
			}
			
			//'electron channel'
			//if(SelectednMu == 0 && SelectednEl >= 1) //not anymore
			else if(ElTrigged)  //or an 'if'??
			{ 
			  //cout<<"**** Electron channel **** "<<endl;
			  if(!useBtagInfoFromFile)
				{
				   //to determine MC b-tag efficiency. Double-check if right place to put it
				   bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
        	 bTool_T->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);  
				}
				
			  //lepton scale factor; apply on MC
				if(!((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			  {
				   if(applyLeptonSF)
					 {
					   //cout<<"(selectedElectrons[0].Eta(), selectedElectrons[0].Pt()) = ("<<selectedElectrons[0].Eta()<<","<< selectedElectrons[0].Pt()<<")"<<endl;
				     scaleFactor = scaleFactor*leptonTools->getElectronSF(selectedElectrons[0].Eta(), selectedElectrons[0].Pt()); //note: this is assuming the leading electron is the trigged one...
					   //cout<<"leptonTools->getElectronSF(selectedElectrons[0].Eta(), selectedElectrons[0].Pt()) = "<<leptonTools->getElectronSF(selectedElectrons[0].Eta(), selectedElectrons[0].Pt())<<endl;
						 for(unsigned int iMu = 0; iMu < selectedMuons.size(); iMu++)
			       {
					     if(selectedMuons[iMu].Pt()>=20)
							 {
							    scaleFactor = scaleFactor*leptonTools->getMuonIdIsoSF(selectedMuons[iMu].Eta(), selectedMuons[iMu].Pt()); //apply id*iso scale factor for muons.
					   			//cout<<"  leptonTools->getMuonIdIsoSF(selectedMuons["<<iMu<<"].Eta(), selectedMuons["<<iMu<<"].Pt()) = "<<leptonTools->getMuonIdIsoSF(selectedMuons[iMu].Eta(), selectedMuons[iMu].Pt())<<endl;
               }
             }
					   for(unsigned int iEl = 0; iEl < selectedElectrons.size(); iEl++)
			       {
					     if(iEl!=0 && (selectedElectrons[iEl].Pt()>=20))
							 {
							    scaleFactor = scaleFactor*leptonTools->getElectronIdIsoSF(selectedElectrons[iEl].Eta(), selectedElectrons[iEl].Pt()); //apply idiso scale factor for non-leading electrons.
               		//cout<<"  leptonTools->getElectronIdIsoSF(selectedElectrons["<<iEl<<"].Eta(), selectedElectrons["<<iEl<<"].Pt()) = "<<leptonTools->getElectronIdIsoSF(selectedElectrons[iEl].Eta(), selectedElectrons[iEl].Pt())<<endl;
							 }
						 }
					 }
				}
				
			  //cout<<"d = "<<d<<endl;
			  MSPlot["MS_nPV_noPUreweighting_el"]->Fill(nPV,datasets[d], true, LuminosityEl);	
			  MSPlot["MS_nPV_el"]->Fill(nPV,datasets[d], true, LuminosityEl*scaleFactor);
				MSPlot["MS_nLeptons_el"]->Fill(SelectednLeptons,datasets[d], true, LuminosityEl*scaleFactor);
	      MSPlot["MS_nMuons_el"]->Fill(SelectednMu,datasets[d], true, LuminosityEl*scaleFactor);
	      MSPlot["MS_nElectrons_el"]->Fill(SelectednEl,datasets[d], true, LuminosityEl*scaleFactor);
	      MSPlot["MS_MET_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
				MSPlot["MS_METoverST_el"]->Fill(met/ST,datasets[d], true, LuminosityEl*scaleFactor);
				MSPlot["MS_HT_el"]->Fill(HT,datasets[d], true, LuminosityEl*scaleFactor);
				
				if(SelectednEl == 1 && SelectednMu == 0)
				{
				   if(selectedJets[0].Pt() >= LeadingJetPtCut*0.5)
					 {
				        MSPlot["MS_MET_singleEl_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
								MSPlot["MS_JetMultiplicity_singleEl_el"]->Fill(selectedJets.size(),datasets[d], true, LuminosityEl*scaleFactor);
					 }
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
			
			float BtagMCWeight = 1.;
			
					
			//define the 'boxes'
			//'muon channel'
			if(MuTrigged)
			{
			 selecTableMu.Fill(d,0,scaleFactor);
			 if(SelectednMu == 1)
			 {
			  if(SelectednEl == 0)
				{
			    if( (selectedJets.size() == 1 || selectedJets.size() == 2) && selectedForwardJets.size() == 1 )
			    {
					  ST = ST + selectedForwardJets[0].Pt();
					  MSPlot["MS_MET_preboxWqq_noleadingjetcut_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
					  if(selectedJets[0].Pt() >= LeadingJetPtCut)
						{
						  //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
						  float MT = sqrt(2*selectedMuons[0].Pt()*met*(1-cos(selectedMuons[0].DeltaPhi(metvector))));
							MSPlot["MS_MTleptonmet_preboxWqq_mu"]->Fill(MT, datasets[d], true, LuminosityMu*scaleFactor);	
							MSPlot["MS_MET_preboxWqq_noMETMTcut_mu"]->Fill(met, datasets[d], true, LuminosityMu*scaleFactor);	
							if((!doMETcutboxWqq || (doMETcutboxWqq && met>METcutboxWqq)) && (!doMTcut || (doMTcut && MT>MTcut)))
							{		
							 MSPlot["MS_MET_preboxWqq_noBveto_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
							 //cout<<"selectedJets.size() = "<<selectedJets.size()<<", selectedAntiBtaggedJets.size() = "<<selectedAntiBtaggedJets.size()<<", selectedOthertaggedJets.size() = "<<selectedOthertaggedJets.size()<<", selectedBtaggedJets.size() = "<<selectedBtaggedJets.size()<<endl;
							 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0)), not entirely equivalent, now I require to be even taggable (discr value >=0). Difference should be marginal.
							 {
							   if(verbosityBox) cout<<"category: Wqq (mu)"<<endl;
							   //apply b-tag SF on MC; taken from LightStopSearch.cc and https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
	               if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			           {
							      BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								    //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								    MSPlot["MS_BtagWeight_boxWqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								    scaleFactor = 	scaleFactor*BtagMCWeight; 
							   }
								 
			 		       MSPlot["MS_MET_boxWqq_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
							   MSPlot["MS_St_boxWqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
								 //float MTQ = sqrt(2*selectedMuons[0].Pt()*met*(1-cos(selectedMuons[0].DeltaPhi(metvector))));			//that's a 2-body decay system of the W boson... not what we want					 
								 float MTQ = sqrt(2*selectedMuons[0].Pt()*met*(1-cos(selectedMuons[0].DeltaPhi(metvector))) + 2*selectedMuons[0].Pt()*selectedJets[0].Pt()*(1-cos(selectedMuons[0].DeltaPhi(selectedJets[0]))) + 2*met*selectedJets[0].Pt()*(1-cos(metvector.DeltaPhi(selectedJets[0]))));
								 MSPlot["MS_MTQ_boxWqq_mu"]->Fill(MTQ,datasets[d], true, LuminosityMu*scaleFactor);
								 MEzCalculator mezcalc;
								 mezcalc.SetMuon(selectedMuons[0]);
								 mezcalc.SetMET(metvector);
								 mezcalc.SetLeadingJet(selectedJets[0]);
								 double neutrinoz = mezcalc.Calculate(3); //type 3 is my hack; taking the neutrino z such that the reconstructed neutrino four-vector has the largest eta difference with the leading jet
								 double neutrinox = metvector.Px();
								 double neutrinoy = metvector.Py();
								 TLorentzVector recneutrino;
								 recneutrino.SetPxPyPzE(neutrinox,neutrinoy,neutrinoz,sqrt(neutrinox*neutrinox + neutrinoy*neutrinoy + neutrinoz*neutrinoz));
								 float MQ = (recneutrino + selectedMuons[0] + selectedJets[0]).M();
							   MSPlot["MS_MQ_boxWqq_mu"]->Fill(MQ,datasets[d], true, LuminosityMu*scaleFactor);
								 
								 if(selectedMuonsCharge[0] == 1)
								 {
								    if(verbosityBox) cout<<"          W+ qq (mu)"<<endl;
								    MSPlot["MS_St_boxWplusqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
										MSPlot["MS_MTQ_boxWplusqq_mu"]->Fill(MTQ,datasets[d], true, LuminosityMu*scaleFactor);
										MSPlot["MS_MQ_boxWplusqq_mu"]->Fill(MQ,datasets[d], true, LuminosityMu*scaleFactor);
								 }
								 else if(selectedMuonsCharge[0] == -1)
								 {
								    if(verbosityBox) cout<<"          W- qq (mu)"<<endl;
								    MSPlot["MS_St_boxWminusqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
										MSPlot["MS_MTQ_boxWminusqq_mu"]->Fill(MTQ,datasets[d], true, LuminosityMu*scaleFactor);
										MSPlot["MS_MQ_boxWminusqq_mu"]->Fill(MQ,datasets[d], true, LuminosityMu*scaleFactor);
								 }
								 else cout<<"WARNING: muon charge "<<selectedMuonsCharge[1]<<" does not make sense!"<<endl;	
							 }	 								 
							}
						}				
				  }
				  else if((doAntiBtagging && selectedAntiBtaggedJets.size()>=2 && selectedJets.size() >= 3) || (!doAntiBtagging && selectedJets.size() >= 3)) //used to be just selectedJets.size() >= 3
			    {
					  if((doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() >= LeadingJetPtCut) && (selectedAntiBtaggedJets[1].Pt() >= SubLeadingJetPtCut)) || (!doAntiBtagging && (selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))) //used to be just if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
						{
						  //if(!useBtagInfoFromFile)
							//{
							//   bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
							//	 bTool_T->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
							//}
														
						  //b jet(s) from Higgs decay... requiring at least 1; not yet a mass constraint...
							//int nBtags = nBjetsPresent(selectedJets_bTagCSV,btagWP);
						  int nBtags = selectedBtaggedJets.size();							
							
							if(!doBtagging || (doBtagging && nBtags>=1))
							{		
							   //if(!doAntiBtagging || (doAntiBtagging && nNonBjetsPresent(selectedJets_bTagCSV,antibtagWP)>=2))   //TO BE FIXED: >=2 non-b jets -> euhm I just changed it now from ==2 to >=2, is it fixed now? Maybe when doing that, the number of selected events gets larger by a considerable amount, also for the background :S
							   //{		
								    if(verbosityBox) cout<<"category: WHqq (mu)"<<endl;
							      if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							           BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								         //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
											   MSPlot["MS_BtagWeight_boxWHqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								         scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }
										MSPlot["MS_St_boxWHqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);	
										
										//cout<<"selectedJets.size() = "<<selectedJets.size()<<", selectedAntiBtaggedJets.size() = "<<selectedAntiBtaggedJets.size()<<" (selectedAntiBtaggedJets_indices.size() = "<<selectedAntiBtaggedJets_indices.size()<<"), selectedOthertaggedJets.size() = "<<selectedOthertaggedJets.size()<<", selectedBtaggedJets.size() = "<<selectedBtaggedJets.size()<<endl;

										vector<pair<int,float> > selectedAvailableJetsindices_sortedbTagCSV;	//available: everything except the two highest-Pt anti-tagged jets, because I'm going to assume they come from the direct decay of the VLQ							    
										for(int iJet = 0; iJet < selectedJets.size(); iJet++)
			              {
										   if(selectedAntiBtaggedJets_indices[0] != iJet && selectedAntiBtaggedJets_indices[1] != iJet) //meaning, if not one of the two leading anti-btagged jets
											 {
											    pair<int,float> jetindex_btagCSV = make_pair(iJet,selectedJets_bTagCSV[iJet]);	
											    selectedAvailableJetsindices_sortedbTagCSV.push_back(jetindex_btagCSV);
											 }										
										}
										std::sort(selectedAvailableJetsindices_sortedbTagCSV.begin(),selectedAvailableJetsindices_sortedbTagCSV.end(),pair_decrease); //after this selectedAvailableJetsindices_sortedbTagCSV is sorted by decreasing btag value										
											
										int indexClosestJettoBjet;
										float minDeltaR = 9999.;
										float mHcandidate; //'reconstructed' Higgs mass
										float deltamH = 30;
										float mHiggs = 125.; //as generated in the signal samples
										TLorentzVector Hcandidate;
										if(selectedJets.size()>=4)
										{
										 for(int iJet = 1; iJet < selectedAvailableJetsindices_sortedbTagCSV.size(); iJet++) //note the loop starts from 1, the highest-btag jet is not considered
										 {
											  if(selectedJets[selectedAvailableJetsindices_sortedbTagCSV[iJet].first].DeltaR(selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first]) < minDeltaR)
												{
														if(nBtags == 1)
														{														    
														    indexClosestJettoBjet = selectedAvailableJetsindices_sortedbTagCSV[iJet].first;
														    minDeltaR = selectedJets[indexClosestJettoBjet].DeltaR(selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first]);
														}
														else if(nBtags > 1 && selectedAvailableJetsindices_sortedbTagCSV[iJet].second > btagWP)
														{														    
																indexClosestJettoBjet = selectedAvailableJetsindices_sortedbTagCSV[iJet].first;
														    minDeltaR = selectedJets[indexClosestJettoBjet].DeltaR(selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first]);
														}
											 }
										 }
										 
										 Hcandidate = selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first] + selectedJets[indexClosestJettoBjet];
										 mHcandidate = Hcandidate.M();
										 MSPlot["MS_MHcandidate_boxWHqq_mu"]->Fill(mHcandidate,datasets[d], true, LuminosityMu*scaleFactor);					 
										}
										else if(selectedJets.size()==3)
										{
										 //cout<<"selectedAvailableJetsindices_sortedbTagCSV[0].first = "<<selectedAvailableJetsindices_sortedbTagCSV[0].first<<endl;
										 Hcandidate = selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first]; //just the b-tagged jet; sort of assuming it's a fat jet... (but if it's realistic to assume a fat Higgs jet would be b-tagged...?)
										}																		
											
										float mHq;
										int myindex; //index (of the vector of anti-b-tagged jets) of the leading or subleading anti-b-tagged jet, furthest from the reconstructed Higgs direction...
									  if(selectedAntiBtaggedJets[0].DeltaR(Hcandidate) > selectedAntiBtaggedJets[1].DeltaR(Hcandidate)) myindex = 0;
										else myindex = 1;
																																				
								    if((selectedJets.size() >= 4 && fabs(mHcandidate - mHiggs) < deltamH) || selectedJets.size() == 3)
										{        
											 float mHq = (Hcandidate + selectedAntiBtaggedJets[myindex]).M();
											 MSPlot["MS_MHq_boxWHqq_mu"]->Fill(mHq,datasets[d], true, LuminosityMu*scaleFactor);
									  }
										
																				  // //trying mass of jet systems
// 										  for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
// 			                {
// 										     for(unsigned int jJet = iJet+1; jJet < selectedJets.size(); jJet++)
// 			                   {
// 											      float dijetmass = (selectedJets[iJet] + selectedJets[jJet]).M();
// 										        MSPlot["MS_Dijetmass_boxWHqq_mu"]->Fill(dijetmass,datasets[d], true, LuminosityMu*scaleFactor);
// 													  if(selectedJets_bTagCSV[iJet]>btagWP && selectedJets_bTagCSV[jJet]>btagWP) MSPlot["MS_Dibjetmass_boxWHqq_mu"]->Fill(dijetmass,datasets[d], true, LuminosityMu*scaleFactor);
// 											   }
// 										  }										
// 										  if(selectedJets.size()>=4)
// 										  {
// 										     float jet34mass = (selectedJets[2] + selectedJets[3]).M();
// 	                       MSPlot["MS_Dijet34mass_boxWHqq_mu"]->Fill(jet34mass,datasets[d], true, LuminosityMu*scaleFactor);
// 										  }
									
										
				
								 //}
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
				  selecTableMu.Fill(d,1,scaleFactor);
				  if( selectedJets.size() >= 2 )
			    {
					  selecTableMu.Fill(d,2,scaleFactor);
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
						   selecTableMu.Fill(d,3,scaleFactor);
						   //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);						
							 if(!doMETcutboxWWqq || (doMETcutboxWWqq && met>METcutboxWWqq))
						   {				
					      if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0)
						    {
							    selecTableMu.Fill(d,4,scaleFactor);
							    if(verbosityBox) cout<<"category: WWqq (mu)"<<endl; //(no need to check the invariant mass of the leptons, they should not come from a Z...)
							    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							         BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								       //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
											 MSPlot["MS_BtagWeight_boxWWqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								       scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }
					        MSPlot["MS_St_boxWWqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);	
									
							/*		///////////Try to construct some mass variable (doesn't seem to improve much w.r.t ST on first sight)
									//first split met vector in two parts which are equal in size, and are orthogonal
									TLorentzVector metvector1, metvector2;
									metvector1.SetPxPyPzE((metvector.Px()+metvector.Py())/2,(metvector.Py()-metvector.Px())/2,0,sqrt((pow(metvector.Px(),2)+pow(metvector.Py(),2))/2));
									metvector2.SetPxPyPzE((metvector.Px()-metvector.Py())/2,(metvector.Py()+metvector.Px())/2,0,sqrt((pow(metvector.Px(),2)+pow(metvector.Py(),2))/2));
									
									TLorentzVector tmpneutrinosQ1[4],tmpneutrinosQ2[4];
									//hypothesis A: assume metvector1 is associated to muon
									MEzCalculator mezcalc1;
								  mezcalc1.SetMuon(selectedMuons[0]);
								  mezcalc1.SetMET(metvector1);
									mezcalc1.Calculate(4);
								  tmpneutrinosQ1[0] = mezcalc1.gettmpneutrino1();
									tmpneutrinosQ1[1] = mezcalc1.gettmpneutrino2();
									//cout<<" tmpneutrinosQ1[0].Eta() = "<<tmpneutrinosQ1[0].Eta()<<endl;
									//cout<<" tmpneutrinosQ1[1].Eta() = "<<tmpneutrinosQ1[1].Eta()<<endl;
									//metvector2 is associated to electron in that case
									MEzCalculator mezcalc2;
								  mezcalc2.SetEl(selectedElectrons[0]);
								  mezcalc2.SetMET(metvector2);
									mezcalc2.Calculate(4);
								  tmpneutrinosQ2[0] = mezcalc2.gettmpneutrino1();
									tmpneutrinosQ2[1] = mezcalc2.gettmpneutrino2();
									//cout<<" tmpneutrinosQ2[0].Eta() = "<<tmpneutrinosQ2[0].Eta()<<endl;
									//cout<<" tmpneutrinosQ2[1].Eta() = "<<tmpneutrinosQ2[1].Eta()<<endl;
									
									//hypothesis B: assume metvector1 is associated to electron
									MEzCalculator mezcalc3;
								  mezcalc3.SetEl(selectedElectrons[0]);
								  mezcalc3.SetMET(metvector1);
									mezcalc3.Calculate(4);
								  tmpneutrinosQ1[2] = mezcalc3.gettmpneutrino1();
									tmpneutrinosQ1[3] = mezcalc3.gettmpneutrino2();
									//cout<<" tmpneutrinosQ1[2].Eta() = "<<tmpneutrinosQ1[2].Eta()<<endl;
									//cout<<" tmpneutrinosQ1[3].Eta() = "<<tmpneutrinosQ1[3].Eta()<<endl;
									//metvector2 is associated to muon in that case
									MEzCalculator mezcalc4;
								  mezcalc4.SetMuon(selectedMuons[0]);
								  mezcalc4.SetMET(metvector2);
									mezcalc4.Calculate(4);
								  tmpneutrinosQ2[2] = mezcalc4.gettmpneutrino1();
									tmpneutrinosQ2[3] = mezcalc4.gettmpneutrino2();
									//cout<<" tmpneutrinosQ2[2].Eta() = "<<tmpneutrinosQ2[2].Eta()<<endl;
									//cout<<" tmpneutrinosQ2[3].Eta() = "<<tmpneutrinosQ2[3].Eta()<<endl;
																		
									float maxaverageDeltaR_hypoA = -9999, maxaverageDeltaR_hypoB = -9999;
									float MQaverage_hypoA = -9999, MQaverage_hypoB = -9999;
									float MQdiff_hypoA = -9999, MQdiff_hypoB = -9999;
									for(unsigned int inQ1=0;inQ1<2;inQ1++)
									{
									   for(unsigned int inQ2=0;inQ2<2;inQ2++)
									   {
										    float MQ1_option1 = (tmpneutrinosQ1[inQ1] + selectedMuons[0] + selectedJets[0]).M();
										    float MQ2_option1 = (tmpneutrinosQ2[inQ2] + selectedElectrons[0] + selectedJets[1]).M();
												float MQ1_option2 = (tmpneutrinosQ1[inQ1] + selectedMuons[0] + selectedJets[1]).M();
										    float MQ2_option2 = (tmpneutrinosQ2[inQ2] + selectedElectrons[0] + selectedJets[0]).M();
												
												float averageDeltaR_option1 = (tmpneutrinosQ1[inQ1].DeltaR(selectedJets[0]) + tmpneutrinosQ2[inQ2].DeltaR(selectedJets[1]))/2;
												float averageDeltaR_option2 = (tmpneutrinosQ1[inQ1].DeltaR(selectedJets[1]) + tmpneutrinosQ2[inQ2].DeltaR(selectedJets[0]))/2;
												
												if(averageDeltaR_option1 >= averageDeltaR_option2)
												{
												   if(averageDeltaR_option1 > maxaverageDeltaR_hypoA)
													 {
													    maxaverageDeltaR_hypoA = averageDeltaR_option1;
															MQaverage_hypoA = (MQ1_option1 + MQ2_option1) / 2;
															MQdiff_hypoA = fabs(MQ1_option1 - MQ2_option1);
													 }
												}
												else
												{
												   if(averageDeltaR_option2 > maxaverageDeltaR_hypoA)
													 {
													    maxaverageDeltaR_hypoA = averageDeltaR_option2;
															MQaverage_hypoA = (MQ1_option2 + MQ2_option2) / 2;
															MQdiff_hypoA = fabs(MQ1_option2 - MQ2_option2);
													 }												
												}
										 }
									}
									for(unsigned int inQ1=2;inQ1<4;inQ1++)
									{
									   for(unsigned int inQ2=2;inQ2<4;inQ2++)
									   {
										    float MQ1_option1 = (tmpneutrinosQ1[inQ1] + selectedElectrons[0] + selectedJets[0]).M();
										    float MQ2_option1 = (tmpneutrinosQ2[inQ2] + selectedMuons[0] + selectedJets[1]).M();
												float MQ1_option2 = (tmpneutrinosQ1[inQ1] + selectedElectrons[0] + selectedJets[1]).M();
										    float MQ2_option2 = (tmpneutrinosQ2[inQ2] + selectedMuons[0] + selectedJets[0]).M();
												
												float averageDeltaR_option1 = (tmpneutrinosQ1[inQ1].DeltaR(selectedJets[0]) + tmpneutrinosQ2[inQ2].DeltaR(selectedJets[1]))/2;
												float averageDeltaR_option2 = (tmpneutrinosQ1[inQ1].DeltaR(selectedJets[1]) + tmpneutrinosQ2[inQ2].DeltaR(selectedJets[0]))/2;
												
												if(averageDeltaR_option1 >= averageDeltaR_option2)
												{
												   if(averageDeltaR_option1 > maxaverageDeltaR_hypoB)
													 {
													    maxaverageDeltaR_hypoB = averageDeltaR_option1;
															MQaverage_hypoB = (MQ1_option1 + MQ2_option1) / 2;
															MQdiff_hypoB = fabs(MQ1_option1 - MQ2_option1);
													}
												}
												else
												{
												   if(averageDeltaR_option2 > maxaverageDeltaR_hypoB)
													 {
													    maxaverageDeltaR_hypoB = averageDeltaR_option2;
															MQaverage_hypoB = (MQ1_option2 + MQ2_option2) / 2;
															MQdiff_hypoB = fabs(MQ1_option2 - MQ2_option2);
													 }												
												}
										 }
									}	
									
									//decide between hypothesis A and B via mass constraint...												
									float MQaverage;
									if(MQdiff_hypoA <= MQdiff_hypoB) MQaverage = MQaverage_hypoA;
									else MQaverage = MQaverage_hypoB;
									
									MSPlot["MS_AveragePseudoMQ_boxWWqq_mu"]->Fill(MQaverage, datasets[d], true, LuminosityMu*scaleFactor);
									*/
						    }
							 }
					  }
					}
				}
				else if(SelectednEl == 2)
				{
				  if( selectedJets.size() >= 2 )
			    {
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
						  //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
					    
							if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0)
						  {
							   if(verbosityBox) cout<<"category: WZqq (mu)"<<endl;
							   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			           {
							         BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								       //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
											 MSPlot["MS_BtagWeight_boxWZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								       scaleFactor = 	scaleFactor*BtagMCWeight; 
							   }
						     MSPlot["MS_nEvts_boxWZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						     MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						  }
					  }
					}
				}
				/*else if(SelectednEl == 3)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
					   if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
						 {
								if(verbosityBox) cout<<"category: ZZqq (mu)"<<endl;
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
							    //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
							    
									if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						      {
									   if(verbosityBox) cout<<"category: Zqq (mu)"<<endl;
									   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			               {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							       }
									   ST = ST + selectedForwardJets[0].Pt();
							       MSPlot["MS_St_boxZqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
									   //calculate mass of (two-lepton and leading jet)-system
									   float massZq = (selectedMuons[0] + selectedMuons[1] + selectedJets[0]).M();
									   MSPlot["MS_MQ_boxZqq_mu"]->Fill(massZq,datasets[d], true, LuminosityMu*scaleFactor);
								  }
							 }
				     }
				     else if((doAntiBtagging && selectedAntiBtaggedJets.size()>=2 && selectedJets.size() >= 3) || (!doAntiBtagging && selectedJets.size() >= 3)) //used to be just selectedJets.size() >= 3
			       {
						   if((doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() >= LeadingJetPtCut) && (selectedAntiBtaggedJets[1].Pt() >= SubLeadingJetPtCut)) || (!doAntiBtagging && (selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))) //used to be just if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				       {
							    //if(!useBtagInfoFromFile)
							    //{
							    //   bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
								  //   bTool_T->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
							    //}
									int nBtags = selectedBtaggedJets.size();	
							    if(!doBtagging || (doBtagging && nBtags>=1)) //used to be nBjetsPresent(selectedJets_bTagCSV,btagWP)>=1
							    {
							      //  if(!doAntiBtagging || (doAntiBtagging && nNonBjetsPresent(selectedJets_bTagCSV,antibtagWP)==2)) //used to be nNonBjetsPresent(selectedJets_bTagCSV,antibtagWP)==2  //TO BE FIXED: >=2 non-b jets -> euhm I just changed it now from ==2 to >=2, is it fixed now?
							      //  {
											   if(verbosityBox) cout<<"category: ZHqq (mu)"<<endl;
											   if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                   {
							              BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								            //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
														MSPlot["MS_BtagWeight_boxZHqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								            scaleFactor = 	scaleFactor*BtagMCWeight; 
							           }
							           MSPlot["MS_St_boxZHqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
												 
												 float massZq;
												 TLorentzVector Zcandidate = (selectedMuons[0] + selectedMuons[1]);
												 if(Zcandidate.DeltaR(selectedAntiBtaggedJets[0]) > Zcandidate.DeltaR(selectedAntiBtaggedJets[1]))
												    massZq = (selectedMuons[0] + selectedMuons[1] + selectedAntiBtaggedJets[0]).M();
												 else
												    massZq = (selectedMuons[0] + selectedMuons[1] + selectedAntiBtaggedJets[1]).M();												 
												 MSPlot["MS_MQ_boxZHqq_mu"]->Fill(massZq,datasets[d], true, LuminosityMu*scaleFactor);												 
												 
										//	}
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
					  selecTableMu.Fill(d,1,scaleFactor);
					  if( selectedJets.size() >= 2 )
			      {
						   selecTableMu.Fill(d,2,scaleFactor);
					     if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				       {
							    selecTableMu.Fill(d,3,scaleFactor);
							    //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
                  if(!doMETcutboxWWqq || (doMETcutboxWWqq && met>METcutboxWWqq))
						      {
									 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						       {
									   selecTableMu.Fill(d,4,scaleFactor);
									   if(verbosityBox) cout<<"category: WWqq (mu)"<<endl;
									   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			               {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWWqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							       }
									   //cout<<"ST = "<<ST<<", met = "<<met<<", selectedJets[0].Pt() = "<<selectedJets[0].Pt()<<", selectedJets.size() = "<<selectedJets.size()<<endl;
						         MSPlot["MS_St_boxWWqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
										 
							/*	  ///////////Try to construct some mass variable
									//first split met vector in two parts which are equal in size, and are orthogonal
									TLorentzVector metvector1, metvector2;
									metvector1.SetPxPyPzE((metvector.Px()+metvector.Py())/2,(metvector.Py()-metvector.Px())/2,0,sqrt((pow(metvector.Px(),2)+pow(metvector.Py(),2))/2));
									metvector2.SetPxPyPzE((metvector.Px()-metvector.Py())/2,(metvector.Py()+metvector.Px())/2,0,sqrt((pow(metvector.Px(),2)+pow(metvector.Py(),2))/2));
									
									TLorentzVector tmpneutrinosQ1[4],tmpneutrinosQ2[4];
									//hypothesis A: assume metvector1 is associated to muon1
									MEzCalculator mezcalc1;
								  mezcalc1.SetMuon(selectedMuons[0]);
								  mezcalc1.SetMET(metvector1);
									mezcalc1.Calculate(4);
								  tmpneutrinosQ1[0] = mezcalc1.gettmpneutrino1();
									tmpneutrinosQ1[1] = mezcalc1.gettmpneutrino2();
									//metvector2 is associated to muon2 in that case
									MEzCalculator mezcalc2;
								  mezcalc2.SetMuon(selectedMuons[1]);
								  mezcalc2.SetMET(metvector2);
									mezcalc2.Calculate(4);
								  tmpneutrinosQ2[0] = mezcalc2.gettmpneutrino1();
									tmpneutrinosQ2[1] = mezcalc2.gettmpneutrino2();
									
									//hypothesis B: assume metvector1 is associated to muon2
									MEzCalculator mezcalc3;
								  mezcalc3.SetMuon(selectedMuons[1]);
								  mezcalc3.SetMET(metvector1);
									mezcalc3.Calculate(4);
								  tmpneutrinosQ1[2] = mezcalc3.gettmpneutrino1();
									tmpneutrinosQ1[3] = mezcalc3.gettmpneutrino2();
									//metvector2 is associated to muon in that case
									MEzCalculator mezcalc4;
								  mezcalc4.SetMuon(selectedMuons[0]);
								  mezcalc4.SetMET(metvector2);
									mezcalc4.Calculate(4);
								  tmpneutrinosQ2[2] = mezcalc4.gettmpneutrino1();
									tmpneutrinosQ2[3] = mezcalc4.gettmpneutrino2();
																		
									float maxaverageDeltaR_hypoA = -9999, maxaverageDeltaR_hypoB = -9999;
									float MQaverage_hypoA = -9999, MQaverage_hypoB = -9999;
									float MQdiff_hypoA = -9999, MQdiff_hypoB = -9999;
									for(unsigned int inQ1=0;inQ1<2;inQ1++)
									{
									   for(unsigned int inQ2=0;inQ2<2;inQ2++)
									   {
										    float MQ1_option1 = (tmpneutrinosQ1[inQ1] + selectedMuons[0] + selectedJets[0]).M();
										    float MQ2_option1 = (tmpneutrinosQ2[inQ2] + selectedMuons[1] + selectedJets[1]).M();
												float MQ1_option2 = (tmpneutrinosQ1[inQ1] + selectedMuons[0] + selectedJets[1]).M();
										    float MQ2_option2 = (tmpneutrinosQ2[inQ2] + selectedMuons[1] + selectedJets[0]).M();
												
												float averageDeltaR_option1 = (tmpneutrinosQ1[inQ1].DeltaR(selectedJets[0]) + tmpneutrinosQ2[inQ2].DeltaR(selectedJets[1]))/2;
												float averageDeltaR_option2 = (tmpneutrinosQ1[inQ1].DeltaR(selectedJets[1]) + tmpneutrinosQ2[inQ2].DeltaR(selectedJets[0]))/2;
												
												if(averageDeltaR_option1 >= averageDeltaR_option2)
												{
												   if(averageDeltaR_option1 > maxaverageDeltaR_hypoA)
													 {
													    maxaverageDeltaR_hypoA = averageDeltaR_option1;
															MQaverage_hypoA = (MQ1_option1 + MQ2_option1) / 2;
															MQdiff_hypoA = fabs(MQ1_option1 - MQ2_option1);
													 }
												}
												else
												{
												   if(averageDeltaR_option2 > maxaverageDeltaR_hypoA)
													 {
													    maxaverageDeltaR_hypoA = averageDeltaR_option2;
															MQaverage_hypoA = (MQ1_option2 + MQ2_option2) / 2;
															MQdiff_hypoA = fabs(MQ1_option2 - MQ2_option2);
													 }												
												}
										 }
									}
									for(unsigned int inQ1=2;inQ1<4;inQ1++)
									{
									   for(unsigned int inQ2=2;inQ2<4;inQ2++)
									   {
										    float MQ1_option1 = (tmpneutrinosQ1[inQ1] + selectedMuons[1] + selectedJets[0]).M();
										    float MQ2_option1 = (tmpneutrinosQ2[inQ2] + selectedMuons[0] + selectedJets[1]).M();
												float MQ1_option2 = (tmpneutrinosQ1[inQ1] + selectedMuons[1] + selectedJets[1]).M();
										    float MQ2_option2 = (tmpneutrinosQ2[inQ2] + selectedMuons[0] + selectedJets[0]).M();
												
												float averageDeltaR_option1 = (tmpneutrinosQ1[inQ1].DeltaR(selectedJets[0]) + tmpneutrinosQ2[inQ2].DeltaR(selectedJets[1]))/2;
												float averageDeltaR_option2 = (tmpneutrinosQ1[inQ1].DeltaR(selectedJets[1]) + tmpneutrinosQ2[inQ2].DeltaR(selectedJets[0]))/2;
												
												if(averageDeltaR_option1 >= averageDeltaR_option2)
												{
												   if(averageDeltaR_option1 > maxaverageDeltaR_hypoB)
													 {
													    maxaverageDeltaR_hypoB = averageDeltaR_option1;
															MQaverage_hypoB = (MQ1_option1 + MQ2_option1) / 2;
															MQdiff_hypoB = fabs(MQ1_option1 - MQ2_option1);
													}
												}
												else
												{
												   if(averageDeltaR_option2 > maxaverageDeltaR_hypoB)
													 {
													    maxaverageDeltaR_hypoB = averageDeltaR_option2;
															MQaverage_hypoB = (MQ1_option2 + MQ2_option2) / 2;
															MQdiff_hypoB = fabs(MQ1_option2 - MQ2_option2);
													 }												
												}
										 }
									}	
									
									//decide between hypothesis A and B via mass constraint...												
									float MQaverage;
									if(MQdiff_hypoA <= MQdiff_hypoB) MQaverage = MQaverage_hypoA;
									else MQaverage = MQaverage_hypoB;
									
									MSPlot["MS_AveragePseudoMQ_boxWWqq_mu"]->Fill(MQaverage, datasets[d], true, LuminosityMu*scaleFactor);		
									*/	
								   }
									}
						   }		
					  }
					}
				}
				else if(SelectednEl == 1)
				{
				  if( selectedJets.size() >= 2 )
			    {
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
						   //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
					     
							 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						   {
							    if(verbosityBox) cout<<"category: WZqq (mu)"<<endl;
							    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }
					        MSPlot["MS_nEvts_boxWZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						      MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						   }
					  }
					}					
				}
				else if(SelectednEl == 2)
				{
				  if( selectedJets.size() >= 2 )
			    {
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
						   //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
					     
							 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						   {
							    if(verbosityBox) cout<<"category: ZZqq (mu)"<<endl;
							    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }
					        MSPlot["MS_nEvts_boxZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						      MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						   }
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
				  if( selectedJets.size() >= 2 )
			    {
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
						   //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
					     
							 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))  //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						   {
							    if(verbosityBox) cout<<"category: WZqq (mu)"<<endl;
							    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }
					        MSPlot["MS_nEvts_boxWZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						      MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						   }
					  }
					}
				}
				/*else if(SelectednEl == 1)
				{
				  if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				  {
					   if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
						 {
								if(verbosityBox) cout<<"category: ZZqq (mu)"<<endl;
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
				  if( selectedJets.size() >= 2 )
			    {
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
						  //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
					    
							if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						  {
							  if(verbosityBox) cout<<"category: ZZqq (mu)"<<endl;
							  if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			          {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							  }
					      MSPlot["MS_nEvts_boxZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						    MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						  }
					  }
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
					   ST = ST + selectedForwardJets[0].Pt();
					   MSPlot["MS_MET_preboxWqq_noleadingjetcut_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
					   if(selectedJets[0].Pt() >= LeadingJetPtCut)
				     { 
						    //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
								float MT = sqrt(2*selectedElectrons[0].Pt()*met*(1-cos(selectedElectrons[0].DeltaPhi(metvector))));
								MSPlot["MS_MTleptonmet_preboxWqq_el"]->Fill(MT, datasets[d], true, LuminosityEl*scaleFactor);
						    MSPlot["MS_MET_preboxWqq_noMETMTcut_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
								if((!doMETcutboxWqq || (doMETcutboxWqq && met>METcutboxWqq)) && (!doMTcut || (doMTcut && MT>MTcut)))
							  {	
								 MSPlot["MS_MET_preboxWqq_noBveto_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
								 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							   {	
								    if(verbosityBox) cout<<"category: Wqq (el)"<<endl;
								    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }
								 
								   MSPlot["MS_MET_boxWqq_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
						       MSPlot["MS_St_boxWqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
									 //float MTQ = sqrt(2*selectedElectrons[0].Pt()*met*(1-cos(selectedElectrons[0].DeltaPhi(metvector))));
									 float MTQ = sqrt(2*selectedElectrons[0].Pt()*met*(1-cos(selectedElectrons[0].DeltaPhi(metvector))) + 2*selectedElectrons[0].Pt()*selectedJets[0].Pt()*(1-cos(selectedElectrons[0].DeltaPhi(selectedJets[0]))) + 2*met*selectedJets[0].Pt()*(1-cos(metvector.DeltaPhi(selectedJets[0]))));								 
								   MSPlot["MS_MTQ_boxWqq_el"]->Fill(MTQ,datasets[d], true, LuminosityEl*scaleFactor);
									 MEzCalculator mezcalc;
								   mezcalc.SetEl(selectedElectrons[0]);
								   mezcalc.SetMET(metvector);
								   mezcalc.SetLeadingJet(selectedJets[0]);
								   double neutrinoz = mezcalc.Calculate(3); //type 3 is my hack; taking the neutrino z such that the reconstructed neutrino four-vector has the largest eta difference with the leading jet
								   double neutrinox = metvector.Px();
								   double neutrinoy = metvector.Py();
								   TLorentzVector recneutrino;
								   recneutrino.SetPxPyPzE(neutrinox,neutrinoy,neutrinoz,sqrt(neutrinox*neutrinox + neutrinoy*neutrinoy + neutrinoz*neutrinoz));
								   float MQ = (recneutrino + selectedElectrons[0] + selectedJets[0]).M();
							     MSPlot["MS_MQ_boxWqq_el"]->Fill(MQ,datasets[d], true, LuminosityEl*scaleFactor);
								 
									 if(selectedElectronsCharge[0] == 1)
								   {
									    if(verbosityBox) cout<<"          W+ qq (el)"<<endl;
								      MSPlot["MS_St_boxWplusqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
											MSPlot["MS_MTQ_boxWplusqq_el"]->Fill(MTQ,datasets[d], true, LuminosityEl*scaleFactor);
											MSPlot["MS_MQ_boxWplusqq_el"]->Fill(MQ,datasets[d], true, LuminosityEl*scaleFactor);
								   }
								   else if(selectedElectronsCharge[0] == -1)
								   {
									    if(verbosityBox) cout<<"          W- qq (el)"<<endl;
								      MSPlot["MS_St_boxWminusqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
											MSPlot["MS_MTQ_boxWminusqq_el"]->Fill(MTQ,datasets[d], true, LuminosityEl*scaleFactor);
											MSPlot["MS_MQ_boxWminusqq_el"]->Fill(MQ,datasets[d], true, LuminosityEl*scaleFactor);
								   }
								   else cout<<"WARNING: electron charge "<<selectedElectronsCharge[1]<<" does not make sense!"<<endl;
								 }
								}
						 }
				  }
				  else if((doAntiBtagging && selectedAntiBtaggedJets.size()>=2 && selectedJets.size() >= 3) || (!doAntiBtagging && selectedJets.size() >= 3))
			    {
					   if((doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() >= LeadingJetPtCut) && (selectedAntiBtaggedJets[1].Pt() >= SubLeadingJetPtCut)) || (!doAntiBtagging && (selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut)))
				     {
						    //if(!useBtagInfoFromFile)
							  //{
							  //     bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
								//     bTool_T->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
							  //}
								
								int nBtags = selectedBtaggedJets.size();	
						    if(!doBtagging || (doBtagging && nBtags>=1))
							  {	
									  if(verbosityBox) cout<<"category: WHqq (el)"<<endl;
									  if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							           BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								         //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
											   MSPlot["MS_BtagWeight_boxWHqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								         scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }
										
							      MSPlot["MS_St_boxWHqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
										
										vector<pair<int,float> > selectedAvailableJetsindices_sortedbTagCSV;	//available: everything except the two highest-Pt anti-tagged jets, because I'm going to assume they come from the direct decay of the VLQ							    
										for(int iJet = 0; iJet < selectedJets.size(); iJet++)
			              {
										   if(selectedAntiBtaggedJets_indices[0] != iJet && selectedAntiBtaggedJets_indices[1] != iJet) //meaning, if not one of the two leading anti-btagged jets
											 {
											    pair<int,float> jetindex_btagCSV = make_pair(iJet,selectedJets_bTagCSV[iJet]);	
											    selectedAvailableJetsindices_sortedbTagCSV.push_back(jetindex_btagCSV);
											 }										
										}
										std::sort(selectedAvailableJetsindices_sortedbTagCSV.begin(),selectedAvailableJetsindices_sortedbTagCSV.end(),pair_decrease); //after this selectedAvailableJetsindices_sortedbTagCSV is sorted by decreasing btag value										
											
										int indexClosestJettoBjet;
										float minDeltaR = 9999.;
										float mHcandidate; //'reconstructed' Higgs mass
										float deltamH = 30;
										float mHiggs = 125.; //as generated in the signal samples
										TLorentzVector Hcandidate;
										if(selectedJets.size()>=4)
										{
										 for(int iJet = 1; iJet < selectedAvailableJetsindices_sortedbTagCSV.size(); iJet++) //note the loop starts from 1, the highest-btag jet is not considered
										 {
											  if(selectedJets[selectedAvailableJetsindices_sortedbTagCSV[iJet].first].DeltaR(selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first]) < minDeltaR)
												{
														if(nBtags == 1)
														{														    
														    indexClosestJettoBjet = selectedAvailableJetsindices_sortedbTagCSV[iJet].first;
														    minDeltaR = selectedJets[indexClosestJettoBjet].DeltaR(selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first]);
														}
														else if(nBtags > 1 && selectedAvailableJetsindices_sortedbTagCSV[iJet].second > btagWP)
														{														    
																indexClosestJettoBjet = selectedAvailableJetsindices_sortedbTagCSV[iJet].first;
														    minDeltaR = selectedJets[indexClosestJettoBjet].DeltaR(selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first]);
														}
											 }
										 }
										 
										 Hcandidate = selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first] + selectedJets[indexClosestJettoBjet];
										 mHcandidate = Hcandidate.M();
										 MSPlot["MS_MHcandidate_boxWHqq_el"]->Fill(mHcandidate,datasets[d], true, LuminosityEl*scaleFactor);									 
										}
										else if(selectedJets.size()==3)
										{
										 //cout<<"selectedAvailableJetsindices_sortedbTagCSV[0].first = "<<selectedAvailableJetsindices_sortedbTagCSV[0].first<<endl;
										 Hcandidate = selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first]; //just the b-tagged jet; sort of assuming it's a fat jet... (but if it's realistic to assume a fat Higgs jet would be b-tagged...?)
										}								
																				
										float mHq;
										int myindex; //index (of the vector of anti-b-tagged jets) of the leading or subleading anti-b-tagged jet, furthest from the reconstructed Higgs direction...
									  if(selectedAntiBtaggedJets[0].DeltaR(Hcandidate) > selectedAntiBtaggedJets[1].DeltaR(Hcandidate)) myindex = 0;
										else myindex = 1;
																																				
								    if((selectedJets.size() >= 4 && fabs(mHcandidate - mHiggs) < deltamH) || selectedJets.size() == 3)
										{                       
											 float mHq = (Hcandidate + selectedAntiBtaggedJets[myindex]).M();
											 MSPlot["MS_MHq_boxWHqq_el"]->Fill(mHq,datasets[d], true, LuminosityEl*scaleFactor);
									  }
											
											//trying mass of jet systems
// 										  for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
// 			                {
// 										     for(unsigned int jJet = iJet+1; jJet < selectedJets.size(); jJet++)
// 			                   {
// 											      float dijetmass = (selectedJets[iJet] + selectedJets[jJet]).M();
// 										        MSPlot["MS_Dijetmass_boxWHqq_el"]->Fill(dijetmass,datasets[d], true, LuminosityEl*scaleFactor);
// 													  if(selectedJets_bTagCSV[iJet]>btagWP && selectedJets_bTagCSV[jJet]>btagWP) MSPlot["MS_Dibjetmass_boxWHqq_el"]->Fill(dijetmass,datasets[d], true, LuminosityEl*scaleFactor);
// 											   }
// 										  }										
// 										  if(selectedJets.size()>=4)
// 										  {
// 										     float jet34mass = (selectedJets[2] + selectedJets[3]).M();
// 	                       MSPlot["MS_Dijet34mass_boxWHqq_el"]->Fill(jet34mass,datasets[d], true, LuminosityEl*scaleFactor);
// 										  }
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
								  //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
								  
									if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						      {
									  if(verbosityBox) cout<<"category: Zqq (el)"<<endl;
									  if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }
							      ST = ST + selectedForwardJets[0].Pt();
							      MSPlot["MS_St_boxZqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
									  //calculate mass of (two-lepton and leading jet)-system
									  float massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedJets[0]).M();
									  MSPlot["MS_MQ_boxZqq_el"]->Fill(massZq,datasets[d], true, LuminosityMu*scaleFactor);
								  }
								}
				     }
				     else if((doAntiBtagging && selectedAntiBtaggedJets.size()>=2 && selectedJets.size() >= 3) || (!doAntiBtagging && selectedJets.size() >= 3))
			       {
						    if((doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() >= LeadingJetPtCut) && (selectedAntiBtaggedJets[1].Pt() >= SubLeadingJetPtCut)) || (!doAntiBtagging && (selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut)))
				        {			
								   //if(!useBtagInfoFromFile)
							     //{
							     //    bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
								   //    bTool_T->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
							     //}					
								   //b jet(s) from Higgs decay... requiring at least 1; not yet a mass constraint...
						       int nBtags = selectedBtaggedJets.size();
									 if(!doBtagging || (doBtagging && nBtags>=1))
							     {	
							        //if(!doAntiBtagging || (doAntiBtagging && nNonBjetsPresent(selectedJets_bTagCSV,antibtagWP)==2))  //TO BE FIXED: >=2 non-b jets
							        //{	
											   if(verbosityBox) cout<<"category: ZHqq (el)"<<endl;
											   if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                   {
							               BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								             //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
														 MSPlot["MS_BtagWeight_boxZHqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								             scaleFactor = 	scaleFactor*BtagMCWeight; 
							           }
							           //calculate mass of (two-lepton and leading jet)-system (CHECK THE JETS ARE SORTED!!!)
							           //float massZq1 = (selectedElectrons[0] + selectedElectrons[1] + selectedJets[0]).M();
							 
							           MSPlot["MS_St_boxZHqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
												 
												 float massZq;
												 TLorentzVector Zcandidate = (selectedElectrons[0] + selectedElectrons[1]);
												 if(Zcandidate.DeltaR(selectedAntiBtaggedJets[0]) > Zcandidate.DeltaR(selectedAntiBtaggedJets[1]))
												    massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedAntiBtaggedJets[0]).M();
												 else
												    massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedAntiBtaggedJets[1]).M();												 
												 MSPlot["MS_MQ_boxZHqq_el"]->Fill(massZq,datasets[d], true, LuminosityEl*scaleFactor);
												 
											//}
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
					   if( selectedJets.size() >= 2 )
			       {
					      if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				        {
								  //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
						      if(!doMETcutboxWWqq || (doMETcutboxWWqq && met>METcutboxWWqq))
									{
									 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						       {
									   if(verbosityBox) cout<<"category: WWqq (el)"<<endl;
									   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			               {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWWqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							       }
						         MSPlot["MS_St_boxWWqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
							     }
									}		
						    }						 
						 }							
					}
				}
				else if(SelectednEl == 3)
				{
				   if( selectedJets.size() >= 2 )
			     {
				     if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				     {
						    //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
								if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						    {
								   if(verbosityBox) cout<<"category: WZqq (el)"<<endl;
								   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			             {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							     }
					         MSPlot["MS_nEvts_boxWZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
							     MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
							  }	
					   }
					}							
				}
				else if(SelectednEl == 4)
				{
				   if( selectedJets.size() >= 2 )
			     {
				     if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				     {
						    //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
					      
								if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						    {
								   if(verbosityBox) cout<<"category: ZZqq (el)"<<endl;
								   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			             {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							     }
					         MSPlot["MS_nEvts_boxZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
							     MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
							  }
					   }
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
				  if( selectedJets.size() >= 2 )
			    {
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
						   //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
					     if(!doMETcutboxWWqq || (doMETcutboxWWqq && met>METcutboxWWqq))
						   {
							  if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						    {
							    if(verbosityBox) cout<<"category: WWqq (el)"<<endl; //(no need to check the invariant mass of the leptons, they should not come from a Z...)
							    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWWqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }
					        MSPlot["MS_St_boxWWqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
						    }
							 }
					  }
					}	
				}
				else if(SelectednEl == 2)
				{
				  if( selectedJets.size() >= 2 )
			    {
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
						   //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
					     
							 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						   {
							    if(verbosityBox) cout<<"category: WZqq (el)"<<endl;
							    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }
					        MSPlot["MS_nEvts_boxWZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						      MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						   }
					  }
					}
				}
// 				else if(SelectednEl == 3)
// 				{
// 				  if(selectedJets[0].Pt() >= LeadingJetPtCut)
// 				  {
// 					   if(!doAntiBtagging || (doAntiBtagging && nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0))
// 						 {
// 								if(verbosityBox) cout<<"category: ZZqq (el)"<<endl;
// 					      MSPlot["MS_nEvts_boxZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
// 						    MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
// 						 }
// 					}
// 				}
				else
				{
				   cout<<"FOUND >= 5-lepton event!"<<endl;		
				}									
			 }
			 else if(SelectednMu == 2)
			 {			  
				if(SelectednEl == 1)
				{
				  if( selectedJets.size() >= 2 )
			    {
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
						  //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
					    
							if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						  {
							   if(verbosityBox) cout<<"category: WZqq (el)"<<endl;
							   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			           {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							   }
					       MSPlot["MS_nEvts_boxWZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);	
						     MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						  }	
				    }
					}			
				}
				else if(SelectednEl == 2)
				{
				  if( selectedJets.size() >= 2 )
			    {
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
						  //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
					    
							if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						  {
							   if(verbosityBox) cout<<"category: ZZqq (el)"<<endl;
							   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			           {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							   }
					       MSPlot["MS_nEvts_boxZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						     MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						  }
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
// 				if(SelectednEl == 1)
// 				{
// 				  if(selectedJets[0].Pt() >= LeadingJetPtCut)
// 				  {
// 						if(verbosityBox) cout<<"category: ZZqq (el)"<<endl;
// 					  MSPlot["MS_nEvts_boxZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
// 						MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
// 					}
// 				}
// 				else 
// 				{
// 				   cout<<"FOUND >= 5-lepton event!"<<endl;		
// 				}
			 }
		 
			  
			}	
			

		 
	//delete selection;
     //cout << "done processing event" << endl;
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
	
	
	
	if (!useBtagInfoFromFile)
	{
	   bTool_L->WriteMCEfficiencyHistos("PlotsForBtagWeights_CSVL_ttbar.root");
		 bTool_T->WriteMCEfficiencyHistos("PlotsForBtagWeights_CSVM_ttbar.root");
	}
	
	
	
	
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
				temp->showNumberEntries(false);
				if(!foundMu && !foundEl)
	          temp->setDataLumi(Luminosity);//if not run on data...	NARNING: I have to reconsider this; when running on muon-data only for example, the  luminosity in the electon plots still needs to be set like this...
        //temp->Draw(false, name, true, true, true, true, true,5,false, true, true, true, true);// old way: (bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal, bool addRatio, bool mergeVV, bool mergeTTV, bool mergeVVV, bool mergeSameSignWW)
        temp->setErrorBandFile("VLQSearch_ErrorBand/ErrorBandFile.root");
				temp->Draw(name,0,false,false,false,5); //label, RatioType, addRatioErrorBand, addErrorBand, ErrorBandAroundTotalInput, scaleNPSignal
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
    selecTableMu.TableCalculator(true, true, true, true, true, true, true,false);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
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

/*unsigned int nBjetsPresent(vector < float > btagvalues, float btagWP)
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
	    if((btagvalues[iJet] < btagWP) && (btagvalues[iJet] >= 0))
			{
					 nNonBjets++; //because the current jet would not be b tagged
			}
	 }
	 return nNonBjets;
};
*/
