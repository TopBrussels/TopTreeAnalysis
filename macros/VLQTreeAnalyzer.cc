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
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysis/VLQSearch/interface/VLQTree.h"
#include "TopTreeAnalysis/VLQSearch/interface/TemplateTools.h"

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

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

bool hasEnding (std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

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
  if( ! (systematic == "Nominal"  || systematic == "JESMinus" || systematic == "JESPlus" || systematic == "METUnClusteredEnergyMinus" || systematic == "METUnClusteredEnergyPlus" || systematic == "JERMinus" || systematic == "JERPlus" || systematic == "bTagMinus" || systematic == "bTagPlus" || systematic == "misTagMinus" || systematic == "misTagPlus" || systematic == "PUMinus" || systematic == "PUPlus" || systematic == "MuonSFMinus" || systematic == "MuonSFPlus" || systematic == "ElectronSFMinus" || systematic == "ElectronSFPlus" ))
  {
    cout << "Unknown systematic!!!" << endl;
    cout << "Possible options are: Nominal JESMinus JESPlus METUnClusteredEnergyMinus METUnClusteredEnergyPlus JERMinus JERPlus bTagMinus bTagPlus misTagMinus misTagPlus PUMinus PUPlus" << endl;
    exit(-1);
  }	  

  string btagger = "CSVM", antibtagger = "CSVL"; //string btagger = "CSVM", antibtagger = "CSVL"
  float btagWP = 9999, antibtagWP = 9999;
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

  cout << "***********************************************" << endl;
  cout << " Beginning of the program for the VLQ search ! " << endl;
  cout << "***********************************************" << endl;

  setTDRStyle();

  string inputpostfixOld = ""; // should be same as postfix in TreeCreator of the trees
	string inputpostfix= inputpostfixOld+"_"+systematic;		

  //string Treespath = "VLQTrees_Summer12_PBS_24June2013"; //v4
	//string Treespath = "VLQTrees_Summer12_PBS_v5_21Aug13"; //v5
	//string Treespath = "VLQTrees_Summer12_PBS_v5_27Jan14"; //v5 with a fix (looser electron veto)
	string Treespath = "VLQTrees_Summer12_PBS_v5_9Mar14"; //v5 with a fix (channel definition)
  Treespath = Treespath + "/"; 		
	bool savePNG = false;
	string outputpostfix = "";
	outputpostfix = outputpostfix+"_"+systematic;
	//string Outputpath = "OutputFiles_VLQAnalyzer_withQCDestimation_28Jan14_signal500overlay_WZjetsFitted_FlatBinning_noModelScaling"; //"OutputFiles_VLQAnalyzer_InvIsoABCD_27Jan14"; // "OutputFiles_VLQAnalyzer_withoutQCDestimation_27Jan14_signal500overlay"; // "OutputFiles_VLQAnalyzer_withoutQCDestimation_16Jan14_signal500overlay_PtTopReweighting_CMStopXS_PtVcuts";
	string Outputpath = "OutputFiles_VLQAnalyzer_withQCDestimation_10Mar14_FlatBinning_NoModelScaling";
	Outputpath = Outputpath + "/";
	mkdir(Outputpath.c_str(),0777);
	//for flat binning
	//string Binningpath = "OutputFiles_VLQAnalyzer_withQCDestimation_28Jan14_signal700overlay_WZjetsFitted_FlatBinning_noModelScaling/";
	string Binningpath = "FlatBinningFiles_VLQAnalyzer_withQCDestimation_10Mar14/";
  mkdir(Binningpath.c_str(),0777);
  
	/////////////////////
  // Configuration
  /////////////////////

	//bool makebinning = false;
  
  //xml file
  string xmlFileName = "";
	xmlFileName = "../config/myVLQConfig.xml";	
  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;   
	
	bool makeflatbinning = false;
	bool useflatbinning = true;
	
	bool doPUreweighting = true;
	bool applyLeptonSF = true;
	bool applyBtagSF = true;
	bool applyPtTopReweighting = true; //MadGraph; see VLQTreeCreator.cc
	
	float LeadingJetPtCut = 200., SubLeadingJetPtCut = 100.; //subleading jet pt cut only for the pair production signal boxes...
	bool doBtagging = true; //Higgs->bbar decay; 'false' maybe not supported anymore due to the complexity of all selections / control regions
	bool doAntiBtagging = true; //all other jets; 'false' maybe not supported anymore due to the complexity of all selections / control regions
	bool useBtagInfoFromFile = true;	
	bool doMETcutboxWqq = true; //would only be in Wqq categories (i.e. single lepton with only 1 high-Pt jet)
	float METcutboxWqq = 60.; //problem: QCD gets mostly rejected with this, but it is also cutting away quite a lot of signal!!
	bool doMTcut = true;//would only be in Wqq categories (i.e. single lepton with only 1 high-Pt jet), instead of MET cut...
	float MTcut = 40.;
	bool doMETcutboxWWqq = true; //would only be in WWqq categories (i.e. double lepton)
	float METcutboxWWqq = 60.;
	bool doMETcutboxWHqq = true; //would only be in WHqq categories (i.e. single lepton with a b-tag)
	float METcutboxWHqq = 60.;	
	bool doMETcutboxWZqq = true; //would only be in WZqq categories (i.e. three lepton (not the semilep box))
	float METcutboxWZqq = 60.;	
	bool doWptcut = true; //in Wqq category, where neutrino is reconstructed
	float Wptcut = 150;
	bool doZptcut = true; //in Zqq category, where neutrino is reconstructed, but also in ZHqq and WZqq (and ZZqq?) categories!
	float Zptcut = 150;
	bool doHptcut = true; //in WHqq category, where higgs is tried to be reconstructed
	float Hptcut = 75; //I saw that 150 was a bit too high for the signal (I think this reflects the fact that the H candidate is not always so well reconstructed...)
	
	bool verbosityBox = false; //to cout in which signal box you end up
	
	//b-tag systematics
	int syst_btag = 0; //-1 for SF up, +1 for SF down -> wait what? it's -1 for down and +1 for up, right??
	if(systematic == "bTagMinus") syst_btag = -1;
	else if(systematic == "bTagPlus") syst_btag = 1;
	int syst_mistag = 0; //-1 for SF up, +1 for SF down -> wait what? it's -1 for down and +1 for up, right??
	if(systematic == "misTagMinus") syst_mistag = -1;
	else if(systematic == "misTagPlus") syst_mistag = 1;
	
	//lepton SF systematics
	string syst_muonSF = "Nominal";
	if(systematic == "MuonSFMinus") syst_muonSF = "Minus";
	else if(systematic == "MuonSFPlus") syst_muonSF = "Plus";
	string syst_electronSF = "Nominal";
	if(systematic == "ElectronSFMinus") syst_electronSF = "Minus";
	else if(systematic == "ElectronSFPlus") syst_electronSF = "Plus";
	
	
	//particular signal model (not only for the overlay mass point, but for all mass points, for consistency. Make sure you don't forget this in the scanning later on!!)
	bool doModelScaling = false; //to reweight according to given Branching Fractions and kappatilde values or not
	float BF_W = 0.5, BF_Z = 0.25;
	float BF_H = 1.0 - BF_W - BF_Z;
	float kappatilde_W = 1.;
	float kappatilde_Z;  
	if(kappatilde_W != 0 && BF_W != 0) kappatilde_Z = kappatilde_W * sqrt( 2 * BF_Z / BF_W); 
	else  kappatilde_Z = sqrt( 2 * BF_Z ); //IS THIS CORRECT??
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
	
	bool doWjetsScaling = true;
	float WjetsScale = 0.84;//0.86 //0.877//result from LF fit in mu channel (fit in el channel maybe not so safe as long as you don't have a QCD estimation: I fit there something like 0.92, a bit higher, consistent with QCD presence pushing the data a bit up?)
  bool doZjetsScaling = true;
	float ZjetsScale = 0.85; //0.86 //result from LF fit in mu channel
  bool doInvIsoScaling = true;
	float InvIsomuScale = 0.704; //0.764;//0.737 //MET_singleMu_Bveto 1.003; //fit of MET_preboxWqq_leadingjet150-200_noForwardJet_noMETMTcut_mu with all other bkgs fixed from other fits (except ttbar, thats just fixed to 1...)
	float InvIsoelScale = 0.423;//0.801;//0.440 //MET_singleEl_Bveto 0.563; //fit of MET_preboxWqq_leadingjet150-200_noForwardJet_noMETMTcut with all other bkgs fixed from other fits (except ttbar, thats just fixed to 1...)


  //MultiSamplePlot options
	bool showEntriesLegend = false;
	bool setCMSPrelim = false;
	int RatioType = 2; 
	bool addRatioErrorBand = false;
	bool addErrorBand = false;
	bool ErrorBandAroundTotalInput = false;
	int scaleNPSignal = 5; //only scaled on canvas, not the template in MSPlot itself!? (TO BE CHECKED!!) (I think I checked...)
	
	
	
  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = anaEnv.Verbose;
  float anaEnvLuminosity = anaEnv.Luminosity;	// in 1/pb 
//  cout << "analysis environment luminosity for rescaling (if not run on data) "<< anaEnvLuminosity << endl;

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
  
  //float Luminosity = anaEnvLuminosity;

  //hardcoded
  float LuminosityMu = 19600.84; //anaEnvLuminosity;
  float LuminosityEl = 19629.525; //anaEnvLuminosity;

  bool foundMu = false;
  bool foundEl = false;
	
	vector<string> inputTrees; //fill with the tree files you want to read!
	for (unsigned int d = 0; d < datasets.size(); d++)
  {				
    cout << "  Dataset " << d << " name : " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;    
    string dataSetName = datasets[d]->Name();
		
		if(dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
      LuminosityMu = datasets[d]->EquivalentLumi();
      foundMu=true;
    }  
    if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) {
      LuminosityEl = datasets[d]->EquivalentLumi();
      foundEl=true;  
    } 
		
		//if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // when only using one channel
    //  Luminosity = datasets[d]->EquivalentLumi();
    //} 
		
		string inputTreeFileName; //should follow convention of TreeFileName in VLQTreeCreator.cc
		if(systematic == "JESMinus" || systematic == "JESPlus"|| systematic == "METUnClusteredEnergyMinus" || systematic == "METUnClusteredEnergyPlus" || systematic == "JERMinus" || systematic == "JERPlus")
		{
				if ( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 || dataSetName.find("InvIso") != string::npos)
				{
					inputTreeFileName = Treespath+"VLQTree_"+dataSetName+inputpostfixOld+".root"; //is actually dummy -> well not really for InvIso
					//cout<<"   Running systematics: data/InvIso will be skipped later on"<<endl;
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
	
	//if(!foundMu && !foundEl && Luminosity != anaEnvLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  //else {
    if(LuminosityMu != anaEnvLuminosity) cout << "Muon PD: changed analysis environment luminosity to "<< LuminosityMu << endl;
    if(LuminosityEl != anaEnvLuminosity) cout << "Electron PD: changed analysis environment luminosity to "<< LuminosityEl << endl;
  //}
	
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
	
	
	

	//////////////////////
	/// for 'flat' binning
	//////////////////////
	//muon channel
	//// Wqq ////
	string xvariable = "MQ";
	int nbinsxvariable_Wqq = 100;
	int nbinsxvariable_Wplusqq = 50;
	int nbinsxvariable_Wminusqq = 50;
	string binningFileName_Wqq_mu = Binningpath+"Binning_VLQSearch_Wqq_mu.root";
	string binningFileName_Wplusqq_mu = Binningpath+"Binning_VLQSearch_Wplusqq_mu.root";
	string binningFileName_Wminusqq_mu = Binningpath+"Binning_VLQSearch_Wminusqq_mu.root";
	TemplateTools VarWqq_mu("Wqq_mu",xvariable,nbinsxvariable_Wqq);
	TemplateTools VarWplusqq_mu("Wplusqq_mu",xvariable,nbinsxvariable_Wplusqq);
  TemplateTools VarWminusqq_mu("Wminusqq_mu",xvariable,nbinsxvariable_Wminusqq);	
	VarWqq_mu.SetDatasets(datasetsMu);
	VarWplusqq_mu.SetDatasets(datasetsMu);
	VarWminusqq_mu.SetDatasets(datasetsMu);
	//// WHqq ////
	xvariable = "MHq";
	int nbinsxvariable_WHqq = 50;
	string binningFileName_WHqq_mu = Binningpath+"Binning_VLQSearch_WHqq_mu.root";
	TemplateTools VarWHqq_mu("WHqq_mu",xvariable,nbinsxvariable_WHqq);
	VarWHqq_mu.SetDatasets(datasetsMu);
	// Zqq ////
	xvariable = "MQ";
	int nbinsxvariable_Zqq = 25;
	string binningFileName_Zqq_mu = Binningpath+"Binning_VLQSearch_Zqq_mu.root";
	TemplateTools VarZqq_mu("Zqq_mu",xvariable,nbinsxvariable_Zqq);
	VarZqq_mu.SetDatasets(datasetsMu);
	//// ZHqq ////
	xvariable = "MQ";
	int nbinsxvariable_ZHqq = 4;
	string binningFileName_ZHqq_mu = Binningpath+"Binning_VLQSearch_ZHqq_mu.root";
	TemplateTools VarZHqq_mu("ZHqq_mu",xvariable,nbinsxvariable_ZHqq);
	VarZHqq_mu.SetDatasets(datasetsMu);
	//// ZZqqWZqq semilep ////
	xvariable = "MQ";
	int nbinsxvariable_SemiLepZZqqWZqq = 10;
	string binningFileName_SemiLepZZqqWZqq_mu = Binningpath+"Binning_VLQSearch_SemiLepZZqqWZqq_mu.root";
	TemplateTools VarSemiLepZZqqWZqq_mu("SemiLepZZqqWZqq_mu",xvariable,nbinsxvariable_SemiLepZZqqWZqq);
	VarSemiLepZZqqWZqq_mu.SetDatasets(datasetsMu);
	//// WWqq ////
	xvariable = "St";
	int nbinsxvariable_WWqq = 11;
	string binningFileName_WWqq_mu = Binningpath+"Binning_VLQSearch_WWqq_mu.root";
	TemplateTools VarWWqq_mu("WWqq_mu",xvariable,nbinsxvariable_WWqq);
	VarWWqq_mu.SetDatasets(datasetsMu);			
	
	//electron channel
	//// Wqq ////
	xvariable = "MQ";
	//int nbinsxvariable_Wqq = 160;
	//int nbinsxvariable_Wplusqq = 60;
	//int nbinsxvariable_Wminusqq = 100;
	string binningFileName_Wqq_el = Binningpath+"Binning_VLQSearch_Wqq_el.root";
	string binningFileName_Wplusqq_el = Binningpath+"Binning_VLQSearch_Wplusqq_el.root";
	string binningFileName_Wminusqq_el = Binningpath+"Binning_VLQSearch_Wminusqq_el.root";
	TemplateTools VarWqq_el("Wqq_el",xvariable,nbinsxvariable_Wqq);
	TemplateTools VarWplusqq_el("Wplusqq_el",xvariable,nbinsxvariable_Wplusqq);
  TemplateTools VarWminusqq_el("Wminusqq_el",xvariable,nbinsxvariable_Wminusqq);	
	VarWqq_el.SetDatasets(datasetsEl);
	VarWplusqq_el.SetDatasets(datasetsEl);
	VarWminusqq_el.SetDatasets(datasetsEl);
	//// WHqq ////
	xvariable = "MHq";
	//int nbinsxvariable_WHqq = 60;
	string binningFileName_WHqq_el = Binningpath+"Binning_VLQSearch_WHqq_el.root";
	TemplateTools VarWHqq_el("WHqq_el",xvariable,nbinsxvariable_WHqq);
	VarWHqq_el.SetDatasets(datasetsEl);
	// Zqq ////
	xvariable = "MQ";
	//int nbinsxvariable_Zqq = 25;
	string binningFileName_Zqq_el = Binningpath+"Binning_VLQSearch_Zqq_el.root";
	TemplateTools VarZqq_el("Zqq_el",xvariable,nbinsxvariable_Zqq);
	VarZqq_el.SetDatasets(datasetsEl);
	//// ZHqq ////
	xvariable = "MQ";
	//int nbinsxvariable_ZHqq = 4;
	string binningFileName_ZHqq_el = Binningpath+"Binning_VLQSearch_ZHqq_el.root";
	TemplateTools VarZHqq_el("ZHqq_el",xvariable,nbinsxvariable_ZHqq);
	VarZHqq_el.SetDatasets(datasetsEl);
	//// ZZqqWZqq semilep ////
	xvariable = "MQ";
	//int nbinsxvariable_SemiLepZZqqWZqq = 10;
	string binningFileName_SemiLepZZqqWZqq_el = Binningpath+"Binning_VLQSearch_SemiLepZZqqWZqq_el.root";
	TemplateTools VarSemiLepZZqqWZqq_el("SemiLepZZqqWZqq_el",xvariable,nbinsxvariable_SemiLepZZqqWZqq);
	VarSemiLepZZqqWZqq_el.SetDatasets(datasetsEl);
	//// WWqq ////
	xvariable = "St";
	//int nbinsxvariable_WWqq = 8;
	string binningFileName_WWqq_el = Binningpath+"Binning_VLQSearch_WWqq_el.root";
	TemplateTools VarWWqq_el("WWqq_el",xvariable,nbinsxvariable_WWqq);
	VarWWqq_el.SetDatasets(datasetsEl);
	
	if(makeflatbinning == false && useflatbinning == true)
  {
	     cout<<" ... The 'flat' binnings will be loaded!"<<endl;
	     // mu channel
       VarWqq_mu.LoadBinning(binningFileName_Wqq_mu);
       VarWplusqq_mu.LoadBinning(binningFileName_Wplusqq_mu);
			 VarWminusqq_mu.LoadBinning(binningFileName_Wminusqq_mu);
			 VarWHqq_mu.LoadBinning(binningFileName_WHqq_mu);
			 VarZqq_mu.LoadBinning(binningFileName_Zqq_mu);
			 VarZHqq_mu.LoadBinning(binningFileName_ZHqq_mu);
			 VarSemiLepZZqqWZqq_mu.LoadBinning(binningFileName_SemiLepZZqqWZqq_mu);
			 VarWWqq_mu.LoadBinning(binningFileName_WWqq_mu);
			 //el channel
			 VarWqq_el.LoadBinning(binningFileName_Wqq_el);
       VarWplusqq_el.LoadBinning(binningFileName_Wplusqq_el);
			 VarWminusqq_el.LoadBinning(binningFileName_Wminusqq_el);
			 VarWHqq_el.LoadBinning(binningFileName_WHqq_el);
			 VarZqq_el.LoadBinning(binningFileName_Zqq_el);
			 VarZHqq_el.LoadBinning(binningFileName_ZHqq_el);
			 VarSemiLepZZqqWZqq_el.LoadBinning(binningFileName_SemiLepZZqqWZqq_el);
			 VarWWqq_el.LoadBinning(binningFileName_WWqq_el);
  }
  else if(makeflatbinning)
  {
       cout<<" ... The 'flat' binnings will be created!"<<endl;
  }
	
	
	
	
	
	
  //Output ROOT file
  string rootFileName (Outputpath+"VLQTreeAnalyzer"+outputpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //Global variable

  string pathPNG = "VLQSearchPlots_TreeAnalyzer"+inputpostfix;
  pathPNG = pathPNG +"/"; 		
  if(savePNG) mkdir(pathPNG.c_str(),0777);
  
  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);

	//'mu channel' = when at least one muon is trigged 
	MSPlot["MS_nPV_noPUreweighting_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices before PU reweighting", 36, -0.5, 35.5, "Nr. of primary vertices");
  MSPlot["MS_nPV_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices", 36, -0.5, 35.5, "Nr. of primary vertices");  //after PU reweighting
	MSPlot["MS_nPV_noPUreweighting_4jets_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices before PU reweighting", 36, -0.5, 35.5, "Nr. of primary vertices");
  MSPlot["MS_nPV_4jets_mu"] = new MultiSamplePlot(datasetsMu, "nPrimaryVertices", 36, -0.5, 35.5, "Nr. of primary vertices");  //after PU reweighting
	MSPlot["MS_nLeptons_mu"] = new MultiSamplePlot(datasetsMu,"number of leptons", 6, -0.5, 5.5, "number of leptons");
	MSPlot["MS_nMuons_mu"] = new MultiSamplePlot(datasetsMu,"number of muons", 6, -0.5, 5.5, "number of muons");
	MSPlot["MS_nElectrons_mu"] = new MultiSamplePlot(datasetsMu,"number of electrons", 6, -0.5, 5.5, "number of electrons");
	MSPlot["MS_MET_mu"] = new MultiSamplePlot(datasetsMu,"MET", 30, 0, 300, "Missing transverse energy (GeV)");
	MSPlot["MS_MET_singleMu_mu"] = new MultiSamplePlot(datasetsMu,"MET", 30, 0, 300, "Missing transverse energy (GeV)");
	MSPlot["MS_MET_singleMu_Bveto_mu"] = new MultiSamplePlot(datasetsMu,"MET", 30, 0, 300, "Missing transverse energy (GeV)");
	MSPlot["MS_METoverST_mu"] = new MultiSamplePlot(datasetsMu,"METoverST", 100, 0, 1, "Missing transverse energy / S_{T}");
	MSPlot["MS_MTleptonmet_preboxWqq_mu"] = new MultiSamplePlot(datasetsMu,"MT(lepton,MET)", 40, 0, 300, "Transverse mass(lepton,MET) (GeV)");
	MSPlot["MS_HT_mu"] = new MultiSamplePlot(datasetsMu,"HT", 200, 0, 2000, "HT (GeV)");
	//jets (not forward)
	MSPlot["MS_JetMultiplicity_mu"] = new MultiSamplePlot(datasetsMu, "JetMultiplicity", 10, -0.5, 9.5, "jet multiplicity");
	MSPlot["MS_JetMultiplicity_singleMu_mu"] = new MultiSamplePlot(datasetsMu, "JetMultiplicity", 10, -0.5, 9.5, "jet multiplicity");
	MSPlot["MS_Pt_alljets_mu"] = new MultiSamplePlot(datasetsMu,"Pt_alljets", 50, 0, 500, "Pt of all jets (GeV)");
	MSPlot["MS_Pt_jet1_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet1", 50, 0, 500, "Pt of jet 1 (GeV)");
	MSPlot["MS_Pt_jet2_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet2", 50, 0, 500, "Pt of jet 2 (GeV)");
	MSPlot["MS_Pt_jet3_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet3", 50, 0, 500, "Pt of jet 3 (GeV)");
	MSPlot["MS_Pt_jet4_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet4", 50, 0, 500, "Pt of jet 4 (GeV)");
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
	MSPlot["MS_Muon1Pt_dimuon_mu"] = new MultiSamplePlot(datasetsMu,"muon pt", 75, 0, 250, "p_{T} muon 1 (GeV)");
	MSPlot["MS_Muon2Pt_dimuon_mu"] = new MultiSamplePlot(datasetsMu,"muon pt", 75, 0, 250, "p_{T} muon 2 (GeV)");
	MSPlot["MS_Muon1Eta_dimuon_mu"] = new MultiSamplePlot(datasetsMu,"muon eta", 50, -2.5 , 2.5, "#eta muon 1");
	MSPlot["MS_Muon2Eta_dimuon_mu"] = new MultiSamplePlot(datasetsMu,"muon eta", 50, -2.5 , 2.5, "#eta muon 2");	
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
	//ttbar-selection plots
	MSPlot["MS_MET_singleMu_ttbarsel1btag_mu"] = new MultiSamplePlot(datasetsMu,"MET", 30, 0, 300, "Missing transverse energy (GeV)");
	MSPlot["MS_MuonPt_singleMu_ttbarsel1btag_mu"] = new MultiSamplePlot(datasetsMu,"MuonPt", 50, 0, 300, "Muon Pt (GeV)");
	MSPlot["MS_MuonEta_singleMu_ttbarsel1btag_mu"] = new MultiSamplePlot(datasetsMu,"MuonEta", 50, -2.5 , 2.5, "Muon Eta");
	MSPlot["MS_Pt_jet1_singleMu_ttbarsel1btag_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet1", 50, 0, 500, "Pt of jet 1 (GeV)");
	MSPlot["MS_Pt_jet2_singleMu_ttbarsel1btag_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet2", 50, 0, 500, "Pt of jet 2 (GeV)");
	MSPlot["MS_Eta_jet1_singleMu_ttbarsel1btag_mu"] = new MultiSamplePlot(datasetsMu,"Eta_jet1", 50, -2.5 , 2.5, "Eta of jet 1");
	MSPlot["MS_Eta_jet2_singleMu_ttbarsel1btag_mu"] = new MultiSamplePlot(datasetsMu,"Eta_jet2", 50, -2.5 , 2.5, "Eta of jet 2");
	MSPlot["MS_HT_singleMu_ttbarsel1btag_mu"] = new MultiSamplePlot(datasetsMu,"HT", 50, 200, 1500, "H_{T} (GeV)");	
	MSPlot["MS_MET_singleMu_ttbarsel2btag_mu"] = new MultiSamplePlot(datasetsMu,"MET", 30, 0, 300, "Missing transverse energy (GeV)");
	MSPlot["MS_MuonPt_singleMu_ttbarsel2btag_mu"] = new MultiSamplePlot(datasetsMu,"MuonPt", 50, 0, 300, "Muon Pt (GeV)");
	MSPlot["MS_MuonEta_singleMu_ttbarsel2btag_mu"] = new MultiSamplePlot(datasetsMu,"MuonEta", 50, -2.5 , 2.5, "Muon Eta");	
	MSPlot["MS_Pt_jet1_singleMu_ttbarsel2btag_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet1", 50, 0, 500, "p_{T} jet 1 (GeV)");
	MSPlot["MS_Pt_jet2_singleMu_ttbarsel2btag_mu"] = new MultiSamplePlot(datasetsMu,"Pt_jet2", 50, 0, 500, "p_{T} jet 2 (GeV)");
	MSPlot["MS_Eta_jet1_singleMu_ttbarsel2btag_mu"] = new MultiSamplePlot(datasetsMu,"Eta_jet1", 50, -2.5 , 2.5, "#eta jet 1");
	MSPlot["MS_Eta_jet2_singleMu_ttbarsel2btag_mu"] = new MultiSamplePlot(datasetsMu,"Eta_jet2", 50, -2.5 , 2.5, "#eta jet 2");
	MSPlot["MS_HT_singleMu_ttbarsel2btag_mu"] = new MultiSamplePlot(datasetsMu,"HT", 50, 200, 1500, "H_{T} (GeV)");
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
	MSPlot["MS_St_boxWWqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box WWqq)", 20, 200, 2200, "S_{T} (GeV)");
	MSPlot["MS_nEvts_boxWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxZZqq_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box ZZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqqZZqq_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq+ZZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_Zwindow_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxZZqq_Zwindow_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box ZZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq+ZZqq)", 1, 0.5, 1.5, "");
	//splitting multilepton boxes
	MSPlot["MS_nEvts_boxWZqq_MuElEl_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_MuMuEl_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_MuMuMu_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_Zwindow_MuElEl_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_Zwindow_MuMuEl_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_Zwindow_MuMuMu_mu"] = new MultiSamplePlot(datasetsMu,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
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
	MSPlot["MS_MET_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"MET (box Wqq)", 30, 0, 300, "Missing transverse energy (GeV)");
 	MSPlot["MS_MTleptonmet_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"MT(lepton,MET)", 40, 0, 300, "Transverse mass(lepton,MET) (GeV)");
	MSPlot["MS_HT_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"HT (box Wqq)", 50, 200, 1500, "H_{T} (GeV)");
	MSPlot["MS_LeadingJetPt_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Pt (box Wqq)", 50, 0, 400, "leading-jet Pt (GeV)");
	MSPlot["MS_LeadingJetEta_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Eta (box Wqq)", 50, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubLeadingJetPt_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Pt (box Wqq)", 50, 0, 400, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetEta_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Eta (box Wqq)", 50, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_ForwardJetPt_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Forward jet Pt (box Wqq)", 50, 0, 300, "forward-jet Pt (GeV)");
	MSPlot["MS_ForwardJetEta_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Forward jet Eta (box Wqq)", 100, -5 , 5, "forward-jet Eta");
	MSPlot["MS_nJets_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"nJets (box Wqq)", 4, -0.5, 3.5, "number of central jets");
	MSPlot["MS_MuonPt_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Muon Pt (box Wqq)", 50, 0, 300, "muon Pt (GeV)");
	MSPlot["MS_MuonEta_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Muon Eta (box Wqq)", 50, -2.5 , 2.5, "muon Eta");
	MSPlot["MS_MCeventWeight_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"MCeventWeight (box Wqq)", 50, 0, 10, "MC event weight");
	MSPlot["MS_MET_preboxWqq_noBveto_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 30, 0, 300, "Missing transverse energy (GeV)");
	MSPlot["MS_MET_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 30, 0, 300, "Missing transverse energy (GeV)");
	MSPlot["MS_HT_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"HT (prebox Wqq)", 50, 200, 1500, "H_{T} (GeV)");
	MSPlot["MS_LeadingJetPt_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Pt (prebox Wqq)", 50, 0, 600, "leading-jet Pt (GeV)");
	MSPlot["MS_LeadingJetEta_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Eta (prebox Wqq)", 50, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubLeadingJetPt_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Pt (prebox Wqq)", 50, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetEta_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Eta (prebox Wqq)", 50, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_ForwardJetPt_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"Forward jet Pt (prebox Wqq)", 50, 0, 300, "forward-jet Pt (GeV)");
	MSPlot["MS_ForwardJetEta_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"Forward jet Eta (prebox Wqq)", 100, -5 , 5, "forward-jet Eta");
	MSPlot["MS_nJets_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"nJets (prebox Wqq)", 4, -0.5, 3.5, "number of central jets");
	MSPlot["MS_MuonPt_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"Muon Pt (prebox Wqq)", 50, 0, 300, "muon Pt (GeV)");	
	MSPlot["MS_MuonEta_preboxWqq_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"Muon Eta (prebox Wqq)", 50, -2.5 , 2.5, "muon Eta");
	MSPlot["MS_MET_preboxWqq_noleadingjetcut_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 30, 0, 300, "Missing transverse energy (GeV)");
	MSPlot["MS_HT_preboxWqq_noleadingjetcut_mu"] = new MultiSamplePlot(datasetsMu,"HT (prebox Wqq)", 50, 0, 1300, "H_{T} (GeV)");
	MSPlot["MS_LeadingJetPt_preboxWqq_noleadingjetcut_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Pt (prebox Wqq)", 50, 0, 600, "leading-jet Pt (GeV)");
	MSPlot["MS_LeadingJetEta_preboxWqq_noleadingjetcut_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Eta (prebox Wqq)", 50, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubLeadingJetPt_preboxWqq_noleadingjetcut_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Pt (prebox Wqq)", 50, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetEta_preboxWqq_noleadingjetcut_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Eta (prebox Wqq)", 50, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_ForwardJetPt_preboxWqq_noleadingjetcut_mu"] = new MultiSamplePlot(datasetsMu,"Forward jet Pt (prebox Wqq)", 50, 0, 300, "forward-jet Pt (GeV)");
	MSPlot["MS_ForwardJetEta_preboxWqq_noleadingjetcut_mu"] = new MultiSamplePlot(datasetsMu,"Forward jet Eta (prebox Wqq)", 100, -5 , 5, "forward-jet Eta");
	//MSPlot["MS_MET_preboxWqq_noforwardjetreq_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
  MSPlot["MS_LeadingJetPt_preboxWqq_leadingjet150-200_METcut_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Pt (prebox Wqq)", 20, 150, 200, "leading-jet Pt (GeV)");
  MSPlot["MS_LeadingJetEta_preboxWqq_leadingjet150-200_METcut_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Eta (prebox Wqq)", 20, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubleadingJetPt_preboxWqq_leadingjet150-200_METcut_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Pt (prebox Wqq)", 20, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubleadingJetEta_preboxWqq_leadingjet150-200_METcut_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Eta (prebox Wqq)", 20, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_HT_preboxWqq_leadingjet150-200_METcut_mu"] = new MultiSamplePlot(datasetsMu,"HT (prebox Wqq)", 20, 150, 1000, "HT (GeV)");
  MSPlot["MS_St_preboxWqq_leadingjet150-200_METcut_mu"] = new MultiSamplePlot(datasetsMu,"ST (prebox Wqq)", 20, 150, 1000, "ST (GeV)");
	MSPlot["MS_MET_preboxWqq_leadingjet150-200_METcut_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 30, 0, 300, "MET (GeV)");
	MSPlot["MS_MuonPt_preboxWqq_leadingjet150-200_METcut_mu"] = new MultiSamplePlot(datasetsMu,"Muon Pt (prebox Wqq)", 20, 0, 300, "muon Pt (GeV)");
	MSPlot["MS_MuonEta_preboxWqq_leadingjet150-200_METcut_mu"] = new MultiSamplePlot(datasetsMu,"Muon Eta (prebox Wqq)", 20, -2.5 , 2.5, "muon Eta");
	MSPlot["MS_MET_preboxWqq_leadingjet150-200_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 30, 0, 300, "MET (GeV)");
  //more box plots	
	MSPlot["MS_DeltaR_WcandLeadingJet_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"DeltaR(W,LeadingJet) (box Wqq)", 30, 0, 6, "#DeltaR(W,LeadingJet)");
	MSPlot["MS_DeltaEta_WcandLeadingJet_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"DeltaEta(W,LeadingJet) (box Wqq)", 30, 0, 3, "#Delta#eta(W,LeadingJet)");
	MSPlot["MS_DeltaR_LeptonNeutrino_boxWqq_mu"] = new MultiSamplePlot(datasetsMu,"DeltaR(lepton,neutrino) (box Wqq)", 30, 0, 3, "#DeltaR(lepton,neutrino)");
	MSPlot["MS_DeltaR_ZcandLeadingJet_boxZqq_mu"] = new MultiSamplePlot(datasetsMu,"DeltaR(Z,LeadingJet) (box Zqq)", 30, 0, 6, "#DeltaR(Z,LeadingJet)");
	MSPlot["MS_DeltaEta_ZcandLeadingJet_boxZqq_mu"] = new MultiSamplePlot(datasetsMu,"DeltaEta(Z,LeadingJet) (box Zqq)", 30, 0, 3, "#Delta#eta(Z,LeadingJet)");
	MSPlot["MS_DeltaR_LeptonLepton_boxZqq_mu"] = new MultiSamplePlot(datasetsMu,"DeltaR(lepton,lepton) (box Zqq)", 30, 0, 3, "#DeltaR(lepton,lepton)");
	//Wqq prebox with all cuts except forward-jet requirement and leading jet pt cut (only 150-200)
	MSPlot["MS_LeadingJetPt_preboxWqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Pt (prebox Wqq)", 20, 150, 200, "leading-jet Pt (GeV)");
  MSPlot["MS_LeadingJetEta_preboxWqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Eta (prebox Wqq)", 20, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubleadingJetPt_preboxWqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Pt (prebox Wqq)", 20, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubleadingJetEta_preboxWqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Eta (prebox Wqq)", 20, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_HT_preboxWqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"HT (prebox Wqq)", 20, 150, 1000, "HT (GeV)");
  MSPlot["MS_St_preboxWqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"ST (prebox Wqq)", 20, 150, 1000, "ST (GeV)");
	MSPlot["MS_MET_preboxWqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 30, 0, 300, "MET (GeV)");
  MSPlot["MS_MTleptonmet_preboxWqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"MT(lepton,MET)", 40, 0, 300, "Transverse mass(lepton,MET) (GeV)");
	MSPlot["MS_MuonPt_preboxWqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Muon Pt (prebox Wqq)", 20, 0, 300, "muon Pt (GeV)");
	MSPlot["MS_MuonEta_preboxWqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Muon Eta (prebox Wqq)", 20, -2.5 , 2.5, "muon Eta");
  MSPlot["MS_nJets_preboxWqq_leadingjet150-200_noForwardJet_noNumberCentralJetCut_mu"] = new MultiSamplePlot(datasetsMu,"nJets (prebox Wqq)", 8, -0.5, 7.5, "number of central jets");
	MSPlot["MS_nJets_preboxWqq_leadingjet150-200_noNumberCentralJetCut_mu"] = new MultiSamplePlot(datasetsMu,"nJets (prebox Wqq)", 8, -0.5, 7.5, "number of central jets");
	MSPlot["MS_MET_preboxWqq_leadingjet150-200_noForwardJet_noMETMTcut_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Wqq)", 30, 0, 300, "MET (GeV)");	
	//other boxes
	MSPlot["MS_Mleptons_preboxZqq_mu"] = new MultiSamplePlot(datasetsMu,"M(leptons) (prebox Zqq/ZHqq)", 100, 0, 200, "M(leptons) (GeV)");
	MSPlot["MS_MET_boxZHqq_mu"] = new MultiSamplePlot(datasetsMu,"MET (box ZHqq)", 30, 0, 300, "MET (GeV)");
	MSPlot["MS_MQ_antiboxZHqq_offZ_mu"] = new MultiSamplePlot(datasetsMu,"Mass Q (antibox ZHqq off-Z)", 15, 0, 1500, "Mass(lepton,lepton,highptjet) (GeV)");
	//Zqq prebox with all cuts except forward-jet requirement
	MSPlot["MS_LeadingJetPt_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Pt (prebox Zqq)", 20, 150, 200, "leading-jet Pt (GeV)");
  MSPlot["MS_LeadingJetEta_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Eta (prebox Zqq)", 20, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubleadingJetPt_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Pt (prebox Zqq)", 20, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubleadingJetEta_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Eta (prebox Zqq)", 20, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_HT_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"HT (prebox Zqq)", 20, 150, 1000, "HT (GeV)");
  MSPlot["MS_St_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"ST (prebox Zqq)", 20, 150, 1000, "ST (GeV)");
	MSPlot["MS_MET_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox Zqq)", 20, 0, 100, "MET (GeV)");
	MSPlot["MS_Muon1Pt_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Muon 1 Pt (prebox Zqq)", 20, 0, 300, "muon 1 Pt (GeV)");
	MSPlot["MS_Muon1Eta_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Muon 1 Eta (prebox Zqq)", 20, -2.5 , 2.5, "muon 1 Eta");
  MSPlot["MS_Muon2Pt_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Muon 2 Pt (prebox Zqq)", 20, 0, 300, "muon 2 Pt (GeV)");
	MSPlot["MS_Muon2Eta_preboxZqq_leadingjet150-200_noForwardJet_mu"] = new MultiSamplePlot(datasetsMu,"Muon 2 Eta (prebox Zqq)", 20, -2.5 , 2.5, "muon 2 Eta");
	MSPlot["MS_nJets_preboxZqq_leadingjet150-200_noForwardJet_noNumberCentralJetCut_mu"] = new MultiSamplePlot(datasetsMu,"nJets (prebox Zqq)", 8, -0.5, 7.5, "number of central jets");
	MSPlot["MS_nJets_preboxZqq_leadingjet150-200_noNumberCentralJetCut_mu"] = new MultiSamplePlot(datasetsMu,"nJets (prebox Zqq)", 8, -0.5, 7.5, "number of central jets");
  //box SemiLepZZqqWZqq
	MSPlot["MS_MQ_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Mass Q (box SemiLepZZqqWZqq)", 40, 200, 1500, "Mass(lepton,lepton,leadingjet) (GeV)");
	MSPlot["MS_MET_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"MET (box SemiLepZZqqWZqq)", 30, 0, 300, "MET (GeV)");
	MSPlot["MS_St_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"ST (box SemiLepZZqqWZqq)", 25, 200, 1700, "S_{T} (GeV)");
	MSPlot["MS_HT_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"HT (box SemiLepZZqqWZqq)", 25, 200, 1500, "H_{T} (GeV)");
	MSPlot["MS_LeadingJetPt_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Pt (box SemiLepZZqqWZqq)", 50, 0, 800, "leading-jet Pt (GeV)");
	MSPlot["MS_LeadingJetEta_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Eta (box SemiLepZZqqWZqq)", 50, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubLeadingJetPt_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Pt (box SemiLepZZqqWZqq)", 50, 0, 400, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetEta_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Eta (box SemiLepZZqqWZqq)", 50, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_nJets_boxSemiLepZZqqWZqq_mu"]  = new MultiSamplePlot(datasetsMu,"nJets (box SemiLepZZqqWZqq)", 8, -0.5, 7.5, "number of central jets");
  MSPlot["MS_Muon1Pt_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Muon 1 Pt (box SemiLepZZqqWZqq)", 20, 0, 300, "muon 1 Pt (GeV)");
	MSPlot["MS_Muon1Eta_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Muon 1 Eta (box SemiLepZZqqWZqq)", 20, -2.5 , 2.5, "muon 1 Eta");
  MSPlot["MS_Muon2Pt_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Muon 2 Pt (box SemiLepZZqqWZqq)", 20, 0, 300, "muon 2 Pt (GeV)");
	MSPlot["MS_Muon2Eta_boxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Muon 2 Eta (box SemiLepZZqqWZqq)", 20, -2.5 , 2.5, "muon 2 Eta");
  //other preboxes
	MSPlot["MS_PtZ_preboxZHqq_mu"] = new MultiSamplePlot(datasetsMu,"Pt Z candidate (prebox ZHqq)", 25, 0, 500, "Pt Z candidate (GeV)");
	MSPlot["MS_PtZ_preboxZqq_mu"] = new MultiSamplePlot(datasetsMu,"Pt Z candidate (prebox Zqq)", 25, 0, 500, "Pt Z candidate (GeV)");
  MSPlot["MS_PtZ_preboxZHqq_LowerJetCuts_LowerZpt_mu"] = new MultiSamplePlot(datasetsMu,"Pt Z candidate (prebox ZHqq)", 15, 0, 250, "Pt Z candidate (GeV)");
  MSPlot["MS_St_preboxZHqq_LowerJetCuts_LowerZpt_mu"] = new MultiSamplePlot(datasetsMu,"ST (prebox ZHqq)", 50, 200, 1700, "S_{T} (GeV)");
	MSPlot["MS_MET_preboxZHqq_LowerJetCuts_LowerZpt_mu"] = new MultiSamplePlot(datasetsMu,"MET (prebox ZHqq)", 30, 0, 300, "MET (GeV)");
	MSPlot["MS_LeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Pt (prebox ZHqq)", 40, 0, 200, "leading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Pt (prebox ZHqq)", 40, 0, 200, "subleading-jet Pt (GeV)");
	MSPlot["MS_LeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_mu"] = new MultiSamplePlot(datasetsMu,"Leading jet Pt (prebox WHqq)", 40, 0, 200, "leading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_mu"] = new MultiSamplePlot(datasetsMu,"Subleading jet Pt (prebox WHqq)", 40, 0, 200, "subleading-jet Pt (GeV)");	
	MSPlot["MS_PtH_preboxWHqq_mu"] = new MultiSamplePlot(datasetsMu,"Pt H candidate (prebox WHqq)", 25, 0, 500, "Pt H candidate (GeV)");
	MSPlot["MS_PtZ_preboxWZqq_Zwindow_mu"] = new MultiSamplePlot(datasetsMu,"Pt Z candidate (prebox WZqq)", 25, 0, 500, "Pt Z candidate (GeV)");
	MSPlot["MS_PtH_preboxWHqq_LowerJetCuts_LowerHpt_mu"] = new MultiSamplePlot(datasetsMu,"Pt H candidate (prebox WHqq)", 15, 0, 250, "Pt H candidate (GeV)");
	MSPlot["MS_St_preboxWHqq_LowerJetCuts_LowerHpt_offH_mu"] = new MultiSamplePlot(datasetsMu,"ST (prebox WHqq)", 50, 200, 1700, "S_{T} (GeV)");
  MSPlot["MS_MV_preboxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"M(dijet) (prebox SemiLepZZqqWZqq)", 50, 0, 300, "M(dijet) (GeV)");
	MSPlot["MS_PtV_preboxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Pt V candidate (prebox SemiLepZZqqWZqq)", 25, 0, 500, "Pt hadronic W/Z candidate (GeV)");
	MSPlot["MS_PtZ_preboxSemiLepZZqqWZqq_mu"] = new MultiSamplePlot(datasetsMu,"Pt Z candidate (prebox SemiLepZZqqWZqq)", 25, 0, 500, "Pt Z candidate (GeV)");
	MSPlot["MS_PtW_preboxWqq_mu"] = new MultiSamplePlot(datasetsMu,"Pt W candidate (prebox Wqq)", 25, 0, 500, "Pt W candidate (GeV)");
	MSPlot["MS_MHcandidate_preboxWHqq_mu"] = new MultiSamplePlot(datasetsMu,"MH (prebox WHqq)", 50, 0, 500, "reconstructed Higgs mass (GeV)");
	
	
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
	MSPlot["MS_MET_singleEl_Bveto_el"] = new MultiSamplePlot(datasetsEl,"MET", 200, 0, 1000, "Missing transverse energy (GeV)");
	MSPlot["MS_METoverST_el"] = new MultiSamplePlot(datasetsEl,"METoverST", 100, 0, 1, "Missing transverse energy / S_{T}");
	MSPlot["MS_MTleptonmet_preboxWqq_el"] = new MultiSamplePlot(datasetsEl,"MT(lepton,MET)", 40, 0, 300, "Transverse mass(lepton,MET) (GeV)");
	MSPlot["MS_HT_el"] = new MultiSamplePlot(datasetsEl,"HT", 200, 0, 2000, "HT (GeV)");
	//jets (not forward)
	MSPlot["MS_JetMultiplicity_el"] = new MultiSamplePlot(datasetsEl, "JetMultiplicity", 10, -0.5, 9.5, "jet multiplicity");
	MSPlot["MS_JetMultiplicity_singleEl_el"] = new MultiSamplePlot(datasetsEl, "JetMultiplicity", 10, -0.5, 9.5, "jet multiplicity");
	MSPlot["MS_Pt_alljets_el"] = new MultiSamplePlot(datasetsEl,"Pt_alljets", 50, 0, 500, "Pt of all jets (GeV)");
	MSPlot["MS_Pt_jet1_el"] = new MultiSamplePlot(datasetsEl,"Pt_jet1", 50, 0, 500, "Pt of jet 1 (GeV)");
	MSPlot["MS_Pt_jet2_el"] = new MultiSamplePlot(datasetsEl,"Pt_jet2", 50, 0, 500, "Pt of jet 2 (GeV)");
	MSPlot["MS_Pt_jet3_el"] = new MultiSamplePlot(datasetsEl,"Pt_jet3", 50, 0, 500, "Pt of jet 3 (GeV)");
	MSPlot["MS_Pt_jet4_el"] = new MultiSamplePlot(datasetsEl,"Pt_jet4", 50, 0, 500, "Pt of jet 4 (GeV)");
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
	MSPlot["MS_Electron1Pt_dielectron_el"] = new MultiSamplePlot(datasetsEl,"electron pt", 75, 0, 250, "p_{T} electron 1 (GeV)");
	MSPlot["MS_Electron2Pt_dielectron_el"] = new MultiSamplePlot(datasetsEl,"electron pt", 75, 0, 250, "p_{T} electron 2 (GeV)");
	MSPlot["MS_Electron1Eta_dielectron_el"] = new MultiSamplePlot(datasetsEl,"electron eta", 50, -2.5 , 2.5, "#eta electron 1");
	MSPlot["MS_Electron2Eta_dielectron_el"] = new MultiSamplePlot(datasetsEl,"electron eta", 50, -2.5 , 2.5, "#eta electron 2");
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
	MSPlot["MS_St_boxWWqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box WWqq)", 20, 200, 2200, "S_{T} (GeV)");
	MSPlot["MS_nEvts_boxWZqq_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxZZqq_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box ZZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqqZZqq_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq+ZZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_Zwindow_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxZZqq_Zwindow_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box ZZqq)", 1, 0.5, 1.5, "");	
	MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq+ZZqq)", 1, 0.5, 1.5, "");
	//splitting multilepton box
	MSPlot["MS_nEvts_boxWZqq_ElElEl_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_ElElMu_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_ElMuMu_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_Zwindow_ElElEl_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_Zwindow_ElElMu_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
	MSPlot["MS_nEvts_boxWZqq_Zwindow_ElMuMu_el"] = new MultiSamplePlot(datasetsEl,"Number of events (box WZqq)", 1, 0.5, 1.5, "");
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
	MSPlot["MS_MET_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"MET (box Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	MSPlot["MS_MTleptonmet_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"MT(lepton,MET)", 40, 0, 300, "Transverse mass(lepton,MET) (GeV)");
	MSPlot["MS_HT_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"HT (box Wqq)", 50, 200, 1500, "H_{T} (GeV)");
	MSPlot["MS_LeadingJetPt_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Pt (box Wqq)", 50, 0, 600, "leading-jet Pt (GeV)");
	MSPlot["MS_LeadingJetEta_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Eta (box Wqq)", 50, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubLeadingJetPt_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Pt (box Wqq)", 50, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetEta_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Eta (box Wqq)", 50, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_ForwardJetPt_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Forward jet Pt (box Wqq)", 50, 0, 300, "forward-jet Pt (GeV)");
	MSPlot["MS_ForwardJetEta_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Forward jet Eta (box Wqq)", 100, -5 , 5, "forward-jet Eta");
	MSPlot["MS_nJets_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"nJets (box Wqq)", 4, -0.5, 3.5, "number of central jets");
	MSPlot["MS_ElectronPt_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Electron Pt (box Wqq)", 50, 0, 300, "electron Pt (GeV)");
	MSPlot["MS_ElectronEta_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"Electron Eta (box Wqq)", 50, -2.5 , 2.5, "electron Eta");
	MSPlot["MS_MCeventWeight_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"MCeventWeight (box Wqq)", 50, 0, 10, "MC event weight");
	MSPlot["MS_MET_preboxWqq_noBveto_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	MSPlot["MS_MET_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	MSPlot["MS_HT_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"HT (prebox Wqq)", 50, 200, 1500, "H_{T} (GeV)");
	MSPlot["MS_LeadingJetPt_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Pt (prebox Wqq)", 50, 0, 600, "leading-jet Pt (GeV)");
	MSPlot["MS_LeadingJetEta_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Eta (prebox Wqq)", 50, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubLeadingJetPt_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Pt (prebox Wqq)", 50, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetEta_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Eta (prebox Wqq)", 50, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_ForwardJetPt_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"Forward jet Pt (prebox Wqq)", 50, 0, 300, "forward-jet Pt (GeV)");
	MSPlot["MS_ForwardJetEta_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"Forward jet Eta (prebox Wqq)", 100, -5 , 5, "forward-jet Eta");
	MSPlot["MS_nJets_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"nJets (prebox Wqq)", 4, -0.5, 3.5, "number of central jets");
	MSPlot["MS_ElectronPt_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"Electron Pt (prebox Wqq)", 50, 0, 300, "electron Pt (GeV)");	
	MSPlot["MS_ElectronEta_preboxWqq_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"Electron Eta (prebox Wqq)", 50, -2.5 , 2.5, "electron Eta");
	MSPlot["MS_MET_preboxWqq_noleadingjetcut_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
	MSPlot["MS_HT_preboxWqq_noleadingjetcut_el"] = new MultiSamplePlot(datasetsEl,"HT (prebox Wqq)", 50, 0, 1300, "H_{T} (GeV)");
	MSPlot["MS_LeadingJetPt_preboxWqq_noleadingjetcut_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Pt (prebox Wqq)", 50, 0, 600, "leading-jet Pt (GeV)");
	MSPlot["MS_LeadingJetEta_preboxWqq_noleadingjetcut_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Eta (prebox Wqq)", 50, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubLeadingJetPt_preboxWqq_noleadingjetcut_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Pt (prebox Wqq)", 50, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetEta_preboxWqq_noleadingjetcut_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Eta (prebox Wqq)", 50, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_ForwardJetPt_preboxWqq_noleadingjetcut_el"] = new MultiSamplePlot(datasetsEl,"Forward jet Pt (prebox Wqq)", 50, 0, 300, "forward-jet Pt (GeV)");
	MSPlot["MS_ForwardJetEta_preboxWqq_noleadingjetcut_el"] = new MultiSamplePlot(datasetsEl,"Forward jet Eta (prebox Wqq)", 100, -5 , 5, "forward-jet Eta");
	//MSPlot["MS_MET_preboxWqq_noforwardjetreq_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 50, 0, 500, "Missing transverse energy (GeV)");
  MSPlot["MS_LeadingJetPt_preboxWqq_leadingjet150-200_METcut_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Pt (prebox Wqq)", 20, 150, 200, "leading-jet Pt (GeV)");
  MSPlot["MS_LeadingJetEta_preboxWqq_leadingjet150-200_METcut_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Eta (prebox Wqq)", 20, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubleadingJetPt_preboxWqq_leadingjet150-200_METcut_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Pt (prebox Wqq)", 20, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubleadingJetEta_preboxWqq_leadingjet150-200_METcut_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Eta (prebox Wqq)", 20, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_HT_preboxWqq_leadingjet150-200_METcut_el"] = new MultiSamplePlot(datasetsEl,"HT (prebox Wqq)", 20, 150, 1000, "HT (GeV)");
  MSPlot["MS_St_preboxWqq_leadingjet150-200_METcut_el"] = new MultiSamplePlot(datasetsEl,"ST (prebox Wqq)", 20, 150, 1000, "ST (GeV)");
	MSPlot["MS_MET_preboxWqq_leadingjet150-200_METcut_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 20, 0, 300, "MET (GeV)");
	MSPlot["MS_ElectronPt_preboxWqq_leadingjet150-200_METcut_el"] = new MultiSamplePlot(datasetsEl,"Electron Pt (prebox Wqq)", 20, 0, 300, "electron Pt (GeV)");
	MSPlot["MS_ElectronEta_preboxWqq_leadingjet150-200_METcut_el"] = new MultiSamplePlot(datasetsEl,"Electron Eta (prebox Wqq)", 20, -2.5 , 2.5, "electron Eta");
  MSPlot["MS_MET_preboxWqq_leadingjet150-200_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 20, 0, 300, "MET (GeV)");
  //more box plots	
	MSPlot["MS_DeltaR_WcandLeadingJet_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"DeltaR(W,LeadingJet) (box Wqq)", 30, 0, 6, "#DeltaR(W,LeadingJet)");
	MSPlot["MS_DeltaEta_WcandLeadingJet_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"DeltaEta(W,LeadingJet) (box Wqq)", 30, 0, 3, "#Delta#eta(W,LeadingJet)");
	MSPlot["MS_DeltaR_LeptonNeutrino_boxWqq_el"] = new MultiSamplePlot(datasetsEl,"DeltaR(lepton,neutrino) (box Wqq)", 30, 0, 3, "#DeltaR(lepton,neutrino)");
	MSPlot["MS_DeltaR_ZcandLeadingJet_boxZqq_el"] = new MultiSamplePlot(datasetsEl,"DeltaR(Z,LeadingJet) (box Zqq)", 30, 0, 6, "#DeltaR(Z,LeadingJet)");
	MSPlot["MS_DeltaEta_ZcandLeadingJet_boxZqq_el"] = new MultiSamplePlot(datasetsEl,"DeltaEta(Z,LeadingJet) (box Zqq)", 30, 0, 3, "#Delta#eta(Z,LeadingJet)");
	MSPlot["MS_DeltaR_LeptonLepton_boxZqq_el"] = new MultiSamplePlot(datasetsEl,"DeltaR(lepton,lepton) (box Zqq)", 30, 0, 3, "#DeltaR(lepton,lepton)");
	//Wqq prebox with all cuts except forward-jet requirement
	MSPlot["MS_LeadingJetPt_preboxWqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Pt (prebox Wqq)", 20, 150, 200, "leading-jet Pt (GeV)");
  MSPlot["MS_LeadingJetEta_preboxWqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Eta (prebox Wqq)", 20, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubleadingJetPt_preboxWqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Pt (prebox Wqq)", 20, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubleadingJetEta_preboxWqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Eta (prebox Wqq)", 20, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_HT_preboxWqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"HT (prebox Wqq)", 20, 150, 1000, "HT (GeV)");
  MSPlot["MS_St_preboxWqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"ST (prebox Wqq)", 20, 150, 1000, "ST (GeV)");
	MSPlot["MS_MET_preboxWqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 20, 0, 300, "MET (GeV)");
  MSPlot["MS_MTleptonmet_preboxWqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"MT(lepton,MET)", 40, 0, 300, "Transverse mass(lepton,MET) (GeV)");
	MSPlot["MS_ElectronPt_preboxWqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Electron Pt (prebox Wqq)", 20, 0, 300, "electron Pt (GeV)");
	MSPlot["MS_ElectronEta_preboxWqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Electron Eta (prebox Wqq)", 20, -2.5 , 2.5, "electron Eta");
	MSPlot["MS_nJets_preboxWqq_leadingjet150-200_noForwardJet_noNumberCentralJetCut_el"] = new MultiSamplePlot(datasetsEl,"nJets (prebox Wqq)", 8, -0.5, 7.5, "number of central jets");
  MSPlot["MS_nJets_preboxWqq_leadingjet150-200_noNumberCentralJetCut_el"] = new MultiSamplePlot(datasetsEl,"nJets (prebox Wqq)", 8, -0.5, 7.5, "number of central jets");
	MSPlot["MS_MET_preboxWqq_leadingjet150-200_noForwardJet_noMETMTcut_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Wqq)", 30, 0, 300, "MET (GeV)");
	//other boxes
	MSPlot["MS_Mleptons_preboxZqq_el"] = new MultiSamplePlot(datasetsEl,"M(leptons) (prebox Zqq/ZHqq)", 100, 0, 200, "M(leptons) (GeV)");
	MSPlot["MS_MET_boxZHqq_el"] = new MultiSamplePlot(datasetsEl,"MET (box ZHqq)", 20, 0, 300, "MET (GeV)");
	MSPlot["MS_MQ_antiboxZHqq_offZ_el"] = new MultiSamplePlot(datasetsEl,"Mass Q (antibox ZHqq off-Z)", 15, 0, 1500, "Mass(lepton,lepton,highptjet) (GeV)");
  //Zqq prebox with all cuts except forward-jet requirement
	MSPlot["MS_LeadingJetPt_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Pt (prebox Zqq)", 20, 150, 200, "leading-jet Pt (GeV)");
  MSPlot["MS_LeadingJetEta_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Eta (prebox Zqq)", 20, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubleadingJetPt_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Pt (prebox Zqq)", 20, 0, 300, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubleadingJetEta_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Eta (prebox Zqq)", 20, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_HT_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"HT (prebox Zqq)", 20, 150, 1000, "HT (GeV)");
  MSPlot["MS_St_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"ST (prebox Zqq)", 20, 150, 1000, "ST (GeV)");
	MSPlot["MS_MET_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox Zqq)", 20, 0, 100, "MET (GeV)");
	MSPlot["MS_Electron1Pt_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Electron 1 Pt (prebox Zqq)", 20, 0, 300, "electron Pt (GeV)");
	MSPlot["MS_Electron1Eta_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Electron 1 Eta (prebox Zqq)", 20, -2.5 , 2.5, "electron Eta");
  MSPlot["MS_Electron2Pt_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Electron 2 Pt (prebox Zqq)", 20, 0, 300, "electron Pt (GeV)");
	MSPlot["MS_Electron2Eta_preboxZqq_leadingjet150-200_noForwardJet_el"] = new MultiSamplePlot(datasetsEl,"Electron 2 Eta (prebox Zqq)", 20, -2.5 , 2.5, "electron Eta");
	MSPlot["MS_nJets_preboxZqq_leadingjet150-200_noForwardJet_noNumberCentralJetCut_el"] = new MultiSamplePlot(datasetsEl,"nJets (prebox Zqq)", 8, -0.5, 7.5, "number of central jets");
	MSPlot["MS_nJets_preboxZqq_leadingjet150-200_noNumberCentralJetCut_el"] = new MultiSamplePlot(datasetsEl,"nJets (prebox Zqq)", 8, -0.5, 7.5, "number of central jets");
  //box SemiLepZZqqWZqq
	MSPlot["MS_MQ_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Mass Q (box SemiLepZZqqWZqq)", 40, 200, 1500, "Mass(lepton,lepton,leadingjet) (GeV)");
	MSPlot["MS_MET_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"MET (box SemiLepZZqqWZqq)", 20, 0, 300, "MET (GeV)");
	MSPlot["MS_St_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"ST (box SemiLepZZqqWZqq)", 25, 200, 1700, "S_{T} (GeV)");
	MSPlot["MS_HT_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"HT (box SemiLepZZqqWZqq)", 25, 200, 1500, "H_{T} (GeV)");
	MSPlot["MS_LeadingJetPt_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Pt (box SemiLepZZqqWZqq)", 50, 0, 800, "leading-jet Pt (GeV)");
	MSPlot["MS_LeadingJetEta_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Eta (box SemiLepZZqqWZqq)", 50, -2.5 , 2.5, "leading-jet Eta");
	MSPlot["MS_SubLeadingJetPt_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Pt (box SemiLepZZqqWZqq)", 50, 0, 400, "subleading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetEta_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Eta (box SemiLepZZqqWZqq)", 50, -2.5 , 2.5, "subleading-jet Eta");
	MSPlot["MS_nJets_boxSemiLepZZqqWZqq_el"]  = new MultiSamplePlot(datasetsEl,"nJets (box SemiLepZZqqWZqq)", 8, -0.5, 7.5, "number of central jets");
  MSPlot["MS_Electron1Pt_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Electron 1 Pt (box SemiLepZZqqWZqq)", 20, 0, 300, "electron 1 Pt (GeV)");
	MSPlot["MS_Electron1Eta_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Electron 1 Eta (box SemiLepZZqqWZqq)", 20, -2.5 , 2.5, "electron 1 Eta");
  MSPlot["MS_Electron2Pt_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Electron 2 Pt (box SemiLepZZqqWZqq)", 20, 0, 300, "electron 2 Pt (GeV)");
	MSPlot["MS_Electron2Eta_boxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Electron 2 Eta (box SemiLepZZqqWZqq)", 20, -2.5 , 2.5, "electron 2 Eta");
  //other preboxes
	MSPlot["MS_PtZ_preboxZHqq_el"] = new MultiSamplePlot(datasetsEl,"Pt Z candidate (prebox ZHqq)", 25, 0, 500, "Pt Z candidate (GeV)");
	MSPlot["MS_PtZ_preboxZqq_el"] = new MultiSamplePlot(datasetsEl,"Pt Z candidate (prebox Zqq)", 25, 0, 500, "Pt Z candidate (GeV)");
  MSPlot["MS_PtZ_preboxZHqq_LowerJetCuts_LowerZpt_el"] = new MultiSamplePlot(datasetsEl,"Pt Z candidate (prebox ZHqq)", 15, 0, 250, "Pt Z candidate (GeV)");
  MSPlot["MS_St_preboxZHqq_LowerJetCuts_LowerZpt_el"] = new MultiSamplePlot(datasetsEl,"ST (prebox ZHqq)", 50, 200, 1700, "S_{T} (GeV)");
	MSPlot["MS_MET_preboxZHqq_LowerJetCuts_LowerZpt_el"] = new MultiSamplePlot(datasetsEl,"MET (prebox ZHqq)", 30, 0, 300, "MET (GeV)");
	MSPlot["MS_LeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Pt (prebox ZHqq)", 40, 0, 200, "leading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Pt (prebox ZHqq)", 40, 0, 200, "subleading-jet Pt (GeV)");
	MSPlot["MS_LeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_el"] = new MultiSamplePlot(datasetsEl,"Leading jet Pt (prebox WHqq)", 40, 0, 200, "leading-jet Pt (GeV)");
	MSPlot["MS_SubLeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_el"] = new MultiSamplePlot(datasetsEl,"Subleading jet Pt (prebox WHqq)", 40, 0, 200, "subleading-jet Pt (GeV)");	
	MSPlot["MS_PtH_preboxWHqq_el"] = new MultiSamplePlot(datasetsEl,"Pt H candidate (prebox WHqq)", 25, 0, 500, "Pt H candidate (GeV)");
	MSPlot["MS_PtZ_preboxWZqq_Zwindow_el"] = new MultiSamplePlot(datasetsEl,"Pt Z candidate (prebox WZqq)", 25, 0, 500, "Pt Z candidate (GeV)");
	MSPlot["MS_PtH_preboxWHqq_LowerJetCuts_LowerHpt_el"] = new MultiSamplePlot(datasetsEl,"Pt H candidate (prebox WHqq)", 15, 0, 250, "Pt H candidate (GeV)");
	MSPlot["MS_St_preboxWHqq_LowerJetCuts_LowerHpt_offH_el"] = new MultiSamplePlot(datasetsEl,"ST (prebox WHqq)", 50, 200, 1700, "S_{T} (GeV)");
  MSPlot["MS_MV_preboxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"M(dijet) (prebox SemiLepZZqqWZqq)", 50, 0, 300, "M(dijet) (GeV)");
	MSPlot["MS_PtV_preboxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Pt V candidate (prebox SemiLepZZqqWZqq)", 25, 0, 500, "Pt hadronic W/Z candidate (GeV)");
	MSPlot["MS_PtZ_preboxSemiLepZZqqWZqq_el"] = new MultiSamplePlot(datasetsEl,"Pt Z candidate (prebox SemiLepZZqqWZqq)", 25, 0, 500, "Pt Z candidate (GeV)");
	MSPlot["MS_PtW_preboxWqq_el"] = new MultiSamplePlot(datasetsEl,"Pt W candidate (prebox Wqq)", 25, 0, 500, "Pt W candidate (GeV)");
	MSPlot["MS_MHcandidate_preboxWHqq_el"] = new MultiSamplePlot(datasetsEl,"MH (prebox WHqq)", 50, 0, 500, "reconstructed Higgs mass (GeV)");
	
	
	cout << " - Declared histograms ..." <<  endl;

  //Configuration and variables for 2D distribution.

	
  ////////////////////////////////////
  /// Selection table   --> to be modified
  ////////////////////////////////////
/*
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
*/
  
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
  //temporary fix for fastsim samples produced with the Spring12 PU, which is an S7 scenario, not S10 like in official Summer12 (and higher-mass private fastsim)
	LumiReWeighting LumiWeights_S7_mu = LumiReWeighting(MCpileuphistofile_S7, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleMu/nominal.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsUp_S7_mu = LumiReWeighting(MCpileuphistofile_S7, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleMu/sys_up.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsDown_S7_mu = LumiReWeighting(MCpileuphistofile_S7, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleMu/sys_down.root", "pileup", "pileup");
 
	//For Run2012A+B_C_D (rereco) SingleElectron
  LumiReWeighting LumiWeights_el = LumiReWeighting(MCpileuphistofile, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleElectron/nominal.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsUp_el = LumiReWeighting(MCpileuphistofile, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleElectron/sys_up.root", "pileup", "pileup");
  LumiReWeighting LumiWeightsDown_el = LumiReWeighting(MCpileuphistofile, "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/Run2012ABCDrereco_SingleElectron/sys_down.root", "pileup", "pileup");
  //temporary fix for fastsim samples produced with the Spring12 PU, which is an S7 scenario, not S10 like in official Summer12 (and higher-mass private fastsim)
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
		 bTool_T = new BTagWeightTools("../../TopTreeAnalysisBase/Calibrations/BTagging/SFb-pt_WITHttbar_payload_EPS13_noCSVSL.txt",btagger,"EPS2013");  //for b-tagging for H->bbar with application of Sfs and created MC eff histos... never needed without antib-tagging? well, in ttbar control regions...?
	   bTool_T->ReadMCEfficiencyHistos("PlotsForBtagWeights_CSVM_ttbar.root");
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
			
		if(makeflatbinning && (dataSetName.find("Data")==0 || dataSetName.find("NP")==0)) //dataSetName.find("InvIso")!=string::npos
			continue; //make binning only using SM MC + data-driven QCD...
			
		if(systematic != "Nominal" && (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0))
		{
		  cout<<"  Running systematics: skipping data"<<endl;  //euhm, why is this? should be possible to plot data on top of systematic shifted samples... TO DO
		  continue;
		}	
			
    if (verbose > 1)
      cout << "file: " << inputTrees[d] << endl;
	 
    TFile* inFile = new TFile(inputTrees[d].c_str(),"READ");
		
		
		//very dirty hack for W+HF/LF splitting
		if(dataSetName.find("WJets_LF_4jets") == 0 || dataSetName.find("WJets_HF_4jets") == 0) inFile = new TFile(TString(Treespath)+"VLQTree_WJets_4jets.root","READ");
		else if(dataSetName.find("WJets_LF_3jets") == 0 || dataSetName.find("WJets_HF_3jets") == 0) inFile = new TFile(TString(Treespath)+"VLQTree_WJets_3jets.root","READ");
		else if(dataSetName.find("WJets_LF_2jets") == 0 || dataSetName.find("WJets_HF_2jets") == 0) inFile = new TFile(TString(Treespath)+"VLQTree_WJets_2jets.root","READ");
		else if(dataSetName.find("WJets_LF_1jets") == 0 || dataSetName.find("WJets_HF_1jets") == 0) inFile = new TFile(TString(Treespath)+"VLQTree_WJets_1jets.root","READ");
		//very dirty hack for Z+HF/LF splitting
		if(dataSetName.find("ZJets_LF_4jets") == 0 || dataSetName.find("ZJets_HF_4jets") == 0) inFile = new TFile(TString(Treespath)+"VLQTree_ZJets_4jets.root","READ");
		else if(dataSetName.find("ZJets_LF_3jets") == 0 || dataSetName.find("ZJets_HF_3jets") == 0) inFile = new TFile(TString(Treespath)+"VLQTree_ZJets_3jets.root","READ");
		else if(dataSetName.find("ZJets_LF_2jets") == 0 || dataSetName.find("ZJets_HF_2jets") == 0) inFile = new TFile(TString(Treespath)+"VLQTree_ZJets_2jets.root","READ");
		else if(dataSetName.find("ZJets_LF_1jets") == 0 || dataSetName.find("ZJets_HF_1jets") == 0) inFile = new TFile(TString(Treespath)+"VLQTree_ZJets_1jets.root","READ");
		
		
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
		//for (int ievt = start; ievt < 2000; ievt++)
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
			bool MuChannel = myBranch_selectedEvents->MuChannel();
			bool ElChannel = myBranch_selectedEvents->ElChannel();
      unsigned int SelectednLeptons = myBranch_selectedEvents->SelectednLeptons();
      unsigned int SelectednMu = myBranch_selectedEvents->SelectednMu();
      unsigned int SelectednEl = myBranch_selectedEvents->SelectednEl();			
						
			TLorentzVector metvector = myBranch_selectedEvents->met();
			float met = (myBranch_selectedEvents->met()).Pt();  //Pt instead of Et, because the Et is not corrected when applying JEC or systematic shifts. update: it seems it is corrected, but not for unclustered MET?
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
						//cout<<"     found antitagged jet: pushing back "<<iJet<<endl;
						selectedAntiBtaggedJets_indices.push_back(iJet);
				 }
				 else  
				    selectedOthertaggedJets.push_back(selectedJets[iJet]); //also non-taggable jets? (discr value -1 or -10)				 
			}
			
			//for 2/3/4-lepton categories
			vector<unsigned int> indicesMuPlus, indicesMuMinus, indicesElPlus, indicesElMinus;
			for(unsigned int iMu=0; iMu<selectedMuonsCharge.size(); iMu++)
			{
			  if(selectedMuonsCharge[iMu]>0) indicesMuPlus.push_back(iMu);
				else if (selectedMuonsCharge[iMu]<0) indicesMuMinus.push_back(iMu);
			}
			for(unsigned int iEl=0; iEl<selectedElectronsCharge.size(); iEl++)
			{
			  if(selectedElectronsCharge[iEl]>0) indicesElPlus.push_back(iEl);
				else if (selectedElectronsCharge[iEl]<0) indicesElMinus.push_back(iEl);
			}



      //W+LF and W+HF, Z+LF and Z+HF splitting (dirty hack)
			bool WHF = false; //HF = heavy flavor, bb and cc (just definition)
			if(dataSetName.find("WJets") == 0)
			{
        for(int iJet = 0; iJet < selectedJets_partonFlavour.size(); iJet++)
			  {
			     if(fabs(selectedJets_partonFlavour[iJet])==5 || fabs(selectedJets_partonFlavour[iJet])==4)
				   {
						  WHF = true;
						  break;
				   } 			
			  }
      }
      bool ZHF = false; //HF = heavy flavor, bb and cc (just definition)
			if(dataSetName.find("ZJets") == 0)
			{
        for(int iJet = 0; iJet < selectedJets_partonFlavour.size(); iJet++)
			  {
			     if(fabs(selectedJets_partonFlavour[iJet])==5 || fabs(selectedJets_partonFlavour[iJet])==4)
				   {
						  ZHF = true;
						  break;
				   } 			
			  }
      }			
			if((dataSetName.find("WJets_LF") == 0 && WHF) || (dataSetName.find("WJets_HF") == 0 && !WHF))
			   continue; //skip event
      if((dataSetName.find("ZJets_LF") == 0 && ZHF) || (dataSetName.find("ZJets_HF") == 0 && !ZHF))
			   continue; //skip event


	    //cout<<"LOADED event with run = "<<myBranch_selectedEvents->runID()<<", lumiblock = "<<myBranch_selectedEvents->lumiBlockID()<<", event id = "<<myBranch_selectedEvents->eventID()<<endl;
	    //cout<<"    => MuTrigged = "<<MuTrigged<<", ElTrigged = "<<ElTrigged<<", SelectednEl = "<<SelectednEl<<", SelectednMu = "<<SelectednMu<<", selectedJets[0].Pt() = "<<selectedJets[0].Pt()<<endl;
	 
	 
      // scale factor for the event
      float scaleFactor = 1.;
      float lumiWeight = 1.;			
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("InvIso") != string::npos) )
			{			
			  bool doS10 = true;
				if(dataSetName.find("NP_") == 0)
				{
					    if(dataSetName.find("QD")== 0 && (dataSetName.find("mass1000")== 0 || dataSetName.find("mass1050")== 0 || dataSetName.find("mass1100")== 0 || dataSetName.find("mass1150")== 0 || dataSetName.find("mass1200")== 0 || dataSetName.find("mass1250")== 0 || dataSetName.find("mass1300")== 0 || dataSetName.find("mass1350")== 0 || dataSetName.find("mass1400")== 0 || dataSetName.find("mass1450")== 0 || dataSetName.find("mass1500")== 0 || dataSetName.find("mass1550")== 0 || dataSetName.find("mass1600")== 0 || dataSetName.find("mass1650")== 0 || dataSetName.find("mass1700")== 0 || dataSetName.find("mass1750")== 0 || dataSetName.find("mass1800")== 0))
							{
							   doS10 = true;
							}
							else if(dataSetName.find("QU")== 0 && (dataSetName.find("mass950")== 0 ||dataSetName.find("mass1000")== 0 || dataSetName.find("mass1050")== 0 || dataSetName.find("mass1100")== 0 || dataSetName.find("mass1150")== 0 || dataSetName.find("mass1200")== 0 || dataSetName.find("mass1250")== 0 || dataSetName.find("mass1300")== 0 || dataSetName.find("mass1350")== 0 || dataSetName.find("mass1400")== 0 || dataSetName.find("mass1450")== 0 || dataSetName.find("mass1500")== 0 || dataSetName.find("mass1550")== 0 || dataSetName.find("mass1600")== 0 || dataSetName.find("mass1650")== 0 || dataSetName.find("mass1700")== 0 || dataSetName.find("mass1750")== 0 || dataSetName.find("mass1800")== 0))
							{
							   doS10 = true;
							}
							else if (dataSetName.find("fullsim") != string::npos || dataSetName.find("NP_TprimeTprime") == 0)
							{
							   doS10 = true;							
							}
							else
							{
							   doS10 = false;
							}
				}		
			
			  //if(MuTrigged)
				if(MuChannel)
				{											
				  //if ( dataSetName.find("NP_") == 0 && dataSetName.find("fullsim") == string::npos && dataSetName.find("NP_TprimeTprime") == string::npos)
					if(!doS10)
					{
					   //cout<<"Identified PU scenario S7"<<endl;
             if(systematic == "PUMinus") lumiWeight = LumiWeightsDown_S7_mu.ITweight( (int) nTruePU );
						 else if(systematic == "PUPlus") lumiWeight = LumiWeightsUp_S7_mu.ITweight( (int) nTruePU );
						 else lumiWeight = LumiWeights_S7_mu.ITweight( (int) nTruePU );
					}
					else
					{			
					   //cout<<"Identified PU scenario S10"<<endl;
					   if(systematic == "PUMinus") lumiWeight = LumiWeightsDown_mu.ITweight( (int) nTruePU );
						 else if(systematic == "PUPlus") lumiWeight = LumiWeightsUp_mu.ITweight( (int) nTruePU );
						 else lumiWeight = LumiWeights_mu.ITweight( (int) nTruePU );
					}
				}
				//else if(ElTrigged)
				else if(ElChannel)
				{
				  //if ( dataSetName.find("NP_") == 0 && dataSetName.find("fullsim") == string::npos && dataSetName.find("NP_TprimeTprime") == string::npos)
					if(!doS10)
					{
					   cout<<"Identified PU scenario S7"<<endl;
					   if(systematic == "PUMinus") lumiWeight = LumiWeightsDown_S7_el.ITweight( (int) nTruePU );
						 else if(systematic == "PUPlus") lumiWeight = LumiWeightsUp_S7_el.ITweight( (int) nTruePU );
						 else lumiWeight = LumiWeights_S7_el.ITweight( (int) nTruePU );
					}
					else
					{
					   cout<<"Identified PU scenario S10"<<endl;
					   if(systematic == "PUMinus") lumiWeight = LumiWeightsDown_el.ITweight( (int) nTruePU );
						 else if(systematic == "PUPlus") lumiWeight = LumiWeightsUp_el.ITweight( (int) nTruePU );
						 else lumiWeight = LumiWeights_el.ITweight( (int) nTruePU );
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

      if(doPUreweighting)
			{
			   scaleFactor = scaleFactor*lumiWeight;
         histo1D["lumiWeights"]->Fill(scaleFactor);
			}
			
			
			if(dataSetName.find("TTbarJets") == 0 && applyPtTopReweighting)
			{
			   float pttopweight = myBranch_selectedEvents->topPtReWeight();	
				 scaleFactor = scaleFactor*pttopweight;
			}


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

      if(doWjetsScaling && dataSetName.find("WJets_LF")==0) scaleFactor = scaleFactor*WjetsScale;
      if(doZjetsScaling && dataSetName.find("ZJets_LF")==0) scaleFactor = scaleFactor*ZjetsScale;
      if(doInvIsoScaling && dataSetName.find("InvIso_Mu")==0) scaleFactor = scaleFactor*InvIsomuScale;
      if(doInvIsoScaling && dataSetName.find("InvIso_El")==0) scaleFactor = scaleFactor*InvIsoelScale;
		 
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
			
			
      float Zmass = 91.1876; //could eventually be replaced by the mean of a fit of the Z peak?
			float Zmasswindow = 7.5; //could eventually be replaced by the resolution of a fit of the Z peak (1 sigma)?


      //'muon channel'
      //if(SelectednMu >= 1) //not anymore
			//if(MuTrigged && selectedJets.size()>=1)
			if(MuChannel && selectedJets.size()>=1)
			{
			  //cout<<"**** Muon channel **** "<<endl;
			  if(!useBtagInfoFromFile)
				{
				   //to determine MC b-tag efficiency. Double-check if right place to put it
				   bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
        	 bTool_T->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);  
				}
				
									
			  //lepton scale factor; apply on MC --> WARNING: this might not be entirely correct, it should be applied on all the leptons that are part of your selection, not all leptons in the event!! -> to be adapted
				if(!((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			  {								
				   if(applyLeptonSF)
					 {
					   //cout<<"(selectedMuons[0].Eta(), selectedMuons[0].Pt()) = ("<<selectedMuons[0].Eta()<<","<< selectedMuons[0].Pt()<<")"<<endl;
					   scaleFactor = scaleFactor*leptonTools->getMuonSF(selectedMuons[0].Eta(), selectedMuons[0].Pt(),syst_muonSF); //note: this is assuming the leading muon is the trigged one...
					   //cout<<"leptonTools->getMuonSF(selectedMuons[0].Eta(), selectedMuons[0].Pt()) = "<<leptonTools->getMuonSF(selectedMuons[0].Eta(), selectedMuons[0].Pt(),syst_muonSF)<<endl;
						 for(unsigned int iMu = 0; iMu < selectedMuons.size(); iMu++)
			       {
					     if(iMu!=0 && (selectedMuons[iMu].Pt()>=20))
							 {
							    scaleFactor = scaleFactor*leptonTools->getMuonIdIsoSF(selectedMuons[iMu].Eta(), selectedMuons[iMu].Pt(),syst_muonSF); //apply id*iso scale factor for non-leading muons.
							    //cout<<"   leptonTools->getMuonIdIsoSF(selectedMuons["<<iMu<<"].Eta(), selectedMuons["<<iMu<<"].Pt()) = "<<leptonTools->getMuonIdIsoSF(selectedMuons[iMu].Eta(), selectedMuons[iMu].Pt(),syst_muonSF)<<endl;
							 }
					   }
					   for(unsigned int iEl = 0; iEl < selectedElectrons.size(); iEl++)
			       {
					     if(selectedElectrons[iEl].Pt()>=20) 
							 {
							   scaleFactor = scaleFactor*leptonTools->getElectronIdIsoSF(selectedElectrons[iEl].Eta(), selectedElectrons[iEl].Pt(),syst_electronSF); //apply idiso scale factor for electrons
					       //cout<<"   leptonTools->getElectronIdIsoSF(selectedElectrons["<<iEl<<"].Eta(), selectedElectrons["<<iEl<<"].Pt()) = "<<leptonTools->getElectronIdIsoSF(selectedElectrons[iEl].Eta(), selectedElectrons[iEl].Pt(),syst_electronSF)<<endl;
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
							if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							{
							      float myBtagMCWeight = 1;
										if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							            myBtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> myBtagMCWeight = "<<myBtagMCWeight<<endl; 
							      }
										MSPlot["MS_MET_singleMu_Bveto_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
							}
							MSPlot["MS_JetMultiplicity_singleMu_mu"]->Fill(selectedJets.size(),datasets[d], true, LuminosityMu*scaleFactor);
					 }
					 if(selectedJets.size()>=4)
					 { 
					    //pt-top reweighting from Michael's https://github.com/TopBrussels/TopTreeAnalysis/blob/CMSSW_53X/macros/BTagTReeCreator.cc:
							/*for (int p=0; p<mcParticles.size(); p++)
							{        
                if (mcParticles[p]->status() != 3 ) continue;
                if (mcParticles[p]->Pt() == 0 ) continue;                 
                double pt = mcParticles[p]->Pt();
                //if (mcParticles[p]->type() == 6) w_pttop = (1.4*pow(10.0,-6.0)*pt*pt)-(2.0*pow(10.0,-3.0)*pt)+1.2;        
                if (mcParticles[p]->type() == 6) w_pttop_top = exp(0.156-(0.00137*pt));
                if (mcParticles[p]->type() == -6) w_pttop_atop = exp(0.156-(0.00137*pt));
              }*/
							
							//if(met > METcutboxWqq && selectedJets[0].Pt() >= LeadingJetPtCut*0.5)
							if(met > METcutboxWqq)
							{
					     if(selectedBtaggedJets.size()>=1)
							 {
							   float myBtagMCWeight = 1;
							   if(applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			           {
							            myBtagMCWeight = bTool_T->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> myBtagMCWeight = "<<myBtagMCWeight<<endl; 
							   }
							   MSPlot["MS_MET_singleMu_ttbarsel1btag_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
								 MSPlot["MS_MuonPt_singleMu_ttbarsel1btag_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
								 MSPlot["MS_MuonEta_singleMu_ttbarsel1btag_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
							   MSPlot["MS_Pt_jet1_singleMu_ttbarsel1btag_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
								 MSPlot["MS_Pt_jet2_singleMu_ttbarsel1btag_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
								 MSPlot["MS_Eta_jet1_singleMu_ttbarsel1btag_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
								 MSPlot["MS_Eta_jet2_singleMu_ttbarsel1btag_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
								 MSPlot["MS_HT_singleMu_ttbarsel1btag_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);
								 //with pt-top reweighting -> outdated, I'm applying it anyway
								 /*MSPlot["MS_MET_singleMu_ttbarsel1btag_pttopreweighting_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
								 MSPlot["MS_MuonPt_singleMu_ttbarsel1btag_pttopreweighting_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								 MSPlot["MS_MuonEta_singleMu_ttbarsel1btag_pttopreweighting_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							   MSPlot["MS_Pt_jet1_singleMu_ttbarsel1btag_pttopreweighting_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								 MSPlot["MS_Pt_jet2_singleMu_ttbarsel1btag_pttopreweighting_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								 MSPlot["MS_Eta_jet1_singleMu_ttbarsel1btag_pttopreweighting_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
								 MSPlot["MS_Eta_jet2_singleMu_ttbarsel1btag_pttopreweighting_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
								 MSPlot["MS_HT_singleMu_ttbarsel1btag_pttopreweighting_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);*/
								 
							 }							
               if(selectedBtaggedJets.size()>=2)
							 {
							   float myBtagMCWeight = 1;
							   if(applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			           {
							            myBtagMCWeight = bTool_T->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> myBtagMCWeight = "<<myBtagMCWeight<<endl; 
							   }
							   MSPlot["MS_MET_singleMu_ttbarsel2btag_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
								 MSPlot["MS_MuonPt_singleMu_ttbarsel2btag_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
								 MSPlot["MS_MuonEta_singleMu_ttbarsel2btag_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
							   MSPlot["MS_Pt_jet1_singleMu_ttbarsel2btag_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								 MSPlot["MS_Pt_jet2_singleMu_ttbarsel2btag_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								 MSPlot["MS_Eta_jet1_singleMu_ttbarsel2btag_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
								 MSPlot["MS_Eta_jet2_singleMu_ttbarsel2btag_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
								 MSPlot["MS_HT_singleMu_ttbarsel2btag_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);
							 }
							}
					 }					 
				}
				
				
				//lepton control plots...
				if(SelectednMu == 2 && SelectednEl == 0)
				{
					 if(selectedJets.size()>=4)
					 {
					    //if(fabs((selectedMuons[0] + selectedMuons[1]).M() - Zmass) < Zmasswindow)
							//{
							MSPlot["MS_Muon1Pt_dimuon_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_Muon1Eta_dimuon_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_Muon2Pt_dimuon_mu"]->Fill(selectedMuons[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_Muon2Eta_dimuon_mu"]->Fill(selectedMuons[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							//}					 
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
			//else if(ElTrigged && selectedJets.size()>=1)  //or an 'if'?? was else if but changed back to if; souble counting was taken care of with the pt requirements on the leptons
			else if(ElChannel && selectedJets.size()>=1)
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
				     scaleFactor = scaleFactor*leptonTools->getElectronSF(selectedElectrons[0].Eta(), selectedElectrons[0].Pt(),syst_electronSF); //note: this is assuming the leading electron is the trigged one...
					   //cout<<"leptonTools->getElectronSF(selectedElectrons[0].Eta(), selectedElectrons[0].Pt()) = "<<leptonTools->getElectronSF(selectedElectrons[0].Eta(), selectedElectrons[0].Pt(),syst_electronSF)<<endl;
						 for(unsigned int iMu = 0; iMu < selectedMuons.size(); iMu++)
			       {
					     if(selectedMuons[iMu].Pt()>=20)
							 {
							    scaleFactor = scaleFactor*leptonTools->getMuonIdIsoSF(selectedMuons[iMu].Eta(), selectedMuons[iMu].Pt(),syst_muonSF); //apply id*iso scale factor for muons.
					   			//cout<<"  leptonTools->getMuonIdIsoSF(selectedMuons["<<iMu<<"].Eta(), selectedMuons["<<iMu<<"].Pt()) = "<<leptonTools->getMuonIdIsoSF(selectedMuons[iMu].Eta(), selectedMuons[iMu].Pt(),syst_muonSF)<<endl;
               }
             }
					   for(unsigned int iEl = 0; iEl < selectedElectrons.size(); iEl++)
			       {
					     if(iEl!=0 && (selectedElectrons[iEl].Pt()>=20))
							 {
							    scaleFactor = scaleFactor*leptonTools->getElectronIdIsoSF(selectedElectrons[iEl].Eta(), selectedElectrons[iEl].Pt(),syst_electronSF); //apply idiso scale factor for non-leading electrons.
               		//cout<<"  leptonTools->getElectronIdIsoSF(selectedElectrons["<<iEl<<"].Eta(), selectedElectrons["<<iEl<<"].Pt()) = "<<leptonTools->getElectronIdIsoSF(selectedElectrons[iEl].Eta(), selectedElectrons[iEl].Pt(),syst_electronSF)<<endl;
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
								//cout<<"   run = "<<myBranch_selectedEvents->runID()<<", lumiblock = "<<myBranch_selectedEvents->lumiBlockID()<<", event id = "<<myBranch_selectedEvents->eventID()<<endl;
				        MSPlot["MS_MET_singleEl_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
								if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							  {
							      float myBtagMCWeight = 1;
										if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							            myBtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> myBtagMCWeight = "<<myBtagMCWeight<<endl; 
							      }
										MSPlot["MS_MET_singleEl_Bveto_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor*myBtagMCWeight);
							  }
								MSPlot["MS_JetMultiplicity_singleEl_el"]->Fill(selectedJets.size(),datasets[d], true, LuminosityEl*scaleFactor);
					 }
				}
				
				
				//lepton control plots...
				if(SelectednEl == 2 && SelectednMu == 0)
				{
					 if(selectedJets.size()>=4)
					 {
					    //if(fabs((selectedElectrons[0] + selectedElectrons[1]).M() - Zmass) < Zmasswindow)
							//{
							MSPlot["MS_Electron1Pt_dielectron_el"]->Fill(selectedElectrons[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							MSPlot["MS_Electron1Eta_dielectron_el"]->Fill(selectedElectrons[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
							MSPlot["MS_Electron2Pt_dielectron_el"]->Fill(selectedElectrons[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							MSPlot["MS_Electron2Eta_dielectron_el"]->Fill(selectedElectrons[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
							//}					 
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
			float BtagMCWeight = 1.;
			
					
			//define the 'boxes'
			//'muon channel'
			//if(MuTrigged)
			if(MuChannel)
			{
//			 selecTableMu.Fill(d,0,scaleFactor);
			 if(SelectednMu == 1)
			 {
			  if(SelectednEl == 0)
				{				
				  //sort of control region
				  if(selectedJets.size() >= 1)
					{
					    if((selectedJets[0].Pt() >= 150) && (selectedJets[0].Pt() <= 200))
						  {
							   float MT = sqrt(2*selectedMuons[0].Pt()*met*(1-cos(selectedMuons[0].DeltaPhi(metvector))));
						     if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							   {	
								    if((!doMETcutboxWqq || (doMETcutboxWqq && met>METcutboxWqq)) && (!doMTcut || (doMTcut && MT>MTcut)))
							      {
										   float myBtagMCWeight = 1;
										   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                 {
							            myBtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
							         }	
											 //cout<<"completely inside: selectedJets.size() = "<<selectedJets.size()<<endl;
											 if(selectedForwardJets.size() == 0) MSPlot["MS_nJets_preboxWqq_leadingjet150-200_noForwardJet_noNumberCentralJetCut_mu"]->Fill(selectedJets.size(), datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);	 		 //when requiring selectedForwardJets.size() == 0 you will only select events with njets>=2 because of the skimming preselection							 
										   else MSPlot["MS_nJets_preboxWqq_leadingjet150-200_noNumberCentralJetCut_mu"]->Fill(selectedJets.size(), datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);										 
										}
								 }
						  }
				  }
				
			    if( (selectedJets.size() == 1 || selectedJets.size() == 2))
			    {
					 if(selectedForwardJets.size() == 1)
					 {
					  ST = ST + selectedForwardJets[0].Pt();
						HT = HT + selectedForwardJets[0].Pt();
					  MSPlot["MS_MET_preboxWqq_noleadingjetcut_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
					  MSPlot["MS_HT_preboxWqq_noleadingjetcut_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);
						MSPlot["MS_LeadingJetPt_preboxWqq_noleadingjetcut_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
						MSPlot["MS_LeadingJetEta_preboxWqq_noleadingjetcut_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
						if(selectedJets.size()>1)
						{
								  MSPlot["MS_SubLeadingJetPt_preboxWqq_noleadingjetcut_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								  MSPlot["MS_SubLeadingJetEta_preboxWqq_noleadingjetcut_mu"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
						}
						
						//some sort of 'control region' (to be checked if not too much signal in it, though)
						if((selectedJets[0].Pt() >= 150) && (selectedJets[0].Pt() <= 200))
						{
						 if(!doMETcutboxWqq || (doMETcutboxWqq && met>METcutboxWqq))
						 {
						  MSPlot["MS_MET_preboxWqq_leadingjet150-200_METcut_mu"]->Fill(met, datasets[d], true, LuminosityMu*scaleFactor);
						  MSPlot["MS_HT_preboxWqq_leadingjet150-200_METcut_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_St_preboxWqq_leadingjet150-200_METcut_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
						  MSPlot["MS_LeadingJetPt_preboxWqq_leadingjet150-200_METcut_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_LeadingJetEta_preboxWqq_leadingjet150-200_METcut_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							if(selectedJets.size()>1)
							{
							   MSPlot["MS_SubleadingJetPt_preboxWqq_leadingjet150-200_METcut_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							   MSPlot["MS_SubleadingJetEta_preboxWqq_leadingjet150-200_METcut_mu"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							}
							MSPlot["MS_MuonPt_preboxWqq_leadingjet150-200_METcut_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_MuonEta_preboxWqq_leadingjet150-200_METcut_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
						 }
						 MSPlot["MS_MET_preboxWqq_leadingjet150-200_mu"]->Fill(met, datasets[d], true, LuminosityMu*scaleFactor);
						}
						
						
						if(selectedJets[0].Pt() >= LeadingJetPtCut)
						{
						  float MT = sqrt(2*selectedMuons[0].Pt()*met*(1-cos(selectedMuons[0].DeltaPhi(metvector))));
							MSPlot["MS_MTleptonmet_preboxWqq_mu"]->Fill(MT, datasets[d], true, LuminosityMu*scaleFactor);	
							MSPlot["MS_MET_preboxWqq_noMETMTcut_mu"]->Fill(met, datasets[d], true, LuminosityMu*scaleFactor);	
							MSPlot["MS_HT_preboxWqq_noMETMTcut_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_LeadingJetPt_preboxWqq_noMETMTcut_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_LeadingJetEta_preboxWqq_noMETMTcut_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							if(selectedJets.size()>1)
							{
								  MSPlot["MS_SubLeadingJetPt_preboxWqq_noMETMTcut_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								  MSPlot["MS_SubLeadingJetEta_preboxWqq_noMETMTcut_mu"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							}
							MSPlot["MS_ForwardJetPt_preboxWqq_noMETMTcut_mu"]->Fill(selectedForwardJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_ForwardJetEta_preboxWqq_noMETMTcut_mu"]->Fill(selectedForwardJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_nJets_preboxWqq_noMETMTcut_mu"]->Fill(selectedJets.size(),datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_MuonPt_preboxWqq_noMETMTcut_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							MSPlot["MS_MuonEta_preboxWqq_noMETMTcut_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							if((!doMETcutboxWqq || (doMETcutboxWqq && met>METcutboxWqq)) && (!doMTcut || (doMTcut && MT>MTcut)))
							{		
							 MSPlot["MS_MET_preboxWqq_noBveto_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
							 //cout<<"selectedJets.size() = "<<selectedJets.size()<<", selectedAntiBtaggedJets.size() = "<<selectedAntiBtaggedJets.size()<<", selectedOthertaggedJets.size() = "<<selectedOthertaggedJets.size()<<", selectedBtaggedJets.size() = "<<selectedBtaggedJets.size()<<endl;
							 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							 {							 
							     if(verbosityBox) cout<<"category: Wqq (mu)"<<endl;							
							     //apply b-tag SF on MC; taken from LightStopSearch.cc and https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
	                 if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			             {
							        BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								      //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								      MSPlot["MS_BtagWeight_boxWqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								      scaleFactor = 	scaleFactor*BtagMCWeight; 
							     }
									 
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
								   TLorentzVector Wcandidate = recneutrino + selectedMuons[0];
									 float PtW = Wcandidate.Pt();
								   
									 MSPlot["MS_PtW_preboxWqq_mu"]->Fill(PtW,datasets[d], true, LuminosityMu*scaleFactor);
									 if((doWptcut && Wcandidate.Pt() > Wptcut) || !doWptcut)
									 {
			 		           MSPlot["MS_MET_boxWqq_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
								     MSPlot["MS_MTleptonmet_boxWqq_mu"]->Fill(MT, datasets[d], true, LuminosityMu*scaleFactor);	
							       MSPlot["MS_St_boxWqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
								     MSPlot["MS_HT_boxWqq_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);
								     MSPlot["MS_LeadingJetPt_boxWqq_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								     MSPlot["MS_LeadingJetEta_boxWqq_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
								     if(selectedJets.size()>1)
								     {
								        MSPlot["MS_SubLeadingJetPt_boxWqq_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
									  	  MSPlot["MS_SubLeadingJetEta_boxWqq_mu"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
								     }
								     MSPlot["MS_ForwardJetPt_boxWqq_mu"]->Fill(selectedForwardJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								     MSPlot["MS_ForwardJetEta_boxWqq_mu"]->Fill(selectedForwardJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
								     MSPlot["MS_nJets_boxWqq_mu"]->Fill(selectedJets.size(),datasets[d], true, LuminosityMu*scaleFactor);
								     MSPlot["MS_MuonPt_boxWqq_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								     MSPlot["MS_MuonEta_boxWqq_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
								     MSPlot["MS_MCeventWeight_boxWqq_mu"]->Fill(LuminosityMu*scaleFactor*datasets[d]->NormFactor(),datasets[d], false, LuminosityMu*scaleFactor); //not scaling the events...								 								 
							       float MQ = (Wcandidate + selectedJets[0]).M();
										 MSPlot["MS_MQ_boxWqq_mu"]->Fill(MQ,datasets[d], true, LuminosityMu*scaleFactor);
								     
										 
										 //for 'flat' binning
										 float fillweight = datasets[d]->NormFactor () * LuminosityMu * scaleFactor;
										 //cout<<"  ---> d = "<<d<<endl;
										 //cout<<"  ---> datasets["<<d<<"] = "<<datasets[d]->Name()<<endl;
										 if(makeflatbinning) VarWqq_mu.Fill_forBinning(MQ,fillweight);
					           else if(useflatbinning) VarWqq_mu.Fill(MQ,fillweight,dataSetName);
						 
								     //other plots
								     float DeltaR_WcandLeadingJet = selectedJets[0].DeltaR(Wcandidate);
								     float DeltaEta_WcandLeadingJet = fabs(selectedJets[0].Eta() - Wcandidate.Eta());
								     float DeltaR_LeptonNeutrino = selectedMuons[0].DeltaR(recneutrino);
								     MSPlot["MS_DeltaR_WcandLeadingJet_boxWqq_mu"]->Fill(DeltaR_WcandLeadingJet,datasets[d], true, LuminosityMu*scaleFactor);								 
								     MSPlot["MS_DeltaEta_WcandLeadingJet_boxWqq_mu"]->Fill(DeltaEta_WcandLeadingJet,datasets[d], true, LuminosityMu*scaleFactor);
								     MSPlot["MS_DeltaR_LeptonNeutrino_boxWqq_mu"]->Fill(DeltaR_LeptonNeutrino,datasets[d], true, LuminosityMu*scaleFactor);
								 
								     if(selectedMuonsCharge[0] == 1)
								     {
								        if(verbosityBox) cout<<"          W+ qq (mu)"<<endl;
								        MSPlot["MS_St_boxWplusqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
									    	MSPlot["MS_MTQ_boxWplusqq_mu"]->Fill(MTQ,datasets[d], true, LuminosityMu*scaleFactor);
									    	MSPlot["MS_MQ_boxWplusqq_mu"]->Fill(MQ,datasets[d], true, LuminosityMu*scaleFactor);
												if(makeflatbinning) VarWplusqq_mu.Fill_forBinning(MQ,fillweight);
												else if(useflatbinning) VarWplusqq_mu.Fill(MQ,fillweight,dataSetName);
								     }
								     else if(selectedMuonsCharge[0] == -1)
								     {
								        if(verbosityBox) cout<<"          W- qq (mu)"<<endl;
								        MSPlot["MS_St_boxWminusqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
									    	MSPlot["MS_MTQ_boxWminusqq_mu"]->Fill(MTQ,datasets[d], true, LuminosityMu*scaleFactor);
									  	  MSPlot["MS_MQ_boxWminusqq_mu"]->Fill(MQ,datasets[d], true, LuminosityMu*scaleFactor);
												if(makeflatbinning) VarWminusqq_mu.Fill_forBinning(MQ,fillweight);
												else if(useflatbinning) VarWminusqq_mu.Fill(MQ,fillweight,dataSetName);
								     }
								     else cout<<"WARNING: muon charge "<<selectedMuonsCharge[1]<<" does not make sense!"<<endl;	
									 }
							 }	 								 
							}
						}
					 }
					 else if(selectedForwardJets.size() == 0)
					 {
					    //sort of control region
							if((selectedJets[0].Pt() >= 150) && (selectedJets[0].Pt() <= 200))
						  {
							   float MT = sqrt(2*selectedMuons[0].Pt()*met*(1-cos(selectedMuons[0].DeltaPhi(metvector))));
						     if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							   {	
								    MSPlot["MS_MET_preboxWqq_leadingjet150-200_noForwardJet_noMETMTcut_mu"]->Fill(met, datasets[d], true, LuminosityMu*scaleFactor);
								    if((!doMETcutboxWqq || (doMETcutboxWqq && met>METcutboxWqq)) && (!doMTcut || (doMTcut && MT>MTcut)))
							      {
										   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                 {
							            BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								          scaleFactor = 	scaleFactor*BtagMCWeight; 
							         }
											 MSPlot["MS_MET_preboxWqq_leadingjet150-200_noForwardJet_mu"]->Fill(met, datasets[d], true, LuminosityMu*scaleFactor);
											 MSPlot["MS_MTleptonmet_preboxWqq_leadingjet150-200_noForwardJet_mu"]->Fill(MT, datasets[d], true, LuminosityMu*scaleFactor);	
						           MSPlot["MS_HT_preboxWqq_leadingjet150-200_noForwardJet_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);
							         MSPlot["MS_St_preboxWqq_leadingjet150-200_noForwardJet_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
						           MSPlot["MS_LeadingJetPt_preboxWqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							         MSPlot["MS_LeadingJetEta_preboxWqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							         if(selectedJets.size()>1)
							         {
							           MSPlot["MS_SubleadingJetPt_preboxWqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
						        	   MSPlot["MS_SubleadingJetEta_preboxWqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
					          	 }
					          	 MSPlot["MS_MuonPt_preboxWqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							         MSPlot["MS_MuonEta_preboxWqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
										}
								 }
						  }					 
					 }				
				  }
				  else if((doAntiBtagging && selectedAntiBtaggedJets.size()>=2 && selectedJets.size() >= 3) || (!doAntiBtagging && selectedJets.size() >= 3))
			    {
					  if((doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() >= LeadingJetPtCut) && (selectedAntiBtaggedJets[1].Pt() >= SubLeadingJetPtCut)) || (!doAntiBtagging && (selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))) //used to be just if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
						{														
						  int nBtags = selectedBtaggedJets.size();														
							if(!doBtagging || (doBtagging && nBtags>=1))
							{
							   if((doMETcutboxWHqq && met > METcutboxWHqq) || !doMETcutboxWHqq)
								 {			
								    if(verbosityBox) cout<<"category: WHqq (mu)"<<endl;
							      if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							           BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								         //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
											   MSPlot["MS_BtagWeight_boxWHqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								         scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }
										
										vector<pair<int,float> > selectedAvailableJetsindices_sortedbTagCSV;	//available: everything except the two highest-Pt anti-tagged jets, because I'm going to assume they come from the direct decay of the VLQ							    
										for(int iJet = 0; iJet < selectedJets.size(); iJet++)
			              {
										  if(doAntiBtagging)
										  {
										    if(selectedAntiBtaggedJets_indices[0] != iJet && selectedAntiBtaggedJets_indices[1] != iJet) //meaning, if not one of the two leading anti-btagged jets
											  {
											    pair<int,float> jetindex_btagCSV = make_pair(iJet,selectedJets_bTagCSV[iJet]);	
											    selectedAvailableJetsindices_sortedbTagCSV.push_back(jetindex_btagCSV);
											  }
										  }
											else
											{
											  if(0 != iJet && 1 != iJet) //meaning, if not one of the two leading jets
											  {
											    pair<int,float> jetindex_btagCSV = make_pair(iJet,selectedJets_bTagCSV[iJet]);	
											    selectedAvailableJetsindices_sortedbTagCSV.push_back(jetindex_btagCSV);
											  }											
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
										 //mHcandidate = Hcandidate.M();
										 //MSPlot["MS_MHcandidate_boxWHqq_mu"]->Fill(mHcandidate,datasets[d], true, LuminosityMu*scaleFactor);
										}																		
											
										float mHq;
										int myindex; //index (of the vector of anti-b-tagged jets) of the leading or subleading anti-b-tagged jet, furthest from the reconstructed Higgs direction...
									  if(doAntiBtagging)
										{
										  if(selectedAntiBtaggedJets[0].DeltaR(Hcandidate) > selectedAntiBtaggedJets[1].DeltaR(Hcandidate)) myindex = 0;
										  else myindex = 1;
										}
										else
										{
										  if(selectedJets[0].DeltaR(Hcandidate) > selectedJets[1].DeltaR(Hcandidate)) myindex = 0;
										  else myindex = 1;
										}
													
													
										MSPlot["MS_PtH_preboxWHqq_mu"]->Fill(Hcandidate.Pt(),datasets[d], true, LuminosityMu*scaleFactor);
										if((doHptcut && Hcandidate.Pt() > Hptcut) || !doHptcut)
									  {																			
								      if((selectedJets.size() >= 4 && fabs(mHcandidate - mHiggs) < deltamH) || selectedJets.size() == 3)
										  //if(fabs(mHcandidate - mHiggs) < deltamH)
										  {        
											   float mHq;
												 if(doAntiBtagging) mHq = (Hcandidate + selectedAntiBtaggedJets[myindex]).M();
												 else mHq = (Hcandidate + selectedJets[myindex]).M();
											   MSPlot["MS_MHq_boxWHqq_mu"]->Fill(mHq,datasets[d], true, LuminosityMu*scaleFactor);
												 //for 'flat' binning
										     float fillweight = datasets[d]->NormFactor () * LuminosityMu * scaleFactor;
												 if(makeflatbinning) VarWHqq_mu.Fill_forBinning(mHq,fillweight);
					               else if(useflatbinning) VarWHqq_mu.Fill(mHq,fillweight,dataSetName);
									    }
										  MSPlot["MS_St_boxWHqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
										}
									}	
							}
						}
						else if(doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() < 150 && selectedAntiBtaggedJets[0].Pt() > 100) && (selectedAntiBtaggedJets[1].Pt() < 100 && selectedAntiBtaggedJets[1].Pt() > 50)) 
						{
						   if((doMETcutboxWHqq && met > METcutboxWHqq) || !doMETcutboxWHqq)
							 {
							    //sort of control region for WHqq
							    int nBtags = selectedBtaggedJets.size();	
							    if(!doBtagging || (doBtagging && nBtags>=1))
							    {									     
										if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							           BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								         //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								         scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }
										
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
										 MSPlot["MS_MHcandidate_preboxWHqq_mu"]->Fill(mHcandidate,datasets[d], true, LuminosityMu*scaleFactor);					 
										}
										else if(selectedJets.size()==3)
										{
										 //cout<<"selectedAvailableJetsindices_sortedbTagCSV[0].first = "<<selectedAvailableJetsindices_sortedbTagCSV[0].first<<endl;
										 Hcandidate = selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first]; //just the b-tagged jet; sort of assuming it's a fat jet... (but if it's realistic to assume a fat Higgs jet would be b-tagged...?)
										 //mHcandidate = Hcandidate.M();
										}																		
											
										float mHq;
										int myindex; //index (of the vector of anti-b-tagged jets) of the leading or subleading anti-b-tagged jet, furthest from the reconstructed Higgs direction...
									  if(doAntiBtagging)
										{
										   if(selectedAntiBtaggedJets[0].DeltaR(Hcandidate) > selectedAntiBtaggedJets[1].DeltaR(Hcandidate)) myindex = 0;
										   else myindex = 1;
										}
										else 
										{
										   if(selectedJets[0].DeltaR(Hcandidate) > selectedJets[1].DeltaR(Hcandidate)) myindex = 0;
										   else myindex = 1;
										}
																										
										MSPlot["MS_PtH_preboxWHqq_LowerJetCuts_LowerHpt_mu"]->Fill(Hcandidate.Pt(),datasets[d], true, LuminosityMu*scaleFactor);
										if(Hcandidate.Pt() < Hptcut)
									  {							
										  //note I reversed the Hwindow cut												
								      if((selectedJets.size() >= 4 && fabs(mHcandidate - mHiggs) > deltamH) || selectedJets.size() == 3)
										  //if(fabs(mHcandidate - mHiggs) < deltamH)
										  {        
												 MSPlot["MS_St_preboxWHqq_LowerJetCuts_LowerHpt_offH_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
											   if(doAntiBtagging)
												 {
												    MSPlot["MS_LeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_mu"]->Fill(selectedAntiBtaggedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
											      MSPlot["MS_SubLeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_mu"]->Fill(selectedAntiBtaggedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
												 }		
												 else 
												 {
												    MSPlot["MS_LeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
											      MSPlot["MS_SubLeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
												 }										 
									    }	 
										}													 
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
//				  selecTableMu.Fill(d,1,scaleFactor);
				  if( selectedJets.size() >= 2 )
			    {
//					  selecTableMu.Fill(d,2,scaleFactor);
				    if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				    {
//						   selecTableMu.Fill(d,3,scaleFactor);
							 if(!doMETcutboxWWqq || (doMETcutboxWWqq && met>METcutboxWWqq))
						   {				
					      if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0)
						    {
//							    selecTableMu.Fill(d,4,scaleFactor);
                  //lepton charges should be opposite
                  if((indicesMuPlus.size()==1 && indicesElMinus.size()==1) || (indicesMuMinus.size()==1 && indicesElPlus.size()==1))
									{
							      if(verbosityBox) cout<<"category: WWqq (mu)"<<endl; //(no need to check the invariant mass of the leptons, they should not come from a Z...)
							      if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							         BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								       //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
											 MSPlot["MS_BtagWeight_boxWWqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								       scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }
					          MSPlot["MS_St_boxWWqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);	
										//for 'flat' binning
										float fillweight = datasets[d]->NormFactor () * LuminosityMu * scaleFactor;
										if(makeflatbinning) VarWWqq_mu.Fill_forBinning(ST,fillweight);
					          else if(useflatbinning) VarWWqq_mu.Fill(ST,fillweight,dataSetName);
								  }
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
						 if((doMETcutboxWZqq && met > METcutboxWZqq) || !doMETcutboxWZqq)
						 {		    
							if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0)
						  {
							 if(indicesElPlus.size()==1 && indicesElMinus.size()==1)
							 {
							   if(verbosityBox) cout<<"category: WZqq (mu)"<<endl;
							   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			           {
							         BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								       //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
											 MSPlot["MS_BtagWeight_boxWZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								       scaleFactor = 	scaleFactor*BtagMCWeight; 
							   }
								 
								 TLorentzVector Zcandidate = (selectedElectrons[0] + selectedElectrons[1]);
								 MSPlot["MS_nEvts_boxWZqq_MuElEl_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						     MSPlot["MS_nEvts_boxWZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						     MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
								 massleptons = (selectedElectrons[0] + selectedElectrons[1]).M();
								 if(fabs(Zmass - massleptons) < Zmasswindow)
								 {
									    MSPlot["MS_PtZ_preboxWZqq_Zwindow_mu"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityMu*scaleFactor);
									    if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
								      {
											   MSPlot["MS_nEvts_boxWZqq_Zwindow_MuElEl_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
											   MSPlot["MS_nEvts_boxWZqq_Zwindow_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
								         MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
									    }
									}
							 }								    
						  }
						 }
					  }
					}
				}
				else
				{
				   cout<<"FOUND >= (incompatible) 4-lepton or 5-lepton event!"<<endl;		
				}									
			 }
			 else if(SelectednMu == 2)
			 {
			  if(SelectednEl == 0)
				{
			    massleptons = (selectedMuons[0] + selectedMuons[1]).M();
					MSPlot["MS_Mleptons_preboxZqq_mu"]->Fill(massleptons,datasets[d], true, LuminosityMu*scaleFactor);	
					
					if( fabs(Zmass - massleptons) < Zmasswindow)
					{				
					   //sort of control region				  
					   if(selectedJets.size() >= 1)
					   {
					     if((selectedJets[0].Pt() >= 150) && (selectedJets[0].Pt() <= 200))
						   {
						     if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							   {	
								   if(indicesMuPlus.size()==1 && indicesMuMinus.size()==1)
									 {
								       float myBtagMCWeight = 1;
										   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                 {
							            myBtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
							         }
											 if(selectedForwardJets.size() == 0) MSPlot["MS_nJets_preboxZqq_leadingjet150-200_noForwardJet_noNumberCentralJetCut_mu"]->Fill(selectedJets.size(), datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
										   else MSPlot["MS_nJets_preboxZqq_leadingjet150-200_noNumberCentralJetCut_mu"]->Fill(selectedJets.size(), datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
								   }
								 }
						   }
				     }
						
				     if( (selectedJets.size() == 1 || selectedJets.size() == 2))
			       {
						  if(selectedForwardJets.size() == 1)
							{
						   if(selectedJets[0].Pt() >= LeadingJetPtCut)
				       {
							    //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
							    
									if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						      {
									 if(indicesMuPlus.size()==1 && indicesMuMinus.size()==1)
									 {
									   if(verbosityBox) cout<<"category: Zqq (mu)"<<endl;
									   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			               {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							       }
									   ST = ST + selectedForwardJets[0].Pt();		
										 								 
										 TLorentzVector Zcandidate = (selectedMuons[0] + selectedMuons[1]);
										 MSPlot["MS_PtZ_preboxZqq_mu"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityMu*scaleFactor);
										 if((doZptcut && (Zcandidate.Pt() > Zptcut)) || !doZptcut)
										 {
							         MSPlot["MS_St_boxZqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
									     //calculate mass of (two-lepton and leading jet)-system
									     float massZq = (Zcandidate + selectedJets[0]).M();
									     MSPlot["MS_MQ_boxZqq_mu"]->Fill(massZq,datasets[d], true, LuminosityMu*scaleFactor);
											 //for 'flat' binning
										   float fillweight = datasets[d]->NormFactor () * LuminosityMu * scaleFactor;
										   if(makeflatbinning) VarZqq_mu.Fill_forBinning(massZq,fillweight);
					             else if(useflatbinning) VarZqq_mu.Fill(massZq,fillweight,dataSetName);
										 
										   //other plots
								       float DeltaR_ZcandLeadingJet = selectedJets[0].DeltaR(Zcandidate);
								       float DeltaEta_ZcandLeadingJet = fabs(selectedJets[0].Eta() - Zcandidate.Eta());
								       float DeltaR_LeptonLepton = selectedMuons[0].DeltaR(selectedMuons[1]);
								       MSPlot["MS_DeltaR_ZcandLeadingJet_boxZqq_mu"]->Fill(DeltaR_ZcandLeadingJet,datasets[d], true, LuminosityMu*scaleFactor);								 
								       MSPlot["MS_DeltaEta_ZcandLeadingJet_boxZqq_mu"]->Fill(DeltaEta_ZcandLeadingJet,datasets[d], true, LuminosityMu*scaleFactor);
								       MSPlot["MS_DeltaR_LeptonLepton_boxZqq_mu"]->Fill(DeltaR_LeptonLepton,datasets[d], true, LuminosityMu*scaleFactor);
									   }
									 }
								  }
							 }
						  }
							else if(selectedForwardJets.size() == 0)
					    {
					      //sort of control region
							  if((selectedJets[0].Pt() >= 150) && (selectedJets[0].Pt() <= 200))
						    {
						     if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							   {	
								     if(indicesMuPlus.size()==1 && indicesMuMinus.size()==1)
										 {
										   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                 {
							            BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								          scaleFactor = 	scaleFactor*BtagMCWeight; 
							         }
											 MSPlot["MS_MET_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(met, datasets[d], true, LuminosityMu*scaleFactor);
						           MSPlot["MS_HT_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);
							         MSPlot["MS_St_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
						           MSPlot["MS_LeadingJetPt_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							         MSPlot["MS_LeadingJetEta_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
							         if(selectedJets.size()>1)
							         {
							           MSPlot["MS_SubleadingJetPt_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
						        	   MSPlot["MS_SubleadingJetEta_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
					          	 }
					          	 MSPlot["MS_Muon1Pt_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							         MSPlot["MS_Muon1Eta_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
											 MSPlot["MS_Muon2Pt_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedMuons[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
							         MSPlot["MS_Muon2Eta_preboxZqq_leadingjet150-200_noForwardJet_mu"]->Fill(selectedMuons[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
										}
								 }
						    }					 
					    } 						 
				     }
				     else if((doAntiBtagging && selectedAntiBtaggedJets.size()>=2 && selectedJets.size() >= 3) || (!doAntiBtagging && selectedJets.size() >= 3)) //used to be just selectedJets.size() >= 3
			       {
						   if((doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() >= LeadingJetPtCut) && (selectedAntiBtaggedJets[1].Pt() >= SubLeadingJetPtCut)) || (!doAntiBtagging && (selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))) //used to be just if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				       {
									int nBtags = selectedBtaggedJets.size();	
							    if(!doBtagging || (doBtagging && nBtags>=1)) //used to be nBjetsPresent(selectedJets_bTagCSV,btagWP)>=1
							    {
									     if(indicesMuPlus.size()==1 && indicesMuMinus.size()==1)
											 {
											   if(verbosityBox) cout<<"category: ZHqq (mu)"<<endl;
												 if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                   {
							                BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								              //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
													  	MSPlot["MS_BtagWeight_boxZHqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								              scaleFactor = 	scaleFactor*BtagMCWeight; 
							           }
												 float massZq;
												 TLorentzVector Zcandidate = (selectedMuons[0] + selectedMuons[1]);
												 MSPlot["MS_PtZ_preboxZHqq_mu"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityMu*scaleFactor);
												 if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
												 {
												   if(doAntiBtagging)
													 {
												     if(Zcandidate.DeltaR(selectedAntiBtaggedJets[0]) > Zcandidate.DeltaR(selectedAntiBtaggedJets[1]))
												      massZq = (selectedMuons[0] + selectedMuons[1] + selectedAntiBtaggedJets[0]).M();
												     else
												      massZq = (selectedMuons[0] + selectedMuons[1] + selectedAntiBtaggedJets[1]).M();
							             }
													 else
													 {
													   if(Zcandidate.DeltaR(selectedJets[0]) > Zcandidate.DeltaR(selectedJets[1]))
												      massZq = (selectedMuons[0] + selectedMuons[1] + selectedJets[0]).M();
												     else
												      massZq = (selectedMuons[0] + selectedMuons[1] + selectedJets[1]).M();													 
													 }
													 MSPlot["MS_St_boxZHqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
												   MSPlot["MS_MET_boxZHqq_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);											 
												   MSPlot["MS_MQ_boxZHqq_mu"]->Fill(massZq,datasets[d], true, LuminosityMu*scaleFactor);
													 //for 'flat' binning
										       float fillweight = datasets[d]->NormFactor () * LuminosityMu * scaleFactor;
										       if(makeflatbinning) VarZHqq_mu.Fill_forBinning(massZq,fillweight);
					                 else if(useflatbinning) VarZHqq_mu.Fill(massZq,fillweight,dataSetName);
												 }
											}												 
									}
									else if(nBtags == 0 && selectedJets.size() >= 4) //not the same as anti-b-tagging every jet
									{
									     if(indicesMuPlus.size()==1 && indicesMuMinus.size()==1)
											 {											   
									       if(verbosityBox) cout<<"category: semilep ZZqq and WZqq (mu)"<<endl; //only one of the Z's decays leptonically, the other Z or W hadronically												 
												 if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                   {
							              BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								            //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								            scaleFactor = 	scaleFactor*BtagMCWeight; 
							           }												 
												 float massZq;
												 TLorentzVector Zcandidate = (selectedMuons[0] + selectedMuons[1]);
												 vector<int> selectedAvailableJetsindices;	//available: everything except the two highest-Pt anti-tagged jets, because I'm going to assume they come from the direct decay of the VLQ							    
										     //cout<<" selectedAntiBtaggedJets.size() = "<<selectedAntiBtaggedJets.size()<<endl;
												 //cout<<" selectedAntiBtaggedJets_indices.size() = "<<selectedAntiBtaggedJets_indices.size()<<endl;
												 for(int iJet = 0; iJet < selectedJets.size(); iJet++)
			                   {
										        if(selectedAntiBtaggedJets_indices[0] != iJet && selectedAntiBtaggedJets_indices[1] != iJet) //meaning, if not one of the two leading anti-btagged jets
											      {
											         selectedAvailableJetsindices.push_back(iJet);
											      }										
										     }
												 //loop over the 'available' jets (here assumed to come from a hadronically decaying W or Z boson) and try to find the two jets closest in deltaR
												 float minDeltaR = 9999.; //doesn't seem to work at all... -> retrial, there was a bug
												 //float dijetMassClosestToZ = 9999;
												 int indexjet1, indexjet2; //two  'available' jets closest together (or with mass closest to Z mass in another version)			 
												 for(int iJet = 0; iJet < selectedAvailableJetsindices.size(); iJet++)
			                   {
												    for(int jJet = iJet+1; jJet < selectedAvailableJetsindices.size(); jJet++)
			                      {
														   //float massdijet = (selectedJets[selectedAvailableJetsindices[iJet]]+selectedJets[selectedAvailableJetsindices[jJet]]).M();
												       if(selectedJets[selectedAvailableJetsindices[iJet]].DeltaR(selectedJets[selectedAvailableJetsindices[jJet]]) < minDeltaR)
															 //if( fabs(Zmass - massdijet) <  fabs(Zmass - dijetMassClosestToZ))
															 {
															    indexjet1 = selectedAvailableJetsindices[iJet];
																	indexjet2 = selectedAvailableJetsindices[jJet];
															    minDeltaR = selectedJets[indexjet1].DeltaR(selectedJets[indexjet2]);
																	//dijetMassClosestToZ = massdijet;
															 }															 															 
														}
												 }
												 TLorentzVector Vcandidate = selectedJets[indexjet1] + selectedJets[indexjet2]; //hadronically decaying Z or W
												 MSPlot["MS_MV_preboxSemiLepZZqqWZqq_mu"]->Fill(Vcandidate.M(),datasets[d], true, LuminosityMu*scaleFactor);
												 MSPlot["MS_PtV_preboxSemiLepZZqqWZqq_mu"]->Fill(Vcandidate.Pt(),datasets[d], true, LuminosityMu*scaleFactor);	//might do a cut on pt or mass window of this too?											 
										     MSPlot["MS_PtZ_preboxSemiLepZZqqWZqq_mu"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityMu*scaleFactor);
												 if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
												 {
												   if(doAntiBtagging)
													 {
												     if(Zcandidate.DeltaR(selectedAntiBtaggedJets[0]) > Zcandidate.DeltaR(selectedAntiBtaggedJets[1]))
												      massZq = (Zcandidate + selectedAntiBtaggedJets[0]).M();
												     else
												      massZq = (Zcandidate + selectedAntiBtaggedJets[1]).M();
												   }
													 else
													 {
													   if(Zcandidate.DeltaR(selectedJets[0]) > Zcandidate.DeltaR(selectedJets[1]))
												      massZq = (Zcandidate + selectedJets[0]).M();
												     else
												      massZq = (Zcandidate + selectedJets[1]).M();													 
													 }
													 MSPlot["MS_MET_boxSemiLepZZqqWZqq_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);
							             MSPlot["MS_St_boxSemiLepZZqqWZqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
								           MSPlot["MS_HT_boxSemiLepZZqqWZqq_mu"]->Fill(HT,datasets[d], true, LuminosityMu*scaleFactor);
								           MSPlot["MS_LeadingJetPt_boxSemiLepZZqqWZqq_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								           MSPlot["MS_LeadingJetEta_boxSemiLepZZqqWZqq_mu"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
								           MSPlot["MS_SubLeadingJetPt_boxSemiLepZZqqWZqq_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								           MSPlot["MS_SubLeadingJetEta_boxSemiLepZZqqWZqq_mu"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
								           MSPlot["MS_nJets_boxSemiLepZZqqWZqq_mu"]->Fill(selectedJets.size(),datasets[d], true, LuminosityMu*scaleFactor);
								           MSPlot["MS_Muon1Pt_boxSemiLepZZqqWZqq_mu"]->Fill(selectedMuons[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								           MSPlot["MS_Muon1Eta_boxSemiLepZZqqWZqq_mu"]->Fill(selectedMuons[0].Eta(),datasets[d], true, LuminosityMu*scaleFactor);
												   MSPlot["MS_Muon2Pt_boxSemiLepZZqqWZqq_mu"]->Fill(selectedMuons[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
								           MSPlot["MS_Muon2Eta_boxSemiLepZZqqWZqq_mu"]->Fill(selectedMuons[1].Eta(),datasets[d], true, LuminosityMu*scaleFactor);												 
												   MSPlot["MS_MQ_boxSemiLepZZqqWZqq_mu"]->Fill(massZq,datasets[d], true, LuminosityMu*scaleFactor);	
													 //for 'flat' binning
										       float fillweight = datasets[d]->NormFactor () * LuminosityMu * scaleFactor;
										       if(makeflatbinning) VarSemiLepZZqqWZqq_mu.Fill_forBinning(massZq,fillweight);
					                 else if(useflatbinning) VarSemiLepZZqqWZqq_mu.Fill(massZq,fillweight,dataSetName);
												 }
										   }											 									
									}						
							 }
							 else if(doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() < 150 && selectedAntiBtaggedJets[0].Pt() > 100) && (selectedAntiBtaggedJets[1].Pt() < 100 && selectedAntiBtaggedJets[1].Pt() > 50)) 
							 {
							    //sort of control region for ZHqq
							    int nBtags = selectedBtaggedJets.size();	
							    if(!doBtagging || (doBtagging && nBtags>=1))
							    {
									     if(indicesMuPlus.size()==1 && indicesMuMinus.size()==1)
											 {
												 if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                   {
							                BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								              //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								              scaleFactor = 	scaleFactor*BtagMCWeight; 
							           }
												 TLorentzVector Zcandidate = (selectedMuons[0] + selectedMuons[1]);
												 MSPlot["MS_PtZ_preboxZHqq_LowerJetCuts_LowerZpt_mu"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityMu*scaleFactor);
												 if(Zcandidate.Pt() < 150)
												 {
							             MSPlot["MS_St_preboxZHqq_LowerJetCuts_LowerZpt_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor);
												   MSPlot["MS_MET_preboxZHqq_LowerJetCuts_LowerZpt_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor);											 
													 if(doAntiBtagging)
													 {
													    MSPlot["MS_LeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_mu"]->Fill(selectedAntiBtaggedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
													    MSPlot["MS_SubLeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_mu"]->Fill(selectedAntiBtaggedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
													 }
													 else
													 {
													    MSPlot["MS_LeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_mu"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityMu*scaleFactor);
													    MSPlot["MS_SubLeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_mu"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityMu*scaleFactor);													 
													 }
												 }
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
//					  selecTableMu.Fill(d,1,scaleFactor);
					  if( selectedJets.size() >= 2 )
			      {
//						   selecTableMu.Fill(d,2,scaleFactor);
					     if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				       {
//							    selecTableMu.Fill(d,3,scaleFactor);
							    //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
                  if(!doMETcutboxWWqq || (doMETcutboxWWqq && met>METcutboxWWqq))
						      {
									 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						       {
//									   selecTableMu.Fill(d,4,scaleFactor);
                     //lepton charges should be opposite
                     if(indicesMuPlus.size()==1 && indicesMuMinus.size()==1)
										 {
									     if(verbosityBox) cout<<"category: WWqq (mu)"<<endl;
										   float myBtagMCWeight = 1;
									     if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                 {
							           myBtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								         //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												 MSPlot["MS_BtagWeight_boxWWqq_mu"]->Fill(myBtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								         //scaleFactor = 	scaleFactor*BtagMCWeight; 
							         }
						           MSPlot["MS_St_boxWWqq_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
											 //for 'flat' binning
										   float fillweight = datasets[d]->NormFactor () * LuminosityMu * scaleFactor;
										   if(makeflatbinning) VarWWqq_mu.Fill_forBinning(ST,fillweight);
					             else if(useflatbinning) VarWWqq_mu.Fill(ST,fillweight,dataSetName);
										 }										 							
								   }
									}
						   }		
					  }
						
						//not an independent category; just some test region (ZHqq with lepton invariant mass outside Z window)
						if((doAntiBtagging && selectedAntiBtaggedJets.size()>=2 && selectedJets.size() >= 3) || (!doAntiBtagging && selectedJets.size() >= 3))
						{
						   if((doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() >= LeadingJetPtCut) && (selectedAntiBtaggedJets[1].Pt() >= SubLeadingJetPtCut)) || (!doAntiBtagging && (selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))) //used to be just if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				       {
							    int nBtags = selectedBtaggedJets.size();	
							    if(!doBtagging || (doBtagging && nBtags>=1)) //used to be nBjetsPresent(selectedJets_bTagCSV,btagWP)>=1
							    {
									     if(indicesMuPlus.size()==1 && indicesMuMinus.size()==1)
											 {
												 float massZq;
												 TLorentzVector Zcandidate = (selectedMuons[0] + selectedMuons[1]);
												 if(doAntiBtagging)
												 {
												   if(Zcandidate.DeltaR(selectedAntiBtaggedJets[0]) > Zcandidate.DeltaR(selectedAntiBtaggedJets[1]))
												    massZq = (selectedMuons[0] + selectedMuons[1] + selectedAntiBtaggedJets[0]).M();
												   else
												    massZq = (selectedMuons[0] + selectedMuons[1] + selectedAntiBtaggedJets[1]).M();
												 }
												 else
												 {
												   if(Zcandidate.DeltaR(selectedJets[0]) > Zcandidate.DeltaR(selectedJets[1]))
												    massZq = (selectedMuons[0] + selectedMuons[1] + selectedJets[0]).M();
												   else
												    massZq = (selectedMuons[0] + selectedMuons[1] + selectedJets[1]).M();												 
												 }
												 float myBtagMCWeight = 1;
											   if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                   {
							              myBtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								            //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
														//MSPlot["MS_St_antiboxZHqq_offZ_mu"]->Fill(ST,datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
												    //MSPlot["MS_MET_antiboxZHqq_offZ_mu"]->Fill(met,datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);											 
							           }
												 
												 if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
												 {
												   MSPlot["MS_MQ_antiboxZHqq_offZ_mu"]->Fill(massZq,datasets[d], true, LuminosityMu*scaleFactor*myBtagMCWeight);
												 }									 
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
						  if((doMETcutboxWZqq && met > METcutboxWZqq) || !doMETcutboxWZqq)
						  {		     
							 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						   {
							  if(indicesMuPlus.size()==1 && indicesMuMinus.size()==1)
								{
							    if(verbosityBox) cout<<"category: WZqq (mu)"<<endl;
							    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }
									
									MSPlot["MS_nEvts_boxWZqq_MuMuEl_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
					        MSPlot["MS_nEvts_boxWZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						      MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
									 
									massleptons = (selectedMuons[indicesMuPlus[0]] + selectedMuons[indicesMuMinus[0]]).M();
									if(fabs(Zmass - massleptons) < Zmasswindow)
									{
										    TLorentzVector Zcandidate = (selectedMuons[0] + selectedMuons[1]);
												MSPlot["MS_PtZ_preboxWZqq_Zwindow_mu"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityMu*scaleFactor);
												if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
												{
												   MSPlot["MS_nEvts_boxWZqq_Zwindow_MuMuEl_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
												   MSPlot["MS_nEvts_boxWZqq_Zwindow_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
								           MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);													 
												}
									}
								}
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
							 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size()))) //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						   {
							  if((indicesMuPlus.size()==1 && indicesMuMinus.size()==1) && (indicesElPlus.size()==1 && indicesElMinus.size()==1))
								{
							    if(verbosityBox) cout<<"category: ZZqq (mu)"<<endl;
							    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }
					        MSPlot["MS_nEvts_boxZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						      MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
									
									float massleptons2 = 9999;
									//probably VERY few events surviving
									if((indicesMuPlus.size()==1 && indicesMuMinus.size()==1) && (indicesElPlus.size()==1 && indicesElMinus.size()==1))
									{
									  TLorentzVector Zcandidate1 = (selectedMuons[indicesMuPlus[0]] + selectedMuons[indicesMuMinus[0]]);
									  TLorentzVector Zcandidate2 = (selectedElectrons[indicesElPlus[0]] + selectedElectrons[indicesElMinus[0]]);
										massleptons = Zcandidate1.M();
										massleptons2 = Zcandidate2.M();
										if(fabs(Zmass - massleptons) < Zmasswindow && fabs(Zmass - massleptons2) < Zmasswindow)
										{
										   if((doZptcut && Zcandidate1.Pt() > Zptcut && Zcandidate2.Pt() > Zptcut) || !doZptcut)
											 {
										     MSPlot["MS_nEvts_boxZZqq_Zwindow_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
								         MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
											 }
										}
									}
									
								}
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
					    if((doMETcutboxWZqq && met > METcutboxWZqq) || !doMETcutboxWZqq)
						  {
							 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))  //used to be nBjetsPresent(selectedJets_bTagCSV,antibtagWP)==0
						   {
							    if(verbosityBox) cout<<"category: WZqq (mu)"<<endl;
							    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }
									MSPlot["MS_nEvts_boxWZqq_MuMuMu_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
					        MSPlot["MS_nEvts_boxWZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						      MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
									
									TLorentzVector Zcandidate;
									if(indicesMuPlus.size()==2 && indicesMuMinus.size()==1)
									{
									  TLorentzVector Zcandidate_opt1 = (selectedMuons[indicesMuPlus[0]] + selectedMuons[indicesMuMinus[0]]);
										TLorentzVector Zcandidate_opt2 = (selectedMuons[indicesMuPlus[1]] + selectedMuons[indicesMuMinus[0]]);
									  float massleptons_opt1 = Zcandidate_opt1.M(); 
									  float massleptons_opt2 = Zcandidate_opt2.M(); 
										if(massleptons_opt1 < massleptons_opt2)
										{
										   massleptons = massleptons_opt1;
											 Zcandidate = Zcandidate_opt1;
										}
										else
										{
										   massleptons = massleptons_opt2;
											 Zcandidate = Zcandidate_opt2;
										}
									}
									else if(indicesMuPlus.size()==1 && indicesMuMinus.size()==2)
									{
									  TLorentzVector Zcandidate_opt1 = (selectedMuons[indicesMuPlus[0]] + selectedMuons[indicesMuMinus[0]]);
										TLorentzVector Zcandidate_opt2 = (selectedMuons[indicesMuPlus[0]] + selectedMuons[indicesMuMinus[1]]);
									  float massleptons_opt1 = Zcandidate_opt1.M(); 
									  float massleptons_opt2 = Zcandidate_opt2.M(); 
										if(massleptons_opt1 < massleptons_opt2)
										{
										   massleptons = massleptons_opt1;
											 Zcandidate = Zcandidate_opt1;
										}
										else
										{
										   massleptons = massleptons_opt2;
											 Zcandidate = Zcandidate_opt2;
										}
									}								
									
								  if(fabs(Zmass - massleptons) < Zmasswindow)
									{
									  MSPlot["MS_PtZ_preboxWZqq_Zwindow_mu"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityMu*scaleFactor);
									  if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
										{
										   MSPlot["MS_nEvts_boxWZqq_Zwindow_MuMuMu_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);										   
									     MSPlot["MS_nEvts_boxWZqq_Zwindow_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
								       MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
										}
									}
						   }
							}
					  }
					}
				}
				else 
				{
				   cout<<"FOUND >= (incompatible) 4-lepton or 5-lepton event!"<<endl;		
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
							 //events with 3 or 4 same-charge leptons not compatible with ZZ signal
							 if((indicesMuPlus.size()==2 && indicesMuMinus.size()==2))
							 {
							   if(verbosityBox) cout<<"category: ZZqq (mu)"<<endl;
							   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			           {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZZqq_mu"]->Fill(BtagMCWeight, datasets[d], true, LuminosityMu*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							   }
					       MSPlot["MS_nEvts_boxZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
						     MSPlot["MS_nEvts_boxWZqqZZqq_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
								
								 TLorentzVector Zcandidate1, Zcandidate2;								 
								 TLorentzVector Zcandidate1_opt1 = (selectedMuons[indicesMuPlus[0]] + selectedMuons[indicesMuMinus[0]]);
								 TLorentzVector Zcandidate1_opt2 = (selectedMuons[indicesMuPlus[1]] + selectedMuons[indicesMuMinus[1]]);
								 TLorentzVector Zcandidate2_opt1 = (selectedMuons[indicesMuPlus[0]] + selectedMuons[indicesMuMinus[1]]);
								 TLorentzVector Zcandidate2_opt2 = (selectedMuons[indicesMuPlus[1]] + selectedMuons[indicesMuMinus[0]]);
								 
								 float massleptons2 = 9999;
								 float massleptons_opt1 = Zcandidate1_opt1.M();
								 float massleptons2_opt1 = Zcandidate2_opt1.M();
								 float massleptons_opt2 = Zcandidate1_opt2.M();
								 float massleptons2_opt2 = Zcandidate2_opt2.M();
									  
								 if(massleptons_opt1 < massleptons_opt2)
								 {
								    Zcandidate1 = Zcandidate1_opt1;
								    massleptons = massleptons_opt1;
								 }
								 else
								 {
								    Zcandidate1 = Zcandidate1_opt2;
								    massleptons = massleptons_opt2;
								 }									
										
								 if(massleptons2_opt1 < massleptons2_opt2)
								 {
								    Zcandidate2 = Zcandidate2_opt1;
								    massleptons2 = massleptons2_opt1;
								 }
								 else
								 {
								    Zcandidate2 = Zcandidate2_opt2;
								    massleptons2 = massleptons2_opt2;		
								 }								
															
								 if(fabs(Zmass - massleptons) < Zmasswindow && fabs(Zmass - massleptons2) < Zmasswindow)
								 {
								    if((doZptcut && Zcandidate1.Pt() > Zptcut && Zcandidate2.Pt() > Zptcut) || !doZptcut)
										{
								       MSPlot["MS_nEvts_boxZZqq_Zwindow_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
								       MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_mu"]->Fill(1,datasets[d], true, LuminosityMu*scaleFactor);
										}
								 }
							 }
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
			//else if(ElTrigged) //or not 'else if' but 'if'???? (then it should be consistent with the TreeCreator...)
			//if(ElTrigged) //changed back to if
			else if(ElChannel)
			{
			 if(SelectednMu == 0)
			 {
			  if(SelectednEl == 1)
				{									
				  //sort of control region
				  if(selectedJets.size() >= 1)
					{
					    if((selectedJets[0].Pt() >= 150) && (selectedJets[0].Pt() <= 200))
						  {
							   float MT = sqrt(2*selectedElectrons[0].Pt()*met*(1-cos(selectedElectrons[0].DeltaPhi(metvector))));
						     if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							   {	
								    if((!doMETcutboxWqq || (doMETcutboxWqq && met>METcutboxWqq)) && (!doMTcut || (doMTcut && MT>MTcut)))
							      {
										   float myBtagMCWeight = 1;
										   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                 {
							            myBtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
							         }	
											 if(selectedForwardJets.size() == 0) MSPlot["MS_nJets_preboxWqq_leadingjet150-200_noForwardJet_noNumberCentralJetCut_el"]->Fill(selectedJets.size(), datasets[d], true, LuminosityEl*scaleFactor*myBtagMCWeight);
								       else MSPlot["MS_nJets_preboxWqq_leadingjet150-200_noNumberCentralJetCut_el"]->Fill(selectedJets.size(), datasets[d], true, LuminosityEl*scaleFactor*myBtagMCWeight);
										}
								 }
						  }
				  }
				
				
				  if(selectedJets.size() == 1 || selectedJets.size() == 2)
			    {
					 if(selectedForwardJets.size() == 1)
					 { 
					   ST = ST + selectedForwardJets[0].Pt();
						 HT = HT + selectedForwardJets[0].Pt();
					   MSPlot["MS_MET_preboxWqq_noleadingjetcut_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
					   MSPlot["MS_HT_preboxWqq_noleadingjetcut_el"]->Fill(HT,datasets[d], true, LuminosityEl*scaleFactor);
						 MSPlot["MS_LeadingJetPt_preboxWqq_noleadingjetcut_el"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
						 MSPlot["MS_LeadingJetEta_preboxWqq_noleadingjetcut_el"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
						 if(selectedJets.size()>1)
						 {
								  MSPlot["MS_SubLeadingJetPt_preboxWqq_noleadingjetcut_el"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
								  MSPlot["MS_SubLeadingJetEta_preboxWqq_noleadingjetcut_el"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
						 }
						
						 //some sort of 'control region' (to be checked if not too much signal in it, though)
						 if((selectedJets[0].Pt() >= 150) && (selectedJets[0].Pt() <= 200))
						 { 
						  if(!doMETcutboxWqq || (doMETcutboxWqq && met>METcutboxWqq))
						  {
						   MSPlot["MS_MET_preboxWqq_leadingjet150-200_METcut_el"]->Fill(met, datasets[d], true, LuminosityEl*scaleFactor);
						   MSPlot["MS_HT_preboxWqq_leadingjet150-200_METcut_el"]->Fill(HT,datasets[d], true, LuminosityEl*scaleFactor);
							 MSPlot["MS_St_preboxWqq_leadingjet150-200_METcut_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
						   MSPlot["MS_LeadingJetPt_preboxWqq_leadingjet150-200_METcut_el"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							 MSPlot["MS_LeadingJetEta_preboxWqq_leadingjet150-200_METcut_el"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
							 if(selectedJets.size()>1)
							 {
							   MSPlot["MS_SubleadingJetPt_preboxWqq_leadingjet150-200_METcut_el"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							   MSPlot["MS_SubleadingJetEta_preboxWqq_leadingjet150-200_METcut_el"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
							 }
							 MSPlot["MS_ElectronPt_preboxWqq_leadingjet150-200_METcut_el"]->Fill(selectedElectrons[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							 MSPlot["MS_ElectronEta_preboxWqq_leadingjet150-200_METcut_el"]->Fill(selectedElectrons[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
							}
							MSPlot["MS_MET_preboxWqq_leadingjet150-200_el"]->Fill(met, datasets[d], true, LuminosityEl*scaleFactor);
						 }

					   if(selectedJets[0].Pt() >= LeadingJetPtCut)
				     { 
						    //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
								float MT = sqrt(2*selectedElectrons[0].Pt()*met*(1-cos(selectedElectrons[0].DeltaPhi(metvector))));
								MSPlot["MS_MTleptonmet_preboxWqq_el"]->Fill(MT, datasets[d], true, LuminosityEl*scaleFactor);
						    MSPlot["MS_MET_preboxWqq_noMETMTcut_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
							  MSPlot["MS_HT_preboxWqq_noMETMTcut_el"]->Fill(HT,datasets[d], true, LuminosityEl*scaleFactor);
							  MSPlot["MS_LeadingJetPt_preboxWqq_noMETMTcut_el"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							  MSPlot["MS_LeadingJetEta_preboxWqq_noMETMTcut_el"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
							  if(selectedJets.size()>1)
							  {
								  MSPlot["MS_SubLeadingJetPt_preboxWqq_noMETMTcut_el"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
								  MSPlot["MS_SubLeadingJetEta_preboxWqq_noMETMTcut_el"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
							  }
							  MSPlot["MS_ForwardJetPt_preboxWqq_noMETMTcut_el"]->Fill(selectedForwardJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							  MSPlot["MS_ForwardJetEta_preboxWqq_noMETMTcut_el"]->Fill(selectedForwardJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
							  MSPlot["MS_nJets_preboxWqq_noMETMTcut_el"]->Fill(selectedJets.size(),datasets[d], true, LuminosityEl*scaleFactor);
							  MSPlot["MS_ElectronPt_preboxWqq_noMETMTcut_el"]->Fill(selectedElectrons[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							  MSPlot["MS_ElectronEta_preboxWqq_noMETMTcut_el"]->Fill(selectedElectrons[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
								if((!doMETcutboxWqq || (doMETcutboxWqq && met>METcutboxWqq)) && (!doMTcut || (doMTcut && MT>MTcut)))
							  {	
								 MSPlot["MS_MET_preboxWqq_noBveto_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
								 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							   {	
								   if(verbosityBox) cout<<"category: Wqq (el)"<<endl;
								   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			             {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							     }
					 
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
								   TLorentzVector Wcandidate = recneutrino + selectedElectrons[0];
									 float PtW = Wcandidate.Pt();
									 
									 MSPlot["MS_PtW_preboxWqq_el"]->Fill(PtW,datasets[d], true, LuminosityEl*scaleFactor);
									 if((doWptcut && Wcandidate.Pt() > Wptcut) || !doWptcut)
									 {
									   MSPlot["MS_MET_boxWqq_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
									   MSPlot["MS_MTleptonmet_boxWqq_el"]->Fill(MT, datasets[d], true, LuminosityEl*scaleFactor);	
						         MSPlot["MS_St_boxWqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
								     MSPlot["MS_HT_boxWqq_el"]->Fill(HT,datasets[d], true, LuminosityEl*scaleFactor);
								     MSPlot["MS_LeadingJetPt_boxWqq_el"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
								     MSPlot["MS_LeadingJetEta_boxWqq_el"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
								     if(selectedJets.size()>1)
								     {
								        MSPlot["MS_SubLeadingJetPt_boxWqq_el"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
										    MSPlot["MS_SubLeadingJetEta_boxWqq_el"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
								     }
								     MSPlot["MS_ForwardJetPt_boxWqq_el"]->Fill(selectedForwardJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
								     MSPlot["MS_ForwardJetEta_boxWqq_el"]->Fill(selectedForwardJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
								     MSPlot["MS_nJets_boxWqq_el"]->Fill(selectedJets.size(),datasets[d], true, LuminosityEl*scaleFactor);
								     MSPlot["MS_ElectronPt_boxWqq_el"]->Fill(selectedElectrons[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
								     MSPlot["MS_ElectronEta_boxWqq_el"]->Fill(selectedElectrons[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
								     MSPlot["MS_MCeventWeight_boxWqq_el"]->Fill(LuminosityEl*scaleFactor*datasets[d]->NormFactor(),datasets[d], false, LuminosityEl*scaleFactor); //not scaling the events...								 								
									   float MQ = (recneutrino + selectedElectrons[0] + selectedJets[0]).M();
							       MSPlot["MS_MQ_boxWqq_el"]->Fill(MQ,datasets[d], true, LuminosityEl*scaleFactor);
								 
								     //for 'flat' binning
										 float fillweight = datasets[d]->NormFactor () * LuminosityEl * scaleFactor;
										 //cout<<"  ---> d = "<<d<<endl;
										 //cout<<"  ---> datasets["<<d<<"] = "<<datasets[d]->Name()<<endl;
										 if(makeflatbinning) VarWqq_el.Fill_forBinning(MQ,fillweight);
					           else if(useflatbinning) VarWqq_el.Fill(MQ,fillweight,dataSetName);
								     
										 //other plots
								     float DeltaR_WcandLeadingJet = selectedJets[0].DeltaR((recneutrino + selectedElectrons[0]));
								     float DeltaEta_WcandLeadingJet = fabs(selectedJets[0].Eta() - (recneutrino + selectedElectrons[0]).Eta());
								     float DeltaR_LeptonNeutrino = selectedElectrons[0].DeltaR(recneutrino);
								     MSPlot["MS_DeltaR_WcandLeadingJet_boxWqq_el"]->Fill(DeltaR_WcandLeadingJet,datasets[d], true, LuminosityEl*scaleFactor);								 
								     MSPlot["MS_DeltaEta_WcandLeadingJet_boxWqq_el"]->Fill(DeltaEta_WcandLeadingJet,datasets[d], true, LuminosityEl*scaleFactor);
								     MSPlot["MS_DeltaR_LeptonNeutrino_boxWqq_el"]->Fill(DeltaR_LeptonNeutrino,datasets[d], true, LuminosityEl*scaleFactor);
									 
									   if(selectedElectronsCharge[0] == 1)
								     {
									    if(verbosityBox) cout<<"          W+ qq (el)"<<endl;
								      MSPlot["MS_St_boxWplusqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
											MSPlot["MS_MTQ_boxWplusqq_el"]->Fill(MTQ,datasets[d], true, LuminosityEl*scaleFactor);
											MSPlot["MS_MQ_boxWplusqq_el"]->Fill(MQ,datasets[d], true, LuminosityEl*scaleFactor);
											//for 'flat' binning
										  float fillweight = datasets[d]->NormFactor () * LuminosityEl * scaleFactor;
										  if(makeflatbinning) VarWplusqq_el.Fill_forBinning(MQ,fillweight);
					            else if(useflatbinning) VarWplusqq_el.Fill(MQ,fillweight,dataSetName);
								     }
								     else if(selectedElectronsCharge[0] == -1)
								     {
									    if(verbosityBox) cout<<"          W- qq (el)"<<endl;
								      MSPlot["MS_St_boxWminusqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
											MSPlot["MS_MTQ_boxWminusqq_el"]->Fill(MTQ,datasets[d], true, LuminosityEl*scaleFactor);
											MSPlot["MS_MQ_boxWminusqq_el"]->Fill(MQ,datasets[d], true, LuminosityEl*scaleFactor);
											//for 'flat' binning
										  float fillweight = datasets[d]->NormFactor () * LuminosityEl * scaleFactor;
										  if(makeflatbinning) VarWminusqq_el.Fill_forBinning(MQ,fillweight);
					            else if(useflatbinning) VarWminusqq_el.Fill(MQ,fillweight,dataSetName);
								     }
								     else cout<<"WARNING: electron charge "<<selectedElectronsCharge[1]<<" does not make sense!"<<endl;
								   }
								 }
								}
						 }
					 }
					 else if(selectedForwardJets.size() == 0)
					 {
					    //sort of control region
							if((selectedJets[0].Pt() >= 150) && (selectedJets[0].Pt() <= 200))
						  {
							   float MT = sqrt(2*selectedElectrons[0].Pt()*met*(1-cos(selectedElectrons[0].DeltaPhi(metvector))));
						     if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							   {	
								    MSPlot["MS_MET_preboxWqq_leadingjet150-200_noForwardJet_noMETMTcut_el"]->Fill(met, datasets[d], true, LuminosityEl*scaleFactor);
								    if((!doMETcutboxWqq || (doMETcutboxWqq && met>METcutboxWqq)) && (!doMTcut || (doMTcut && MT>MTcut)))
							      {
										   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                 {
							            BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								          scaleFactor = 	scaleFactor*BtagMCWeight; 
							         }
											 MSPlot["MS_MET_preboxWqq_leadingjet150-200_noForwardJet_el"]->Fill(met, datasets[d], true, LuminosityEl*scaleFactor);
											 MSPlot["MS_MTleptonmet_preboxWqq_leadingjet150-200_noForwardJet_el"]->Fill(MT, datasets[d], true, LuminosityEl*scaleFactor);	
						           MSPlot["MS_HT_preboxWqq_leadingjet150-200_noForwardJet_el"]->Fill(HT,datasets[d], true, LuminosityEl*scaleFactor);
							         MSPlot["MS_St_preboxWqq_leadingjet150-200_noForwardJet_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
						           MSPlot["MS_LeadingJetPt_preboxWqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							         MSPlot["MS_LeadingJetEta_preboxWqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
							         if(selectedJets.size()>1)
							         {
							           MSPlot["MS_SubleadingJetPt_preboxWqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
						        	   MSPlot["MS_SubleadingJetEta_preboxWqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
					          	 }
					          	 MSPlot["MS_ElectronPt_preboxWqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedElectrons[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							         MSPlot["MS_ElectronEta_preboxWqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedElectrons[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
										}
								 }
						  }					 
					 }				 
				  }
				  else if((doAntiBtagging && selectedAntiBtaggedJets.size()>=2 && selectedJets.size() >= 3) || (!doAntiBtagging && selectedJets.size() >= 3))
			    {
					   if((doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() >= LeadingJetPtCut) && (selectedAntiBtaggedJets[1].Pt() >= SubLeadingJetPtCut)) || (!doAntiBtagging && (selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut)))
				     {								
								int nBtags = selectedBtaggedJets.size();	
						    if(!doBtagging || (doBtagging && nBtags>=1))
							  {	
								  if((doMETcutboxWHqq && met > METcutboxWHqq) || !doMETcutboxWHqq)
									{
									  if(verbosityBox) cout<<"category: WHqq (el)"<<endl;
									  if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							           BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								         //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
											   MSPlot["MS_BtagWeight_boxWHqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								         scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }

										vector<pair<int,float> > selectedAvailableJetsindices_sortedbTagCSV;	//available: everything except the two highest-Pt anti-tagged jets, because I'm going to assume they come from the direct decay of the VLQ							    
										for(int iJet = 0; iJet < selectedJets.size(); iJet++)
			              {
										  if(doAntiBtagging)
											{
										    if(selectedAntiBtaggedJets_indices[0] != iJet && selectedAntiBtaggedJets_indices[1] != iJet) //meaning, if not one of the two leading anti-btagged jets
											  {
											    pair<int,float> jetindex_btagCSV = make_pair(iJet,selectedJets_bTagCSV[iJet]);	
											    selectedAvailableJetsindices_sortedbTagCSV.push_back(jetindex_btagCSV);
											  }
											}	
											else
											{
											  if(0 != iJet && 1 != iJet) //meaning, if not one of the two leading jets
											  {
											    pair<int,float> jetindex_btagCSV = make_pair(iJet,selectedJets_bTagCSV[iJet]);	
											    selectedAvailableJetsindices_sortedbTagCSV.push_back(jetindex_btagCSV);
											  }
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
										 //mHcandidate = Hcandidate.M();
										 //MSPlot["MS_MHcandidate_boxWHqq_el"]->Fill(mHcandidate,datasets[d], true, LuminosityEl*scaleFactor);
										}								
																				
										float mHq;
										int myindex; //index (of the vector of anti-b-tagged jets) of the leading or subleading anti-b-tagged jet, furthest from the reconstructed Higgs direction...
									  if(doAntiBtagging)
										{
										   if(selectedAntiBtaggedJets[0].DeltaR(Hcandidate) > selectedAntiBtaggedJets[1].DeltaR(Hcandidate)) myindex = 0;
										   else myindex = 1;
										}
										else
										{
										   if(selectedJets[0].DeltaR(Hcandidate) > selectedJets[1].DeltaR(Hcandidate)) myindex = 0;
										   else myindex = 1;
										}
												
										MSPlot["MS_PtH_preboxWHqq_el"]->Fill(Hcandidate.Pt(),datasets[d], true, LuminosityEl*scaleFactor);
										if((doHptcut && Hcandidate.Pt() > Hptcut) || !doHptcut)
									  {																									
								      if((selectedJets.size() >= 4 && fabs(mHcandidate - mHiggs) < deltamH) || selectedJets.size() == 3)
										  //if(fabs(mHcandidate - mHiggs) < deltamH)
										  {                       
											  float mHq;
												if(doAntiBtagging) mHq = (Hcandidate + selectedAntiBtaggedJets[myindex]).M();
											  else mHq = (Hcandidate + selectedJets[myindex]).M();
												MSPlot["MS_MHq_boxWHqq_el"]->Fill(mHq,datasets[d], true, LuminosityEl*scaleFactor);
												//for 'flat' binning
										    float fillweight = datasets[d]->NormFactor () * LuminosityEl * scaleFactor;
										    if(makeflatbinning) VarWHqq_el.Fill_forBinning(mHq,fillweight);
					              else if(useflatbinning) VarWHqq_el.Fill(mHq,fillweight,dataSetName);
									    }
											MSPlot["MS_St_boxWHqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
										}
									}
								}
						 }
						 else if(doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() < 150 && selectedAntiBtaggedJets[0].Pt() > 100) && (selectedAntiBtaggedJets[1].Pt() < 100 && selectedAntiBtaggedJets[1].Pt() > 50)) 
						 {
							    //sort of control region for WHqq
							    int nBtags = selectedBtaggedJets.size();	
							    if(!doBtagging || (doBtagging && nBtags>=1))
							    {
									 if((doMETcutboxWHqq && met > METcutboxWHqq) || !doMETcutboxWHqq)
									 {					     
										if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							           BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								         //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								         scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }
										
										vector<pair<int,float> > selectedAvailableJetsindices_sortedbTagCSV;	//available: everything except the two highest-Pt anti-tagged jets, because I'm going to assume they come from the direct decay of the VLQ							    
										for(int iJet = 0; iJet < selectedJets.size(); iJet++)
			              {
										  if(doAntiBtagging)
											{
										    if(selectedAntiBtaggedJets_indices[0] != iJet && selectedAntiBtaggedJets_indices[1] != iJet) //meaning, if not one of the two leading anti-btagged jets
											  {
											    pair<int,float> jetindex_btagCSV = make_pair(iJet,selectedJets_bTagCSV[iJet]);	
											    selectedAvailableJetsindices_sortedbTagCSV.push_back(jetindex_btagCSV);
											  }		
											}		
											else
											{
											  if(0 != iJet && 1 != iJet) //meaning, if not one of the two leading jets
											  {
											    pair<int,float> jetindex_btagCSV = make_pair(iJet,selectedJets_bTagCSV[iJet]);	
											    selectedAvailableJetsindices_sortedbTagCSV.push_back(jetindex_btagCSV);
											  }												
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
										 MSPlot["MS_MHcandidate_preboxWHqq_el"]->Fill(mHcandidate,datasets[d], true, LuminosityEl*scaleFactor);					 
										}
										else if(selectedJets.size()==3)
										{
										 //cout<<"selectedAvailableJetsindices_sortedbTagCSV[0].first = "<<selectedAvailableJetsindices_sortedbTagCSV[0].first<<endl;
										 Hcandidate = selectedJets[selectedAvailableJetsindices_sortedbTagCSV[0].first]; //just the b-tagged jet; sort of assuming it's a fat jet... (but if it's realistic to assume a fat Higgs jet would be b-tagged...?)
										 //mHcandidate = Hcandidate.M();
										}																		
											
										float mHq;
										int myindex; //index (of the vector of anti-b-tagged jets) of the leading or subleading anti-b-tagged jet, furthest from the reconstructed Higgs direction...
									  if(doAntiBtagging)
										{
										   if(selectedAntiBtaggedJets[0].DeltaR(Hcandidate) > selectedAntiBtaggedJets[1].DeltaR(Hcandidate)) myindex = 0;
										   else myindex = 1;
										}
										else
										{
										   if(selectedJets[0].DeltaR(Hcandidate) > selectedJets[1].DeltaR(Hcandidate)) myindex = 0;
										   else myindex = 1;
										} 
																										
										MSPlot["MS_PtH_preboxWHqq_LowerJetCuts_LowerHpt_el"]->Fill(Hcandidate.Pt(),datasets[d], true, LuminosityEl*scaleFactor);
										if(Hcandidate.Pt() < Hptcut)
									  {							
										  //note I reversed the Hwindow cut												
								      if((selectedJets.size() >= 4 && fabs(mHcandidate - mHiggs) > deltamH) || selectedJets.size() == 3)
										  //if(fabs(mHcandidate - mHiggs) < deltamH)
										  {        
												 MSPlot["MS_St_preboxWHqq_LowerJetCuts_LowerHpt_offH_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
											   if(doAntiBtagging)
										     {
												    MSPlot["MS_LeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_el"]->Fill(selectedAntiBtaggedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
											      MSPlot["MS_SubLeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_el"]->Fill(selectedAntiBtaggedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
												 }		
												 else
												 {
												    MSPlot["MS_LeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_el"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
											      MSPlot["MS_SubLeadingJetPt_preboxWHqq_LowerJetCuts_LowerHpt_offH_el"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
												 }										 
									    }	 
										}
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
					MSPlot["MS_Mleptons_preboxZqq_el"]->Fill(massleptons,datasets[d], true, LuminosityEl*scaleFactor);	
					
					if( fabs(Zmass - massleptons) < Zmasswindow)
					{								   
						 //sort of control region
				     if(selectedJets.size() >= 1)
					   {
					    if((selectedJets[0].Pt() >= 150) && (selectedJets[0].Pt() <= 200))
						  {
						     if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							   {	
								   if((indicesElPlus.size()==1 && indicesElMinus.size()==1))
							     {
								       float myBtagMCWeight = 1;
										   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                 {
							            BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
							         }
											 if(selectedForwardJets.size() == 0) MSPlot["MS_nJets_preboxZqq_leadingjet150-200_noForwardJet_noNumberCentralJetCut_el"]->Fill(selectedJets.size(), datasets[d], true, LuminosityEl*scaleFactor*myBtagMCWeight);
									     else MSPlot["MS_nJets_preboxZqq_leadingjet150-200_noNumberCentralJetCut_el"]->Fill(selectedJets.size(), datasets[d], true, LuminosityEl*scaleFactor*myBtagMCWeight);
								   }
								 }
						  }
				     }
					
							
				     if( (selectedJets.size() == 1 || selectedJets.size() == 2))
			       {
						  if(selectedForwardJets.size() == 1 )
						 	{					 
						    if(selectedJets[0].Pt() >= LeadingJetPtCut)
				        {
								  //if(!useBtagInfoFromFile) bTool_L->FillMCEfficiencyHistos(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV);
								  
									if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						      {
									 if((indicesElPlus.size()==1 && indicesElMinus.size()==1))
							     {
									  if(verbosityBox) cout<<"category: Zqq (el)"<<endl;
									  if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }
							      ST = ST + selectedForwardJets[0].Pt();
							       
										TLorentzVector Zcandidate = (selectedElectrons[0] + selectedElectrons[1]);
										MSPlot["MS_PtZ_preboxZqq_el"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityEl*scaleFactor);
										if((doZptcut && (Zcandidate.Pt() > Zptcut)) || !doZptcut)
										{
										  MSPlot["MS_St_boxZqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
									    //calculate mass of (two-lepton and leading jet)-system
									    float massZq = (Zcandidate + selectedJets[0]).M();
									    MSPlot["MS_MQ_boxZqq_el"]->Fill(massZq,datasets[d], true, LuminosityEl*scaleFactor);
										
										  //for 'flat' binning
										  float fillweight = datasets[d]->NormFactor () * LuminosityEl * scaleFactor;
										  if(makeflatbinning) VarZqq_el.Fill_forBinning(massZq,fillweight);
					            else if(useflatbinning) VarZqq_el.Fill(massZq,fillweight,dataSetName);
									  	
											//other plots
								      float DeltaR_ZcandLeadingJet = selectedJets[0].DeltaR((selectedElectrons[0] + selectedElectrons[1]));
								      float DeltaEta_ZcandLeadingJet = fabs(selectedJets[0].Eta() - (selectedElectrons[0] + selectedElectrons[1]).Eta());
								      float DeltaR_LeptonLepton = selectedElectrons[0].DeltaR(selectedElectrons[1]);
								      MSPlot["MS_DeltaR_ZcandLeadingJet_boxZqq_el"]->Fill(DeltaR_ZcandLeadingJet,datasets[d], true, LuminosityEl*scaleFactor);								 
								      MSPlot["MS_DeltaEta_ZcandLeadingJet_boxZqq_el"]->Fill(DeltaEta_ZcandLeadingJet,datasets[d], true, LuminosityEl*scaleFactor);
								      MSPlot["MS_DeltaR_LeptonLepton_boxZqq_el"]->Fill(DeltaR_LeptonLepton,datasets[d], true, LuminosityEl*scaleFactor);
									  }
									 }
								  }
								}
							}
							else if(selectedForwardJets.size() == 0)
					    {
					      //sort of control region
							  if((selectedJets[0].Pt() >= 150) && (selectedJets[0].Pt() <= 200))
						    {
						     if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
							   {	
								   if((indicesElPlus.size()==1 && indicesElMinus.size()==1))
							     {
										   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                 {
							            BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								          scaleFactor = 	scaleFactor*BtagMCWeight; 
							         }
											 MSPlot["MS_MET_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(met, datasets[d], true, LuminosityEl*scaleFactor);
						           MSPlot["MS_HT_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(HT,datasets[d], true, LuminosityEl*scaleFactor);
							         MSPlot["MS_St_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
						           MSPlot["MS_LeadingJetPt_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							         MSPlot["MS_LeadingJetEta_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
							         if(selectedJets.size()>1)
							         {
							           MSPlot["MS_SubleadingJetPt_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
						        	   MSPlot["MS_SubleadingJetEta_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
					          	 }
					          	 MSPlot["MS_Electron1Pt_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedElectrons[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							         MSPlot["MS_Electron1Eta_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedElectrons[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
											 MSPlot["MS_Electron2Pt_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedElectrons[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
							         MSPlot["MS_Electron2Eta_preboxZqq_leadingjet150-200_noForwardJet_el"]->Fill(selectedElectrons[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
								   }
								 }
						    }					 
					    }		
				     }
				     else if((doAntiBtagging && selectedAntiBtaggedJets.size()>=2 && selectedJets.size() >= 3) || (!doAntiBtagging && selectedJets.size() >= 3))
			       {
						    if((doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() >= LeadingJetPtCut) && (selectedAntiBtaggedJets[1].Pt() >= SubLeadingJetPtCut)) || (!doAntiBtagging && (selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut)))
				        {										
								   int nBtags = selectedBtaggedJets.size();
									 if(!doBtagging || (doBtagging && nBtags>=1))
							     {	
							        if((indicesElPlus.size()==1 && indicesElMinus.size()==1))
							        {
											   if(verbosityBox) cout<<"category: ZHqq (el)"<<endl;
												 if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                   {
							               BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								             //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
														 MSPlot["MS_BtagWeight_boxZHqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								             scaleFactor = 	scaleFactor*BtagMCWeight; 
							           }
                         float massZq;
												 TLorentzVector Zcandidate = (selectedElectrons[0] + selectedElectrons[1]);
												 MSPlot["MS_PtZ_preboxZHqq_el"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityEl*scaleFactor);
												 if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
												 {
												   if(doAntiBtagging)
													 {
												     if(Zcandidate.DeltaR(selectedAntiBtaggedJets[0]) > Zcandidate.DeltaR(selectedAntiBtaggedJets[1]))
												      massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedAntiBtaggedJets[0]).M();
												     else
												      massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedAntiBtaggedJets[1]).M();													 
							             }
													 else
													 {
													   if(Zcandidate.DeltaR(selectedJets[0]) > Zcandidate.DeltaR(selectedJets[1]))
												      massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedJets[0]).M();
												     else
												      massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedJets[1]).M();													 
													 }
													 MSPlot["MS_St_boxZHqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
												   MSPlot["MS_MET_boxZHqq_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);											 
												   MSPlot["MS_MQ_boxZHqq_el"]->Fill(massZq,datasets[d], true, LuminosityEl*scaleFactor);
													 //for 'flat' binning
										       float fillweight = datasets[d]->NormFactor () * LuminosityEl * scaleFactor;
										       if(makeflatbinning) VarZHqq_el.Fill_forBinning(massZq,fillweight);
					                 else if(useflatbinning) VarZHqq_el.Fill(massZq,fillweight,dataSetName);												 
											   }
											}
									 }
									 else if(nBtags == 0 && selectedJets.size() >= 4) //not the same as anti-b-tagging every jet
									 {
									     if(indicesElPlus.size()==1 && indicesElMinus.size()==1)
											 {
									       if(verbosityBox) cout<<"category: semilep ZZqq and WZqq (el)"<<endl; //only one of the Z's decays leptonically, the other Z or W hadronically												
												 if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                   {
							                BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								            //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								            scaleFactor = 	scaleFactor*BtagMCWeight; 
							           }												
												 float massZq;
												 TLorentzVector Zcandidate = (selectedElectrons[0] + selectedElectrons[1]);									 
												 vector<int> selectedAvailableJetsindices;	//available: everything except the two highest-Pt anti-tagged jets, because I'm going to assume they come from the direct decay of the VLQ							    
										     for(int iJet = 0; iJet < selectedJets.size(); iJet++)
			                   {
												    if(doAntiBtagging)
														{
										          if(selectedAntiBtaggedJets_indices[0] != iJet && selectedAntiBtaggedJets_indices[1] != iJet) //meaning, if not one of the two leading anti-btagged jets
											        {
											          selectedAvailableJetsindices.push_back(iJet);
											        }
														}
														else
														{
														  if(0 != iJet && 1 != iJet) //meaning, if not one of the two leading jets
											        {
											          selectedAvailableJetsindices.push_back(iJet);
											        }
														}									
										     }
												 //loop over the 'available' jets (here assumed to come from a hadronically decaying W or Z boson) and try to find the two jets closest in deltaR
												 float minDeltaR = 9999.; //doesn't seem to work at all... -> retial, there was a bug
												 //float dijetMassClosestToZ = 9999;
												 int indexjet1, indexjet2; //two  'available' jets closest together (or with mass closest to Z mass in another version)		 
												 for(int iJet = 0; iJet < selectedAvailableJetsindices.size(); iJet++)
			                   {
												    for(int jJet = iJet+1; jJet < selectedAvailableJetsindices.size(); jJet++)
			                      {
														   //float massdijet = (selectedJets[selectedAvailableJetsindices[iJet]]+selectedJets[selectedAvailableJetsindices[jJet]]).M();
												       if(selectedJets[selectedAvailableJetsindices[iJet]].DeltaR(selectedJets[selectedAvailableJetsindices[jJet]]) < minDeltaR)
															 //if( fabs(Zmass - massdijet) <  fabs(Zmass - dijetMassClosestToZ))
															 {
															    indexjet1 = selectedAvailableJetsindices[iJet];
																	indexjet2 = selectedAvailableJetsindices[jJet];
															    minDeltaR = selectedJets[indexjet1].DeltaR(selectedJets[indexjet2]);
																	//dijetMassClosestToZ = massdijet;
															 }															 															 
														}
												 }
												 TLorentzVector Vcandidate = selectedJets[indexjet1] + selectedJets[indexjet2]; //hadronically decaying Z or W
												 MSPlot["MS_MV_preboxSemiLepZZqqWZqq_el"]->Fill(Vcandidate.M(),datasets[d], true, LuminosityEl*scaleFactor);
												 MSPlot["MS_PtV_preboxSemiLepZZqqWZqq_el"]->Fill(Vcandidate.Pt(),datasets[d], true, LuminosityEl*scaleFactor);	//might do a cut on pt or mass window of this too?											 
										     MSPlot["MS_PtZ_preboxSemiLepZZqqWZqq_el"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityEl*scaleFactor);
												 if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
												 {				
												   if(doAntiBtagging)
													 { 
												     if(Zcandidate.DeltaR(selectedAntiBtaggedJets[0]) > Zcandidate.DeltaR(selectedAntiBtaggedJets[1]))
												      massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedAntiBtaggedJets[0]).M();
												     else
												      massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedAntiBtaggedJets[1]).M();							
												   }
													 else
													 {
													   if(Zcandidate.DeltaR(selectedJets[0]) > Zcandidate.DeltaR(selectedJets[1]))
												      massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedJets[0]).M();
												     else
												      massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedJets[1]).M();
													 }
													 MSPlot["MS_MET_boxSemiLepZZqqWZqq_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);
							             MSPlot["MS_St_boxSemiLepZZqqWZqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
								           MSPlot["MS_HT_boxSemiLepZZqqWZqq_el"]->Fill(HT,datasets[d], true, LuminosityEl*scaleFactor);
								           MSPlot["MS_LeadingJetPt_boxSemiLepZZqqWZqq_el"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
								           MSPlot["MS_LeadingJetEta_boxSemiLepZZqqWZqq_el"]->Fill(selectedJets[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
								           MSPlot["MS_SubLeadingJetPt_boxSemiLepZZqqWZqq_el"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
								           MSPlot["MS_SubLeadingJetEta_boxSemiLepZZqqWZqq_el"]->Fill(selectedJets[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
								           MSPlot["MS_nJets_boxSemiLepZZqqWZqq_el"]->Fill(selectedJets.size(),datasets[d], true, LuminosityEl*scaleFactor);
								           MSPlot["MS_Electron1Pt_boxSemiLepZZqqWZqq_el"]->Fill(selectedElectrons[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
								           MSPlot["MS_Electron1Eta_boxSemiLepZZqqWZqq_el"]->Fill(selectedElectrons[0].Eta(),datasets[d], true, LuminosityEl*scaleFactor);
												   MSPlot["MS_Electron2Pt_boxSemiLepZZqqWZqq_el"]->Fill(selectedElectrons[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
								           MSPlot["MS_Electron2Eta_boxSemiLepZZqqWZqq_el"]->Fill(selectedElectrons[1].Eta(),datasets[d], true, LuminosityEl*scaleFactor);											 
												   MSPlot["MS_MQ_boxSemiLepZZqqWZqq_el"]->Fill(massZq,datasets[d], true, LuminosityEl*scaleFactor);	
													 //for 'flat' binning
										       float fillweight = datasets[d]->NormFactor () * LuminosityEl * scaleFactor;
										       if(makeflatbinning) VarSemiLepZZqqWZqq_el.Fill_forBinning(massZq,fillweight);
					                 else if(useflatbinning) VarSemiLepZZqqWZqq_el.Fill(massZq,fillweight,dataSetName);												 									
									      }
									    }
								  }
								}
								else if(doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() < 150 && selectedAntiBtaggedJets[0].Pt() > 100) && (selectedAntiBtaggedJets[1].Pt() < 100 && selectedAntiBtaggedJets[1].Pt() > 50)) 
							  {
							    //sort of control region for ZHqq
							    int nBtags = selectedBtaggedJets.size();	
							    if(!doBtagging || (doBtagging && nBtags>=1))
							    {
									     if(indicesElPlus.size()==1 && indicesElMinus.size()==1)
											 {
												 if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                   {
							                BtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								              //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
								              scaleFactor = 	scaleFactor*BtagMCWeight; 
							           }
												 TLorentzVector Zcandidate = (selectedElectrons[0] + selectedElectrons[1]);
												 MSPlot["MS_PtZ_preboxZHqq_LowerJetCuts_LowerZpt_el"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityEl*scaleFactor);
												 if(Zcandidate.Pt() < 150)
												 {
							             MSPlot["MS_St_preboxZHqq_LowerJetCuts_LowerZpt_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
												   MSPlot["MS_MET_preboxZHqq_LowerJetCuts_LowerZpt_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor);											 
													 if(doAntiBtagging)
													 {
													   MSPlot["MS_LeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_el"]->Fill(selectedAntiBtaggedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
													   MSPlot["MS_SubLeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_el"]->Fill(selectedAntiBtaggedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
												   }
													 else
													 {
													   MSPlot["MS_LeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_el"]->Fill(selectedJets[0].Pt(),datasets[d], true, LuminosityEl*scaleFactor);
													   MSPlot["MS_SubLeadingJetPt_preboxZHqq_LowerJetCuts_LowerZpt_el"]->Fill(selectedJets[1].Pt(),datasets[d], true, LuminosityEl*scaleFactor);													 
													 }
												 }
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
					  if( selectedJets.size() >= 2 )
			      {
					      if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				        {
						      if(!doMETcutboxWWqq || (doMETcutboxWWqq && met>METcutboxWWqq))
									{
									 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						       {
									   //lepton charges should be opposite
                     if(indicesElPlus.size()==1 && indicesElMinus.size()==1)
										 {
									     if(verbosityBox) cout<<"category: WWqq (el)"<<endl;
									     float myBtagMCWeight = 1;
										   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                 {
							            myBtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								          //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
										  		MSPlot["MS_BtagWeight_boxWWqq_el"]->Fill(myBtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								          //scaleFactor = 	scaleFactor*BtagMCWeight; 
							         }
						           MSPlot["MS_St_boxWWqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor*myBtagMCWeight);
											 //for 'flat' binning
										   float fillweight = datasets[d]->NormFactor () * LuminosityEl * scaleFactor;
										   if(makeflatbinning) VarWWqq_el.Fill_forBinning(ST,fillweight);
					             else if(useflatbinning) VarWWqq_el.Fill(ST,fillweight,dataSetName);
							       }
									 }
									}		
						    }						 
						}
						 
						//not an independent category; just some test region (ZHqq with lepton invariant mass outside Z window)
						if((doAntiBtagging && selectedAntiBtaggedJets.size()>=2 && selectedJets.size() >= 3) || (!doAntiBtagging && selectedJets.size() >= 3))
						{
						   if((doAntiBtagging && (selectedAntiBtaggedJets[0].Pt() >= LeadingJetPtCut) && (selectedAntiBtaggedJets[1].Pt() >= SubLeadingJetPtCut)) || (!doAntiBtagging && (selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))) //used to be just if((selectedJets[0].Pt() >= LeadingJetPtCut) && (selectedJets[1].Pt() >= SubLeadingJetPtCut))
				       {
							    int nBtags = selectedBtaggedJets.size();	
							    if(!doBtagging || (doBtagging && nBtags>=1)) //used to be nBjetsPresent(selectedJets_bTagCSV,btagWP)>=1
							    {
									     if(indicesElPlus.size()==1 && indicesElMinus.size()==1)
											 {
												 float massZq;
												 TLorentzVector Zcandidate = (selectedElectrons[0] + selectedElectrons[1]);
												 if(doAntiBtagging)
												 {
												   if(Zcandidate.DeltaR(selectedAntiBtaggedJets[0]) > Zcandidate.DeltaR(selectedAntiBtaggedJets[1]))
												    massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedAntiBtaggedJets[0]).M();
												   else
												    massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedAntiBtaggedJets[1]).M();
											   }
												 else
												 {
												   if(Zcandidate.DeltaR(selectedJets[0]) > Zcandidate.DeltaR(selectedJets[1]))
												    massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedJets[0]).M();
												   else
												    massZq = (selectedElectrons[0] + selectedElectrons[1] + selectedJets[1]).M();
												 }
												 float myBtagMCWeight = 1;
												 if(doBtagging && doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			                   {
							              myBtagMCWeight = bTool_LT->getMCEventWeight_LT(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								            //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
														//MSPlot["MS_St_antiboxZHqq_offZ_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor*myBtagMCWeight);
												    //MSPlot["MS_MET_antiboxZHqq_offZ_el"]->Fill(met,datasets[d], true, LuminosityEl*scaleFactor*myBtagMCWeight);											 
							           }	
												 
												 if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
												 {
												   MSPlot["MS_MQ_antiboxZHqq_offZ_el"]->Fill(massZq,datasets[d], true, LuminosityEl*scaleFactor*myBtagMCWeight);											 
											   }
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
							if((doMETcutboxWZqq && met > METcutboxWZqq) || !doMETcutboxWZqq)
						  {
								if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						    {
								  if(verbosityBox) cout<<"category: WZqq (el)"<<endl;
								  if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }
									MSPlot["MS_nEvts_boxWZqq_ElElEl_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
					        MSPlot["MS_nEvts_boxWZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
							    MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
									 
									TLorentzVector Zcandidate;
									if(indicesElPlus.size()==2 && indicesElMinus.size()==1)
									{
									  TLorentzVector Zcandidate_opt1 = (selectedElectrons[indicesElPlus[0]] + selectedElectrons[indicesElMinus[0]]);
										TLorentzVector Zcandidate_opt2 = (selectedElectrons[indicesElPlus[1]] + selectedElectrons[indicesElMinus[0]]);
									  float massleptons_opt1 = Zcandidate_opt1.M(); 
									  float massleptons_opt2 = Zcandidate_opt2.M(); 
										if(massleptons_opt1 < massleptons_opt2)
										{
										   massleptons = massleptons_opt1;
											 Zcandidate = Zcandidate_opt1;
										}
										else
										{
										   massleptons = massleptons_opt2;
											 Zcandidate = Zcandidate_opt2;
										}
									}
									else if(indicesElPlus.size()==1 && indicesElMinus.size()==2)
									{
									  TLorentzVector Zcandidate_opt1 = (selectedElectrons[indicesElPlus[0]] + selectedElectrons[indicesElMinus[0]]);
										TLorentzVector Zcandidate_opt2 = (selectedElectrons[indicesElPlus[0]] + selectedElectrons[indicesElMinus[1]]);
									  float massleptons_opt1 = Zcandidate_opt1.M(); 
									  float massleptons_opt2 = Zcandidate_opt2.M(); 
										if(massleptons_opt1 < massleptons_opt2)
										{
										   massleptons = massleptons_opt1;
											 Zcandidate = Zcandidate_opt1;
										}
										else
										{
										   massleptons = massleptons_opt2;
											 Zcandidate = Zcandidate_opt2;
										}
									}								
									
								  if(fabs(Zmass - massleptons) < Zmasswindow)
									{
									  MSPlot["MS_PtZ_preboxWZqq_Zwindow_el"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityEl*scaleFactor);
									  if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
										{
										   MSPlot["MS_nEvts_boxWZqq_Zwindow_ElElEl_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
									     MSPlot["MS_nEvts_boxWZqq_Zwindow_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
								       MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
										}
									}
							  }
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
								if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						    {
								  //events with 3 or 4 same-charge leptons not compatible with ZZ signal
								  if(indicesElPlus.size()==2 && indicesElMinus.size()==2)
								  {
								   if(verbosityBox) cout<<"category: ZZqq (el)"<<endl;
								   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			             {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							     }
					         MSPlot["MS_nEvts_boxZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
							     MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
									 
									 TLorentzVector Zcandidate1, Zcandidate2;								 
								   TLorentzVector Zcandidate1_opt1 = (selectedElectrons[indicesElPlus[0]] + selectedElectrons[indicesElMinus[0]]);
								   TLorentzVector Zcandidate1_opt2 = (selectedElectrons[indicesElPlus[1]] + selectedElectrons[indicesElMinus[1]]);
								   TLorentzVector Zcandidate2_opt1 = (selectedElectrons[indicesElPlus[0]] + selectedElectrons[indicesElMinus[1]]);
								   TLorentzVector Zcandidate2_opt2 = (selectedElectrons[indicesElPlus[1]] + selectedElectrons[indicesElMinus[0]]);
								 
								   float massleptons2 = 9999;
								   float massleptons_opt1 = Zcandidate1_opt1.M();
								   float massleptons2_opt1 = Zcandidate2_opt1.M();
								   float massleptons_opt2 = Zcandidate1_opt2.M();
								   float massleptons2_opt2 = Zcandidate2_opt2.M();
									
									 if(massleptons_opt1 < massleptons_opt2)
								   {
								      Zcandidate1 = Zcandidate1_opt1;
								      massleptons = massleptons_opt1;
								   }
								   else
								   {
								      Zcandidate1 = Zcandidate1_opt2;
								      massleptons = massleptons_opt2;
								   }									
										
								   if(massleptons2_opt1 < massleptons2_opt2)
								   {
								      Zcandidate2 = Zcandidate2_opt1;
								      massleptons2 = massleptons2_opt1;
								   }
								   else
							  	 {
								      Zcandidate2 = Zcandidate2_opt2;
								      massleptons2 = massleptons2_opt2;		
							  	 }								
															
								   if(fabs(Zmass - massleptons) < Zmasswindow && fabs(Zmass - massleptons2) < Zmasswindow)
								   {
								      if((doZptcut && Zcandidate1.Pt() > Zptcut && Zcandidate2.Pt() > Zptcut) || !doZptcut)
										  {
								         MSPlot["MS_nEvts_boxZZqq_Zwindow_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
								         MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
										  }
								   }
									 
									}
									 
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
					     if(!doMETcutboxWWqq || (doMETcutboxWWqq && met>METcutboxWWqq))
						   {
							  if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						    {
								  //lepton charges should be opposite
									if((indicesMuPlus.size()==1 && indicesElMinus.size()==1) || (indicesMuMinus.size()==1 && indicesElPlus.size()==1))
									{
							      if(verbosityBox) cout<<"category: WWqq (el)"<<endl; //(no need to check the invariant mass of the leptons, they should not come from a Z...)
							      if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			              {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWWqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							      }
					          MSPlot["MS_St_boxWWqq_el"]->Fill(ST,datasets[d], true, LuminosityEl*scaleFactor);
										//for 'flat' binning
										float fillweight = datasets[d]->NormFactor () * LuminosityEl * scaleFactor;
										if(makeflatbinning) VarWWqq_el.Fill_forBinning(ST,fillweight);
					          else if(useflatbinning) VarWWqq_el.Fill(ST,fillweight,dataSetName);
						      }
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
					    if((doMETcutboxWZqq && met > METcutboxWZqq) || !doMETcutboxWZqq)
						  {			     
							 if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						   {
							  if(indicesElPlus.size()==1 && indicesElMinus.size()==1)
								{
							    if(verbosityBox) cout<<"category: WZqq (el)"<<endl;
							    if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			            {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							    }								 
									
									TLorentzVector Zcandidate = (selectedElectrons[0] + selectedElectrons[1]);
					        MSPlot["MS_nEvts_boxWZqq_ElElMu_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
									MSPlot["MS_nEvts_boxWZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						      MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);									
									massleptons = (selectedElectrons[indicesElPlus[0]] + selectedElectrons[indicesElMinus[0]]).M();
								  if(fabs(Zmass - massleptons) < Zmasswindow)
								  {
									    MSPlot["MS_PtZ_preboxWZqq_Zwindow_el"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityEl*scaleFactor);
									    if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
								      {
											   MSPlot["MS_nEvts_boxWZqq_Zwindow_ElElMu_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
											   MSPlot["MS_nEvts_boxWZqq_Zwindow_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
								         MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
									    }
									}
								}
						   }
							}
					  }
					}
				}
				else
				{
				   	cout<<"FOUND >= (incompatible) 4-lepton or 5-lepton event!"<<endl;			
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
						 if((doMETcutboxWZqq && met > METcutboxWZqq) || !doMETcutboxWZqq)
						 {				    
							if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						  {
							  if(indicesMuPlus.size()==1 && indicesMuMinus.size()==1)
								{
							   if(verbosityBox) cout<<"category: WZqq (el)"<<endl;
							   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			           {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxWZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							   }
								 MSPlot["MS_nEvts_boxWZqq_ElMuMu_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);	
					       MSPlot["MS_nEvts_boxWZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);	
						     MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
								 
								 massleptons = (selectedMuons[indicesMuPlus[0]] + selectedMuons[indicesMuMinus[0]]).M();
								 if(fabs(Zmass - massleptons) < Zmasswindow)
								 {
										    TLorentzVector Zcandidate = (selectedMuons[0] + selectedMuons[1]);
												MSPlot["MS_PtZ_preboxWZqq_Zwindow_el"]->Fill(Zcandidate.Pt(),datasets[d], true, LuminosityEl*scaleFactor);
												if((doZptcut && Zcandidate.Pt() > Zptcut) || !doZptcut)
												{
												   MSPlot["MS_nEvts_boxWZqq_Zwindow_ElMuMu_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
												   MSPlot["MS_nEvts_boxWZqq_Zwindow_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
								           MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);													 
												}
								 }
								
								}
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
							if(!doAntiBtagging || (doAntiBtagging && (selectedJets.size()==selectedAntiBtaggedJets.size())))
						  {
							 if((indicesMuPlus.size()==1 && indicesMuMinus.size()==1) && (indicesElPlus.size()==1 && indicesElMinus.size()==1))
							 {
							   if(verbosityBox) cout<<"category: ZZqq (el)"<<endl;
							   if(doAntiBtagging && applyBtagSF && useBtagInfoFromFile && !((dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) || (dataSetName.find("InvIso") != string::npos)))
			           {
							          BtagMCWeight = bTool_L->getMCEventWeight(selectedJets,selectedJets_partonFlavour,selectedJets_bTagCSV,syst_btag,syst_mistag);
								        //cout << "   ---> BtagMCWeight = "<<BtagMCWeight<<endl;
												MSPlot["MS_BtagWeight_boxZZqq_el"]->Fill(BtagMCWeight, datasets[d], true, LuminosityEl*scaleFactor);
								        scaleFactor = 	scaleFactor*BtagMCWeight; 
							   }
					       MSPlot["MS_nEvts_boxZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
						     MSPlot["MS_nEvts_boxWZqqZZqq_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
								 
								 float massleptons2 = 9999;
								 //probably VERY few events surviving
								 if((indicesMuPlus.size()==1 && indicesMuMinus.size()==1) && (indicesElPlus.size()==1 && indicesElMinus.size()==1))
								 {
									  TLorentzVector Zcandidate1 = (selectedMuons[indicesMuPlus[0]] + selectedMuons[indicesMuMinus[0]]);
									  TLorentzVector Zcandidate2 = (selectedElectrons[indicesElPlus[0]] + selectedElectrons[indicesElMinus[0]]);
										massleptons = Zcandidate1.M();
										massleptons2 = Zcandidate2.M();
										if(fabs(Zmass - massleptons) < Zmasswindow && fabs(Zmass - massleptons2) < Zmasswindow)
										{
										   if((doZptcut && Zcandidate1.Pt() > Zptcut && Zcandidate2.Pt() > Zptcut) || !doZptcut)
											 {
										     MSPlot["MS_nEvts_boxZZqq_Zwindow_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
								         MSPlot["MS_nEvts_boxWZqqZZqq_Zwindow_el"]->Fill(1,datasets[d], true, LuminosityEl*scaleFactor);
											 }
										}
								 }
							 }
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
		 
	   //delete selection;
     //cout << "done processing event" << endl;
		} //loop on events

    
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
				temp->showNumberEntries(showEntriesLegend);
				temp->setPreliminary(setCMSPrelim);
				if(!foundMu && !foundEl) //warning, to be considered: what when you only run on e.g. mu data but not electron data
				{
				    if(hasEnding(name,"_mu")) temp->setDataLumi(LuminosityMu); //if not run on data...
	          if(hasEnding(name,"_el")) temp->setDataLumi(LuminosityEl); //if not run on data...
        }
				//temp->Draw(false, name, true, true, true, true, true,5,false, true, true, true, true);// old way: (bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal, bool addRatio, bool mergeVV, bool mergeTTV, bool mergeVVV, bool mergeSameSignWW)
        temp->setErrorBandFile("VLQSearch_ErrorBand/ErrorBandFile_7Jan14.root");
				//temp->Draw(name,2,true,true,false,5); //label, RatioType, addRatioErrorBand, addErrorBand, ErrorBandAroundTotalInput, scaleNPSignal
				temp->Draw(name,RatioType, addRatioErrorBand, addErrorBand, ErrorBandAroundTotalInput, scaleNPSignal); 
				//temp->Draw(name,0,false,false,false,5); //for inviso
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
    fout->cd();
    
    //Selection tables
/*    selecTableMu.TableCalculator(true, true, true, true, true, true, true,false);//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST, bool mergeVV, bool mergettV, bool NP_mass)
    string selectiontableMu = "VLQSearch_SelectionTable_TreeAnalyzer_Mu";
    selectiontableMu = Outputpath+selectiontableMu+outputpostfix+".tex"; 	
    selecTableMu.Write(selectiontableMu.c_str(),false, true, false, false, false, false, false); //(filename, error, merged, lines, unscaled, eff, totaleff, landscape)
*/
   /* 
		selecTableMultiLep.TableCalculator(true, true, true, true, true, false, true, true);
    string selectiontableMultiLep = "InclFourthGenSearch_SelectionTable_TreeAnalyzer_MultiLepton";
    selectiontableMultiLep = Outputpath+selectiontableMultiLep+outputpostfix+".tex"; 	
    selecTableMultiLep.Write(selectiontableMultiLep.c_str(),false, true, false, false, false, false, false);
		*/
		
		//for 'flat' binning
    if (makeflatbinning)
    {	
		    //mu channel
        VarWqq_mu.Write_forBinning(binningFileName_Wqq_mu);
				VarWplusqq_mu.Write_forBinning(binningFileName_Wplusqq_mu);
				VarWminusqq_mu.Write_forBinning(binningFileName_Wminusqq_mu);
				VarWHqq_mu.Write_forBinning(binningFileName_WHqq_mu);
			  VarZqq_mu.Write_forBinning(binningFileName_Zqq_mu);
			  VarZHqq_mu.Write_forBinning(binningFileName_ZHqq_mu);
			  VarSemiLepZZqqWZqq_mu.Write_forBinning(binningFileName_SemiLepZZqqWZqq_mu);
			  VarWWqq_mu.Write_forBinning(binningFileName_WWqq_mu);
				//el channel
				VarWqq_el.Write_forBinning(binningFileName_Wqq_el);
				VarWplusqq_el.Write_forBinning(binningFileName_Wplusqq_el);
				VarWminusqq_el.Write_forBinning(binningFileName_Wminusqq_el);
				VarWHqq_el.Write_forBinning(binningFileName_WHqq_el);
			  VarZqq_el.Write_forBinning(binningFileName_Zqq_el);
			  VarZHqq_el.Write_forBinning(binningFileName_ZHqq_el);
			  VarSemiLepZZqqWZqq_el.Write_forBinning(binningFileName_SemiLepZZqqWZqq_el);
			  VarWWqq_el.Write_forBinning(binningFileName_WWqq_el);
    }
		else if (useflatbinning)
		{
		    //mu channel
		    VarWqq_mu.Write(fout,th1dir);
				VarWplusqq_mu.Write(fout,th1dir);
				VarWminusqq_mu.Write(fout,th1dir);
				VarWHqq_mu.Write(fout,th1dir);
			  VarZqq_mu.Write(fout,th1dir);
			  VarZHqq_mu.Write(fout,th1dir);
			  VarSemiLepZZqqWZqq_mu.Write(fout,th1dir);
			  VarWWqq_mu.Write(fout,th1dir);
				//el channel
				VarWqq_el.Write(fout,th1dir);
				VarWplusqq_el.Write(fout,th1dir);
				VarWminusqq_el.Write(fout,th1dir);
				VarWHqq_el.Write(fout,th1dir);
			  VarZqq_el.Write(fout,th1dir);
			  VarZHqq_el.Write(fout,th1dir);
			  VarSemiLepZZqqWZqq_el.Write(fout,th1dir);
			  VarWWqq_el.Write(fout,th1dir);				
		}
		 
		 
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
