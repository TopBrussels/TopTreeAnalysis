// isis.marina.van.parijs@cern.ch
//  based on SingleTop-tW.cc of rebeca  
//
// Summer12 round
// 8 Tev

// This is a program that runs over the toptrees

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/JetTools.h"
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
#include "../MCInformation/interface/ResolutionFit.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../MCInformation/interface/LumiReWeighting.h"
#include "../macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

int main(int argc, char* argv[]) {
    
     
    /////////////////////////////////////////////
    ///               SYSTEMATICS             ///
    /////////////////////////////////////////////
    
    bool JESPlus = false; 
    bool JESMinus = false; 
    
    bool JERPlus = false; 
    bool JERMinus = false; 
    
    bool SFplus = false; 
    bool SFminus = false; 
    
    bool unclusteredUp = false; 
    bool unclusteredDown = false; 
    
    bool PUsysUp = false; 
    bool PUsysDown = false; 
    
    /////////////////////////////////////////////
    ///                 B-tag SF              ///
    /////////////////////////////////////////////
    bool scaleFactor = false; 

    /////////////////////////////////////////////
    ///              Naked Option             ///
    /////////////////////////////////////////////
    bool isRAW = false; 

    /////////////////////////////////////////////
    ///                  Run HLT              ///
    /////////////////////////////////////////////    
    bool runHLT = false; 
    
    /////////////////////////////////////////////
    ///                 PU reweighting        ///
    /////////////////////////////////////////////    
    
    // use by default the official PU reweighting, so not the 3d kind of 2011
    bool reweightPU = true; 
    bool Pu3D = false;    
    
    
    
    
    /////////////////////////////////////////////
    ///       SINGLE TOP TW ANALYSIS          ///
    /////////////////////////////////////////////
    
    //Make the plots nicer: colors, style,... 
    setMyStyle(); 
    
    // To see how long the program takes to run, I include a clock
    clock_t start = clock(); 
    
    cout << "******************************************************************************" << endl; 
    cout << " Welcome to the single top tW analysis: Summer12 53X at 8 TeV with MC " << endl;   
    cout << "******************************************************************************" << endl; 
    
    // Definition of the different working modes: 
    //     0 = emu 
    //     1 = mumu 
    //     2 = ee 
    
    int mode = 0; 
    double lumi = 1000; //by default is the luminosity 1000
    string xmlfile = "twemu.xml" ;  // by default, use mode 0
    
    // Set as default the emu mode
    if(mode != 0 && mode != 1 && mode != 2){
    	mode = 0; 
    }
    
    

    // Give the different options for executing 
    // argc is the number of arguments you give, your name of the code is already 1 argument, if you put space --something then you have 2 arguments so if argc >1, something happens
    // argv is the array that contains all the possible arguments
    // argval is the argument on place iarg in the array argv
    
    for(int iarg = 0; iarg < argc && argc>1 ; iarg++){
    	std::string argval=argv[iarg];
    	if(argval=="--help" || argval =="--h"){
      		cout << "--ee:   Di-Electron" << endl;
      		cout << "--emu:  Electron-Muon" << endl;
      		cout << "--mumu: Di-Muon" << endl;
      		cout << "--JESplus: JES sys +1 sigma MET included" << endl;
      		cout << "--JESminus: JES sys -1 sigma MET included" << endl;
      		cout << "--JERplus: JER +" << endl;
      		cout << "--JERminus: JER -" << endl;
      		cout << "--SFplus: SF up +10% syst" << endl;
      		cout << "--SFminus: SF down -10% syst" << endl;
      		cout << "--PUup: PU reweghting scaled up " << endl;
      		cout << "--PUdown: PU reweghting scaled down " << endl;
      		cout << "--uncMETup: Unclustered MET syst. Up " << endl;
      		cout << "--uncMETdown: Unclustered MET syst. Down " << endl;
      		cout << "--NoPU: Do not apply pileup re-weighting" << endl;
      		cout << "--NoSF: Do not apply b-tag scale factor" << endl;
      		cout << "--RAW: Do not apply pileup re-weighting or b-tag scale factor" << endl;
      		cout << "--3D: 3D Pileup reweighting" << endl;
      		return 0;
    	}
    	if (argval=="--ee"){
     		mode = 2;
    	}
    	if (argval=="--emu"){
     		mode = 0;
    	}
    	if (argval=="--mumu"){
    		mode = 1;
    	}
    	if (argval=="--uncMETup"){
       		unclusteredUp = true;
    	}
    	if (argval=="--uncMETdown"){
       		unclusteredDown = true;
    	}
    	if (argval=="--PUup" ){
    		PUsysUp = true;
    	}
    	if (argval=="--PUdown" ){
    		PUsysDown = true;
	}
    	if (argval=="--JESplus") {
    		JESPlus = true;
	}
    	if (argval=="--JESminus") {
    		JESMinus = true;
	}
    	if (argval=="--JERplus") {
    		JERPlus = true;
	}
    	if (argval=="--JERminus") {
    		JERMinus = true;
	}
    	if (argval=="--SFplus") {
    		SFplus = true;
	}
    	if (argval=="--SFminus"){
    		SFminus = true;
	}
    	if (argval=="--NoPU") {
    		reweightPU = false;
	}
    	if (argval=="--NoSF") {
    		scaleFactor = false;
	}
    	if (argval=="--RAW") {
    		reweightPU = false; 
		scaleFactor = false; 
		isRAW = true;
	}
    	if (argval=="--3D") {
    		Pu3D = true;
	}
    
    }   
    
    // Give the modes on the screen 
    cout << "You have chosen: ";
    if(mode == 0){
    	cout << " Electron - Muon channel " << endl;
    }
    if(mode == 1){
    	cout << " Di-Electron channel " << endl;
    }
    if(mode == 2){
    	cout << " Di-Muon channel " << endl;
    }
    cout << "******************************************************************************" << endl; 
    
    
    /////////////////////////////////
    ///       CONFIGURATION       ///
    /////////////////////////////////
    
    // xml file, this file contains the rootfiles of the skimmed toptrees
    if( mode == 0){
    	lumi = 5085.246;   // Still to check! 
	xmlfile = "twemu.xml";
    }
    if( mode == 1){
    	lumi = 1000;        // Still to check! 
	xmlfile = "twmumu.xml";
    }
    if( mode == 2){
    	lumi = 5103.58; // Still to check! 
	xmlfile = "twee.xml";
    }
    
    
    // Configuration of the output format this takes the configurations of the xml file 
       // Make a top tree with the name: set configurations, you can call this tree with configTree
       //   A TTree object has a header with a name and a title. It consists of a list of independent branches (TBranch). Each branch has its own definition and list of buffers.
    TTree *configTree = new TTree("configTree", "configuration tree"); 
    
       // Make an array of identical objects that expands automatically when objects are addded, the lower bound is 1000, and you call this with Datasets
    TClonesArray* tcdatasets = new TClonesArray("Datasets", 1000);
    
       //Make a branch in the toptrees
       // For every element of Datasets a toplevel is created. If this element is an collection on its own, then for every element of this collection a toplevel is created, etc this continues untill the number &tcdatasets is zero. (it decreases with one
       //  everytime elements are taken. The class of datasets is TClonesArray at adress &tcdatasets  (IS DIT CORRECT ???)
    configTree->Branch("Datasets","TClonesArray", &tcdatasets);
   
       // Make another array
    TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000); 
    
       // Make a branch for previous array 
    configTree->Branch("AnalysisEnvironment", "TClonesArray", &tcAnaEnv); 
    
    
      
    /////////////////////////////////
    ///    Analysis environment   ///     
    /////////////////////////////////    
    
    //Create an analysisenvironment
    AnalysisEnvironment anaEnv; 
    AnalysisEnvironmentLoader anaLoad(anaEnv, xmlfile.c_str()); 
    new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv); 
    
    
    /////////////////////////////////
    ///  LOADING OF THE DATASETS  ///  
    /////////////////////////////////    
    // Load a toptree
    TTreeLoader treeLoader; 
    
    //make a vector containing all the datasets
    vector <Dataset*> datasets;
    
    //Load the xmlfiles in dataset, each sample in the xml file is an element of the vector datasets (IS THIS CORRECT ??? )
    treeLoader.LoadDatasets (datasets, xmlfile.c_str());
    
    
    
    /////////////////////////////////
    ///  START OF THE ANALYSIS    ///  
    ///////////////////////////////// 
    
    // Do the analysis for every dataset in the vector with name datasets containing the datasets as elements
    for(unsigned int d = 0; d < datasets.size(); d++){
    	// Take the dataset on place d in the vector
	treeLoader.LoadDataset (datasets[d], anaEnv); 
	
	// Take the name of the chosen dataset
	string dataSetName = datasets[d]->Name(); 	
	
	/////////////////////////////////
	///        DECLARATIONS       ///
	///////////////////////////////// 
	bool isData = false;    // To make the division between data an MC
	bool isTop = false;     // To make the division between top pair MC and single top MC, this is needed for the btagging SF
	double xlweight;        // To define the reweighting of the MC compared to the data, if this is 1, then no reweighting is applied. 
	                        // The reweighting is done as follows: 
				//       xlweight = (cross-section x luminosity)/number of events in the toptree before the skimming
				//This number of events is given in the mail that you receive with the urls of the skimmed toptrees (so you have to save these numbers)
	
	
	// Define the cross sections and weights for every data set
	// sprintf makes the string you call name contain data 
	char name[100];
	
	if (dataSetName == "data"){		sprintf(name, "data");  	xlweight = 1; 				isData = true;}
        else if (dataSetName == "tt"){          sprintf(name, "tt");            xlweight = lumi*225.197/6709118; 	isTop = true;} 
        else if (dataSetName == "twdr"){        sprintf(name, "tw_dr");         xlweight = lumi*11.1/497657; 		} 
        else if (dataSetName == "atwdr"){       sprintf(name, "atw_dr");        xlweight = lumi*11.1/493460; 		} 
        else if (dataSetName == "t"){         	sprintf(name, "t");         	xlweight = lumi*56.4/23777; 		} 
        else if (dataSetName == "at"){         	sprintf(name, "at");         	xlweight = lumi*30.7/1935071; 		} 
        else if (dataSetName == "ww"){         	sprintf(name, "ww");         	xlweight = lumi*57.07/9969958; 		} 
        else if (dataSetName == "wz"){         	sprintf(name, "wz");         	xlweight = lumi*22.44/8080197; 		} 
        else if (dataSetName == "zz"){         	sprintf(name, "zz");         	xlweight = lumi*9.03/9799902; 		} 
        else if (dataSetName == "zjets"){       sprintf(name, "zjets");         xlweight = lumi*3532.8/16080506; 	} 
        else if (dataSetName == "zjets_lowmll"){sprintf(name, "zjets_lowmll");  xlweight = lumi*860.5/7132214; 	        } 
        else if (dataSetName == "wjets"){  	sprintf(name, "wjets");  	xlweight = lumi*37509/18036994; 	}  
	
	// Define the output rootfiles 
	//  ==>  sprintf(rootFileName,"outputs/naked_%d_%s.root", mode, name)
        //       This makes the rootfile to be called outputs/naked_modeName_sampleName.root
	//       eg: if you are looking at emu (0) and data, it will be called naked_0_data.root and placed in the directory 
	char rootFileName[100];
      
        if  (isRAW){                         sprintf(rootFileName,"outputs/naked_%d_%s.root", mode, name);}
        else if (!isData && JESPlus){        sprintf(rootFileName,"outputs/JESsysUp_%d_%s.root", mode, name);}
        else if (!isData && JESMinus){       sprintf(rootFileName,"outputs/JESsysDown_%d_%s.root", mode, name);}
        else if (!isData && JERMinus){       sprintf(rootFileName,"outputs/JERsysDown_%d_%s.root", mode, name);}
        else if (!isData && JERPlus){        sprintf(rootFileName,"outputs/JERsysUp_%d_%s.root", mode, name);}
        else if (!isData && SFplus){         sprintf(rootFileName,"outputs/SFsysUp_%d_%s.root", mode, name);}
        else if (!isData && SFminus){        sprintf(rootFileName,"outputs/SFsysDown_%d_%s.root", mode, name);}
        else if (!isData && unclusteredUp){  sprintf(rootFileName,"outputs/METsysUp_%d_%s.root", mode, name);}
        else if (!isData && unclusteredDown){sprintf(rootFileName,"outputs/METsysDown_%d_%s.root", mode, name);}
        else if (!isData && PUsysUp){        sprintf(rootFileName,"outputs/PUsysUp_%d_%s.root", mode, name);}
        else if (!isData && PUsysDown){      sprintf(rootFileName,"outputs/PUsysDown_%d_%s.root", mode, name);}
        else if (!isData && !reweightPU){    sprintf(rootFileName,"outputs/out_noPU_%d_%s.root", mode, name);}
        else if (!isData && Pu3D){           sprintf(rootFileName,"outputs/out_3D_%d_%s.root", mode, name);}
        else{                                sprintf(rootFileName,"outputs/out_%d_%s.root", mode, name);}
      


        // Define information output files 
	//    ==> ofstream output(myTexFile)
    	//        This writes the output file stream to myTexFile, these files are stored in the directory information
	//     WHAT IS WRITTEN IN THESE FILES???
        char myTexFile[300];
        sprintf(myTexFile,"information/lepsel_info_run_lumi_event_%d_%s.txt", mode, name);
        ofstream OutPut(myTexFile);
        sprintf(myTexFile,"information/lepveto_run_lumi_event_%d_%s.txt", mode, name);
        ofstream OutPut2(myTexFile);
        sprintf(myTexFile,"information/met_run_lumi_event_%d_%s.txt", mode, name);
        ofstream OutPut3(myTexFile);
        sprintf(myTexFile,"information/jet_run_lumi_event_%d_%s.txt", mode, name);
        ofstream salida4(myTexFile);
        sprintf(myTexFile,"information/bt_run_lumi_event_%d_%s.txt", mode, name);
        ofstream OutPut5(myTexFile);
	
	
	// Define the objects
	// ==> vector <kind> V
	//     This defines a vector V with elements of sort kind
		
	vector <TRootVertex*>   vertex; 
	vector <TRootMuon*>     initial_muons; 
	vector <TRootElectron*> initial_electrons; 
	vector <TRootJet*>      initial_jets; 
	vector <TRootJet*>      initial_jets_corrected; 
	vector <TRootMET*>      mets; 
	vector <TRootGenJet*>   genjets;
	
	
	/////////////////////////////////
	///    Pile Up reweighting    ///     HOW DOES THIS REWEIGHTING WORK? I know you divide MC histo with data Histo for resulting normalizing factor, but how to get these histos? 
	////////////////////////////////      where are they given? 
	
	 
	LumiReWeighting LumiWeights; 
	LumiWeights = LumiReWeighting("pileupHistos/Summer12.root","pileupHistos/Run2012AB_new.root","pileup","pileup");  // THIS WRITES STUFF (2x), what is this??
	
	//systematics    WHAT DOES THIS DO? 
	reweight::PoissonMeanShifter PShiftDown_= reweight::PoissonMeanShifter(-0.6); 
	reweight::PoissonMeanShifter PShiftUp_= reweight::PoissonMeanShifter(0.6); 
	
	
	
	////////////////////////////////////
	/// INITALISATION JEC/JER Factors //    This isn't applied for JEC, but needed for JER  ==> HOW DOES THIS SECTION WORK? 
	////////////////////////////////////
	vector<JetCorrectorParameters> vCorrParam;
      
      	// Create the JetCorrectorParameter objects, the order does not matter.
      	// YYYY is the first part of the txt files: usually the global tag from which they are retrieved
      	JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L3Absolute_AK5PFchs.txt");
      	JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L2Relative_AK5PFchs.txt");
      	JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("JECFiles/Summer12_V3_MC_L1FastJet_AK5PFchs.txt");

      	//  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
      	vCorrParam.push_back(*L1JetPar);
      	vCorrParam.push_back(*L2JetPar);
      	vCorrParam.push_back(*L3JetPar);

      	if(!isData) { 
		JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Summer12_V3_DATA_L2L3Residual_AK5PFchs.txt");
		vCorrParam.push_back(*ResJetCorPar);	
      	} 
      
      	//I think this is not used!
      		JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Summer12_V3_MC_Uncertainty_AK5PFchs.txt");
      	// if (isData) JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Summer12_V3_DATA_Uncertainty_AK5PFchs.txt");
    
      	// true means redo also the L1
      		JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); 
       
	
	
	
	////////////////////////////////
	///        ROOTSTUFF          //    
	////////////////////////////////
	// Open the created rootfiles and RECREATE:  create a new file, if the file already exists it will be overwritten. 
	// rootFileName is created before with sprintf(rootFileName,"outputs/naked_%d_%s.root", mode, name), so rootfiles were created with the name outputs/naked_modeName_sampleName.root
	// eg: if you are looking at emu (0) and data, it was called naked_0_data.root and placed in the directory  outputs
	TFile *fout = new TFile (rootFileName, "RECREATE"); 
	
	//Make an event
	TRootEvent* event = 0; 
	
	
	
	////////////////////////////////
	///        HISTOGRAMS         //    
	////////////////////////////////
	map <string,TH1F> histo1D; 
	map <string, TH2F> histo2D; 
	
	// Definition of the histograms
	TH1F* cutflow = new TH1F("cutflow", "The cutflow", 31, -0.5,30.5);     // A 1 dimensional histo with reference cutflow, title "The cutflow", 31 bins, starting from -0.5 till 30.5
										// This shows the influence of every selection the number of events, reweighted
	TH1F* cutflow_raw = new TH1F("cutflow_raw", "The raw cutflow", 31,-0.5, 30.5); // same but not reweighted 
	TH1F* Regions = new TH1F ("Regions" ," The different regions", 40,0,40); // for example  according to number of jets and b-tagged jets
	
	TH1F* pileup_weights = new TH1F("pileup_weights", "The calculated PU weights", 1000,0,10); 
	TH1F* pileup_weights3D = new TH1F("pileup_weights3D", "The calculated 3DPU weights", 1000,0,10); 
	
	
	
	 //Create structure to store sum of squares of weights:   if histogram is already filled, the sum of squares of weights is filled with the existing bin contents. 
	 //							  The error per bin will be computed as sqrt(sum of squares of weight) for each bin
	cutflow->Sumw2(); 
	cutflow_raw->Sumw2(); 
	Regions->Sumw2();  
	
	pileup_weights->Sumw2(); 
	pileup_weights3D->Sumw2(); 
	
	
	
	////////////////////////////////
	///    CREATION OUTPUT TREE   //    
	////////////////////////////////
	// Variables to put as branches in the toptree
	double xlWeight; 
	double puweight; 
	double rawWeight; 
	
	double lum; 
	
	int npu; 
	int nvertex; 
	
	double metPt; 
	double metPx; 
	double metPy; 
	
	std::vector<double> *ptLepton; 
	std::vector<double> *pxLepton; 
	std::vector<double> *pyLepton; 
	std::vector<double> *pzLepton; 
	std::vector<double> *eLepton; 
	std::vector<double> *qLepton; 
	
	
	std::vector<double> *ptJet; 
	std::vector<double> *pxJet; 
	std::vector<double> *pyJet; 
	std::vector<double> *pzJet; 
	std::vector<double> *eJet; 
	std::vector<double> *qJet;  
	std::vector<double> *btJPBJet; 
	std::vector<double> *btBJPBJet; 
	std::vector<double> *btCSVBJet; 
	std::vector<double> *btCSVBmvaJet; 
	
	
	// Make the output tree 
	TTree* myTree  = new TTree("myTree", " "); 
	
	// Set the branches for the doubles
	myTree -> Branch("xlWeight", &xlWeight, "xlWeight/D");      // Make a branch with name xlWeight, on location xlWeight, with WHAT IS XlWeight/D???
	myTree -> Branch("puweight", &puweight, "puweight/D"); 
	myTree -> Branch("rawWeight", &rawWeight, "rawWeight/D"); 
	
	myTree -> Branch("lum", &lum, "lum/D");
	
	myTree -> Branch("npu", &npu, "npu/D");
	myTree -> Branch("nvertex", &nvertex, "nvertex/D");
	
	myTree -> Branch("metPt", &metPt, "metPt/D");
	myTree -> Branch("metPx", &metPx, "metPx/D");
	myTree -> Branch("metPy", &metPy, "metPy/D");
	
	// Set the branches for the vectors 
	myTree->Branch("ptLepton","std::vector<double>",&ptLepton);   // Make a branch with name ptLepton, from type vector(double), on loaction ptLepton
        myTree->Branch("pxLepton","std::vector<double>",&pxLepton);
        myTree->Branch("pyLepton","std::vector<double>",&pyLepton);
        myTree->Branch("pzLepton","std::vector<double>",&pzLepton);
        myTree->Branch("eLepton","std::vector<double>",&eLepton);
        myTree->Branch("qLepton","std::vector<double>",&qLepton);
      
        myTree->Branch("ptJet","std::vector<double>",&ptJet);
        myTree->Branch("pxJet","std::vector<double>",&pxJet);
        myTree->Branch("pyJet","std::vector<double>",&pyJet);
        myTree->Branch("pzJet","std::vector<double>",&pzJet);
        myTree->Branch("eJet","std::vector<double>",&eJet);
        myTree->Branch("qJet","std::vector<double>",&qJet);
        myTree->Branch("btJPBJet","std::vector<double>",&btJPBJet);
        myTree->Branch("btBJPBJet","std::vector<double>",&btBJPBJet);
        myTree->Branch("btCSVBJet","std::vector<double>",&btCSVBJet);
        myTree->Branch("btCSVBmvaJet","std::vector<double>",&btCSVBmvaJet);
	
	
	////////////////////////////////
	///    INFORMATION FOR USER  ///    
	////////////////////////////////
	cout << "[Info:] output rootfile named " << rootFileName << endl; 
      	cout << "[Info:] mode = " << mode << ", lumi: " <<  lumi << " pb, sample: " << name << ", base weight: " << xlweight << endl;
      
      	if (JERPlus ||JERMinus || JESPlus || JESMinus ||  SFplus || SFminus || unclusteredUp || unclusteredDown 
	  || !reweightPU || !scaleFactor || PUsysUp || PUsysDown || Pu3D) {
		cout << "[Warning:] Non-standard options, ignore if you did it conciously" << endl;
		
		if (JERPlus) cout << "[Warning:] JER systematics on, plus" << endl;
		if (JERMinus) cout << "[Warning:] JER systematics on, minus" << endl;
		if (JESPlus) cout << "[Warning:] JES systematics on, plus" << endl;
		if (JESMinus) cout << "[Warning:] JES systematics on, minus" << endl;
		if (SFplus) cout <<"[Warning:] SF up 10% " << endl;
		if (SFminus) cout <<"[Warning:]  SF down 10% " << endl;
		if (unclusteredUp) cout <<"[Warning:] unclustered MET up 10% " << endl;
		if (unclusteredDown) cout <<"[Warning:] unclustered MET down 10% " << endl;
		if (!reweightPU && !isData) cout << "[Warning:] You are NOT applying PU re-weighting " << endl;
		if (!scaleFactor && !isData) cout << "[Warning:] You are NOT applying the b-tagging SF " << endl;
		if (PUsysUp) cout <<"[Warning:] PU up " << endl;
		if (PUsysDown) cout <<"[Warning:] PU down " << endl;
	
      	} 
	else{
		cout << "[Info:] Standard setup " << endl;
	}
	
      	cout << "[Info:] " << datasets[d]->NofEvtsToRunOver() << " total events" << endl;
      	if (runHLT) cout << "[Info:] You have the HLT activated, this might be slower than the usual. " << endl;
	
	////////////////////////////////
	///    LOOP OVER THE EVENTS  ///    
	////////////////////////////////
	for(int ievent = 0; ievent < datasets[d]->NofEvtsToRunOver(); ievent++){
		if(ievent%500==0){
			std::cout << "Processing the " << ievent << "th event" << flush << "\r";        // << flush << "\r" means this line will be overwritten next time 
		}
		
			
		//Load the event from the toptree
		event = treeLoader.LoadEvent (ievent, vertex, initial_muons, initial_electrons, initial_jets_corrected, mets); 
		
		
		
		//For the jet energy resolution of the simulation, we need generated jets. So we only need to add to the MC datasets
		if(!isData){
			genjets =treeLoader.LoadGenJet(ievent, false); 
			sort(genjets.begin(),genjets.end(),HighestPt());
		}
		
		////////////////////////////////
		//APPLYING THE PU reweighting///    
		////////////////////////////////
		//Declare the reweighting
		double weight = xlweight;
		double lumiWeight = 1.0; 
		
		//The PU reweighting is only done for MC datasets 
		if(!isData){
			lumiWeight = LumiWeights.ITweight( (int) event -> nTruePU() );
			
			// If one wants the PU reweighting scaled down   ==> WHY WOULD SOMEONE WANT THIS?
			if(PUsysDown){
				lumiWeight*PShiftDown_.ShiftWeight(event ->nPu(0));     // Why not event->nTruePU() like in michael's code? 
			}
			
			// If one wants the PU reweighting scaled up
			if(PUsysUp){
				lumiWeight*PShiftUp_.ShiftWeight(event ->nPu(0)); 
			}
			
			//If PU reweighting is turned on
			if(reweightPU){
				weight *= lumiWeight;      // HOW COME THIS IS NOT WITH HISTO'S? 
			}
		} // closing loop for PU reweighting for MC
		
		////////////////////////////////
		///  JES (Jet Energy Scale)  ///     
		///        SYSTEMATICS       ///      Isn't applied, ready for use later   ==> commented
		//////////////////////////////// 
		/*
		// If the JES systematics are on, minus
		if(JESPlus){
			jetTools->correctJetJESUnc(initial_jets_corrected, "minus", 1); 
		}
		// If the JES systematics are on, plus
		else if(JESMinus){
			jetTools->correctJetJESUnc(initial_jets_corrected, "plus", 1); 
		}
		*/
		
		////////////////////////////////
		///    JER(JE Resolution)    ///     The jet energy resolution in the MC is better than the one from the real data.
		///       SMEARING           ///     Smear the MC energy resolution in order to mimic the one from data.       
		////////////////////////////////
		if(JERMinus){
			jetTools->correctJetJER(initial_jets_corrected,genjets,mets[0], "minus", false); // false means don't use old numbers, but new ones
		}
		else if(JERPlus){
			jetTools->correctJetJER(initial_jets_corrected,genjets,mets[0], "plus", false); // false means don't use old numbers, but new ones
		}
		else{
			jetTools->correctJetJER(initial_jets_corrected,genjets,mets[0], "nominal", false); // false means don't use old numbers, but new ones
		}
		
		////////////////////////////////
		///      TRIGGER SETUP       ///
		////////////////////////////////
		bool trigged = false; 
		
		//If the HLT is applied 
		if(runHLT){
			int currentRun = event->runId(); 
			bool itrigger = false; 
			bool isecondtrigger = false; 
			
			//The HLT is only used for data
			if(isData){
				//The HLT path is dependent of the mode, these paths are the several steps or software modules. Each module performs a well defined task 
				// such as reconstruction of physics objects, making intermediate decisions, triggering more refined reconstructions in subsequent modules, 
				// or calculating the final decision for that trigger path.
				if(mode == 0){
					itrigger = treeLoader.iTrigger ("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", currentRun);
					isecondtrigger = treeLoader.iTrigger ("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", currentRun);
	      			} 
				else if (mode == 1){
					itrigger = treeLoader.iTrigger ("HLT_Mu17_Mu8_v16", currentRun);
					isecondtrigger = treeLoader.iTrigger ("HLT_Mu17_TkMu8_v9", currentRun);
	      			} 
				else if (mode == 2){
					itrigger = treeLoader.iTrigger ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17", currentRun);
	      			}
	    		} // closing the HLT for data loop
			
			//For the MC, there is no triggerpath
			else {
				itrigger = true;      // WHY?
	      			isecondtrigger = true;
	    		} // closing the HLT for MC
			
			//if one of the triggers is used, the trigging is set on true   ==> WHAT IS HAPPENING? THESE STATEMENTS MAKE NO SENSE
	    		if (itrigger || isecondtrigger){ 
				trigged = true;
			 } 
			 else{
			 	trigged = true;
			}	
		} // closing the HLT run loop
		
		
		////////////////////////////////
		///        SELECTION         ///
		////////////////////////////////
		Selection selection(initial_jets_corrected, initial_muons,initial_electrons,mets); 
		
		//----------------------------------------------------------------
		//Cut on primary vertex ==> WHY IS THIS USELESS FOR SINGLE TOP? 
		//-----------------------------------------------------------------
		bool isGoodPV = false; 
		isGoodPV = selection.isPVSelected(vertex,4,24,2.);   // Look in class Selection for explanation: Here the function isPVselected is defined as
									// 	bool Selection::isPVSelected(const std::vector<TRootVertex*>& vertex, int NdofCut, float Zcut, float RhoCut){
 	 								//		if(vertex.size()==0) return false;
  									//		float Rho = sqrt(vertex[0]->x()*vertex[0]->x()+vertex[0]->y()*vertex[0]->y());
  									//		if(!vertex[0]->isFake() && vertex[0]->ndof()>NdofCut && abs(vertex[0]->z())<Zcut && Rho<RhoCut) return true;
  									//		return false;
									//	}
	
	
		//is useless for single top  ==> WHY? 
		isGoodPV = true; 
		
		//--------------------------------------------------------------
		// Set up of the unclustered MET systematic    ==> WHAT IS THIS?
		//---------------------------------------------------------------
		
		
		
		
		
		
	
	} // Closing the loop over the events in one dataset
      
    
    
    
    
	////////////////////////////////
	///   ROOTSTUFF:output files ///    
	////////////////////////////////
	// Write the recreated rootfiles 
	// rootFileName is created before with sprintf(rootFileName,"outputs/naked_%d_%s.root", mode, name), so rootfiles were created with the name outputs/naked_modeName_sampleName.root
	// eg: if you are looking at emu (0) and data, it was called naked_0_data.root and placed in the directory  outputs    
    	fout->Write(); 
	fout->Close(); 
    }  // closing the loop over the datasets 
    
    
    
    
    // To display how long it took for the program to run
    cout << "*******************************************************************************************" << endl; 
    cout << "It took you " << ((double)clock() - start) /CLOCKS_PER_SEC << " to run the program" << endl; 
  
}

