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
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
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
#include "TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/MCInformation/interface/Lumi3DReWeighting.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
#include "TopTreeAnalysis/Isis-SingleTop-tW/BtagFiles/Tools/BTagSFUtil.h"
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
    bool scaleFactor = true; 
    // load btag SF
    BTagWeightTools *bTool = new BTagWeightTools("BtagFiles/SFb-pt_payload_Moriond13.txt","CSVM");
    
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
    string xmlfile = "config/twemu.xml" ;  // by default, use mode 0
    
    // Set as default the emu mode
    if(mode != 0 && mode != 1 && mode != 2){
        mode = 0; 
    }
    
    

    // Give the different options for executing 
    // argc is the number of arguments you give, your name of the code is already 1 argument, if you put space --something then you have 2 arguments so if argc >1, something happens
    // argv is the array that contains all the possible arguments
    // argval is the argument on place iarg in the array argv
   
    std::string tempxml; 
    bool foundxml = false; 
    
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
		cout << "--xml myxml.xml Xml file" << endl; 
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
       if (argval=="--xml") {
                iarg++;
		tempxml = argv[iarg];
		foundxml = true; 
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
    
      // Luminosity and xml files
  // Only runs A, A recover, B, C (24, v2)
     if      (mode == 0){ 	 lumi = 11966.617;  	xmlfile ="config/twemu.xml";}
     else if (mode == 1){	 lumi = 12067.294;  	xmlfile = "config/twmumu.xml";}
     else if (mode == 2){	 lumi = 12093.792;  	xmlfile = "config/twee.xml";}
 
    
    // xml file, this file contains the rootfiles of the skimmed toptrees
  /*  if( mode == 0){
        lumi = 4399;   // only run B 13July2012
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
  */  
    if (foundxml){
      xmlfile = tempxml; 
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
        bool isSingleTop = false; 
	double xlweight;        // To define the reweighting of the MC compared to the data, if this is 1, then no reweighting is applied. 
                                // The reweighting is done as follows: 
                                //       xlweight = (cross-section x luminosity)/number of events in the toptree before the skimming
                                //This number of events is given in the mail that you receive with the urls of the skimmed toptrees (so you have to save these numbers)
        
        
        // Define the cross sections and weights for every data set
        // sprintf makes the string you call name contain data 
        char name[100];
	
        if (dataSetName == "data"){             sprintf(name, "data");          xlweight = 1;                           isData = true;}
        else if (dataSetName == "tt"){          sprintf(name, "tt");            xlweight = lumi*225.197/6830443;        isTop = true;} 
        else if (dataSetName == "twdr"){        sprintf(name, "tw_dr");         xlweight =lumi*11.1/497657;             isSingleTop = true;} 
        else if (dataSetName == "atwdr"){       sprintf(name, "atw_dr");        xlweight = lumi*11.1/481071;            isSingleTop = true;} 
        else if (dataSetName == "t"){           sprintf(name, "t");             xlweight = lumi*56.4/3748832;             } 
        else if (dataSetName == "at"){          sprintf(name, "at");            xlweight = lumi*30.7/180719;           }
	else if (dataSetName == "s"){           sprintf(name, "s");             xlweight = lumi*3.79/259960;             } 
        else if (dataSetName == "as"){          sprintf(name, "as");            xlweight = lumi*1.76/139974;           } 
        else if (dataSetName == "ww"){          sprintf(name, "ww");            xlweight = lumi*54.838/10000413;          } 
        else if (dataSetName == "wz"){          sprintf(name, "wz");            xlweight = lumi*22.44/9900267;          } 
        else if (dataSetName == "zz"){          sprintf(name, "zz");            xlweight = lumi*9.03/9799891 ;           } 
        else if (dataSetName == "zjets"){       sprintf(name, "zjets");         xlweight = lumi*3532.8/30364599;        } 
        else if (dataSetName == "zjets_lowmll"){sprintf(name, "zjets_lowmll");  xlweight = lumi*860.5/7059426;          } 
        else if (dataSetName == "wjets"){       sprintf(name, "wjets");         xlweight = lumi*36257.2/57411352;         }  

        
	
	
	

	
	
	
	
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
        //     In these files, is the information about: 
	//          1: the run ID   (so run A,B,C, D)
	//          2: the lumiblock ID ( so which luminosity used in that run (every run has different parts with different lumi)
	// 	    3: the eventID ( so which event in that lumiblock)
        char myTexFile[300];
        sprintf(myTexFile,"information/lepsel_info_run_lumi_event_%d_%s.txt", mode, name);
        ofstream OutPut(myTexFile);
        sprintf(myTexFile,"information/lepveto_run_lumi_event_%d_%s.txt", mode, name);
        ofstream OutPut2(myTexFile);
        sprintf(myTexFile,"information/met_run_lumi_event_%d_%s.txt", mode, name);
        ofstream OutPut3(myTexFile);
        sprintf(myTexFile,"information/jet_run_lumi_event_%d_%s.txt", mode, name);
        ofstream OutPut4(myTexFile);
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
        
	//////////////////////////////////
	//// Btagging eff             /////
	////////////////////////////////////
	std::vector<double> btag_eff_tt; 
	std::vector<double> btag_eff_twdr; 	 
 	std::vector<double> btag_eff_atwdr;  
 	std::vector<double> btag_eff_t; 
 	std::vector<double> btag_eff_at;
	std::vector<double> btag_eff_s;    
 	std::vector<double> btag_eff_as;
	std::vector<double> btag_eff_ww;  
 	std::vector<double> btag_eff_wz;
	std::vector<double> btag_eff_zz;  
 	std::vector<double> btag_eff_zjets;  
	std::vector<double> btag_eff_zjets_lowmll; 
 	std::vector<double> btag_eff_wjets;  
	
	/////////////////////////////////
	//// Mistag eff  - fakerate - non btag eff /////
	////////////////////////////////////
	std::vector<double> fake_eff_tt; 
	std::vector<double> fake_eff_twdr; 	 
 	std::vector<double> fake_eff_atwdr;  
 	std::vector<double> fake_eff_t; 
 	std::vector<double> fake_eff_at;
	std::vector<double> fake_eff_s;    
 	std::vector<double> fake_eff_as;
	std::vector<double> fake_eff_ww;  
 	std::vector<double> fake_eff_wz;
	std::vector<double> fake_eff_zz;  
 	std::vector<double> fake_eff_zjets;  
	std::vector<double> fake_eff_zjets_lowmll; 
 	std::vector<double> fake_eff_wjets; 


	
	
	
	
        /////////////////////////////////
        ///    Pile Up reweighting    /// 
	
	    
        
         
        LumiReWeighting LumiWeights; 
	
	    bool PUsysUp = false; 
    bool PUsysDown = false; 
    
      if(PUsysUp){
           LumiWeights = LumiReWeighting("pileupHistos/pileup_MC_Summer12.root","pileupHistos/pileup_tW_2012Data53X_UpToRun203002/sys_up.root", "pileup", "pileup"); 
      }else if(PUsysDown){
          LumiWeights = LumiReWeighting("pileupHistos/pileup_MC_Summer12.root","pileupHistos/pileup_tW_2012Data53X_UpToRun203002/sys_down.root", "pileup", "pileup"); 
      }else{
	  LumiWeights = LumiReWeighting("pileupHistos/pileup_MC_Summer12.root","pileupHistos/pileup_tW_2012Data53X_UpToRun203002/nominal.root", "pileup", "pileup"); 
	   }
	  
	//LumiWeights = LumiReWeighting("pileupHistos/toptree_id_2126_canvas_36.root", "pileupHistos/pileup_2012Data53X_UpToRun196531.root", "pileup", "pileup");  
       //  LumiWeights = LumiReWeighting("pileupHistos/Summer12.root","pileupHistos/Run2012AB_new.root","pileup","pileup");  // gives PU weights per bin
        /* OLD 7 TeV analysis
        //systematics   : je maakt dezelfde  distributie met mean 0.6 naar links or rechts
        reweight::PoissonMeanShifter PShiftDown_= reweight::PoissonMeanShifter(-0.6); 
        reweight::PoissonMeanShifter PShiftUp_= reweight::PoissonMeanShifter(0.6); 
        */
        
        
        ////////////////////////////////////
        /// INITALISATION JEC/JER Factors //    This isn't applied for JEC, but needed for JER  ==> HOW DOES THIS SECTION WORK? 
        ////////////////////////////////////l
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
      
        //With the new JEC correction uncertainties , you only apply these onces for systematics because the correction itself is already applied in the MC Samples themselves, 
	// each new tag comes with a new correction, in brussels they have  way of avoiding the problem that you have to skim your samples again, but in this case it isnt necessary 
	// Load the JEC corrections uncertainties                                    
                JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Summer12_V3_MC_Uncertainty_AK5PFchs.txt");
         //if (isData) JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt");
    
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
	
	TH1F* nvertex_beforePU = new TH1F("nvertex_beforePU", "The #vertices before PU reweighting", 1000,0,1000); 
        TH1F* nvertex_afterPU = new TH1F("nvertex_afterPU", "The #vertices after PU reweighting", 1000,0,1000); 
        
        
        
         //Create structure to store sum of squares of weights:   if histogram is already filled, the sum of squares of weights is filled with the existing bin contents. 
         //                                                       The error per bin will be computed as sqrt(sum of squares of weight) for each bin
        cutflow->Sumw2(); 
        cutflow_raw->Sumw2(); 
        Regions->Sumw2();  
        
        pileup_weights->Sumw2(); 
        pileup_weights3D->Sumw2(); 
        
	//////////////////////////////
	/// BTAG AND nonBTAG EFF   //	
	/////////////////////////////////



	// define histograms
 	char titlePlot[100];
	
	sprintf(titlePlot,"ptBTagB_%d_%s",mode,name);
  	TH1F* histo_ptBTagB = new TH1F(titlePlot, " pt of btagged jets for bjets", 100,  0, 200);
	
	sprintf(titlePlot,"ptBTagC_%d_%s",mode,name);
  	TH1F* histo_ptBTagC = new TH1F(titlePlot, " pt of btagged jets for cjets", 100,  0, 200);
	
  	sprintf(titlePlot,"ptBTagL_%d_%s",mode,name);
  	TH1F* histo_ptBTagL = new TH1F(titlePlot, " pt of btagged jets for ljets", 100,  0, 200);
  
  	sprintf(titlePlot,"ptBjet_%d_%s",mode,name);
  	TH1F* histo_ptBjet = new TH1F(titlePlot, " pt of jets coming from a b quark ", 100,  0, 200);
  
 	sprintf(titlePlot,"ptCjet_%d_%s",mode,name);
 	TH1F* histo_ptCjet = new TH1F(titlePlot, " pt of jets coming from a c quark ", 100 , 0, 200);
 
  	sprintf(titlePlot,"ptLjet_%d_%s",mode,name);
  	TH1F* histo_ptLjet = new TH1F(titlePlot, " pt of jets coming from light quarks ", 100 , 0, 200);   
	
	sprintf(titlePlot,"btagged_jets_%d_%s",mode,name);
  	TH1F* histo_btagged_jets = new TH1F(titlePlot, " Btagged jets",  10,  -0.5, 9.5 ); 
  	
  
 	
 	
 
  	
  	    
        
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
	std::vector<double> *etaJet;
        std::vector<double> *qJet;  
	std::vector<double> *SFjet;
        std::vector<double> *btJPBJet; 
        std::vector<double> *btBJPBJet; 
        std::vector<double> *btCSVBJet; 
        std::vector<double> *btCSVBmvaJet; 
        
        
        // Make the output tree (home made ;) ) 
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
	myTree->Branch("etaJet","std::vector<double>",&etaJet);
        myTree->Branch("qJet","std::vector<double>",&qJet);
	myTree->Branch("SFjet", "std::vector<double>",&SFjet);
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
                        
			/* OLD 7 tEV 
                        // If one wants the PU reweighting scaled down   ==> WHY WOULD SOMEONE WANT THIS?
                        if(PUsysDown){
                                lumiWeight*PShiftDown_.ShiftWeight(event ->nPu(0));     // Why not event->nTruePU() like in michael's code? 
                        }
                        
                        // If one wants the PU reweighting scaled up
                        if(PUsysUp){
                                lumiWeight*PShiftUp_.ShiftWeight(event ->nPu(0)); 
                        }
                        */
			
                        //If PU reweighting is turned on
                        if(reweightPU){
                                weight *= lumiWeight;      
                        }
                } // closing loop for PU reweighting for MC
                
                ////////////////////////////////
                ///  JES (Jet Energy Scale)  ///     
                ///        SYSTEMATICS       ///      
                //////////////////////////////// 
                
                // If the JES systematics are on, minus
                if(JESPlus){
                        jetTools->correctJetJESUnc(initial_jets_corrected, "minus", 1); 
                }
                // If the JES systematics are on, plus
                else if(JESMinus){
                        jetTools->correctJetJESUnc(initial_jets_corrected, "plus", 1); 
                }
                
                
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
                        
                        //if one of the triggers is used, the trigging is set on true   ==> WHAT IS HAPPENING? THESE STATEMENTS MAKE NO SENSE, what is this trigged? 
			// only for trigged = true, are cuts made
                        if (itrigger || isecondtrigger){ 
                                trigged = true;
                         } 
                         else{
                                trigged = true;
                        }       
                } // closing the HLT run loop
                else{               
			trigged = true; 
		}// HLT makes no difference (CHECK PREVIOUS LOOP)
		
                ////////////////////////////////
                /// SELECTION & CUTFLOW      ///
                ////////////////////////////////
                Selection selection(initial_jets_corrected, initial_muons,initial_electrons,mets); 
                
                //----------------------------------------------------------------
                //Cut on primary vertex ==> WHY IS THIS USELESS FOR SINGLE TOP? 
                //-----------------------------------------------------------------
                bool isGoodPV = false; 
                isGoodPV = selection.isPVSelected(vertex,4,24,2.);   // Look in class Selection for explanation: Here the function isPVselected is defined as
                                                                        //      bool Selection::isPVSelected(const std::vector<TRootVertex*>& vertex, int NdofCut, float Zcut, float RhoCut){
                                                                        //              if(vertex.size()==0) return false;
                                                                        //              float Rho = sqrt(vertex[0]->x()*vertex[0]->x()+vertex[0]->y()*vertex[0]->y());
                                                                        //              if(!vertex[0]->isFake() && vertex[0]->ndof()>NdofCut && abs(vertex[0]->z())<Zcut && Rho<RhoCut) return true;
                                                                        //              return false;
                                                                        //      }
        
        
                //is useless for single top  ==> WHY? 
                isGoodPV = true; 
                
                //--------------------------------------------------------------
                // Set up of the unclustered MET systematic    ==> WHAT IS THIS?
                //---------------------------------------------------------------
	  	double uncmet_px = mets[0]->Px();
	  	double uncmet_py = mets[0]->Py();
	  	for(unsigned int i=0; i<initial_jets.size(); i++){
	    		uncmet_px += initial_jets[i]->Px();
	    		uncmet_py += initial_jets[i]->Py();
	  	}
	  	for(unsigned int i=0; i<initial_muons.size(); i++){
	    		uncmet_px += initial_muons[i]->Px();
	    		uncmet_py += initial_muons[i]->Py();
	  	}	
	  	for(unsigned int i=0; i<initial_electrons.size(); i++){
	    		uncmet_px += initial_electrons[i]->Px();
	    		uncmet_py += initial_electrons[i]->Py();
	  	}	
	    
	  	double met_px = mets[0]->Px();
	  	double met_py = mets[0]->Py();
	    
	  	if(unclusteredUp){
	    		met_px += uncmet_px*0.1;
	    		met_py += uncmet_py*0.1;
	  	} 
		if(unclusteredDown){
	    		met_px -= uncmet_px*0.1;
	    		met_py -= uncmet_py*0.1;
	  	}
	    
	  	double met_pt = sqrt(met_px*met_px + met_py*met_py);
                
                //--------------------------------------------------------------
                // START OF CUTFLOW
                //---------------------------------------------------------------
		
		// Fill the histograms cutflow and cutflow_raw in the first bin with the number of events before any cut
                cutflow->Fill(1, weight);
	  	cutflow_raw->Fill(1);
		
		//If the trigged , the cutflow is started
	  	if(trigged){
			//fill the histos after triggering
	    		cutflow->Fill(2, weight);
	    		cutflow_raw->Fill(2);
			
			//If there is good primary vertex==> MAKES NO DIFFERENCE BECAUSE ISN'T APPLIED
	    		if(isGoodPV){
				//fill the histos after the goodPV criteria 
	      			cutflow->Fill(3, weight);
	      			cutflow_raw->Fill(3);
	
	      			// Select Objects -> Cuts
				// From the class Selection the following functions are used: 
				// 	void Selection::setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon)
				//	void Selection::setDiElectronCuts(float Et, float Eta, float RelIso, float d0, float MVAId, float DistVzPVz, float DRJets, int MaxMissingHits) 
				//	void Selection::setLooseDiElectronCuts(float Et, float Eta, float RelIso) 
				//	void Selection::setDiMuonCuts(float Pt, float Eta, float RelIso, float d0) 
				// 	void Selection::setLooseMuonCuts(float Pt, float Eta, float RelIso) 
				//
	      			selection.setJetCuts(20.,5.,0.01,1.,0.98,0.3,0.1);
             		 	selection.setDiMuonCuts(20.,2.4,0.20,999.);
              			selection.setDiElectronCuts(20.,2.5,0.15,0.04,0.,1,0.3,1);
              			selection.setLooseMuonCuts(10.,2.5,0.2);
              			selection.setLooseDiElectronCuts(15.0,2.5,0.2,0.5); 
		
		
		
	      			//Select Objects and put them in a vector 
	      			vector<TRootJet*> selectedJets = selection.GetSelectedJets(true);
	      			vector<TRootMuon*> selectedMuons = selection.GetSelectedDiMuons();
	      			vector<TRootMuon*> looseMuons = selection.GetSelectedLooseMuons();
	      			vector<TRootElectron*> selectedElectrons = selection.GetSelectedDiElectrons();
	      			vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseDiElectrons();
				
					
				
	      			// Tight lepton selection: make three modes with a very tight cut
	      			bool leptonSelection = false;
	      			if 		(mode == 0 && selectedElectrons.size()== 1 && selectedMuons.size()== 1) leptonSelection = true;
	      			else if 	(mode == 1 && selectedElectrons.size()== 0 && selectedMuons.size()== 2) leptonSelection = true;
	      			else if 	(mode == 2 && selectedElectrons.size()== 2 && selectedMuons.size()== 0) leptonSelection = true;

				// With this tight lepton selection, define the lepton charges and leptons
	      			if (leptonSelection) {
		  
					bool charge = false;
					double q0, q1;
					TLorentzVector lepton0, lepton1;
					
					//for the emu mode 
					if (mode == 0){
						// Take the first selected electron
		  				TRootElectron* electron = (TRootElectron*) selectedElectrons[0];
						
						// Take the first selected muon
		  				TRootMuon* muon = (TRootMuon*) selectedMuons[0];
						
						// if the leptons have an opposite sign, set the boolean on true
		  				if (electron->charge()*muon->charge() < 0){
						 	charge = true;
						}
						
						//set lepton0 as the lepton with the highest Pt
		  				if (electron->Pt() > muon->Pt()){
		    					lepton0.SetPxPyPzE(electron->Px(), electron->Py(), electron->Pz(), electron->Energy());
		    					lepton1.SetPxPyPzE(muon->Px(), muon->Py(), muon->Pz(), muon->Energy());
		    					q0 = electron->charge();
		    					q1 = muon->charge();
		  				} 
						else {
		    					lepton0.SetPxPyPzE(muon->Px(), muon->Py(), muon->Pz(), muon->Energy());
		    					lepton1.SetPxPyPzE(electron->Px(), electron->Py(), electron->Pz(), electron->Energy());
		    					q0 = muon->charge();
		    					q1 = electron->charge();
		  				}
					} 
					// for the mumu mode
					else if (mode == 1){
		  				TRootMuon* muon0 = (TRootMuon*) selectedMuons[0];
		  				TRootMuon* muon1 = (TRootMuon*) selectedMuons[1];
		  				if (muon0->charge()*muon1->charge() < 0) charge = true;
		  				lepton0.SetPxPyPzE(muon0->Px(), muon0->Py(), muon0->Pz(), muon0->Energy());
		  				lepton1.SetPxPyPzE(muon1->Px(), muon1->Py(), muon1->Pz(), muon1->Energy());
		  				q0 = muon0->charge();
		  				q1 = muon1->charge();
					} 
					// for the ee mode 
					else {
		  				TRootElectron* electron0 = (TRootElectron*) selectedElectrons[0];
		  				TRootElectron* electron1 = (TRootElectron*) selectedElectrons[1];
		  				if (electron0->charge()*electron1->charge() < 0) charge = true;
		  				lepton0.SetPxPyPzE(electron0->Px(), electron0->Py(), electron0->Pz(), electron0->Energy());
		  				lepton1.SetPxPyPzE(electron1->Px(), electron1->Py(), electron1->Pz(), electron1->Energy());
		  				q0 = electron0->charge();
		  				q1 = electron1->charge();
					}
		  			
					
					//If the event has two leptons with an opposite sign 
					if (charge){
						// fill the histos with the events with 2 leptons with an opposite sign
		  				cutflow->Fill(4, weight);
		  				cutflow_raw->Fill(4);
		  				// Loose lepton veto: there are no almost good leptons
		  				bool leptonVeto = false;
		  				if 	  	(mode == 0 && looseMuons.size()== 1 && looseElectrons.size() == 1) leptonVeto = true;
		  				else if 	(mode == 1 && looseMuons.size()== 2 && looseElectrons.size() == 0) leptonVeto = true;
		  				else if 	(mode == 2 && looseMuons.size()== 0 && looseElectrons.size() == 2) leptonVeto = true;
		  				//Write the lepton information in the previous defined output file
						// Define information output files 
        					//    ==> ofstream output(myTexFile)
        					//        This writes the output file stream to myTexFile, these files are stored in the directory information
        					//     In these files, is the information about: 
						//          1: the run ID   (so run A,B,C, D)
						//          2: the lumiblock ID ( so which luminosity used in that run (every run has different parts with different lumi)
						// 	    3: the eventID ( so which event in that lumiblock)
		 				OutPut << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;
					        
						
						//////////////////////////////////////////////////
						// CALCULATION BTAG EFFICIENCY  - fakerate     ///
						//////////////////////////////////////////////////

						double btag_eff = 0.; 
						double btag_eff_c = 0.;
						double fake_eff = 0.; 
						
						if(!isData){					
							sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //Sort the selected jets based on Pt
							
							vector <TLorentzVector> mcParticlesTLV,selectedJetsTLV;
							TLorentzVector bQuark_vector;
							
							mcParticlesTLV.clear(); //make sure nothing is inside this vector
							selectedJetsTLV.clear();


						        int pdgID; 
							for(unsigned int iJet=0;iJet<selectedJets.size(); iJet++){
	    							selectedJetsTLV.push_back(*selectedJets[iJet]);
								
								TRootJet* Jet_temp = (TRootJet*) selectedJets[iJet];
		    						pdgID = Jet_temp->partonFlavour(); 
								
								mcParticlesTLV.push_back(pdgID);
								
							} // closing for loop over selectedJets.size();
							

	
							vector< pair<unsigned int, unsigned int> > JetPartonPair;
							vector< pair<unsigned int, unsigned int> > JetPartonPair_c;
							vector< pair<unsigned int, unsigned int> > JetPartonPair_fake;
							
							
							for(unsigned int iJ=0; iJ<selectedJetsTLV.size(); iJ++) {
	    							if(fabs(pdgID) == 5)
								{
	      								JetPartonPair.push_back( pair<unsigned int, unsigned int> (pdgID, iJ) );
									
	 							 } // END JETPARTONPAIR FILLING
								if(pdgID == 4 )
								{
	      								JetPartonPair_c.push_back( pair<unsigned int, unsigned int> (pdgID, iJ) );
	 							 } // END JETPARTONPAIR FILLING
								if(pdgID == 21 || fabs(pdgID) < 4) // <5 light quarks, 21 gluon no fabs because electrically neutral
								{
	      								JetPartonPair_fake.push_back( pair<unsigned int, unsigned int> (pdgID, iJ) );
	 							 } // END JETPARTONPAIR FILLING
							} // end for loop 
							
							for(unsigned int iJ=0; iJ<JetPartonPair.size(); iJ++) {
								unsigned int jetnumber = JetPartonPair[iJ].second;
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								
								double jetPT = Jet_tempo->Pt(); 
								histo_ptBjet->Fill(jetPT);
							
							}
							
							for(unsigned int iJ=0; iJ<JetPartonPair_c.size(); iJ++) {
								unsigned int jetnumber = JetPartonPair_c[iJ].second;
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								
								double jetPT = Jet_tempo->Pt(); 
								histo_ptCjet->Fill(jetPT);
							
							}
							
							for(unsigned int iJ=0; iJ<JetPartonPair_fake.size(); iJ++) {
								unsigned int jetnumber = JetPartonPair_fake[iJ].second;
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								
								double jetPT = Jet_tempo->Pt(); 
								histo_ptLjet->Fill(jetPT);
							
							}
							
							double Number_matched_bquarks =JetPartonPair.size(); 
							double Number_matched_cquarks =JetPartonPair_c.size(); 
							double Number_matched_lquarks =JetPartonPair_fake.size();
							double Number_btagged_jets = 0; 
							double Number_btagged_jets_c = 0; 
							double Number_btagged_jets_fake = 0;
							
							for(unsigned int iPair=0; iPair<JetPartonPair.size(); iPair++) 
							{
	    							unsigned int jetnumber = JetPartonPair[iPair].second;
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								double jetPT = Jet_tempo->Pt();
								
								if (Jet_tempo->btag_combinedSecondaryVertexBJetTags()> 0.679){
			    						Number_btagged_jets++;
									histo_ptBTagB->Fill(jetPT);
									
								} // closing loop counting btagged jets
							} // closing for loop over jetpartonpair.size()
							
							for(unsigned int iPair=0; iPair<JetPartonPair_c.size(); iPair++) 
							{
	    							unsigned int jetnumber = JetPartonPair_c[iPair].second;
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								double jetPT = Jet_tempo->Pt();
								
								if (Jet_tempo->btag_combinedSecondaryVertexBJetTags()> 0.679){
			    						Number_btagged_jets_c++;
									histo_ptBTagC->Fill(jetPT);
								} // closing loop counting btagged jets
							} // closing for loop over jetpartonpair.size()
							
							for(unsigned int iPair=0; iPair<JetPartonPair_fake.size(); iPair++) 
							{
	    							unsigned int jetnumber = JetPartonPair_fake[iPair].second;
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								double jetPT = Jet_tempo->Pt();
								if (Jet_tempo->btag_combinedSecondaryVertexBJetTags()> 0.679){
			    						Number_btagged_jets_fake++;
									histo_ptBTagL->Fill(jetPT);
								} // closing loop counting btagged jets
							} // closing for loop over jetpartonpair.size()
							
							btag_eff = Number_btagged_jets/Number_matched_bquarks;
							btag_eff_c = Number_btagged_jets_c/Number_matched_cquarks;
							fake_eff = Number_btagged_jets_fake/Number_matched_lquarks;
							
												
							
							if(Number_matched_bquarks != 0){
							
		 						if (dataSetName == "tt"){ 	
								   btag_eff_tt.push_back(btag_eff ); }
 								else if (dataSetName == "twdr"){    btag_eff_twdr.push_back(btag_eff ); }	
 								else if (dataSetName == "atwdr"){   btag_eff_atwdr.push_back(btag_eff ); }	
 								else if (dataSetName == "t"){  	    btag_eff_t.push_back(btag_eff ); }
 								else if (dataSetName == "at"){ 	   btag_eff_at.push_back(btag_eff ); }
								else if (dataSetName == "s"){      btag_eff_s.push_back(btag_eff ); }	
 								else if (dataSetName == "as"){    btag_eff_as.push_back(btag_eff ); }
								else if (dataSetName == "ww"){    btag_eff_ww.push_back(btag_eff ); }	
 								else if (dataSetName == "wz"){ 	    btag_eff_wz.push_back(btag_eff ); }
								else if (dataSetName == "zz"){   btag_eff_zz.push_back(btag_eff ); }	
 								else if (dataSetName == "zjets"){     btag_eff_zjets.push_back(btag_eff ); }	
								else if (dataSetName == "zjets_lowmll"){    btag_eff_zjets_lowmll.push_back(btag_eff ); }
 								else if (dataSetName == "wjets"){   btag_eff_wjets.push_back(btag_eff ); }	
								// cout << "BTAGEFFICIENTIE CALCULATION: Btagged jets = " << Number_btagged_jets << " Number bquarks = " <<  Number_matched_bquarks << " Btag_eff = " << btag_eff << endl; 
						 	}
							
							if(Number_matched_lquarks != 0){
							
		 						if (dataSetName == "tt"){ 	
								   fake_eff_tt.push_back(fake_eff ); }
 								else if (dataSetName == "twdr"){    fake_eff_twdr.push_back(fake_eff ); }	
 								else if (dataSetName == "atwdr"){   fake_eff_atwdr.push_back(fake_eff ); }	
 								else if (dataSetName == "t"){  	    fake_eff_t.push_back(fake_eff ); }
 								else if (dataSetName == "at"){ 	   fake_eff_at.push_back(fake_eff ); }
								else if (dataSetName == "s"){      fake_eff_s.push_back(fake_eff ); }	
 								else if (dataSetName == "as"){    fake_eff_as.push_back(fake_eff ); }
								else if (dataSetName == "ww"){    fake_eff_ww.push_back(fake_eff ); }	
 								else if (dataSetName == "wz"){ 	    fake_eff_wz.push_back(fake_eff ); }
								else if (dataSetName == "zz"){   fake_eff_zz.push_back(fake_eff ); }	
 								else if (dataSetName == "zjets"){     fake_eff_zjets.push_back(fake_eff ); }	
								else if (dataSetName == "zjets_lowmll"){    fake_eff_zjets_lowmll.push_back(fake_eff ); }
 								else if (dataSetName == "wjets"){   fake_eff_wjets.push_back(fake_eff ); }	
								// cout << "NON BTAGEFFICIENTIE CALCULATION: fake Btagged jets = " << Number_btagged_jets_fake << " Number lquarks = " <<  Number_matched_lquarks << " fake_eff = " << fake_eff << endl; 
						 	}
								
						} // closing if(!isData)
						
					
	
						
						
						
						
						
						
						
						//If the event has no extra loose leptons
						if (leptonVeto) {
							// fill the histos in bin 5 with events after the loose lepton veto
		    					cutflow->Fill(5, weight);
		    					cutflow_raw->Fill(5);

		    					OutPut2 << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;
		    					
							// Low mll cut (all final states), in order to remove low invariant mass Z/gamma* events
		    					TLorentzVector pair = lepton0 + lepton1;   
		    					if (pair.M() > 20){
							       //////////////////////////////////////////////
							       ///                 TAKING BTAG SF        ////
							       //////////////////////////////////////////////
								int SFsys = 0;  
								if (SFminus) 	SFsys = -1;   // systematics down
		      						if (SFplus) 	SFsys = +1;    // systematics up
								
		      						int nJetsBT = 0;
		      						int nTightJetsBT = 0;
		     						int nJets = 0;
								int iJet = -5;

								/////////////////////////////////////////
								//  CALCULATION BTAGGING SF  adption /////// you have to correct for the fact that the btagging efficiency is not 100%
								/////////////////////////////////////////  and thus adapt your SF calculated with BtagWeight
								double bTag_SF_b;
								
								double BTagSF; 
								double BTagEff; 
								double LightJetSF; 
								double LightJetEff; 
								
								int jet_flavorSF; 
								double jet_phiSF; 
								double jet_etaSF; 
								bool bTagged = false; 
								
								
								
								for (unsigned int iJ =0; iJ < selectedJets.size(); iJ ++){
									TRootJet* tempJet = (TRootJet*) selectedJets[iJ];
									TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
									
									float TempEta = tempJet->Eta(); 
									if(fabs(TempEta)>2.4){ TempEta = 2.4;}  // temporarly fix
									
									// Use for the moment default values (from my calculations)
									if (dataSetName == "tt"){  BTagEff = 0.52;     LightJetEff = 0.31; }
									else if (dataSetName == "twdr"){   BTagEff = 0.52;     LightJetEff = 0.26; }	
									else if (dataSetName == "atwdr"){   BTagEff = 0.53;     LightJetEff = 0.25; }	
									else if (dataSetName == "t"){  	    BTagEff = 0.40;     LightJetEff = 0.14; }
									else if (dataSetName == "at"){ 	  BTagEff = 0.37;     LightJetEff = 0.14; }
									else if (dataSetName == "s"){      BTagEff = 0.54;     LightJetEff = 0.30; }	
									else if (dataSetName == "as"){     BTagEff = 0.33;     LightJetEff = 0.15; }
									else if (dataSetName == "ww"){     BTagEff = 0.40;     LightJetEff = 0.01;}	
									else if (dataSetName == "wz"){ 	     BTagEff = 0.34;     LightJetEff = 0.01; }
									else if (dataSetName == "zz"){   BTagEff = 0.48;     LightJetEff = 0.03; }	
									else if (dataSetName == "zjets"){      BTagEff = 0.40;     LightJetEff = 0.02;}	
									else if (dataSetName == "zjets_lowmll"){     BTagEff = 0.5;     LightJetEff = 0.04; }
									else if (dataSetName == "wjets"){   BTagEff = 0.5;     LightJetEff = 0.02; }	
										


									
									if(!isData){
										jet_flavorSF = tempJet->partonFlavour();
									}  

									jet_phiSF = tempJet->Phi(); 
									jet_etaSF = tempJet-> Eta(); 
									
									
									if (tempJet->Pt() > 30 && fabs(tempJet->Eta()) < 2.5 && TMath::Min(fabs(lepton0.DeltaR(tJet)),fabs(lepton1.DeltaR(tJet))) > 0.3) { // if proper jet 
			 							nJets++;
										iJet = iJ;

										
										if (tempJet->btag_combinedSecondaryVertexBJetTags()> 0.679){  // then btagged
										
											
											
											bTagged = true; 
											
											if (!isData){
												if(fabs(jet_flavorSF) == 5){											
											   		BTagSF = bTool->getWeight(tempJet->Pt(), TempEta,5,"CSVM",SFsys);
												}else if(fabs(jet_flavorSF) == 4){
													BTagSF = bTool->getWeight(tempJet->Pt(), TempEta,4,"CSVM",SFsys);
												}else if( fabs(jet_flavorSF)< 4 || fabs(jet_flavorSF)==21)  {
													LightJetSF = bTool->getWeight(tempJet->Pt(), TempEta,1,"CSVM",SFsys);
												}
												
												if(BTagSF < 0){
												
													cout << "WARNING: negative SF" << "Jetpt: " << tempJet->Pt() << " Eta: " << TempEta << endl; 
													cout << "BtagSF: " << BTagSF << " LightJetSF: " << LightJetSF << endl; 
													
													
													BTagSF = 1;  // temporarly fix
												} else if (LightJetSF < 0){
												
													cout << "WARNING: negative SF" << "Jetpt: " << tempJet->Pt() << " Eta: " << TempEta << endl; 
													cout << "BtagSF: " << BTagSF << " LightJetSF: " << LightJetSF << endl;
												
													LightJetSF = 1;    // temporarly fix 
												}
											}
										}
										
																	
										//set a unique seed 
										double phi = jet_phiSF; 
										double sin_phi = sin(phi*1000000);
										double seed = fabs(static_cast<int>(sin_phi*100000));
										if(!isData){
											//Initialize class
											BTagSFUtil* btsfutil = new BTagSFUtil(seed);
										
											//cout << "BTAGGING SF ADAPTION:: Before modification: bTagged = " << bTagged << endl; 
									
											//modify tags 
										
											btsfutil->modifyBTagsWithSF(bTagged, jet_flavorSF, BTagSF, BTagEff, LightJetSF, LightJetEff);
										}	
										//cout << "BTAGGING SF ADAPTION:: After modification: bTagged = " << bTagged << endl; 
										
										if(bTagged){
											nJetsBT++;
			      								nTightJetsBT++;	
										}
									} // end proper jet demand
									else if (tempJet->btag_combinedSecondaryVertexBJetTags()> 0.679 && fabs(tempJet->Eta()) < 2.5){
										bTagged = true; 
											
										if (!isData ){
											if(fabs(jet_flavorSF) == 5){											
										   		BTagSF = bTool->getWeight(tempJet->Pt(), TempEta,5,SFsys);
											}else if(fabs(jet_flavorSF) == 4){
												BTagSF = bTool->getWeight(tempJet->Pt(), TempEta,4,SFsys);
											}else if( fabs(jet_flavorSF)< 4 || fabs(jet_flavorSF)==21) {
												LightJetSF = bTool->getWeight(tempJet->Pt(), TempEta,1,SFsys);
										
											}
											
											if(BTagSF < 0){
												
												cout << "WARNING: negative SF" << "Jetpt: " << tempJet->Pt() << " Eta: " << TempEta << endl; 
												cout << "BtagSF: " << BTagSF << " LightJetSF: " << LightJetSF << endl; 
													
													
												BTagSF = 1;  // temporarly fix
											} else if (LightJetSF < 0){
												
												cout << "WARNING: negative SF" << "Jetpt: " << tempJet->Pt() << " Eta: " << TempEta << endl; 
												cout << "BtagSF: " << BTagSF << " LightJetSF: " << LightJetSF << endl;
											
												LightJetSF = 1;    // temporarly fix 
											}
										}
									} // closing else if 

									//set a unique seed 
									double phi = jet_phiSF; 
									double sin_phi = sin(phi*1000000);
									double seed = fabs(static_cast<int>(sin_phi*100000));
									
									if(!isData){
										//Initialize class
										BTagSFUtil* btsfutil = new BTagSFUtil(seed);
									
										//cout << "BTAGGING SF ADAPTION:: Before modification: bTagged = " << bTagged << endl; 
									
										//modify tags 
										btsfutil->modifyBTagsWithSF(bTagged, jet_flavorSF, BTagSF, BTagEff, LightJetSF, LightJetEff);
									}
									//cout << "BTAGGING SF ADAPTION:: After modification: bTagged = " << bTagged << endl; 
										
									if(bTagged){
										nJetsBT++;
									}
	

								} // closing loop over jets 

								histo_btagged_jets->Fill(nJetsBT); 

							
							
							
						      		//Filling the Tree (at pre-selection level, leptons and mll)
		      						lum = lumi;
			
		      						xlWeight = weight;
			
		      						puweight =  lumiWeight;
		      						rawWeight = xlweight;
					
		      						npu = event->nPu(0);
		      						nvertex = vertex.size();
			
		      						metPt = mets[0]->Pt();
		      						metPx = mets[0]->Px();
		      						metPy = mets[0]->Py();
		      
		      						ptLepton = new std::vector<double>; 
		      						pxLepton = new std::vector<double>; 
		      						pyLepton = new std::vector<double>; 
		      						pzLepton = new std::vector<double>; 
		      						eLepton = new std::vector<double>; 
		      						qLepton = new std::vector<double>;
			
		      						ptJet = new std::vector<double>; 
		      						pxJet = new std::vector<double>; 
		     						pyJet = new std::vector<double>; 
		      						pzJet = new std::vector<double>; 
		      						eJet = new std::vector<double>; 
								etaJet = new std::vector<double>; 
		      						qJet = new std::vector<double>; 
								SFjet = new std::vector<double>;
		      						btJPBJet = new std::vector<double>; 
		     						btBJPBJet = new std::vector<double>; 
		      						btCSVBJet = new std::vector<double>;
		      						btCSVBmvaJet = new std::vector<double>; 
			
		      						ptLepton->push_back(lepton0.Pt());
		      						ptLepton->push_back(lepton1.Pt());
			
		      						pxLepton->push_back(lepton0.Px());
		      						pxLepton->push_back(lepton1.Px());
			
		      						pyLepton->push_back(lepton0.Py());
		      						pyLepton->push_back(lepton1.Py());
			
		      						pzLepton->push_back(lepton0.Pz());
		      						pzLepton->push_back(lepton1.Pz());
			
		      						eLepton->push_back(lepton0.Energy());
		      						eLepton->push_back(lepton1.Energy());
			
		      						qLepton->push_back(q0);
		      						qLepton->push_back(q1);
			
		      						for (unsigned int i =0; i < selectedJets.size(); i ++){
									TRootJet* tempJet = (TRootJet*) selectedJets[i];
									ptJet->push_back(tempJet->Pt());
									pxJet->push_back(tempJet->Px());
									pyJet->push_back(tempJet->Py());
									pzJet->push_back(tempJet->Pz());
									eJet->push_back(tempJet->Energy());
									etaJet->push_back(tempJet->Eta());
									qJet->push_back(tempJet->charge());
									double bTag_SF_b ;
									if ( fabs(tempJet->Eta()) < 2.4){
										if (isData || !scaleFactor){
											bTag_SF_b= 1;	
										}else{									
									   		bTag_SF_b = bTool->getWeight(tempJet->Pt(), tempJet->Eta(),5,"CSVM",0);
										}
									} else { 
										bTag_SF_b = 1;
									}
									SFjet->push_back(bTag_SF_b);
									btJPBJet->push_back(tempJet->btag_jetProbabilityBJetTags() );
									btBJPBJet->push_back(tempJet->btag_jetBProbabilityBJetTags());
									btCSVBJet->push_back(tempJet->btag_combinedSecondaryVertexBJetTags() );
									btCSVBmvaJet->push_back(tempJet->btag_combinedSecondaryVertexMVABJetTags());  
		      						}
				
	
		      						myTree->Fill();
			
		      						delete ptLepton;
		      						delete pxLepton;
		      						delete pyLepton;
		      						delete pzLepton;
		      						delete eLepton;
		      						delete qLepton;
			
		     				 		delete ptJet;
		      						delete pxJet;
		      						delete pyJet;
		      						delete pzJet;
		      						delete eJet;
								delete etaJet;
		      						delete qJet;
								delete SFjet;
		      						delete btJPBJet;
		      						delete btBJPBJet;
		      						delete btCSVBJet;
		      						delete btCSVBmvaJet;
								
								
								

								// --> Invariant mass cut for ee and mumu such that these are outside the Z mass window 
		      						if (pair.M() > 101 || pair.M() < 81 || mode == 0){
									// Fill the histos in bin 6 with events after the Z mass window cut
									cutflow->Fill(6, weight);
									cutflow_raw->Fill(6);
									
									// --> MET cut in ee and mumu such that events without genuine MET are reduced
									if (met_pt > 30 || mode == 0){
			 							cutflow->Fill(7, weight);
			  							cutflow_raw->Fill(7);
			  							OutPut3 << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;
										
										// Filling all the regions, so just regions with jets
			  							if (nJets !=0){
			    								TRootJet* jet = (TRootJet*) selectedJets[iJet];
			    								double Ht = lepton0.Pt() + lepton1.Pt() + jet->Pt() + met_pt; 
											
											// --> Ht cut for the emu mode in order to remove additional drell yann background
			    								if (Ht > 160 || mode != 0){
			      									if (nJets == 1 && nTightJetsBT == 1 && nJetsBT == 1 && bTagged)Regions->Fill(1, weight);
			      									if (nJets == 1 && nTightJetsBT == 2)  Regions->Fill(2, weight);
			      									if (nJets == 1 && nTightJetsBT > 0)  Regions->Fill(3, weight);
			      									if (nJets == 1 && nTightJetsBT > 1)  Regions->Fill(4, weight);
			      									if (nJets == 2 && nTightJetsBT == 0)  Regions->Fill(5, weight);
			      									if (nJets == 2 && nTightJetsBT == 1)  Regions->Fill(6, weight);
			      									if (nJets == 2 && nTightJetsBT == 2)  Regions->Fill(7, weight);
			      									if (nJets == 2 && nTightJetsBT > 0)  Regions->Fill(8, weight);
			      									if (nJets == 2 && nTightJetsBT > 1)  Regions->Fill(9, weight);
			      									if (nJets > 1 && nTightJetsBT == 0)  Regions->Fill(10, weight);
			      									if (nJets > 1 && nTightJetsBT == 1)  Regions->Fill(11, weight);
			      									if (nJets > 1 && nTightJetsBT == 2)  Regions->Fill(12, weight);
			      									if (nJets > 1 && nTightJetsBT !=0 )  Regions->Fill(13, weight);
			      									if (nJets > 1 && nTightJetsBT > 1 )  Regions->Fill(14, weight);
			      									if (nJets == 3 && nTightJetsBT ==0 )  Regions->Fill(15, weight);
			      									if (nJets == 3 && nTightJetsBT ==1 )  Regions->Fill(16, weight);
			      									if (nJets == 3 && nTightJetsBT ==2 )  Regions->Fill(17, weight);
			      									if (nJets == 3 && nTightJetsBT ==3 )  Regions->Fill(18, weight);
			    								}// closing the loop for Ht cut in emu
			 							} // closing the loop for filling all regions
			    
			  							// --> (cut) Filling the signal region, regions where there is only 1 jet
			  							if(nJets == 1){
			   								OutPut4 << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;

			    								TRootJet* jet = (TRootJet*) selectedJets[iJet];
											
											// Fill the histos after claim of exactly 1 jet
			    								cutflow->Fill(8, weight);
			    								cutflow_raw->Fill(8);
											
											// --> (cut) exactly one btagged jet
			    								if (nJets == 1 && nTightJetsBT == 1 && nJetsBT == 1 && bTagged){
			      									OutPut5 << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;

			      									cutflow->Fill(9,weight);
			      									cutflow_raw->Fill(9);
												
												// --> Ht cut for emu mode 
			      									double Ht = lepton0.Pt() + lepton1.Pt() + jet->Pt() + met_pt; 
			      									if (Ht > 160 || mode != 0){
													cutflow->Fill(10, weight);
													cutflow_raw->Fill(10);	
												} // closing loop for Ht cut in signal region
											} // closing  loop for exactly 1 btagged jet 
										} // Closing loop for exatly one jet
									} // closing the loop for met cut in ee and mumu
								} // closing the  loop for the mll cut in ee and mumu 
							} // closing the loop over the events with mll > 20 GeV
						} // closing the loop over the events with no extra loose leptons				
					} // closing the loop over the events with two opposite charge leptons
				} //closing the loop over the events with a tight leptons selection	
			} // Closing the loop over the events with good reconstructed primary vertices 
		} //Closing the loop over the triggered events
        } // Closing the loop over the events in one dataset
	
	delete jecUnc;
      	delete jetTools;

        ////////////////////////////////
        ///      INFORMATION EVENTS  ///  ==> WHY THE USE OF SCALER 1 AND 2
	///          IN 1 DATASET    ///   
        ////////////////////////////////
      	double scaler1 = cutflow->GetBinContent(2) ;
	
      	if(scaler1<=0.0)		scaler1=1.;
	
      	cout << "--------------------------------------------------" << endl;
      	cout << "[Results Normalized:] " <<  endl;
      	cout << "All:       " <<  cutflow->GetBinContent(2) << " +/- "  << cutflow->GetBinError(2) << "\t = " << 100.*cutflow->GetBinContent(2)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(2)/scaler1 << "%" << endl;
      	cout << "HLT:       " <<  cutflow->GetBinContent(3) << " +/- "  << cutflow->GetBinError(3) <<  "\t = " << 100.*cutflow->GetBinContent(3)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(3)/scaler1 << "%" <<endl;
      	//cout << "PV:        " <<  cutflow->GetBinContent(4) << " +/- "  << cutflow->GetBinError(4) <<  "\t = " << 100.*cutflow->GetBinContent(4)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(4)/scaler1 << "%" <<endl;
      	cout << "Lep. Sel:  " <<  cutflow->GetBinContent(5) << " +/- "  << cutflow->GetBinError(5) <<  "\t = " << 100.*cutflow->GetBinContent(5)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(5)/scaler1 << "%" <<endl;
      	cout << "Lep. Veto: " <<  cutflow->GetBinContent(6) << " +/- "  << cutflow->GetBinError(6) <<  "\t = " << 100.*cutflow->GetBinContent(6)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(6)/scaler1 << "%" <<endl;
      	cout << "mll:       " <<  cutflow->GetBinContent(7) << " +/- "  << cutflow->GetBinError(7) <<  "\t = " << 100.*cutflow->GetBinContent(7)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(7)/scaler1 << "%" << endl;
      	cout << "MET:       " <<  cutflow->GetBinContent(8) << " +/- "  << cutflow->GetBinError(8) <<  "\t = " << 100.*cutflow->GetBinContent(8)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(8)/scaler1 << "%" << endl;
      	cout << "1 jet:     " <<  cutflow->GetBinContent(9) << " +/- "  << cutflow->GetBinError(9) <<  "\t = " << 100.*cutflow->GetBinContent(9)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(9)/scaler1 << "%" <<endl;
      	cout << "1 jet BT:  " <<  cutflow->GetBinContent(10) << " +/- " << cutflow->GetBinError(10) <<  "\t = " << 100.*cutflow->GetBinContent(10)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(10)/scaler1 << "%" <<endl;
      	cout << "Ht tW:     " <<  cutflow->GetBinContent(11) << " +/- " << cutflow->GetBinError(11) <<  "\t = " << 100.*cutflow->GetBinContent(11)/scaler1 << " +/- "  << 100.*cutflow->GetBinError(11)/scaler1 << "%" <<endl;
    
      
      	cout << "--------------------------------------------------" << endl;
      	cout << "[Results Raw:] " <<  endl;
	
      	double scaler2 =  cutflow_raw->GetBinContent(2);
	
      	if(scaler2 <=0.0) scaler2=1.;
	
      	cout << "All:       " <<  cutflow_raw->GetBinContent(2) << " +/- "  << cutflow_raw->GetBinError(2) <<  "\t = " << 100.*cutflow_raw->GetBinContent(2)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(2)/scaler2 << "%" << endl;
      	cout << "HLT:       " <<  cutflow_raw->GetBinContent(3) << " +/- "  << cutflow_raw->GetBinError(3) <<  "\t = " << 100.*cutflow_raw->GetBinContent(3)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(3)/scaler2 << "%" << endl;
      	// cout << "PV:        " <<  cutflow_raw->GetBinContent(4) << " +/- "  << cutflow_raw->GetBinError(4) <<  "\t = " << 100.*cutflow_raw->GetBinContent(4)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(4)/scaler2 << "%" <<endl;
      	cout << "Lep. Sel:  " <<  cutflow_raw->GetBinContent(5) << " +/- "  << cutflow_raw->GetBinError(5) <<  "\t = " << 100.*cutflow_raw->GetBinContent(5)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(5)/scaler2 << "%" <<endl;
      	cout << "Lep. Veto: " <<  cutflow_raw->GetBinContent(6) << " +/- "  << cutflow_raw->GetBinError(6) <<  "\t = " << 100.*cutflow_raw->GetBinContent(6)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(6)/scaler2 << "%" <<endl;
      	cout << "mll:       " <<  cutflow_raw->GetBinContent(7) << " +/- "  << cutflow_raw->GetBinError(7) <<  "\t = " << 100.*cutflow_raw->GetBinContent(7)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(7)/scaler2 << "%" <<endl;
      	cout << "MET:       " <<  cutflow_raw->GetBinContent(8) << " +/- "  << cutflow_raw->GetBinError(8) <<  "\t = " << 100.*cutflow_raw->GetBinContent(8)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(8)/scaler2 << "%" <<endl;
      	cout << "1 jet:     " <<  cutflow_raw->GetBinContent(9) << " +/- "  << cutflow_raw->GetBinError(9) <<  "\t = " << 100.*cutflow_raw->GetBinContent(9)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(9)/scaler2 << "%" <<endl;
      	cout << "1 jet BT:  " <<  cutflow_raw->GetBinContent(10) << " +/- " << cutflow_raw->GetBinError(10) <<  "\t = " << 100.*cutflow_raw->GetBinContent(10)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(10)/scaler2 << "%" <<endl;
      	cout << "Ht:        " <<  cutflow_raw->GetBinContent(11) << " +/- " << cutflow_raw->GetBinError(11) <<  "\t = " << 100.*cutflow_raw->GetBinContent(11)/scaler2 << " +/- "  << 100.*cutflow_raw->GetBinError(11)/scaler2 << "%" <<endl;
      
      	cout << "--------------------------------------------------" << endl;
      	cout << "[Jet Multiplicity Check:]" << endl;
      	cout << "1 jet 1 tag: " << Regions->GetBinContent(2) << " +/- " << Regions->GetBinError(2) << "\t = " << 100.*Regions->GetBinContent(2)/scaler1 << " +/- "  << 100.*Regions->GetBinError(2)/scaler1 << "%" << endl;
      	cout << "2 jet 1 tag: " << Regions->GetBinContent(7) << " +/- " << Regions->GetBinError(2) << "\t = " << 100.*Regions->GetBinContent(7)/scaler1 << " +/- "  << 100.*Regions->GetBinError(7)/scaler1 << "%" << endl;
      	cout << "2 jet 2 tag: " << Regions->GetBinContent(8) << " +/- " << Regions->GetBinError(2) << "\t = " << 100.*Regions->GetBinContent(8)/scaler1 << " +/- "  << 100.*Regions->GetBinError(8)/scaler1 << "%" <<endl;
     	cout << "--------------------------------------------------" << endl;
/*      
        ////////////////////////////
	//   info btag eff       //
	///////////////////////////

     	if (dataSetName == "tt"){ 	
	        double  eff_tt = 0.; 
		double  mean_eff_tt = 0.;  
		for(int i =0; i <  btag_eff_tt.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_tt[i] << endl; 
			eff_tt = eff_tt + btag_eff_tt[i];
		}
		mean_eff_tt = eff_tt/btag_eff_tt.size();
		cout << "Mean Btag Eff for tt is: " << mean_eff_tt << endl; 
	}
 	else if (dataSetName == "twdr"){    
		double  eff_twdr = 0.; 
		double  mean_eff_twdr = 0.;  
		for(int i =0; i <  btag_eff_twdr.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_twdr[i] << endl; 
			eff_twdr = eff_twdr + btag_eff_twdr[i];
		}
		mean_eff_twdr = eff_twdr/btag_eff_twdr.size();
		cout << "Mean Btag Eff for twdr is: " << mean_eff_twdr << endl; 
	}	
 	else if (dataSetName == "atwdr"){   
		double  eff_atwdr = 0.; 
		double  mean_eff_atwdr = 0.;  
		for(int i =0; i <  btag_eff_atwdr.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_atwdr[i] << endl; 
			eff_atwdr = eff_atwdr + btag_eff_atwdr[i];
		}
		mean_eff_atwdr = eff_atwdr/btag_eff_atwdr.size();
		cout << "Mean Btag Eff for atwdr is: " << mean_eff_atwdr << endl; 
	}	
 	else if (dataSetName == "t"){  	    
	        double  eff_t = 0.; 
		double  mean_eff_t = 0.;  
		for(int i =0; i <  btag_eff_t.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_t[i] << endl; 
			eff_t = eff_t + btag_eff_t[i];
		}
		mean_eff_t = eff_t/btag_eff_t.size();
		cout << "Mean Btag Eff for t is: " << mean_eff_t << endl; 
	}
 	else if (dataSetName == "at"){ 	   
		double  eff_at = 0.; 
		double  mean_eff_at = 0.;  
		for(int i =0; i <  btag_eff_at.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_at[i] << endl; 
			eff_at = eff_at + btag_eff_at[i];
		}
		mean_eff_at = eff_at/btag_eff_at.size();
		cout << "Mean Btag Eff for at is: " << mean_eff_at << endl; 
	}
	else if (dataSetName == "s"){     
		double  eff_s = 0.; 
		double  mean_eff_s = 0.;  
		for(int i =0; i <  btag_eff_s.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_s[i] << endl; 
			eff_s = eff_s + btag_eff_s[i];
		}
		mean_eff_s = eff_s/btag_eff_s.size();
		cout << "Mean Btag Eff for s is: " << mean_eff_s << endl; 
	}	
 	else if (dataSetName == "as"){   
		double  eff_as = 0.; 
		double  mean_eff_as = 0.;  
		for(int i =0; i <  btag_eff_as.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_as[i] << endl; 
			eff_as = eff_as + btag_eff_as[i];
		}
		mean_eff_as = eff_as/btag_eff_as.size();
		cout << "Mean Btag Eff for as is: " << mean_eff_as << endl; 
	}
	else if (dataSetName == "ww"){    
		double  eff_ww = 0.; 
		double  mean_eff_ww = 0.;  
		for(int i =0; i <  btag_eff_ww.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_ww[i] << endl; 
			eff_ww = eff_ww + btag_eff_ww[i];
		}
		mean_eff_ww = eff_ww/btag_eff_ww.size();
		cout << "Mean Btag Eff for ww is: " << mean_eff_ww << endl; 
	}	
 	else if (dataSetName == "wz"){ 	   
		double  eff_wz = 0.; 
		double  mean_eff_wz = 0.;  
		for(int i =0; i <  btag_eff_wz.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_wz[i] << endl; 
			eff_wz = eff_wz + btag_eff_wz[i];
		}
		mean_eff_wz = eff_wz/btag_eff_wz.size();
		cout << "Mean Btag Eff for wz is: " << mean_eff_wz << endl; 
	}
	else if (dataSetName == "zz"){  
		double  eff_zz = 0.; 
		double  mean_eff_zz = 0.;  
		for(int i =0; i <  btag_eff_zz.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_zz[i] << endl; 
			eff_zz = eff_zz + btag_eff_zz[i];
		}
		mean_eff_zz = eff_zz/btag_eff_zz.size();
		cout << "Mean Btag Eff for zz is: " << mean_eff_zz << endl; 
	}	
 	else if (dataSetName == "zjets"){    
	        double  eff_zjets = 0.; 
		double  mean_eff_zjets = 0.;  
		for(int i =0; i <  btag_eff_zjets.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_zjets[i] << endl; 
			eff_zjets = eff_zjets + btag_eff_zjets[i];
		}
		mean_eff_zjets = eff_zjets/btag_eff_zjets.size();
		cout << "Mean Btag Eff for zjets is: " << mean_eff_zjets << endl; 
	}	
	else if (dataSetName == "zjets_lowmll"){    
		double  eff_zjets_lowmll = 0.; 
		double  mean_eff_zjets_lowmll = 0.;  
		for(int i =0; i <  btag_eff_zjets_lowmll.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_zjets_lowmll[i] << endl; 
			eff_zjets_lowmll = eff_zjets_lowmll + btag_eff_zjets_lowmll[i];
		}
		mean_eff_zjets_lowmll = eff_zjets_lowmll/btag_eff_zjets_lowmll.size();
		cout << "Mean Btag Eff for zjets_lowmll is: " << mean_eff_zjets_lowmll << endl; 
	}
 	else if (dataSetName == "wjets"){  
		double  eff_wjets = 0.; 
		double  mean_eff_wjets = 0.;  
		for(int i =0; i <  btag_eff_wjets.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  btag_eff_wjets[i] << endl; 
			eff_wjets = eff_wjets + btag_eff_wjets[i];
		}
		mean_eff_wjets = eff_wjets/btag_eff_wjets.size();
		cout << "Mean Btag Eff for wjets is: " << mean_eff_wjets << endl; 
	}    
	
	cout << "-----------------------------------------------------------------------------------------------------------" << endl; 
     
        ////////////////////////////
	//   info btag eff       //
	///////////////////////////

     	if (dataSetName == "tt"){ 	
	        double  eff_tt = 0.; 
		double  mean_eff_tt = 0.;  
		for(int i =0; i <  fake_eff_tt.size(); i++){
			//cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_tt[i] << endl; 
			eff_tt = eff_tt + fake_eff_tt[i];
		}
		mean_eff_tt = eff_tt/fake_eff_tt.size();
		cout << "Mean fake Eff for tt is: " << mean_eff_tt << endl; 
	}
 	else if (dataSetName == "twdr"){    
		double  eff_twdr = 0.; 
		double  mean_eff_twdr = 0.;  
		for(int i =0; i <  fake_eff_twdr.size(); i++){
		//	cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_twdr[i] << endl; 
			eff_twdr = eff_twdr + fake_eff_twdr[i];
		}
		mean_eff_twdr = eff_twdr/fake_eff_twdr.size();
		cout << "Mean fake Eff for twdr is: " << mean_eff_twdr << endl; 
	}	
 	else if (dataSetName == "atwdr"){   
		double  eff_atwdr = 0.; 
		double  mean_eff_atwdr = 0.;  
		for(int i =0; i <  fake_eff_atwdr.size(); i++){
			//cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_atwdr[i] << endl; 
			eff_atwdr = eff_atwdr + fake_eff_atwdr[i];
		}
		mean_eff_atwdr = eff_atwdr/fake_eff_atwdr.size();
		cout << "Mean fake Eff for atwdr is: " << mean_eff_atwdr << endl; 
	}	
 	else if (dataSetName == "t"){  	    
	        double  eff_t = 0.; 
		double  mean_eff_t = 0.;  
		for(int i =0; i <  fake_eff_t.size(); i++){
			//cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_t[i] << endl; 
			eff_t = eff_t + fake_eff_t[i];
		}
		mean_eff_t = eff_t/fake_eff_t.size();
		cout << "Mean fake Eff for t is: " << mean_eff_t << endl; 
	}
 	else if (dataSetName == "at"){ 	   
		double  eff_at = 0.; 
		double  mean_eff_at = 0.;  
		for(int i =0; i <  fake_eff_at.size(); i++){
			//cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_at[i] << endl; 
			eff_at = eff_at + fake_eff_at[i];
		}
		mean_eff_at = eff_at/fake_eff_at.size();
		cout << "Mean fake Eff for at is: " << mean_eff_at << endl; 
	}
	else if (dataSetName == "s"){     
		double  eff_s = 0.; 
		double  mean_eff_s = 0.;  
		for(int i =0; i <  fake_eff_s.size(); i++){
			//cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_s[i] << endl; 
			eff_s = eff_s + fake_eff_s[i];
		}
		mean_eff_s = eff_s/fake_eff_s.size();
		cout << "Mean fake Eff for s is: " << mean_eff_s << endl; 
	}	
 	else if (dataSetName == "as"){   
		double  eff_as = 0.; 
		double  mean_eff_as = 0.;  
		for(int i =0; i <  fake_eff_as.size(); i++){
			//cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_as[i] << endl; 
			eff_as = eff_as + fake_eff_as[i];
		}
		mean_eff_as = eff_as/fake_eff_as.size();
		cout << "Mean fake Eff for as is: " << mean_eff_as << endl; 
	}
	else if (dataSetName == "ww"){    
		double  eff_ww = 0.; 
		double  mean_eff_ww = 0.;  
		for(int i =0; i <  fake_eff_ww.size(); i++){
		//	cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_ww[i] << endl; 
			eff_ww = eff_ww + fake_eff_ww[i];
		}
		mean_eff_ww = eff_ww/fake_eff_ww.size();
		cout << "Mean fake Eff for ww is: " << mean_eff_ww << endl; 
	}	
 	else if (dataSetName == "wz"){ 	   
		double  eff_wz = 0.; 
		double  mean_eff_wz = 0.;  
		for(int i =0; i <  fake_eff_wz.size(); i++){
		//	cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_wz[i] << endl; 
			eff_wz = eff_wz + fake_eff_wz[i];
		}
		mean_eff_wz = eff_wz/fake_eff_wz.size();
		cout << "Mean fake Eff for wz is: " << mean_eff_wz << endl; 
	}
	else if (dataSetName == "zz"){  
		double  eff_zz = 0.; 
		double  mean_eff_zz = 0.;  
		for(int i =0; i <  fake_eff_zz.size(); i++){
		//	cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_zz[i] << endl; 
			eff_zz = eff_zz + fake_eff_zz[i];
		}
		mean_eff_zz = eff_zz/fake_eff_zz.size();
		cout << "Mean fake Eff for zz is: " << mean_eff_zz << endl; 
	}	
 	else if (dataSetName == "zjets"){    
	        double  eff_zjets = 0.; 
		double  mean_eff_zjets = 0.;  
		for(int i =0; i <  fake_eff_zjets.size(); i++){
		//	cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_zjets[i] << endl; 
			eff_zjets = eff_zjets + fake_eff_zjets[i];
		}
		mean_eff_zjets = eff_zjets/fake_eff_zjets.size();
		cout << "Mean fake Eff for zjets is: " << mean_eff_zjets << endl; 
	}	
	else if (dataSetName == "zjets_lowmll"){    
		double  eff_zjets_lowmll = 0.; 
		double  mean_eff_zjets_lowmll = 0.;  
		for(int i =0; i <  fake_eff_zjets_lowmll.size(); i++){
		//	cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_zjets_lowmll[i] << endl; 
			eff_zjets_lowmll = eff_zjets_lowmll + fake_eff_zjets_lowmll[i];
		}
		mean_eff_zjets_lowmll = eff_zjets_lowmll/fake_eff_zjets_lowmll.size();
		cout << "Mean fake Eff for zjets_lowmll is: " << mean_eff_zjets_lowmll << endl; 
	}
 	else if (dataSetName == "wjets"){  
		double  eff_wjets = 0.; 
		double  mean_eff_wjets = 0.;  
		for(int i =0; i <  fake_eff_wjets.size(); i++){
		//	cout << "fakeefficiency for " << i << "th event: " <<  fake_eff_wjets[i] << endl; 
			eff_wjets = eff_wjets + fake_eff_wjets[i];
		}
		mean_eff_wjets = eff_wjets/fake_eff_wjets.size();
		cout << "Mean fake Eff for wjets is: " << mean_eff_wjets << endl; 
	}    
	
	cout << "-----------------------------------------------------------------------------------------------------------" << endl; 
*/    





       /////////////////////////////////////////////
       /// Writing the efficiency plots in a file //
       /////////////////////////////////////////////
        char eff_file[100];
        sprintf(eff_file,"BtagFiles/btag_eff_%d_%s.root", mode, name);
        
	TFile *fileEf = new TFile(eff_file,"RECREATE"); 
	
 	char titlePlotEff[100];
	
	sprintf(titlePlotEff,"eff_Bjet_%d_%s",mode,name);	
	TH1F* histo_eff_Bjet = (TH1F*) histo_ptBTagB->Clone("histo_eff_Bjet"); 
	histo_eff_Bjet ->SetNameTitle(titlePlotEff, " Efficiency of btagging for jets from b quark");
	histo_eff_Bjet->Divide(histo_ptBjet);
	
	sprintf(titlePlotEff,"eff_Cjet_%d_%s",mode,name);	
	TH1F* histo_eff_Cjet = (TH1F*) histo_ptBTagC->Clone("histo_eff_Cjet"); 
	histo_eff_Cjet ->SetNameTitle(titlePlotEff, " Efficiency of Btagging for jets from C quark");
	histo_eff_Cjet->Divide(histo_ptCjet);
	
	sprintf(titlePlotEff,"eff_Ljet_%d_%s",mode,name);	
	TH1F* histo_eff_Ljet = (TH1F*) histo_ptBTagL->Clone("histo_eff_Ljet"); 
	histo_eff_Ljet ->SetNameTitle(titlePlotEff, " Efficiency of Btagging for jets from L quark");
	histo_eff_Ljet->Divide(histo_ptLjet);
	
	histo_btagged_jets->SetDirectory(fileEf); 
	
	
	// Drawing them
	
	TCanvas *c1_Eff = new TCanvas();
	histo_eff_Bjet->SetMaximum(1);
	//histo_eff_Bjet->Rebin(2);
        histo_eff_Bjet->Draw();
        histo_eff_Bjet->GetYaxis()->SetTitle("Eff");
        histo_eff_Bjet->GetXaxis()->SetTitle("pt");
	
	histo_eff_Cjet->SetMaximum(1);
	//histo_eff_Cjet->Rebin(2);
	histo_eff_Cjet ->SetLineColor(kBlue);
	histo_eff_Cjet->Draw("sames");
    
      char plotNameEff[100]; 
      sprintf(plotNameEff,"eff_B_C_jet_%d_%s",mode,name);
      TString string = plotNameEff;   
      c1_Eff->SaveAs("BtagFiles/plots/" + string  + ".png");


	TCanvas *c2_Eff = new TCanvas();
	histo_eff_Ljet ->SetLineColor(kGreen);
	histo_eff_Ljet->SetMaximum(1);
	//histo_eff_Ljet->Rebin(2);
      	histo_eff_Ljet->Draw();
      	histo_eff_Ljet->GetYaxis()->SetTitle("Eff");
      	histo_eff_Ljet->GetXaxis()->SetTitle("pt");

      char plotNameMisEff[100]; 
      sprintf(plotNameMisEff,"fake_eff_Bjet_%d_%s",mode,name);
      TString stringMis = plotNameMisEff;  
      c2_Eff->SaveAs("BtagFiles/plots/" + stringMis  + ".png");
      
      
     

      TCanvas *c3_Eff = new TCanvas();
      histo_btagged_jets->Draw();
      histo_btagged_jets->GetYaxis()->SetTitle("#evt");
      histo_btagged_jets->GetXaxis()->SetTitle("#btagged jets");

      char plotNameCheckEff[100]; 
      sprintf(plotNameCheckEff,"Number_Btagged_jets_%d_%s",mode,name);
      TString stringCheck = plotNameCheckEff;  
      c3_Eff->SaveAs("BtagFiles/plots/" + stringCheck + ".png");
      
      
      

	
      fileEf->Write(); 
      fileEf->Close(); 
	
	

	
    }  // closing the loop over the datasets 
    
    
    

    
    // To display how long it took for the program to run
    cout << "*******************************************************************************************" << endl; 
    cout << "It took you " << ((double)clock() - start) /CLOCKS_PER_SEC << " to run the program" << endl; 
  
}

