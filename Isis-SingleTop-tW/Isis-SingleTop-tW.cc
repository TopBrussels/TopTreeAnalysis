// isis.marina.van.parijs@cern.ch
//  based on SingleTop-tW.cc of rebeca  
//
// Summer12 round
// 8 Tev

// This is a program that runs over the toptrees

#include "TStyle.h"
#include "TH3F.h"
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
    bool pdf = false;
    
    bool JESPlus = false; 
    bool JESMinus = false; 
    
    bool JERPlus = false; 
    bool JERMinus = false; 
    
    bool SFplus = false; 
    bool SFminus = false; 
    bool SFplus_c = false; 
    bool SFminus_c = false; 
    bool SFplus_l = false; 
    bool SFminus_l = false; 
    
    bool unclusteredUp = false; 
    bool unclusteredDown = false; 
    
    bool PUsysUp = false; 
    bool PUsysDown = false; 
    
    bool eleSFsysUp = false; 
    bool eleSFsysDown = false; 
    
    bool topmass_plus = false; 
    bool topmass_minus = false; 
    
    bool Q2_plus = false; 
    bool Q2_minus = false; 
    
    bool matching_plus = false; 
    bool matching_minus = false; 
    
    bool ZSFplus = false; 
    bool ZSFminus = false; 
    
    
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
    double leptonEff_Trig_ID = 1; // default value
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
                cout << "--SFplus: SF up +10% syst for b quarks" << endl;
                cout << "--SFminus: SF down -10% syst for b quarks" << endl;
		cout << "--SFplus_c: SF up +10% systfor c quarks" << endl;
                cout << "--SFminus_c: SF down -10% syst for c quarks" << endl;
		cout << "--SFplus_l: SF up +10% syst for light quarks" << endl;
                cout << "--SFminus_l: SF down -10% syst for light quarks" << endl;
                cout << "--PUup: PU reweghting scaled up " << endl;
                cout << "--PUdown: PU reweghting scaled down " << endl;
                cout << "--uncMETup: Unclustered MET syst. Up " << endl;
                cout << "--uncMETdown: Unclustered MET syst. Down " << endl;
                cout << "--NoPU: Do not apply pileup re-weighting" << endl;
                cout << "--NoSF: Do not apply b-tag scale factor" << endl;
                cout << "--RAW: Do not apply pileup re-weighting or b-tag scale factor" << endl;
                cout << "--3D: 3D Pileup reweighting" << endl;
		cout << "--xml myxml.xml Xml file" << endl; 
		cout << "--eleSFplus   lepton ID/trigger eff scaled up" << endl; 
		cout << "--eleSFminus:  lepton ID/trigger eff scaled down" << endl; 
		cout << "--ZSFplus: Z/gamma reweighing up, twice the scale factors" << endl; 
		cout << "--ZSFminus: Z/gamma reweighing down, no scale factor" << endl; 
                return 0;
        }
	if (argval == "--ZSFplus"){
		ZSFplus = true; 
	}
	if (argval == "--ZSFminus"){
		ZSFminus = true; 
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
        if (argval=="--eleSFplus") {
                eleSFsysUp = true;
        }
        if (argval=="--eleSFminus"){
                eleSFsysDown = true;
        }	
	if (argval=="--SFplus_c") {
                SFplus_c = true;
        }
        if (argval=="--SFminus_c"){
                SFminus_c = true;
        }
	if (argval=="--SFplus_l") {
                SFplus_l = true;
        }
        if (argval=="--SFminus_l"){
                SFminus_l = true;
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
    if(mode == 2){
        cout << " Di-Electron channel " << endl;
    }
    if(mode == 1){
        cout << " Di-Muon channel " << endl;
    }
    cout << "******************************************************************************" << endl; 
    
    
    /////////////////////////////////
    /// LOAD BTAG EFF  rootfiles  ///
    /////////////////////////////////
    char Load_Eff_File_B[100];
    sprintf(Load_Eff_File_B,"BtagFiles/rootfiles/EFFbtag_0_Bjets.root");
    
    char Load_Eff_File_C[100];
    sprintf(Load_Eff_File_C,"BtagFiles/rootfiles/EFFbtag_0_Cjets.root");
    
    char Load_Eff_File_L[100];
    sprintf(Load_Eff_File_L,"BtagFiles/rootfiles/EFFbtag_0_Ljets.root");
    
    TFile *Efile_B = new TFile(Load_Eff_File_B, "read"); 
    TFile *Efile_C = new TFile(Load_Eff_File_C, "read"); 
    TFile *Efile_L = new TFile(Load_Eff_File_L, "read");  
    

    char TotitlePlot[100];
    char GetPlot[100];
    
    sprintf(TotitlePlot,"Eff_pt_btagB_%d_tt",mode);
    sprintf(GetPlot,"histo_tt_ptEffB");
    TH1F* histo_Eff_pt_btagB_tt = (TH1F*) Efile_B->Get(GetPlot); 
    if(histo_Eff_pt_btagB_tt == 0){
    	cout << " ERROR: nullpointer for btag eff plots Bjets" <<endl; 
	//return; 
    
    }
    histo_Eff_pt_btagB_tt->SetName(TotitlePlot); 
    
    sprintf(TotitlePlot,"Eff_pt_btagB_%d_twdr",mode);
    sprintf(GetPlot,"histo_twdr_ptEffB");
    TH1F* histo_Eff_pt_btagB_twdr = (TH1F*) Efile_B->Get(GetPlot); 
    histo_Eff_pt_btagB_twdr->SetName(TotitlePlot); 
    
     sprintf(TotitlePlot,"Eff_pt_btagC_%d_tt",mode);
    sprintf(GetPlot,"histo_tt_ptEffC");
    TH1F* histo_Eff_pt_btagC_tt = (TH1F*) Efile_C->Get(GetPlot); 
    histo_Eff_pt_btagC_tt->SetName(TotitlePlot); 
    
    sprintf(TotitlePlot,"Eff_pt_btagC_%d_twdr",mode);
    sprintf(GetPlot,"histo_twdr_ptEffC");
    TH1F* histo_Eff_pt_btagC_twdr = (TH1F*) Efile_C->Get(GetPlot); 
    histo_Eff_pt_btagC_twdr->SetName(TotitlePlot); 
    
    sprintf(TotitlePlot,"Eff_pt_btagL_%d_tt",mode);
    sprintf(GetPlot,"histo_tt_ptEffL");
    TH1F* histo_Eff_pt_btagL_tt = (TH1F*) Efile_L->Get(GetPlot); 
    histo_Eff_pt_btagL_tt->SetName(TotitlePlot); 
    
    sprintf(TotitlePlot,"Eff_pt_btagL_%d_twdr",mode);
    sprintf(GetPlot,"histo_twdr_ptEffL");
    TH1F* histo_Eff_pt_btagL_twdr = (TH1F*) Efile_L->Get(GetPlot); 
    histo_Eff_pt_btagL_twdr->SetName(TotitlePlot); 
    
    
    
        
    /////////////////////////////////
    ///       CONFIGURATION       ///  leptoneff should be checked!!! (12/3/2013)
    /////////////////////////////////
    
    
      // Luminosity and xml files
  // Only runs A, A recover, B, C (24, v2)
     if      (mode == 0){  // emu	 
     	lumi = 11966.617;  
	xmlfile ="config/twemu.xml";    

	leptonEff_Trig_ID = 0.915;
	if(eleSFsysUp){
		leptonEff_Trig_ID += 0.017;
	}else if(eleSFsysDown){
		leptonEff_Trig_ID = leptonEff_Trig_ID - 0.017;
	}
      }
     else if (mode == 1){ 
     	 lumi = 12067.294;  	
	 xmlfile = "config/twmumu.xml";
	 leptonEff_Trig_ID = 0.963;
	 if(eleSFsysUp){
		leptonEff_Trig_ID += 0.022;
	}else if(eleSFsysDown){
		leptonEff_Trig_ID = leptonEff_Trig_ID - 0.022;
	}
     }
     else if (mode == 2){	 
     	lumi = 12093.792;  	
     	xmlfile = "config/twee.xml";
	leptonEff_Trig_ID = 0.903;
	if(eleSFsysUp){
		leptonEff_Trig_ID += 0.021;
	}else if(eleSFsysDown){
		leptonEff_Trig_ID = leptonEff_Trig_ID - 0.021;
	}
     }
 
    

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
	bool isZjets = false; 
	bool isQ2up = false;
	bool isQ2down = false; 
	bool isTopmassup = false; 
	bool isTopmassdown = false; 
	bool isMatchingup = false; 
	bool isMatchingdown = false; 
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
	else if (dataSetName == "twds"){        sprintf(name, "tw_ds"); 	xlweight =lumi*11.1*0.1/2970009;             isSingleTop = true;} 
        else if (dataSetName == "atwds"){       sprintf(name, "atw_ds");     	xlweight = lumi*11.1*0.1/2940588;            isSingleTop = true;} 
        else if (dataSetName == "t"){           sprintf(name, "t");             xlweight = lumi*56.4/3748832;             } 
        else if (dataSetName == "at"){          sprintf(name, "at");            xlweight = lumi*30.7/180719;           }
	else if (dataSetName == "s"){           sprintf(name, "s");             xlweight = lumi*3.79/259960;             } 
        else if (dataSetName == "as"){          sprintf(name, "as");            xlweight = lumi*1.76/139974;           } 
        else if (dataSetName == "ww"){          sprintf(name, "ww");            xlweight = lumi*54.838/10000413;          } 
        else if (dataSetName == "wz"){          sprintf(name, "wz");            xlweight = lumi*22.44/9900267;          } 
        else if (dataSetName == "zz"){          sprintf(name, "zz");            xlweight = lumi*9.03/9799891 ;           } 
        else if (dataSetName == "zjets"){       sprintf(name, "zjets");         xlweight = lumi*3532.8/30364599;    isZjets =true;    } 
        else if (dataSetName == "zjets_lowmll"){sprintf(name, "zjets_lowmll");  xlweight = lumi*860.5/7059426;       isZjets = true;   } 
        else if (dataSetName == "wjets"){       sprintf(name, "wjets");         xlweight = lumi*36257.2/57411352;         }  
	
	// systematics
	else if (dataSetName == "tt_Q2_up"){          sprintf(name, "tt");       xlweight = lumi*225.197/4703202;        isTop = true;       isQ2up = true;} 
        else if (dataSetName == "twdr_Q2_up"){        sprintf(name, "tw_dr");        xlweight =lumi*11.1*0.1/1493129;             isSingleTop = true;  isQ2up = true;} //10% BR, onlydilepton
        else if (dataSetName == "atwdr_Q2_up"){       sprintf(name, "atw_dr");  	     xlweight = lumi*11.1*0.1/1492532;            isSingleTop = true;   isQ2up = true;} 
	else if (dataSetName == "tt_Q2_down"){          sprintf(name, "tt");           xlweight = lumi*225.197/5346767;       isTop = true;       isQ2down = true;} 
        else if (dataSetName == "twdr_Q2_down"){        sprintf(name,	"tw_dr");         xlweight =lumi*11.1*0.1/14923129;             isSingleTop = true;  isQ2down = true;} 
        else if (dataSetName == "atwdr_Q2_down"){       sprintf(name,	"atw_dr");        xlweight = lumi*11.1*0.1/1493099;            isSingleTop = true;   isQ2down = true;}
	
	else if (dataSetName == "tt_Topmass_up"){          sprintf(name,"tt");            xlweight = lumi*225.197/4733472;        isTop = true;       isTopmassup = true;} 
        else if (dataSetName == "twdr_Topmass_up"){ sprintf(name,"tw_dr");                xlweight =lumi*11.1*0.1/1493427;             isSingleTop = true;  isTopmassup = true;} 
        else if (dataSetName == "atwdr_Topmass_up"){       sprintf(name,	"atw_dr");        xlweight = lumi*11.1*0.1/1493387;            isSingleTop = true;   isTopmassup = true;} 
	else if (dataSetName == "tt_Topmass_down"){          sprintf(name,"tt");            xlweight = lumi*225.197/4358130;        isTop = true;       isTopmassdown = true;} 
        else if (dataSetName == "twdr_Topmass_down"){        sprintf(name,"tw_dr");         xlweight =lumi*11.1*0.1/1489878;             isSingleTop = true;  isTopmassdown = true;} 
        else if (dataSetName == "atwdr_Topmass_down"){       sprintf(name, "atw_dr");        xlweight = lumi*11.1*0.1/1478196;            isSingleTop = true;   isTopmassdown = true;}
	
	else if (dataSetName == "tt_Matching_up"){          sprintf(name,"tt");            xlweight = lumi*225.197/5415003;        isTop = true;       isMatchingup = true;} 
	else if (dataSetName == "tt_Matching_down"){          sprintf(name,"tt");            xlweight = lumi*225.197/5456715;        isTop = true;       isMatchingdown = true;} 
        
	
	
	

	
	
	
	
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
	else if (!isData && isTopmassup){   sprintf(rootFileName,"outputs/TopMassUp_%d_%s.root", mode, name);}
	else if (!isData && isTopmassdown){  sprintf(rootFileName,"outputs/TopMassDown_%d_%s.root", mode, name);}
	else if (!isData && isQ2up)	 {   sprintf(rootFileName,"outputs/Q2Up_%d_%s.root", mode, name);}
	else if (!isData && isQ2down){       sprintf(rootFileName,"outputs/Q2Down_%d_%s.root", mode, name);}
	else if (!isData && isMatchingup){  sprintf(rootFileName,"outputs/matchingUp_%d_%s.root", mode, name);}
	else if (!isData && isMatchingdown){ sprintf(rootFileName,"outputs/matchingDown_%d_%s.root", mode, name);}
        else if (!isData && eleSFsysUp){     sprintf(rootFileName,"outputs/eleSFsysUp_%d_%s.root", mode, name);}
        else if (!isData && eleSFsysDown){   sprintf(rootFileName,"outputs/eleSFsysDown_%d_%s.root", mode, name);}
        else if (!isData && unclusteredUp){  sprintf(rootFileName,"outputs/METsysUp_%d_%s.root", mode, name);}
        else if (!isData && unclusteredDown){sprintf(rootFileName,"outputs/METsysDown_%d_%s.root", mode, name);}
        else if (!isData && PUsysUp){        sprintf(rootFileName,"outputs/PUsysUp_%d_%s.root", mode, name);}
        else if (!isData && PUsysDown){      sprintf(rootFileName,"outputs/PUsysDown_%d_%s.root", mode, name);}
        else if (!isData && !reweightPU){    sprintf(rootFileName,"outputs/out_noPU_%d_%s.root", mode, name);}
	else if (!isData && ZSFplus){	     sprintf(rootFileName,"outputs/out_ZSFsysUp_%d_%s.root", mode, name);}
	else if (!isData && ZSFminus){	     sprintf(rootFileName,"outputs/out_ZSFsysDown_%d_%s.root", mode, name);}
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
        
	
	
	char myFile[300];
        sprintf(myFile,"pdf_unc/forPDF/pdf_signal_%d_%s.txt", mode, name);
        ofstream salida(myFile); 
	
	sprintf(myFile,"pdf_unc/forPDF/pdf_2j1t_%d_%s.txt", mode, name);
        ofstream salida2j1t(myFile);
	
	sprintf(myFile,"pdf_unc/forPDF/pdf_2j2t_%d_%s.txt", mode, name);
        ofstream salida2j2t(myFile);
        
        // Define the objects
        // ==> vector <kind> V
        //     This defines a vector V with elements of sort kind
                
        vector <TRootVertex*>   vertex; 
        vector <TRootMuon*>     initial_muons; 
        vector <TRootElectron*> initial_electrons; 
        vector <TRootJet*>      initial_jets; 
        vector <TRootJet*>      initial_jets_corrected; 
	vector <TRootJet*>      initial_jets_corrected_beforeJES_plus; 
	vector <TRootJet*>      initial_jets_corrected_beforeJES_minus;
	vector <TRootJet*>      initial_jets_corrected_beforeJER_plus;
	vector <TRootJet*>      initial_jets_corrected_beforeJER_minus;
        vector <TRootMET*>      mets; 
        vector <TRootGenJet*>   genjets;
/*        
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
	//// non btag eff for c jets /////
	////////////////////////////////////
	std::vector<double> cfake_eff_tt; 
	std::vector<double> cfake_eff_twdr; 	 
 	std::vector<double> cfake_eff_atwdr;  
 	std::vector<double> cfake_eff_t; 
 	std::vector<double> cfake_eff_at;
	std::vector<double> cfake_eff_s;    
 	std::vector<double> cfake_eff_as;
	std::vector<double> cfake_eff_ww;  
 	std::vector<double> cfake_eff_wz;
	std::vector<double> cfake_eff_zz;  
 	std::vector<double> cfake_eff_zjets;  
	std::vector<double> cfake_eff_zjets_lowmll; 
 	std::vector<double> cfake_eff_wjets; 

	
	
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
*/
	
	
	
	// for checking btagging
	std::vector<double> btag_eff_check1; 
	std::vector<double> btag_eff_check2; 
	
	
	
	
	
	
	
        /////////////////////////////////
        ///    Pile Up reweighting    /// 
	
	    
        
         
        LumiReWeighting LumiWeights; 
	
       //bool PUsysUp = false; 
       //bool PUsysDown = false; 
    
      if(PUsysUp){
           LumiWeights = LumiReWeighting("pileupHistos/pileup_MC_Summer12.root","pileupHistos/pileup_tW_2012Data53X_UpToRun203002/sys_up.root", "pileup", "pileup"); 
      }else if(PUsysDown){
          LumiWeights = LumiReWeighting("pileupHistos/pileup_MC_Summer12.root","pileupHistos/pileup_tW_2012Data53X_UpToRun203002/sys_down.root", "pileup", "pileup"); 
      }else{
	  LumiWeights = LumiReWeighting("pileupHistos/pileup_MC_Summer12.root","pileupHistos/pileup_tW_2012Data53X_UpToRun203002/nominal.root", "pileup", "pileup"); 
	   }
	  

        
        
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
        TH1F* Regions = new TH1F ("Regions" ," The different regions", 40,0.5,40.5); // for example  according to number of jets and b-tagged jets
        
        TH1F* pileup_weights = new TH1F("pileup_weights", "The calculated PU weights", 1000,0,10); 
        TH1F* pileup_weights3D = new TH1F("pileup_weights3D", "The calculated 3DPU weights", 1000,0,10); 
	
	TH1F* nvertex_beforePU = new TH1F("nvertex_beforePU", "The #vertices before PU reweighting", 1000,0,1000); 
        TH1F* nvertex_afterPU = new TH1F("nvertex_afterPU", "The #vertices after PU reweighting", 1000,0,1000); 
        
        TH1F* met_check_nom = new TH1F("met_check_nom", "The checking of nominal uncMet", 1000,0,2); 
        TH1F* met_check_up = new TH1F("met_check_up", "The checking of up uncMet", 1000,0,2);
	TH1F* met_check_down = new TH1F("met_check_down", "The checking of down uncMet", 1000,0,2); 
	
	TH1F* ptsys_check = new TH1F("ptsys_check", "The checking of pt sys for uncMET", 1000,0,2);
	TH1F* ptsys_check_1j1t = new TH1F("ptsys_check_1j1t", "The checking of pt sys for uncMET", 1000,0,2);
        
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
	
	sprintf(titlePlot,"btagged_jets_%d_%s",mode,name);
  	TH1F* histo_btagged_jets = new TH1F(titlePlot, " Btagged all jets",  10,  0, 10 ); 
	
	sprintf(titlePlot,"btagged_tightjets_%d_%s",mode,name);
  	TH1F* histo_btagged_tightjets = new TH1F(titlePlot, " Btagged tight jets",  10,  0, 10 ); 
	
	sprintf(titlePlot,"jets_%d_%s",mode,name);
  	TH1F* histo_jets = new TH1F(titlePlot, " tight jets",  10,  0, 10 ); 
	
	sprintf(titlePlot,"3d_btagged_tightjets_%d_%s",mode,name);
  	TH3F* histo_3d_btagged_tightjets = new TH3F(titlePlot, " Btagged  jets",  10,  -0.5, 9.5,   10,  -0.5, 9.5,10,  -0.5, 9.5  ); 
	
	sprintf(titlePlot,"2d_btagged_tightjets_%d_%s",mode,name);
  	TH2F* histo_2d_btagged_tightjets = new TH2F(titlePlot, " Btagged  jets",  10,  -0.5, 9.5,10,  -0.5, 9.5  );
	
	sprintf(titlePlot,"2d_btagged_tightjets_noHt_%d_%s",mode,name);
  	TH2F* histo_2d_btagged_tightjets_noHt = new TH2F(titlePlot, " Btagged  jets - no Ht cut",  10,  -0.5, 9.5,10,  -0.5, 9.5  );
/*	
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
	
	sprintf(titlePlot,"etaBTagB_%d_%s",mode,name);
  	TH1F* histo_etaBTagB = new TH1F(titlePlot, " eta of btagged jets for bjets", 100,  -3, 3);
	
	sprintf(titlePlot,"etaBTagC_%d_%s",mode,name);
  	TH1F* histo_etaBTagC = new TH1F(titlePlot, " eta of btagged jets for cjets", 100,  -3, 3);
	
  	sprintf(titlePlot,"etaBTagL_%d_%s",mode,name);
  	TH1F* histo_etaBTagL = new TH1F(titlePlot, " eta of btagged jets for ljets", 100,  -3, 3);
  
  	sprintf(titlePlot,"etaBjet_%d_%s",mode,name);
  	TH1F* histo_etaBjet = new TH1F(titlePlot, " eta of jets coming from a b quark ", 100,  -3, 3);
  
 	sprintf(titlePlot,"etaCjet_%d_%s",mode,name);
 	TH1F* histo_etaCjet = new TH1F(titlePlot, " eta of jets coming from a c quark ", 100 , -3, 3);
 
  	sprintf(titlePlot,"etaLjet_%d_%s",mode,name);
  	TH1F* histo_etaLjet = new TH1F(titlePlot, " eta of jets coming from light quarks ", 100 , -3, 3); 
	

	
        

 	histo_ptBTagB ->Sumw2(); 
   	histo_ptBTagC ->Sumw2(); 
	histo_ptBTagL ->Sumw2(); 
	histo_ptBjet ->Sumw2(); 
	histo_ptCjet ->Sumw2(); 
	histo_ptLjet ->Sumw2(); 
	histo_etaBTagB ->Sumw2(); 
	histo_etaBTagC ->Sumw2(); 
	histo_etaBTagL ->Sumw2(); 
	histo_etaBjet ->Sumw2(); 
	histo_etaCjet ->Sumw2(); 
	histo_etaLjet ->Sumw2(); 
*/	histo_btagged_jets->Sumw2(); 
 	histo_btagged_tightjets->Sumw2();
 	histo_jets->Sumw2();
	histo_3d_btagged_tightjets->Sumw2();
  	histo_2d_btagged_tightjets->Sumw2();
  	histo_2d_btagged_tightjets_noHt->Sumw2();    
        
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
	std::vector<bool> *Btagjet;
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
	myTree->Branch("Btagjet", "std::vector<bool>",&Btagjet);
        myTree->Branch("btJPBJet","std::vector<double>",&btJPBJet);
        myTree->Branch("btBJPBJet","std::vector<double>",&btBJPBJet);
        myTree->Branch("btCSVBJet","std::vector<double>",&btCSVBJet);
        myTree->Branch("btCSVBmvaJet","std::vector<double>",&btCSVBmvaJet);
        
        
        ////////////////////////////////
        ///    INFORMATION FOR USER  ///    
        ////////////////////////////////
        cout << "[Info:] output rootfile named " << rootFileName << endl; 
        cout << "[Info:] mode = " << mode << ", lumi: " <<  lumi << " pb, sample: " << name << ", base weight: " << xlweight << " , xml file: " << xmlfile << endl;
      
        if (JERPlus ||JERMinus || JESPlus || JESMinus ||  SFplus || SFminus ||SFplus_c || SFminus_c ||SFplus_l || SFminus_l || unclusteredUp || unclusteredDown 
          || !reweightPU || !scaleFactor || PUsysUp || PUsysDown || Pu3D || topmass_plus ||topmass_minus || ZSFplus || ZSFminus) {
                cout << "[Warning:] Non-standard options, ignore if you did it conciously" << endl;
                
                if (JERPlus) cout << "[Warning:] JER systematics on, plus. Note that for btagging the nominal values are used" << endl;
                if (JERMinus) cout << "[Warning:] JER systematics on, minus. Note that for btagging the nominal values are used" << endl;
                if (JESPlus) cout << "[Warning:] JES systematics on, plus. Note that for btagging the nominal values are used" << endl;
                if (JESMinus) cout << "[Warning:] JES systematics on, minus. Note that for btagging the nominal values are used" << endl;
		if ( eleSFsysUp) cout <<"[Warning:] SF up with one sigma for leptons " << endl;
                if (eleSFsysDown) cout <<"[Warning:]  SF down with one sigma for leptons " << endl;
                if (SFplus) cout <<"[Warning:] SF up 10% for b quarks " << endl;
                if (SFminus) cout <<"[Warning:]  SF down 10% for b quarks " << endl;
		if (SFplus_c) cout <<"[Warning:] SF up 10% for c quarks" << endl;
                if (SFminus_c) cout <<"[Warning:]  SF down 10% for c quarks " << endl;
		if (SFplus_l) cout <<"[Warning:] SF up 10% for light quarks " << endl;
                if (SFminus_l) cout <<"[Warning:]  SF down 10% for light quarks" << endl;
                if (unclusteredUp) cout <<"[Warning:] unclustered MET up 10% " << endl;
                if (unclusteredDown) cout <<"[Warning:] unclustered MET down 10% " << endl;
                if (!reweightPU && !isData) cout << "[Warning:] You are NOT applying PU re-weighting " << endl;
                if (!scaleFactor && !isData) cout << "[Warning:] You are NOT applying the b-tagging SF " << endl;
                if (PUsysUp) cout <<"[Warning:] PU up " << endl;
                if (PUsysDown) cout <<"[Warning:] PU down " << endl;
        	if (topmass_plus) cout << " Warning: topmass is 178.5 GeV " << endl; 
		if (topmass_minus) cout << "Warning: topmass is 166.5 GeV" << endl; 
		if (ZSFplus) cout << "Warning: Z/gamma scalefactors doubled" << endl; 
		if (ZSFminus) cout << "Warning: no Z/gamma scalefactors applied" << endl; 
        } 
        else{
                cout << "[Info:] Standard setup " << endl;
		pdf = true; 
        }
        
        cout << "[Info:] " << datasets[d]->NofEvtsToRunOver() << " total events" << endl;
        if (runHLT) cout << "[Info:] You have the HLT activated, this might be slower than the usual. " << endl;
	
	
        if(SFplus && !SFplus_c){
		cout << "[WARNING: have you forgotten the c jets? ]" << endl; 
	}else if(!SFplus && SFplus_c){
		cout << "[WARNING: have you forgotten the b jets? ]" << endl;
	}
	
	

        ////////////////////////////////
        ///    LOOP OVER THE EVENTS  ///    
        ////////////////////////////////
        for(int ievent = 0; ievent < datasets[d]->NofEvtsToRunOver(); ievent++){
		// HT DEFINING 
		double Ht_temp = 0;
		double ptsysX_temp = 0; 
		double ptsysY_temp = 0; 
	
	        // For setting tbranch 
	        std::vector<double> btag_booleans;
		
		// for zjets sf 
		double ZjetsSF = 1.0;
		
		 
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
                        
			
			
                        //If PU reweighting is turned on
                        if(reweightPU){
			       
                                weight *= lumiWeight ;      
                        }
			
                } // closing loop for PU reweighting for MC
                
                ////////////////////////////////
                ///  JES (Jet Energy Scale)  ///     
                ///        SYSTEMATICS       ///      
                //////////////////////////////// 
                
                // If the JES systematics are on, minus
                if(JESPlus){
			initial_jets_corrected_beforeJES_plus = initial_jets_corrected; 
                        jetTools->correctJetJESUnc(initial_jets_corrected, "plus", 1); 
                }
                // If the JES systematics are on, plus
                else if(JESMinus){
		 	initial_jets_corrected_beforeJES_minus = initial_jets_corrected;
                        jetTools->correctJetJESUnc(initial_jets_corrected, "minus", 1); 
                }
                
                
                ////////////////////////////////
                ///    JER(JE Resolution)    ///     The jet energy resolution in the MC is better than the one from the real data.
                ///       SMEARING           ///     Smear the MC energy resolution in order to mimic the one from data.       
                ////////////////////////////////
                if(JERMinus){
		        initial_jets_corrected_beforeJER_minus = initial_jets_corrected;
                        jetTools->correctJetJER(initial_jets_corrected,genjets,mets[0], "minus", false); // false means don't use old numbers, but new ones
                }
                else if(JERPlus){
			initial_jets_corrected_beforeJER_plus = initial_jets_corrected;
                        jetTools->correctJetJER(initial_jets_corrected,genjets,mets[0], "plus", false); // false means don't use old numbers, but new ones
                }
                else{
			if(JESMinus) jetTools->correctJetJER(initial_jets_corrected_beforeJES_minus,genjets,mets[0], "nominal", false); // false means don't use old numbers, but new ones
			if(JESPlus) jetTools->correctJetJER(initial_jets_corrected_beforeJES_plus,genjets,mets[0], "nominal", false); // false means don't use old numbers, but new ones
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
                                        itrigger = treeLoader.iTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", currentRun);
                                        isecondtrigger = treeLoader.iTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", currentRun);
                                } 
                                else if (mode == 1){
                                        itrigger = treeLoader.iTrigger("HLT_Mu17_Mu8_v*", currentRun);
                                        isecondtrigger = treeLoader.iTrigger("HLT_Mu17_TkMu8_v*", currentRun);
                                } 
                                else if (mode == 2){
                                        itrigger = treeLoader.iTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", currentRun);
                                }
                        } // closing the HLT for data loop
                        
                        //For the MC, there is no triggerpath
                        else {
                                itrigger = true;     
                                isecondtrigger = true;
                        } // closing the HLT for MC
                        
                        //if one of the triggers is used, the trigging is set on true   ==> WHAT IS HAPPENING? THESE STATEMENTS MAKE NO SENSE, what is this trigged? 
			// only for trigged = true, are cuts made
                        if (itrigger || isecondtrigger){ 
                                trigged = true;
				if(!isData){weight *= leptonEff_Trig_ID;} // so the event weight is the PUeff * Trig/ID eff * weight from lumi 
                         } 
                         else{
                                trigged = false;
                        }       
                } // closing the HLT run loop
                else{               
	      		trigged = true; 
			if(!isData){weight *= leptonEff_Trig_ID;} 
		}// if HLT makes no difference --> it doesn't (is checked, the event count stays the same)
		
                ////////////////////////////////
                /// SELECTION & CUTFLOW      ///
                ////////////////////////////////
                Selection selection(initial_jets_corrected, initial_muons,initial_electrons,mets); 
		
		
			Selection selection_beforeJES_plus(initial_jets_corrected_beforeJES_plus, initial_muons,initial_electrons,mets);
		
			Selection selection_beforeJES_minus(initial_jets_corrected_beforeJES_minus, initial_muons,initial_electrons,mets);
		
			Selection selection_beforeJER_plus(initial_jets_corrected_beforeJER_plus, initial_muons,initial_electrons,mets);
		
			Selection selection_beforeJER_minus(initial_jets_corrected_beforeJER_minus, initial_muons,initial_electrons,mets);
		

                
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
        
        
                //is already applied at skim level so won't make any difference 
                isGoodPV = true; 
                
                //--------------------------------------------------------------
                // Set up of the unclustered MET systematic    
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
		
		double met_px_0 = mets[0]->Px();
	  	double met_py_0 = mets[0]->Py();
	    
	  	if(unclusteredUp){
	    		met_px += uncmet_px*0.1;
	    		met_py += uncmet_py*0.1;
	  	} 
		if(unclusteredDown){
	    		met_px -= uncmet_px*0.1;
	    		met_py -= uncmet_py*0.1;
	  	}
	    
	  	double met_pt = sqrt(met_px*met_px + met_py*met_py);
		double met_pt_0 = sqrt(met_px_0*met_px_0 + met_py_0*met_py_0);
		
		double met_diff = met_pt_0/met_pt; 
		
		if(unclusteredUp){
			met_check_up ->Fill(met_diff);
		} else if(unclusteredDown){
			met_check_down->Fill(met_diff);
		} else {
			met_check_nom->Fill(met_diff);
 		}
                
		//---------------------------------------------
		// Zjets sf only for MC zjets (from danny meeting single top 21/3)
		//------------------------------------------------
		if(isZjets && !ZSFminus && !ZSFplus){
			if(mode == 0){
				if(met_pt < 10){ ZjetsSF = 0.9028; }
				else if(met_pt > 10 && met_pt < 20){ ZjetsSF = 0.9497; }
				else if(met_pt >= 20 && met_pt < 30){ ZjetsSF = 1.0189; }
				else if(met_pt >= 30 && met_pt < 40){ ZjetsSF = 1.0988; }
				else if(met_pt >= 40 && met_pt < 50){ ZjetsSF = 1.17415; }
				else if(met_pt >= 50 && met_pt < 60){ ZjetsSF = 1.25145; }	
				else{ ZjetsSF = 1.26325; }
				
			}
			else if (mode == 1){
				if(met_pt < 10){ ZjetsSF = 0.8841; }
				else if(met_pt > 10 && met_pt < 20){ ZjetsSF = 0.9386; }
				else if(met_pt >= 20 && met_pt < 30){ ZjetsSF = 1.0131; }
				else if(met_pt >= 30 && met_pt < 40){ ZjetsSF = 1.1012; }
				else if(met_pt >= 40 && met_pt < 50){ ZjetsSF = 1.1850; }
				else if(met_pt >= 50 && met_pt < 60){ ZjetsSF = 1.2500; }	
				else{ ZjetsSF = 1.3071; }			
			}
			else{
				if(met_pt < 10){ ZjetsSF = 0.9215; }
				else if(met_pt > 10 && met_pt < 20){ ZjetsSF = 0.9608; }
				else if(met_pt >= 20 && met_pt < 30){ ZjetsSF = 1.0247; }
				else if(met_pt >= 30 && met_pt < 40){ ZjetsSF = 1.0964; }
				else if(met_pt >= 40 && met_pt < 50){ ZjetsSF = 1.1633; }
				else if(met_pt >= 50 && met_pt < 60){ ZjetsSF = 1.2529; }	
				else{ ZjetsSF = 1.2194; }			
			}
		
		
		
		
			weight *= ZjetsSF;
		
		}
		if(isZjets && !ZSFminus && ZSFplus){
			if(mode == 0){
				if(met_pt < 10){ ZjetsSF = 2*0.9028; }
				else if(met_pt > 10 && met_pt < 20){ ZjetsSF = 2*0.9497; }
				else if(met_pt >= 20 && met_pt < 30){ ZjetsSF = 2*1.0189; }
				else if(met_pt >= 30 && met_pt < 40){ ZjetsSF = 2*1.0988; }
				else if(met_pt >= 40 && met_pt < 50){ ZjetsSF = 2*1.17415; }
				else if(met_pt >= 50 && met_pt < 60){ ZjetsSF = 2*1.25145; }	
				else{ ZjetsSF = 2*1.26325; }
				
			}
			else if (mode == 1){
				if(met_pt < 10){ ZjetsSF = 2*0.8841; }
				else if(met_pt > 10 && met_pt < 20){ ZjetsSF = 2*0.9386; }
				else if(met_pt >= 20 && met_pt < 30){ ZjetsSF = 2*1.0131; }
				else if(met_pt >= 30 && met_pt < 40){ ZjetsSF = 2*1.1012; }
				else if(met_pt >= 40 && met_pt < 50){ ZjetsSF = 2*1.1850; }
				else if(met_pt >= 50 && met_pt < 60){ ZjetsSF = 2*1.2500; }	
				else{ ZjetsSF = 2*1.3071; }			
			}
			else{
				if(met_pt < 10){ ZjetsSF = 0.9215; }
				else if(met_pt > 10 && met_pt < 20){ ZjetsSF = 2*0.9608; }
				else if(met_pt >= 20 && met_pt < 30){ ZjetsSF = 2*1.0247; }
				else if(met_pt >= 30 && met_pt < 40){ ZjetsSF = 2*1.0964; }
				else if(met_pt >= 40 && met_pt < 50){ ZjetsSF = 2*1.1633; }
				else if(met_pt >= 50 && met_pt < 60){ ZjetsSF = 2*1.2529; }	
				else{ ZjetsSF = 2*1.2194; }			
			}
		
		
		
		
			weight *= ZjetsSF;
		
		}		
		
		
		
		
                //--------------------------------------------------------------
                // START OF CUTFLOW
                //---------------------------------------------------------------
		
		// Fill the histograms cutflow and cutflow_raw in the first bin with the number of events before any cut
                cutflow->Fill(1, weight);
	  	cutflow_raw->Fill(1);
		
		if(!trigged){
			cout << " not trigged" << endl; 
		}
		
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
				if(JESPlus) selection_beforeJES_plus.setJetCuts(20.,5.,0.01,1.,0.98,0.3,0.1);
				if(JERPlus) selection_beforeJER_plus.setJetCuts(20.,5.,0.01,1.,0.98,0.3,0.1);
				if(JESMinus) selection_beforeJES_minus.setJetCuts(20.,5.,0.01,1.,0.98,0.3,0.1);
				if(JERMinus) selection_beforeJER_minus.setJetCuts(20.,5.,0.01,1.,0.98,0.3,0.1);
             		 	selection.setDiMuonCuts(20.,2.4,0.20,999.);
              			selection.setDiElectronCuts(20.,2.5,0.15,0.04,0.5,1,0.3,1); //(20.,2.5,0.15,0.04,0.,1,0.3,1)
              			selection.setLooseMuonCuts(10.,2.5,0.2);
              			selection.setLooseDiElectronCuts(15.0,2.5,0.2,0.5); 
		
		
		
	      			//Select Objects and put them in a vector 
	      			vector<TRootJet*> selectedJets = selection.GetSelectedJets(true);
				vector<TRootJet*> selectedJets_beforeJES_plus = selection_beforeJES_plus.GetSelectedJets(true);
				vector<TRootJet*> selectedJets_beforeJER_plus = selection_beforeJER_plus.GetSelectedJets(true);
				vector<TRootJet*> selectedJets_beforeJES_minus = selection_beforeJES_minus.GetSelectedJets(true);
				vector<TRootJet*> selectedJets_beforeJER_minus = selection_beforeJER_minus.GetSelectedJets(true);
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
					        
/*						
						//////////////////////////////////////////////////
						// CALCULATION BTAG EFFICIENCY  - fakerate     ///
						//////////////////////////////////////////////////

						double btag_eff = 0.; 
						double btag_eff_c = 0.;
						double fake_eff = 0.; 
						
						if(!isData){					
							sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //Sort the selected jets based on Pt
							
							vector <int> mcParticlesTLV,selectedJetsTLV;
							TLorentzVector bQuark_vector;
							
							mcParticlesTLV.clear(); //make sure nothing is inside this vector
							selectedJetsTLV.clear();
							
      
							for(unsigned int iJet=0;iJet<selectedJets.size(); iJet++){
								TRootJet* tempJet = (TRootJet*) selectedJets[iJet];
								TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
								
								int pdgID = tempJet->partonFlavour();
								
								if(tempJet->Pt() > 30 && fabs(tempJet->Eta()) < 2.5 && TMath::Min(fabs(lepton0.DeltaR(tJet)),fabs(lepton1.DeltaR(tJet))) > 0.3){
	    								selectedJetsTLV.push_back(iJet);
									mcParticlesTLV.push_back(pdgID);
								} // closing proper jet statement
								
							} // closing for loop over selectedJets.size();
							

						//	for(unsigned int iJ=0; iJ<selectedJetsTLV.size(); iJ++) {
					//			int jet_Number = selectedJetsTLV[iJ];
					//			int pdgID = mcParticlesTLV[iJ];
								
						//		cout << "JETnr: " << jet_Number << " pdgID: " << pdgID << endl; 
					//		}
						
	
							vector< pair<unsigned int, unsigned int> > JetPartonPair;
							vector< pair<unsigned int, unsigned int> > JetPartonPair_c;
							vector< pair<unsigned int, unsigned int> > JetPartonPair_fake;
							
								
							for(unsigned int iJ=0; iJ<selectedJetsTLV.size(); iJ++) {
								int jet_Number =selectedJetsTLV[iJ];
								int pdgID = mcParticlesTLV[iJ];
								
	    							if(fabs(pdgID) == 5)
								{
	      								JetPartonPair.push_back( pair<unsigned int, unsigned int> (pdgID, jet_Number) );
										
	 							 } // END JETPARTONPAIR FILLING
								else if(pdgID == 4 )
								{
	      								JetPartonPair_c.push_back( pair<unsigned int, unsigned int> (pdgID, jet_Number) );
	 							 } // END JETPARTONPAIR FILLING
								else //if(pdgID == 21 || fabs(pdgID) < 4) // <5 light quarks, 21 gluon no fabs because electrically neutral
								{
	      								JetPartonPair_fake.push_back( pair<unsigned int, unsigned int> (pdgID, jet_Number) );
	 							 } // END JETPARTONPAIR FILLING
							} // end for loop 
								
								
							for(unsigned int iJ=0; iJ<JetPartonPair.size(); iJ++) {
								unsigned int jetnumber = JetPartonPair[iJ].second;
								//cout << "JETNR " << jetnumber << " JET " << endl; 
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								
								double jetPT = Jet_tempo->Pt(); 
								histo_ptBjet->Fill(jetPT);
								
								double jetETA = Jet_tempo->Eta(); 
								histo_etaBjet->Fill(jetETA);
							
							}
								
							for(unsigned int iJ=0; iJ<JetPartonPair_c.size(); iJ++) {
								unsigned int jetnumber = JetPartonPair_c[iJ].second;
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								
								double jetPT = Jet_tempo->Pt(); 
								histo_ptCjet->Fill(jetPT);
								
								double jetETA = Jet_tempo->Eta(); 
								histo_etaCjet->Fill(jetETA);
							
							}
							
							for(unsigned int iJ=0; iJ<JetPartonPair_fake.size(); iJ++) {
								unsigned int jetnumber = JetPartonPair_fake[iJ].second;
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								
								double jetPT = Jet_tempo->Pt(); 
								histo_ptLjet->Fill(jetPT);
								
								double jetETA = Jet_tempo->Eta(); 
								histo_etaLjet->Fill(jetETA);
							
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
								double jetETA = Jet_tempo->Eta();
							
								if (Jet_tempo->btag_combinedSecondaryVertexBJetTags()> 0.679){
			    						Number_btagged_jets++;
									histo_ptBTagB->Fill(jetPT);
									histo_etaBTagB->Fill(jetETA);
									
								} // closing loop counting btagged jets
							} // closing for loop over jetpartonpair.size()
							
							for(unsigned int iPair=0; iPair<JetPartonPair_c.size(); iPair++) 
							{
	    							unsigned int jetnumber = JetPartonPair_c[iPair].second;
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								double jetPT = Jet_tempo->Pt();
								double jetETA = Jet_tempo->Eta();
							
								if (Jet_tempo->btag_combinedSecondaryVertexBJetTags()> 0.679){
			    						Number_btagged_jets_c++;
									histo_ptBTagC->Fill(jetPT);
									histo_etaBTagC->Fill(jetETA);
								} // closing loop counting btagged jets
							} // closing for loop over jetpartonpair.size()
							
							for(unsigned int iPair=0; iPair<JetPartonPair_fake.size(); iPair++) 
							{
	    							unsigned int jetnumber = JetPartonPair_fake[iPair].second;
								TRootJet* Jet_tempo = (TRootJet*) selectedJets[jetnumber];
								double jetPT = Jet_tempo->Pt();
								double jetETA = Jet_tempo->Eta();
								if (Jet_tempo->btag_combinedSecondaryVertexBJetTags()> 0.679){
			    						Number_btagged_jets_fake++;
									histo_ptBTagL->Fill(jetPT);
									histo_etaBTagL->Fill(jetETA);
								} // closing loop counting btagged jets
							} // closing for loop over jetpartonpair.size()
							
								btag_eff = Number_btagged_jets/Number_matched_bquarks;
								btag_eff_c = Number_btagged_jets_c/Number_matched_cquarks;
								fake_eff = Number_btagged_jets_fake/Number_matched_lquarks;
							
							 					
							
							if(Number_matched_bquarks != 0){
							
		 						if (dataSetName == "tt"){ btag_eff_tt.push_back(btag_eff ); }
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
							
							if(Number_matched_cquarks != 0){
							
		 						if (dataSetName == "tt"){ 	
								   cfake_eff_tt.push_back(btag_eff_c ); }
 								else if (dataSetName == "twdr"){    cfake_eff_twdr.push_back(btag_eff_c ); }	
 								else if (dataSetName == "atwdr"){   cfake_eff_atwdr.push_back(btag_eff_c ); }	
 								else if (dataSetName == "t"){  	    cfake_eff_t.push_back(btag_eff_c ); }
 								else if (dataSetName == "at"){ 	   cfake_eff_at.push_back(btag_eff_c ); }
								else if (dataSetName == "s"){      cfake_eff_s.push_back(btag_eff_c ); }	
 								else if (dataSetName == "as"){    cfake_eff_as.push_back(btag_eff_c ); }
								else if (dataSetName == "ww"){    cfake_eff_ww.push_back(btag_eff_c ); }	
 								else if (dataSetName == "wz"){ 	    cfake_eff_wz.push_back(btag_eff_c ); }
								else if (dataSetName == "zz"){   cfake_eff_zz.push_back(btag_eff_c ); }	
 								else if (dataSetName == "zjets"){     cfake_eff_zjets.push_back(btag_eff_c ); }	
								else if (dataSetName == "zjets_lowmll"){    cfake_eff_zjets_lowmll.push_back(btag_eff_c ); }
 								else if (dataSetName == "wjets"){   cfake_eff_wjets.push_back(btag_eff_c ); }	
												 		
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
								// cout << "NON BTAGEFFICIENTIE CALCULATION: fake Btagged jets = " << Number_btagged_jets_fake << " Number lquarks = " <<  Number_matched_lquarks << " fake_eff = " << fake_eff << endl; 						 		}
							}
								
						} // closing if(!isData)
						
*/					
	
						
						
						
						
						
						
						
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
								int SFsys_b = 0;  
								if (SFminus) 	SFsys_b = -1;   // systematics down
		      						if (SFplus) 	SFsys_b = +1;    // systematics up
								
								int SFsys_c = 0;  
								if (SFminus_c) 	SFsys_c = -1;   // systematics down
		      						if (SFplus_c) 	SFsys_c = +1;    // systematics up
								
								int SFsys_l= 0;  
								if (SFminus_l) 	SFsys_l = -1;   // systematics down
		      						if (SFplus_l) 	SFsys_l = +1;    // systematics up
								
		      						int nJetsBT = 0;
		      						int nTightJetsBT = 0;
		     						int nJets = 0;
								int iJet = -5;

								/////////////////////////////////////////
								//  CALCULATION BTAGGING SF  adption /////// you have to correct for the fact that the btagging efficiency is not 100%
								/////////////////////////////////////////  and thus adapt your SF calculated with BtagWeight
								//double bTag_SF_b;
								
								double BTagSF; 
								double BTagEff; 
								double BTagSF_c; 
								double BTagEff_c; 
								double LightJetSF; 
								double LightJetEff; 
								
								double BTagEff_new;
								double BTagEff_c_new;
								double LightJetEff_new;
							
								int jet_flavorSF;
								double jet_phiSF; 
								double jet_etaSF; 
								bool bTagged = false; 	
								
								double check_btagged=0.; 
								double check_notbtagged=0.;
								double check_eff = 0.; 						
								
								btag_booleans.clear(); 
								
								for (unsigned int iJ =0; iJ < selectedJets.size(); iJ ++){
									bTagged = false;
									
									TRootJet* tempJet = (TRootJet*) selectedJets[iJ];
									if(JESPlus) tempJet = (TRootJet*) selectedJets_beforeJES_plus[iJ]; 
									if(JERPlus) tempJet = (TRootJet*) selectedJets_beforeJER_plus[iJ];
									if(JESMinus) tempJet = (TRootJet*) selectedJets_beforeJES_minus[iJ];
									if(JERMinus) tempJet = (TRootJet*) selectedJets_beforeJER_minus[iJ];
									TLorentzVector tJet(tempJet->Px(), tempJet->Py(), tempJet->Pz(), tempJet->Energy());
									
									if (tempJet->Pt() > 20 && fabs(tempJet->Eta()) < 2.4 ){	
										float TempEta = tempJet->Eta(); 
													
										float Temp_jetPt = tempJet->Pt();
										BTagEff = -1; 
										LightJetEff = -1; 
									
									
										if (tempJet->btag_combinedSecondaryVertexBJetTags()> 0.679){ 
											bTagged = true; 
										}
										
									
																	
										// Getting efficiencies out of histograms 
										if(!isData &&!isSingleTop){
											double xvalue = Temp_jetPt;
										 		
											int xbinB = histo_Eff_pt_btagB_tt->GetXaxis()->FindBin(xvalue);
											float btageff_B = histo_Eff_pt_btagB_tt->GetBinContent(xbinB); 							
											//cout << "Btagging eff (b): trial: jetpt = " << xvalue << " xbin = " << xbinB << " eff = " << btageff_B << endl;
											BTagEff = btageff_B;
											
											int xbinC = histo_Eff_pt_btagC_tt->GetXaxis()->FindBin(xvalue);
											float btageff_C = histo_Eff_pt_btagC_tt->GetBinContent(xbinC); 							
											//cout << "Btagging eff (c): trial: jetpt = " << xvalue << " xbin = " << xbinC << " eff = "<< btageff_C << endl;
											BTagEff_c = btageff_C; 
												
											int xbinL = histo_Eff_pt_btagL_tt->GetXaxis()->FindBin(xvalue);
											float btageff_L = histo_Eff_pt_btagL_tt->GetBinContent(xbinL); 							
											//cout << "Btagging eff (L): trial: jetpt = " << xvalue << " xbin = " << xbinL << " eff = "<<btageff_L << endl;
											LightJetEff = btageff_L;  
											
										}else if(isSingleTop && !isData){
											double xvalue = Temp_jetPt;
										 		
											int xbinB = histo_Eff_pt_btagB_twdr->GetXaxis()->FindBin(xvalue);
											float btageff_B = histo_Eff_pt_btagB_twdr->GetBinContent(xbinB); 							
											//cout << "Btagging eff (b): trial: jetpt = " << xvalue << " xbin = " << xbinB << " eff = " << btageff_B << endl;
											BTagEff = btageff_B;
												
											int xbinC = histo_Eff_pt_btagC_twdr->GetXaxis()->FindBin(xvalue);
											float btageff_C = histo_Eff_pt_btagC_twdr->GetBinContent(xbinC); 							
											//cout << "Btagging eff (c): trial: jetpt = " << xvalue << " xbin = " << xbinC << " eff = "<< btageff_C << endl;
											BTagEff_c = btageff_C; 
											
											int xbinL = histo_Eff_pt_btagL_twdr->GetXaxis()->FindBin(xvalue);
											float btageff_L = histo_Eff_pt_btagL_twdr->GetBinContent(xbinL); 							
											//cout << "Btagging eff (L): trial: jetpt = " << xvalue << " xbin = " << xbinL << " eff = "<<btageff_L << endl;
											LightJetEff = btageff_L;
										}
															
										
									

									

										if(!isData){	
											// Get scale factors from btag POG (pt and eta dependent)
											if(fabs(jet_flavorSF) == 5){											
												BTagSF = bTool->getWeight(tempJet->Pt(), TempEta,5,"CSVM",SFsys_b);
											}else if(fabs(jet_flavorSF) == 4){
												BTagSF_c = bTool->getWeight(tempJet->Pt(), TempEta,4,"CSVM",SFsys_c);
											}else //if( fabs(jet_flavorSF)< 4 || fabs(jet_flavorSF)==21) 
											{
												LightJetSF = bTool->getWeight(tempJet->Pt(), TempEta,1,"CSVM",SFsys_l);
											}
										
											// warnings	
											if(BTagSF < 0){
												cout << "WARNING: negative SF" << "Jetpt: " << tempJet->Pt() << " Eta: " << TempEta << endl; 
												cout << "BtagSF: " << BTagSF << " LightJetSF: " << LightJetSF << endl; 
												
												BTagSF = 1;  // temporarly fix
											} else if (LightJetSF < 0){
												
												cout << "WARNING: negative SF" << "Jetpt: " << tempJet->Pt() << " Eta: " << TempEta << endl; 
												cout << "BtagSF: " << BTagSF << " LightJetSF: " << LightJetSF << endl;
											
											LightJetSF = 1;    // temporarly fix 
											} else if (BTagSF_c < 0){
												
												cout << "WARNING: negative SF" << "Jetpt: " << tempJet->Pt() << " Eta: " << TempEta << endl; 
												cout << "BtagSF: " << BTagSF << " LightJetSF: " << LightJetSF << " c " <<BTagSF_c << endl;
												BTagSF_c = 1;    // temporarly fix 
											}
										}
											
										//set a unique seed 
										jet_flavorSF = tempJet->partonFlavour();  
										jet_phiSF = tempJet->Phi(); 
										
										
										double phi = jet_phiSF; 
										double sin_phi = sin(phi*1000000);
										int seed =(int) fabs(static_cast<int>(sin_phi*100000));
										if(!isData && scaleFactor){
											//Initialize class
											BTagSFUtil* btsfutil = new BTagSFUtil(seed);
										
											//cout << "BTAGGING SF ADAPTION:: Before modification: bTagged = " << bTagged << endl; 
									
											//modify tags 
										
											btsfutil->modifyBTagsWithSF(bTagged, jet_flavorSF, BTagSF, BTagEff, LightJetSF, LightJetEff,BTagSF_c, BTagEff_c);
										
											//cout << "BTAGGING SF ADAPTION:: After modification: bTagged = " << bTagged << endl; 
										} // closing MC adaption loop		
										
										
 
									
										
	
										if(!isData && jet_flavorSF == 5){
											if(bTagged){ 
												check_btagged++;
												//cout << "btag after is 0" << endl; 
											}
												check_notbtagged++;
												//cout << "btag after is 1" << endl;
										
											
										}

										if (tempJet->Pt() > 30  && TMath::Min(fabs(lepton0.DeltaR(tJet)),fabs(lepton1.DeltaR(tJet))) > 0.3) { // if proper jet 
				 							nJets++;
											iJet = iJ;
											
											Ht_temp += tempJet->Pt();
											ptsysX_temp += tempJet->Px();
											ptsysY_temp += tempJet->Py();
										
											if(bTagged){
												nTightJetsBT++;	
												nJetsBT++;
											}
											
										} // end proper jet statement
										else if(bTagged)
										{
											nJetsBT++;
										}
									} // pt > 20   eta < 2.4  such that btagging can be used
									else{
										bTagged = false; 
									
									} // jets with pt < 20 and eta > 2.4 aren't btagged
									// btag booleans for looper
									btag_booleans.push_back(bTagged);	
								} // closing loop over jets
								
								if(btag_booleans.size() != selectedJets.size()){
								
									cout << "WARNING tree for btagging is not correctly filled " << endl; 
								}
								
								 double check_eff_unc;
								
								if(!isData && check_notbtagged > 0){
									check_eff = check_btagged/check_notbtagged;
									double term1 = check_btagged/(check_notbtagged*check_notbtagged);
									double term2 = (check_btagged*check_btagged)/check_notbtagged;
									check_eff_unc = sqrt(term1+term2); 
									 
									//cout << "number of btagged MC: " << check_btagged << " - number of jets MC " << check_notbtagged << endl; 
									btag_eff_check1.push_back(check_eff);
									btag_eff_check2.push_back(check_eff_unc);
									//cout << "eff: " << check_eff << endl;  
								}
								        //cout << "number of btagged (all): " << nJetsBT << " tight ones: " << nTightJetsBT << endl; 
								
								
							     //   cout << "Tight btagged jets: " << nTightJetsBT << " Btagged ones: " << nJetsBT << endl; 
							
							
						      		//Filling the Tree (at pre-selection level, leptons and mll)
		      						lum = lumi;
			
		      						xlWeight = weight;
			
		      						puweight =  lumiWeight;
		      						rawWeight = xlweight;
					
		      						npu = event->nPu(0);
		      						nvertex = vertex.size();
			                                        
								
								
		      						metPt = met_pt; //mets[0]->Pt();
		      						metPx = met_px; // mets[0]->Px();
		      						metPy = met_py; //mets[0]->Py();
		      
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
								Btagjet = new std::vector<bool>;
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
									Btagjet->push_back(btag_booleans[i]);
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
								delete Btagjet;
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
									if (met_pt > 50 || mode == 0){
			 							cutflow->Fill(7, weight);
			  							cutflow_raw->Fill(7);
			  							OutPut3 << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;
										
										// Filling all the regions, so just regions with jets
			  							if (nJets !=0){
			    								TRootJet* jet = (TRootJet*) selectedJets[iJet];
			    								double Ht = lepton0.Pt() + lepton1.Pt() + Ht_temp + met_pt; 
											double ptsysX = lepton0.Px() + lepton1.Px() + ptsysX_temp + met_px;
											double ptsysY = lepton0.Py() + lepton1.Py() + ptsysY_temp + met_py;
											double ptsystem = sqrt(ptsysX *ptsysX + ptsysY *ptsysY);
											
											double ptsysX_0 = lepton0.Px() + lepton1.Px() + ptsysX_temp + met_px_0;
											double ptsysY_0 = lepton0.Py() + lepton1.Py() + ptsysY_temp + met_py_0;
											double ptsystem_0 = sqrt(ptsysX_0 *ptsysX_0 + ptsysY_0 *ptsysY_0);
											
											double ptdiff = ptsystem_0/ptsystem; 
											histo_2d_btagged_tightjets_noHt->Fill(nTightJetsBT,nJets, weight);
											
											
											ptsys_check->Fill(ptdiff);
											
											
											
											// --> Ht cut for the emu mode in order to remove additional drell yann background
			    								if (Ht > 160 || mode != 0){
											
												if (nJets == 1 && nTightJetsBT == 1 && nJetsBT == 1){
													ptsys_check_1j1t->Fill(ptdiff);
													xlWeight = weight;
													int id1 = event->idParton1();
													int id2 = event->idParton2();
													float x1 = event->xParton1();
													float x2 = event->xParton2();
													float q = event->factorizationScale();
													float ptsys = ptsystem;
													float ht = Ht;
				
													if (pdf) salida << weight << " " << x1 << " " << x2 << " " << q << " " << id1 << " " << id2 << " " << ptsys << " " << ht << " " << name << " " << mode << endl;
				
												}
												if (nJets == 2 && nTightJetsBT == 1 && nJetsBT == 1){
				
													xlWeight = weight;
													int id1 = event->idParton1();
													int id2 = event->idParton2();
													float x1 = event->xParton1();
													float x2 = event->xParton2();
													float q = event->factorizationScale();
													float ptsys = ptsystem;
													float ht = Ht;
				
													if (pdf) salida2j1t << weight << " " << x1 << " " << x2 << " " << q << " " << id1 << " " << id2 << " " << ptsys << " " << ht << " " << name << " " << mode << endl;
				
												}
												if (nJets == 2 && nTightJetsBT == 2 && nJetsBT == 2){
				
													xlWeight = weight;
													int id1 = event->idParton1();
													int id2 = event->idParton2();
													float x1 = event->xParton1();
													float x2 = event->xParton2();
													float q = event->factorizationScale();
													float ptsys = ptsystem;
													float ht = Ht;
				
													if (pdf) salida2j2t << weight << " " << x1 << " " << x2 << " " << q << " " << id1 << " " << id2 << " " << ptsys << " " << ht << " " << name << " " << mode << endl;
				
												}
												histo_jets->Fill(nJets,weight); 
								                        	histo_btagged_jets->Fill(nJetsBT,weight); 
								                        	histo_btagged_tightjets->Fill(nTightJetsBT,weight); 
												histo_3d_btagged_tightjets->Fill(nJetsBT,nTightJetsBT,nJets, weight);
												histo_2d_btagged_tightjets->Fill(nTightJetsBT,nJets, weight);
											
		      									        if (nJets == 1 && nTightJetsBT == 1 && nJetsBT == 1)Regions->Fill(1, weight);
			      									if (nJets == 1 && nJetsBT == 2)  Regions->Fill(2, weight);
			      									if (nJets == 1 && nTightJetsBT > 0)  Regions->Fill(3, weight);
			      									if (nJets == 1 && nTightJetsBT > 1)  Regions->Fill(4, weight);
			      									if (nJets == 2 && nTightJetsBT == 0)  Regions->Fill(5, weight);
			      									if (nJets == 2 && nTightJetsBT == 1 && nJetsBT ==1)  Regions->Fill(6, weight);
			      									if (nJets == 2 && nTightJetsBT == 2 && nJetsBT == 2)  Regions->Fill(7, weight);
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
											
											//cout << "EXACTLY ONE JET " << endl; 
											
											// --> (cut) exactly one btagged jet
			    								if (nJets == 1 && nTightJetsBT == 1 && nJetsBT == 1){
											        
											        //cout << "EXACTLY ONE BTAGGED JET " << endl; 
			      									OutPut5 << event->runId() << "\t" << event->lumiBlockId() << "\t" << event->eventId() << endl;

			      									cutflow->Fill(9,weight);
			      									cutflow_raw->Fill(9);
												
												// --> Ht cut for emu mode 
			      									double Ht = lepton0.Pt() + lepton1.Pt() + Ht_temp + met_pt; 
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
      	cout << "1 jet 1 tag: " << Regions->GetBinContent(1) << " +/- " << Regions->GetBinError(1) << "\t = " << 100.*Regions->GetBinContent(1)/scaler1 << " +/- "  << 100.*Regions->GetBinError(1)/scaler1 << "%" << endl;
      	cout << "2 jet 1 tag: " << Regions->GetBinContent(6) << " +/- " << Regions->GetBinError(6) << "\t = " << 100.*Regions->GetBinContent(6)/scaler1 << " +/- "  << 100.*Regions->GetBinError(6)/scaler1 << "%" << endl;
      	cout << "2 jet 2 tag: " << Regions->GetBinContent(7) << " +/- " << Regions->GetBinError(7) << "\t = " << 100.*Regions->GetBinContent(7)/scaler1 << " +/- "  << 100.*Regions->GetBinError(7)/scaler1 << "%" <<endl;
     	cout << "--------------------------------------------------" << endl;
	
	
	cout << "--------------------------------------------------" << endl;
	cout << " btag check (no logistics) " << endl; 
	cout << "--------------------------------------------------" << endl;
	
	cout << "VETOING LOOSE JETS: "<< endl;
	cout << "1jet 1tag: " << histo_3d_btagged_tightjets->GetBinContent(2,2,2) << "+/-" << histo_3d_btagged_tightjets->GetBinError(2,2,2) << endl; 
	cout << "2jet 1tag: " << histo_3d_btagged_tightjets->GetBinContent(2,2,3) << "+/-" << histo_3d_btagged_tightjets->GetBinError(2,2,3) << endl;
	cout << "2jet 2tag: " << histo_3d_btagged_tightjets->GetBinContent(3,3,3) << "+/-" << histo_3d_btagged_tightjets->GetBinError(3,3,3) << endl;
	cout << "" << endl;
	cout << " NO VETO ON LOOSE JETS: " << endl; 
	cout << "1jet 1tag: " << histo_2d_btagged_tightjets->GetBinContent(2,2) << "+/-" << histo_2d_btagged_tightjets->GetBinError(2,2) << endl; 
	cout << "2jet 1tag: " << histo_2d_btagged_tightjets->GetBinContent(2,3) << "+/-" << histo_2d_btagged_tightjets->GetBinError(2,3) << endl;
	cout << "2jet 2tag: " << histo_2d_btagged_tightjets->GetBinContent(3,3) << "+/-" << histo_2d_btagged_tightjets->GetBinError(3,3) << endl;
	cout << " " << endl; 
	cout << " NO VETO ON LOOSE JETS - NO HT CUT " << endl; 
	cout << "1jet 1tag: " << histo_2d_btagged_tightjets_noHt->GetBinContent(2,2) << "+/-" << histo_2d_btagged_tightjets_noHt->GetBinError(2,2) << endl; 
	cout << "2jet 1tag: " << histo_2d_btagged_tightjets_noHt->GetBinContent(2,3) << "+/-" << histo_2d_btagged_tightjets_noHt->GetBinError(2,3) << endl;
	cout << "2jet 2tag: " << histo_2d_btagged_tightjets_noHt->GetBinContent(3,3) << "+/-" << histo_2d_btagged_tightjets_noHt->GetBinError(3,3) << endl;
	
	cout << "--------------------------------------------------" << endl;
	cout << "BTAG CHECKING " << endl; 
	cout << "--------------------------------------------------" << endl;
	cout << "ALL (tight) JETS : " << endl; 
	cout << "Number of events with 0 jets: " << histo_jets->GetBinContent(1+0)  << " +/- " <<histo_jets->GetBinError(1+0) << endl; 
	cout << "Number of events with 1 jets: " << histo_jets->GetBinContent(1+1)  << " +/- " <<histo_jets->GetBinError(1+1) << endl; 
	cout << "Number of events with 2 jets: " << histo_jets->GetBinContent(1+2)  << " +/- " <<histo_jets->GetBinError(1+2) << endl; 
	cout << "Number of events with 3 jets: " << histo_jets->GetBinContent(1+3)  << " +/- " <<histo_jets->GetBinError(1+3) << endl; 
	cout << "Number of events with 4 jets: " << histo_jets->GetBinContent(1+4)  << " +/- " <<histo_jets->GetBinError(1+4) << endl; 
	cout << "Number of events with 5 jets: " << histo_jets->GetBinContent(1+5)  << " +/- " <<histo_jets->GetBinError(1+5) << endl; 
	cout << "--------------------------------------------------" << endl;
	cout << "ALL btagged JETS : " << endl; 
	cout << "Number of events with 0 bjets: " << histo_btagged_jets->GetBinContent(1+0)  << " +/- " <<histo_btagged_jets->GetBinError(1+0) << endl; 
	cout << "Number of events with 1 bjets: " << histo_btagged_jets->GetBinContent(1+1)  << " +/- " <<histo_btagged_jets->GetBinError(1+1) << endl; 
	cout << "Number of events with 2 bjets: " << histo_btagged_jets->GetBinContent(1+2)  << " +/- " <<histo_btagged_jets->GetBinError(1+2) << endl; 
	cout << "Number of events with 3 bjets: " << histo_btagged_jets->GetBinContent(1+3)  << " +/- " <<histo_btagged_jets->GetBinError(1+3) << endl; 
	cout << "Number of events with 4 bjets: " << histo_btagged_jets->GetBinContent(1+4)  << " +/- " <<histo_btagged_jets->GetBinError(1+4) << endl; 
	cout << "Number of events with 5 bjets: " << histo_btagged_jets->GetBinContent(1+5)  << " +/- " <<histo_btagged_jets->GetBinError(1+5) << endl; 
	cout << "--------------------------------------------------" << endl;
	cout << "TIGHT btagged JETS : " << endl; 
	cout << "Number of events with 0 bjets: " << histo_btagged_tightjets->GetBinContent(1+0)  << " +/- " <<histo_btagged_tightjets->GetBinError(1+0) << endl; 
	cout << "Number of events with 1 bjets: " << histo_btagged_tightjets->GetBinContent(1+1)  << " +/- " <<histo_btagged_tightjets->GetBinError(1+1) << endl; 
	cout << "Number of events with 2 bjets: " << histo_btagged_tightjets->GetBinContent(1+2)  << " +/- " <<histo_btagged_tightjets->GetBinError(1+2) << endl; 
	cout << "Number of events with 3 bjets: " << histo_btagged_tightjets->GetBinContent(1+3)  << " +/- " <<histo_btagged_tightjets->GetBinError(1+3) << endl; 
	cout << "Number of events with 4 bjets: " << histo_btagged_tightjets->GetBinContent(1+4)  << " +/- " <<histo_btagged_tightjets->GetBinError(1+4) << endl; 
	cout << "Number of events with 5 bjets: " << histo_btagged_tightjets->GetBinContent(1+5)  << " +/- "<< histo_btagged_tightjets->GetBinError(1+5) << endl;
	cout << "--------------------------------------------------" << endl;
        double ef = 0;
	double ef_unc = 0; 
        for(int i = 0; i< btag_eff_check1.size(); i++){
		ef = ef + btag_eff_check1[i];
		ef_unc = ef_unc + (btag_eff_check2[i] * btag_eff_check2[i]);
	}
	double eff = ef /btag_eff_check1.size(); 
	double eff_unc = sqrt(ef_unc)/btag_eff_check1.size();
	cout << " mean btag eff for bjets (MC): " << eff  << endl;
	cout << "--------------------------------------------------" << endl; 




/*      
        ////////////////////////////
	//   info btag eff: b       //
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
	//   info btag eff: c       //
	///////////////////////////

     	if (dataSetName == "tt"){ 	
	        double  eff_tt = 0.; 
		double  mean_eff_tt = 0.;  
		for(int i =0; i <  cfake_eff_tt.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_tt[i] << endl; 
			eff_tt = eff_tt + cfake_eff_tt[i];
		}
		mean_eff_tt = eff_tt/cfake_eff_tt.size();
		cout << "Mean Btag Eff c for tt is: " << mean_eff_tt << endl; 
	}
 	else if (dataSetName == "twdr"){    
		double  eff_twdr = 0.; 
		double  mean_eff_twdr = 0.;  
		for(int i =0; i <  cfake_eff_twdr.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_twdr[i] << endl; 
			eff_twdr = eff_twdr + cfake_eff_twdr[i];
		}
		mean_eff_twdr = eff_twdr/cfake_eff_twdr.size();
		cout << "Mean Btag Eff c for twdr is: " << mean_eff_twdr << endl; 
	}	
 	else if (dataSetName == "atwdr"){   
		double  eff_atwdr = 0.; 
		double  mean_eff_atwdr = 0.;  
		for(int i =0; i <  cfake_eff_atwdr.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_atwdr[i] << endl; 
			eff_atwdr = eff_atwdr + cfake_eff_atwdr[i];
		}
		mean_eff_atwdr = eff_atwdr/cfake_eff_atwdr.size();
		cout << "Mean Btag Eff c for atwdr is: " << mean_eff_atwdr << endl; 
	}	
 	else if (dataSetName == "t"){  	    
	        double  eff_t = 0.; 
		double  mean_eff_t = 0.;  
		for(int i =0; i <  cfake_eff_t.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_t[i] << endl; 
			eff_t = eff_t + cfake_eff_t[i];
		}
		mean_eff_t = eff_t/cfake_eff_t.size();
		cout << "Mean Btag Eff c for t is: " << mean_eff_t << endl; 
	}
 	else if (dataSetName == "at"){ 	   
		double  eff_at = 0.; 
		double  mean_eff_at = 0.;  
		for(int i =0; i <  cfake_eff_at.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_at[i] << endl; 
			eff_at = eff_at + cfake_eff_at[i];
		}
		mean_eff_at = eff_at/cfake_eff_at.size();
		cout << "Mean Btag Eff c for at is: " << mean_eff_at << endl; 
	}
	else if (dataSetName == "s"){     
		double  eff_s = 0.; 
		double  mean_eff_s = 0.;  
		for(int i =0; i <  cfake_eff_s.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_s[i] << endl; 
			eff_s = eff_s + cfake_eff_s[i];
		}
		mean_eff_s = eff_s/cfake_eff_s.size();
		cout << "Mean Btag Eff c for s is: " << mean_eff_s << endl; 
	}	
 	else if (dataSetName == "as"){   
		double  eff_as = 0.; 
		double  mean_eff_as = 0.;  
		for(int i =0; i <  cfake_eff_as.size(); i++){
			//cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_as[i] << endl; 
			eff_as = eff_as + cfake_eff_as[i];
		}
		mean_eff_as = eff_as/cfake_eff_as.size();
		cout << "Mean Btag Eff c for as is: " << mean_eff_as << endl; 
	}
	else if (dataSetName == "ww"){    
		double  eff_ww = 0.; 
		double  mean_eff_ww = 0.;  
		for(int i =0; i <  cfake_eff_ww.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_ww[i] << endl; 
			eff_ww = eff_ww + cfake_eff_ww[i];
		}
		mean_eff_ww = eff_ww/cfake_eff_ww.size();
		cout << "Mean Btag Eff c for ww is: " << mean_eff_ww << endl; 
	}	
 	else if (dataSetName == "wz"){ 	   
		double  eff_wz = 0.; 
		double  mean_eff_wz = 0.;  
		for(int i =0; i <  cfake_eff_wz.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_wz[i] << endl; 
			eff_wz = eff_wz + cfake_eff_wz[i];
		}
		mean_eff_wz = eff_wz/cfake_eff_wz.size();
		cout << "Mean Btag Eff c for wz is: " << mean_eff_wz << endl; 
	}
	else if (dataSetName == "zz"){  
		double  eff_zz = 0.; 
		double  mean_eff_zz = 0.;  
		for(int i =0; i <  cfake_eff_zz.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_zz[i] << endl; 
			eff_zz = eff_zz + cfake_eff_zz[i];
		}
		mean_eff_zz = eff_zz/cfake_eff_zz.size();
		cout << "Mean Btag Eff c for zz is: " << mean_eff_zz << endl; 
	}	
 	else if (dataSetName == "zjets"){    
	        double  eff_zjets = 0.; 
		double  mean_eff_zjets = 0.;  
		for(int i =0; i <  cfake_eff_zjets.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_zjets[i] << endl; 
			eff_zjets = eff_zjets + cfake_eff_zjets[i];
		}
		mean_eff_zjets = eff_zjets/cfake_eff_zjets.size();
		cout << "Mean Btag Eff c for zjets is: " << mean_eff_zjets << endl; 
	}	
	else if (dataSetName == "zjets_lowmll"){    
		double  eff_zjets_lowmll = 0.; 
		double  mean_eff_zjets_lowmll = 0.;  
		for(int i =0; i <  cfake_eff_zjets_lowmll.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_zjets_lowmll[i] << endl; 
			eff_zjets_lowmll = eff_zjets_lowmll + cfake_eff_zjets_lowmll[i];
		}
		mean_eff_zjets_lowmll = eff_zjets_lowmll/cfake_eff_zjets_lowmll.size();
		cout << "Mean Btag Eff c for zjets_lowmll is: " << mean_eff_zjets_lowmll << endl; 
	}
 	else if (dataSetName == "wjets"){  
		double  eff_wjets = 0.; 
		double  mean_eff_wjets = 0.;  
		for(int i =0; i <  cfake_eff_wjets.size(); i++){
		//	cout << "Btagefficiency for " << i << "th event: " <<  cfake_eff_wjets[i] << endl; 
			eff_wjets = eff_wjets + cfake_eff_wjets[i];
		}
		mean_eff_wjets = eff_wjets/cfake_eff_wjets.size();
		cout << "Mean Btag c Eff for wjets is: " << mean_eff_wjets << endl; 
	}    
	
	cout << "-----------------------------------------------------------------------------------------------------------" << endl;      
        ////////////////////////////
	//   info btag eff : L      //
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
    





      
       /////////////////////////////////////////////
       /// Writing the efficiency plots in a file //
       /////////////////////////////////////////////
        char eff_file[100];
        sprintf(eff_file,"BtagFiles/btag_eff_%d_%s.root", mode, name);
        
	TFile *fileEf = new TFile(eff_file,"RECREATE"); 
	
 	char titlePlotEff[100];
	
	
	histo_btagged_jets->SetDirectory(fileEf); 
	histo_ptBjet->SetDirectory(fileEf);
	histo_ptCjet->SetDirectory(fileEf);
	histo_ptLjet->SetDirectory(fileEf);
	histo_etaBjet->SetDirectory(fileEf); 
	histo_etaCjet->SetDirectory(fileEf);
	histo_etaLjet->SetDirectory(fileEf);
	histo_ptBTagB->SetDirectory(fileEf);
	histo_ptBTagC->SetDirectory(fileEf);
	histo_ptBTagL->SetDirectory(fileEf);
	histo_etaBTagB->SetDirectory(fileEf);
	histo_etaBTagC->SetDirectory(fileEf);
	histo_etaBTagL->SetDirectory(fileEf);
	


      fileEf->Write(); 
      fileEf->Close(); 
*/      
     
     fout->Write(); 
      fout->Close(); 
	

	
    }  // closing the loop over the datasets 
    
    
    

    
    // To display how long it took for the program to run
    cout << "*******************************************************************************************" << endl; 
    cout << "It took you " << ((double)clock() - start) /CLOCKS_PER_SEC << " to run the program" << endl; 
  
}

