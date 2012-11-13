// isis.marina.van.parijs@cern.ch
//
// Summer12 round
// 8 Tev

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
    
    
    // Configuration of the output format 
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
    ///    Analysis environment   ///     HOW DOES THIS WORK ??? 
    /////////////////////////////////    
    
    //Create an analysisenvironment
    AnalysisEnvironment anaEnv; 
    AnalysisEnvironmentLoader anaLoad(anaEnv, xmlfile.c_str()); 
    new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv); 
    
    
    
    // To display how long it took for the program to run
    cout << "It took you " << ((double)clock() - start) /CLOCKS_PER_SEC << " to run the program" << endl; 
  
}

