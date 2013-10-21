// isis.marina.van.parijs@cern.ch

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
    ///       SINGLE TOP TW ANALYSIS          ///
    /////////////////////////////////////////////
    
    //Make the plots nicer: colors, style,... 
    setMyStyle(); 
    
    // To see how long the program takes to run, I include a clock
    clock_t start = clock(); 
    
    cout << "******************************************************************************" << endl; 
    cout << " Welcome to the single top tW analysis: Summer12 53X at 8 TeV with MC " << endl;   
    cout << "******************************************************************************" << endl; 
    
    
    int mode = 0; 
    double lumi = 1000; //by default is the luminosity 1000
    double leptonEff_Trig_ID = 1; // default value
    string xmlfile = "config/twemu.xml" ;  // by default, use mode 0
    
   

    
    
    
    
    
        
    /////////////////////////////////
    ///       CONFIGURATION       ///  
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
        
	double xlweight;        // To define the reweighting of the MC compared to the data, if this is 1, then no reweighting is applied. 
                                // The reweighting is done as follows: 
                                //       xlweight = (cross-section x luminosity)/number of events in the toptree before the skimming
                                //This number of events is given in the mail that you receive with the urls of the skimmed toptrees (so you have to save these numbers)
        
        
        // Define the cross sections and weights for every data set
        // sprintf makes the string you call name contain data 
        char name[100];
	
        if (dataSetName == "data"){             sprintf(name, "data");          xlweight = 1;                           isData = true;}
        else if (dataSetName == "tt"){          sprintf(name, "tt");            xlweight = lumi*245.8/6830443;        isTop = true;} 
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
        else if (dataSetName == "zjets"){       sprintf(name, "zjets");         xlweight = lumi*3532.8/30364599;    isZjets = true;    } 
        else if (dataSetName == "zjets_lowmll"){sprintf(name, "zjets_lowmll");  xlweight = lumi*860.5/7059426;       isZjets = true;   } 
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
      	
	int idParton_1; 
        int idParton_2; 
        double xParton_1; 
        double xParton_2;
	double Q2_event;
	
				
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
        
        myTree -> Branch("npu", &npu, "npu/I");
        myTree -> Branch("nvertex", &nvertex, "nvertex/D");
        
        myTree -> Branch("metPt", &metPt, "metPt/D");
        myTree -> Branch("metPx", &metPx, "metPx/D");
        myTree -> Branch("metPy", &metPy, "metPy/D");
        
        // Set the branches for the vectors 
	myTree->Branch("idParton_1", &idParton_1,"idParton_1/I");   
        myTree->Branch("idParton_2", &idParton_2,"idParton_2/I"); 
	myTree->Branch("xParton_1", &xParton_1,"xParton_1/D");   
        myTree->Branch("xParton_2", &xParton_2,"xParton_2/D"); 
	myTree->Branch("Q2_event", &Q2_event,"Q2_event/D"); 
	
	
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
      
                 cout << "[Info:] Standard setup " << endl;
		
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
	    
	    
	  	double met_pt = sqrt(met_px*met_px + met_py*met_py);
		double met_pt_0 = sqrt(met_px_0*met_px_0 + met_py_0*met_py_0);
		
		double met_diff = met_pt_0/met_pt; 
		
		
                
		
		
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
								
								idParton_1 = event->idParton1();
								idParton_2 = event->idParton2();
								xParton_1 =  event->xParton1();
								xParton_2 =  event->xParton2();
								Q2_event = event->factorizationScale();
								
//								std::cout << "test:  " << idParton_1 << " " << idParton_2 << " " << xParton_1 << " " << xParton_2 << " " << Q2_event << std::endl;
		      
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
								
								
								


							} // closing the loop over the events with mll > 20 GeV
						} // closing the loop over the events with no extra loose leptons				
					} // closing the loop over the events with two opposite charge leptons
				} //closing the loop over the events with a tight leptons selection	
			} // Closing the loop over the events with good reconstructed primary vertices 
		} //Closing the loop over the triggered events
        } // Closing the loop over the events in one dataset
	
	delete jecUnc;
      	delete jetTools;

    
     fout->Write(); 
      fout->Close(); 
	

	
    }  // closing the loop over the datasets 
    
    
    

    
    // To display how long it took for the program to run
    cout << "*******************************************************************************************" << endl; 
    cout << "It took you " << ((double)clock() - start) /CLOCKS_PER_SEC << " to run the program" << endl; 
  
}

