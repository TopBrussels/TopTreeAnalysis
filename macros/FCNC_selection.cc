// isis.marina.van.parijs@cern.ch 
// 2013
// This is a program that runs over the toptrees

#include "TStyle.h"
#include "TH3F.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

#include "../../TopTreeProducer/interface/TRootRun.h"
#include "../../TopTreeProducer/interface/TRootEvent.h"

#include "../../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"

#include "../../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../../TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "../../TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "../../TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"

#include "../../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../../TopTreeAnalysisBase/Content/interface/Dataset.h"

#include "../../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

#include "../../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"



#include "Style.C"

using namespace std;	//needed for cout and stuff
using namespace TopTree;	//needed for TT
using namespace reweight;  //needed for PUreweighting




int main(int argc, char *argv[]){
	//Make plots nicer: color, style, ... 
	setMyStyle();

	//see how long the program takes to run with a clock 
	clock_t start = clock(); 
	std::cout << "******************************************" << std::endl; 
	std::cout << " Starting clock" << endl; 
	std::cout << "******************************************"<<std::endl; 
	std::cout << " Beginning of the program for the FCNC selection " << std::endl; 
	std::cout << "******************************************"<<std::endl; 




        //set the xml file
	string xmlfile = "FCNC_config.xml";     //place of the xml file 
	
	
	//set a default luminosity in pb^-1
	float luminosity = 100000; 

	//Load the analysisenvironment
	AnalysisEnvironment anaEnv; 
	std::cout << "[PROCES]	Loading the analysisenvironment" << endl; 
	AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());    //load via the xml file the environment
	//int verbose = 2; // why do we need this? 
	
	//Load the datasets
	TTreeLoader treeLoader; 
	vector <Dataset*> datasets; //vector that will contain all datasets
	std::cout << "[PROCES]	Loading the datasets " << endl; 
	treeLoader.LoadDatasets(datasets, xmlfile.c_str()); //put datasets via xmlfile in the dataset vector


	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//             systematics   booleans                    //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	// - PU reweighing
	bool reweighPU = false; 
	
	// - naked option for raw data 
	bool isRAW = false; 
	
	
	// - HLT 
	bool runHLT = false; 
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//   different options for executing this macro          //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	std::string tempxml; 
	bool foundxml = false; 
    
    	for(int iarg = 0; iarg < argc && argc>1 ; iarg++)
	{
        	std::string argval=argv[iarg];
		
        	if(argval=="--help" || argval =="--h")
		{
                	cout << "--NoPU: Do not apply pileup re-weighting" << endl;
                	// cout << "--RAW: Do not apply pileup re-weighting or b-tag scale factor" << endl;
			cout << "--xml myxml.xml: change Xml file" << endl; 
			cout << "--noHLT: Do not apply HLT" << endl;
                	return 0;
        	}
        	if (argval=="--NoPU") {
                	reweighPU = false;
        	}
                //if (argval=="--RAW") {
                //	reweighPU = false;  
                //	isRAW = true;
        	//}
       		if (argval=="--xml") {
                	iarg++;
			tempxml = argv[iarg];
			foundxml = true; 
        	}
		if (argval=="--noHLT") {
                	runHLT = false; 
        	}
    	}   
    
    	if (foundxml)
	{
		xmlfile = tempxml; 
	}
	std::cout << "[INFO]	Used configuration file: " << xmlfile << endl; 
	
	if(!reweighPU ||isRAW ||!runHLT)
	{	
		cout << "[INFO]	You are using NON standard options:" << endl;
		if(!reweighPU)	cout << "[INFO]	- You are NOT applying PU reweighting" << endl; 
		if(isRAW)	cout << "[INFO]	- You are using raw data" << endl; 
		if(!runHLT)	cout << "[INFO]	- You are NOT applying HLT " << endl; 
	
	}
	else
	{
		cout << "[INFO]	Using standard set up" << endl; 
	}

	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//   end different options for executing this macro      //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////




	//Initialize PUreweighting
	LumiReWeighting LumiWeights; 
	cout << "[PROCES]	Initialized PU reweighting" << endl; 








	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//            BEGIN REWEIGHTING THE DATASETS            //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	for(unsigned int d = 0; d < datasets.size(); d++){
		//Take the dataset on place d in the datasetvector with the given
		//analysis environment
		treeLoader.LoadDataset(datasets[d], anaEnv); 
		string datasetName = datasets[d]->Name(); 
		
		//define a reweighting weight for the MC sample, 1 means no reweighting 
		// the equivalent luminosity is given with the dataset in the xml file, where the equiv lumi is defined as the #evts before skimming (this you get from the TT) divided 
		// by the cross-section of that event. For data, the equivalent lumi is the real lumi of the experiment. 
		std::cout<< "[INFO]	Found dataset " << datasetName << " with an equivalent luminosity of " << datasets[d]->EquivalentLumi() << " pb^-1"<< std::endl;  
		
		
		if(datasetName.find("data")!=string::npos || datasetName.find("Data")!=string::npos || datasetName.find("DATA")!=string::npos)
		{
			luminosity = datasets[d]->EquivalentLumi();
			reweighPU = false; 
			break; 
		}
		
		std::cout<<"[INFO]	Rescaling to an integrated luminosity of " << luminosity << " pb^-1" << endl; 
	

		

	}
	std::cout<<"[INFO]	Rescaled to an integrated luminosity of " << luminosity << " pb^-1"<< std::endl; 
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//            END REWEIGHTING THE DATASETS               //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	
	
	// Set an output rootfile
	string OutputRootFile("Output_FCNC_Selection.root"); 
	// Open the created rootfile and RECREATE:  create a new file, if the file already exists it will be overwritten.
	TFile *outputFile = new TFile(OutputRootFile.c_str(),"RECREATE"); 
	
	
	// Declare variables: 
	cout << "[PROCES]	Variable declaration  "<< endl;
  	vector < TRootVertex* >   vertex;
  	vector < TRootMuon* >     init_muons;
  	vector < TRootElectron* > init_electrons;
  	vector < TRootJet* >      init_jets;
  	vector < TRootMET* >      mets;
	
	//Define an event (global variable)
	TRootEvent* event = 0;
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//            DEFINING THE MULTISAMPLEPLOTS             //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////	
	//cout << "[PROCES]	MSplot declaration  "<< endl;
	
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//      END OF DEFINING THE MULTISAMPLEPLOTS             //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////	
	
	

	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                   DEFINING THE 1D PLOTS               //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////	
	//cout << "[PROCES]	1D plot declaration  "<< endl;
	
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//          END OF DEFINING THE 1D PLOTS                 //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                    SELECTIONTABLES                    //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	cout << "[PROCES]	Selection table declaration  "<< endl;	
	//Define a table where all kinematic cuts are stored eg > 2jets  to make at the end a cutflow table)
	vector<string> CutsSelectionTable;
	CutsSelectionTable.push_back(string("initial"));
	
	//Define a selection table with the previously defined kin. cuts and the given datasets
	SelectionTable Selectiontable(CutsSelectionTable, datasets);
	// give it the rescaled luminosity such that the number of events corresponds with the right luminosity 
	Selectiontable.SetLuminosity(luminosity); 
	// set the precision at 1 decimal 
	Selectiontable.SetPrecision(1);
	
	 
	 
	
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                END OF SELECTIONTABLES                 //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                START LOOPING OVER THE DATASETS        //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	cout << "[PROCES]	Looping over the datasets:  " << datasets.size()<< " datasets" << endl;
	for(unsigned int d = 0; d < datasets.size();d++)
	{
		//Load datasets
		treeLoader.LoadDataset(datasets[d], anaEnv); 
		string datasetName = datasets[d]->Name(); 
		
		//if(verbose > 1)
		//{
			cout << "[INFO]	Dataset " << d << " name : " << datasetName << " / title : " << datasets[d]->Title() << endl;
      			cout << "[INFO]	Cross section = " << datasets[d]->Xsection() << endl;
      			cout << "[INFO]	IntLumi = " << datasets[d]->EquivalentLumi() << "  NormFactor = " << datasets[d]->NormFactor() << endl;
      			cout << "[INFO]	Nb of events : " << datasets[d]->NofEvtsToRunOver() << endl;
		
		//}
		
		
		///////////////////////////////////////////////////////////
		//                START LOOPING OVER THE EVENTS          //
		///////////////////////////////////////////////////////////
		
		//if(verbose > 1)
		//{
			cout << "[PROCES]	looping over the events: " << datasets[d]->NofEvtsToRunOver() << " events." << endl; 
		//}
		
		for(int ievent = 0; ievent < datasets[d]->NofEvtsToRunOver(); ievent++)
		{
			if(ievent%1000 == 0)
			{
				// << flush << "\r" means this line will be overwritten next time 
				std::cout << "[PROCES]	Processing the " << ievent << "th event" << flush << "\r";  
			}      
			
			//Load the event 
			event = treeLoader.LoadEvent(ievent, vertex, init_muons, init_electrons, init_jets, mets);
			
			
			//Fill the selection table
			Selectiontable.SetLuminosity(luminosity);
			Selectiontable.Fill(d,0,1);  // Fill the initial number of events in the cutflow table on row 0
			
			
			//Make a selection 
			Selection selection(init_jets, init_muons,init_electrons,mets);
			//define selection cuts --> have to be validated!!!
			// From the class Selection the following functions are used: 
				// 	void Selection::setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon)
				//	void Selection::setDiElectronCuts(float Et, float Eta, float RelIso, float d0, float MVAId, float DistVzPVz, float DRJets, int MaxMissingHits) 
				//	void Selection::setLooseDiElectronCuts(float Et, float Eta, float RelIso) 
				//	void Selection::setDiMuonCuts(float Pt, float Eta, float RelIso, float d0) 
				// 	void Selection::setLooseMuonCuts(float Pt, float Eta, float RelIso) 
			selection.setJetCuts(20.,5.,0.01,1.,0.98,0.3,0.1);
			selection.setDiMuonCuts(20.,2.4,0.20,999.);
              		selection.setDiElectronCuts(20.,2.5,0.15,0.04,0.5,1,0.3,1); 
              		selection.setLooseMuonCuts(10.,2.5,0.2);
              		selection.setLooseDiElectronCuts(15.0,2.5,0.2,0.5); 
			
			//select the right objects and put them in a vector
			vector<TRootJet*> selectedJets = selection.GetSelectedJets(true);
			vector<TRootMuon*> selectedMuons = selection.GetSelectedDiMuons();
	      		vector<TRootMuon*> looseMuons = selection.GetSelectedLooseMuons();
	      		vector<TRootElectron*> selectedElectrons = selection.GetSelectedDiElectrons();
	      		vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseDiElectrons();
			
			
			//order the jets according to the Pt 
			sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order jets wrt Pt.                                                                    
    			
	
	
		}
		
		///////////////////////////////////////////////////////////
		//                END LOOPING OVER THE EVENTS            //
		///////////////////////////////////////////////////////////
		
	
	
	}
	cout << "[PROCES]	End of looping over the datasets:  " << datasets.size()<< " datasets" << endl;
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                END LOOPING OVER THE DATASETS          //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	
	std::cout << "******************************************"<<std::endl; 
	std::cout << " End of the program for the FCNC selection " << std::endl; 
	std::cout << "******************************************"<<std::endl;

	// display how long it took for the program to run
	std::cout << "******************************************" << std::endl; 
	std::cout << " It took the program " << ((double)clock() - start) /CLOCKS_PER_SEC << " to run. " << endl; 
	std::cout << "******************************************" << std::endl; 
	return 0; 
	
}
