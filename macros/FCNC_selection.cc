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
#include "../../TopTreeAnalysisBase/MCInformation/interface/Lumi3DReWeighting.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

#include "../../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"



#include "Style.C"

using namespace std;
using namespace TopTree;




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
	std::cout << "[INFO]	Used configuration file: " << xmlfile << endl; 
	
	//set a default luminosity in pb^-1
	float luminosity = 100000; 

	//Load the analysisenvironment
	AnalysisEnvironment anaEnv; 
	std::cout << "[PROCES]	Loading the analysisenvironment" << endl; 
	AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());    //load via the xml file the environment
	int verbose = 2;
	
	//Load the datasets
	TTreeLoader treeLoader; 
	vector <Dataset*> datasets; //vector that will contain all datasets
	std::cout << "[PROCES]	Loading the datasets " << endl; 
	treeLoader.LoadDatasets(datasets, xmlfile.c_str()); //put datasets via xmlfile in the dataset vector



	




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
		std::cout<< "[INFO]	Found dataset " << datasetName << " with an equivalent luminosity of " << datasets[d]->EquivalentLumi() << std::endl;  
		
		
		if(datasetName.find("data")!=string::npos || datasetName.find("Data")!=string::npos || datasetName.find("DATA")!=string::npos)
		{
			luminosity = datasets[d]->EquivalentLumi();
			break; 
		}
		
		std::cout<<"[INFO]	Rescaling to an integrated luminosity of " << luminosity << " pb^-1" << endl; 
	

		

	}
	std::cout<<"[INFO]	Rescaled to an integrated luminosity of " << luminosity << std::endl; 
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
	cout << "[PROCES]	MSplot declaration  "<< endl;
	
	
	
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
	cout << "[PROCES]	1D plot declaration  "<< endl;
	
	
	
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
	//Define a table where all kinematic cuts are stored eg > 2jets)
	vector<string> CutsSelectionTable;
	CutsSelectionTable.push_back(string("initial"));
	
	//Define a selection table with the previously defined kin. cuts and the given datasets
	SelectionTable Selectiontable(CutsSelectionTable, datasets);
	// give it the rescaled luminosity 
	Selectiontable.SetLuminosity(luminosity); 
	
	 
	 
	
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                END OF SELECTIONTABLES                 //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                START LOOPING OVER THE EVENTS          //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	cout << "[PROCES]	Looping over the datasets:  " << datasets.size()<< " datasets" << endl;
		



	std::cout << "******************************************"<<std::endl; 
	std::cout << " End of the program for the FCNC selection " << std::endl; 
	std::cout << "******************************************"<<std::endl;

	// display how long it took for the program to run
	std::cout << "******************************************" << std::endl; 
	std::cout << " It took the program " << ((double)clock() - start) /CLOCKS_PER_SEC << " to run. " << endl; 
	std::cout << "******************************************" << std::endl; 
	return 0; 
	
}
