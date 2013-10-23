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


	//see how long the program takes to run with a clock 
	clock_t start = clock(); 

	std::cout << "******************************************"<<std::endl; 
	std::cout << " Beginning of the program for the FCNC selection " << std::endl; 
	std::cout << "******************************************"<<std::endl; 

        //set the xml file
	string xmlfile = "FCNC_config.xml";     //place of the xml file 
	std::cout << "Used configuration file: " << xmlfile << endl; 

	//Load the analysisenvironment
	AnalysisEnvironment anaEnv; 
	std::cout << "Loading the analysisenvironment" << endl; 
	AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());    //load via the xml file the environment
	//new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);     //what does this do? 
	
	//Load the datasets
	TTreeLoader treeLoader; 
	vector <Dataset*> datasets; //vector that will contain all datasets
	std::cout << "Loading the datasets " << endl; 
	treeLoader.LoadDatasets(datasets, xmlfile.c_str()); //put datasets via xmlfile in the dataset vector


	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//            BEGIN RUNNING OVER THE DATASETS            //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	for(unsigned int d = 0; d < datasets.size(); d++){
	
	
	
	
	
	
	
	
	}
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//            END RUNNING OVER THE DATASETS            //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////

	std::cout << "******************************************"<<std::endl; 
	std::cout << " End of the program for the FCNC selection " << std::endl; 
	std::cout << "******************************************"<<std::endl;

	// display how long it took for the program to run
	std::cout << " ******************************************" << std::endl; 
	std::cout << " It took the program " << ((double)clock() - start) /CLOCKS_PER_SEC << " to run. " << endl; 
	
	return 0; 
	
}
