// isis.marina.van.parijs@cern.ch 
// 2013
// This is a program that runs over the toptrees and calculates the 
// efficiencies of certain cuts in the datasamples. 

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

	// bool for debugging
	bool debug = false; 
	bool warnings = true; 
	bool information = true; 

        //set the xml file
	string xmlfile = "FCNC_config.xml";     //place of the xml file 
	
	//set the channel 
	string channel = "undefined";
	
	
	//set a default luminosity in pb^-1
	float luminosity = 100000; 
	float NofEvts = 100000;

	

	
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
			cout << "--xml myxml.xml: change Xml file" << endl; 
			cout << "--3gamma: use the 3 gamma channel" << endl; 
			cout << "--1L3B: use the 1 lepton + 3 b-tags channel" << endl; 
			cout << "--SSdilepton: use the same sign dilepton channel" << endl; 
			cout << "--3L: use the channel with at least 3 leptons" << endl;
                	return 0;
        	}
		if (argval=="--xml") {
			iarg++;
			tempxml = argv[iarg];
			foundxml = true; 
		}
		if (argval=="--3gamma") {
                	channel = "3gamma";
			xmlfile = "FCNC_3gamma_config.xml";
        	}
		if (argval=="--1L3B") {
                	channel = "1L3B";
			xmlfile = "FCNC_1L3B_config.xml";
        	}
		if (argval=="--SSdilepton") {
                	channel = "SSdilepton";
			xmlfile = "FCNC_SSdilepton_config.xml";
        	}
		if (argval=="--3L") {
                	channel = "3L";
			xmlfile = "FCNC_3L_config.xml";
        	}

    	} 
	
    	if (foundxml)
	{
		xmlfile = tempxml; 
	}
	if(information)	std::cout << "[INFO]	Used configuration file: " << xmlfile << endl;
	if(information)	std::cout << "[INFO]	Used channel: " << channel << endl;
	if(channel.find("undefined")!=string::npos && warnings) std:cout << "[WARNING]	No channel was defined" << endl; 
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//   end different options for executing this macro      //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////

	//Load the analysisenvironment
	AnalysisEnvironment anaEnv; 
	if(debug) std::cout << "[PROCES]	Loading the analysisenvironment" << endl; 
	AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());    //load via the xml file the environment
	
	
	//Load the datasets
	TTreeLoader treeLoader; 
	vector <Dataset*> datasets; //vector that will contain all datasets
	if(debug) std::cout << "[PROCES]	Loading the datasets " << endl; 
	treeLoader.LoadDatasets(datasets, xmlfile.c_str()); //put datasets via xmlfile in the dataset vector

	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                       output stuff                    //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	
	// Set an output rootfile
	char rootFileName[900];
	char channelchar[900];
	if(channel.find("3gamma")!=string::npos)	sprintf(channelchar, "3gamma");
	if(channel.find("3L")!=string::npos)		sprintf(channelchar, "3L");	
	if(channel.find("1L3B")!=string::npos)		sprintf(channelchar, "1L3B");
	if(channel.find("SSdilepton")!=string::npos)	sprintf(channelchar, "SSdilepton");
	
	sprintf(rootFileName,"Output_FCNC_selection_%s.root",channelchar);
	TFile *fout = new TFile (rootFileName, "RECREATE");
	if(debug) cout << "[PROCES]	Declared output rootfiles  "<< endl;
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                        end output stuff               //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	// Declare variables: 
	if(debug) cout << "[PROCES]	Variable declaration  "<< endl;
  	vector < TRootVertex* >   vertex;
  	vector < TRootMuon* >     init_muons;
  	vector < TRootElectron* > init_electrons;
  	vector < TRootJet* >      init_jets;
  	vector < TRootMET* >      mets;
	
	//Define an event (global variable)
	TRootEvent* event = 0;
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                Cut flow histograms		        //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
		
	char plotTitle_total[900];
	sprintf(plotTitle_total,"The total cutflow for %s channel",channelchar); 
	
	TH1F* cutflow_total = new TH1F("cutflow_total", plotTitle_total, 6, -0.5,5.5); 
	cutflow_total->Sumw2();
	
	if(debug) cout << "[PROCES]	Declared total cutflow histogram  "<< endl;
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                START LOOPING OVER THE DATASETS        //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	if(information)	cout << "[PROCES]	Looping over the datasets:  " << datasets.size()<< " datasets" << endl;
	for(unsigned int d = 0; d < datasets.size();d++)
	{
		//Load datasets
		treeLoader.LoadDataset(datasets[d], anaEnv); 
		string datasetName = datasets[d]->Name(); 
		
		char datasetNamechar[900];
		if(datasetName.find("ttbar")!=string::npos) {sprintf(datasetNamechar,"ttbar");}
		if(datasetName.find("Wjets")!=string::npos) {sprintf(datasetNamechar,"wjets");}
		if(datasetName.find("ttt")!=string::npos) {sprintf(datasetNamechar,"ttt");}
		if(datasetName.find("ttW")!=string::npos) {sprintf(datasetNamechar,"ttw");}
		if(datasetName.find("WZ")!=string::npos) {sprintf(datasetNamechar,"wz");}
		if(datasetName.find("ZZ")!=string::npos) {sprintf(datasetNamechar,"zz");}
		if(datasetName.find("ttZ")!=string::npos) {sprintf(datasetNamechar,"ttz");}
		
		
		if(information) cout << "[INFO]	Dataset " << d << " name : " << datasetName << " / title : " << datasets[d]->Title() << endl;
      			
		
		
		// Define different plots for each channel and dataset
		char plotTitle[900];
		char NamePlot[900];
		sprintf(plotTitle,"The cutflow for %s channel: %s dataset",channelchar,datasetNamechar); 
		sprintf(NamePlot,"cutflow_%s",datasetNamechar);
	
		TH1F* cutflow = new TH1F(NamePlot, plotTitle, 6, -0.5,5.5); 
		cutflow->Sumw2();
	
		if(debug) cout << "[PROCES]	Declared cutflow histogram  "<< endl;
		
		
		///////////////////////////////////////////////////////////
		//                START LOOPING OVER THE EVENTS          //
		///////////////////////////////////////////////////////////

		float NofRuns = 0; 
		if( NofEvts > datasets[d]->NofEvtsToRunOver()) 
		{
			NofRuns = datasets[d]->NofEvtsToRunOver(); 
		}
		else
		{
			NofRuns = NofEvts; 
		} 
		
		if(information) cout << "[PROCES]	looping over " << NofRuns <<" events "<< endl;
		
		for(int ievent = 0; ievent <NofRuns; ievent++)
		{
			if(ievent%1000 == 0 && information)
			{
				// << flush << "\r" means this line will be overwritten next time 
				std::cout << "[PROCES]	Processing the " << ievent << "th event" << flush << "\r";  
			}    
			
			
			//Load the event 
			event = treeLoader.LoadEvent(ievent, vertex, init_muons, init_electrons, init_jets, mets);
			
			
			cutflow->Fill(1);
			cutflow_total->Fill(1);
			
			cutflow_total->GetXaxis()->SetBinLabel(2, "initial");
			cutflow->GetXaxis()->SetBinLabel(2, "initial");
			
			//Make a selection 
			Selection selection(init_jets, init_muons,init_electrons,mets);
			//define selection cuts --> have to be validated!!!
			// From the class Selection the following functions are used: 
				// 	void Selection::setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon)
				//	void Selection::setDiElectronCuts(float Et, float Eta, float RelIso, float d0, float MVAId, float DistVzPVz, float DRJets, int MaxMissingHits) 
				//	void Selection::setLooseDiElectronCuts(float ptt, float Eta, float RelIso, MVAid) 
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
			// vector<TRootPhoton*> selectedPhotons = selection.GetSelecetedPhotons(); Photons not yet included in the selection class!!!!
			
			
			//order the jets according to the Pt 
			sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order jets wrt Pt.  
			
			
			
			//Create a vector containing all the bjets (since it is simulation, I can know this)
			vector <int> mcParticlesTLV,selectedJetsTLV;
			float bjets = 0;
			
					
			mcParticlesTLV.clear(); //make sure nothing is inside this vector
			selectedJetsTLV.clear();
			
					
      			for(unsigned int iJet=0;iJet<selectedJets.size(); iJet++){
				TRootJet* tempJet = (TRootJet*) selectedJets[iJet];
						
				int pdgID = tempJet->partonFlavour();
						
				selectedJetsTLV.push_back(iJet);
				mcParticlesTLV.push_back(pdgID);
				
				if(fabs(pdgID) == 5)
				{
	      				bjets++;
								
	 			}
			
			} 

			
			if(debug) cout << "looseElectrons.size() = " << looseElectrons.size() << endl; 
			if(debug) cout << "looseMuons.size() = " << looseMuons.size() << endl; 
			
			//more than 3 leptons
			if(channel.find("3L")!=string::npos)
			{
				if(debug) cout << "in 3L channel" << endl;
				if(looseElectrons.size() > 2  || looseMuons.size() > 2)
				{ 
					if(debug) cout << "fill 3L" << endl;
					cutflow_total->Fill(2);
					cutflow->Fill(2);
					cutflow_total->GetXaxis()->SetBinLabel(3, "3L");
					cutflow->GetXaxis()->SetBinLabel(3, "3L");
					if(debug) cout << "filled 3L" << endl;
				}
				if(debug)	cout << "out fill 3L loop" << endl; 
			}
			//1 lepton + 3 b-jets
			if(channel.find("1L3B")!=string::npos)
			{
				if(debug) cout << "in 1L3B channel" << endl;
				if(looseElectrons.size() +  looseMuons.size() == 1)
				{
					if(debug) cout << "in fill 1l3b loop" << endl;
					cutflow_total->Fill(2);
					cutflow->Fill(2);
					cutflow_total->GetXaxis()->SetBinLabel(3, "1L");
					cutflow->GetXaxis()->SetBinLabel(3, "1L");
					if(debug) cout << "selectedJets.size() = " << selectedJets.size() << endl;
					if(selectedJets.size() == 3)
					{
						if(debug) cout << "in fill 1l3b loop: 3jets" << endl;
						cutflow_total->Fill(3);
						cutflow->Fill(3);
						cutflow_total->GetXaxis()->SetBinLabel(4, "3jets");
						cutflow->GetXaxis()->SetBinLabel(4, "3jets");
						if(bjets == 3)
						{
							if(debug) cout << "in fill 1l3b loop: 3bjets" << endl;
							cutflow_total->Fill(4);
							cutflow->Fill(4);
							cutflow_total->GetXaxis()->SetBinLabel(5, "3 bjets");
							cutflow->GetXaxis()->SetBinLabel(5, "3 bjets");
						}
					}
				
					if(debug) cout << "out fill 1l3b loop" << endl;
				}
			}
			if(channel.find("SSdilepton")!=string::npos)
			{
				if(debug) cout << "in SSdilepton channel" << endl;
				if(looseElectrons.size() + looseMuons.size() == 2)
				{
					if(debug) cout << "in fill SS dilepton " << endl; 
					cutflow_total->Fill(2);
					cutflow->Fill(2);
					
					cutflow_total->GetXaxis()->SetBinLabel(3, "2L");
					cutflow->GetXaxis()->SetBinLabel(3, "2L");
					
					bool electron = false; 
					bool muon = false; 
		  			TRootElectron* electron0 = 0;
					TRootElectron* electron1 = 0;
					TRootMuon* muon0 = 0;
					TRootMuon* muon1 = 0;
				
					if(looseElectrons.size() > 1)
					{
						for(unsigned int i = 0; i<looseElectrons.size()-1; i++)
						{
							for(unsigned int j = i+1;j<looseElectrons.size();j++)
							{
								if(debug) cout << "in fill SS dilepton: electronloop " << endl;
								electron0 = (TRootElectron*) looseElectrons[i];
								electron1 = (TRootElectron*) looseElectrons[j];
							}
						}
						if(electron0->charge()*electron1->charge()>0) electron = true; 
						if(debug) cout << "Electron boolean defined" << endl; 
					}
				
					if(looseMuons.size() > 1)
					{
						for(unsigned int i = 0; i<looseMuons.size()-1; i++)
						{
							for(unsigned int j = i+1;j<looseMuons.size();j++)
							{
								if(debug) cout << "in fill SS dilepton: muonloop " << endl;
								muon0 = (TRootMuon*) looseMuons[i];
								muon1 = (TRootMuon*) looseMuons[j];
							}
						}
						if(muon0->charge()== muon1->charge()) muon = true; 
					}
				
					if(muon || electron)
					{
						if(debug) cout << "in fill SS dilepton: same sign " << endl;
						cutflow_total->Fill(3);
						cutflow->Fill(3);
						cutflow_total->GetXaxis()->SetBinLabel(4, "2 SS L");
						cutflow->GetXaxis()->SetBinLabel(4, "2 SS L");
					}
					if(debug) cout << "out fill SS dilepton " << endl;
				}
			}
			if(channel.find("3gamma")!=string::npos)
			{
				cout << "Photons are not yet included in selection class !!! " << endl; 
			/*	if(debug) cout << "in fill 3gamma " << endl;
				TRootPhoton* photon0 = 0;
				TRootPhoton* photon1 = 0;
				TRootPhoton* photon3 = 0;
				
				if(selectedPhotons.size() > 2)
				{
					if(debug) cout << "in fill 3 gamma: selected photons loop " << endl;
					cutflow->Fill(2);
					cutflow_total->Fill(2)
					
					cutflow_total->GetXaxis()->SetBinLabel(3, "3 photons");
					cutflow->GetXaxis()->SetBinLabel(3, "3 photons");
				}	
				if(debug) cout << "out fill 3gamma " << endl;
			*/
			}
			                                                                  
    			
	
	
		}
		
		///////////////////////////////////////////////////////////
		//                END LOOPING OVER THE EVENTS            //
		///////////////////////////////////////////////////////////
		
		
	
	}
	if(information)	cout << "[PROCES]	End of looping over the datasets:  " << datasets.size()<< " datasets" << endl;
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                END LOOPING OVER THE DATASETS          //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	fout ->Write(); 
	fout->Close();
	
	
	
	
	 
	
	std::cout << "******************************************"<<std::endl; 
	std::cout << " End of the program for the FCNC selection " << std::endl; 
	std::cout << "******************************************"<<std::endl;

	// display how long it took for the program to run
	std::cout << "******************************************" << std::endl; 
	std::cout << " It took the program " << ((double)clock() - start) /CLOCKS_PER_SEC << " to run. " << endl; 
	std::cout << "******************************************" << std::endl; 
	return 0; 
	
}
