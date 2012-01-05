///////////////////////////
///// TODO & COMMENTS /////
/////////////////////////// 

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
#include "../JESMeasurement/interface/FullKinFit.h"
#include "../JESMeasurement/interface/JetCombiner.h"
#include "Style.C"

using namespace std;
using namespace TopTree;

int main() {

  clock_t start = clock();

  cout << "*************************************************" << endl;
  cout << " Running The Lepton+Jets RefSel for the SyncEx ! " << endl;
  cout << "*************************************************" << endl;
  
  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

    unsigned int lepton = 1; // 0 = muon+jets 1= electron +jets

  //xml file
  string xmlfile ="../config/myRefSelconfig.xml";

  //Output ROOT file
  string rootFileName ("Esel.root");

  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);

  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<"Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  
  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile.c_str());
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);

  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
  float Luminosity = oldLuminosity;
  for (unsigned int d = 0; d < datasets.size (); d++)
    {
      //cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
      if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
      string dataSetName = datasets[d]->Name();
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	{
	  Luminosity = datasets[d]->EquivalentLumi();
	}
    }
  if(Luminosity != oldLuminosity) cout << "changed analysis environment luminosity to "<< Luminosity << endl;
  

  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* > vertex;
  vector < TRootMuon* > init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* > init_jets;
  vector < TRootMET* > mets;

  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  //Global variable
  TRootEvent* event = 0;
  
  ////////////////////////////////////
  /// Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////

  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;

   ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTable;

  if (lepton == 0) {
    cout << " - Preparing the cutflow table for mu+jets refSelV4" << endl;
    CutsSelecTable.push_back(string("initial"));
    CutsSelecTable.push_back(string("trigged"));
    CutsSelecTable.push_back(string("Good PV"));
    CutsSelecTable.push_back(string("1 selected muon"));
    CutsSelecTable.push_back(string("Veto 2nd muon"));
    CutsSelecTable.push_back(string("Veto electron"));
    char LabelNJets[100];
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
    CutsSelecTable.push_back(string(LabelNJets));
  } else if (lepton == 1){
    cout << " - Preparing the cutflow table for e+jets refSelV4" << endl;
    CutsSelecTable.push_back(string("initial"));
    CutsSelecTable.push_back(string("trigged"));
    CutsSelecTable.push_back(string("Good PV"));
    CutsSelecTable.push_back(string("1 isolated electron"));
    CutsSelecTable.push_back(string("Loose muon veto"));
    CutsSelecTable.push_back(string("Z veto"));
    CutsSelecTable.push_back(string("Conversion rejection"));
    CutsSelecTable.push_back(string("Partnertrack veto"));
    char LabelNJets[100];
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
    CutsSelecTable.push_back(string(LabelNJets));
    sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
    CutsSelecTable.push_back(string(LabelNJets));

  } else {

    cout << "!!!!Lepton should equal 0 or 1, exiting" << endl;

    exit(1);

  }

  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ... size: " << CutsSelecTable.size() <<endl;
  SelectionTable selecTable(CutsSelecTable, datasets);
  //selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;


  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
		
		//open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);

    //selecTable.Fill(d,0, datasets[d]->Xsection() * datasets[d]->EquivalentLumi() );

    string dataSetName = datasets[d]->Name();

    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    //nEvents[d] = 0;
    int itrigger = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
//    for (unsigned int ievt = 0; ievt < 50000; ievt++)
    {

      selecTable.Fill(d,0,1.);

      
      //nEvents[d]++;

      //cout << "anaEnv.JetType -> " << anaEnv.JetType << " --- " <<  endl;
      //cout << "anaEnv.METType -> " << anaEnv.METType << " --- " <<  endl;


      if(ievt%500 == 0)
        std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
        
      //load event

      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

      /////////////
      // TRIGGER //
      /////////////
      
      bool trigged=false;
      
      if (lepton == 0) {
	//if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
	
	int currentRun = event->runId();
	//if(previousRun != currentRun) {
	//  previousRun = currentRun;
	//cout << currentRun << endl;
	itrigger = treeLoader.iTrigger ("HLT_Mu9", currentRun);
	//}
	
	trigged = treeLoader.EventTrigged (itrigger);
	
	//} else 
	//trigged=true;
      }
      
      else if (lepton == 1) {

	trigged=false;
	
	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
	  
	  // SHOULD BE CHECKED ?!
	  std::map<std::string,std::vector<double> > filters= event->getTriggerFilters();
	  
	  if (event->runId() < 140041 && filters.find("hltL1NonIsoHLTNonIsoSingleElectronLWEt10PixelMatchFilter") != filters.end())
	    trigged=true;
	  else if (event->runId() >= 140041 && event->runId() <= 143962 && filters.find("hltL1NonIsoHLTNonIsoSingleElectronEt15PixelMatchFilter") != filters.end())
	    trigged=true;
	  else if (event->runId() >= 143963 && event->runId() <= 146427 && filters.find("hltL1NonIsoHLTNonIsoSingleElectronEt15CaloEleIdPixelMatchFilter") != filters.end())
	    trigged=true;
	  else if (event->runId() >= 146428 && event->runId() <= 147116 && filters.find("hltL1NonIsoHLTNonIsoSingleElectronEt17CaloEleIdPixelMatchFilter") != filters.end())
	    trigged=true;
	  else if (event->runId() >= 147117 && filters.find("hltL1NonIsoHLTNonIsoSingleElectronEt17CaloEleIdPixelMatchFilter") != filters.end())
	    trigged=true;
	  else
	    trigged=false;
	  
	}
	else
	  trigged=true;
      }

      /////////////////////////////
      //   Selection
      /////////////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets);
      
      bool isGoodPV = isGoodPV = selection.isPVSelected(vertex, 4,24,2.);  

      //*******************************//
      //***** Muon+Jets RefSel V4 *****//
      //*******************************//

      if (lepton == 0) {
    
	// FILL THE SELECTION TABLE //
	
	if(trigged){
	  selecTable.Fill(d,1,1.);
	  
	  if(isGoodPV){
	    selecTable.Fill(d,2,1.);
	    
	    vector<TRootJet*> selectedJets;
	    vector<TRootMuon*> selectedMuons;
	   
	    bool doPF2PAT = false; // not supported/synced atm...

	    if (init_jets.size() > 0) {
	      if (init_jets[0]->jetType() == 1 || doPF2PAT) { // calojets
		//cout << "Selecting for caloJets" << endl;
		selectedJets = selection.GetSelectedJets(true);
		selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
	      }	else {
		//cout << "Selecting for PF/JPT jets" << endl;
		vector<TRootMuon*> overlapMuons = selection.GetSelectedMuons(vertex[0]);
		//selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1); // refSelV4 values
		selectedJets = selection.GetSelectedJets(overlapMuons,true);
		selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);	
	      }
	    }

	    vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
	    vector<TRootElectron*> vetoElectrons = selection.GetSelectedLooseElectrons(false);

	    if(selectedMuons.size()==1){
	      selecTable.Fill(d,3,1.);
	      
	      if(vetoMuons.size()==1){
		selecTable.Fill(d,4,1.);
		
		if(vetoElectrons.size()==0) {
		  selecTable.Fill(d,5,1.);			    
		  
		  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3) {
		    selecTable.Fill(d,6,1.);
		    
		    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2) {
		      selecTable.Fill(d,7,1.);
		      
		      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1)
			{
			  selecTable.Fill(d,8,1.);
			  
			  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets)
			    {
			      selecTable.Fill(d,9,1.);

			      // EVENT IS SELECTED!!!!
 
			    }
			}
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      //***********************************//
      //***** ELECTRON+Jets RefSel V4 *****//
      //***********************************//
      
      else if (lepton == 1) {
	
	// FILL THE SELECTION TABLE //
	
	if(trigged){
	  selecTable.Fill(d,1,1.);
	  
	  if(isGoodPV){
	    selecTable.Fill(d,2,1.);
	    
	    vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(vertex[0]);
	    vector<TRootJet*> selectedJets = selection.GetSelectedJets(selectedElectrons,true);
	    vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseElectrons(20,2.5,1.,true);
	    vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
	    
	    if (selectedElectrons.size()==1) {
	      selecTable.Fill(d,3,1.);
	      
	      TRootElectron* electron = (TRootElectron*) selectedElectrons[0];
	      
	      if (histo1D.find("Et") == histo1D.end())
		histo1D["Et"] = new TH1F("Et","Et",300,0,300);
	      
	      histo1D["Et"]->Fill(electron->Et());
	      
	      if (histo1D.find("Eta") == histo1D.end())
		histo1D["Eta"] = new TH1F("Eta","Eta",60,0,6);
	      
	      histo1D["Eta"]->Fill(fabs(electron->Eta()));

	      if (histo1D.find("d0") == histo1D.end())
		histo1D["d0"] = new TH1F("d0","d0",50,0,0.5);
		
	      histo1D["d0"]->Fill(fabs(electron->d0()));
	      
	      if (histo1D.find("RelIso") == histo1D.end())
		histo1D["RelIso"] = new TH1F("RelIso","RelIso",50,0,0.5);
	      
	      float relISO = (electron->caloIso(3)+electron->trackerIso(3)) / electron->Et();
	      histo1D["RelIso"]->Fill(relISO);
	      
	      if (vetoMuons.size()==0) {
		selecTable.Fill(d,4,1.);

		bool passZVeto = true;
		for (unsigned int e=0;e<looseElectrons.size();e++) {
		  
		  TRootElectron* el = (TRootElectron*) looseElectrons[e];
		  
		  if (fabs(el->superClusterEta()) > 1.5660 || fabs(el->superClusterEta()) < 1.4442) {
		    TLorentzVector Zcand = *looseElectrons[e]+*selectedElectrons[0];
		    
		    //cout << Zcand.M() << endl;
		    if (Zcand.M() > 76 && Zcand.M() < 106)
		      passZVeto = false;
		  }
		}
		
		if (passZVeto) {
		  selecTable.Fill(d,5,1.);
		  
		  if (selectedElectrons[0]->missingHits() == 0) {
		    selecTable.Fill(d,6,1.);
		    		    
		    if (fabs(selectedElectrons[0]->Dist()) >= 0.02 || fabs(selectedElectrons[0]->DCot()) >= 0.02) {
		      selecTable.Fill(d,7,1.);
			
		      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3) {
			selecTable.Fill(d,8,1.);
			
			if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2) {
			  selecTable.Fill(d,9,1.);
			  
			  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1) {
			    selecTable.Fill(d,10,1.);
			    
			    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets) {
			      selecTable.Fill(d,11,1.);

			      // EVENT IS SELECTED

			    }
			  }
			}
		      }
		    }
		  }
		}	
	      }
	    }
	  }
	}
      }
      
    } // event loop
    
  } // dataset loop

  //Selection tables
  selecTable.TableCalculator(false, true, true, true);
  string selectiontable = "SelectionTable_SyncEx";
  selectiontable = selectiontable +".tex"; 	
  selecTable.Write(selectiontable.c_str());


  //TFile* fout = new TFile("out.root","RECREATE");

  fout->cd();

  TDirectory* th1dir = fout->mkdir("1D_histograms");

  th1dir->cd();

  fout->cd();

  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++) {
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    //temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    //temp->SetBinContent(N+1,0);
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( ("plots/"+it->first+".png").c_str() );
  }
  
  fout->Close();

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

}
 
