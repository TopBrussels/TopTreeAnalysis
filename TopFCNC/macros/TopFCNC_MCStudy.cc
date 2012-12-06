#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../Tools/interface/JetTools.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../Reconstruction/interface/MakeBinning.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../TopFCNC/interface/TopFCNC_GenEvt.h"

#include "../../macros/Style.C"

using namespace std;
using namespace TopTree;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;


int main (int argc, char *argv[])
{

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FCNC MC study ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
  //setMyStyle();

  string postfix = "_MCStudy"; // to relabel the names of the output file

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string xmlFileName = "../config/myTopFCNCconfig_MCStudy.xml";
  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////// AnalysisEnvironment /////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = 2;//anaEnv.Verbose;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  

  //Output ROOT file
  string rootFileName ("TopFCNC"+postfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //vector of objects
  cout << " - Variable declaration ..." << endl;

  vector < TRootVertex* >   vertex;
  vector < TRootMuon* >     init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* >      init_jets;
  vector < TRootMET* >      mets;

  //Global variable
  TRootEvent* event = 0;
  TopFCNC_GenEvt* MyTopFCNC_GenEvtCand = 0;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  histo1D["DiLeptonMass"] = new TH1F("DiLeptonMass",";m_{ll};#events",320,50,130);

  histo1D["TopMass_FCNC_Decay_DiLep"] = new TH1F("TopMass_FCNC_Decay_DiLep",";m_{Zq};#events",500,100,250);
  histo1D["TopMass_Had_Decay_DiLep"]  = new TH1F("TopMass_Had_Decay_DiLep",";m_{bqq};#events",500,100,250);
  histo1D["WMass_Had_Decay_DiLep"]    = new TH1F("WMass_Had_Decay_DiLep",";m_{qq};#events",480,0,160);

  histo1D["TopMass_FCNC_Decay_TriLep"] = new TH1F("TopMass_FCNC_Decay_TriLep",";m_{Zq};#events",500,100,250);
  histo1D["TopMass_Had_Decay_TriLep"]  = new TH1F("TopMass_Had_Decay_TriLep",";m_{blv};#events",500,100,250);
  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "TopFCNC"+postfix;
  pathPNG += "_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Loop on datasets ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

  for (unsigned int d = 0; d < datasets.size(); d++) //d < datasets.size()
  {
    if (verbose > 1){
      cout << "   Dataset " << d << " name : " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;
      cout << " - Cross section = " << datasets[d]->Xsection() << endl;
      cout << " - IntLumi = " << datasets[d]->EquivalentLumi() << "  NormFactor = " << datasets[d]->NormFactor() << endl;
      cout << " - Nb of events : " << datasets[d]->NofEvtsToRunOver() << endl;
    }
    //open files and load
    cout<<"Load Dataset"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
		
    string dataSetName = datasets[d]->Name();	
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////// Loop on events //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();
     
    if (verbose > 1) cout << " - Loop over events " << endl;      
    
    for (unsigned int ievt = start; ievt < end; ievt++)
    {        

      if(ievt%1000 == 0)
		  std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

	    //load event
	    event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
	    //vector<TRootGenJet*>     genjets     = treeLoader.LoadGenJet(ievt);
	    vector<TRootMCParticle*> mcParticles   = treeLoader.LoadMCPart(ievt);

	    Selection selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...
	    selection.setJetCuts(10.,2.4,0.01,1.,0.98,0.3,0.1);
	   	selection.setDiMuonCuts(15,2.4,0.20,999.); //Et,Eta,RelIso,d0
  	  selection.setDiElectronCuts(15,2.5,0.15,0.04,0.,1,0.3,1); //Et,Eta,RelIso,d0,MVAId,DistVzPVz,DRJets,MaxMissHits

	    vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId
	    vector<TRootMuon*>     selectedMuons       = selection.GetSelectedDiMuons();

	    if(selectedMuons.size()<2) continue;

	    MyTopFCNC_GenEvtCand = new TopFCNC_GenEvt();
	    MyTopFCNC_GenEvtCand->ReconstructEvt(mcParticles);
	    MyTopFCNC_GenEvtCand->MatchJetsToPartons(selectedJets);
	    MyTopFCNC_GenEvtCand->MatchLeptonsToZ(selectedMuons);
      
	    histo1D["DiLeptonMass"]->Fill(MyTopFCNC_GenEvtCand->matchedZ().M());
	    if(MyTopFCNC_GenEvtCand->isDiLeptonic()){
		    histo1D["TopMass_FCNC_Decay_DiLep"]->Fill((MyTopFCNC_GenEvtCand->matchedZ()+MyTopFCNC_GenEvtCand->matchedQ()).M());
		    if(MyTopFCNC_GenEvtCand->matchedQuark1FromW().E()>0 && MyTopFCNC_GenEvtCand->matchedQuark2FromW().E()>0){
			    histo1D["WMass_Had_Decay_DiLep"]->Fill((MyTopFCNC_GenEvtCand->matchedQuark1FromW()+MyTopFCNC_GenEvtCand->matchedQuark2FromW()).M());
			    if(MyTopFCNC_GenEvtCand->matchedB().E()>0)
				    histo1D["TopMass_Had_Decay_DiLep"]->Fill((MyTopFCNC_GenEvtCand->matchedQuark1FromW()+MyTopFCNC_GenEvtCand->matchedQuark2FromW()+MyTopFCNC_GenEvtCand->matchedB()).M());
		    }
	    }
	    else{
		    histo1D["TopMass_FCNC_Decay_TriLep"]->Fill((MyTopFCNC_GenEvtCand->matchedZ()+MyTopFCNC_GenEvtCand->matchedQ()).M());
	    }

	    //delete selection;
	    if(MyTopFCNC_GenEvtCand) delete MyTopFCNC_GenEvtCand;
    }//loop on events
    
    cout<<endl;
    
    //important: free memory
    treeLoader.UnLoadDataset();

  } //loop on datasets
    
  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;
  fout->cd();
  TDirectory* th1dir = fout->mkdir("Histos1D");
  th1dir->cd();
  for(map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
	TH1F *temp = it->second;
	temp->Write();
	//TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	//tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  TDirectory* th2dir = fout->mkdir("Histos2D");
  th2dir->cd();
   for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {
	TH2F *temp = it->second;
	temp->Write();
	//TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	//tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  //delete  
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

