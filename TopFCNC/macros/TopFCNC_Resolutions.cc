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
  cout << " Beginning of the program for the FCNC resolution study ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
  //setMyStyle();


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string xmlFileName = "../config/myTopFCNCconfig_Resolutions.xml";
  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    

  string postfix = "_Resolutions"; // to relabel the names of the output file

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
  TopFCNC_GenEvt *GenEvent = 0;

  ResolutionFit *resFitMuons     = new ResolutionFit("Muon");
  ResolutionFit *resFitElectrons = new ResolutionFit("Electron");
  ResolutionFit *resFitBJets     = new ResolutionFit("BJet");
  ResolutionFit *resFitQJets     = new ResolutionFit("BJet");
  ResolutionFit *resFitLightJets = new ResolutionFit("LightJet");
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
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  histo1D["DR_muon_MC_vs_RECO"] = new TH1F("DR_muon_MC_vs_RECO",";#Delta R;#events",100,0,0.005);
  histo1D["DR_elec_MC_vs_RECO"] = new TH1F("DR_elec_MC_vs_RECO",";#Delta R;#events",100,0,0.005);
  histo1D["DR_qjets_MC_vs_RECO"] = new TH1F("DR_qjets_MC_vs_RECO",";#Delta R;#events",100,0,0.005);
  histo1D["DR_ljets_MC_vs_RECO"] = new TH1F("DR_ljets_MC_vs_RECO",";#Delta R;#events",100,0,0.005);

  histo1D["DiMuonMass"] = new TH1F("DiMuonMass",";m_{ll};#events",320,50,130);
  histo1D["DiElecMass"] = new TH1F("DiElecMass",";m_{ll};#events",320,50,130);

  histo1D["WMass_Had_Decay"]    = new TH1F("WMass_Had_Decay",";m_{qq};#events",480,0,160);
  histo1D["TopMass_Had_Decay"]  = new TH1F("TopMass_Had_Decay",";m_{bqq};#events",500,100,250);
  histo1D["TopMass_Fcnc_Decay"] = new TH1F("TopMass_Fcnc_Decay",";m_{llq};#events",500,100,250);

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
    //unsigned int end = 10000;
     
    if (verbose > 1) cout << " - Loop over events " << endl;      
    
    for (unsigned int ievt = start; ievt < end; ievt++)
    {        

      if(ievt%1000 == 0)
		  std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

	    // load Event
	    event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      // load MCParticles
	    vector<TRootMCParticle*> mcParticles = treeLoader.LoadMCPart(ievt);

	    Selection selection(init_jets, init_muons, init_electrons, mets); //mets can also be corrected...
	    selection.setJetCuts(20.,2.5,0.01,1.,0.98,0.3,0.1); //Pt,Eta,EMF,n90Hits,fHPD,dRJetElectron,dRJetMuon
    	selection.setDiMuonCuts(15,2.4,0.20,999.); //Et,Eta,RelIso,d0
  	  selection.setDiElectronCuts(15,2.5,0.15,0.04,0.,1,0.3,1); //Et,Eta,RelIso,d0,MVAId,DistVzPVz,DRJets,MaxMissHits

	    vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId
	    vector<TRootMuon*>     selectedMuons       = selection.GetSelectedDiMuons();
	    vector<TRootElectron*> selectedElectrons   = selection.GetSelectedDiElectrons();

      // instanciate topfcnc object
      GenEvent = new TopFCNC_GenEvt();
      // reconstruct top fcnc event
      GenEvent->ReconstructEvt(mcParticles);
      // match jets to top fcnc partons
      GenEvent->MatchJetsToPartons(selectedJets);

      if(GenEvent->zLeptonicChannel() == TopFCNC_GenEvt::kMuon){
        // match muons to top fcnc Z decay products
        GenEvent->MatchLeptonsToZ(selectedMuons);
        // fill corresponding histograms
        histo1D["DR_muon_MC_vs_RECO"]->Fill(GenEvent->DR_MatchLepton1FromZ());
        histo1D["DR_muon_MC_vs_RECO"]->Fill(GenEvent->DR_MatchLepton2FromZ());
        histo1D["DiMuonMass"]->Fill((GenEvent->matchedLepton1FromZ()+GenEvent->matchedLepton2FromZ()).M());
      }
      else if(GenEvent->zLeptonicChannel() == TopFCNC_GenEvt::kElec){
        // match electrons to top fcnc Z decay products
        GenEvent->MatchLeptonsToZ(selectedElectrons);
        // fill corresponding histograms
        histo1D["DR_elec_MC_vs_RECO"]->Fill(GenEvent->DR_MatchLepton1FromZ());
        histo1D["DR_elec_MC_vs_RECO"]->Fill(GenEvent->DR_MatchLepton2FromZ());
        histo1D["DiElecMass"]->Fill((GenEvent->matchedLepton1FromZ()+GenEvent->matchedLepton2FromZ()).M());
      }

      histo1D["DR_bjets_MC_vs_RECO"]->Fill(GenEvent->DR_MatchBFromTop());
      histo1D["DR_qjets_MC_vs_RECO"]->Fill(GenEvent->DR_MatchQFromTop());

      histo1D["DR_ljets_MC_vs_RECO"]->Fill(GenEvent->DR_MatchQuark1FromW());
      histo1D["DR_ljets_MC_vs_RECO"]->Fill(GenEvent->DR_MatchQuark2FromW());

      histo1D["WMass_Had_Decay"]   ->Fill((GenEvent->matchedQuark1FromW()+GenEvent->matchedQuark2FromW()).M());
      histo1D["TopMass_Had_Decay"] ->Fill((GenEvent->matchedQuark1FromW()+GenEvent->matchedQuark2FromW()+GenEvent->matchedB()).M());
      histo1D["TopMass_Fcnc_Decay"]->Fill((GenEvent->matchedLepton1FromZ()+GenEvent->matchedLepton2FromZ()+GenEvent->matchedQ()).M());

      GenEvent->FillResolutions(resFitMuons, resFitElectrons, resFitBJets, resFitQJets, resFitLightJets);

	    //delete selection;
	    delete GenEvent;

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

  resFitMuons->WritePlots(fout, true, pathPNG+"resFit_Muon/");
  resFitMuons->WriteResolutions("muonReso.root");
  resFitElectrons->WritePlots(fout, true, pathPNG+"resFit_Electron/");
  resFitElectrons->WriteResolutions("electronReso.root");

  resFitBJets->WritePlots(fout, true, pathPNG+"resFit_BJet/");
  resFitBJets->WriteResolutions("bJetReso.root");
  resFitQJets->WritePlots(fout, true, pathPNG+"resFit_QJet/");
  resFitQJets->WriteResolutions("qJetReso.root");
  resFitLightJets->WritePlots(fout, true, pathPNG+"resFit_LightJet/");
  resFitLightJets->WriteResolutions("lightJetReso.root");

  //delete  
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

