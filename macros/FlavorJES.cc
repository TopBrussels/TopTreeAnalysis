#include "TStyle.h"
#include "TF2.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysis/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysis/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysis/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysis/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysis/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysis/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysis/Content/interface/Dataset.h"
#include "TopTreeAnalysis/Tools/interface/JetTools.h"
#include "TopTreeAnalysis/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysis/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysis/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysis/MCInformation/interface/JetPartonMatching.h"
#include "Style.C"

using namespace std;
using namespace TopTree;

/// TGraphAsymmErrors
map<string,TGraphAsymmErrors*> graphAsymmErr;
map<string,TGraphErrors*> graphErr;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

int main (int argc, char *argv[])
{
  bool doPF2PAT = false;

  clock_t start = clock();

  cout << "*****************************************************" << endl;
  cout << " Beginning of the program for Flavor dependent JES ! " << endl;
  cout << "*****************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////
  
  //xml file
  string xmlFileName = "../config/myFlavorJESconfig.xml";
  if (argc >= 2)
    xmlFileName=string(argv[1]);
  const char *xmlfile = xmlFileName.c_str();
  cout << "Used config file:  " << xmlfile << endl;
  
  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<"Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = anaEnv.Verbose;
  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
 
 	cout << "analysis environment luminosity for rescaling "<< oldLuminosity << endl;
 	
  /////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  treeLoader.LoadDatasets (datasets, xmlfile);
  
  float Luminosity = oldLuminosity;
	for (unsigned int d = 0; d < datasets.size (); d++)
  {
		//cout << "luminosity of dataset "<< d << " is " << datasets[d]->EquivalentLumi() << endl;
		if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
    string dataSetName = datasets[d]->Name();
		if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
		{
		  Luminosity = datasets[d]->EquivalentLumi();
		  break;
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

  TFile *fout = new TFile ("FlavorJES.root", "RECREATE");
  //Global variable
  TRootEvent* event = 0;
  
  string pathPNG = "PlotsFlavorJES/"; 	
  mkdir(pathPNG.c_str(),0777);
	
	//nof selected events
  float NEvtsData = 0;
  
  TH2F* KinFitPNegDeriv = 0;
  
  // Define the plots
  MSPlot["nSelectedJets"] = new MultiSamplePlot(datasets, "nSelectedJets", 20, -0.5, 19.5, "Nr. of Selected Jets");
  MSPlot["JetPt"] = new MultiSamplePlot(datasets, "JetPt", 100, 0, 2000, "P_{T} of the Selected Jets");
  
  histo1D["JetPt"] = new TH1F("JetPt","JetPt",100,0,100);
  
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTable;
  CutsSelecTable.push_back(string("Initial"));
  CutsSelecTable.push_back(string("trigged"));
  CutsSelecTable.push_back(string("Good PV"));
  CutsSelecTable.push_back(string("1 good jet"));
  
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  
  /////////////////////////////
  /// ResolutionFit Stuff
  /////////////////////////////
  
  bool CalculateResolutions = true; // If false, the resolutions will be loaded from a previous calculation
  bool ResolutionsClosure = false;
  
  ResolutionFit *resFitLightJets = 0, *resFitBJets = 0, *resFitLightJetsL7 = 0, *resFitBJetsL7 = 0, *resFitBJets_B = 0, *resFitBJets_Bbar = 0;
  if(CalculateResolutions)
  {
    resFitLightJets = new ResolutionFit("LightJet");
    resFitBJets = new ResolutionFit("BJet");
    resFitBJets_B = new ResolutionFit("BJet_B");
    resFitBJets_Bbar = new ResolutionFit("BJet_Bbar");
    if(ResolutionsClosure)
    {
      resFitLightJets->LoadResolutions("resolutions/lightJetReso.root");
      resFitBJets->LoadResolutions("resolutions/bJetReso.root");
      resFitBJets_B->LoadResolutions("resolutions/bJetReso.root");
      resFitBJets_Bbar->LoadResolutions("resolutions/bJetReso.root");
    }
  }
  else
  {
    resFitLightJets = new ResolutionFit("LightJet");
    resFitLightJets->LoadResolutions("resolutions/lightJetReso.root");
    resFitLightJetsL7 = new ResolutionFit("LightJetL7");
    resFitLightJetsL7->LoadResolutions("resolutions/lightJetReso_AfterL7.root");
    resFitBJets = new ResolutionFit("BJet");
    resFitBJets->LoadResolutions("resolutions/bJetReso.root");
    resFitBJetsL7 = new ResolutionFit("BJetL7");
    resFitBJetsL7->LoadResolutions("resolutions/bJetReso_AfterL7.root");
  }
  if (verbose > 0)
    cout << " - ResolutionFit instantiated ..." << endl;
  
  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;
    if (verbose > 1)
      std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
		
		//open files and load
    cout<<"LoadEvent"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout<<"LoadEvent"<<endl;
    
    string dataSetName = datasets[d]->Name();
    string previousFilename = "";
    int iFile = -1;
    
    /////////////////////////////////////
   	/// Initialize JEC factors
   	/////////////////////////////////////
   	
    vector<JetCorrectorParameters> vCorrParam;
//    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) // Data!
//    {
//      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Jec11V2_db_AK5PFchs_L2L3Residual.txt");
//      vCorrParam.push_back(*ResJetCorPar);
//    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/Jec11V2_db_AK5PFchs_Uncertainty.txt");
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false);
    
    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    int itrigger = -1, itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;

//    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    for (unsigned int ievt = 0; ievt < 5000; ievt++)
    {
      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
      
      //load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      vector<TRootGenJet*> genjets;
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
      {
        genjets = treeLoader.LoadGenJet(ievt, false);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
//      cout << "run: " << event->runId() << "  lumi: " << event->lumiBlockId() << "  event: " << event->eventId() << endl;
      
      // scale factors for the event
      float scaleFactor = 1.;
      float lumiWeight = 1.;
      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
        lumiWeight = 1;
      
      // check which file in the dataset it is to have the HLTInfo right
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename)
      {
      	previousFilename = currentFilename;
	      iFile++;
	      cout<<"File changed!!! => iFile = "<<iFile<<endl;
      }
      
      int currentRun = event->runId();
      if(previousRun != currentRun)
      {
        previousRun = currentRun;
        
        // semi-muon
        if(dataSetName.find("Data_MuHad") == 0 || dataSetName.find("data_MuHad") == 0 || dataSetName.find("DATA_MuHad") == 0 )
        {
          if( event->runId() <= 161176 )
            itrigger = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v1"), currentRun, iFile);
        }
        else
        {
          itrigger = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v2"), currentRun);
      	  if (itrigger == 9999)
      	    itrigger = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v3"), currentRun); // Fall11 44X Chamonix
        }
      }
      // Apply Jet Corrections on-the-fly
//      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
//        jetTools->correctJets(init_jets, vertex);
      
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
      {
        jetTools->correctJetJER(init_jets, genjets, "nominal");
        
        //Scale jets with a certain factor
//        jetTools->scaleJets(init_jets, 1.);
      }
      /////////////////////////////
      //   Selection
      /////////////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets);
      selection.setJetCuts(30.,2.5,0.01,1.,0.98,0.3,0.1);
      
      bool trigged = false;
//      trigged = treeLoader.EventTrigged (itrigger);
      trigged = true;
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);
      
      vector<TRootJet*> selectedJets = selection.GetSelectedJets(true);
      
      bool eventSelected = false;
      
      selecTable.Fill(d,0,scaleFactor*lumiWeight);
      if( trigged )
      {
        selecTable.Fill(d,1,scaleFactor*lumiWeight);
        if( isGoodPV )
        {
          selecTable.Fill(d,2,scaleFactor*lumiWeight);
          if( selectedJets.size() >= 1 )
          {
            eventSelected = true;
          }
        }
      }
      
      if( ! eventSelected ) continue;
      
      vector<TRootMCParticle*> mcParticles;
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
      {
        mcParticles = treeLoader.LoadMCPart(ievt, false);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      // Now do something with those events!
      MSPlot["nSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      
      vector<TLorentzVector> mcPartTLV, jetsTLV, genJetsTLV;
      for(size_t iJet=0; iJet<selectedJets.size(); iJet++)
      {
        TRootJet* jet = selectedJets[iJet];
        histo1D["JetPt"]->Fill(jet->Pt());
        MSPlot["JetPt"]->Fill(jet->Pt(), datasets[d], true, Luminosity*scaleFactor);
        jetsTLV.push_back(*jet);
      }
      
      for(size_t iGenJet=0; iGenJet<genjets.size(); iGenJet++)
        genJetsTLV.push_back(*genjets[iGenJet]);
      
      for(size_t iPart=0; iPart<mcParticles.size(); iPart++)
      {
        TRootMCParticle* part = mcParticles[iPart];
        if( part->status() != 3 ) continue;
        if( part->motherType() == 2212 ) continue;
//        cout << part->Eta() << " " << part->Pt() << " | " << part->status() << " | " << part->type() << " | granny: " << part->grannyType() << " mother: " << part->motherType() << " | nDau: " << part->nDau() << "  " << part->dauOneId() << " " << part->dauTwoId() << " " << part->dauThreeId() << " " << part->dauFourId() << endl;
        mcPartTLV.push_back(*part);
      }
      
      JetPartonMatching matchingJetParton = JetPartonMatching(mcPartTLV, jetsTLV, 2, true, true, 0.3);
      JetPartonMatching matchingJetGenJet = JetPartonMatching(genJetsTLV, jetsTLV, 2, true, true, 0.3);
      JetPartonMatching matchingGenJetParton = JetPartonMatching(mcPartTLV, genJetsTLV, 2, true, true, 0.3);
      
      vector< pair<unsigned int, unsigned int> > JetPartonPair, JetGenJetPair, GenJetPartonPair; // First one is jet number, second one is mcParticle number
      for(unsigned int i=0; i<mcPartTLV.size(); i++)
      {
        int matchedJetNumber = matchingJetParton.getMatchForParton(i, 0);
        if(matchedJetNumber != -1)
          JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
        matchedJetNumber = matchingGenJetParton.getMatchForParton(i, 0);
        if(matchedJetNumber != -1)
          GenJetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
      }
      for(unsigned int i=0; i<genJetsTLV.size(); i++)
      {
        int matchedJetNumber = matchingJetGenJet.getMatchForParton(i, 0);
        if(matchedJetNumber != -1)
          JetGenJetPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
      }

      
      
      
      
      
      
      
      
    }				//loop on events
    cout<<endl;

    if(jetTools) delete jetTools;
    
    treeLoader.UnLoadDataset(); //important: free memory
  }				//loop on datasets
  
  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;
  
  //Selection tables
  selecTable.TableCalculator(false, true, true, true, true);
	selecTable.Write("SelectionTable_FlavorJES.tex");
  
  // Do some special things with certain plots (normalize, BayesDivide, ... )
  if (verbose > 0)
    cout << "Treating the special plots." << endl;
  
  ///////////////////
  // Writing
  //////////////////
  if (verbose > 1)
  	cout << " - Writing outputs to the files ..." << endl;
  // Fill the resolution histograms and calculate the resolutions
  if(CalculateResolutions)
  {
    mkdir((pathPNG+"resFit_LightJet/").c_str(),0777);
    mkdir((pathPNG+"resFit_BJet/").c_str(),0777);
    mkdir((pathPNG+"resFit_BJet_B/").c_str(),0777);
    mkdir((pathPNG+"resFit_BJet_Bbar/").c_str(),0777);
  
//    resFitLightJets->WritePlots(fout, true, pathPNG+"resFit_LightJet/");
//    resFitLightJets->WriteResolutions("lightJetReso.root");
//    resFitBJets->WritePlots(fout, true, pathPNG+"resFit_BJet/");
//    resFitBJets->WriteResolutions("bJetReso.root");
//    resFitBJets_B->WritePlots(fout, true, pathPNG+"resFit_BJet_B/");
//    resFitBJets_B->WriteResolutions("bJetReso_B.root");
//    resFitBJets_Bbar->WritePlots(fout, true, pathPNG+"resFit_BJet_Bbar/");
//    resFitBJets_Bbar->WriteResolutions("bJetReso_Bbar.root");
  }
  
  cout << "Writing the histograms..." << endl;
  // 1D 
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  mkdir((pathPNG+"MSPlot/").c_str(),0777);
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(false, name, true, false, true, true, true);
    temp->Write(fout, name, true, pathPNG+"MSPlot/");
  }

  //Write histograms
  fout->cd();
  th1dir->cd();

  fout->cd();

	for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
	{
		TH1F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}

  // 2D
  TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
  th2dir->cd();
	for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
	{
		TH2F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}

	//Write TGraphAsymmErrors
	for(map<string,TGraphAsymmErrors*>::const_iterator it = graphAsymmErr.begin(); it != graphAsymmErr.end(); it++)
	{
	  TGraphAsymmErrors *temp = it->second;
	  temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
	}
	
  fout->cd();
  //add configuration info
  fout->cd();

  //Write TGraphErrors
  fout->cd();
  for(map<string,TGraphErrors*>::const_iterator it = graphErr.begin(); it != graphErr.end(); it++)
  {
    TGraphErrors *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }

  //
  if (verbose > 1)
    cout << " - Done with writing the module outputs in the ouput file ..." << endl;
  cout << " - Closing the output file now..." << endl;
//  fout->Write();
  fout->Close();

  //delete
//  if(resFitLightJets) delete resFitLightJets;
//  if(resFitBJets) delete resFitBJets;
//  if(resFitBJets_B) delete resFitBJets_B;
//  if(resFitBJets_Bbar) delete resFitBJets_Bbar;
  
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;

  return 0;
}
