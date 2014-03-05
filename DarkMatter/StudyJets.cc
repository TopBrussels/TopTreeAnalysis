#define _USE_MATH_DEFINES

// standard library
#include <iostream>
#include <map>
#include <cstdlib> 
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// root library
#include "TStyle.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TRandom.h"
#include "TMath.h"

// toptree library
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiCutPlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
//#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"

// custom cms plot style (tdr-like ?)
#include "TopTreeAnalysis/macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

int analysis(string, TString, int, float, float, float, int, float, float);

int main(int argc, char *argv[])
{

  string  xmlfile="configLite.xml";
  TString path="test";

  int cut_nJets=0;
  float cut_pt1=100;
  float cut_pt2=100;
  float cut_deltaphi=2;
  int puweight=0;

  float magnify=1.3;
  float magnifyLog=1.5;

  if(argc>1) xmlfile   = argv[1];
  if(argc>2) path      = argv[2];
  if(argc>3) cut_nJets = (int)atof(argv[3]);
  if(argc>4) cut_pt1   = atof(argv[4]);
  if(argc>5) cut_pt2   = atof(argv[5]);
  if(argc>6) cut_deltaphi = atof(argv[6]);
  if(argc>7) puweight  = (int)atof(argv[7]);

  analysis(xmlfile, path, cut_nJets, cut_pt1, cut_pt2, cut_deltaphi, puweight, magnify, magnifyLog);

  return 0;

}

int analysis(string xmlfile, TString path, 
	     int cut_nJets=2, float cut_pt1=100, float cut_pt2=100, float cut_deltaphi=2, int puweight=0, 
	     float magnify=1.3, float magnifyLog=1.5)
{

  /////////////////
  // ENVIRONMENT //
  /////////////////
  
  // Generic stuff
  clock_t start = clock();
  int verbose = 2;
  string pathPNG = path.Data();
  mkdir(pathPNG.c_str(),0777);

  // TLegend coordinates
  float x1=0.55; float y1=0.66; float x2=0.88; float y2=0.88;
  
  // Output ROOT file
  TString foutname = path+"/results.root";
  cout << "foutname=" << foutname << endl;
  TFile *fout = new TFile (foutname, "RECREATE");

  // Analysis environment
  AnalysisEnvironment anaEnv;
  cout << "- Loading environment : xmlfile=" << xmlfile << endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());
  cout << "-- env loaded !" << endl;

  // Load datasets
  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << "- Loading datasets" << endl;
  treeLoader.LoadDatasets(datasets, xmlfile.c_str());
  float Luminosity = 20000;

  // PU Reweighting
  LumiReWeighting LumiWeights = LumiReWeighting("PUReweighting/pileup_MC_S10.root", "PUReweighting/PURewtoptree_id_2014_1350565897_PileupHistogram.root", "pileup", "pileup");
  double lumiWeight=1;
  double lumiWeightOLD=1;

  // Apply JES
  int doJESShift = 0; // 0: off 1: minus 2: plus
  cout << "doJESShift: " << doJESShift << endl;

  // Global variable
  TRootEvent* event = 0;

  // Vectors of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* >   vertex;
  vector < TRootMuon* >     init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* >      init_jets;
  vector < TRootPFJet* >    init_jets_PF;
  vector < TRootMET* >      mets;
  vector<TRootJet*>         selectedJets;
  vector<TRootPFJet*>       selectedJetsPF;

  // jet variables
  const int nJet=2;
  Float_t tot_neutral_E_frac[nJet]={0,0};
  Float_t tot_charged_E_frac[nJet]={0,0};
  Float_t tot_chargedMult[nJet]={0,0};
  Float_t jet_E[  nJet]={0,0};
  Float_t jet_pt[ nJet]={0,0};
  Float_t jet_phi[nJet]={0,0};
  Float_t DeltaPhi=0;

  ///////////////////
  // DECLARE PLOTS //
  ///////////////////

  map<string,MultiSamplePlot*> MSPlot;

  // Event
  MSPlot["RhoCorrection"]         = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
  MSPlot["NbOfVertices"]          = new MultiSamplePlot(datasets, "NbOfVertices", 40, 0, 40, "Nb. of vertices");

  // Jet Numbers
  MSPlot["NbOfSelectedJets"]      = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedLightJets"] = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
  MSPlot["NbOfSelectedBJets"]     = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0, 8, "Nb. of jets");

  // Jet Variables
  string nameJet="";
  string titleJet[nJet]={"Leading Jet", "Sub-Leading Jet"};

  MSPlot["charged_E_frac_dijet"] = new MultiSamplePlot(datasets, "charged_E_frac_dijet", 100, 0, 1, "dijet charged energy fraction");
  MSPlot["chargedMult_dijet"   ] = new MultiSamplePlot(datasets, "chargedMult_dijet",    100, 0, 50,"dijet charged multiplicity");
  MSPlot["mass_dijet"          ] = new MultiSamplePlot(datasets, "mass_dijet",           6000, 0, 3000,"dijet mass (GeV)");

  for(int iJ=0 ; iJ<nJet ; iJ++) {
    nameJet=Form("_Jet%d",iJ+1);
    cout << "--- nameJet=" << nameJet << endl;
    
    MSPlot["Eta"+nameJet] = new MultiSamplePlot(datasets, "Eta"+nameJet, 30,  -5, 5,     titleJet[iJ]+" #eta");
    MSPlot["Phi"+nameJet] = new MultiSamplePlot(datasets, "Phi"+nameJet, 50,  -4, 4,     titleJet[iJ]+" #phi");
    MSPlot["Pt" +nameJet] = new MultiSamplePlot(datasets, "Pt" +nameJet, 600, 0,  3000,  titleJet[iJ]+" p_{T} (GeV)");
    MSPlot["E" +nameJet]  = new MultiSamplePlot(datasets, "E"  +nameJet, 600, 0,  3000,  titleJet[iJ]+" E (GeV)");

    MSPlot["chargedHad_E_frac"+nameJet] = new MultiSamplePlot(datasets, "chargedHad_E_frac"+nameJet, 100, 0, 1,   titleJet[iJ]+" charged had energy fraction");
    MSPlot["chargedEm_E_frac" +nameJet] = new MultiSamplePlot(datasets, "chargedEm_E_frac" +nameJet, 100, 0, 1,   titleJet[iJ]+" charged EM  energy fraction");
    MSPlot["chargedMu_E_frac" +nameJet] = new MultiSamplePlot(datasets, "chargedMu_E_frac" +nameJet, 100, 0, 1,   titleJet[iJ]+" charged Mu  energy fraction");
    MSPlot["neutralEm_E_frac" +nameJet] = new MultiSamplePlot(datasets, "neutralEm_E_frac" +nameJet, 100, 0, 1,   titleJet[iJ]+" neutral EM  energy fraction");
    MSPlot["neutralHad_E_frac"+nameJet] = new MultiSamplePlot(datasets, "neutralHad_E_frac"+nameJet, 100, 0, 1,   titleJet[iJ]+" neutral had energy fraction");    

    MSPlot["tot_charged_E_frac"+nameJet]= new MultiSamplePlot(datasets, "tot_charged_E_frac"+nameJet, 100, 0, 1,   titleJet[iJ]+" tot charged energy fraction");
    MSPlot["tot_neutral_E_frac"+nameJet]= new MultiSamplePlot(datasets, "tot_neutral_E_frac"+nameJet, 100, 0, 1,   titleJet[iJ]+" tot neutral energy fraction");

    MSPlot["chargedMult"      +nameJet] = new MultiSamplePlot(datasets, "chargedMult"      +nameJet, 100, 0, 50,  titleJet[iJ]+" charged multiplicity");
    MSPlot["neutralMult"      +nameJet] = new MultiSamplePlot(datasets, "neutralMult"      +nameJet, 100, 0, 50,  titleJet[iJ]+" neutral multiplicity");
    MSPlot["muonMult"         +nameJet] = new MultiSamplePlot(datasets, "muonMult"         +nameJet, 100, 0, 50,  titleJet[iJ]+" muon multiplicity");

    MSPlot["tot_chargedMult"  +nameJet] = new MultiSamplePlot(datasets, "tot_chargedMult"  +nameJet, 100, 0, 50,  titleJet[iJ]+" total charged multiplicity");
  }
  
  ///////////////////
  // Histograms
  ///////////////////

  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;
  
  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);

  ////////////////////////
  // LOOP OVER DATASETS //
  ////////////////////////

  u_int nDS=datasets.size();
  cout << "NUMBER OF DATASETS : " << nDS << endl;
  
  for (unsigned int d = 0; d < nDS; d++) {
    if (verbose > 1){
      cout << "   Dataset " << d << " name : " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;
      cout << " - Cross section = " << datasets[d]->Xsection() << endl;
      cout << " - IntLumi = " << datasets[d]->EquivalentLumi() << "  NormFactor = " << datasets[d]->NormFactor() << endl;
      cout << " - Nb of events : " << datasets[d]->NofEvtsToRunOver() << endl;
    }
    // open files and load
    cout << "- Load Dataset" << endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
    cout << "- loaded!" << endl;

    string previousFilename = "";
    int iFile = -1;
  
    string dataSetName = datasets[d]->Name();	
    
    /// Initialize JEC factors
    vector<JetCorrectorParameters> vCorrParam;
    /*    
    // Create the JetCorrectorParameter objects, the order does not matter.
    // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L1FastJet.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L2Relative.txt");
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L3Absolute.txt");
  
    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vCorrParam.push_back(*L1JetPar);
    vCorrParam.push_back(*L2JetPar);
    vCorrParam.push_back(*L3JetPar);
    
    if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") // Data!
    {
    JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("../macros/JECFiles/START42_V17_AK5PFchs_L2L3Residual.txt");
    vCorrParam.push_back(*ResJetCorPar);
    }
    */
    JetCorrectionUncertainty *jecUnc;// = new JetCorrectionUncertainty("../macros/JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false); // last boolean ('startFromRaw') = false!
    
    //////////////////////
    // LOOP OVER EVENTS //
    //////////////////////

    int itrigger = -1, previousRun = -1;
    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();
    cout <<"Number of events = "<<  end  <<endl;    

    bool debug = true;
    int event_start=0;
    double currentfrac =0.;
    double end_d = end;
    
    // Start loop over events
    for (unsigned int ievt = event_start; ievt <end_d ; ievt++) {
      
      double ievt_d = ievt;
      currentfrac = ievt_d/end_d;
      if(ievt%1000 == 0)
	std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      float rho = event->kt6PFJets_rho();
      MSPlot["RhoCorrection"]->Fill(rho, datasets[d], true, Luminosity);

      vector<TRootMCParticle*> mcParticles_flav;
      TRootGenEvent* genEvt_flav = 0;
      genEvt_flav = treeLoader.LoadGenEvent(ievt,false);
      treeLoader.LoadMCEvent(ievt, genEvt_flav, 0, mcParticles_flav,false); 

      vector<TRootGenJet*> genjets;
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
	  // loading GenJets as I need them for JER
	  genjets = treeLoader.LoadGenJet(ievt);
	}
      
      // check which file in the dataset it is to have the HLTInfo right
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename)
	{
	  previousFilename = currentFilename;
	  iFile++;
	  cout<<"File changed!!! => iFile = "<<iFile<<endl;
	}

      ///////////////////////////////////////
      //  Beam scraping and PU reweighting //
      ///////////////////////////////////////

      if(verbose>2) cout << "-- enter beam scraping and PU reweighting" << endl;

      // scale factor for the event
      float scaleFactor = 1.;
      
      if(dataSetName.find("Data") != 0 || dataSetName.find("data") != 0 || dataSetName.find("DATA") != 0) {
	// Apply the scraping veto. (Is it still needed?)
	bool isBeamBG = true;
	if(event->nTracks() > 10) {
	  if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
	    isBeamBG = false;
	}
	if(isBeamBG) continue;
      }
      else if(puweight!=0) {
	lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
	scaleFactor = scaleFactor*lumiWeight;	
      }
      histo1D["lumiWeights"]->Fill(scaleFactor);	

      ///////////////////////////////////////////////////////////
      //   Event selection
      ///////////////////////////////////////////////////////////

      if(verbose>2) cout << "-- enter the event selection" << endl;

      // Trigger information
      /*
	bool trigged = false;
	//std::string filterName = "";
	int currentRun = event->runId();
	if(previousRun != currentRun) {
	previousRun = currentRun;
	itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);
	//filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";
	}				  
	trigged = treeLoader.EventTrigged (itrigger);
	if(verbose>2) cout << "triggered? Y/N?  " << trigged  << endl;
	// if(!trigged) continue;
	*/

      //////////////////////
      // OBJECT SELECTION //
      //////////////////////

      // Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets);

      // CaloJets selection (we should implement a PFJet selection in the Selection class)
      // Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets
      //selection.setJetCuts(40.,2.5,0.01,1.,0.98,0,0); // Jet ID
      selection.setJetCuts(0.,5.,0.,0.,0.,0.,0.);//ND : No Jet ID

      // get jet collection 
      selectedJets   = selection.GetSelectedJets(false); // false : do not apply hardcoded ID
      selectedJetsPF = jetTools->convertToPFJets(selectedJets); // convert TRootJet to TRootPFJet
      sort(selectedJetsPF.begin(),selectedJetsPF.end(),HighestPt()); // order in pT

      /////////////////////

      /////////////////////
      // EVENT SELECTION //
      /////////////////////

      // Primary vertex
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
      if(!isGoodPV) continue;

      // JETS

      // number of jets
      if(selectedJets.size()<cut_nJets) continue;
      if(verbose>2) cout << "-- selectedJets.size()=" << selectedJets.size() << endl;

      // leading jet pt
      jet_pt[0]=selectedJets[0]->Pt();
      if(jet_pt[0]<cut_pt1) continue;

      // subleading jet pt
      jet_pt[1]=selectedJets[1]->Pt();
      if(jet_pt[1]<cut_pt2) continue;

      // delta phi
      jet_phi[0]=selectedJets[0]->Phi();
      jet_phi[1]=selectedJets[1]->Phi();
      DeltaPhi = TMath::Abs(selectedJets[0]->Phi() - selectedJets[1]->Phi());
      if(DeltaPhi<cut_deltaphi) continue;

      /////////////////////

      //////////////////////
      /// FILL HISTOGRAMS //
      //////////////////////

      // EVENT
      MSPlot["NbOfVertices"]         ->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

      // JETS
      MSPlot["NbOfSelectedJets"]     ->Fill(selectedJetsPF.size(), datasets[d], true, Luminosity*scaleFactor);

      for(int iJ=0 ; iJ<nJet ; iJ++) { // loop over jets

	if(selectedJetsPF.size()<=iJ) break;
	
	nameJet=Form("_Jet%d",iJ+1);

	jet_E[iJ]=selectedJetsPF[iJ]->E();

	tot_neutral_E_frac[iJ] = selectedJetsPF[iJ]->neutralHadronEnergyFraction() 
	  + selectedJetsPF[iJ]->neutralEmEnergyFraction();

	tot_charged_E_frac[iJ] = selectedJetsPF[iJ]->chargedHadronEnergyFraction() 
	  + selectedJetsPF[iJ]->chargedEmEnergyFraction() 
	  + selectedJetsPF[iJ]->chargedMuEnergyFraction();

	tot_chargedMult[iJ] = selectedJetsPF[iJ]->chargedMultiplicity() + selectedJetsPF[iJ]->muonMultiplicity();

	MSPlot["Eta"+nameJet]->Fill(selectedJetsPF[iJ]->Eta(), datasets[d], true, Luminosity*scaleFactor);

	if(iJ>1) {
	  MSPlot["Phi"+nameJet]->Fill(selectedJetsPF[iJ]->Phi(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["Pt"+nameJet] ->Fill(selectedJetsPF[iJ]->Pt(),  datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["E"+nameJet]  ->Fill(selectedJetsPF[iJ]->E(),   datasets[d], true, Luminosity*scaleFactor);
	}
	else {
	  MSPlot["Phi"+nameJet]->Fill(jet_phi[iJ], datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["Pt"+nameJet] ->Fill(jet_pt[iJ],  datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["E"+nameJet]  ->Fill(jet_E[iJ],   datasets[d], true, Luminosity*scaleFactor);
	}

	MSPlot["chargedHad_E_frac"+nameJet]->Fill(selectedJetsPF[iJ]->chargedHadronEnergyFraction(),  datasets[d], true, Luminosity*scaleFactor);
	MSPlot["chargedEm_E_frac" +nameJet]->Fill(selectedJetsPF[iJ]->chargedEmEnergyFraction(),      datasets[d], true, Luminosity*scaleFactor);
	MSPlot["chargedMu_E_frac" +nameJet]->Fill(selectedJetsPF[iJ]->chargedMuEnergyFraction(),      datasets[d], true, Luminosity*scaleFactor);
	MSPlot["neutralEm_E_frac" +nameJet]->Fill(selectedJetsPF[iJ]->neutralEmEnergyFraction(),      datasets[d], true, Luminosity*scaleFactor);
	MSPlot["neutralHad_E_frac"+nameJet]->Fill(selectedJetsPF[iJ]->neutralHadronEnergyFraction(),  datasets[d], true, Luminosity*scaleFactor);

	MSPlot["tot_neutral_E_frac"+nameJet]->Fill(tot_neutral_E_frac[iJ],  datasets[d], true, Luminosity*scaleFactor);
	MSPlot["tot_charged_E_frac"+nameJet]->Fill(tot_charged_E_frac[iJ],  datasets[d], true, Luminosity*scaleFactor);

	MSPlot["chargedMult"      +nameJet]->Fill(selectedJetsPF[iJ]->chargedMultiplicity(),          datasets[d], true, Luminosity*scaleFactor);
	MSPlot["neutralMult"      +nameJet]->Fill(selectedJetsPF[iJ]->neutralMultiplicity(),          datasets[d], true, Luminosity*scaleFactor);
	MSPlot["muonMult"         +nameJet]->Fill(selectedJetsPF[iJ]->muonMultiplicity(),             datasets[d], true, Luminosity*scaleFactor);
	MSPlot["tot_chargedMult"  +nameJet]->Fill(tot_chargedMult[iJ], datasets[d], true, Luminosity*scaleFactor);
      }

      Float_t charged_E_frac_dijet=0, chargedMult_dijet=0, mass_dijet=0;
      TRootPFJet dijet;
      if(selectedJetsPF.size()>=2) {
	charged_E_frac_dijet = (jet_E[0]+jet_E[1]!=0) ? (jet_E[0]*tot_charged_E_frac[0]+jet_E[1]*tot_charged_E_frac[1]) / (jet_E[0]+jet_E[1]) : -999;
	chargedMult_dijet = tot_chargedMult[0] + tot_chargedMult[1];
	dijet = (*selectedJetsPF[0]) + (*selectedJetsPF[1]) ;
	mass_dijet = dijet.M();
      }

      MSPlot["charged_E_frac_dijet"]->Fill(charged_E_frac_dijet, datasets[d], true, Luminosity*scaleFactor);
      MSPlot["chargedMult_dijet"   ]->Fill(chargedMult_dijet,    datasets[d], true, Luminosity*scaleFactor);
      MSPlot["mass_dijet"          ]->Fill(mass_dijet,           datasets[d], true, Luminosity*scaleFactor);


    } // end loop over events

    // free memory
    //if(jetTools) delete jetTools;
    treeLoader.UnLoadDataset();

  } // end loop over datasets
  
  //////////////////
  // WRITE OUTPUT //
  //////////////////

  cout << "- Start writing output" << endl;
  fout->cd();

  cout << "-- Loop over MSPlots" << endl;
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++) {
    string name = it->first;
    cout << "--- name=" << name << endl;
    
    MultiSamplePlot *temp = it->second;
    if(!temp) {
      cout <<"--- error: no such MSPlot!" << endl;
      continue;
    }

    TH1F *tempHisto_data;
    TH1F *tempHisto_TTTT;

    //Option to write ROOT files containing histograms with systematics varied +/-
    //TFile * tempErrorFile;
        
    if (doJESShift == 1){
      string filename = "ErrorBands/Error_"+name+".root";
      TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
      //TH1F * tempHisto = temp->getTH1F("TTJets");
      //tempHisto->Write("Minus");
      // tempErrorFile->Write();
      //tempErrorFile->Close();
    }
    else if  (doJESShift ==2){
            
      string filename = "ErrorBands/Error_"+name+".root";
      TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
      //TH1F * tempHisto = temp->getTH1F("TTJets");
      //tempHisto->Write("Plus");
      //tempErrorFile->Write();
      //tempErrorFile->Close();

      if(verbose>2) cout <<"--- JES sys down"<<endl;
            
    }
    else if  (doJESShift ==0){
      if(verbose>2) cout <<"--- JES sys off "<<endl;

      // MultiSamplePlot::Draw(std::string, unsigned int, bool, bool, bool, int)
      // void MultiSamplePlot::Draw(string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal)

      // MultiSamplePlot::Write(TFile*, std::string, bool, std::string, std::string)
      // void MultiSamplePlot::Write(TFile* fout, string label, bool savePNG, string pathPNG, string ext) 

      // temp->addText("CMS preliminary");

      // ND
      if(verbose>2) cout << "--- draw" << endl;
      temp->Draw(name, 0, false, false, false, 0, x1, y1, x2, y2, magnify);
      //temp->Draw(name, 0, false, false, false, 0);
      if(verbose>2) cout << "--- drawn!" << endl;

      if(verbose>2) cout << "--- write in pathPNG=" << pathPNG << endl;
      float maxY = temp->getMaxY();
      cout << "###############################" << endl 
	   << "########## maxY=" << maxY << " ##########" << endl 
	   << "###############################" << endl;
      //temp->setMaxY(1000000*maxY);

      temp->Write(fout, name, true, pathPNG+"/", "png", magnifyLog); //ND true => SaveAs the Canvas as image => seg fault probably caused by empty plots !
      //temp->Write(fout, name, false, pathPNG, "png");
      if(verbose>2) cout << "--- written!" << endl;
    }
  } // end loop over MSPlots
  cout << "-- end loop over MSPlots" << endl;
  
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
    
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << "s to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;


  return 0;
}
