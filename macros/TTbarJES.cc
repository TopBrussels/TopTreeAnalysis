///////////////////////////
///// TODO & COMMENTS /////
/////////////////////////// 

#include "TStyle.h"
#include "TF2.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../Selection/interface/SelectionTable.h"
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../Tools/interface/MVATrainer.h"
#include "../Tools/interface/MVAComputer.h"
#include "../Tools/interface/JetTools.h"
#include "../MCInformation/interface/ResolutionFit.h"
#include "../JESMeasurement/interface/FullKinFit.h"
#include "../JESMeasurement/interface/JetCombiner.h"
#include "../Reconstruction/interface/JetCorrectorParameters.h"
#include "../Reconstruction/interface/JetCorrectionUncertainty.h"
#include "../JESMeasurement/interface/LightMonster.h"
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

  cout << "******************************************" << endl;
  cout << " Beginning of the program for TTbar JES ! " << endl;
  cout << "******************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

  //which systematic to run?
  string systematic = "Nominal";
  if (argc >= 2)
		systematic = string(argv[1]);
  cout << "Systematic to be used:  " << systematic << endl;
  if( ! (systematic == "Nominal" || systematic == "InvertedIso" || systematic == "JESPlus" || systematic == "JESMinus" || systematic == "JERPlus" || systematic == "JERMinus" || systematic == "AlignPlus" || systematic == "AlignMinus") )
  {
    cout << "Unknown systematic!!!" << endl;
    cout << "Possible options are: Nominal, InvertedIso, JESPlus, JESMinus, JERPlus, JERMinus" << endl;
    exit(-1);
  }
  
  //xml file
  string xmlFileName = "../config/myJESconfig.xml";
  if (argc >= 3)
    xmlFileName=string(argv[2]);
  const char *xmlfile = xmlFileName.c_str();
  cout << "Used config file:  " << xmlfile << endl;
  
  //Output ROOT file
  string rootFileName ("TTbarJES.root");
 	const char *rootfile = rootFileName.c_str();  
  
  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  
  /////////////////////////
  // Which decay channel //
  /////////////////////////
  
  bool semiElectron = true; // use semiElectron channel?
  bool semiMuon = true; // use semiMuon channel?
  if(semiElectron && semiMuon) cout << "  --> Using semiMuon and semiElectron channel..." << endl;
  else
  {
    if(semiMuon) cout << " --> Using the semiMuon channel..." << endl;
    else cout << " --> Using the semiElectron channel..." << endl;
  }
  
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
  for(unsigned int i=0;i<datasets.size();i++) // Rename datasets, taking systematic name into account
  {
    string dataSetName = datasets[i]->Name();
    datasets[i]->SetName( dataSetName+"_"+systematic );
  }
  
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
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
  vector < TRootJet* > init_jets_corrected;
  vector < TRootMET* > mets;

  TFile *fout = new TFile (rootfile, "RECREATE");
  //Global variable
  TRootEvent* event = 0;
  
  string pathPNG = "PlotsJES/"; 	
  mkdir(pathPNG.c_str(),0777);
	
	//nof selected events
  float NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
  TH2F* KinFitPNegDeriv = 0;
  
  // Define the plots
  MSPlot["AllMuonsRelPFIso"] = new MultiSamplePlot(datasets, "AllMuonsRelPFIso", 100, 0, 5, "Muons PF Iso");
  MSPlot["AllElectronsRelPFIso"] = new MultiSamplePlot(datasets, "AllElectronsRelPFIso", 100, 0, 5, "Electrons PF Iso");
  MSPlot["SelectedEventsMuonsRelPFIso"] = new MultiSamplePlot(datasets, "SelectedEventsMuonsRelPFIso", 100, 0, 5, "Muons PF Iso");
  MSPlot["SelectedEventsElectronsRelPFIso"] = new MultiSamplePlot(datasets, "SelectedEventsElectronsRelPFIso", 100, 0, 5, "Electrons PF Iso");
  MSPlot["nEventsAfterCutsSemiMu"] = new MultiSamplePlot(datasets, "nEventsAfterCutsSemiMu", 23, -0.5, 22.5, "Nr. of events after each cut, SemiMu");
  MSPlot["nEventsAfterCutsSemiEl"] = new MultiSamplePlot(datasets, "nEventsAfterCutsSemiEl", 23, -0.5, 22.5, "Nr. of events after each cut, SemiEl");
	MSPlot["nPrimaryVertices"] = new MultiSamplePlot(datasets, "nPrimaryVertices", 12, -0.5, 11.5, "Nr. of primary vertices");
	MSPlot["nGoodPrimaryVertices"] = new MultiSamplePlot(datasets, "nGoodPrimaryVertices", 12, -0.5, 11.5, "Nr. of good primary vertices");
  MSPlot["nSelectedMuons"] = new MultiSamplePlot(datasets, "nSelectedMuons", 5, -0.5, 4.5, "Nr. of selected Muons");
  MSPlot["nSelectedElectrons"] = new MultiSamplePlot(datasets, "nSelectedElectrons", 5, -0.5, 4.5, "Nr. of selected Electrons");
  MSPlot["nLooseOtherMuons"] = new MultiSamplePlot(datasets, "nLooseOtherMuons", 5, -0.5, 4.5, "Nr. of Loose Other Muons");
  MSPlot["nLooseElectronsSemiMu"] = new MultiSamplePlot(datasets, "nLooseElectronsSemiMu", 5, -0.5, 4.5, "Nr. of Loose Electrons, SemiMu");
  MSPlot["nLooseElectronsSemiEl"] = new MultiSamplePlot(datasets, "nLooseElectronsSemiEl", 5, -0.5, 4.5, "Nr. of Loose Electrons, SemiEl");
  MSPlot["nSelectedJets"] = new MultiSamplePlot(datasets, "nSelectedJets", 10, -0.5, 9.5, "Nr. of Selected Jets");
  
  histo1D["FourthJetPt"] = new TH1F("FourthJetPt","FourthJetPt",100,0,100);
  histo1D["FourthJetPtTriggered"] = new TH1F("FourthJetPtTriggered","FourthJetPtTriggered",100,0,100);
  histo1D["AlignSystSF"] = new TH1F("AlignSystSF","AlignSystSF",200,-.001,.001);
  
  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////

  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("preselected"));
  CutsSelecTableSemiMu.push_back(string("trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
  CutsSelecTableSemiMu.push_back(string("Veto electron"));
  
  vector<string> CutsSelecTableSemiEl;
  CutsSelecTableSemiEl.push_back(string("preselected"));
  CutsSelecTableSemiEl.push_back(string("trigged"));
  CutsSelecTableSemiEl.push_back(string("Good PV"));
  CutsSelecTableSemiEl.push_back(string("1 selected electron"));
  CutsSelecTableSemiEl.push_back(string("Veto muon"));
  CutsSelecTableSemiEl.push_back(string("Veto 2nd electron from Z-decay"));
  CutsSelecTableSemiEl.push_back(string("Conversion veto"));
  
  char LabelNJets[100];
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  selecTableSemiMu.SetLuminosity(Luminosity);
  SelectionTable selecTableSemiEl(CutsSelecTableSemiEl, datasets);
  selecTableSemiEl.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  
  /////////////////////////////
  /// ResolutionFit Stuff
  /////////////////////////////
  
  bool CalculateResolutions = false; // If false, the resolutions will be loaded from a previous calculation
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
  
  /////////////////////////////////////
  // Initialize JetCombination Stuff //
  /////////////////////////////////////
	bool Tprime = false; // If false, regular variables are used in MVA

  bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff
  string MVAmethod = "Likelihood"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)
    
  JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, Tprime);
  if (verbose > 0)
    cout << " - JetCombiner instantiated ..." << endl;
  
  /////////////////////////////////////////////////
  // Which analysis to execute? Top Mass or JES? //
  /////////////////////////////////////////////////
  
  bool measureTopMass = false;
  bool measureTopMassDifference = true;
  
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
    
    histo1D["mTop_LeptonPlus_"+dataSetName] = new TH1F(("mTop_LeptonPlus_"+dataSetName).c_str(),("mTop_LeptonPlus_"+dataSetName).c_str(),150,100,250);
    histo1D["mTop_LeptonMinus_"+dataSetName] = new TH1F(("mTop_LeptonMinus_"+dataSetName).c_str(),("mTop_LeptonMinus_"+dataSetName).c_str(),150,100,250);
    
    /////////////////////////////////////
   	/// Initialize JEC factors
   	/////////////////////////////////////
   	
    vector<JetCorrectorParameters> vCorrParam;
//    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) // Data!
//    {
//      JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("JECFiles/Jec11V2_db_AK5PFchs_L2L3Residual.txt");
//      vCorrParam.push_back(*ResJetCorPar);
//    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("JECFiles/START42_V17_AK5PFchs_Uncertainty.txt");
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, false);
    
    histo2D["jesUncPlus"] = new TH2F("jesUncPlus", "jesUncPlus", 51, 29.5, 80.5, 51, -2.55, 2.55);
    histo2D["jesUncMinus"] = new TH2F("jesUncMinus", "jesUncMinus", 51, 29.5, 80.5, 51, -2.55, 2.55);
    for(float jetPt=30; jetPt<=80; jetPt++)
    {
      for(float jetEta=-2.4; jetEta<=2.4; jetEta+=0.1)
      {
        jecUnc->setJetEta(jetEta);
        jecUnc->setJetPt(jetPt);
        histo2D["jesUncPlus"]->Fill(jetPt, jetEta, jecUnc->getUncertainty(true));
        jecUnc->setJetEta(jetEta);
        jecUnc->setJetPt(jetPt);
        histo2D["jesUncMinus"]->Fill(jetPt, jetEta, jecUnc->getUncertainty(false));
      }
    }
    
    ////////////////////////////////
    // LOAD THE FULLKINFIT OBJECT //
    ////////////////////////////////
    
    FullKinFit* kinFit = NULL;
    if( ! CalculateResolutions )
    {
      kinFit = new FullKinFit(datasets[d], resFitLightJets, resFitBJets, measureTopMass, measureTopMassDifference);
      kinFit->SetResFitL7(resFitLightJetsL7, resFitBJetsL7);
    }
    
    ////////////////////////////////////////////////////////////
    // CREATE OUTPUT FILE AND TTREE FOR STORAGE OF THE NTUPLE //
    ////////////////////////////////////////////////////////////
    
    string decayChannel;
    if(semiElectron && semiMuon) decayChannel = "SemiLep";
    else
    {
      if(semiMuon) decayChannel = "SemiMu";
      if(semiElectron) decayChannel = "SemiEl";
    }
    
    string monsterFileTitle = "Monsters/KinFit_Monsters_"+dataSetName+"_"+decayChannel+".root";
    if(measureTopMass)
      monsterFileTitle = "Monsters/KinFit_Monsters_TopMass_"+dataSetName+"_"+decayChannel+".root";
    else if(measureTopMassDifference)
      monsterFileTitle = "Monsters/KinFit_LightMonsters_TopMassDiff_"+dataSetName+"_"+decayChannel+".root";
    
    cout << "INFO: creating Monsters file "+monsterFileTitle << endl;
    
    TFile* MonsterFile = new TFile(monsterFileTitle.c_str(),"RECREATE");
    
    LightMonster* lightMonster = 0;
    
    TTree* MonsterTree = new TTree("MonsterTree","Tree containing the monsters");
    
    if(measureTopMassDifference)
      MonsterTree->Branch("TheLightMonster","LightMonster",&lightMonster);

    ////////////////////////////////////
    //	Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1, itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
//    for (unsigned int ievt = 0; ievt < 10; ievt++)
    {
      nEvents[d]++;
      if(ievt%1000 == 0)
        std::cout<<"Processing the "<<ievt<<"th event" <<flush<<"\r";
      
      //load event
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      vector<TRootGenJet*> genjets;
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
      {
        genjets = treeLoader.LoadGenJet(ievt);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
//      cout << "run: " << event->runId() << "  lumi: " << event->lumiBlockId() << "  event: " << event->eventId() << endl;
      
      // scale factors for the event
      float scaleFactor = 1.;
      float lumiWeight = 1.;
      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
        lumiWeight = 1;
      
      TRootGenEvent* genEvt = 0;
      // Load the GenEvent and calculate the branching ratio correction
      if(dataSetName.find("TTbarJets") == 0)
      {
        genEvt = treeLoader.LoadGenEvent(ievt);
        if( genEvt->isSemiLeptonic() )
					scaleFactor *= (0.108*9.)*(0.676*1.5);
			  else if( genEvt->isFullHadronic() )
			    scaleFactor *= (0.676*1.5)*(0.676*1.5);
			  else if( genEvt->isFullLeptonic() )
  		    scaleFactor *= (0.108*9.)*(0.108*9.);
  		  scaleFactor = 1.;
/*  		  if(genEvt->isSemiLeptonic( TRootGenEvent::kMuon ) || genEvt->isSemiLeptonic( TRootGenEvent::kElec ) )
  		  {
          vector<TRootMCParticle*> mcParticles = treeLoader.LoadMCPart(ievt);
          bool leptonPlus = false;
          for(unsigned int i=0; i<mcParticles.size(); i++)
          {
            if( mcParticles[i]->status() != 3) continue;
            
            if( mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 )
              leptonPlus = false;
            else if( mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 )
              leptonPlus = true;
            else if( mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 )
              leptonPlus = false;
            else if( mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 )
              leptonPlus = true;
          }
          TLorentzVector hadrB = genEvt->hadronicDecayB();
          TLorentzVector hadrQ1 = genEvt->hadronicDecayQuark();
          TLorentzVector hadrQ2 = genEvt->hadronicDecayQuarkBar();
          TLorentzVector hadrBGenJet, hadrQ1GenJet, hadrQ2GenJet;
          float minDRhadrB = 9999, minDRhadrQ1 = 9999, minDRhadrQ2 = 9999;
          for(unsigned int i=0; i<genjets.size(); i++)
          {
            if(hadrB.DeltaR(*genjets[i]) < minDRhadrB)
            {
              minDRhadrB = hadrB.DeltaR(*genjets[i]);
              hadrBGenJet = *genjets[i];
            }
            if(hadrQ1.DeltaR(*genjets[i]) < minDRhadrQ1)
            {
              minDRhadrQ1 = hadrQ1.DeltaR(*genjets[i]);
              hadrQ1GenJet = *genjets[i];
            }
            if(hadrQ2.DeltaR(*genjets[i]) < minDRhadrQ2)
            {
              minDRhadrQ2 = hadrQ2.DeltaR(*genjets[i]);
              hadrQ2GenJet = *genjets[i];
            }
          }
          if(hadrBGenJet.Pt() > 20. && hadrQ1GenJet.Pt() > 20. && hadrQ2GenJet.Pt() > 20. && ( hadrBGenJet.Pt() != hadrQ1GenJet.Pt() ) && ( hadrBGenJet.Pt() != hadrQ2GenJet.Pt() ) && ( hadrQ1GenJet.Pt() != hadrQ2GenJet.Pt() ))
          { // fill the plots!!
            float mTop = (hadrBGenJet+hadrQ1GenJet+hadrQ2GenJet).M();
            if(leptonPlus) histo1D["mTop_LeptonPlus_"+dataSetName]->Fill( mTop );
            else histo1D["mTop_LeptonMinus_"+dataSetName]->Fill( mTop );
          }
  		  }*/
      }
      
      // Clone the init_jets vector, otherwise the corrections will be removed
      for(unsigned int i=0; i<init_jets_corrected.size(); i++)
        if(init_jets_corrected[i]) delete init_jets_corrected[i];
      init_jets_corrected.clear();

      for(unsigned int i=0; i<init_jets.size(); i++)
        init_jets_corrected.push_back( (TRootJet*) init_jets[i]->Clone() );
      
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
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v1"), currentRun, iFile);
          else if( event->runId() >= 161217 && event->runId() <= 163261 )
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v2"), currentRun, iFile);
          else if( event->runId() >= 163270 && event->runId() <= 163869 )
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v4"), currentRun, iFile);
          else if( event->runId() >= 165088 && event->runId() <= 165633 )
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_Mu17_TriCentralJet30_v5"), currentRun, iFile);
          else if( event->runId() == 166346 )
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_TriCentralJet30_v2"), currentRun, iFile);
          else if( event->runId() >= 165970 && event->runId() <= 167043 )
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_TriCentralJet30_v1"), currentRun, iFile);
          else if( event->runId() >= 167078 && event->runId() <= 167913 )
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_TriCentralJet30_v3"), currentRun, iFile);
          else if( event->runId() >= 170826 && event->runId() <= 173198 )
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_TriCentralJet30_v5"), currentRun, iFile);
          else if( event->runId() >= 173236 && event->runId() <= 178380 )
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralJet30_v1"), currentRun, iFile);
          else if( event->runId() >= 178381 && event->runId() <= 179889 )
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v2"), currentRun, iFile);
          else if( event->runId() >= 179959 && event->runId() <= 180252 )
            itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v3"), currentRun, iFile);
          else
            cout << "Unknown run for SemiMu HLTpath selection: " << event->runId() << endl;
          if( itriggerSemiMu == 9999 )
          {
            cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
            exit(-1);
          }
        }
        else
        {
          itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralJet30_v1"), currentRun);
//      	  if (itriggerSemiMu == 9999)
//      	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v3"), currentRun); // Fall11 44X Chamonix
        }
        
        // semi-electron
        if(dataSetName.find("Data_ElectronHad") == 0 || dataSetName.find("data_ElectronHad") == 0 || dataSetName.find("DATA_ElectronHad") == 0 )
        {
          if( event->runId() <= 161176 )
            itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v1"), currentRun, iFile);
          else if( event->runId() >= 161177 && event->runId() <= 163261 )
            itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v2"), currentRun, iFile);
          else if( event->runId() >= 163262 && event->runId() <= 165633 )
            itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30_v3"), currentRun, iFile);
          else if( event->runId() >= 165970 && event->runId() <= 166967 )
            itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v1"), currentRun, iFile);
          else if( event->runId() >= 167039 && event->runId() <= 167913 )
            itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v2"), currentRun, iFile);
          else if( event->runId() >= 170826 && event->runId() <= 173198 )
            itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v4"), currentRun, iFile);
          else if( event->runId() >= 173236 && event->runId() <= 178380 )
            itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v5"), currentRun, iFile);
          else if( event->runId() >= 178381 && event->runId() <= 179889 )
            itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v2"), currentRun, iFile);
          else if( event->runId() >= 179959 && event->runId() <= 180252 )
            itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v3"), currentRun, iFile);
          else
            cout << "Unknown run for SemiEl HLTpath selection: " << event->runId() << endl;
          if( itriggerSemiEl == 9999 )
          {
            cout << "itriggerSemiEl == 9999 for SemiEl HLTpath selection: " << event->runId() << endl;
            exit(-1);
          }
        }
        else
        {
          itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v5"), currentRun);
//      	  if (itriggerSemiEl == 9999)
//      	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v3"), currentRun); // Fall11 44X Chamonix
        }
      }
      // Apply Jet Corrections on-the-fly
//      if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
//        jetTools->correctJets(init_jets_corrected, vertex);
      
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
      {
        if(systematic == "JERPlus")
          jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus", false);
        else if(systematic == "JERMinus")
          jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus", false);
        else
          jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal", false);
        
        // Correct jets for JES uncertainy systematics
        if(systematic == "JESPlus")
          jetTools->correctJetJESUnc(init_jets_corrected, mets[0], "plus");
        else if(systematic == "JESMinus")
          jetTools->correctJetJESUnc(init_jets_corrected, mets[0], "minus");
        
        //Scale jets with a certain factor
//        jetTools->scaleJets(init_jets_corrected, 1.);
      }
      
      for(unsigned i=0; i<init_muons.size(); i++)
        MSPlot["AllMuonsRelPFIso"]->Fill((init_muons[i]->chargedHadronIso()+init_muons[i]->neutralHadronIso()+init_muons[i]->photonIso())/init_muons[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      for(unsigned i=0; i<init_electrons.size(); i++)
        MSPlot["AllElectronsRelPFIso"]->Fill((init_electrons[i]->chargedHadronIso()+init_electrons[i]->neutralHadronIso()+init_electrons[i]->photonIso())/init_electrons[i]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      
      if(systematic == "AlignPlus" || systematic == "AlignMinus")
      {
        // Scale PFJets up/down according to their charge (for tracker misalignment studies)
        for(unsigned int iJet=0; iJet<init_jets_corrected.size(); iJet++)
        {
          TRootPFJet* jet = jetTools->convertToPFJets(init_jets_corrected[iJet]);
          if(jet->chargedMultiplicity() > 0)
          {
            float chargedFraction = jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction() + jet->chargedMuEnergyFraction();
            float chargedAveragePt = jet->Pt() * chargedFraction / jet->chargedMultiplicity();
            float charge = jet->charge();
            float deltaPtFraction = ( chargedAveragePt * chargedAveragePt * charge * 0.0001 ) / jet->Pt();
            if( systematic == "AlignMinus" ) deltaPtFraction *= -1.;
            init_jets_corrected[iJet]->SetPxPyPzE(jet->Px()*(1+deltaPtFraction), jet->Py()*(1+deltaPtFraction), jet->Pz()*(1+deltaPtFraction), jet->E()*(1+deltaPtFraction));
            histo1D["AlignSystSF"]->Fill(deltaPtFraction);
          }
        }
      }
      
      /////////////////////////////
      //   Selection
      /////////////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
      selection.setJetCuts(30.,2.4,0.01,1.,0.98,0.3,0.1);
      selection.setMuonCuts(20,2.1,0.125,10,0.02,0.3,1,1,1);
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setElectronCuts(30,2.5,0.1,0.02,1,0.3); // FIXME: RelIso to 0.1
      selection.setLooseElectronCuts(15,2.5,0.2); // semiMu looseElectron cuts
      
      bool triggedSemiMu = false;
      if( ! (dataSetName.find("Data_ElectronHad") == 0 || dataSetName.find("data_ElectronHad") == 0 || dataSetName.find("DATA_ElectronHad") == 0) )
        triggedSemiMu = treeLoader.EventTrigged (itriggerSemiMu);
      bool triggedSemiEl = false;
      if( ! (dataSetName.find("Data_MuHad") == 0 || dataSetName.find("data_MuHad") == 0 || dataSetName.find("DATA_MuHad") == 0) )
        triggedSemiEl = treeLoader.EventTrigged (itriggerSemiEl);
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);
      
      vector<TRootJet*> selectedJets, selectedJetsNoMu, selectedJetsNoEl;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons();
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(false);
      vector<TRootElectron*> vetoElectronsSemiEl = selection.GetSelectedLooseElectrons(20,2.5,0.2,true);
      if( systematic == "InvertedIso" )
      {
        vector<TRootMuon*> overlapMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0]);
        vector<TRootElectron*> overlapElectrons = selection.GetSelectedElectronsInvIso(0.2, vertex[0]);
        selectedJetsNoMu = selection.GetSelectedJets(overlapMuons,true);
        selectedJetsNoEl = selection.GetSelectedJets(overlapElectrons,true);
        selectedMuons = selection.GetSelectedMuonsInvIso(0.2, vertex[0], selectedJetsNoMu);
        selectedElectrons = selection.GetSelectedElectronsInvIso(0.2, vertex[0],selectedJetsNoEl);
      }
      else // Normal selection
      {
        selectedJets = selection.GetSelectedJets(true);
        selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
        selectedElectrons = selection.GetSelectedElectrons(vertex[0],selectedJets);
      }
      
      MSPlot["nPrimaryVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
   		MSPlot["nGoodPrimaryVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nSelectedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nSelectedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLooseOtherMuons"]->Fill(vetoMuons.size()-1, datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLooseElectronsSemiMu"]->Fill(vetoElectronsSemiMu.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nLooseElectronsSemiEl"]->Fill(vetoElectronsSemiEl.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["nSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      
      bool eventSelectedSemiMu = false;
      bool eventSelectedSemiEl = false;
      
      MSPlot["nEventsAfterCutsSemiMu"]->Fill(0, datasets[d], true, Luminosity*scaleFactor);
      selecTableSemiMu.Fill(d,0,scaleFactor*lumiWeight);
      if( triggedSemiMu && semiMuon )
      {
        MSPlot["nEventsAfterCutsSemiMu"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
        selecTableSemiMu.Fill(d,1,scaleFactor*lumiWeight);
        if( isGoodPV )
        {
          MSPlot["nEventsAfterCutsSemiMu"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
          selecTableSemiMu.Fill(d,2,scaleFactor*lumiWeight);
          if( selectedMuons.size() == 1 )
          {
            MSPlot["nEventsAfterCutsSemiMu"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
  		      selecTableSemiMu.Fill(d,3,scaleFactor*lumiWeight);
        		if( vetoMuons.size() == 1 || ( systematic == "InvertedIso" && vetoMuons.size() == 0 ) ) // if InvertedIso, selected muon not part of vetoMuons vector!
        		{
        		  MSPlot["nEventsAfterCutsSemiMu"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
              selecTableSemiMu.Fill(d,4,scaleFactor*lumiWeight);
              if( vetoElectronsSemiMu.size() == 0 )
              {
                MSPlot["nEventsAfterCutsSemiMu"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
                selecTableSemiMu.Fill(d,5,scaleFactor*lumiWeight);
                if(systematic == "InvertedIso") selectedJets = selectedJetsNoMu;
                if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-3)
                {
                  MSPlot["nEventsAfterCutsSemiMu"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
                  selecTableSemiMu.Fill(d,6,scaleFactor*lumiWeight);
                  if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-2)
                  {
                    MSPlot["nEventsAfterCutsSemiMu"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
                    selecTableSemiMu.Fill(d,7,scaleFactor*lumiWeight);
                    if(selectedJets.size()>=(unsigned int)anaEnv.NofJets-1)
                    {
                      MSPlot["nEventsAfterCutsSemiMu"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
                      selecTableSemiMu.Fill(d,8,scaleFactor*lumiWeight);
                      if(selectedJets.size()>=(unsigned int)anaEnv.NofJets)
                      {
                        if(selectedJets[2]->Pt() > 30)
                        {
                          histo1D["FourthJetPt"]->Fill(selectedJets[3]->Pt());
                          if(triggedSemiMu) histo1D["FourthJetPtTriggered"]->Fill(selectedJets[3]->Pt());
                        }
                        MSPlot["nEventsAfterCutsSemiMu"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
                        selecTableSemiMu.Fill(d,9,scaleFactor*lumiWeight);
                        eventSelectedSemiMu = true;
                        float reliso = (selectedMuons[0]->chargedHadronIso()+selectedMuons[0]->neutralHadronIso()+selectedMuons[0]->photonIso())/selectedMuons[0]->Pt();
                        MSPlot["SelectedEventsMuonsRelPFIso"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      
      MSPlot["nEventsAfterCutsSemiEl"]->Fill(0, datasets[d], true, Luminosity*scaleFactor);
      selecTableSemiEl.Fill(d,0,scaleFactor*lumiWeight);
      if( semiElectron && triggedSemiEl )
      {
        MSPlot["nEventsAfterCutsSemiEl"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
        selecTableSemiEl.Fill(d,1,scaleFactor*lumiWeight);
        if( isGoodPV )
        {
          MSPlot["nEventsAfterCutsSemiEl"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
          selecTableSemiEl.Fill(d,2,scaleFactor*lumiWeight);
          if( selectedElectrons.size() == 1 )
          {
            MSPlot["nEventsAfterCutsSemiEl"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
            selecTableSemiEl.Fill(d,3,scaleFactor*lumiWeight);
            if( vetoMuons.size() == 0 )
            {
              MSPlot["nEventsAfterCutsSemiEl"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
              selecTableSemiEl.Fill(d,4,scaleFactor*lumiWeight);
              if( !selection.foundZCandidate(selectedElectrons[0], vetoElectronsSemiEl) )
              {
                MSPlot["nEventsAfterCutsSemiEl"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
                selecTableSemiEl.Fill(d,5,scaleFactor*lumiWeight);
                if( selection.passConversionRejection(selectedElectrons[0]) )
                {
                  MSPlot["nEventsAfterCutsSemiEl"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
                  selecTableSemiEl.Fill(d,6,scaleFactor*lumiWeight);
                  if(systematic == "InvertedIso") selectedJets = selectedJetsNoEl;
                  if( selectedJets.size()>=(unsigned int)anaEnv.NofJets-3 )
                  {
                    MSPlot["nEventsAfterCutsSemiEl"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
                    selecTableSemiEl.Fill(d,7,scaleFactor*lumiWeight);
                    if( selectedJets.size()>=(unsigned int)anaEnv.NofJets-2 )
                    {
                      MSPlot["nEventsAfterCutsSemiEl"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
                      selecTableSemiEl.Fill(d,8,scaleFactor*lumiWeight);
                      if( selectedJets.size()>=(unsigned int)anaEnv.NofJets-1 )
                      {
                        MSPlot["nEventsAfterCutsSemiEl"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
                        selecTableSemiEl.Fill(d,9,scaleFactor*lumiWeight);
                        if( selectedJets.size()>=(unsigned int)anaEnv.NofJets )
                        {
                          if(selectedJets[2]->Pt() > 30)
                          {
                            histo1D["FourthJetPt"]->Fill(selectedJets[3]->Pt());
                            if(triggedSemiEl) histo1D["FourthJetPtTriggered"]->Fill(selectedJets[3]->Pt());
                          }
                          MSPlot["nEventsAfterCutsSemiEl"]->Fill(10, datasets[d], true, Luminosity*scaleFactor);
                          selecTableSemiEl.Fill(d,10,scaleFactor*lumiWeight);
                          eventSelectedSemiEl = true;
                          float reliso = (selectedElectrons[0]->chargedHadronIso()+selectedElectrons[0]->neutralHadronIso()+selectedElectrons[0]->photonIso())/selectedElectrons[0]->Pt();
                          MSPlot["SelectedEventsElectronsRelPFIso"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
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
      
      if( ! ( eventSelectedSemiEl || eventSelectedSemiMu ) ) continue;
      
      if( eventSelectedSemiEl && eventSelectedSemiMu )
        cout << "Event selected in semiEl and semiMu channel???" << endl;
      
//      continue;
        
      vector<TRootMCParticle*> mcParticles;
      if( dataSetName.find("TTbarJets") == 0 || dataSetName.find("TT_") == 0 )
      {
        mcParticles = treeLoader.LoadMCPart(ievt);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      if( eventSelectedSemiMu )
        jetCombiner->ProcessEvent(datasets[d], mcParticles, selectedJets, selectedMuons[0], init_electrons, init_muons, genEvt, scaleFactor);
      else
        jetCombiner->ProcessEvent(datasets[d], mcParticles, selectedJets, selectedElectrons[0], init_electrons, init_muons, genEvt, scaleFactor);
      
      if(CalculateResolutions && dataSetName.find("TTbarJets") == 0)
        jetCombiner->FillResolutions(resFitLightJets, resFitBJets, resFitBJets_B, resFitBJets_Bbar);
      
      if ( !TrainMVA && !CalculateResolutions )
      {
        //get the MC matched jet combination, not the MVA best matched
	      vector<unsigned int> goodCombi = jetCombiner->GetGoodJetCombination();
        
				///////////////
	      // KINFITTER //
	      ///////////////

   	    kinFit->SetJets(selectedJets);
   	    
   	    vector<TH2F> monsterVector;
   	    vector< float > mvaValsVector;
   	    vector< vector<unsigned int> > mvaResultsVector;
   	    vector< vector< float > > topMassVector; // mTopFit, sigmaMTopFit, chi2MTopFit
   	    for(unsigned int iCombi=0; iCombi<12; iCombi++)
   	    {
   	      pair<float, vector<unsigned int> > tmpMvaVals = jetCombiner->getMVAValue(MVAmethod, iCombi+1);
   	      mvaResultsVector.push_back(tmpMvaVals.second);
   	      mvaValsVector.push_back(tmpMvaVals.first);
          kinFit->SetMVAStuff(tmpMvaVals);
          
          if(measureTopMassDifference)
          {
            vector<float> tmp;
            float* res = 0;
//              if( goodCombi[2] == tmpMvaVals.second[2] && ( ( goodCombi[1] == tmpMvaVals.second[1] && goodCombi[0] == tmpMvaVals.second[0] )
//                 || ( goodCombi[0] == tmpMvaVals.second[1] && goodCombi[1] == tmpMvaVals.second[0] ) ) )
            res = kinFit->EstimateTopMass(event, 80.4, false, iCombi);
//              else
//              {
//                res = new float[3];
//                res[0] = 99.;
//                res[1] = 99.;
//                res[2] = 99.;
//              }
            tmp.push_back(res[0]);
            tmp.push_back(res[1]);
            tmp.push_back(res[2]);
            topMassVector.push_back(tmp);
            delete res;
          }
          else
          {
        	  TH2F* histo = 0;
        	  if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
        	    histo = kinFit->FitEvent(event, 80.4, 173.3, false, iCombi); // As measured by the Tevatron //switch first boolean to true to save all monsters
            else
              histo = kinFit->FitEvent(event, 80.4, 172.5, false, iCombi); // As used in the MC
            monsterVector.push_back(*histo);
            delete histo;
          }
        }
//          continue;
        vector<unsigned int> mcJetCombi = jetCombiner->GetGoodJetCombination();
        int hadrBJetIndex = mcJetCombi[2], lightJet1Index = mcJetCombi[0], lightJet2Index = mcJetCombi[1], leptBJetIndex = mcJetCombi[3];

        vector<float> bTagTCHE, bTagTCHP, bTagSSVHE, bTagSSVHP;
        vector<TLorentzVector> otherSelectedJets;
        for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++)
        {
          otherSelectedJets.push_back( *selectedJets[iJet] );
          bTagTCHE.push_back(selectedJets[iJet]->btag_trackCountingHighEffBJetTags());
          bTagTCHP.push_back(selectedJets[iJet]->btag_trackCountingHighPurBJetTags());
          bTagSSVHE.push_back(selectedJets[iJet]->btag_simpleSecondaryVertexHighEffBJetTags());
          bTagSSVHP.push_back(selectedJets[iJet]->btag_simpleSecondaryVertexHighPurBJetTags());
        }

        if(measureTopMassDifference)
        {
          // check top quark masses and which top is decaying hadronically
          float topMass = -1., antiTopMass = -1.;
          bool topDecayedLept = false;
          
          if(dataSetName.find("TTbarJets") == 0 || dataSetName.find("TT_") == 0)
          {
            for(unsigned int iPart=0; iPart<mcParticles.size(); iPart++)
            {
              if( mcParticles[iPart]->status() != 3) continue;
//                cout << "type: " << mcParticles[iPart]->type() << endl;
              if( mcParticles[iPart]->type() == 6 )
                topMass = mcParticles[iPart]->M();
              else if( mcParticles[iPart]->type() == -6 )
                antiTopMass = mcParticles[iPart]->M();
              else if( mcParticles[iPart]->type() == -13 && mcParticles[iPart]->motherType() == 24 && mcParticles[iPart]->grannyType() == 6 )
                topDecayedLept = true;
              else if( mcParticles[iPart]->type() == -11 && mcParticles[iPart]->motherType() == 24 && mcParticles[iPart]->grannyType() == 6 )
                topDecayedLept = true;
            }
          }
          
          lightMonster = new LightMonster();
          lightMonster->setEventID( event->eventId() );
          lightMonster->setRunID( event->runId() );
          lightMonster->setLumiBlockID( event->lumiBlockId() );
          lightMonster->setIdParton1( event->idParton1() );
          lightMonster->setXParton1( event->xParton1() );
          lightMonster->setIdParton2( event->idParton2() );
          lightMonster->setXParton2( event->xParton2() );
          lightMonster->setFactorizationScale( event->factorizationScale() );
          lightMonster->setNPV(vertex.size());
          lightMonster->setNPUBXm1(event->nPu(-1));
          lightMonster->setNPU(event->nPu(0));
          lightMonster->setNPUBXp1(event->nPu(1));
          lightMonster->setTopMass(topMass);
          lightMonster->setAntiTopMass(antiTopMass);
          lightMonster->setSelectedSemiMu(eventSelectedSemiMu);
          if(dataSetName.find("TTbarJets") == 0 || dataSetName.find("TT_") == 0)
          {
            lightMonster->setSemiMuDecay(genEvt->isSemiLeptonic( TRootGenEvent::kMuon ));
            lightMonster->setSemiElDecay(genEvt->isSemiLeptonic( TRootGenEvent::kElec ));
          }
          lightMonster->setTopDecayedLept(topDecayedLept);
          lightMonster->setAll4JetsMCMatched( jetCombiner->All4JetsMatched_MCdef() );
          lightMonster->setAllHadronicJetsMCMatched( jetCombiner->HadronicTopJetsMatched_MCdef() );
          lightMonster->setMvaVals(mvaValsVector);
          lightMonster->setMvaResults(mvaResultsVector);
          lightMonster->setEventWeight(scaleFactor);
          lightMonster->setMTopFitResults(topMassVector);
          lightMonster->setHadrBJet( hadrBJetIndex );
          lightMonster->setHadrLJet1( lightJet1Index );
          lightMonster->setHadrLJet2( lightJet2Index );
          lightMonster->setLeptBJet( leptBJetIndex );
          lightMonster->setMET( *mets[0] );
          lightMonster->setSelectedJets( otherSelectedJets );
          lightMonster->setBTagTCHE(bTagTCHE);
          lightMonster->setBTagTCHP(bTagTCHP);
          lightMonster->setBTagSSVHE(bTagSSVHE);
          lightMonster->setBTagSSVHP(bTagSSVHP);
          if( eventSelectedSemiMu )
          {
            lightMonster->setLepton( *selectedMuons[0] );
            lightMonster->setLeptonCharge( selectedMuons[0]->charge() );
            lightMonster->setLeptonPFRelIso( (selectedMuons[0]->chargedHadronIso()+selectedMuons[0]->neutralHadronIso()+selectedMuons[0]->photonIso())/selectedMuons[0]->Pt() );
          }
          else
          {
            lightMonster->setLepton( *selectedElectrons[0] );
            lightMonster->setLeptonCharge( selectedElectrons[0]->charge() );
            lightMonster->setLeptonPFRelIso( (selectedElectrons[0]->chargedHadronIso()+selectedElectrons[0]->neutralHadronIso()+selectedElectrons[0]->photonIso())/selectedElectrons[0]->Pt() );
          }
          lightMonster->setHadrBQuark( jetCombiner->GetHadrBQuark() );
          lightMonster->setHadrLQuark1( jetCombiner->GetLightQuark1() );
          lightMonster->setHadrLQuark2( jetCombiner->GetLightQuark2() );
          lightMonster->setLeptBQuark( jetCombiner->GetLeptBQuark() );
          
          MonsterTree->Fill();
          delete lightMonster;
	      }
      }// end !TrainMVA && eventSelected
    }				//loop on events
    cout<<endl;

    if( !CalculateResolutions ) kinFit->Write(fout, true, pathPNG+"FullKinFit/");
    delete kinFit;
    
    MonsterFile->cd();
    
    TTree *configTreeMonsterFile = new TTree("configTreeMonsterFile","configuration Tree in Monster File");
    TClonesArray* tcdatasetmonsterfile = new TClonesArray("Dataset",1);
    configTreeMonsterFile->Branch("Dataset","TClonesArray",&tcdatasetmonsterfile);
    new ((*tcdatasetmonsterfile)[0]) Dataset(*datasets[d]);
    
    configTreeMonsterFile->Fill();
    configTreeMonsterFile->Write();
    MonsterTree->Write();
    MonsterFile->Close();
    delete MonsterFile;
    
    if(jetTools) delete jetTools;
    
    treeLoader.UnLoadDataset(); //important: free memory
    
/*    if(!CalculateResolutions)
    {
      histo1D["mTop_LeptonPlus_"+dataSetName]->Fit("gaus", "QR","",160,180);
      histo1D["mTop_LeptonMinus_"+dataSetName]->Fit("gaus", "QR","",160,180);
    
      TF1* funcPlus = histo1D["mTop_LeptonPlus_"+dataSetName]->GetFunction("gaus");
      TF1* funcMinus = histo1D["mTop_LeptonMinus_"+dataSetName]->GetFunction("gaus");
      cout << "mTop:\nleptonPlus: " << funcPlus->GetParameter(1) << " +/- " <<  funcPlus->GetParError(1);
      cout << "  leptonMinus: " << funcMinus->GetParameter(1) << " +/- " <<  funcMinus->GetParError(1) << endl;
      float error = sqrt(funcPlus->GetParError(1)*funcPlus->GetParError(1) + funcMinus->GetParError(1)*funcMinus->GetParError(1));
      cout << "mTopDiff:  " << funcPlus->GetParameter(1) - funcMinus->GetParameter(1) << " +/- " << error << endl;
    }*/
  }				//loop on datasets
  
  //Once everything is filled ...
  if (verbose > 0)
    cout << " We ran over all the data ;-)" << endl;
  
  //Selection tables
  selecTableSemiMu.TableCalculator(false, true, true, true, true);
  selecTableSemiEl.TableCalculator(false, true, true, true, true);
  string selectiontableSemiMu = "SelectionTable_SemiMu_JES.tex";
  string selectiontableSemiEl = "SelectionTable_SemiEl_JES.tex";
	selecTableSemiMu.Write(selectiontableSemiMu.c_str());
	selecTableSemiEl.Write(selectiontableSemiEl.c_str());

  // Do some special things with certain plots (normalize, BayesDivide, ... )
  if (verbose > 0)
    cout << "Treating the special plots." << endl;
  
  ///////////////////
  // Writing
  //////////////////
  if (verbose > 1)
  	cout << " - Writing outputs to the files ..." << endl;

	string pathPNGJetCombi = pathPNG+"JetCombination/";

  mkdir((pathPNG+"MSPlot/").c_str(),0777);
  mkdir(pathPNGJetCombi.c_str(),0777);
  
  jetCombiner->Write(fout, true, pathPNGJetCombi, false);
  
  // Fill the resolution histograms and calculate the resolutions
  if(CalculateResolutions)
  {
    cout << "Fit the resolution stuff..." << endl;
    mkdir((pathPNG+"resFit_LightJet/").c_str(),0777);
    mkdir((pathPNG+"resFit_BJet/").c_str(),0777);
    mkdir((pathPNG+"resFit_BJet_B/").c_str(),0777);
    mkdir((pathPNG+"resFit_BJet_Bbar/").c_str(),0777);
  
    resFitLightJets->WritePlots(fout, true, pathPNG+"resFit_LightJet/");
    resFitLightJets->WriteResolutions("lightJetReso.root");
    resFitBJets->WritePlots(fout, true, pathPNG+"resFit_BJet/");
    resFitBJets->WriteResolutions("bJetReso.root");
    resFitBJets_B->WritePlots(fout, true, pathPNG+"resFit_BJet_B/");
    resFitBJets_B->WriteResolutions("bJetReso_B.root");
    resFitBJets_Bbar->WritePlots(fout, true, pathPNG+"resFit_BJet_Bbar/");
    resFitBJets_Bbar->WriteResolutions("bJetReso_Bbar.root");
  }
  
  cout << "Writing the histograms..." << endl;
  // 1D 
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(false, name, true, true, true, true, true);
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
  configTree->Fill();
  configTree->Write();

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
  if(jetCombiner) delete jetCombiner;
  
  delete fout;
  delete tcdatasets;
  delete configTree;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;

  return 0;
}
