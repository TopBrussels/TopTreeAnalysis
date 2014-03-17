#include "Headers.h"

using namespace std;
using namespace TopTree;
using namespace reweight;

//typedef pair< string , vector< pair<string,float> > > MSIntegral;
typedef map<string,MultiSamplePlot*> MapMSP;
typedef map<string,MapMSP>           MapMapMSP;
typedef map<string,bool>             MapAnaCut;

///////////////////////
// DECLARE FUNCTIONS //
///////////////////////

int analysis(string, TString, int, float, float, float, int, float, float);

//MSIntegral eval_msi(MultiSamplePlot*);

int SelectCompute(TRootEvent* event, Selection selection, int verbose,
		  vector<TRootVertex*> vertex , vector<TRootPFJet*> selectedJetsPF, 
		  int cut_nJets, float cut_pt1, float cut_pt2, float cut_deltaphi,
		  int m_nJet, float* jet_E, float* jet_pt, float* jet_phi, float* jet_eta,  
		  float* neutral_had_E_frac, float* neutral_em_E_frac, float* charged_em_E_frac, float* charged_had_E_frac, float* charged_mu_E_frac,
		  float* charged_mult, float* neutral_mult, float* muon_mult,
		  float* tot_em_E_frac, float* tot_neutral_E_frac, float* tot_charged_E_frac, float* tot_chargedMult,
		  float& charged_E_frac_dijet, float& chargedMult_dijet, float& mass_dijet, float& em_E_frac_dijet,
		  float& rho, int& nVertices, int& nJetInEvt, float& DeltaPhi);

int FillHistos(MapMSP MSPlot, string name_MSP,
	       Dataset* dataset,
	       int verbose, float Luminosity, float scaleFactor,
	       int m_nJet, float* jet_E, float* jet_pt, float* jet_phi, float* jet_eta,  
	       float* neutral_had_E_frac, float* neutral_em_E_frac, float* charged_em_E_frac, float* charged_had_E_frac, float* charged_mu_E_frac,
	       float* charged_mult, float* neutral_mult, float* muon_mult,
	       float* tot_em_E_frac, float* tot_neutral_E_frac, float* tot_charged_E_frac, float* tot_chargedMult,
	       float charged_E_frac_dijet, float chargedMult_dijet, float mass_dijet, float em_E_frac_dijet,
	       float rho, int nVertices, int nJetInEvt, float DeltaPhi);

int LocalFill(MapMSP MSPlot, string name, float value, Dataset* dataset, bool scale, float Luminosity, float scaleFactor);

//////////
// MAIN //
//////////

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

//////////////
// ANALYSIS //
//////////////

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

  // Outlog
  ofstream outlog(path+"/log.txt",ios::out);

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
  float rho=0;
  int nVertices=0;
  int nJetInEvt=0;

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
  const int m_nJet=2;

  float jet_E[  m_nJet]={0,0};
  float jet_pt[ m_nJet]={0,0};
  float jet_phi[m_nJet]={0,0};
  float jet_eta[m_nJet]={0,0};
  float DeltaPhi=0;
  //
  float neutral_had_E_frac[m_nJet]={0,0};
  float neutral_em_E_frac[ m_nJet]={0,0};
  float charged_had_E_frac[m_nJet]={0,0};
  float charged_em_E_frac[ m_nJet]={0,0};
  float charged_mu_E_frac[ m_nJet]={0,0};
  // zia
  float charged_mult[m_nJet]={0,0};
  float neutral_mult[m_nJet]={0,0};
  float muon_mult[   m_nJet]={0,0};
  //
  float tot_em_E_frac[     m_nJet]={0,0};
  float tot_neutral_E_frac[m_nJet]={0,0};
  float tot_charged_E_frac[m_nJet]={0,0};
  float tot_chargedMult[   m_nJet]={0,0};
  //
  float charged_E_frac_dijet=0, chargedMult_dijet=0, mass_dijet=0, em_E_frac_dijet=0;

  // SelectCompute output
  int outSelectCompute=-1;
  int outfill=-1;

  ///////////////////
  // DECLARE PLOTS //
  ///////////////////

  MapMapMSP MSPlot;
  MapAnaCut AnaCuts;

  const int nMSP=18;
  string name_MSP[nMSP]={"AnaCut_none","AnaCut_2jets", "AnaCut_dphi", "AnaCut_pt1", "AnaCut_pt2", "AnaCut_kine",
			 "AnaCut_jet1_emf0p2","AnaCut_jet2_emf0p2","AnaCut_pair_emf0p2","AnaCut_both_emf0p2",
			 "AnaCut_jet1_chf0p1","AnaCut_jet2_chf0p1","AnaCut_pair_chf0p1","AnaCut_both_chf0p1",
			 "AnaCut_jet1_chmult","AnaCut_jet2_chmult","AnaCut_pair_chmult","AnaCut_both_chmult"};

  // Objects to compute background rejection
  MultiSamplePlot::MSIntegral msi[nMSP];
  float val_before=0, val_after=0, reduction=0;
  string name_before="", name_after="";
  //vector<float> reductions;
  
  for(u_int iMSP=0 ; iMSP<nMSP ; iMSP++) {
    // analysis cuts within boolean map
    AnaCuts[name_MSP[iMSP]] = false;

    // Event
    MSPlot[name_MSP[iMSP]]["RhoCorrection"]         = new MultiSamplePlot(datasets, "RhoCorrection"+name_MSP[iMSP], 100, 0, 100, "#rho");
    MSPlot[name_MSP[iMSP]]["NbOfVertices"]          = new MultiSamplePlot(datasets, "NbOfVertices"+name_MSP[iMSP], 40, 0, 40, "Nb. of vertices");
    //
    MSPlot[name_MSP[iMSP]]["NbOfSelectedJets"]      = new MultiSamplePlot(datasets, "NbOfSelectedJets"+name_MSP[iMSP], 15, 0, 15, "Nb. of jets");

    // Jet Variables
    string nameJet="";
    string nameJetH="";
    string titleJet[m_nJet]={"Leading Jet", "Sub-Leading Jet"};

    MSPlot[name_MSP[iMSP]]["em_E_frac_dijet"     ] = new MultiSamplePlot(datasets, "em_E_frac_dijet"+name_MSP[iMSP]     , 100, 0, 1, "dijet em energy fraction");
    MSPlot[name_MSP[iMSP]]["charged_E_frac_dijet"] = new MultiSamplePlot(datasets, "charged_E_frac_dijet"+name_MSP[iMSP], 100, 0, 1, "dijet charged energy fraction");
    MSPlot[name_MSP[iMSP]]["chargedMult_dijet"   ] = new MultiSamplePlot(datasets, "chargedMult_dijet"+name_MSP[iMSP],    100, 0, 100,"dijet charged multiplicity");
    MSPlot[name_MSP[iMSP]]["mass_dijet"          ] = new MultiSamplePlot(datasets, "mass_dijet"+name_MSP[iMSP],           6000, 0, 3000,"dijet mass (GeV)");
    
    for(int iJ=0 ; iJ<m_nJet ; iJ++) {

      nameJet  = Form("_Jet%d",iJ+1);
      nameJetH = nameJet + "_" + name_MSP[iMSP];

      cout << "--- nameJet=" << nameJet << endl;
    
      MSPlot[name_MSP[iMSP]]["Eta"+nameJet] = new MultiSamplePlot(datasets, "Eta"+nameJetH, 100,  -5, 5,     titleJet[iJ]+" #eta");
      MSPlot[name_MSP[iMSP]]["Phi"+nameJet] = new MultiSamplePlot(datasets, "Phi"+nameJetH, 33,  -3.15, 3.15,titleJet[iJ]+" #phi");
      MSPlot[name_MSP[iMSP]]["Pt" +nameJet] = new MultiSamplePlot(datasets, "Pt" +nameJetH, 600, 0,  3000,  titleJet[iJ]+" p_{T} (GeV)");
      MSPlot[name_MSP[iMSP]]["E" +nameJet]  = new MultiSamplePlot(datasets, "E"  +nameJetH, 600, 0,  3000,  titleJet[iJ]+" E (GeV)");

      MSPlot[name_MSP[iMSP]]["chargedHad_E_frac"+nameJet] = new MultiSamplePlot(datasets, "chargedHad_E_frac"+nameJetH, 100, 0, 1,   titleJet[iJ]+" charged had energy fraction");
      MSPlot[name_MSP[iMSP]]["chargedEm_E_frac" +nameJet] = new MultiSamplePlot(datasets, "chargedEm_E_frac" +nameJetH, 100, 0, 1,   titleJet[iJ]+" charged EM  energy fraction");
      MSPlot[name_MSP[iMSP]]["chargedMu_E_frac" +nameJet] = new MultiSamplePlot(datasets, "chargedMu_E_frac" +nameJetH, 100, 0, 1,   titleJet[iJ]+" charged Mu  energy fraction");
      MSPlot[name_MSP[iMSP]]["neutralEm_E_frac" +nameJet] = new MultiSamplePlot(datasets, "neutralEm_E_frac" +nameJetH, 100, 0, 1,   titleJet[iJ]+" neutral EM  energy fraction");
      MSPlot[name_MSP[iMSP]]["neutralHad_E_frac"+nameJet] = new MultiSamplePlot(datasets, "neutralHad_E_frac"+nameJetH, 100, 0, 1,   titleJet[iJ]+" neutral had energy fraction");    

      MSPlot[name_MSP[iMSP]]["em_E_frac"     +nameJet]= new MultiSamplePlot(datasets, "em_E_frac"     +nameJetH, 100, 0, 1,   titleJet[iJ]+" tot EM energy fraction");
      MSPlot[name_MSP[iMSP]]["charged_E_frac"+nameJet]= new MultiSamplePlot(datasets, "charged_E_frac"+nameJetH, 100, 0, 1,   titleJet[iJ]+" tot charged energy fraction");
      MSPlot[name_MSP[iMSP]]["neutral_E_frac"+nameJet]= new MultiSamplePlot(datasets, "neutral_E_frac"+nameJetH, 100, 0, 1,   titleJet[iJ]+" tot neutral energy fraction");

      MSPlot[name_MSP[iMSP]]["chargedMult"      +nameJet] = new MultiSamplePlot(datasets, "chargedMult"      +nameJetH, 100, 0, 100,  titleJet[iJ]+" charged multiplicity");
      MSPlot[name_MSP[iMSP]]["neutralMult"      +nameJet] = new MultiSamplePlot(datasets, "neutralMult"      +nameJetH, 100, 0, 100,  titleJet[iJ]+" neutral multiplicity");
      MSPlot[name_MSP[iMSP]]["muonMult"         +nameJet] = new MultiSamplePlot(datasets, "muonMult"         +nameJetH, 100, 0, 100,  titleJet[iJ]+" muon multiplicity");

      MSPlot[name_MSP[iMSP]]["chargedMultTot"  +nameJet] = new MultiSamplePlot(datasets, "chargedMultTot"    +nameJetH, 100, 0, 100,  titleJet[iJ]+" total charged multiplicity");
    }
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
//       float rho = event->kt6PFJets_rho();
//       MSPlot["RhoCorrection"]->Fill(rho, datasets[d], true, Luminosity);

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

      if(verbose>2) cout << "-- enter object selection" << endl;

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
      
      // Apply inclusive selection
      // and use TRootPFJet methods to get members and put them within variables
      if(verbose>2) cout << "-- Call SelectCompute on the map : " << "noAnaCut" << endl;

//       outSelectCompute = 
// 	SelectCompute(event, selection, verbose,
// 		      vertex , selectedJetsPF, 
// 		      cut_nJets, cut_pt1, cut_pt2, cut_deltaphi,
// 		      m_nJet, jet_E, jet_pt, jet_phi, jet_eta,  
// 		      neutral_had_E_frac, neutral_em_E_frac, 
// 		      charged_em_E_frac, charged_had_E_frac, charged_mu_E_frac,
// 		      charged_mult, neutral_mult, muon_mult,
// 		      tot_em_E_frac, tot_neutral_E_frac, tot_charged_E_frac, tot_chargedMult,
// 		      charged_E_frac_dijet, chargedMult_dijet, mass_dijet, em_E_frac_dijet,
// 		      rho, nVertices, nJetInEvt, DeltaPhi);

      outSelectCompute = 
	SelectCompute(event, selection, verbose,
		      vertex , selectedJetsPF, 
		      -999, -999, -999, -999, // do no apply ANY cuts yet
		      m_nJet, jet_E, jet_pt, jet_phi, jet_eta,  
		      neutral_had_E_frac, neutral_em_E_frac, 
		      charged_em_E_frac, charged_had_E_frac, charged_mu_E_frac,
		      charged_mult, neutral_mult, muon_mult,
		      tot_em_E_frac, tot_neutral_E_frac, tot_charged_E_frac, tot_chargedMult,
		      charged_E_frac_dijet, chargedMult_dijet, mass_dijet, em_E_frac_dijet,
		      rho, nVertices, nJetInEvt, DeltaPhi);

      if(verbose>2 && ievt%1000==0) 
	cout << "%% DIJET : " << charged_E_frac_dijet 
	     << "  "          << chargedMult_dijet 
	     << "  "          << mass_dijet 
	     << "  "          << em_E_frac_dijet 
	     << endl;

      if(outSelectCompute==2) continue;
      
      // Analysis selections
      AnaCuts["AnaCut_none"] = true;
      //
      AnaCuts["AnaCut_2jets"] = nJetInEvt >= cut_nJets ;
      AnaCuts["AnaCut_dphi" ] = DeltaPhi  >= cut_deltaphi ;
      AnaCuts["AnaCut_pt1"  ] = jet_pt[0] >= cut_pt1 ;
      AnaCuts["AnaCut_pt2"  ] = jet_pt[1] >= cut_pt2 ;
      AnaCuts["AnaCut_kine" ] = AnaCuts["AnaCut_2jets"] && AnaCuts["AnaCut_dphi"] 
	&& AnaCuts["AnaCut_pt1"] && AnaCuts["AnaCut_pt2"];
      //
      AnaCuts["AnaCut_jet1_emf0p2"] = tot_em_E_frac[0] < 0.2 ;
      AnaCuts["AnaCut_jet2_emf0p2"] = tot_em_E_frac[1] < 0.2 ;
      AnaCuts["AnaCut_pair_emf0p2"] = em_E_frac_dijet  < 0.2 ;
      AnaCuts["AnaCut_both_emf0p2"] = AnaCuts["AnaCut_jet1_emf0p2"] && AnaCuts["AnaCut_jet2_emf0p2"] ;
      //
      AnaCuts["AnaCut_jet1_chf0p1"] = tot_charged_E_frac[0] < 0.1 ;
      AnaCuts["AnaCut_jet2_chf0p1"] = tot_charged_E_frac[1] < 0.1 ;
      AnaCuts["AnaCut_pair_chf0p1"] = charged_E_frac_dijet  < 0.1 ;
      AnaCuts["AnaCut_both_chf0p1"] = AnaCuts["AnaCut_jet1_chf0p1"] && AnaCuts["AnaCut_jet2_chf0p1"] ;
      //
      AnaCuts["AnaCut_jet1_chmult"] = tot_chargedMult[0]==0 ;
      AnaCuts["AnaCut_jet2_chmult"] = tot_chargedMult[1]==0 ;
      AnaCuts["AnaCut_pair_chmult"] = chargedMult_dijet ==0;
      AnaCuts["AnaCut_both_chmult"] = AnaCuts["AnaCut_jet1_chmult"] && AnaCuts["AnaCut_jet2_chmult"] ;
      
      // Fill various MapMSPlots depending on various selections
      for(MapMapMSP::const_iterator iMSP=MSPlot.begin() ; iMSP!=MSPlot.end() ; iMSP++) {

	if(AnaCuts[iMSP->first]) outfill = 
	  FillHistos(iMSP->second, iMSP->first, 
		     datasets[d],
		     verbose, Luminosity, scaleFactor,
		     m_nJet, jet_E, jet_pt, jet_phi, jet_eta,  
		     neutral_had_E_frac, neutral_em_E_frac, charged_em_E_frac, charged_had_E_frac, charged_mu_E_frac,
		     charged_mult, neutral_mult, muon_mult,
		     tot_em_E_frac, tot_neutral_E_frac, tot_charged_E_frac, tot_chargedMult,
		     charged_E_frac_dijet, chargedMult_dijet, mass_dijet, em_E_frac_dijet,
		     rho, nVertices, nJetInEvt, DeltaPhi);      
      }

    } // end loop over events

    // free memory
    //if(jetTools) delete jetTools;
    treeLoader.UnLoadDataset();

  } // end loop over datasets
 
  //////////////////////////////////
  // COMPUTE BACKGROUND REDUCTION //
  //////////////////////////////////

  // Reminder: typedef pair< string , vector< pair<string,float> > > MSIntegral;

  for(u_int iMSP=0 ; iMSP<nMSP ; iMSP++) {
    //msi[iMSP] = eval_msi( MSPlot[name_MSP[iMSP]]["NbOfVertices"] );
    msi[iMSP] = MSPlot[name_MSP[iMSP]]["NbOfVertices"] -> Integrate();
  }

  outlog << endl
	 << "$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl
	 << "$$ BACKGROUND REDUCTION $$" << endl;
  
  for(u_int iS=0 ; iS<(msi[0].second).size() ; iS++) {

    name_before = name_after = "";
    val_before = val_after = reduction = 0;

    name_before = ((msi[0].second)[iS]).first;  // name of the sample
    val_before  = ((msi[0].second)[iS]).second; // integral of the sample

    outlog << "$$ Sample #"   << iS << " " << name_before << " $$" << endl;
    
    for(u_int iMSP=0 ; iMSP<nMSP ; iMSP++) {
      
      if( (msi[iMSP].second).size() <= iS ) {
	outlog << "Error: missing AFTER sample #" << iS << endl;
	break;
      }
      
      name_after  = ((msi[iMSP] .second)[iS]).first;
      val_after   = ((msi[iMSP] .second)[iS]).second;
      
      reduction = val_after !=0 ? val_before / val_after : -999 ;

      outlog << "Selection : " << name_MSP[iMSP]
	     << " : before="   << val_before
	     << " ; after="    << val_after
	     << " ; reduction = " << reduction << " $$" << endl;
    }
  }
  outlog   << "$$$$$$$$$$$$$$$$$$$$$$$$$$"
	   << endl << endl;

  //////////////////
  // WRITE OUTPUT //
  //////////////////

  cout << "- Start writing output" << endl;
  fout->cd();

  cout << "-- Loop over MSPlots" << endl;

  for(MapMapMSP::const_iterator itMap = MSPlot.begin() ; itMap != MSPlot.end(); itMap++) {

    cout << "--- MAP : " << itMap->first << endl;

    for( MapMSP::const_iterator it    = (itMap->second).begin() ; it    != (itMap->second).end(); it++) {

      string name = (it->first)+"_"+(itMap->first);
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
	
	// ND
	if(verbose>2) cout << "--- draw" << endl;
	temp->Draw(name, 0, false, false, false, 0, x1, y1, x2, y2, magnify);

	if(verbose>2) cout << "--- write in pathPNG=" << pathPNG << endl;
	float maxY = temp->getMaxY();
	cout << "########## maxY=" << maxY << " ###" << endl;

	temp->Write(fout, name, true, pathPNG+"/", "png", magnifyLog); 
	//ND true => SaveAs the Canvas as image => seg fault probably caused by empty plots !
      }
    } // end loop over MSPlots
  }
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

//////////////
// EVAL_MSI //
//////////////

/** evaluate integrals of individual components in a MSPlot */
/*
MSIntegral eval_msi(MultiSamplePlot* plot)
{
  MSIntegral msi;
  vector<float> integrals;
  vector<string> names;
  float integral=0;
  float total=0;
  TH1F* h=0;

  // Error case
  if(!plot) {
    msi.first = "nothing";
    cout << "ERROR in eval_msi : no input plot" << endl;
    return msi;
  }

  // Fill the MSIntegral
  msi.first = plot->getplotName();
  names = plot->getTH1FNames();

  for(u_int i=0 ; i<names.size() ; i++) {
    h = plot->getTH1F(names[i]);
    integral = h!=0 ? h->Integral() : 0;
    (msi.second).push_back( make_pair(names[i],integral) );
    total += integral;
  }

  (msi.second).push_back( make_pair("total",total) );

  return msi;
}
*/
/////////////
// SelectCompute //
/////////////

/** Apply event and object selection then fill MSPlots. 
    Called within a loop over events. 
    Returns integers : 0 for correct execution, 
    2 for sending a "continue" instruction within the loop over events.
*/

int SelectCompute(TRootEvent* event, Selection selection, int verbose,
		  vector<TRootVertex*> vertex , vector<TRootPFJet*> selectedJetsPF, 
		  int cut_nJets, float cut_pt1, float cut_pt2, float cut_deltaphi,
		  int m_nJet, float* jet_E, float* jet_pt, float* jet_phi, float* jet_eta,  
		  float* neutral_had_E_frac, float* neutral_em_E_frac, float* charged_em_E_frac, float* charged_had_E_frac, float* charged_mu_E_frac,
		  float* charged_mult, float* neutral_mult, float* muon_mult,
		  float* tot_em_E_frac, float* tot_neutral_E_frac, float* tot_charged_E_frac, float* tot_chargedMult,
		  float& charged_E_frac_dijet, float& chargedMult_dijet, float& mass_dijet, float& em_E_frac_dijet,
		  float& rho, int& nVertices, int& nJetInEvt, float& DeltaPhi)
{

  if(verbose>2) cout << "--- enter the function" << endl;

  bool exist_dijet=false;

  rho       = event->kt6PFJets_rho();
  nVertices = vertex.size();

  /////////////////////
  // EVENT SELECTION //
  /////////////////////

  // Primary vertex
  bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
  if(!isGoodPV) return 2;

  // JETS

  // number of jets
  nJetInEvt = selectedJetsPF.size();
  if(nJetInEvt==0)        return 2;
  if(nJetInEvt<cut_nJets) return 2;
  if(nJetInEvt>=2)        exist_dijet=true;
  if(verbose>2) cout << "-- nJetInEvt=" << nJetInEvt << endl;

  // leading jet pt
  jet_pt[0]=selectedJetsPF[0]->Pt();
  if(jet_pt[0]<cut_pt1) return 2;

  // subleading jet pt
  if(exist_dijet) {
    jet_pt[1]=selectedJetsPF[1]->Pt();
    if(jet_pt[1]<cut_pt2) return 2;
  }
  else jet_pt[1]=-999;

  // delta phi
  jet_phi[0]=selectedJetsPF[0]->Phi();
  if(exist_dijet) {
    jet_phi[1]=selectedJetsPF[1]->Phi();
    DeltaPhi = TMath::Abs(jet_phi[0] - jet_phi[1]);
    if(DeltaPhi<cut_deltaphi) return 2;
  }
  else {
    jet_phi[2]=-999;
    DeltaPhi  =-999;
  }

  //////////////////////
  /// FILL HISTOGRAMS //
  //////////////////////

  if(verbose>2) cout << "--- look at jet pairs" << endl;

  for(int iJ=0 ; iJ<m_nJet ; iJ++) { // loop over jets

    if(nJetInEvt<=iJ) break;
	
    jet_E[  iJ]=selectedJetsPF[iJ]->E();
    jet_eta[iJ]=selectedJetsPF[iJ]->Eta();

    if(iJ>1) {
      jet_pt[ iJ] = selectedJetsPF[iJ]->Pt();
      jet_phi[iJ] = selectedJetsPF[iJ]->Phi();
    }

    // subgroups
    neutral_had_E_frac[iJ] = selectedJetsPF[iJ]->neutralHadronEnergyFraction();
    neutral_em_E_frac[ iJ] = selectedJetsPF[iJ]->neutralEmEnergyFraction();
    charged_had_E_frac[iJ] = selectedJetsPF[iJ]->chargedHadronEnergyFraction();
    charged_em_E_frac[ iJ] = selectedJetsPF[iJ]->chargedEmEnergyFraction();
    charged_mu_E_frac[ iJ] = selectedJetsPF[iJ]->chargedMuEnergyFraction();
    charged_mult[      iJ] = selectedJetsPF[iJ]->chargedMultiplicity();
    muon_mult[         iJ] = selectedJetsPF[iJ]->muonMultiplicity();
    neutral_mult[      iJ] = selectedJetsPF[iJ]->neutralMultiplicity();

    // sums
    tot_neutral_E_frac[iJ] = neutral_had_E_frac[iJ] + neutral_em_E_frac[iJ];
    tot_charged_E_frac[iJ] = charged_had_E_frac[iJ] + charged_em_E_frac[iJ] + charged_mu_E_frac[iJ];
    tot_em_E_frac[     iJ] = neutral_em_E_frac[ iJ] + charged_em_E_frac[iJ];
    tot_chargedMult[   iJ] = charged_mult[iJ] + muon_mult[iJ];

  } // end loop over jets

  if(verbose>2) cout << "--- end loop over jets" << endl
		     << "--- start building dijet quantities" << endl;

  // dijet
  TRootPFJet dijet;
  charged_E_frac_dijet=0; chargedMult_dijet=0; mass_dijet=0; em_E_frac_dijet=0;

  if(nJetInEvt>=2) {
    charged_E_frac_dijet = (jet_E[0]+jet_E[1]!=0) ? (jet_E[0]*tot_charged_E_frac[0]+jet_E[1]*tot_charged_E_frac[1]) / (jet_E[0]+jet_E[1]) : -999;
    em_E_frac_dijet      = (jet_E[0]+jet_E[1]!=0) ? (jet_E[0]*tot_em_E_frac[0]+jet_E[1]*tot_em_E_frac[1])           / (jet_E[0]+jet_E[1]) : -999;
    chargedMult_dijet = tot_chargedMult[0] + tot_chargedMult[1];
    dijet = (*selectedJetsPF[0]) + (*selectedJetsPF[1]) ;
    mass_dijet = dijet.M();
  }

  return 0;
}

int FillHistos(MapMSP MSPlot, string name_MSP,
	       Dataset* dataset,
	       int verbose, float Luminosity, float scaleFactor,
	       int m_nJet, float* jet_E, float* jet_pt, float* jet_phi, float* jet_eta,  
	       float* neutral_had_E_frac, float* neutral_em_E_frac, float* charged_em_E_frac, float* charged_had_E_frac, float* charged_mu_E_frac,
	       float* charged_mult, float* neutral_mult, float* muon_mult,
	       float* tot_em_E_frac, float* tot_neutral_E_frac, float* tot_charged_E_frac, float* tot_chargedMult,
	       float charged_E_frac_dijet, float chargedMult_dijet, float mass_dijet, float em_E_frac_dijet,
	       float rho, int nVertices, int nJetInEvt, float DeltaPhi)
{

  if(verbose>2) cout << "--- fill histograms" << endl;
  int outlocal=0;

  // EVENT
  outlocal=LocalFill(MSPlot, "NbOfVertices" , nVertices, dataset, true, Luminosity, scaleFactor);
  outlocal=LocalFill(MSPlot, "RhoCorrection", rho,       dataset, true, Luminosity, scaleFactor);

  // JETS
  outlocal=LocalFill(MSPlot, "NbOfSelectedJets", nJetInEvt, dataset, true, Luminosity, scaleFactor);

  outlocal=LocalFill(MSPlot, "em_E_frac_dijet"     , em_E_frac_dijet,      dataset, true, Luminosity, scaleFactor);
  outlocal=LocalFill(MSPlot, "charged_E_frac_dijet", charged_E_frac_dijet, dataset, true, Luminosity, scaleFactor);
  outlocal=LocalFill(MSPlot, "chargedMult_dijet"   , chargedMult_dijet,    dataset, true, Luminosity, scaleFactor);
  outlocal=LocalFill(MSPlot, "mass_dijet"          , mass_dijet,           dataset, true, Luminosity, scaleFactor);

  string nameJet="";

  for(int iJ=0 ; iJ<m_nJet ; iJ++) { // loop over jets

    if(nJetInEvt<=iJ) break;

    nameJet=Form("_Jet%d",iJ+1);
    //nameJet = nameJet + "_" + name_MSP;

    outlocal=LocalFill(MSPlot, "Eta"+nameJet, jet_eta[iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "E"+  nameJet, jet_E[iJ],   dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "Phi"+nameJet, jet_phi[iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "Pt"+nameJet , jet_pt[iJ],  dataset, true, Luminosity, scaleFactor);

    outlocal=LocalFill(MSPlot, "chargedHad_E_frac"+nameJet, charged_had_E_frac[iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "chargedEm_E_frac" +nameJet, charged_em_E_frac[ iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "chargedMu_E_frac" +nameJet, charged_mu_E_frac[ iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "neutralEm_E_frac" +nameJet, neutral_em_E_frac[ iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "neutralHad_E_frac"+nameJet, neutral_had_E_frac[iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "neutral_E_frac"   +nameJet, tot_neutral_E_frac[iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "charged_E_frac"   +nameJet, tot_charged_E_frac[iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "em_E_frac"        +nameJet, tot_em_E_frac[     iJ], dataset, true, Luminosity, scaleFactor);

    outlocal=LocalFill(MSPlot, "chargedMult"      +nameJet, charged_mult[      iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "muonMult"         +nameJet, muon_mult[         iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "chargedMultTot"   +nameJet, tot_chargedMult[   iJ], dataset, true, Luminosity, scaleFactor);
    outlocal=LocalFill(MSPlot, "neutralMult"      +nameJet, neutral_mult[      iJ], dataset, true, Luminosity, scaleFactor);
  }

  return 0;
}

int LocalFill(MapMSP MSPlot, string name, float value, Dataset* dataset, bool scale, float Luminosity, float scaleFactor)
{
  
  if(MSPlot[name]) {
    MSPlot[name]->Fill(value, dataset, true, Luminosity*scaleFactor);
    return 0;
  }
  else return -1;

}
