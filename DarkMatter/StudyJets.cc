#define _USE_MATH_DEFINES

#include "TStyle.h"
#include "TPaveText.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"
#include "TRandom.h"
#include <iostream>
#include <map>
#include <cstdlib> 

//user code
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

//This header file is taken directly from the BTV wiki. It contains
// to correctly apply an event level Btag SF. It is not yet on CVS
// as I hope to merge the functionality into BTagWeigtTools.h
//#include "TopTreeAnalysisBase/Tools/interface/BTagSFUtil.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"

//#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"

#include "TopTreeAnalysis/macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

int analysis(string, string, string);

int main(int argc, char *argv[])
{

  string outname="tree";
  string xmlfile="configLite.xml";
  string pathpng="Results/";

  if(argc>1) outname=argv[1];
  if(argc>2) xmlfile=argv[2];
  if(argc>3) pathpng=argv[3];
  
  analysis(outname,xmlfile,pathpng);

  return 0;

}

int analysis(string outname, string xmlfile, string pathPNG)
{

  /////////////////////////////
  /// AnalysisEnvironment
  /////////////////////////////
  
  //string pathPNG="Results/";

  AnalysisEnvironment anaEnv;
  cout << "- Loading environment : xmlfile=" << xmlfile << endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());
  cout << "-- env loaded !" << endl;
  int verbose = 2;//anaEnv.Verbose;
  clock_t start = clock();

  // Load datasets
  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << "- Loading datasets" << endl;
  treeLoader.LoadDatasets(datasets, xmlfile.c_str());
  float Luminosity = 20000;

  // PU Reweighting
  LumiReWeighting LumiWeights = LumiReWeighting("PUReweighting/pileup_MC_S10.root", "PUReweighting/PURewtoptree_id_2014_1350565897_PileupHistogram.root", "pileup", "pileup");

  // Apply JES
  int doJESShift = 0; // 0: off 1: minus 2: plus
  cout << "doJESShift: " << doJESShift << endl;

  // Output ROOT file
  string foutname = pathPNG+"/"+outname+".root";
  cout << "foutname=" << foutname << endl;
  TFile *fout = new TFile (foutname.c_str(), "RECREATE");

  // Vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* >   vertex;
  vector < TRootMuon* >     init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* >      init_jets;
  vector < TRootPFJet* >    init_jets_PF;
  vector < TRootMET* >      mets;

  // Global variable
  TRootEvent* event = 0;

  ////////////////////////////
  /// MultiSample plots  /////
  ////////////////////////////

  map<string,MultiSamplePlot*> MSPlot;

  // Event
  MSPlot["RhoCorrection"]         = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
  MSPlot["NbOfVertices"]          = new MultiSamplePlot(datasets, "NbOfVertices", 40, 0, 40, "Nb. of vertices");
  // B-tagging discriminators
  //MSPlot["BdiscBJetCand_CSV"]     = new MultiSamplePlot(datasets, "HighestBdisc_CSV", 75, 0, 1, "CSV b-disc.");
  // Jets
  MSPlot["NbOfSelectedJets"]      = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0, 15, "Nb. of jets");
  MSPlot["NbOfSelectedLightJets"] = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
  MSPlot["NbOfSelectedBJets"]     = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0, 8, "Nb. of jets");
  // Jet Kinematics
  MSPlot["EtaJet1"]               = new MultiSamplePlot(datasets, "EtaJet1", 30,-4, 4, "Leading Jet #eta");
  MSPlot["PhiJet1"]               = new MultiSamplePlot(datasets, "PhiJet1", 50, -4,4 ,"Leading Jet #phi");
  MSPlot["PtJet1"]                = new MultiSamplePlot(datasets, "PtJet1",  200, 0, 1000 ,"Leading Jet p_{T}");
  MSPlot["EtaJet2"]               = new MultiSamplePlot(datasets, "EtaJet2", 30,-4, 4, "Sub-Leading Jet #eta");
  MSPlot["PhiJet2"]               = new MultiSamplePlot(datasets, "PhiJet2", 50, -4,4 ,"Sub-Leading Jet #phi");
  MSPlot["PtJet2"]                = new MultiSamplePlot(datasets, "PtJet2",  200, 0, 1000 ,"Sub-Leading Jet p_{T}");

  ///////////////////
  // Histograms
  ///////////////////

  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;
  
  histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);

  ////////////////////////////
  /// Loop over datasets  ////
  ////////////////////////////
  
  for (unsigned int d = 0; d < datasets.size(); d++) {
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
    
    // LOOP OVER EVENTS
    int itrigger = -1, previousRun = -1;
    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();
    cout <<"Number of events = "<<  end  <<endl;    

    bool debug = true;
    int event_start=0;
    double currentfrac =0.;
    double end_d = end;

    for (unsigned int ievt = event_start; ievt <end_d ; ievt++) { // start loop
      
      double ievt_d = ievt;
      currentfrac = ievt_d/end_d;
      if(ievt%1000 == 0)
	std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);
      float rho = event->kt6PFJets_rho();

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
      

      ////////////////////////////////////////////////////////////////////////////////////
      //  JES Corrections: The nominal corrections are already applied at PAT level     //
      //    so these tools should only be used for studies of the effect of systematics //
      ////////////////////////////////////////////////////////////////////////////////////
      
      // Apply Jet Corrections on-the-fly
      /*
	coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
      	if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
	jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
      	else
	jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
	coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");
      */      
      
      // JER smearing
      /*
	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
        {
	//JER
	doJERShift == 0;
	if(doJERShift == 1)
	jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
	else if(doJERShift == 2)
	jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
	else
	jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
	  
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");
          
	// JES systematic!
	if (doJESShift == 1)
	jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
	else if (doJESShift == 2)
	jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
            
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");
	  
        }
      */

      // Type I MET corrections: (Only for |eta| <=4.7 )
      /*
	coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before MET type I correction:");      
	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
      	jetTools->correctMETTypeOne(init_jets,mets[0],true);
	else
      	jetTools->correctMETTypeOne(init_jets,mets[0],false);
	coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After MET type I correction:");
      */      

      ////////////////////////////////////////
      //  Beam scraping and PU reweighting
      ////////////////////////////////////////

      if(verbose>2) cout << "-- enter beam scraping and PU reweighting" << endl;

      // scale factor for the event
      float scaleFactor = 1.;
      
      if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
	// Apply the scraping veto. (Is it still needed?)
	bool isBeamBG = true;
	if(event->nTracks() > 10) {
	  if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
	    isBeamBG = false;
	}
	if(isBeamBG) continue;
      }
      else{
	double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
	double lumiWeightOLD=lumiWeight;
        if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0){
	  lumiWeight=1;
        }
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

      // Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets);

      // Define object selection cuts
      //selection.setJetCuts(40.,2.5,0.01,1.,0.98,0,0);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets
      selection.setJetCuts(0.,5.,0.,0.,0.,0.,0.);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets //ND
      
      vector<TRootJet*> selectedJets   = selection.GetSelectedJets(true); // ApplyJetId

      vector<TRootJet*>      selectedBJets; // B-Jets
      vector<TRootJet*>      selectedLightJets; // light-Jets
      vector<TRootJet*>      selectedCSVOrderedJets     = selection.GetSelectedJets(true); //CSV ordered Jet collection added by JH

      // ND CHECKS
      if(verbose>2) cout << "-- selectedJets.size()=" << selectedJets.size() << endl;

      // order jets wrt to Pt, then set bool corresponding to RefSel cuts.
      sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt Pt.                                                                    
      //sort(selectedCSVOrderedJets.begin(), selectedCSVOrderedJets.end(), HighestCVSBtag()); //order Jets wrt CSVtag

      // get PFJets from the sorted selected Jets
      vector<TRootPFJet*> selectedJetsPF = jetTools->convertToPFJets(selectedJets);
      
      int JetCut =0;
      //int nMu = selectedMuons.size();
      //int nEl = selectedElectrons.size();

      // Apply primary vertex selection
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
      if(!isGoodPV) continue;

      // Applying baseline offline event selection here
      //if (!(selectedJets.size() >= 2)) continue;

      //////////////////////////////////////
      /// Filling histograms
      //////////////////////////////////////

      // Event
      MSPlot["NbOfVertices"]         ->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
      // Jets
      MSPlot["NbOfSelectedJets"]     ->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfSelectedBJets"]    ->Fill(selectedBJets.size(), datasets[d], true, Luminosity*scaleFactor);
      // Jet Kinematics
      if(selectedJetsPF.size()>0) {
	MSPlot["EtaJet1"]->Fill(selectedJetsPF[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["PhiJet1"]->Fill(selectedJetsPF[0]->Phi(), datasets[d], true, Luminosity*scaleFactor);
	MSPlot["PtJet1"] ->Fill(selectedJetsPF[0]->Pt(),  datasets[d], true, Luminosity*scaleFactor);
	if(selectedJetsPF.size()>1) {
	   MSPlot["EtaJet2"]->Fill(selectedJetsPF[1]->Eta(), datasets[d], true, Luminosity*scaleFactor);
	   MSPlot["PhiJet2"]->Fill(selectedJetsPF[1]->Phi(), datasets[d], true, Luminosity*scaleFactor);
	   MSPlot["PtJet2"] ->Fill(selectedJetsPF[1]->Pt(),  datasets[d], true, Luminosity*scaleFactor);
	}
      }
      //MSPlot["RhoCorrection"]->Fill(event->kt6PFJetsPF2PAT_rho(), datasets[d], true, Luminosity*scaleFactor);
      // Jets kinematics
      //MSPlot["JetEta"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      //MSPlot["JetPhi"]
      
      // plots to to inspire staggered Jet Pt selection
      /*
	if (selectedJets.size()>=4) MSPlot["4thJetPt"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	if (selectedJets.size()>=5) MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	if (selectedJets.size()>=6) MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	if (selectedJets.size()>=7) MSPlot["7thJetPt"]->Fill(selectedJets[6]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      */

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
  //  TString pathPNGJetCombi = pathPNG + "JetCombination/";
  //  mkdir(pathPNGJetCombi,0777);

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
      temp->Draw(name, 0, false, false, false, 0);
      if(verbose>2) cout << "--- drawn!" << endl;

      if(verbose>2) cout << "--- write in pathPNG=" << pathPNG << endl;
      temp->Write(fout, name, true, pathPNG, "png"); //ND true => SaveAs the Canvas as image => seg fault probably caused by empty plots !
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

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;


  return 0;
}
